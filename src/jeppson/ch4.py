#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This reimplements the pipe flow analysis program given in chapter 4 of *Steady
Flow Analysis of Pipe Networks: An Instructional Manual* (1974). *Reports.*
Paper 300.  http://digitalcommons.usu.edu/water_rep/300 and *Analysis of Flow
in Pipe Networks* (1976). Ann Arbor Science Publishers, Inc.
http://www.worldcat.org/title/analysis-of-flow-in-pipe-networks/oclc/927534147
by Roland W. Jeppson.

Then run `python setup.py install` which will install the command `jeppson_ch4`
inside your current environment.
"""
from __future__ import division, print_function, absolute_import

import argparse
from collections import namedtuple, OrderedDict
import logging
from math import isnan, log
import sys

# import iapws
import scipy.constants as sc
from fluids.friction import friction_factor

from . import _logger, ureg, Q_
# from jeppson.pipe import Pipe
from jeppson.input import InputLine
from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"

Nomenclature = namedtuple('Nomenclature',
                          ['tag', 'description', 'mks_units', 'us_units'])

idataprop = OrderedDict([
    ('vol_flow',   Nomenclature('vol_flow',
                                'volumetric flowrate',
                                ureg.meter**3 / ureg.second,
                                ureg.foot**3 / ureg.second)),
    ('idiameter',  Nomenclature('idiameter',
                                'pipe inner diameter',
                                ureg.meter,
                                ureg.foot)),
    ('lpipe',      Nomenclature('lpipe',
                                'pipe length',
                                ureg.meter,
                                ureg.foot)),
    ('froughness', Nomenclature('froughness',
                                'absolute pipe roughness',
                                ureg.meter,
                                ureg.foot)),
    ('fvol_flow',  Nomenclature('fvol_flow',
                                'fractional volumetric flow',
                                ureg.parse_expression('dimensionless'),
                                ureg.parse_expression('dimensionless'))),
    ('kin_visc',   Nomenclature('kin_visc',
                                'kinematic viscosity',
                                ureg.meter**2 / ureg.second,
                                ureg.foot**2 / ureg.second)),
    ('flow_area',  Nomenclature('flow_area',
                                'pipe flow area',
                                ureg.meter**2,
                                ureg.foot**2)),
    ('arl',        Nomenclature('arl',
                                'power law model constant term',
                                ureg.second**2 / ureg.meter**5,
                                ureg.second**2 / ureg.feet**5
                                )),
    ('dvol_flow',  Nomenclature('dvol_flow',
                                'volumentric flowrate deviation',
                                ureg.meter**3 / ureg.second,
                                ureg.foot**3 / ureg.second)),
    ('eroughness', Nomenclature('eroughness',
                                'relative pipe roughness',
                                ureg.parse_expression('dimensionless'),
                                ureg.parse_expression('dimensionless'))),
    ('friction',   Nomenclature('friction',
                                'darcy-weisback friction factor',
                                ureg.parse_expression('dimensionless'),
                                ureg.parse_expression('dimensionless')))
])

inputconv = OrderedDict([
    ('vol_flow',   'ft**3/s'),
    ('idiameter',  'in'),
    ('lpipe',      'ft'),
    ('froughness', 'in'),
    ('fvol_flow',  ''),
    ('kin_visc',   'ft**2/s')])

ikeys = inputconv.keys()


def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Incompressible flow calculator")
    parser.add_argument(
        '--version',
        action='version',
        version='jeppson-python {ver}'.format(ver=__version__))
    parser.add_argument(
        dest="file",
        help="input files (STDIN if not specified)",
        type=argparse.FileType('r'),
        nargs='*',
        default=[sys.stdin],
        metavar="FILE")
    parser.add_argument(
        '-v',
        '--verbose',
        dest="loglevel",
        help="set loglevel to INFO",
        action='store_const',
        const=logging.INFO)
    parser.add_argument(
        '-vv',
        '--very-verbose',
        dest="loglevel",
        help="set loglevel to DEBUG",
        action='store_const',
        const=logging.DEBUG)
    return parser.parse_args(args)


def setup_logging(loglevel):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s: %(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout,
                        format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def extract_case_input(iline, force_units=None):
    """Extract and validate case input from pre-processed data line

        Args:
            iline (InputLine): Pre-processed input data line
            force_units (str): Unit system of input data. Accepted values are
              'SI' (mks) and 'Traditional' (foot-pound-second). If left
              undefined, units will be inferred from gravitational
              acceleration.

        Returns:
            (dict): Case input data and metadata"""
    _logger.debug('Extracting input data for single case')

    # Set default results
    results = {
        'status': 'undefined',
        'msg':    'No results generated yet',
        'idata':  {}
    }

    idata = results['idata']

    for ikey in ikeys:
        idata[ikey] = float('NaN')

    # Check number of tokens
    mintok = len(ikeys)
    if iline.ntok < mintok:
        results['status'] = 'error'
        results['msg'] = 'Too few tokens ({} found, {} expected)' \
                         .format(iline.ntok, mintok)
        return results

    elif iline.ntok > mintok:
        results['status'] = 'warning'
        results['msg'] = 'Too many tokens ({} found, {} expected)' \
                         .format(iline.ntok, mintok)
    else:
        results['status'] = 'ok'
        results['msg'] = 'Proper token count ({})'.format(mintok)

    # Convert tokens to floats; check for parsing errors
    _logger.debug('On read:')
    for i, ikey in enumerate(ikeys):
        desc = idataprop[ikey].description
        try:
            idata[ikey] = Q_(float(iline.token[i]), inputconv[ikey])
            _logger.debug('{0:s} "{1:s}" is {2:0.4E~}'
                          .format(ikey, desc,
                                  idata[ikey].to(inputconv[ikey])))
        except ValueError as err:
            _logger.error('Numeric parse failure, '
                          'field {0:d}, {1:s} "{2:s}": {3:s}'
                          .format(i, ikey, desc, str(err)))
            results['status'] = 'error'
            results['msg'] = 'Cannot parse values from input line, field ' \
                             '{0:d}, "{1:s}"'.format(i, ikey)
            return results

        if isnan(idata[ikey].magnitude):
            _logger.error('Numeric parse failure, field {0:d}, {1:s}, "{2:s}"'
                          .format(i, ikey, desc))
            results['status'] = 'error'
            results['msg'] = 'Cannot parse values from input line, field ' \
                '{0:d}, {1:s}, "{2:s}"'.format(i, ikey, desc)
            return results

    # Set units (Pint)
    _logger.debug('As Pint quantities:')
    for ikey in ikeys:
        desc = idataprop[ikey].description
        _logger.debug('  {0:s} = {1:0.4E~P}'
                      .format(desc, idata[ikey].to_base_units()))

    _logger.info(results['status'] + ': ' + results['msg'])

    return results


def pipe_friction(vol_flow, idiameter, kin_visc, flow_area, eroughness):
    """Calculate friction factor from """
    vflow = vol_flow / flow_area
    Re = vflow * idiameter / kin_visc
    friction = friction_factor(Re=Re, eD=eroughness)
    _logger.debug('{0:16s}{1:12.4E}'.format('vflow', vflow))
    _logger.debug('{0:16s}{1:12.4E}'.format('Re', Re))

    return friction


def generate_results(idata):
    """ Generate intermediate quantities and modeling
    coeffcients for incompressible pipe flow cases"""
    ugrav = Q_(sc.g, "m/s**2")

    # Results
    odata = {}

    ddata = {}
    ddata['flow_area'] = sc.pi / 4.0 * idata['idiameter']**2
    ddata['dvol_flow'] = idata['fvol_flow'] * idata['vol_flow']
    ddata['eroughness'] = idata['froughness'] / idata['idiameter']
    ddata['arl'] = idata['lpipe'] \
        / (2.0 * ugrav * ddata['flow_area']**2 * idata['idiameter'])

    # Logging
    _logger.debug('Input data:')
    for tag in ikeys:
        _logger.debug('  {0:s} is {1:0.4E~}'
                      .format(tag, idata[tag].to(idataprop[tag].us_units)))

    _logger.debug('Derived data:')
    for tag in ddata:
        _logger.debug('  {0:s} is {1:0.4E~}'
                      .format(tag, ddata[tag].to(idataprop[tag].us_units)))

    vol_flow1 = idata['vol_flow'] - ddata['dvol_flow']
    friction1 = Q_(pipe_friction(vol_flow1.to('m**3/s').magnitude,
                                 idata['idiameter'].to('m').magnitude,
                                 idata['kin_visc'].to('m**2/s').magnitude,
                                 ddata['flow_area'].to('m**2').magnitude,
                                 ddata['eroughness'].to('').magnitude), '')

    vol_flow2 = idata['vol_flow'] + ddata['dvol_flow']
    friction2 = Q_(pipe_friction(vol_flow2.to('m**3/s').magnitude,
                                 idata['idiameter'].to('m').magnitude,
                                 idata['kin_visc'].to('m**2/s').magnitude,
                                 ddata['flow_area'].to('m**2').magnitude,
                                 ddata['eroughness'].to('').magnitude), '')

    # Logging
    _logger.debug('{0:16s}{1:12.4E~} {2:12.4E~}'
                  .format('vol_flow', vol_flow1, vol_flow2))
    _logger.debug('{0:16s}{1:12.4E~}             {2:12.4E~}'
                  .format('friction', friction1, friction2))

    odata['be'] = ((log(friction1.magnitude) - log(friction2.magnitude))
                   / (log(vol_flow2.magnitude) - log(vol_flow1.magnitude)))
    odata['ae'] = friction1.magnitude \
        * (vol_flow1.to('ft**3/s').magnitude)**odata['be']
    odata['ep'] = 2.0 - odata['be']
    odata['ck'] = odata['ae'] * ddata['arl'].to('s**2/ft**5').magnitude

    # Logging
    _logger.debug('Output data:')
    for ikey in ('ae', 'be', 'ck', 'ep'):
        _logger.debug('{0:s} = {1:0.4E}'.format(ikey, odata[ikey]))

    return ddata, odata


def generate_legacy_output(idata, odata):
    """Return a string containing volumetric flowrate, pipe diameter,
    and flow correlation parameters in the format of the original Jeppson
    code.

    Args:
        idata (dict): Input to incompressible flow calculations (values are
          Pint Quantity objects)
        odata (dict): Correlation parameters (values are floats)

    Returns:
        (str): FORTRAN-formatted head loss results
    """

    outstr = ''

    outstr = '{:12.5f}{:12.5f}{:12.5f}{:12.5f}{:12.5f}{:16.6E}' \
             .format(idata['vol_flow'].to('ft**3/s').magnitude,
                     idata['idiameter'].to('ft').magnitude,
                     odata['be'], odata['ae'],
                     odata['ep'], odata['ck'])

    return outstr


def main(args):
    """Main entry point allowing external calls

    Args:
        args ([str]): command line parameter list
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.info("Starting jeppson_ch4")

    for fh in args.file:
        msg = 'Processing file: {0:s}'.format(fh.name)
        _logger.info(msg)

        for ict, rawline in enumerate(fh):
            iline = InputLine(line=rawline, ipos=ict+1)
            _logger.debug(iline.as_log())

            if iline.typecode == 'D':
                _logger.info('Processing data line:')
                _logger.info(iline.as_log())

                results = extract_case_input(iline)

                # Log parse status;
                if results['status'] == 'ok':
                    _logger.debug('Case processed successfully')
                elif results['status'] == 'warning':
                    _logger.warning('Case processed with warning: {0:s}'
                                    .format(results['msg']))
                elif results['status'] == 'error':
                    _logger.error('Case not processed due to input '
                                  'error: {0:s}'.format(results['msg']))
                    _logger.debug('Skipping to next case')
                    continue
                else:
                    _logger.error('Unknown status processing case: '
                                  '"{0:s}", {1:s}'
                                  .format(results['status'], results['msg']))
                    _logger.debug('Skipping to next case')
                    continue

                idata = results['idata']

                ddata, odata = generate_results(idata)

                ostr = generate_legacy_output(idata, odata)

                print(ostr)

    _logger.info("Ending jeppson_ch4")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
