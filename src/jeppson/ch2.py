#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This reimplements the pipe flow analysis program given in chapter 2 of *Steady
Flow Analysis of Pipe Networks: An Instructional Manual* (1974). *Reports.*
Paper 300.  http://digitalcommons.usu.edu/water_rep/300 and *Analysis of Flow
in Pipe Networks* (1976). Ann Arbor Science Publishers, Inc.
http://www.worldcat.org/title/analysis-of-flow-in-pipe-networks/oclc/927534147
by Roland W. Jeppson.

Then run `python setup.py install` which will install the command `jeppson_ch2`
inside your current environment.
"""
from __future__ import division, print_function, absolute_import

import argparse
from collections import namedtuple, OrderedDict
import logging
from math import isnan
import sys

# import iapws
import scipy.constants as sc
from fluids.friction import friction_factor

from . import _logger, ureg, Q_
from jeppson.input import InputLine
from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"

Nomenclature = namedtuple('Nomenclature',
                          ['tag', 'description', 'mks_units', 'us_units'])

idataprop = OrderedDict([
    ('idiameter',  Nomenclature('idiameter',
                                'pipe inner diameter',
                                ureg.meter,
                                ureg.foot)),
    ('vol_flow',   Nomenclature('vol_flow',
                                'volumetric flowrate',
                                ureg.meter**3 / ureg.second,
                                ureg.foot**3 / ureg.second)),
    ('lpipe',      Nomenclature('lpipe',
                                'pipe length',
                                ureg.meter,
                                ureg.foot)),
    ('kin_visc',   Nomenclature('kin_visc',
                                'kinematic viscosity',
                                ureg.meter**2 / ureg.second,
                                ureg.foot**2 / ureg.second)),
    ('froughness', Nomenclature('froughness',
                                'absolute pipe roughness',
                                ureg.meter,
                                ureg.foot)),
    ('grav',       Nomenclature('grav',
                                'gravitational acceleration',
                                ureg.meter/ureg.second**2,
                                ureg.foot/ureg.second**2)),
    ('flow_area',  Nomenclature('flow_area',
                                'pipe flow area',
                                ureg.meter**2,
                                ureg.foot**2)),
    ('vflow',      Nomenclature('vflow',
                                'flow velocity',
                                ureg.meter/ureg.second,
                                ureg.foot/ureg.second)),
    ('Re',         Nomenclature('Re',
                                'reynolds number',
                                ureg.parse_expression('dimensionless'),
                                ureg.parse_expression('dimensionless'))),
    ('eroughness', Nomenclature('eroughness',
                                'relative pipe roughness',
                                ureg.parse_expression('dimensionless'),
                                ureg.parse_expression('dimensionless'))),
    ('head_loss',  Nomenclature('head_loss',
                                'head loss',
                                ureg.meter,
                                ureg.foot)),
    ('friction',   Nomenclature('friction',
                                'darcy-weisback friction factor',
                                ureg.parse_expression('dimensionless'),
                                ureg.parse_expression('dimensionless')))
])

inputconv_us = OrderedDict([
    ('idiameter',  'ft'),
    ('vol_flow',   'ft**3/s'),
    ('lpipe',      'ft'),
    ('kin_visc',   'ft**2/s'),
    ('froughness', 'ft'),
    ('grav',       'ft/s**2')
])

inputconv_mks = OrderedDict([
    ('idiameter',  'm'),
    ('vol_flow',   'm**3/s'),
    ('lpipe',      'm'),
    ('kin_visc',   'm**2/s'),
    ('froughness', 'm'),
    ('grav',       'm/s**2')
])

ikeys = inputconv_us.keys()


def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Frictional head loss calculator")
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
        'input':  {}
    }

    idata = results['input']

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

    tmpdata = {}
    # Convert tokens to floats; check for parsing errors
    _logger.debug('On read:')
    for i, ikey in enumerate(ikeys):
        desc = idataprop[ikey].description
        try:
            tmpdata[ikey] = float(iline.token[i])
            _logger.debug('{0:s} "{1:s}" is {2:0.4E}'
                          .format(ikey, desc, tmpdata[ikey]))
        except ValueError as err:
            _logger.error('Numeric parse failure, '
                          'field {0:d}, {1:s} "{2:s}": {3:s}'
                          .format(i, ikey, desc, str(err)))
            results['status'] = 'error'
            results['msg'] = 'Cannot parse values from input line, field ' \
                             '{0:d}, "{1:s}"'.format(i, ikey)
            return results

        if isnan(tmpdata[ikey]):
            _logger.error('Numeric parse failure, field {0:d}, {1:s}, "{2:s}"'
                          .format(i, ikey, desc))
            results['status'] = 'error'
            results['msg'] = 'Cannot parse values from input line, field ' \
                '{0:d}, {1:s}, "{2:s}"'.format(i, ikey, desc)
            return results

    # Select units via heuristic (within 10% of SI gravitational acceleration)
    if abs(tmpdata['grav'] - sc.g) / sc.g < 0.1:
        _logger.info("Assuming SI units")
        results['units'] = 'SI'
    else:
        _logger.info("Assuming traditional units (US/English)")
        results['units'] = 'Traditional'

    # Attempt to coerce units
    if force_units:
        if force_units.trim().upper() == 'SI':
            _logger.info("Forcing use of SI units")
            results['units'] = 'SI'
        elif force_units.trim().lower() == 'traditional':
            _logger.info("Forcing use of traditional units")
            results['units'] = 'Traditional'
        else:
            msg = 'Cannot force units of measure to "{0:s}"' \
                  .format(force_units)
            _logger.warning(msg)
            if results['status'] == 'ok':
                results['status'] = 'warning'
                results['msg'] = msg
            # pass through; use inferred units

    # Set units (Pint)
    _logger.debug('As Pint quantities:')
    for ikey in ikeys:
        desc = idataprop[ikey].description
        if results['units'] == 'SI':
            idata[ikey] = Q_(tmpdata[ikey], idataprop[ikey].mks_units)
        else:
            idata[ikey] = Q_(tmpdata[ikey], idataprop[ikey].us_units)

        _logger.debug('  {0:s} = {1:0.4E~P}'
                      .format(desc, idata[ikey].to_base_units()))

    _logger.info(results['status'] + ': ' + results['msg'])

    return results


def calculate_headloss(vol_flow, flow_area, lpipe, idiameter, eroughness,
                       kin_visc, grav):
    """Calculate Darcy-Weisbach friction factor and head loss

       Calculates head loss (m) and Darcy-Wesibach friction factor and
       intermediate quantities flow velocity (m/s), and Reynolds number. These
       values are returned in a named tuple with key names ``head_loss``,
       ``friction``, ``vflow``, and ``Re`` respectively.

    Args:
        vol_flow (float): volumetric flow, in cubic meters per second
        flow_area (float): pipe flow area, in square meters
        lpipe (float): pipe length, in meters
        idiameter (float): pipe inner diameter, in meters
        eroughness (float): pipe relative roughness, dimensionless
        kin_visc (float): kinematic viscosity, in square meters per second
        grav (float): gravitational acceleration, in meters per second squared

    Returns:
        (namedtuple): Results and intermediate quantities

    """
    HeadLoss = namedtuple('HeadLoss', ['head_loss', 'friction', 'vflow', 'Re'])

    flow_vel = vol_flow / flow_area
    Re = flow_vel * idiameter / kin_visc
    friction = friction_factor(Re=Re, eD=eroughness)
    head_loss = (friction * lpipe * flow_vel**2) / (2.0 * grav * idiameter)

    hldat = HeadLoss(head_loss, friction, flow_vel, Re)

    return hldat


def generate_results(idata):
    """Generate head loss and Darcy-Weisbach friction factor from processed
    input. Intermediate and final results will be returned with metadata

    Args:
        idata (dict): A dict with the following keys with values as Pint
          objects: ``vflow`` (volumetric flowrate), ``idiameter`` (pipe inner
          diameter), ``lpipe`` (pipe length), ``kin_visc`` (kinematic
          viscosity), ``froughness`` (absolute pipe roughness), and ``grav``
          (gravitational acceleration)

    Returns:
        (dict): Calculation results and diagnostic info
    """
    dfmt = '  {0:s} = {1:0.4E~}'
#    dufmt = '  {0:s} = {1:0.4E~P}'

    dkeys = (
        'flow_area',
        'eroughness',
        'vflow',
        'Re'
    )

    okeys = (
        'friction',
        'head_loss'
    )

    _logger.debug('Calculating head loss and friction factor')

    # Set up results structure
    results = {
        'status':  'undefined',
        'msg':     'No results generated yet',
        'derived': {},
        'output':  {}
    }

    # Alias the parameter dicts
    ddata = results['derived']
    odata = results['output']

#    _logger.debug('Input object echo')
#     for kk in idata:
#         _logger.debug(dfmt.format(kk, idata[kk]))

    for kk in dkeys:
        ddata[kk] = float('NaN')

    for kk in okeys:
        odata[kk] = float('NaN')

    for kk in ikeys:
        if kk not in idata:
            _logger.debug('Cannot find {0:s} in input'.format(kk))
            # Record the first error
            if results['status'] != 'error':
                results['status'] = 'error'
                results['msg'] = 'Required input "{0:s}" not specified' \
                    .format(kk)

    if results['status'] == 'error':
        return results

    ddata['flow_area'] = sc.pi / 4.0 * idata['idiameter']**2
    ddata['eroughness'] = idata['froughness'] / idata['idiameter']

    # Calculate results
    hldat = calculate_headloss(
        vol_flow=idata['vol_flow'].to('m**3/s').magnitude,
        flow_area=ddata['flow_area'].to('m**2').magnitude,
        lpipe=idata['lpipe'].to('m').magnitude,
        idiameter=idata['idiameter'].to('m').magnitude,
        eroughness=ddata['eroughness'].magnitude,
        kin_visc=idata['kin_visc'].to('m**2/s').magnitude,
        grav=idata['grav'].to('m/s**2').magnitude)

    ddata['vflow'] = Q_(hldat.vflow, 'm/s')
    ddata['Re'] = Q_(hldat.Re, '')

    odata['head_loss'] = Q_(hldat.head_loss, 'm')
    odata['friction'] = Q_(hldat.friction, '')

    results['status'] = 'ok'
    results['msg'] = 'Calculation complete'

    # Diagnostics/audit trail
    _logger.debug('Input values:')
    for kk in ikeys:
        _logger.debug(dfmt.format(kk, idata[kk]))

    _logger.debug('Derived values:')
    for kk in dkeys:
        _logger.debug(dfmt.format(kk, ddata[kk]))

    _logger.debug('Output values:')
    for kk in okeys:
        _logger.debug(dfmt.format(kk, odata[kk]))

    _logger.info(results['status'] + ': ' + results['msg'])

    return results


def generate_legacy_output(idata, odata, units):
    """Return a head loss, friction factor, and initial conditions in the
    format of the original Jeppson code

    Args:
        idata (dict): Input to head loss calculations
        odata (dict): Results of head loss calculations
        units (str): Display units for reporting results - either 'SI'
          or 'Traditional'

    Returns:
        (str): FORTRAN-formatted head loss results
    """
    unitmap = {
        'vol_flow': {'SI': 'm**3/s', 'Traditional': 'ft**3/s'},
        'idiameter': {'SI': 'm', 'Traditional': 'ft'},
        'lpipe': {'SI': 'm', 'Traditional': 'ft'},
        'friction': {'SI': '', 'Traditional': ''},
        'head_loss': {'SI': 'm', 'Traditional': 'ft'}
    }

    if units in ('SI', 'Traditional'):
        dunits = units
    else:
        dunits = 'SI'

    ofmt = 'Q={0:10.4f}D={1:10.4f}L={2:10.2f}F={3:10.5f}HEADLOSS={4:10.4f}'

    _logger.debug('Display units are {0:s}'.format(units))
    outstr = ofmt.format(
        idata['vol_flow'].to(unitmap['vol_flow'][dunits]).magnitude,
        idata['idiameter'].to(unitmap['idiameter'][dunits]).magnitude,
        idata['lpipe'].to(unitmap['lpipe'][dunits]).magnitude,
        odata['friction'].magnitude,
        odata['head_loss'].to(unitmap['head_loss'][dunits]).magnitude)

    return outstr


def main(args):
    """Main entry point allowing external calls

    Args:
        args ([str]): command line parameter list
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.info("Starting jeppson_ch2")

    for fh in args.file:
        msg = 'Processing file: {0:s}'.format(fh.name)
        _logger.info(msg)

        for ict, rawline in enumerate(fh):
            iline = InputLine(line=rawline, ipos=ict+1)
            _logger.debug(iline.as_log())

            if iline.typecode == 'D':
                _logger.info('Processing data line:')
                _logger.info(iline.as_log())

                indat = extract_case_input(iline)
                idata = indat['input']
                results = generate_results(idata)
#                results = generate_results(indat['input'], indat['uinput'])

                # Log calcularion status;
                if results['status'] == 'ok':
                    _logger.debug('Case processed successfully')
                elif results['status'] == 'warning':
                    _logger.warning('Case processed with warning: {0:s}'
                                    .format(results['msg']))
                elif results['status'] == 'error':
                    _logger.error('Case not processed due to input '
                                  'error: {0:s}'.format(results['msg']))
                else:
                    _logger.error('Unknown status processing case: '
                                  '"{0:s}", {1:s}'
                                  .format(results['status'], results['msg']))

                # Finally, print results as original code; infer output units
                # from input units
                if results['status'] in ('ok', 'warning'):
                    ostr = generate_legacy_output(idata,
                                                  results['output'],
                                                  indat['units'])
                    print(ostr)
                else:
                    # Differs from original code; broken input would crash and
                    # all input was provided via STDIN. Here multiple files may
                    # be procssed, each potentially containing multiple cases.
                    # Report errors with line and source file name
                    print('! Case defined on line {0:d} of {1:s} failed'
                          .format(ict, fh.name))

    _logger.info("Ending jeppson_ch2")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
