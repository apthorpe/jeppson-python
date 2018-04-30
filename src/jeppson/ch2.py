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
import math
import sys

# import iapws
import scipy.constants as sc
from fluids.friction import friction_factor

from jeppson.pipe import Pipe
from jeppson.input import InputLine
from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"

_logger = logging.getLogger(__name__)


def fib(n):
    """Fibonacci example function

    Args:
      n (int): integer

    Returns:
      int: n-th Fibonacci number
    """
    assert n > 0
    a, b = 1, 1
    for i in range(n-1):
        a, b = b, a+b
    return a


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
    mintok = 6

    FlowInput = namedtuple('FlowInput', ['tag', 'description',
                                         'siunits', 'tunits', 'T2SI_conv'])

#     ikeys = (
#         'pipe inner diameter',
#         'volumetric flowrate',
#         'pipe length',
#         'kinematic viscosity',
#         'absolute pipe roughness',
#         'gravitational acceleration'
#     )
#
#     unitconv = (
#         sc.foot,
#         sc.foot**3,
#         sc.foot,
#         sc.foot**2,
#         sc.foot,
#         sc.foot
#     )

    idataprop = OrderedDict([
        ('idiameter',  FlowInput('idiameter',  'pipe inner diameter',
                                 'm',    'ft',    sc.foot)),
        ('vol_flow',   FlowInput('vol_flow',   'volumetric flowrate',
                                 'm/s',  'ft/s',  sc.foot**3)),
        ('lpipe',      FlowInput('lpipe',      'pipe length',
                                 'm',    'ft',    sc.foot)),
        ('kin_visc',   FlowInput('kin_visc',   'kinematic viscosity',
                                 'm2/s', 'ft/s2', sc.foot**2)),
        ('eroughness', FlowInput('eroughness', 'absolute pipe roughness',
                                 '-',    '-',     sc.foot)),
        ('grav',       FlowInput('grav',       'gravitational acceleration',
                                 'm/s2', 'ft/s2', sc.foot))
    ])

    # Set default results
    results = {
        'status': 'undefined',
        'msg':    'No results generated yet',
        'units':  '',
        'input':  {}
    }

    for kk in idataprop:
        desc = idataprop[kk].description
        results['input'][desc] = float('NaN')

    # Check number of tokens
    if iline.ntok < mintok:
        results['status'] = 'error'
        results['msg'] = 'Too few tokens ({} found, {} expected)' \
                         .format(iline.ntok, mintok)
        return results

    if iline.ntok > mintok:
        results['status'] = 'warning'
        results['msg'] = 'Too many tokens ({} found, {} expected)' \
                         .format(iline.ntok, mintok)
    else:
        results['status'] = 'ok'
        results['msg'] = 'Proper token count ({})'.format(mintok)

    # Convert tokens to floats; check for parsing errors
#    for i, desc in enumerate(ikeys):
    for i, kk in enumerate(idataprop.keys()):
        desc = idataprop[kk].description
        try:
            results['input'][desc] = float(iline.token[i])
        except ValueError as err:
            _logger.error('Numeric parse failure, field {0:d}, "{1:s}": {2:s}'
                          .format(i, desc, str(err)))
            results['status'] = 'error'
            results['msg'] = 'Cannot parse values from input line, field ' \
                             '{0:d}, "{1:s}"'.format(i, desc)
            return results

        if math.isnan(results['input'][desc]):
            _logger.error('Numeric parse failure, field {0:d}, "{1:s}"'
                          .format(i, desc))
            results['status'] = 'error'
            results['msg'] = 'Cannot parse values from input line, field ' \
                             '{0:d}, "{1:s}"'.format(i, desc)
            return results

    # Select units via heuristic (within 10% of SI gravitational acceleration)
    if abs(results['input']['gravitational acceleration'] - sc.g) / sc.g < 0.1:
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

    # Convert units
    if results['units'] == 'Traditional':
#        for i, fconv in enumerate(unitconv):
#            results['input'][ikeys[i]] *= fconv
        for kk in idataprop:
            desc = idataprop[kk].description
            results['input'][desc] *= idataprop[kk].T2SI_conv

    for i, kk in enumerate(idataprop.keys()):
        desc = idataprop[kk].description
        _logger.debug('{0:s} = {1:0.4E} {2:s}'
                      .format(desc, results['input'][desc],
                              idataprop[kk].siunits))

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
            vol_flow (float): Z, in cubic meters per second
            flow_area (float): Z, in square meters
            lpipe (float): Z, in meters
            idiameter (float): Z, in meters
            eroughness (float): Z, dimensionless
            kin_visc (float): Z, in square meters per second
            grav (float): Z, in meters per second squared

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


def generate_results(kwinput):
    """Generate head loss and Darcy-Weisbach friction factor from processed
    input. Intermediate and final results will be returned with metadata

    Args:
      kwinput (dict): A dict with the following keys with values in SI units:
        ``volumetric flowrate``, ``pipe inner diameter``, ``pipe length``,
        ``kinematic viscosity``, ``absolute pipe roughness``,
        and ``gravitational acceleration``

    Returns:
      (dict): Calculation results in SI units and diagnostic info
    """
    dfmt = '  {0:s} = {1:0.4E}'

    ikeys = (
        'volumetric flowrate',
        'pipe inner diameter',
        'pipe length',
        'kinematic viscosity',
        'absolute pipe roughness',
        'gravitational acceleration'
    )

    dkeys = (
        'flow area',
        'relative pipe roughness',
        'flow velocity',
        'reynolds number'
    )

    okeys = ('darcy friction factor', 'head loss')

    _logger.debug('Calculating head loss and friction factor')

    # Set up results structure
    results = {
        'status':  'undefined',
        'msg':     'No results generated yet',
        'input':   {},
        'derived': {},
        'output':  {}
    }

#    _logger.debug('Input object echo')
#     for kk in kwinput:
#         _logger.debug(dfmt.format(kk, kwinput[kk]))

    for kk in dkeys:
        results['derived'][kk] = float('NaN')

    for kk in okeys:
        results['output'][kk] = float('NaN')

    for kk in ikeys:
        if kk in kwinput:
            results['input'][kk] = kwinput[kk]
            _logger.debug(dfmt.format(kk, results['input'][kk]))
        else:
            _logger.debug('Cannot find {0:s} in input'.format(kk))
            # Record the first error
            if results['status'] != 'error':
                results['status'] = 'error'
                results['msg'] = 'Required input "{0:s}" not specified' \
                    .format(kk)

    if results['status'] == 'error':
        return results

    # Alias the parameter dicts
    idata = results['input']
    ddata = results['derived']
    odata = results['output']

    # Arbitrary wall thickness to ensure complete pipe object definition
    twall = 0.1 * idata['pipe inner diameter']
    try:
        pipe = Pipe(label='Example pipe', length=idata['pipe length'],
                    idiameter=idata['pipe inner diameter'],
                    twall=twall, froughness=idata['absolute pipe roughness'])
    except ValueError as err:
        results['status'] = 'error'
        results['msg'] = 'Cannot model pipe: {0:s}'.format(str(err))
        return results

    ddata['flow area'] = pipe.flow_area
    ddata['relative pipe roughness'] = pipe.eroughness

    # Calculate results
    hldat = calculate_headloss(vol_flow=idata['volumetric flowrate'],
                               flow_area=ddata['flow area'],
                               lpipe=idata['pipe length'],
                               idiameter=idata['pipe inner diameter'],
                               eroughness=ddata['relative pipe roughness'],
                               kin_visc=idata['kinematic viscosity'],
                               grav=idata['gravitational acceleration'])

    results['status'] = 'ok'
    results['msg'] = 'Calculation complete'

    ddata['flow velocity'] = hldat.vflow
    ddata['reynolds number'] = hldat.Re

    odata['head loss'] = hldat.head_loss
    odata['darcy friction factor'] = hldat.friction

    # Diagnostics/audit trail
    _logger.debug('Input values:')
    for kk in ikeys:
        _logger.debug(dfmt.format(kk, results['input'][kk]))

    _logger.debug('Derived values:')
    for kk in dkeys:
        _logger.debug(dfmt.format(kk, results['derived'][kk]))

    _logger.debug('Output values:')
    for kk in okeys:
        _logger.debug(dfmt.format(kk, results['output'][kk]))

    _logger.info(results['status'] + ': ' + results['msg'])

    return results


def main(args):
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.info("Starting jeppson_ch2")

    ofmt = 'Q={0:10.4f}D={1:10.4f}L={2:10.2f}F={3:10.5f}HEADLOSS={4:10.4f}'

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
                results = generate_results(indat['input'])

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

# Finally, print results as original code
                if results['status'] in ('ok', 'warning'):
                    if indat['units'] == 'SI':
                        _logger.debug('Display units are SI')
                        print(ofmt.format(
                            results['input']['volumetric flowrate'],
                            results['input']['pipe inner diameter'],
                            results['input']['pipe length'],
                            results['output']['darcy friction factor'],
                            results['output']['head loss']))
                    else:
                        _logger.debug('Display units are Traditional')
                        print(ofmt.format(
                            results['input']['volumetric flowrate']
                            / sc.foot**3,
                            results['input']['pipe inner diameter'] / sc.foot,
                            results['input']['pipe length'] / sc.foot,
                            results['output']['darcy friction factor'],
                            results['output']['head loss'] / sc.foot))
                else:
                    print('! Case defined on line {0:d} of {1:s} failed'
                          .format(ict, fh.name))

    _logger.info("Ending jeppson_ch2")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
