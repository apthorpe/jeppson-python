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
# from collections import namedtuple
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
    mintok = 6

    ikeys = (
        'pipe inner diameter',
        'volumetric flowrate',
        'pipe length',
        'kinematic viscosity',
        'absolute pipe roughness',
        'gravitational acceleration'
    )

    unitconv = (
        sc.foot,
        sc.foot**3,
        sc.foot,
        sc.foot**2,
        sc.foot,
        sc.foot
    )

    # Set default results
    results = {
        'status': 'undefined',
        'msg':    'No results generated yet',
        'units':  '',
        'input':  {}
    }

    for kk in ikeys:
        results['input'][kk] = float('NaN')

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
    for i, kk in enumerate(ikeys):
        try:
            results['input'][kk] = float(iline.token[i])
        except ValueError as err:
            _logger.error('Numeric parse failure, field {0:d}, "{1:s}": {2:s}'
                          .format(i, kk, str(err)))
            results['status'] = 'error'
            results['msg'] = 'Cannot parse values from input line, field ' \
                             '{0:d}, "{1:s}"'.format(i, kk)
            return results

        if math.isnan(results['input'][kk]):
            _logger.error('Numeric parse failure, field {0:d}, "{1:s}"'
                          .format(i, kk))
            results['status'] = 'error'
            results['msg'] = 'Cannot parse values from input line, field ' \
                             '{0:d}, "{1:s}"'.format(i, kk)
            return results

    # Select units via heuristic (within 10% of SI gravitational acceleration)
    if abs(results['input']['gravitational acceleration'] - sc.g) / sc.g < 0.1:
        _logger.info("Assuming SI units")
        results['units'] = 'SI'
    else:
        _logger.info("Assuming traditional units (US/English)")
        results['units'] = 'traditional'

    # Attempt to coerce units
    if force_units:
        if force_units.trim().upper() == 'SI':
            _logger.info("Forcing use of SI units")
            results['units'] = 'SI'
        elif force_units.trim().lower() == 'traditional':
            _logger.info("Forcing use of traditional units")
            results['units'] = 'traditional'
        else:
            msg = 'Cannot force units of measure to "{0:s}"' \
                  .format(force_units)
            _logger.warning(msg)
            if results['status'] == 'ok':
                results['status'] = 'warning'
                results['msg'] = msg
            # pass through; use inferred units

    # Convert units
    if results['units'] == 'traditional':
        for i, fconv in enumerate(unitconv):
            results['input'][ikeys[i]] *= fconv

    for kk in ikeys:
        _logger.debug('{0:s} = {1:0.4E}'.format(kk, results['input'][kk]))

    return results


def calculate_headloss(kwinput):
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
#         _logger.debug('{0:s} is {1:0.4E}'.format(kk, kwinput[kk]))

    for kk in dkeys:
        results['derived'][kk] = float('NaN')

    for kk in okeys:
        results['output'][kk] = float('NaN')

    for kk in ikeys:
        if kk in kwinput:
            results['input'][kk] = kwinput[kk]
            _logger.debug('{0:s} is {1:0.4E}'
                          .format(kk, results['input'][kk]))
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

    # Calculate results
    ddata['flow area'] = pipe.flow_area
    ddata['relative pipe roughness'] = pipe.eroughness
    ddata['flow velocity'] = idata['volumetric flowrate'] / ddata['flow area']
    ddata['reynolds number'] = \
        ddata['flow velocity'] * idata['pipe inner diameter'] \
        / idata['kinematic viscosity']
    odata['darcy friction factor'] = \
        friction_factor(Re=ddata['reynolds number'],
                        eD=ddata['relative pipe roughness'])
    odata['head loss'] = \
        (odata['darcy friction factor'] * idata['pipe length']
         * ddata['flow velocity']**2) \
        / (2.0 * idata['gravitational acceleration']
           * idata['pipe inner diameter'])

    results['status'] = 'ok'
    results['msg'] = 'Calculation complete'

    # Diagnostics/audit trail
    _logger.debug('Input values:')
    for kk in ikeys:
        _logger.debug('{0:s} = {1:0.4E}'
                      .format(kk, results['input'][kk]))

    _logger.debug('Derived values:')
    for kk in dkeys:
        _logger.debug('{0:s} = {1:0.4E}'
                      .format(kk, results['derived'][kk]))

    _logger.debug('Output values:')
    for kk in okeys:
        _logger.debug('{0:s} = {1:0.4E}'
                      .format(kk, results['output'][kk]))

    _logger.debug(results['status'] + ': ' + results['msg'])

    return results


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
        print(msg)
        _logger.info(msg)

        for ict, rawline in enumerate(fh):
            iline = InputLine(line=rawline, ipos=ict+1)
            _logger.debug(iline.as_log())

            if iline.typecode == 'D':
                _logger.debug('Processing data line')
                indat = extract_case_input(iline)
                results = calculate_headloss(indat['input'])

#                results = process_ch2_case(iline.ntok, iline.token)

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

#     do
#         ! 1) Read flow conditions and pipe geometry
#         ! D   - Pipe diameter, ft
#         ! Q   - Flow rate, cfs
#         ! FL  - Length of pipe, ft
#         ! VIS - Kinematic viscosity of fluid (nu)
#         ! E   - Absolute roughness of pipe, ft
#         ! G   - Acceleration of gravity, ft/s**2
#         read(STDIN, 100, end=99) D, Q, FL, VIS, E, G
#
#         ! 2) Calculate friction factor
#         F = f_darcy_weisbach(Q, D, E, VIS, PARLIM, MAXITER, MAXDIF)
#
#         ! 3) Calculate bulk flow velocity
#         V = Q / circ_area(D)
#
#         ! 4) Calculate head loss
#         HL = F * FL * V * V / (2.0 * G * D)
#
#         ! 5) Display results
#         ! Q  - Flow rate, cfs
#         ! D  - Pipe diameter, ft
#         ! FL - Length of pipe, ft
#         ! F  - Darcy-Weisbach friction factor
#         ! HL - Head loss along pipe, ft
#         write(STDOUT, 101) Q, D, FL, F, HL
#
#     end do

    _logger.info("Ending jeppson_ch2")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
