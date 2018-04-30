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
from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"

_logger = logging.getLogger(__name__)


class InputLine(dict):
    """ Categorized and tokenized user input for Jeppson Ch. 2 friction factor
        and head loss calculator.

        All attributes are immutable except ``ipos``.

        Attributes:
            line (str): Original line of input text with newline(s) removed
            type (str): One of 'blank', 'comment', or 'data'
            typecode (str): One of 'B', 'C', or 'D', corresponding to type
            ipos (int): Line number in original file (count starts at 1)
            ntok (int): Number of tokens found (0 except for 'data' lines)
            token ([str]): Tokens parsed from line ('data' lines only,
              otherwise empty)

        Args:
            line (str): Original line of input text with newline(s) removed
            ipos (int): Line number in original file
            commentchar (str): leading character/string denoting that a line is
              a comment
    """

    def __init__(self, line, ipos=0, commentchar='#'):
        """Constructor """
        self.ipos = ipos

        self._line = line.rstrip('\r\n')
        self._type = 'unknown'
        self._ntok = 0
        self._token = ()

        tline = self._line.strip()
        ltline = len(tline)

        if ltline == 0:
            self._type = 'blank'
        else:
            if commentchar and self._line.startswith(commentchar):
                self._type = 'comment'
            else:
                self._type = 'data'
                self._token = tline.split()
                self._ntok = len(self.token)

    @property
    def line(self):
        """Line read accessor

            Returns:
                (str): Original input line stripped of line terminators
        """
        return self._line

#    @line.setter
#    def line(self, line):
#        """Line write accessor
#            Raises:
#                ValueError
#        """
#        return self._line

    @property
    def type(self):
        """Type read accessor

            Returns:
                (str): Type of input line; one of 'blank', 'comment', or 'data'
        """
        return self._type

    @property
    def typecode(self):
        """Type code read accessor

            Returns:
                (str): Type code of input line; one of 'B', 'C', or 'D',
                  corresponding to 'blank', 'comment', or 'data', respectively
        """
        return self._type[0].upper()

    @property
    def ntok(self):
        """Token count read accessor

            Returns:
                (int): Number of tokens found (0 except for 'data' lines)
        """
        return self._ntok

    @property
    def token(self):
        """Token list read accessor

            Returns:
                ([str]): Tokens parsed from line ('data' lines only, otherwise
                  empty)
        """
        return self._token

    def as_log(self, logfmt='{0:-6d} [{1:s}] {2:s}'):
        """Return line number (ipos), type code, and original line to assist
        in finding input errors. Type code is 'B', 'C', or 'D', corresponding
        to 'blank', 'comment', and 'data', respectively

            Args:
                logfmt (str): Format string for producing log output. Field 0
                  is the `ipos` attribute, field 1 is the type code, and field
                  2 is the `line` attribute

            Returns:
                (str): Formatted line with metadata
        """
        return logfmt.format(self.ipos, self.typecode, self._line)


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


# def process_ch2_case(ntok, token):
#     """Generate results from tokenized data
#
#     Args:
#       ntok (int): number of data tokens
#       token ([str]): tokens parsed from data line.
#         Contains at least one token
#
#     Returns:
#       (dict): Calculation results and diagnostic info
#     """
#     mintok = 6
#
#     ikeys = (
#         'volumetric flowrate', 'pipe inner diameter', 'pipe length',
#         'kinematic viscosity', 'absolute pipe roughness',
#         'gravitational acceleration', 'flow velocity', 'reynolds number'
#     )
#     okeys = (
#         'darcy friction factor', 'head loss'
#     )
#
#     results = {
#         'status': 'undefined',
#         'msg':    'No results generated yet',
#         'units':  '',
#         'input':  {},
#         'output': {}
#     }
#
#     for kk in ikeys:
#         results['input'][kk] = float('NaN')
#
#     for kk in okeys:
#         results['output'][kk] = float('NaN')
#
#     if ntok < mintok:
#         results['status'] = 'error'
#         results['msg'] = 'Too few tokens ({} found, {} expected)' \
#                          .format(ntok, mintok)
#         return results
#
#     if ntok > mintok:
#         results['status'] = 'warning'
#         results['msg'] = 'Too many tokens ({} found, {} expected)' \
#                          .format(ntok, mintok)
#     else:
#         results['status'] = 'ok'
#         results['msg'] = 'Proper token count ({})'.format(mintok)
#
#     try:
#         # D   - Pipe diameter, ft
#         idiameter = float(token[0])
#         # Q   - Flow rate, cfs
#         vol_flowrate = float(token[1])
#         # FL  - Length of pipe, ft
#         length = float(token[2])
#         # VIS - Kinematic viscosity of fluid (nu)
#         kin_visc = float(token[3])
#         # E   - Absolute roughness of pipe, ft
#         froughness = float(token[4])
#         # G   - Acceleration of gravity, ft/s**2
#         agrav = float(token[5])
#     except ValueError as err:
#         _logger.error("Numeric parse failure: {0:s}".format(str(err)))
#         results['status'] = 'error'
#         results['msg'] = 'Cannot parse values from input line'
#         return results
#
#     if abs(agrav - sc.g) / sc.g < 0.1:
#         _logger.info("Assuming SI units")
#         results['units'] = 'SI'
#     else:
#         _logger.info("Assuming traditional units (US/English)")
#         results['units'] = 'traditional'
#         idiameter *= sc.foot
#         vol_flowrate *= sc.foot**3
#         length *= sc.foot
#         kin_visc *= sc.foot**2
#         froughness *= sc.foot
#         agrav *= sc.foot
#
#     _logger.debug('idiameter = {0:0.4E} m'.format(idiameter))
#     _logger.debug('vol_flowrate = {0:0.4E} m3/s'.format(vol_flowrate))
#     _logger.debug('length = {0:0.4E} m'.format(length))
#     _logger.debug('kin_visc = {0:0.4E} m2/s'.format(kin_visc))
#     _logger.debug('froughness = {0:0.4E} m'.format(froughness))
#     _logger.debug('agrav = {0:0.4E} m/s2'.format(agrav))
#
#     twall = 0.1 * idiameter
#     try:
#         pipe = Pipe(label='Example pipe', length=length, idiameter=idiameter,
#                     twall=twall, froughness=froughness)
#     except ValueError as err:
#         results['status'] = 'error'
#         results['msg'] = 'Cannot model pipe: {0:s}'.format(str(err))
#         return results
#
#     # Calculate results
#     vflow = vol_flowrate / pipe.flow_area
#     Re = vflow * pipe.idiameter / kin_visc
#     friction = friction_factor(Re, eD=pipe.eroughness)
#     head_loss = friction * pipe.length * vflow**2 \
#         / (2.0 * agrav * pipe.idiameter)
#
#     # Package results
#     results['input']['volumetric flowrate'] = vol_flowrate
#     results['input']['pipe inner diameter'] = pipe.idiameter
#     results['input']['pipe length'] = pipe.length
#     results['input']['kinematic viscosity'] = kin_visc
#     results['input']['absolute pipe roughness'] = froughness
#     results['input']['gravitational acceleration'] = agrav
#     results['input']['flow velocity'] = vflow
#     results['input']['reynolds number'] = Re
#     results['output']['darcy friction factor'] = friction
#     results['output']['head loss'] = head_loss
#
#     for kk in ikeys:
#         _logger.debug('{0:s} = {1:0.4E}'
#                       .format(kk, results['input'][kk]))
#
#     for kk in okeys:
#         _logger.debug('{0:s} = {1:0.4E}'
#                       .format(kk, results['output'][kk]))
#
#     return results


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
#            line = rawline.rstrip('\r\n')
#            ltype, ntok, token = categorize_line(line)
#
#            lcode = ltype[0].upper()
#            msg = '{0:-6d} [{1:s}] {2:s}' \
#                  .format(ict+1, ltype[0].upper(), line)
#            _logger.debug(msg)
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
