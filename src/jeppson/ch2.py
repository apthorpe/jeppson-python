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
import logging
import sys

import iapws

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
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout,
                        format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def main(args):
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.info("Starting jeppson_ch2")

    for fh in args.file:
        print("Processing file: {0:s}".format(fh.name))
#    print("The {}-th Fibonacci number is {}".format(args.n, fib(args.n)))

        for ict, rawline in enumerate(fh):
            line = rawline.strip()
            token = line.split()
            if token:
                if token[0].startswith('#'):
                    print("Line {}: [C] {}".format(ict, line))
                else:
                    print("Line {}: [D] {}".format(ict, line))
                    if len(token) >= 6:
                        print("  Enough parts.")
                    else:
                        print("  TOO FEW PARTS.")
            else:
                print("Line {}: [B]".format(ict))

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
