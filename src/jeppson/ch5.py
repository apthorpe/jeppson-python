#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This reimplements the pipe flow analysis program given in chapter 5 of *Steady
Flow Analysis of Pipe Networks: An Instructional Manual* (1974). *Reports.*
Paper 300.  http://digitalcommons.usu.edu/water_rep/300 and *Analysis of Flow
in Pipe Networks* (1976). Ann Arbor Science Publishers, Inc.
http://www.worldcat.org/title/analysis-of-flow-in-pipe-networks/oclc/927534147
by Roland W. Jeppson.

Then run `python setup.py install` which will install the command `jeppson_ch5`
inside your current environment.
"""
from __future__ import division, print_function, absolute_import

import argparse
from collections import namedtuple, OrderedDict
import logging
# from math import isnan, log
import sys

# import iapws
# import scipy.constants as sc
# from fluids.friction import friction_factor

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
        description="Linear method flow solver")
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


def extract_case_parameters(deck, iptr):
    """Parse out network characteristics

    Args:
        deck [(InputLine)]: List of parsed input data lines
        iptr (int): Starting pointer for reading case parameters

    Returns:
        (dict): Case parameters
    """
    tags = ('npipes', 'njunctions', 'nloops', 'maxiter', 'unitcode',
            'errlimit', 'kin_visc', 'fvol_flow')
    mintok = len(tags)
    CaseParameters = namedtuple("CaseParameters", tags)

    results = {
        'status': 'unknown',
        'msg': 'CaseParameters are not initialized',
        'iread': 1
    }

    cpline = deck[iptr]
    if cpline.ntok < mintok:
        results['status'] = 'error'
        results['msg'] = 'Too few entries found for case parameters'
        results['case_parameters'] = ()
    elif cpline.ntok > mintok:
        results['status'] = 'warning'
        results['msg'] = 'More entries found than expected for case parameters'
    else:
        results['status'] = 'ok'
        results['msg'] = 'warning'

    if results['status'] in ('ok', 'warning'):
        npipes = int(cpline.token[0])
        njunctions = int(cpline.token[1])
        nloops = int(cpline.token[2])
        maxiter = int(cpline.token[3])
        unitcode = int(cpline.token[4])
        errlimit = float(cpline.token[5])
        if unitcode <= 1:
            kinvisc = Q_(float(cpline.token[6]), 'ft**2/s')
        else:
            kinvisc = Q_(float(cpline.token[6]), 'm**2/s')
        fvol_flow = float(cpline.token[7])

        results['case_parameters'] = CaseParameters(
            npipes, njunctions, nloops, maxiter,
            unitcode, errlimit, kinvisc, fvol_flow)

    return results


def main(args):
    """Main entry point allowing external calls

    Args:
        args ([str]): command line parameter list
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.info("Starting jeppson_ch5")

    for fh in args.file:
        msg = 'Processing file: {0:s}'.format(fh.name)
        _logger.info(msg)

        deck = []

        for ict, rawline in enumerate(fh):
            iline = InputLine(line=rawline, ipos=ict+1)
            _logger.debug(iline.as_log())

            if iline.typecode == 'D':
                deck.append(iline)

# Step 1. Extract case parameters
        iptr = 0
        case_info = extract_case_parameters(deck, iptr)
        if case_info['status'] == 'error':
            _logger.error('Error reading case parameters: {0:s}'
                          .format(case_info['msg']))
            _logger.info('Skipping remainder of {0:s}'.format(fh.file))
            continue

        iptr += case_info['iread']

# Step 2. Read pipe data
#        pipe = extract_pipe_definitions(deck, iptr)
#        iptr += len(pipe)
# Step 3. Read junction data
#        pipemap, iread = extract_junctions(deck, iptr, pipe)
#        iptr += iread
# Step 4. Read loop data
#        pipeloop, iread = extract_loops(deck, iptr, pipe)

# Done reading input; assert (iptr + iread - 1) == len(deck)
# Step 5. Assemble matrix
# Step 6. Solve matrix
# Step 7. Display results
    _logger.info("Ending jeppson_ch5")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
