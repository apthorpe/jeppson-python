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
from math import copysign, log
import sys

# import iapws
import scipy.constants as sc
# from fluids.core import Reynolds
from fluids.friction import friction_factor
import numpy as np

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

ugrav = Q_(sc.g, 'm/s**2')

pipe_inputconv = [
    {'idiameter': 'in', 'lpipe': 'ft', 'froughness': 'in'},
    {'idiameter': 'ft', 'lpipe': 'ft', 'froughness': 'ft'},
    {'idiameter': 'm', 'lpipe': 'm', 'froughness': 'm'},
    {'idiameter': 'cm', 'lpipe': 'm', 'froughness': 'cm'}
]

flow_inputconv = [
    'skip',
    'gallon / minute',
    'ft**3 / sec',
    'm**3 / sec'
]


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
            'tolerance', 'kin_visc', 'fvol_flow')
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
        tolerance = float(cpline.token[5])
        if unitcode <= 1:
            kinvisc = Q_(float(cpline.token[6]), 'ft**2/s')
        else:
            kinvisc = Q_(float(cpline.token[6]), 'm**2/s')
        fvol_flow = float(cpline.token[7])

        results['case_parameters'] = CaseParameters(
            npipes, njunctions, nloops, maxiter,
            unitcode, tolerance, kinvisc, fvol_flow)

        _logger.debug('Successfully read case parameters')
        _logger.debug('  npipes = {0:d}'
                      .format(results['case_parameters'].npipes))
        _logger.debug('  njunctions = {0:d}'
                      .format(results['case_parameters'].njunctions))
        _logger.debug('  nloops = {0:d}'
                      .format(results['case_parameters'].nloops))
        _logger.debug('  maxiter = {0:d}'
                      .format(results['case_parameters'].maxiter))
        _logger.debug('  unitcode = {0:d}'
                      .format(results['case_parameters'].unitcode))
        _logger.debug('  tolerance = {0:0.4E}'
                      .format(results['case_parameters'].tolerance))
        _logger.debug('  kin_visc = {0:0.4E~}'
                      .format(results['case_parameters'].kin_visc))
        _logger.debug('  fvol_flow = {0:0.4f}'
                      .format(results['case_parameters'].fvol_flow))

    return results


def extract_pipe_definitions(deck, iptr, npipes, npipecards, unitcode):
    """Extract pipe definitions

    Args:
        deck [(InputLine)]: List of InputLine objects; user input lines
        iptr (int): Starting pointer for reading case parameters
        npipes (int): Number of pipes expected in model
        npipecards (int): Number of input lines needed to define a single
          characteristic of all pipes.
        unitcode (int): Code from input file selecting pipe dimension units

    Returns:
        (dict): Pipe dimensions and metadata
    """

    inunit = pipe_inputconv[unitcode]

    tmppipe = {
        'idiameter': [],
        'lpipe': [],
        'froughness': []
    }

    results = {
        'status': 'unknown',
        'msg': 'Pipe info results not initialized',
        'iread': 3 * npipecards,
        'pipe_info': []
    }

    for ict in range(npipecards):
        idiamctr = iptr + ict
        ilpipectr = idiamctr + npipecards
        ifroughctr = idiamctr + 2 * npipecards
        for tok in deck[idiamctr].token:
            tmppipe['idiameter'].append(Q_(float(tok), inunit['idiameter']))
        for tok in deck[ilpipectr].token:
            tmppipe['lpipe'].append(Q_(float(tok), inunit['lpipe']))
        for tok in deck[ifroughctr].token:
            tmppipe['froughness'].append(Q_(float(tok), inunit['froughness']))

    _logger.debug('Model contains {0:d} pipes:'.format(npipes))
    for ict in range(npipes):
        results['pipe_info'].append({'id': ict,
                                     'idiameter': tmppipe['idiameter'][ict],
                                     'lpipe': tmppipe['lpipe'][ict],
                                     'froughness': tmppipe['froughness'][ict]})
        currpipe = results['pipe_info'][ict]
        _logger.debug('  id={0:d}  D= {1:16.4E~} L={2:9.1f~} e={3:16.4E~}'
                      .format(currpipe['id'],
                              currpipe['idiameter'],
                              currpipe['lpipe'],
                              currpipe['froughness']))
    results['status'] = 'ok'
    results['message'] = 'Read diameter, length and roughness ' \
                         'for {:d} pipes'.format(npipes)

    return results


def extract_junctions(deck, iptr, njunctions):
    """Extract junction information from user input

    Args:
        deck [(InputLine)]: List of parsed lines of user input
        iptr (int): Pointer to first line of junction input
        njunctions (int): Number of junctions expected in model

    Returns:
        (dict): Junction info and metadata
    """
#    PipeRoute = namedtuple('PipeRoute', ['pipe_id', 'from', 'to', 'inflow'])
    results = {
        'status': 'unknown',
        'msg': 'Junction info results not initialized',
        'iread': 0,
        'junc_info': {
            'inflows': {},
            'pipe_map': {},
        }
    }

    ji = results['junc_info']

    jptr = iptr
    iread = 0
    ictr = 0
    while ictr < njunctions:
        ictr += 1
        junc_id = ictr - 1
        iread += 1
        inflow_units = int(deck[jptr].token[0])
        njpipes = int(deck[jptr].token[1])
        _logger.debug('Processing: {0:s}'.format(deck[jptr].as_log()))
        for itok in range(njpipes):
            pipe_id = int(deck[jptr].token[itok+2])
            pipe_in = pipe_id > 0
            pipe_read_id = abs(pipe_id)
            pipe_stored_id = pipe_read_id - 1
            if pipe_stored_id not in ji['pipe_map']:
                ji['pipe_map'][pipe_stored_id] = {'id': pipe_stored_id}

            if pipe_read_id == 0:
                results['status'] = 'error'
                results['msg'] = 'Junction {0:d} member {1:d} has pipe id ' \
                                 'out of range: (0):'.format(junc_id, itok+1)
                _logger.error(results['msg'])
            else:
                # Set to zero-index
                _logger.debug('Note: Pipe id adjusted from {0:d} to {1:d} due '
                              'to zero-indexing'.format(pipe_read_id,
                                                        pipe_stored_id))
                if pipe_in:
                    pipe_dir = 'inflow'
                    ji['pipe_map'][pipe_stored_id]['to'] = junc_id
                else:
                    pipe_dir = 'outflow'
                    ji['pipe_map'][pipe_stored_id]['from'] = junc_id
                _logger.debug('Junction {0:d} member {1:d} is pipe {2:d}, '
                              '{3:s}'.format(junc_id, itok+1, pipe_stored_id,
                                             pipe_dir))

        qinflow = Q_(0.0, 'm**3 / sec')
        if inflow_units < 1:
            _logger.debug('No inflow')
        elif inflow_units < 4:
            jptr += 1
            iread += 1
            _logger.debug('Processing inflow: {0:s}'
                          .format(deck[jptr].as_log()))
            qinflow = Q_(float(deck[jptr].token[0]),
                         flow_inputconv[inflow_units])
            _logger.debug('Inflow of {0:0.4E~}'.format(qinflow))
        else:
            _logger.error('Junction inflow unit specifier out of range')
            results['status'] = 'error'
            results['msg'] = 'Pipe {0:d} junction inflow unit specifier ' \
                             'out of range'.format(ictr)
        ji['inflows'][junc_id] = qinflow

        jptr += 1

    _logger.debug('jptr-iptr={0:d} vs iread={1:d}'.format(jptr-iptr, iread))
    _logger.debug('next line to read is #{0:d}: {1:s}'
                  .format(jptr, deck[jptr].as_log()))

    results['iread'] = iread
    if results['status'] == 'unknown':
        results['status'] = 'ok'
        results['msg'] = 'Read {0:d} junctions over {1:d} input lines.' \
                         .format(ictr, iread)

    _logger.debug(results['status'] + ': ' + results['msg'])

    _logger.debug('Pipe network topology:')
    for idx in sorted(ji['pipe_map'].keys()):
        _logger.debug('  Pipe {0:d} connects junction {1:d} to junction {2:d}'
                      .format(ji['pipe_map'][idx]['id'],
                              ji['pipe_map'][idx]['from'],
                              ji['pipe_map'][idx]['to']))

    _logger.debug('Junction inflows:')
    for idx in sorted(ji['inflows'].keys()):
        _logger.debug('  Inflow to junction {0:d} is {1:0.4E~}'
                      .format(idx, ji['inflows'][idx]))

    return results


def extract_loops(deck, iptr, nloops):
    """Extract pipe loop data for continuity
    Args:
        deck [(InputLine)]: List of parsed lines of user input
        iptr (int): Pointer to first line of loop input
        nloops (int): Number of loops expected in case

    Returns:
        (dict): Pipe loop info and metadata
    """

    results = {
        'status': 'unknown',
        'msg': 'Pipe loop results not initialized',
        'iread': 0,
        'loop_info': {}
    }
    LoopElement = namedtuple("LoopElement", ['pipe_id', 'flow_dir'])

    iloop = 0
    iread = 0
    jptr = iptr
    while iread < nloops:
        _logger.debug('Processing {0:s}'.format(deck[jptr].as_log()))
        nlooppipe = int(deck[jptr].token[0])
        assert nlooppipe == deck[jptr].ntok - 1
        results['loop_info'][iloop] = []
        for ipos in range(nlooppipe):
            pipe_read_id = int(deck[jptr].token[1+ipos])
            flow_conv = copysign(1.0, pipe_read_id)
            pipe_read_id = abs(pipe_read_id)
            pipe_stored_id = pipe_read_id - 1
            looppart = LoopElement(pipe_stored_id, flow_conv)
            results['loop_info'][iloop].append(looppart)
            if pipe_read_id == 0:
                results['status'] = 'error'
                results['msg'] = 'Element {0:d} of loop {1:d} has an ' \
                                 'invalid pipe id (0)'.format(ipos+1, iloop+1)

        iread += 1
        iloop += 1
        jptr += 1

    results['iread'] = iread
    if results['status'] == 'unknown':
        results['status'] = 'ok'
        results['msg'] = '{0:d} loops have been read'.format(iloop)

    _logger.debug(results['status'] + ': ' + results['msg'])

    assert nloops == iloop

    for idx in range(iloop):
        nparts = len(results['loop_info'][idx])
        _logger.debug('Loop {0:d} contains {1:d} elements:'
                      .format(idx, nparts))
        for jdx in range(nparts):
            if results['loop_info'][idx][jdx].flow_dir > 0.0:
                flow_dir = 'clockwise (forward)'
            else:
                flow_dir = 'counter-clockwise (reverse)'
            _logger.debug('  Element {0:d} is pipe {1:d}, flow is {2:s}'
                          .format(jdx+1,
                                  results['loop_info'][idx][jdx].pipe_id,
                                  flow_dir))

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

        case_dom = {}
        deck = []

        for ict, rawline in enumerate(fh):
            iline = InputLine(line=rawline, ipos=ict+1)
            _logger.debug(iline.as_log())

            if iline.typecode == 'D':
                deck.append(iline)

        iptr = 0
        while iptr < len(deck):

            # Step 1. Extract case parameters
            _logger.debug('1. Extracting case parameters')
            case_info = extract_case_parameters(deck, iptr)
            if case_info['status'] == 'error':
                _logger.error('Error reading case parameters: {0:s}'
                              .format(case_info['msg']))
                _logger.info('Skipping remainder of {0:s}'.format(fh.name))
                break

            case_dom['params'] = case_info['case_parameters']

            iptr += case_info['iread']

            npipes = case_info['case_parameters'].npipes
            npipecards = 1
            pipect = deck[iptr].ntok
            while pipect < npipes:
                npipecards += 1
                pipect += deck[iptr+npipecards-1].ntok

            _logger.debug('Found {0:d} of {1:d} pipes defined in {2:d} lines'
                          .format(pipect, npipes, npipecards))
            # Step 2. Read pipe data
            _logger.debug('2. Reading pipe data')
            unitcode = case_info['case_parameters'].unitcode
            pipe_info = extract_pipe_definitions(deck, iptr, npipes,
                                                 npipecards, unitcode)

            if pipe_info['status'] == 'error':
                _logger.error('Error reading pipe definitions: {0:s}'
                              .format(pipe_info['msg']))
                _logger.info('Skipping remainder of {0:s}'.format(fh.name))
                break

            case_dom['pipe'] = pipe_info['pipe_info']

            iptr += pipe_info['iread']

            # Step 3. Read junction data
            _logger.debug('3. Reading junction inflows and pipe network '
                          'topology')
            njunctions = case_info['case_parameters'].njunctions
            pipemap_info = extract_junctions(deck, iptr, njunctions)

            if pipemap_info['status'] == 'error':
                _logger.error('Error reading junction data: {0:s}'
                              .format(pipemap_info['msg']))
                _logger.info('Skipping remainder of {0:s}'.format(fh.name))
                break

            case_dom['junc'] = pipemap_info['junc_info']

            iptr += pipemap_info['iread']

            # Step 4. Read loop data
            _logger.debug('4. Reading loop continuity data')
            nloops = case_info['case_parameters'].nloops
            pipeloop_info = extract_loops(deck, iptr, nloops)

            if pipeloop_info['status'] == 'error':
                _logger.error('Error reading loop continuity data: {0:s}'
                              .format(pipeloop_info['msg']))
                _logger.info('Skipping remainder of {0:s}'.format(fh.name))
                break

            case_dom['loop'] = pipeloop_info['loop_info']

            iptr += pipeloop_info['iread']
            # Done reading input; assert (iptr + iread - 1) == len(deck) for a
            # single case

# #######################################################################
#        ! Calculate loss coefficient KP based on length and diameter
#        ! unit of measure

            # Note: Use SI coefficient since base units are SI (mks)
            kpcoeff = 2.12E-3
#            unitcode = case_dom['params'].unitcode
#            if unitcode in (0, 1):
#                # Coefficient for traditional (US) units
#                kpcoeff = 9.517E-4
#            elif unitcode in (2, 3):
#                # Coefficient for mks (SI) units
#                kpcoeff = 2.12E-3
#            else:
#                raise ValueError('Unknown unit code {0:d}, expected (0..3)'
#                                 .format(unitcode))

            for currpipe in case_dom['pipe']:
                ld = (currpipe['lpipe'] / currpipe['idiameter']
                      ).to_base_units().magnitude
                currpipe['flow_area'] = (
                    sc.pi / 4.0 * currpipe['idiameter']**2).to_base_units()
                currpipe['arl'] = (
                    currpipe['lpipe'] / (2.0 * ugrav * currpipe['idiameter']
                                         * currpipe['flow_area']**2)
                ).to_base_units()
#                _logger.debug('ARL is in {0:s}'.format(currpipe['arl'].units))
                currpipe['eroughness'] = (
                    currpipe['froughness']
                    / currpipe['idiameter']
                ).to_base_units().magnitude

                currpipe['expp'] = 0.0
                currpipe['kp'] = kpcoeff * ld ** 4.87

                _logger.debug('Pipe {0:d}: Kp = {1:0.4E}, Aflow = {2:0.4E~}'
                              .format(currpipe['id'],
                                      currpipe['kp'],
                                      currpipe['flow_area']))

# +        select case(NUNIT)
# +          case(0, 1)
# +            KP(1:NP) = 9.517E-4 * L(1:NP) / D(1:NP)**4.87
# +          case(2, 3)
# +            KP(1:NP) = 2.12E-3 * L(1:NP) / D(1:NP)**4.87
# +        end select
            nct = 0
            ssum = 100.0
            done = False
            converged = False
            qpredict = np.zeros(case_dom['params'].npipes)

# +        NCT = 0
# +        SSUM = 100.0
# +        CONVERGED = .false.
# +        DONE = .false.

# The goal of this project is to reimplement the JEPPSON_CH5 code in Python,
# not simply translate the original Fortran into Python. For this reason, pipe,
# junction, and loop indices are zero-based, default NumPy matrix storage is
# used. This complicates the comparison between the original Fortran and the
# Python implementations, but effectively illusrates the differences in
# implementing the solution in each language.

# The matrix a is constructed in the same manner as in the original Fortran
# application with a few modifications. While NumPy multidimensional array
# storage can be set to either C-style (row-major) or Fortran-style
# (column-major), the default NumPy style is used - C-style.

# Matrix columns represent flow through pipes. The first (njunctions-1) rows
# contain 1.0 in a column for outflow via that column's pipe or -1.0 for inflow
# from that column's pipe. The corresponding term in the first (njunctions - 1)
# rows of the column vector b contains the fixed flow into the junction from
# outside the system, positive for inflow and negative for outflow.

# The remaining nloops rows in the a matrix contain the flow resistance of each
# pipe, positive if the flow convention is clockwise/forward in the loop,
# negative if the flow convention is counter-clockwise/reverse in the loop.
# Flow direction/convention is heuristically assigned by the modeler; the
# actual flow direction will be determined from the problem solution. The
# corresponding entries in the b column vector are zero since the flow around a
# loop is conservative - no net increase or decrease.

            flow_units = ''
            npipes = case_dom['params'].npipes
            njunctions = case_dom['params'].njunctions
            while not done:
                # Step 5. Assemble matrix
                _logger.debug('5. Assemble matrix')
                pmap = case_dom['junc']['pipe_map']
                b = np.zeros((npipes))
                a = np.zeros((npipes, npipes))
                for idx in pmap:
                    pipe_id = pmap[idx]['id']
                    jfrom = pmap[idx]['from']
                    jto = pmap[idx]['to']
                    _logger.debug('Pipe {0:d} goes from {1:d} to {2:d}'
                                  .format(pipe_id, jfrom, jto))
                    if jfrom < njunctions - 1:
                        a[jfrom, pipe_id] = -1.0
                    if jto < njunctions - 1:
                        a[jto, pipe_id] = 1.0

                for idx in case_dom['junc']['inflows']:
                    # Use base units of first non-zero flow for result
                    # conversion.
                    if flow_units == '' and case_dom['junc']['inflows'][idx] \
                                            .magnitude != 0.0:
                        flow_units = case_dom['junc']['inflows'][idx] \
                                     .to_base_units().units
                    if idx < njunctions - 1:
                        b[idx] = case_dom['junc']['inflows'][idx] \
                                 .to_base_units().magnitude

                for iloop in case_dom['loop']:
                    row_id = njunctions - 1 + iloop
                    _logger.debug('Row id is {0:d} = njunctions + iloop = '
                                  '{1:d} + {2:d}'.format(row_id, njunctions,
                                                         iloop))
                    for pipe in case_dom['loop'][iloop]:
                        col_id = pipe.pipe_id
                        resistance = (pipe.flow_dir
                                      * case_dom['pipe'][pipe.pipe_id]['kp'])
                        _logger.debug('  Col id is {0:d}; resistance = '
                                      '{1:0.4E}'.format(col_id, resistance))
                        a[row_id, col_id] = resistance

                _logger.debug('Resultant flows are in units of {0:s}'
                              .format(flow_units))

                print(repr(a))
                print(repr(b))

                # Call matrix solver
                # Step 6. Solve matrix
                _logger.debug('6. Solve matrix')
                try:
                    x = np.linalg.solve(a, b)
                except np.linalg.LinAlgError as err:
                    msg = 'Cannot solve matrix: {0:s}'.format(str(err))
                    _logger.error(msg)
                    print('Error: ' + msg)
                    converged = False
                    break

                print(repr(x))

                _logger.debug('7. Adjust matrix')
                if nct > 0:
                    ssum = 0.0

# +L20:    do while (.not. DONE)
# +            II = 1
# +            do I = 1, NJ1
# +                A(I, 1:NPP) = 0.0
# +                NNJ = NN(I)        

# +                do J = 1, NNJ
# +                    A(I, abs(JN(I,J))) = PJDIR(I, J)
# +                end do

# +                if (IFLOW(I) /= 0) then
# +                    A(I,NPP) = QJ(II)
# +                    II = II + 1
# +                end if
# +            end do

# +            do I = NJ, NP
# +!?                A(I, 1:NPP) = 0.0
# +                A(I, 1:NP) = 0.0
# +                II = I - NJ1
# +!                NNJ = LP(NPLMAX+1, II)
# +                do J = 1, LP(NPLMAX+1, II)
# +                    IJ = LP(J ,II)
# +                    IIJ = abs(IJ)
# +                    if (IJ < 0) then
# +                        A(I, IIJ) = -KP(IIJ)
# +                    else
# +                        A(I, IIJ) = KP(IIJ)
# +                    end if
# +                end do
# +                A(I,NPP) = 0.0
# +            end do

# +            V(1) = 4.0
# +! System subroutine from UNIVAC MATH-PACK to solve linear system of
# +! equations
# +! Sperry alternate return point syntax not supported
# +!      CALL GJR(A, 51, 50, NP, NPP, $98, JC, V)
# +            call GJR(A, NPMAX+1, NPMAX, NP, NPP, 98, JC, V)

# +            if (JC(1) /= NP) then
# +                write(STDOUT, "('Matrix failure:', /,"                     &
# +                    // "'JC(1): ', I3, ' <> NP: ', I3)") JC(1), NP
# +                goto 98
# +            end if

# +            if (NCT > 0) then
# +                SSUM = 0.0
# +            end if

                for ipipe, currpipe in enumerate(case_dom['pipe']):
                    if nct > 0:
                        qm = 0.5 * (qpredict[ipipe] + x[ipipe])
                        ssum += abs(qpredict[ipipe] - x[ipipe])
                    else:
                        qm = x[ipipe]

                    qpredict[ipipe] = qm
                    dq = Q_(qm * case_dom['params'].fvol_flow, flow_units)
                    qmu = Q_(abs(qm), flow_units)
#                    vflowe = qmu / currpipe['flow_area']

                    qq_lo = qmu - dq
                    vflowv_lo = qq_lo / currpipe['flow_area']

                    qq_hi = qmu + dq
                    vflowv_hi = qq_hi / currpipe['flow_area']

                    if vflowv_lo.magnitude < 0.001:
                        vflowv_lo = Q_(0.002, vflowv_lo.units)
                        _logger.info('  Flow velocity low endpoint for pipe '
                                     '{0:d} increased to {1:0.4E~}'
                                     .format(ipipe, vflowv_lo))

                    re_lo = (vflowv_lo * currpipe['idiameter']
                             / case_dom['params'].kin_visc).to_base_units()

                    re_hi = (vflowv_hi * currpipe['idiameter']
                             / case_dom['params'].kin_visc).to_base_units()

                    _logger.debug('  Pipe {0:d} Re varies from {1:0.4E~} to '
                                  '{1:0.4E~}'.format(ipipe, re_lo, re_hi))

                    friction_lo = friction_factor(
                        Re=re_lo, eD=currpipe['eroughness'])

                    friction_hi = friction_factor(
                        Re=re_hi, eD=currpipe['eroughness'])

                    _logger.debug('  Pipe {0:d} f varies from {1:0.4E~} to '
                                  '{1:0.4E~}'.format(ipipe, friction_lo,
                                                     friction_hi))

                    # Calculate new Kp for each pipe based on flow regime and
                    # friction factor

                    # Note: Flow is laminar for Reynolds number less than 2050
                    if re_lo < 2050.0:
                        currpipe['expp'] = 1.0
                        tmp_kp = (
                            2.0 * ugrav * case_dom['params'].kin_visc
                            * currpipe['arl'] / currpipe['idiameter']
                        )
                        _logger.debug('  tmp_kp is in units of {0:s}'
                                .format(tmp_kp.units))
                        currpipe['kp'] = tmp_kp.to('1/ft**3/s').magnitude
                        _logger.debug('  Pipe {0:d} flow in laminar region'
                                      .format(ipipe))
                    else:
                        # Consider only transition regime, not transition-rough
                        be = ((log(friction_lo) - log(friction_hi))
                              / (log(qq_lo.to('ft**3/s').magnitude)
                                 - log(qq_hi.to('ft**3/s').magnitude)))
                        ae = friction_lo * qq_lo.to('ft**3/s').magnitude**be
                        ep = 1.0 - be
                        currpipe['expp'] = 2.0 - be
                        currpipe['kp'] = (
                            ae * currpipe['arl'].to('s**2/ft**5').magnitude
                            * qmu.to('ft**3/s').magnitude**ep).magnitude
                        _logger.debug('  arl is in units of {0:s}'
                                      .format(currpipe['arl'].units))
                        _logger.debug('  Pipe {0:d} flow in transition / '
                                      'turbulent region'.format(ipipe))
                    _logger.debug('  Pipe {0:d} Kp is updated to {1:0.4E}'
                                  .format(ipipe, currpipe['kp']))

#             do I = 1, NP
#                 BB = A(I,NPP)

#                 if (NCT > 0) then
#                     QM = 0.5 * (Q(I) + BB)
#                     SSUM = SSUM + abs(Q(I) - BB)
#                 else
#                     QM = BB
#                 end if

#                 Q(I) = QM
#                 DELQ = QM * DELQ1
#                 QM = abs(QM)
#                 VE = QM / AR(I)

#                 QQ(1) = QM - DELQ
#                 QQ(2) = QM + DELQ

#                 VV = QQ / AR(I)
#                 if (VV(1) < 0.001) then
#                     write(STDOUT, "('Note: V1 adjusted from ', ES12.4," &
#                                   // "' to ', ES12.4)") VV(1), 0.002
#                     VV(1) = 0.002
#                 end if
#                 RE = VV * D(I) / VIS

# !                if (RE(2) <= 2.1E3) then
#                 if (maxval(RE) <= 2.1E3) then
#                     ! Laminar flow
#                     F = f_laminar(RE)
#                     EXPP(I) = 1.0
# !                    KP(I) = 64.4 * VIS * ARL(I) / D(I)
# !orig                KP(I) = 2.0 * GRAV_ES * VIS * ARL(I) / D(I)
#                     KP(I) = G2 * VIS * ARL(I) / D(I)
#                 else
#                     F = f_rough(E(I), D(I))
#                     PAR = transition_rough_metric1(F(1), E(I), VE, VIS)
# !                    PAR = VE * sqrt(0.125 * F) * D(I) * E(I) / VIS
#                     if (PAR > 65.0) then
#                         ! Fully-rough-pipe flow
#                         KP(I) = F(1) * ARL(I) * QM**2
#                         EXPP(I) = 2.0
#                     else
#                         ! Newton's method root finder to calculate
#                         ! Darcy-Weisback friction factor in transition
#                         ! region

#                         F(1) = f_transition(F(1), QQ(1), D(I), E(I),    &
#                                             VIS, MAXMITER, MAXDIF)
#                         F(2) = f_transition(F(1), QQ(2), D(I), E(I),    &
#                                             VIS, MAXMITER, MAXDIF)

#                         BE = (log(F(1)) - log(F(2)))                    &
#                            / (log(QQ(1)) - log(QQ(2)))
#                         AE = F(1) * QQ(1)**BE
#                         EP = 1.0 - BE
#                         EXPP(I) = 2.0 - BE
#                         KP(I) = AE * ARL(I) * QM ** EP
#                     end if
#                 end if
#             end do

#             NCT = NCT + 1

# ! The next three cards can be removed
#             write(STDOUT,157) NCT, SSUM, Q(1:NP)
#             write(STDOUT,344) EXPP(1:NP)
#             write(STDOUT,344) KP(1:NP)

#             CONVERGED = (SSUM <= TOL)
#             DONE = (CONVERGED .or. (NCT >= MAXITER))
#         end do L20
                _logger.debug('Iteration {0:d}'.format(nct))
                _logger.debug('Deviation {0:0.4E}'.format(ssum))
                for ipipe, currpipe in enumerate(case_dom['pipe']):
                    _logger.debug('Pipe {0:d}: Kp = {1:0.4E}, expp = '
                        '{2:0.4E}, Q = {3:0.4E~}'
                        .format(ipipe,
                                currpipe['kp'], 
                                currpipe['expp'],
                                Q_(qpredict[ipipe], 'm**3/s')
                                .to('ft**3/s')))

                for iflow, xflow in enumerate(x):
                    qfinal = Q_(xflow, 'm**3/s')
                    _logger.debug('Pipe {0:d}: {1:12.4E~}    {2:12.4E~}'
                                  '    {3:12.4E~}'
                                  .format(iflow,
                                          qfinal.to('m**3/s'),
                                          qfinal.to('ft**3/s'), 
                                          qfinal.to('gallon/minute')))

                nct += 1

                _logger.debug('8. Check convergence')

                converged = ssum <= case_dom['params'].tolerance
                done = converged or (nct >= case_dom['params'].maxiter)
            # End iteration
# #######################################################################
            if not converged:
                _logger.info('Case not converged: ssum = {0:0.4E}'
                             .format(ssum))

            # Step 7. Display results
            _logger.debug('9. Display results')

        _logger.info('Done processing {0:s}'.format(fh.name))
    _logger.info("Ending jeppson_ch5")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
