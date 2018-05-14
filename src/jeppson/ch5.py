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
from __future__ import absolute_import, division, print_function

import argparse
# from collections import namedtuple, OrderedDict
import logging
from math import copysign, log
from os.path import abspath, splitext
import sys

# import iapws
import scipy.constants as sc
# from fluids.core import Reynolds
from fluids.friction import friction_factor
import numpy as np
import pygraphviz as pgv

from . import _logger, Q_
# from jeppson.pipe import Pipe
from jeppson.input import InputLine
from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"


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
        deck ([InputLine]): List of parsed input data lines
        iptr (int): Starting pointer for reading case parameters

    Returns:
        (dict): Case parameters

    Raises:
        ValueError: Too few case parameters detected
    """
    tags = ('npipes', 'njunctions', 'nloops', 'maxiter', 'unitcode',
            'tolerance', 'kin_visc', 'fvol_flow')
    mintok = len(tags)

    case_params = {
        '_iread': 1
    }

    cpline = deck[iptr]
    if cpline.ntok < mintok:
        msg = 'Too few entries found for case parameters ({0:d} expected)' \
              .format(mintok)
        raise ValueError(msg)
    elif cpline.ntok > mintok:
        msg = 'Too many entries found for case parameters ({0:d} found, ' \
              '{1:d} expected)'.format(cpline.ntok, mintok)
        _logger.warning(msg)

    case_params['npipes'] = int(cpline.token[0])
    case_params['njunctions'] = int(cpline.token[1])
    case_params['nloops'] = int(cpline.token[2])
    case_params['maxiter'] = int(cpline.token[3])
    case_params['unitcode'] = int(cpline.token[4])
    case_params['tolerance'] = float(cpline.token[5])
    if case_params['unitcode'] <= 1:
        case_params['kin_visc'] = Q_(float(cpline.token[6]), 'ft**2/s')
    else:
        case_params['kin_visc'] = Q_(float(cpline.token[6]), 'm**2/s')
    case_params['fvol_flow'] = float(cpline.token[7])

    _logger.debug('Successfully read case parameters')
    _logger.debug('  npipes = {0:d}'
                  .format(case_params['npipes']))
    _logger.debug('  njunctions = {0:d}'
                  .format(case_params['njunctions']))
    _logger.debug('  nloops = {0:d}'
                  .format(case_params['nloops']))
    _logger.debug('  maxiter = {0:d}'
                  .format(case_params['maxiter']))
    _logger.debug('  unitcode = {0:d}'
                  .format(case_params['unitcode']))
    _logger.debug('  tolerance = {0:0.4E}'
                  .format(case_params['tolerance']))
    _logger.debug('  kin_visc = {0:0.4E~}'
                  .format(case_params['kin_visc']))
    _logger.debug('  fvol_flow = {0:0.4f}'
                  .format(case_params['fvol_flow']))

    return case_params


def extract_pipe_definitions(deck, iptr, npipes, npipecards, unitcode):
    """Extract pipe definitions

    Args:
        deck ([InputLine]): List of InputLine objects; user input lines
        iptr (int): Starting pointer for reading case parameters
        npipes (int): Number of pipes expected in model
        npipecards (int): Number of input lines needed to define a single
          characteristic of all pipes.
        unitcode (int): Code from input file selecting pipe dimension units

    Returns:
        (list): List of dicts containing pipe dimensions
    """
    pipe_inputconv = [
        {'idiameter': 'in', 'lpipe': 'ft', 'froughness': 'in'},
        {'idiameter': 'ft', 'lpipe': 'ft', 'froughness': 'ft'},
        {'idiameter': 'm', 'lpipe': 'm', 'froughness': 'm'},
        {'idiameter': 'cm', 'lpipe': 'm', 'froughness': 'cm'}
    ]

    inunit = pipe_inputconv[unitcode]

    tmppipe = {
        'idiameter': [],
        'lpipe': [],
        'froughness': []
    }

    pipe_info = []

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

    ndiam = len(tmppipe['idiameter'])
    nlpipe = len(tmppipe['lpipe'])
    nrough = len(tmppipe['froughness'])

    if ndiam != nlpipe or ndiam != nrough:
        raise ValueError('Mismatched number of pipe attributes: {0:d} '
                         'diameters, {1:d} lengths, {2:d} roughnesses - all '
                         'should be equal'.format(ndiam, nlpipe, nrough))
    elif ndiam != npipes:
        raise ValueError('{0:d} of each pipe attribute found, {1:d} expected'
                         .format(ndiam, npipes))
    else:
        _logger.debug('Model contains {0:d} pipes ({1:d} expected):'
                      .format(ndiam, npipes))

    for ict in range(npipes):
        pipe_info.append({'id': ict,
                          'idiameter': tmppipe['idiameter'][ict],
                          'lpipe': tmppipe['lpipe'][ict],
                          'froughness': tmppipe['froughness'][ict]})
        currpipe = pipe_info[ict]
        _logger.debug('  id={0:d}  D= {1:16.4E~} L={2:9.1f~} e={3:16.4E~}'
                      .format(currpipe['id'],
                              currpipe['idiameter'],
                              currpipe['lpipe'],
                              currpipe['froughness']))

    return pipe_info


def extract_junctions(deck, iptr, njunctions, npipes):
    """Extract junction information from user input

    Args:
        deck ([InputLine]): List of parsed lines of user input
        iptr (int): Pointer to first line of junction input
        njunctions (int): Number of junctions expected in model

    Returns:
        (dict): Junction info and metadata

    Raises:
        ValueError: Pipe ID out of range (<1) or junction flow unit specifier
          out of range (not in [0..3])
    """

    flow_inputconv = [
        'skip',
        'gallon / minute',
        'ft**3 / sec',
        'm**3 / sec'
    ]

    results = {
        '_iread': 0,
        'junc_info': {
            'inflows': [],
            'pipe_map': []
        }
    }

    ji = results['junc_info']
    for idx in range(njunctions):
        ji['inflows'].append(Q_(0.0, 'm**3/s'))

    for idx in range(npipes):
        ji['pipe_map'].append({'id': idx, 'to': None, 'from': None})

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

            if pipe_read_id == 0:
                msg = 'Junction {0:d} member {1:d} has pipe id ' \
                      'out of range: (0):'.format(junc_id, itok+1)
                _logger.error(msg)
                raise ValueError(msg)
            else:
                # Set to zero-index
                _logger.debug('Note: Pipe id adjusted from {0:d} to {1:d} due '
                              'to zero-indexing'.format(pipe_read_id,
                                                        pipe_stored_id))
                pipespan = ji['pipe_map'][pipe_stored_id]
                if 'id' not in pipespan:
                    pipespan['id'] = pipe_stored_id

                if pipe_in:
                    pipe_dir = 'inflow'
                    pipespan['to'] = junc_id
                else:
                    pipe_dir = 'outflow'
                    pipespan['from'] = junc_id

                _logger.debug('Junction {0:d} member {1:d} is pipe {2:d}, '
                              '{3:s}'.format(junc_id, itok+1, pipe_stored_id,
                                             pipe_dir))

        if inflow_units < 1:
            _logger.debug('No inflow')
        elif inflow_units < 4:
            jptr += 1
            iread += 1
            _logger.debug('Processing junction {0:d} inflow: {1:s}'
                          .format(junc_id, deck[jptr].as_log()))
            ji['inflows'][junc_id] = Q_(float(deck[jptr].token[0]),
                                        flow_inputconv[inflow_units])
            _logger.debug('Inflow of {0:0.4E~}'
                          .format(ji['inflows'][junc_id]))
        else:
            msg = 'Pipe {0:d} junction inflow unit specifier out of range' \
                  .format(ictr)
            _logger.error(msg)
            raise ValueError(msg)

        jptr += 1

    _logger.debug('jptr-iptr={0:d} vs iread={1:d}'.format(jptr-iptr, iread))
    _logger.debug('next line to read is #{0:d}: {1:s}'
                  .format(jptr, deck[jptr].as_log()))

    results['_iread'] = iread
    _logger.debug('Read {0:d} junctions over {1:d} input lines.'
                  .format(ictr, iread))

    _logger.debug('Pipe network topology:')
    for idx, currpipe in enumerate(ji['pipe_map']):
        _logger.debug('  Pipe {0:d} connects junction {1:d} to junction {2:d}'
                      .format(currpipe['id'],
                              currpipe['from'],
                              currpipe['to']))

    _logger.debug('Junction inflows:')
    for idx, qinflow in enumerate(ji['inflows']):
        _logger.debug('  Inflow to junction {0:d} is {1:0.4E~}'
                      .format(idx, qinflow))

    return results


def extract_loops(deck, iptr, nloops):
    """Extract pipe loop data for continuity
    Args:
        deck ([InputLine]): List of parsed lines of user input
        iptr (int): Pointer to first line of loop input
        nloops (int): Number of loops expected in case

    Returns:
        (dict): Pipe loop info and metadata

    Raises:
        ValueError: Pipe ID out of range (<1) in loop definition
    """

    results = {
        '_iread': 0,
        'loop_info': []
    }

    iloop = 0
    iread = 0
    jptr = iptr
    while iread < nloops:
        _logger.debug('Processing {0:s}'.format(deck[jptr].as_log()))
        nlooppipe = int(deck[jptr].token[0])
        assert nlooppipe == deck[jptr].ntok - 1
        results['loop_info'].append([])
        for ipos in range(nlooppipe):
            pipe_read_id = int(deck[jptr].token[1+ipos])
            flow_conv = copysign(1.0, pipe_read_id)
            pipe_read_id = abs(pipe_read_id)
            pipe_stored_id = pipe_read_id - 1
            results['loop_info'][iloop].append(
                {'pipe_id': pipe_stored_id, 'flow_dir': flow_conv})
            if pipe_read_id == 0:
                msg = 'Element {0:d} of loop {1:d} has an invalid pipe id ' \
                      '(0)'.format(ipos+1, iloop+1)
                _logger.error(msg)
                raise ValueError(msg)

        iread += 1
        iloop += 1
        jptr += 1

    results['_iread'] = iread

    assert nloops == iloop

    for idx in range(iloop):
        nparts = len(results['loop_info'][idx])
        _logger.debug('Loop {0:d} contains {1:d} elements:'
                      .format(idx, nparts))
        for jdx in range(nparts):
            if results['loop_info'][idx][jdx]['flow_dir'] > 0.0:
                flow_dir = 'clockwise (forward)'
            else:
                flow_dir = 'counter-clockwise (reverse)'
            _logger.debug('  Element {0:d} is pipe {1:d}, flow is {2:s}'
                          .format(jdx+1,
                                  results['loop_info'][idx][jdx]['pipe_id'],
                                  flow_dir))

    return results


def extract_case(iptr, deck):
    """Extract complete pipe flow network analysis case from user input

    Args:
        iptr (int): Pointer to first unread line of user input in deck
        deck ([InputLine]): List of tokenized lines of user input

    Returns:
        (dict): Pipe flow network object model
        """

    case_dom = {}

    # Step 1. Extract case parameters
    _logger.debug('1. Extracting case parameters')

    case_dom['params'] = extract_case_parameters(deck, iptr)

    iptr += case_dom['params']['_iread']

    del(case_dom['params']['_iread'])

    npipes = case_dom['params']['npipes']
    npipecards = 1
    pipect = deck[iptr].ntok
    while pipect < npipes:
        npipecards += 1
        pipect += deck[iptr+npipecards-1].ntok

    _logger.debug('Found {0:d} of {1:d} pipes defined in {2:d} lines'
                  .format(pipect, npipes, npipecards))

    # Step 2. Read pipe data
    _logger.debug('2. Reading pipe data')
    unitcode = case_dom['params']['unitcode']
    case_dom['pipe'] = extract_pipe_definitions(deck, iptr, npipes,
                                                npipecards, unitcode)

    iptr += 3 * npipecards

    # Step 3. Read junction data
    _logger.debug('3. Reading junction inflows and pipe network '
                  'topology')
    njunctions = case_dom['params']['njunctions']
    pipemap_info = extract_junctions(deck, iptr, njunctions, npipes)

    # Migrate case_dom['junc']['pipe_map'] info into case_dom['pipe']
    for pmap in pipemap_info['junc_info']['pipe_map']:
        pipe_id = pmap['id']
        case_dom['pipe'][pipe_id]['to'] = pmap['to']
        case_dom['pipe'][pipe_id]['from'] = pmap['from']

    case_dom['inflows'] = pipemap_info['junc_info']['inflows']

    iptr += pipemap_info['_iread']

    # Step 4. Read loop data
    _logger.debug('4. Reading loop continuity data')
    nloops = case_dom['params']['nloops']
    pipeloop_info = extract_loops(deck, iptr, nloops)

    case_dom['loop'] = pipeloop_info['loop_info']

    iptr += pipeloop_info['_iread']

    # Done reading input; assert (iptr + iread - 1) == len(deck) for a
    # single case
    case_dom['params']['next_iptr'] = iptr

    return case_dom


def solve_network_flows(case_dom):
    """Find the volumetric flow and head loss for the piping network defined in
    the case_dom structure by using the linear method as described in Chapter 5
    of Jeppson.

    Args:
        case_dom (dict): Pipe flow network object model

    Raises:
        ValueError: Network solution matrix is singular or does not converge.
    """

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

    ugrav = Q_(sc.g, 'm/s**2')

    nct = 0
    ssum = 100.0
    flow_units = ''
    npipes = case_dom['params']['npipes']
    njunctions = case_dom['params']['njunctions']
    qpredict = np.zeros(npipes)

    # Set (njunctions-1) independent, conservative junction equations.
    # Note that the remaining junction equation can be derived from the
    # junction equations specified so far.

    _logger.debug('5a. Assemble constant portions of matrix and RHS vector')

    a = np.zeros((npipes, npipes))
    # This portion of the A matrix is constant
    for pipe in case_dom['pipe']:
        pipe_id = pipe['id']
        jfrom = pipe['from']
        jto = pipe['to']
        _logger.debug('Pipe {0:d} goes from {1:d} to {2:d}'
                      .format(pipe_id, jfrom, jto))
        if jfrom < njunctions - 1:
            a[jfrom, pipe_id] = -1.0
        if jto < njunctions - 1:
            a[jto, pipe_id] = 1.0

    # The B vector is constant
    b = np.zeros((npipes))
    for idx, inflow in enumerate(case_dom['inflows']):
        # Use base units of first non-zero flow for result
        # conversion.
        if flow_units == '' and inflow.magnitude != 0.0:
            flow_units = inflow.to_base_units().units
        if idx < njunctions - 1:
            b[idx] = inflow.to_base_units().magnitude

    done = False
    converged = False

    while not done:
        # Step 5. Assemble matrix
        _logger.debug('5b. Assemble matrix rows of loop equations')
        for iloop, looppipe in enumerate(case_dom['loop']):
            row_id = njunctions - 1 + iloop
            _logger.debug('Row id is {0:d} = njunctions + iloop = '
                          '{1:d} + {2:d}'.format(row_id, njunctions,
                                                 iloop))
            for pipe in looppipe:
                pid = pipe['pipe_id']
                col_id = pid
                resistance = (pipe['flow_dir']
                              * case_dom['pipe'][pid]['kp'])
                _logger.debug('  Col id is {0:d}; resistance = {1:0.4E}'
                              .format(col_id, resistance))
                a[row_id, col_id] = resistance

        _logger.debug('Resultant flows are in units of {0:s}'
                      .format(flow_units))

#        print(repr(a))
#        print(repr(b))

        # Call matrix solver
        # Step 6. Solve matrix
        _logger.debug('6. Solve matrix')
        try:
            x = np.linalg.solve(a, b)
        except np.linalg.LinAlgError as err:
            msg = 'Cannot solve matrix: {0:s}'.format(str(err))
            _logger.error(msg)
            print('Error: ' + msg)
            print('A Matrix:\n{:s}\n'.format(repr(a)))
            print('B Vector:\n{:s}\n'.format(repr(b)))
            converged = False
            # force-exit iteration loop
            break

#        print(repr(x))

        _logger.debug('7. Adjust matrix')
        if nct > 0:
            ssum = 0.0

        for ipipe, currpipe in enumerate(case_dom['pipe']):
            if nct > 0:
                qm = 0.5 * (qpredict[ipipe] + x[ipipe])
                ssum += abs(qpredict[ipipe] - x[ipipe])
            else:
                qm = x[ipipe]

            qpredict[ipipe] = qm
            dq = Q_(qm * case_dom['params']['fvol_flow'], flow_units)
            qmu = Q_(abs(qm), flow_units)
#            vflowe = qmu / currpipe['flow_area']

            qq_lo = qmu - dq
            vflowv_lo = qq_lo / currpipe['flow_area']

            qq_hi = qmu + dq
            vflowv_hi = qq_hi / currpipe['flow_area']

            if vflowv_lo.magnitude < 0.001:
                vflowv_lo = Q_(0.002, vflowv_lo.units)
                _logger.info('  Flow velocity low endpoint for pipe {0:d} '
                             'increased to {1:0.4E~}'.format(ipipe, vflowv_lo))

            re_lo = (vflowv_lo * currpipe['idiameter']
                     / case_dom['params']['kin_visc']).to_base_units()

            re_hi = (vflowv_hi * currpipe['idiameter']
                     / case_dom['params']['kin_visc']).to_base_units()

            _logger.debug('  Pipe {0:d} Re varies from {1:0.4E~} to '
                          '{1:0.4E~}'.format(ipipe, re_lo, re_hi))

            friction_lo = friction_factor(Re=re_lo, eD=currpipe['eroughness'])

            friction_hi = friction_factor(Re=re_hi, eD=currpipe['eroughness'])

            _logger.debug('  Pipe {0:d} f varies from {1:0.4E~} to {1:0.4E~}'
                          .format(ipipe, friction_lo, friction_hi))

            # Calculate new Kp for each pipe based on flow regime and
            # friction factor

            # Note: Flow is laminar for Reynolds number less than 2050
            if re_lo < 2050.0:
                currpipe['expp'] = 1.0
                tmp_kp = (2.0 * ugrav * case_dom['params']['kin_visc']
                          * currpipe['arl'] / currpipe['idiameter'])
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

        _logger.debug('8. Display interim results')
        print('Iteration {0:d}'.format(nct))
        print('Deviation {0:0.4E} (Tolerance {1:0.4E}'
              .format(ssum, case_dom['params']['tolerance']))

        print()
        print('Pipe   Kp            expp          Qcurrent                  '
              'Qpredict')
        for ipipe, currpipe in enumerate(case_dom['pipe']):
            print('{0:3d}    {1:0.4E}    {2:0.4E}    {3:0.4E~}    {4:0.4E~}'
                  .format(ipipe,
                          currpipe['kp'],
                          currpipe['expp'],
                          Q_(x[ipipe], 'm**3/s').to('ft**3/s'),
                          Q_(qpredict[ipipe], 'm**3/s').to('ft**3/s')))
        print()

        nct += 1

        _logger.debug('9. Check convergence')

        converged = ssum <= case_dom['params']['tolerance']
        done = converged or (nct >= case_dom['params']['maxiter'])

    # End iteration
# ########################################################################
    if not converged:
        msg = 'Case not converged: ssum = {0:0.4E} > tolerance {1:0.4E}' \
              .format(ssum, case_dom['params']['tolerance'])
        # Advance to next case
        raise ValueError(msg)

    # Add final results to case_dom
    flow_disp_units = 'm**3/s'
    for qext in case_dom['inflows']:
        if qext != 0.0:
            flow_disp_units = qext.units
            break

    for ipipe, currpipe in enumerate(case_dom['pipe']):
        currpipe['vol_flow'] = Q_(x[ipipe], 'm**3/s').to(flow_disp_units)
        currpipe['head_loss'] = (Q_(currpipe['kp']
                                    * currpipe['vol_flow']
                                    .to('ft**3/s').magnitude, 'ft'))

    return


def set_pipe_derived_properties(pipelist):
    """Set derived properties of pipe - flow area, arl, eroughness, and initial
    kp and expp

    Args:
        pipelist ([dicr]): List of dicts containing pipe dimensions and
          attributes"""
    ugrav = Q_(sc.g, 'm/s**2')

    for currpipe in pipelist:
        currpipe['LD'] = (
            currpipe['lpipe'] / currpipe['idiameter']
        ).to_base_units().magnitude

        currpipe['flow_area'] = (
            sc.pi / 4.0 * currpipe['idiameter']**2
        ).to_base_units()

        currpipe['arl'] = (
            currpipe['lpipe'] / (2.0 * ugrav * currpipe['idiameter']
                                 * currpipe['flow_area']**2)
        ).to_base_units()

        currpipe['eroughness'] = (
            currpipe['froughness']
            / currpipe['idiameter']
        ).to_base_units().magnitude

        # For the initial pass, the value of kpcoeff should not matter.
        # kpcoeff is a constant multiplier on all non-zero terms in each loop
        # equation. Loop equations are homogeneous (B[iloop] = 0) so the
        # results are only affected by the (L/D) exponential term for each
        # pipe; only the relative magnitudes of the coefficient in the loop
        # equations determine the calculated flows.

        # Later iterations of the solver loop use laminar or turbulent
        # correlations to calculate kp; kpcoeff is not used after this point.

# Method 1:
#        # Note: Use SI coefficient since base units are SI (mks)
#        kpcoeff = 2.12E-3.

# Method 2:
        # Note: Use US coefficient since dP/dQ calculations are done using
        # Q in ft**3/s
        kpcoeff = 9.517E-4

# Method 3:
#        unitcode = case_dom['params'].unitcode
#        if unitcode in (0, 1):
#            # Coefficient for traditional (US) units
#            kpcoeff = 9.517E-4
#        elif unitcode in (2, 3):
#            # Coefficient for mks (SI) units
#            kpcoeff = 2.12E-3
#        else:
#            raise ValueError('Unknown unit code {0:d}, expected (0..3)'
#                             .format(unitcode))

        currpipe['expp'] = 0.0
        currpipe['kp'] = kpcoeff * currpipe['LD'] ** 4.87

        _logger.debug('Pipe {:d}: '
                      'Aflow = {:0.4E~}, '
                      'arl = {:0.4E~}, '
                      'eD = {:0.4E}, '
                      'LD = {:0.4E}, '
                      'Kp = {:0.4E}, '
                      'expp = {:0.4E}'
                      .format(currpipe['id'],
                              currpipe['flow_area'],
                              currpipe['arl'],
                              currpipe['eroughness'],
                              currpipe['LD'],
                              currpipe['kp'],
                              currpipe['expp']))

    return


def pipe_dimension_table(pipelist):
    """Print pipe dimensions in a tabular text format

    Args:
        pipelist ([dict]): List of dicts containing pipe data

    Returns:
        (str): Text table of pipe dimensions"""

    result = 'Pipe    Diameter         Length           Rel. Roughness\n'
    for currpipe in pipelist:
        result += '{0:3d}     {1:12.4E~}    {2:12.4E~}  {3:12.4E}\n' \
                  .format(currpipe['id'], currpipe['idiameter'],
                          currpipe['lpipe'], currpipe['eroughness'])
    return result


def flow_and_head_loss_report(case_dom):
    """Return a string containing the final calculated volumetric flow and head
    loss in each pipe

    Args:
        case_dom (dict): Pipe network object model

    Returns:
        (str): Table containing the final calculated volumetric flow and head
          loss in each pipe
    """
    outstr = 'Pipe  Flow                     Flow                      Flow' \
             '                    Head Loss       Head Loss\n'
    for ipipe, currpipe in enumerate(case_dom['pipe']):
        outstr += '{0:-3d}   {1:12.4E~}    {2:12.4E~}    {3:12.4E~}    ' \
                  '{4:12.4E~}    {5:12.4E~}\n' \
                  .format(ipipe,
                          currpipe['vol_flow'].to('m**3/s'),
                          currpipe['vol_flow'].to('ft**3/s'),
                          currpipe['vol_flow'].to('gallon/minute'),
                          currpipe['head_loss'].to('m'),
                          currpipe['head_loss'].to('ft'))

    return outstr


def create_topology_dotfile(case_dom, filepath='tmp.gv'):
    """Create directed graph in GraphViz ``dot`` format

    Args:
        case_dom (dict): Pipe network object model
        filepath (str): Absolute path of dotfile"""

    jfmt = 'J{:d}'
    jxfmt = 'JX{:d}'
#    pfmt = 'P{:d}'
    pqfmt = 'P{:d}: {:0.1f~}'

    G = pgv.AGraph(directed=True, splines=False, ratio='fill', overlap=False)

    for idx, inflow in enumerate(case_dom['inflows']):
        jtag = jfmt.format(idx)
        jxtag = jxfmt.format(idx)
        G.add_node(jtag, shape='circle', label=jtag)
        if inflow.magnitude > 0.0:
            G.add_node(jxtag, shape='none', label='')
            G.add_edge(jxfmt.format(idx), jtag, color='blue',
                       label='{:0.1f~}'.format(inflow))
        elif inflow.magnitude < 0.0:
            G.add_node(jxtag, shape='none', label='')
            G.add_edge(jtag, jxfmt.format(idx), color='red',
                       label='{:0.1f~}'.format(abs(inflow)))

    for idx, link in enumerate(case_dom['pipe']):
        pqtag = pqfmt.format(idx, link['vol_flow'])
        jfrom = jfmt.format(link['from'])
        jto = jfmt.format(link['to'])
        G.add_edge(jfrom, jto, label=pqtag)

    dotfn = abspath(filepath)
    basepath, ext = splitext(dotfn)
    pngpath = abspath(basepath + '.png')
    _logger.debug('Writing png to {:s}'.format(pngpath))
    G.write(dotfn)
    G.draw(path=pngpath, prog='dot')

    return


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

        icase = 0
        deck = []

        for ict, rawline in enumerate(fh):
            iline = InputLine(line=rawline, ipos=ict+1)
            _logger.debug(iline.as_log())

            if iline.typecode == 'D':
                deck.append(iline)

        iptr = 0
        while iptr < len(deck):
            icase += 1

            _logger.info('Reading case {0:d} from {1:s}'
                         .format(icase, fh.name))
            try:
                case_dom = extract_case(iptr, deck)
                iptr = case_dom['params']['next_iptr']
                del case_dom['params']['next_iptr']
            except ValueError as err:
                _logger.error('Failed to read case {0:d} from {1:s}: {2:s}'
                              .format(icase, fh.name, str(err)))
                _logger.info('Advancing to next file.')
                break

            set_pipe_derived_properties(case_dom['pipe'])

            print(pipe_dimension_table(case_dom['pipe']))
            print()

            try:
                solve_network_flows(case_dom)
            except ValueError as err:
                _logger.error('Failed to solve case {0:d} from {1:s}: {2:s}'
                              .format(icase, fh.name, str(err)))
                _logger.info('Advancing to next case.')
                continue

            # Step 10. Display results
            _logger.debug('10. Display final results')

            print(flow_and_head_loss_report(case_dom))

            dotfn = abspath((splitext(fh.name))[0] + '_{:d}.gv'.format(icase))
            create_topology_dotfile(case_dom, dotfn)

            _logger.info('Done processing case {:d}'.format(icase))
        _logger.info('Done processing {0:s}'.format(fh.name))
    _logger.info("Ending jeppson_ch5")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
