#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This reimplements the second pipe flow analysis program given in chapter 6 of
*Steady Flow Analysis of Pipe Networks: An Instructional Manual* (1974).
*Reports.* Paper 300.  http://digitalcommons.usu.edu/water_rep/300 and
*Analysis of Flow in Pipe Networks* (1976). Ann Arbor Science Publishers, Inc.
http://www.worldcat.org/title/analysis-of-flow-in-pipe-networks/oclc/927534147
by Roland W. Jeppson.

Then run `python setup.py install` which will install the command
`jeppson_ch6b` inside your current environment.
"""
from __future__ import division, print_function, absolute_import

import argparse
# from collections import namedtuple, OrderedDict
import logging
import pprint
from math import copysign
from os.path import abspath, splitext
import sys

# import iapws
import scipy.constants as sc
# from fluids.core import Reynolds
import numpy as np
import pygraphviz as pgv

from . import _logger, Q_
from jeppson.input import InputLine
from jeppson.constants import ahws_us, echw, edhw

from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"

_pp = pprint.PrettyPrinter(indent=4)


def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Corrective flow Newton-Raphson flow network solver")
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
    tags = ('npipes', 'nloops', 'maxiter', 'npumps', 'npseudoloops')
    mintok = len(tags)

    case_params = {
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

    for itok, tag in enumerate(tags):
        case_params[tag] = int(cpline.token[itok])

    _logger.debug('Successfully read case parameters')
    for tag in tags:
        _logger.debug('  {0:s} = {1:d}'.format(tag, case_params[tag]))

    return case_params


def extract_pipe_definitions(deck, iptr, npipes):
    """Extract pipe definitions

    Args:
        deck ([InputLine]): List of InputLine objects; user input lines
        iptr (int): Starting pointer for reading case parameters
        npipes (int): Number of pipes expected in model

    Returns:
        (list): List of dicts containing pipe dimensions
    """

    tags = ('id', 'idiameter', 'lpipe', 'froughness', 'init_vol_flow')
    mintok = len(tags)

    pipe_info = []
    for ipipe in range(npipes):
        pipe_info.append({'id': ipipe})

    for ipipe in range(npipes):
        currline = deck[iptr + ipipe]
        if currline.ntok < mintok:
            msg = 'Too few tokens found in pipe definition line ({0:d} ' \
                'found, {1:d} expected)'.format(currline.ntok, mintok)
            _logger.error(msg)
            _logger.error(currline.as_log())
            raise ValueError(msg)
        if currline.ntok > mintok:
            msg = 'Too many tokens found in pipe definition line ({0:d} ' \
                'found, {1:d} expected)'.format(currline.ntok, mintok)
            _logger.warning(msg)

        pipe_id = int(currline.token[0]) - 1

        pipe_info[pipe_id]['idiameter'] = Q_(float(currline.token[1]), 'in')
        pipe_info[pipe_id]['lpipe'] = Q_(float(currline.token[2]), 'foot')
        pipe_info[pipe_id]['chw'] = float(currline.token[3])
        pipe_info[pipe_id]['init_vol_flow'] = Q_(float(currline.token[4]),
                                                 'ft**3/sec')

    for currpipe in pipe_info:
        _logger.debug('  Pipe {:d} D={:4.1f~} L={:8.1E~} CHW={:9.1f} '
                      'Qinit={:8.2f~}'
                      .format(currpipe['id'],
                              currpipe['idiameter'],
                              currpipe['lpipe'],
                              currpipe['chw'],
                              currpipe['init_vol_flow']))

    return pipe_info


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

    loop_info = []
    for iloop in range(nloops):
        currline = deck[iptr+iloop]
        _logger.debug('Processing {0:s}'.format(currline.as_log()))
        nlooppipe = int(currline.token[0])
        assert nlooppipe == currline.ntok - 1

        loop_info.append([])

        for ipos in range(nlooppipe):
            pipe_read_id = int(currline.token[1+ipos])
            flow_conv = copysign(1.0, pipe_read_id)
            pipe_read_id = abs(pipe_read_id)
            pipe_stored_id = pipe_read_id - 1
            loop_info[iloop].append(
                {'pipe_id': pipe_stored_id, 'flow_dir': flow_conv})
            if pipe_read_id == 0:
                msg = 'Element {0:d} of loop {1:d} has an invalid pipe id ' \
                      '(0)'.format(ipos+1, iloop+1)
                _logger.error(msg)
                raise ValueError(msg)

    for idx in range(nloops):
        nparts = len(loop_info[idx])
        _logger.debug('Loop {0:d} contains {1:d} elements:'
                      .format(idx, nparts))
        for jdx in range(nparts):
            if loop_info[idx][jdx]['flow_dir'] > 0.0:
                flow_dir = 'clockwise (forward)'
            else:
                flow_dir = 'counter-clockwise (reverse)'
            _logger.debug('  Element {0:d} is pipe {1:d}, flow is {2:s}'
                          .format(jdx+1,
                                  loop_info[idx][jdx]['pipe_id'],
                                  flow_dir))

    return loop_info


def extract_pumps(deck, iptr, npumps):
    """Extract pump curve, pipe ID, and flow direction
    Args:
        deck ([InputLine]): List of parsed lines of user input
        iptr (int): Pointer to first line of pump input
        npumps (int): Number of pumps expected in case

    Returns:
        (list): Pump info

    Raises:
        ValueError: Pipe ID out of range (<1) in pump definition
    """

    mintok = 4
    pump_info = []
    for ipump in range(npumps):
        currline = deck[iptr+ipump]
        _logger.debug('Processing {0:s}'.format(currline.as_log()))

        if currline.ntok < mintok:
            msg = 'Too few entries found for pump ({0:d} found, {1:d} ' \
                  'expected)'.format(currline.ntok, mintok)
            raise ValueError(msg)
        elif currline.ntok > mintok:
            msg = 'Too many entries found for pump ({0:d} found, {1:d} ' \
                  'expected)'.format(currline.ntok, mintok)
            _logger.warning(msg)

        ipipe = int(currline.token[0])
        flow_dir = copysign(1.0, ipipe)
        ipipe = abs(ipipe) - 1

        pump_info.append({
            'pipe_id': ipipe,
            'a': Q_(float(currline.token[1]), 's**2/ft**5'),
            'b': Q_(float(currline.token[2]), 's/ft**2'),
            'h0': Q_(float(currline.token[3]), 'ft'),
            'flow_dir': flow_dir
        })

    for ipump, currpump in enumerate(pump_info):
        if currpump['flow_dir'] < 0.0:
            dir_tag = 'reverse'
        else:
            dir_tag = 'forward'

        _logger.debug('Pump {0:d} curve is {1:0.4E~} Q**2 + ({2:0.4E~}) Q + '
                      '{3:0.4E~}, {4:s} flow along pipe {5:d}'
                      .format(ipump,
                              currpump['a'],
                              currpump['b'],
                              currpump['h0'],
                              dir_tag,
                              currpump['pipe_id']))

    return pump_info


def extract_pseudoloops(deck, iptr, npseudoloops):
    """Extract pipe pseudoloop data for continuity
    Args:
        deck ([InputLine]): List of parsed lines of user input
        iptr (int): Pointer to first line of pseudoloop input
        npseudoloops (int): Number of pseudoloops expected in case

    Returns:
        (list): Pseudoloop info (loop ID and head difference)

    Raises:
        ValueError: Too few user input entries
    """

    mintok = 2
    pseudoloop_info = []
    for iploop in range(npseudoloops):
        currline = deck[iptr+iploop]
        _logger.debug('Processing {0:s}'.format(currline.as_log()))

        if currline.ntok < mintok:
            msg = 'Too few entries found for pseudoloop ({0:d} found, {1:d} ' \
                  'expected)'.format(currline.ntok, mintok)
            raise ValueError(msg)
        elif currline.ntok > mintok:
            msg = 'Too many entries found for pseudoloop ({0:d} found, ' \
                  '{1:d} expected)'.format(currline.ntok, mintok)
            _logger.warning(msg)

        iloop = int(currline.token[0]) - 1
        head_diff = Q_(float(currline.token[1]), 'ft')
        pseudoloop_info.append({'loop_id': iloop, 'head_diff': head_diff})

    _logger.debug('Head differences for each pseudoloop:')
    for currploop in pseudoloop_info:
        _logger.debug('Loop {0:d}: {1:0.4E~}'
                      .format(currploop['loop_id'], currploop['head_diff']))

    return pseudoloop_info


def extract_case(iptr, deck):
    """Extract complete pipe flow network analysis case from user input

    Args:
        iptr (int): Pointer to first unread line of user input in deck
        deck ([InputLine]): List of tokenized lines of user input

    Returns:
        (dict): Pipe flow network data model
        """

    case_dom = {}

    # Step 1. Extract case parameters
    _logger.debug('1. Extracting case parameters')

    case_dom['params'] = extract_case_parameters(deck, iptr)

    iptr += 1

    # Step 2. Read pipe data
    _logger.debug('2. Reading pipe data')
    npipes = case_dom['params']['npipes']
    case_dom['pipe'] = extract_pipe_definitions(deck, iptr, npipes)

    iptr += npipes

    # Step 3. Read loop data
    _logger.debug('3. Reading loop continuity data')
    nloops = case_dom['params']['nloops']

    case_dom['loop'] = extract_loops(deck, iptr, nloops)

    iptr += nloops

    # Step 4. Read pump data
    npumps = case_dom['params']['npumps']
    if npumps > 0:
        pump_info = extract_pumps(deck, iptr, npumps)

        case_dom['pump'] = pump_info
        iptr += npumps

    # Step 5. Read pseudoloop data
    npseudoloops = case_dom['params']['npseudoloops']
    if npseudoloops > 0:
        pseudoloop_info = extract_pseudoloops(deck, iptr, npseudoloops)

        case_dom['pseudoloop'] = pseudoloop_info
        iptr += npseudoloops

    # Done reading input
    case_dom['params']['next_iptr'] = iptr

    return case_dom


def update_vol_flows(case_dom, dq):
    """Add corrective flows to initial flow to get current flow in each pipe

    Args:
        case_dom (dict): Pipe flow network data model
        dq (numpy.ndarray): List of corrective flows in cubic feet per second
    """

    for currpipe in case_dom['pipe']:
        currpipe['vol_flow'] = currpipe['init_vol_flow']
        for iloop, flow_dir in currpipe['loopdir'].items():
            currpipe['vol_flow'] += Q_(flow_dir * dq[iloop], 'ft**3/s')


def set_loop_head_loss(case_dom, dq):
    """Find total head loss around a pipe loop

    Args:
        case_dom (dict): Pipe flow network data model

    Returns:
        (numpy.ndarray): List of loop head losses in units of feet
    """

    b = np.zeros((case_dom['params']['nloops']))
    for iloop, loop in enumerate(case_dom['loop']):
        for ilpipe, looppipe in enumerate(loop):
            currpipe = case_dom['pipe'][looppipe['pipe_id']]
            qabs = abs(currpipe['vol_flow'].to('ft**3/s').magnitude)
            qsgn = copysign(1.0, currpipe['vol_flow'].magnitude)

            b[iloop] += (looppipe['flow_dir'] * currpipe['kp']
                         * qsgn * qabs**echw)

    if case_dom['params']['npseudoloops'] > 0:
        for psloop in case_dom['pseudoloop']:
            iloop = psloop['loop_id']
            b[iloop] -= psloop['head_diff'].to('ft').magnitude
            pipelist = []
            for looppipe in case_dom['loop'][iloop]:
                pipelist.append(looppipe['pipe_id'])

            if case_dom['params']['npumps'] > 0:
                for pump in case_dom['pump']:
                    pid = pump['pipe_id']
                    if pid not in pipelist:
                        continue
#                    pipedir = 0
#                    for looppipe in case_dom['loop'][iloop]:
#                        if looppipe['pipe_id'] == pid:
#                            pipedir = looppipe['flow_dir']
                    currpipe = case_dom['pipe'][pid]
                    qi = currpipe['init_vol_flow'].to('ft**3/s').magnitude
                    qcfs = abs(qi * pump['flow_dir'] + dq[iloop])
                    q = Q_(qcfs, 'ft**3/s')
                    hp = ((pump['a'] * q) + pump['b']) * q + pump['h0']
#                    dhp = 2.0 * pump['a'] * q + pump['b']
                    b[iloop] -= pump['flow_dir'] * hp.to('ft').magnitude

    return b


def set_flow_jacobian(case_dom, dq):
    """Find total head loss around a pipe loop

    Args:
        case_dom (dict): Pipe flow network data model

    Returns:
        (numpy.ndarray): Matrix of partial derivativese
    """

    nloops = case_dom['params']['nloops']
    a = np.zeros((nloops, nloops))
    for iloop, loop in enumerate(case_dom['loop']):
        for looppipe in loop:
            pid = looppipe['pipe_id']
            currpipe = case_dom['pipe'][pid]
            q = currpipe['init_vol_flow'].to('ft**3/s').magnitude
            for jloop, jflow_dir in currpipe['loopdir'].items():
                q += jflow_dir * dq[jloop]

            for jloop in currpipe['loopdir']:
                a[iloop, jloop] += (currpipe['loopdir'][jloop]
                                    * echw * looppipe['flow_dir']
                                    * currpipe['kp'] * abs(q)**(echw - 1.0))

    if case_dom['params']['npseudoloops'] > 0:
        for psloop in case_dom['pseudoloop']:
            iloop = psloop['loop_id']
            pipelist = []
            for looppipe in case_dom['loop'][iloop]:
                pipelist.append(looppipe['pipe_id'])

            if case_dom['params']['npumps'] > 0:
                for pump in case_dom['pump']:
                    pid = pump['pipe_id']
                    if pid not in pipelist:
                        continue
#                    pipedir = 0
#                    for looppipe in case_dom['loop'][iloop]:
#                        if looppipe['pipe_id'] == pid:
#                            pipedir = looppipe['flow_dir']
                    currpipe = case_dom['pipe'][pid]
#                    pipedir = currpipe['loopdir'][iloop]
                    qi = currpipe['init_vol_flow'].to('ft**3/s').magnitude
                    qcfs = abs(qi * pump['flow_dir'] + dq[iloop])
                    q = Q_(qcfs, 'ft**3/s')
#                    hp = ((pump['a'] * q) + pump['b']) * q + pump['h0']
                    dhp = (2.0 * pump['a'] * q + pump['b']) \
                        .to('sec/ft**2').magnitude
#                    b[iloop] -= pump['flow_dir'] * hp.to('ft').magnitude
#                    LLP is a pipe ID, not a loop ID - how does this work?
#                    a[iloop,LLP(IK)] += pump['flow_dir'] * dhp
                    a[iloop, pid] += pump['flow_dir'] * dhp

    return a


def solve_network_flows(case_dom):
    """Find the volumetric flow and head loss for the piping network defined in
    the case_dom structure by using the first method described in Chapter 6 of
    Jeppson.

    Args:
        case_dom (dict): Pipe flow network data model

    Raises:
        ValueError: Network solution matrix is singular or does not converge.
    """

#    ugrav = Q_(sc.g, 'm/s**2')
    nct = 0
    ssum = 0.0
    case_dom['params']['tolerance'] = 1.0E-3
#    npipes = case_dom['params']['npipes']
    nloops = case_dom['params']['nloops']

    # Set (npipes) independent equations.

    _logger.debug('5a. Prepare to enter iteration loop')

    done = False
    converged = False

    a = np.zeros((nloops, nloops))
    b = np.zeros((nloops))
    dq = np.zeros((nloops))

# Set initialize current volumetric flow in each pipe
    for currpipe in case_dom['pipe']:
        currpipe['vol_flow'] = currpipe['init_vol_flow']

    while not done:
        # Step 5. Assemble matrix
        _logger.debug('5b. Assemble a matrix and b vector')

        a = set_flow_jacobian(case_dom, dq)
        b = set_loop_head_loss(case_dom, dq)

        print('A:')
        _pp.pprint(a)
        print('B:')
        _pp.pprint(b)

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

        print('X Vector:\n' + repr(x))
        dq -= x
        ssum = sum(abs(x))

        _logger.debug('7. Adjust matrix')

#        for ipipe, currpipe in enumerate(case_dom['pipe']):
#            if nct > 0:
#                qm = 0.5 * (qpredict[ipipe] + x[ipipe])
#                ssum += abs(qpredict[ipipe] - x[ipipe])
#            else:
#                qm = x[ipipe]
#
#            _logger.debug('  Pipe {0:d} Kp is updated to {1:0.4E}'
#                          .format(ipipe, currpipe['kp']))

        _logger.debug('8. Display interim results')
        print('Iteration {0:d}'.format(nct))
        print('Deviation {0:0.4E} (Tolerance {1:0.4E}'
              .format(ssum, case_dom['params']['tolerance']))

#        print()
#        print('Pipe   Kp            expp          Qcurrent                  '
#              'Qpredict')
#        for ipipe, currpipe in enumerate(case_dom['pipe']):
#            print('{0:3d}    {1:0.4E}    {2:0.4E}    {3:0.4E~}    {4:0.4E~}'
#                  .format(ipipe,
#                          currpipe['kp'],
#                          currpipe['expp'],
#                          Q_(x[ipipe], 'm**3/s').to('ft**3/s'),
#                          Q_(qpredict[ipipe], 'm**3/s').to('ft**3/s')))
#        print()

        update_vol_flows(case_dom, dq)
        set_pipe_head_loss(case_dom, dq)

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
    return


def set_pipe_derived_properties(case_dom):
    """Set derived properties of pipe - flow area, and initial kp and expp

    Args:
        case_dom (dict): Pipe network data model"""
#    ugrav = Q_(sc.g, 'm/s**2')

    for currpipe in case_dom['pipe']:
        currpipe['loopdir'] = {}
        currpipe['kp'] = (ahws_us * currpipe['lpipe'].to('ft').magnitude
                          / (currpipe['chw']**echw
                             * currpipe['idiameter'].to('ft').magnitude**edhw))

        currpipe['flow_area'] = (
            sc.pi / 4.0 * currpipe['idiameter']**2
        ).to('ft**2')

        currpipe['expp'] = 0.0

        _logger.debug('Pipe {:d}: '
                      'Aflow = {:0.4E~}, '
                      'Kp = {:0.4E}, '
                      'expp = {:0.4E}'
                      .format(currpipe['id'],
                              currpipe['flow_area'],
                              currpipe['kp'],
                              currpipe['expp']))

    for iloop, currloop in enumerate(case_dom['loop']):
        for pipeloop in currloop:
            pid = pipeloop['pipe_id']
            case_dom['pipe'][pid]['loopdir'][iloop] = pipeloop['flow_dir']

    return


def pipe_dimension_table(pipelist):
    """Print pipe dimensions in a tabular text format

    Args:
        pipelist ([dict]): List of dicts containing pipe data

    Returns:
        (str): Text table of pipe dimensions"""

    result = 'Pipe    Diameter         Length           CHW\n'
    for currpipe in pipelist:
        result += '{0:3d}     {1:12.4E~}    {2:12.4E~}  {3:8.1f}\n' \
                  .format(currpipe['id'], currpipe['idiameter'],
                          currpipe['lpipe'], currpipe['chw'])
    return result


def set_pipe_head_loss(case_dom, dq):
    """Find head loss along a pipe

    Args:
        case_dom (dict): Pipe flow network data model
        dq (numpy.ndarray): Array of loop flow corrections
    """
    for currpipe in case_dom['pipe']:
        currpipe['head_loss'] = Q_(currpipe['kp']
                                   * abs(currpipe['vol_flow'].to('ft**3/s')
                                         .magnitude)**echw, 'ft')


def flow_and_head_loss_report(case_dom):
    """Return a string containing the final calculated volumetric flow and head
    loss in each pipe

    Args:
        case_dom (dict): Pipe network data model

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
        case_dom (dict): Pipe network data model
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
    _logger.info("Starting jeppson_ch6b")

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

            set_pipe_derived_properties(case_dom)

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
#
#            dotfn = abspath((splitext(fh.name))[0] + '_{:d}.gv'.format(icase))
#            create_topology_dotfile(case_dom, dotfn)

#            print('case_dom:')
#            _pp.pprint(case_dom)

            _logger.info('Done processing case {:d}'.format(icase))
        _logger.info('Done processing {0:s}'.format(fh.name))
    _logger.info("Ending jeppson_ch6b")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
