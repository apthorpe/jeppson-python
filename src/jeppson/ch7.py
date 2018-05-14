#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This reimplements the Hardy Cross pipe flow analysis program given in chapter 7
of *Steady Flow Analysis of Pipe Networks: An Instructional Manual* (1974).
*Reports.* Paper 300.  http://digitalcommons.usu.edu/water_rep/300 and
*Analysis of Flow in Pipe Networks* (1976). Ann Arbor Science Publishers, Inc.
http://www.worldcat.org/title/analysis-of-flow-in-pipe-networks/oclc/927534147
by Roland W. Jeppson.

Then run `python setup.py install` which will install the command `jeppson_ch7`
inside your current environment.
"""
from __future__ import division, print_function, absolute_import

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
import pygraphviz as pgv

from . import _logger, ureg, Q_
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
    tags = ('npipes', 'nloops', 'maxiter', 'kin_visc', 'tolerance', 'grav2')
    mintok = len(tags)

    case_params = {
    }

    cpline = deck[iptr]
    if cpline.ntok < mintok:
        msg = 'Too few entries found for case parameters ({0:d} found, )' \
              '{1:d} expected)'.format(cpline.ntok, mintok)
        raise ValueError(msg)
    elif cpline.ntok > mintok:
        msg = 'Too many entries found for case parameters ({0:d} found, ' \
              '{1:d} expected)'.format(cpline.ntok, mintok)
        _logger.warning(msg)

    case_params['npipes'] = int(cpline.token[0])
    case_params['nloops'] = int(cpline.token[1])
    case_params['maxiter'] = int(cpline.token[2])
    case_params['kin_visc'] = Q_(float(cpline.token[3]), 'ft**2/s')
    case_params['tolerance'] = float(cpline.token[4])
    case_params['grav2'] = Q_(float(cpline.token[5]), 'ft/s**2')

    _logger.debug('Successfully read case parameters')
    _logger.debug('  {0:s} = {1:d}'
                  .format('npipes', case_params['npipes']))
    _logger.debug('  {0:s} = {1:d}'
                  .format('nloops', case_params['nloops']))
    _logger.debug('  {0:s} = {1:d}'
                  .format('maxiter', case_params['maxiter']))
    _logger.debug('  {0:s} = {1:0.4E~}'
                  .format('kin_visc', case_params['kin_visc']))
    _logger.debug('  {0:s} = {1:0.4E}'
                  .format('tolerance', case_params['tolerance']))
    _logger.debug('  {0:s} = {1:0.4E~}'
                  .format('grav2', case_params['grav2']))

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
        pipe_info[pipe_id]['froughness'] = Q_(float(currline.token[3]), 'in')
        pipe_info[pipe_id]['init_vol_flow'] = Q_(float(currline.token[4]),
                                                 'ft**3/sec')

    for currpipe in pipe_info:
        _logger.debug('  Pipe {0:d} D={1:4.1f~} L={2:8.1E~} '
                      'froughness={3:0.4E~} Qinit={4:8.2f~}'
                      .format(currpipe['id'],
                              currpipe['idiameter'],
                              currpipe['lpipe'],
                              currpipe['froughness'],
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

    # Hardcoded; should be user input
    case_dom['params']['fvol_flow'] = 0.1

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

    # Done reading input
    case_dom['params']['next_iptr'] = iptr

    return case_dom


def solve_network_flows(case_dom):
    """Find the volumetric flow and head loss for the piping network defined in
    the case_dom structure by using the first method described in Chapter 6 of
    Jeppson.

    Args:
        case_dom (dict): Pipe flow network object model

    Raises:
        ValueError: Network solution matrix is singular or does not converge.
    """

# The goal of this project is to reimplement the JEPPSON_CH7 code in Python,
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

    nct = 0
    ssum = 100.0
    kin_visc = case_dom['params']['kin_visc']

    # Set kp and expp for each pipe
    for pipe in case_dom['pipe']:
        qm = abs(pipe['init_vol_flow'])
        pipe['vol_flow'] = qm
        dq = case_dom['params']['fvol_flow'] * qm

        q1 = qm - dq
        q2 = qm + dq
        v1 = q1 / pipe['flow_area'].to('ft**2')
        v2 = q2 / pipe['flow_area'].to('ft**2')
        Re1 = (v1 * pipe['idiameter'].to('ft') / kin_visc) \
            .to_base_units().magnitude
        Re2 = (v2 * pipe['idiameter'].to('ft') / kin_visc) \
            .to_base_units().magnitude
        f1 = friction_factor(Re=Re1, eD=pipe['eroughness'])
        f2 = friction_factor(Re=Re2, eD=pipe['eroughness'])

        _logger.debug('Pipe {0:d} Re = {1:0.4E}'.format(pipe['id'], Re1))
        if Re1 < 2050.0:
            # laminar
            pipe['expp'] = 1.0
            pipe['kp'] = (case_dom['params']['grav2'] * kin_visc
                          * pipe['arl'] / pipe['idiameter'].to('ft')
                          ).magnitude

        else:
            # Transition/tubulent
            be = ((log(f1) - log(f2))
                  / (log(q2.to('ft**3/s').magnitude)
                     - log(q1.to('ft**3/s').magnitude)))
            ae = f1 * (q1.to('ft**3/s').magnitude)**be

            pipe['expp'] = 2.0 - be
            pipe['kp'] = (ae * pipe['arl']).magnitude

        _logger.debug('Pipe {0:d} kp = {1:0.4E} expp = {2:0.4f}'
                      .format(pipe['id'], pipe['kp'], pipe['expp']))

        _logger.debug('Pipe {0:d} ae = {1:0.4f} be = {2:0.4f}'
                      .format(pipe['id'], ae, be))

    done = False
    converged = False

    while not done:
        # Step 5. Iteratively correct flows
        _logger.debug('Iteratively correct flows')
        ssum = 0.0
        for iloop, looppipe in enumerate(case_dom['loop']):
            sum1 = 0.0
            sum2 = 0.0
            for pipe in looppipe:
                pid = pipe['pipe_id']
                currpipe = case_dom['pipe'][pid]
                hl = ((pipe['flow_dir'] * currpipe['kp']
                       * currpipe['vol_flow']
                       .to('ft**3/s').magnitude**currpipe['expp']))
                currpipe['head_loss'] = Q_(abs(hl), 'ft')
                sum1 += hl
                _logger.debug('q = {0:0.4E~}'.format(currpipe['vol_flow']))
                _logger.debug('expp = {0:0.4E}'.format(currpipe['expp']))
                _logger.debug('hl = {0:0.4E~}'.format(currpipe['head_loss']))
                sum2 += (currpipe['expp'] * abs(hl)
                         / (currpipe['vol_flow'].to('ft**3/s').magnitude))

            dq = sum1 / sum2

            # 5c) Increment the deviation accumulator with the net flow
            # correction

            ssum = ssum + abs(dq)

            for pipe in looppipe:
                pid = pipe['pipe_id']
                currpipe = case_dom['pipe'][pid]
                currpipe['vol_flow'] -= Q_(pipe['flow_dir'] * dq, 'ft**3/s')

        converged = ssum <= case_dom['params']['tolerance']
        done = converged or (nct >= case_dom['params']['maxiter'])

    # End iteration
# ########################################################################
    if not converged:
        msg = 'Case not converged: ssum = {0:0.4E} > tolerance {1:0.4E}' \
              .format(ssum, case_dom['params']['tolerance'])
        # Advance to next case
        raise ValueError(msg)

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
        ).to('ft**2')

        currpipe['arl'] = (
            currpipe['lpipe'] / (2.0 * ugrav * currpipe['idiameter']
                                 * currpipe['flow_area']**2)
        ).to('s**2/ft**5')

        currpipe['eroughness'] = (
            currpipe['froughness']
            / currpipe['idiameter']
        ).to_base_units().magnitude

        currpipe['expp'] = 0.0
        currpipe['kp'] = 0.0

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
                  .format(currpipe['id'], currpipe['idiameter'].to('in'),
                          currpipe['lpipe'].to('ft'), currpipe['eroughness'])
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
    ureg.default_system = 'US'

    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.info("Starting jeppson_ch7")

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

#           dotfn = abspath((splitext(fh.name))[0] + '_{:d}.gv'.format(icase))
#           create_topology_dotfile(case_dom, dotfn)

#            print(repr(case_dom))

            _logger.info('Done processing case {:d}'.format(icase))
        _logger.info('Done processing {0:s}'.format(fh.name))
    _logger.info("Ending jeppson_ch7")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
