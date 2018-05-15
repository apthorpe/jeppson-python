#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This reimplements the first pipe flow analysis program given in chapter 6 of
*Steady Flow Analysis of Pipe Networks: An Instructional Manual* (1974).
*Reports.* Paper 300.  http://digitalcommons.usu.edu/water_rep/300 and
*Analysis of Flow in Pipe Networks* (1976). Ann Arbor Science Publishers, Inc.
http://www.worldcat.org/title/analysis-of-flow-in-pipe-networks/oclc/927534147
by Roland W. Jeppson.

Then run `python setup.py install` which will install the command
`jeppson_ch6a` inside your current environment.
"""
from __future__ import absolute_import, division, print_function

import argparse
from collections import OrderedDict
import logging
# import pprint
from math import copysign
from os.path import abspath, splitext
import sys

import numpy as np
import pygraphviz as pgv

from . import _logger, Q_
from jeppson.input import InputLine
from jeppson.constants import ahws_us, eshw, echw, edhw

from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"

# _pp = pprint.PrettyPrinter(indent=4)


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
    tags = ('npipes', 'njunctions', 'jfixed', 'maxiter', 'tolerance')
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
        if tag == 'tolerance':
            case_params[tag] = float(cpline.token[itok])
        else:
            case_params[tag] = int(cpline.token[itok])

    if (case_params['jfixed'] < 1
            or case_params['jfixed'] > case_params['njunctions']):
        msg = 'ID of fixed junction out of range; found {0:d}, expected ' \
              '[1 .. {1:d}]'.format(case_params['jfixed'],
                                    case_params['njunctions'])
        raise ValueError(msg)
    else:
        # Reduce junction ID by one (zero-indexed)
        case_params['jfixed'] -= 1

    _logger.debug('Successfully read case parameters')
    _logger.debug('  npipes = {0:d}'
                  .format(case_params['npipes']))
    _logger.debug('  njunctions = {0:d}'
                  .format(case_params['njunctions']))
    _logger.debug('  jfixed = {0:d} (zero-indexed)'
                  .format(case_params['jfixed']))
    _logger.debug('  maxiter = {0:d}'
                  .format(case_params['maxiter']))
    _logger.debug('  tolerance = {0:0.4E}'
                  .format(case_params['tolerance']))

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

    tags = ('from', 'to', 'idiameter', 'chw', 'lpipe')
    mintok = len(tags)

    pipe_info = []

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

        pipe_info.append({
            'id': ipipe,
            'from': int(currline.token[0]) - 1,
            'to': int(currline.token[1]) - 1,
            'idiameter': Q_(float(currline.token[2]), 'in'),
            'chw': float(currline.token[3]),
            'lpipe': Q_(float(currline.token[4]), 'foot')
        })

    for currpipe in pipe_info:
        _logger.debug('  Pipe {0:d} from {1:d} to {2:d} D={3:4.1f~} '
                      'L={4:8.1E~} CHW={5:9.1f}'
                      .format(currpipe['id'],
                              currpipe['from'],
                              currpipe['to'],
                              currpipe['idiameter'],
                              currpipe['lpipe'],
                              currpipe['chw']))

    return pipe_info


def extract_junctions(deck, iptr, njunctions):
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

    tags = ('id', 'inflow', 'init_head')
    mintok = len(tags)
    junc_info = []

    for ijunc in range(njunctions):
        junc_info.append({'id': ijunc})

    for idx in range(njunctions):
        currline = deck[iptr + idx]

        if currline.ntok < mintok:
            msg = 'Too few tokens found in junction definition line ({0:d} ' \
                'found, {1:d} expected)'.format(currline.ntok, mintok)
            _logger.error(msg)
            _logger.error(currline.as_log())
            raise ValueError(msg)

        if currline.ntok > mintok:
            msg = 'Too many tokens found in junction definition line ' \
                '({0:d} found, {1:d} expected)'.format(currline.ntok, mintok)
            _logger.warning(msg)

        ijunc = int(currline.token[0]) - 1
        if ijunc < 0 or ijunc >= njunctions:
            msg = 'Junction identifier out of range (found {0:s}, ' \
                  'expected [1 .. {1:d}])' \
                  .format(currline.token[0], njunctions)
            _logger.error(msg)
            _logger.error(currline.as_log())

            raise ValueError(msg)
        else:
            junc_info[ijunc]['inflow'] = Q_(float(currline.token[1]),
                                            'ft**3/sec')
            junc_info[ijunc]['init_head'] = Q_(float(currline.token[2]), 'ft')

    for currjunc in junc_info:
        _logger.debug('  Junction {0:d} Qin={1:8.3f~} H={2:7.2f~}'
                      .format(currjunc['id'],
                              currjunc['inflow'],
                              currjunc['init_head']))

    return junc_info


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

    # Step 3. Read junction data
    _logger.debug('3. Reading junction inflow and head')
    njunctions = case_dom['params']['njunctions']
    case_dom['junc'] = extract_junctions(deck, iptr, njunctions)

    iptr += njunctions

    # Done reading input; assert (iptr + iread - 1) == len(deck) for a
    # single case
    case_dom['params']['next_iptr'] = iptr

    return case_dom


def calculate_head_loss(case_dom):
    """Calculate head loss along pipes based on junction heads

    Args:
        case_dom (dict): Pipe network data model"""

    for pipe in case_dom['pipe']:
        pipe['head_loss'] = (case_dom['junc'][pipe['from']]['head']
                             - case_dom['junc'][pipe['to']]['head'])

        abs_flow = (abs(pipe['head_loss'].to('ft').magnitude)
                    / pipe['kp']**eshw)
        pipe['vol_flow'] = Q_(copysign(abs_flow, pipe['head_loss'].magnitude),
                              'ft**3/s')


def solve_network_flows(case_dom):
    """Find the volumetric flow and head loss for the piping network defined in
    the case_dom structure by using the first method described in Chapter 6 of
    Jeppson.

    Args:
        case_dom (dict): Pipe flow network data model

    Raises:
        ValueError: Network solution matrix is singular or does not converge.
    """

# The goal of this project is to reimplement the JEPPSON_CH6A code in Python,
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
    njunctions = case_dom['params']['njunctions']

    for junc in case_dom['junc']:
        junc['head'] = junc['init_head']

    calculate_head_loss(case_dom)

# In some universe, the line below is clearer than the four which follow it.
# Then again some people think Perl's map() syntax is clear too. The difference
# is that few people in the Perl community will give you grief for avoiding the
# use of map().
#
# For the purpose of satisfying those who demand at least one list
# comprehension per Python project (warranted or not), the following is
# provided. It lacks the clarity and robustness of the commented equivalent
# code (does not support diagnostics or error-handling), however it does uphold
# the Perl tradition of being inscrutable, ridden with punctuation, and best of
# all, it's crammed into a single line.

# Ref: https://google.github.io/styleguide/pyguide.html#List_Comprehensions

    jfixed = case_dom['params']['jfixed']

    # Junction ID from matrix index
    livejunc = [i for i in range(njunctions) if i != jfixed]

    # Matrix index from junction ID
    rcid_junc = OrderedDict()
    for idx, ijunc in enumerate(livejunc):
        rcid_junc[ijunc] = idx

    assert jfixed not in rcid_junc

    done = False
    converged = False

    while not done:
        # Step 5. Assemble matrix
        _logger.debug('5a. Assemble matrix rows of junction equations')

        a = np.zeros((njunctions-1, njunctions-1))
        b = np.zeros((njunctions-1))

        for pipe in case_dom['pipe']:
            ifrom = pipe['from']
            ito = pipe['to']

            arg = ((case_dom['junc'][ifrom]['head']
                    - case_dom['junc'][ito]['head']).to('ft').magnitude
                   / pipe['kp'])

            if ifrom in rcid_junc:
                idx = rcid_junc[ifrom]
                arge = arg**eshw
                b[idx] += arge
                z = arge * eshw / (pipe['kp'] * arg)
                a[idx, idx] += z
                if ito in rcid_junc:
                    jdx = rcid_junc[ito]
                    a[jdx, idx] -= z

            if ito in rcid_junc:
                idx = rcid_junc[ito]
                arge = -(arg**eshw)
                b[idx] += arge
                z = arge * eshw / (pipe['kp'] * arg)
                a[idx, idx] -= z
                if ifrom in rcid_junc:
                    jdx = rcid_junc[ifrom]
                    a[jdx, idx] += z

        for idx, ijunc in enumerate(livejunc):
            b[idx] -= case_dom['junc'][ijunc]['inflow'] \
                      .to('ft**3/s').magnitude

#        print('a: ' + repr(a))
#        print('b: ' + repr(b))

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

            raise ValueError(msg)

#        print('x: ' + repr(x))

        _logger.debug('7. Adjust matrix')

        ssum = 0.0
        for idx, ijunc in enumerate(livejunc):
            ssum += abs(x[idx])
            # Add corrective pressures
            case_dom['junc'][ijunc]['head'] -= Q_(x[idx], 'ft')

        calculate_head_loss(case_dom)

        recalc_head = False
        for ipipe, currpipe in enumerate(case_dom['pipe']):
            if currpipe['head_loss'].magnitude < 0.0:
                recalc_head = True
#                currpipe['head_loss'] = -currpipe['head_loss']
                ifrom = currpipe['from']
                ito = currpipe['to']
                print('Flow reversal detected in pipe {:d} from {:d} to {:d}'
                      .format(ipipe, ifrom, ito))
                currpipe['to'] = ifrom
                currpipe['from'] = ito

        if recalc_head:
            calculate_head_loss(case_dom)
            print()

        _logger.debug('8. Display interim results')
        print('Iteration {0:d}'.format(nct))
        print('Deviation {0:0.4E} (Tolerance {1:0.4E}'
              .format(ssum, case_dom['params']['tolerance']))

        print()
        print('Pipe   Kp            Head loss        Volumetric Flow')
        for ipipe, currpipe in enumerate(case_dom['pipe']):
            print('{0:3d}    {1:0.4E}    {2:0.4E~}    {3:0.4E~}'
                  .format(ipipe,
                          currpipe['kp'],
                          currpipe['head_loss'].to('ft'),
                          currpipe['vol_flow'].to('ft**3/s')))
        print()
        print('Junction    Head')
        for ijunc, currjunc in enumerate(case_dom['junc']):
            print('{0:3d}         {1:0.4E~}'
                  .format(ijunc, currjunc['head'].to('ft')))
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

    return


def set_pipe_derived_properties(pipelist):
    """Set derived properties of pipe - flow area, arl, eroughness, and initial
    kp and expp

    Args:
        pipelist ([dicr]): List of dicts containing pipe dimensions and
          attributes"""
    for currpipe in pipelist:
        currpipe['LD'] = (
            currpipe['lpipe'] / currpipe['idiameter']
        ).to_base_units().magnitude

        currpipe['kp'] = (ahws_us * currpipe['lpipe'].to('ft').magnitude
                          / (currpipe['chw']**echw
                             * currpipe['idiameter'].to('ft').magnitude**edhw))

        _logger.debug('Pipe {:d}: kp = {:0.4E}'
                      .format(currpipe['id'], currpipe['kp']))

    return


def pipe_dimension_table(pipelist):
    """Print pipe dimensions in a tabular text format

    Args:
        pipelist ([dict]): List of dicts containing pipe data

    Returns:
        (str): Text table of pipe dimensions"""

    result = 'Pipe    Diameter         Length         CHW\n'
    for currpipe in pipelist:
        result += '{0:3d}     {1:12.4E~}    {2:12.4E~}  {3:7.3f}\n' \
                  .format(currpipe['id'], currpipe['idiameter'],
                          currpipe['lpipe'], currpipe['chw'])
    return result


def pipe_results_table(case_dom):
    """Print pipe dimensions in a tabular text format

    Args:
        case_dom (dict): Pipe flow network data model

    Returns:
        (str): Text table of pipe dimensions, routing, and results"""

    result = '  Pipe  From    To  Diameter         Length         CHW' \
             '      Volumetric Flow    Head Loss      From Head      To Head\n'
    for currpipe in case_dom['pipe']:
        ifrom = currpipe['from']
        ito = currpipe['to']
        result += '{0:3d}   {1:3d}     {2:3d}   {3:12.4E~}    {4:12.4E~}  ' \
                  '{5:7.3f}  {6:7.3f~}  {7:12.4E~}  {8:12.4E~}  {9:12.4E~}\n' \
                  .format(currpipe['id'],
                          currpipe['from'],
                          currpipe['to'],
                          currpipe['idiameter'],
                          currpipe['lpipe'],
                          currpipe['chw'],
                          currpipe['vol_flow'],
                          currpipe['head_loss'],
                          case_dom['junc'][ifrom]['head'],
                          case_dom['junc'][ito]['head'])
    return result


def junction_results_table(case_dom):
    """Generate a text table of junction results

    Args:
        case_dom (dict): Pipe flow network data model

    Returns:
        (str): Table of junction results"""

    results = 'Junction  Inflow                  Head\n'
    for ijunc, currjunc in enumerate(case_dom['junc']):
        if ijunc == case_dom['params']['jfixed']:
            note = 'Fixed head'
        else:
            note = ''

        results += '{0:3d}       {1:12.4E~}  {2:-12.4E~}  {3:s}\n' \
                   .format(currjunc['id'],
                           currjunc['inflow'],
                           currjunc['init_head'],
                           note)

    return results


def create_topology_dotfile(case_dom, filepath='tmp.gv'):
    """Create directed graph in GraphViz ``dot`` format

    Args:
        case_dom (dict): Pipe flow network data model
        filepath (str): Absolute path of dotfile"""

    jfmt = 'J{:d}'
    jhfmt = 'J{:d}\n{:0.3f~}'
    jxfmt = 'JX{:d}'
#    pfmt = 'P{:d}'
    pqfmt = 'P{:d}: {:0.3f~}'

    G = pgv.AGraph(directed=True, splines=False, ratio='fill', overlap=False)

    for idx, junc in enumerate(case_dom['junc']):
        jtag = jfmt.format(idx)
        jhtag = jhfmt.format(idx, junc['head'].to('ft'))
        jxtag = jxfmt.format(idx)
        G.add_node(jtag, shape='circle', label=jhtag)

        if idx == case_dom['params']['jfixed']:
            ftag = '(fixed)'
            style = 'dashed'
        else:
            ftag = ''
            style = 'solid'

        if junc['inflow'].magnitude > 0.0:
            G.add_node(jxtag, shape='none', label='')
            G.add_edge(jxfmt.format(idx), jtag, color='blue', style=style,
                       label='{:0.3f~} {:s}'.format(junc['inflow'], ftag))
        elif junc['inflow'].magnitude < 0.0:
            G.add_node(jxtag, shape='none', label='')
            G.add_edge(jtag, jxfmt.format(idx), color='red', style=style,
                       label='{:0.3f~} {:s}'.format(abs(junc['inflow']), ftag))

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
    _logger.info("Starting jeppson_ch6a")

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

            failed = False
            try:
                solve_network_flows(case_dom)
            except ValueError as err:
                failed = True
                if str(err).startswith('Cannot solve'):
                    _logger.error('Failed to solve case {0:d} '
                                  'from {1:s}: {2:s}'
                                  .format(icase, fh.name, str(err)))
                    _logger.info('Advancing to next case.')
                    continue

            # Step 10. Display results
            _logger.debug('10. Display final results')

            print('Final results:\n')
            if failed:
                print('WARNING: Case {0:d} did not converge.\n'
                      .format(icase))

            print(pipe_results_table(case_dom))
            print()
            print(junction_results_table(case_dom))

            dotfn = abspath((splitext(fh.name))[0] + '_{:d}.gv'.format(icase))
            create_topology_dotfile(case_dom, dotfn)

#            print('case_dom:')
#            _pp.pprint(case_dom)

            _logger.info('Done processing case {:d}'.format(icase))
        _logger.info('Done processing {0:s}'.format(fh.name))
    _logger.info("Ending jeppson_ch6a")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
