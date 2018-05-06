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
from math import copysign
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
        _logger.debug('  errlimit = {0:0.4E}'
                      .format(results['case_parameters'].errlimit))
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

    inputconv = [
        {'idiameter': 'in', 'lpipe': 'ft', 'froughness': 'in'},
        {'idiameter': 'ft', 'lpipe': 'ft', 'froughness': 'ft'},
        {'idiameter': 'm', 'lpipe': 'm', 'froughness': 'm'},
        {'idiameter': 'cm', 'lpipe': 'm', 'froughness': 'cm'}
    ]
    inunit = inputconv[unitcode]

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
        results['pipe_info'].append({'idiameter': tmppipe['idiameter'][ict],
                                     'lpipe': tmppipe['lpipe'][ict],
                                     'froughness': tmppipe['froughness'][ict]})
        curpipe = results['pipe_info'][ict]
        _logger.debug('  id={0:d}  D= {1:16.4E~} L={2:9.1f~} e={3:16.4E~}'
                      .format(ict,
                              curpipe['idiameter'],
                              curpipe['lpipe'],
                              curpipe['froughness']))
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

    inputconv = [
        'skip',
        'gallon / minute',
        'ft**3 / sec',
        'm**3 / sec'
    ]

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
            qinflow = Q_(float(deck[jptr].token[0]), inputconv[inflow_units])
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

            iptr += pipemap_info['iread']

            # Step 4. Read loop data
            _logger.debug('4. Reading loop continuity data')
            nloops = case_info['case_parameters'].nloops
            pipeloop_info = extract_loops(deck, iptr, nloops)

            if pipemap_info['status'] == 'error':
                _logger.error('Error reading loop continuity data: {0:s}'
                              .format(pipeloop_info['msg']))
                _logger.info('Skipping remainder of {0:s}'.format(fh.name))
                break

            iptr += pipeloop_info['iread']
            # Done reading input; assert (iptr + iread - 1) == len(deck) for a
            # single case

            # Step 5. Assemble matrix
            _logger.debug('5. Assemble matrix')
            # Step 6. Solve matrix
            _logger.debug('6. Solve matrix')
            # Step 7. Display results
            _logger.debug('7. Display results')

        _logger.info('Done processing {0:s}'.format(fh.name))
    _logger.info("Ending jeppson_ch5")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
