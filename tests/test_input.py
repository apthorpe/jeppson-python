#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pytest import approx, raises
from jeppson.input import InputLine

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"


def test_input():
    blank = (
        '',
        '\r',
        '\r\n',
        '                                        '
            '                                        \r\n',
        '\n',
        '\r\n\r\n',
        ' \n',
        '\t\n',
        '\t   \v\n',
    )

    for line in blank:
        iline = InputLine(line)
        tline = line.rstrip('\n\r')
        assert iline.type == 'blank'
        assert iline.typecode == 'B'
        assert iline.ntok == 0
        assert len(iline.token) == 0
        assert iline.line == tline
        assert iline.as_log()
        
    # Comments:
    comment = (
        '#\r\n',
        '                    #                   '
            '                                        \r',
        '                                        '
            '                                       #\n',
        '! Alternate comment character\r\n',
        '##### Just a normal comment\n',
        'an abnormal # non-comment line\n'
    )

    icmt = 0
    line = comment[icmt]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'comment'
    assert iline.typecode == 'C'
    assert iline.ntok == 0
    assert len(iline.token) == 0
    assert iline.line == tline
    assert iline.as_log()

    icmt += 1
    line = comment[icmt]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'comment'
    assert iline.typecode == 'C'
    assert iline.ntok == 0
    assert len(iline.token) == 0
    assert iline.line == tline
    assert iline.as_log()

    icmt += 1
    line = comment[icmt]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'comment'
    assert iline.typecode == 'C'
    assert iline.ntok == 0
    assert len(iline.token) == 0
    assert iline.line == tline
    assert iline.as_log()

    icmt += 1
    line = comment[icmt]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type != 'comment'
    assert iline.typecode != 'C'
    assert iline.ntok != 0
    assert len(iline.token) != 0
    assert iline.line == tline
    assert iline.as_log()

    line = comment[icmt]
    iline = InputLine(line, commentchar='!')
    tline = line.rstrip('\n\r')
    assert iline.type == 'comment'
    assert iline.typecode == 'C'
    assert iline.ntok == 0
    assert len(iline.token) == 0
    assert iline.line == tline
    assert iline.as_log()

    icmt += 1
    line = comment[icmt]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'comment'
    assert iline.typecode == 'C'
    assert iline.ntok == 0
    assert len(iline.token) == 0
    assert iline.line == tline
    assert iline.as_log()

    icmt += 1
    line = comment[icmt]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type != 'comment'
    assert iline.typecode != 'C'
    assert iline.ntok != 0
    assert len(iline.token) != 0
    assert iline.line == tline
    assert iline.as_log()

# Data
    data = (
        '0    3    -1   -2    5',
        '   0.33333   0.497   150.0      1.217E-5    7.0E-6  32.2    ',
        '-1000.0',
        '    1   0.0     200.0    ',
        '   10   8.     1200.0     130.        2.0    ',
        '1106.0    751.0     1000.0    500.0     1200.0     600.0     800.0',
        '  The wheels 4.0 turn slowly -1.9E-17 '
            'but they grind 8 exceedingly fine.#'
    )

    idat = 0
    line = data[idat]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'data'
    assert iline.typecode == 'D'
    assert iline.ntok == 5
    assert len(iline.token) == 5
    assert iline.line == tline
    assert iline.as_log()

    idat += 1
    line = data[idat]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'data'
    assert iline.typecode == 'D'
    assert iline.ntok == 6
    assert len(iline.token) == 6
    assert iline.line == tline
    assert iline.as_log()

    idat += 1
    line = data[idat]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'data'
    assert iline.typecode == 'D'
    assert iline.ntok == 1
    assert len(iline.token) == 1
    assert iline.line == tline
    assert iline.as_log()

    idat += 1
    line = data[idat]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'data'
    assert iline.typecode == 'D'
    assert iline.ntok == 3
    assert len(iline.token) == 3
    assert iline.line == tline
    assert iline.as_log()

    idat += 1
    line = data[idat]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'data'
    assert iline.typecode == 'D'
    assert iline.ntok == 5
    assert len(iline.token) == 5
    assert iline.line == tline
    assert iline.as_log()

    idat += 1
    line = data[idat]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'data'
    assert iline.typecode == 'D'
    assert iline.ntok == 7
    assert len(iline.token) == 7
    assert iline.line == tline
    assert iline.as_log()

    idat += 1
    line = data[idat]
    iline = InputLine(line)
    tline = line.rstrip('\n\r')
    assert iline.type == 'data'
    assert iline.typecode == 'D'
    assert iline.ntok == 12
    assert len(iline.token) == 12
    assert iline.line == tline
    assert iline.as_log()

    return
