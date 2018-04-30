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
        assert len(iline.token) == 0
        assert iline.line == tline

    return
