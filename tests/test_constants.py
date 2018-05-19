#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pytest import approx, raises
import jeppson.constants as jc

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"


def test_constants():
    assert jc.eshw == approx(0.54)

    assert jc.erhw == approx(0.63)

    assert jc.ahwq_us == approx(1.318)

    assert jc.ahwq_si == approx(0.849)

    assert jc.echw == approx(1.85185185)

    assert jc.edhw == approx(4.87037037)

    assert jc.ahws_us == approx(4.73)

    assert jc.ahws_si == approx(10.67)

    return
