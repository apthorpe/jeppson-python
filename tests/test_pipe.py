#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pytest import approx, raises
from jeppson.pipe import Pipe
import scipy.constants as sc

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"


def test_pipe():
    p1 = Pipe(label='Pipe 1', length=100.0 * sc.foot, 
              idiameter=12.0 * sc.inch, eroughness=8.5E-4, 
              zbottom=0.0, ztop=0.0)
    assert p1.label == 'Pipe 1'
    assert p1.length == approx(30.48)
    assert p1.idiameter == approx(0.3048)
    assert p1.eroughness == approx(8.5E-4)

    p2 = Pipe(label='Pipe 2', length=100.0 * sc.foot, 
              idiameter=12.0 * sc.inch, eroughness=8.5E-4, 
              zbottom=0.0, ztop=0.0)

    with raises(ValueError):
        p2.idiameter = 1.0E-4

    with raises(ValueError):
        p2.idiameter = 11.0

    with raises(ValueError):
        p2.eroughness = -1.0E-6

    with raises(ValueError):
        p2.eroughness = 0.2

    with raises(ValueError):
        p2.length = 1.0E-4

    with raises(ValueError):
        p2.length = 1001.0

    assert p2.idiameter == approx(0.3048)
    assert p2.flow_area == approx(7.2965877E-2)

    with raises(ValueError):
        p2.flow_area = 7.2965877E-2
        
    p2.nearest_material_roughness('cast iron', 'clean')
    assert p2.eroughness == approx(2.59E-4)

    p2.nearest_dimensions_from_schedule('80')
    # The schedule 80 pipe size with an inner diameter closest to 12" is
    # 14" with Di = 12.75", Do = 14.0, and twall = 0.75"
    assert p2._twall == approx(0.75 * sc.inch)
    assert p2._odiameter == approx(14.0 * sc.inch)
    assert p2.idiameter == approx(12.5 * sc.inch)

    p2.nearest_dimensions_from_schedule(schedule='80', dnominal=12)
    assert p2.idiameter == approx(11.38 * sc.inch, abs=1.0E-3)
    assert p2._odiameter == approx(12.75 * sc.inch, abs=1.0E-3)
    assert p2._twall == approx((p2._odiameter - p2.idiameter) / 2.0)

    # 12" Schedule 80 pipe has a Di of 11.38", Do = 12.75", and twall = 0.68"
    p3 = Pipe(label='Pipe 3', length=10.0 * sc.foot, 
              idiameter=11.37 * sc.inch, schedule='80',
              eroughness=8.5E-4, zbottom=0.0, ztop=0.0)

    assert p3.idiameter == approx(11.38 * sc.inch, abs=1.0E-3)
    assert p3._odiameter == approx(12.75 * sc.inch, abs=1.0E-3)
    assert p3._twall == approx((p3._odiameter - p3.idiameter) / 2.0)

#    print(p1.as_table(headers=['length', 'idiameter', 'eroughness']))
#    assert fib(1) == 1
#    assert fib(2) == 1
#    assert fib(7) == 13
#    with raises(AssertionError):
#        fib(-10)
    return
