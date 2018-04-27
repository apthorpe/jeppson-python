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
              idiameter=12.0 * sc.inch, eroughness=8.5E-4)
    assert p1.label == 'Pipe 1'
    assert p1.length == approx(30.48)
    assert p1.idiameter == approx(0.3048)
    assert p1.eroughness == approx(8.5E-4)

    assert p1.as_table(headers=['length', 'idiameter', 'eroughness'])

    p2 = Pipe(label='Pipe 2', length=100.0 * sc.foot, 
              idiameter=12.0 * sc.inch, eroughness=8.5E-4)

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
        
    with raises(ValueError):
        p2.nearest_material_roughness('cheese', is_clean=False)

    p2.nearest_material_roughness('cast iron', is_clean=True)
    assert p2.eroughness == approx(2.59E-4)

    p2.nearest_dimensions_from_schedule('80')
    # The schedule 80 pipe size with an inner diameter closest to 12" is
    # 14" with Di = 12.75", Do = 14.0, and twall = 0.75"
    assert p2._twall == approx(0.75 * sc.inch)
    assert p2._odiameter == approx(14.0 * sc.inch)
    assert p2.idiameter == approx(12.5 * sc.inch)

    with raises(ValueError):
        p2.nearest_dimensions_from_schedule(schedule='80', dnominal=120)

    p2.nearest_dimensions_from_schedule(schedule='80', dnominal=12)
    assert p2.idiameter == approx(11.38 * sc.inch, abs=1.0E-3)
    assert p2._odiameter == approx(12.75 * sc.inch, abs=1.0E-3)
    assert p2._twall == approx((p2._odiameter - p2.idiameter) / 2.0)

    p2.idiameter = 36.0 * sc.inch
    with raises(ValueError):
        p2.nearest_dimensions_from_schedule(schedule='80')

    with raises(ValueError):
        p2.clean = 'manky'

    with raises(ValueError):
        p2.surface = 'smooth'

    # 12" Schedule 80 pipe has a Di of 11.38", Do = 12.75", and twall = 0.68"
    p3 = Pipe(label='Pipe 3', length=10.0 * sc.foot, 
              idiameter=11.37 * sc.inch, schedule='80',
              eroughness=8.5E-4)

    assert p3.idiameter == approx(11.38 * sc.inch, abs=1.0E-3)
    assert p3._odiameter == approx(12.75 * sc.inch, abs=1.0E-3)
    assert p3._twall == approx((p3._odiameter - p3.idiameter) / 2.0)

    # 12" Schedule 80 pipe has a Di of 11.38", Do = 12.75", and twall = 0.68"
    p4 = Pipe(label='Pipe 4', length=10.0 * sc.foot, 
              idiameter=11.37 * sc.inch, 
              odiameter=12.75 * sc.inch, 
              eroughness=8.5E-4)

    assert p4._twall == approx((p4._odiameter - p4.idiameter) / 2.0)

    # 12" Schedule 80 pipe has a Di of 11.38", Do = 12.75", and twall = 0.68"
    p5 = Pipe(label='Pipe 5', length=10.0 * sc.foot, 
              idiameter=11.37 * sc.inch, 
              twall=0.68 * sc.inch, 
              eroughness=8.5E-4)

    assert p5._odiameter == approx(12.75 * sc.inch, abs=1.0E-3)

    # 12" Schedule 80 pipe has a Di of 11.38", Do = 12.75", and twall = 0.68"
    p6 = Pipe(label='Pipe 6', length=10.0 * sc.foot, 
              schedule='80', nps=12)

    assert p6.idiameter == approx(11.38 * sc.inch, abs=1.0E-3)
    assert p6._odiameter == approx(12.75 * sc.inch, abs=1.0E-3)
    assert p6._twall == approx((p6._odiameter - p6.idiameter) / 2.0)
    assert p6.eroughness == 0.0
    assert p6.as_table(headers=['length', 'idiameter', 'eroughness'])

    # 12" Schedule 80 pipe has a Di of 11.38", Do = 12.75", and twall = 0.68"
    p7 = Pipe(label='Pipe 7', length=10.0 * sc.foot, 
              schedule='80', nps=12, froughness=2.46850408E-4)

    assert p7.idiameter == approx(11.38 * sc.inch, abs=1.0E-3)
    assert p7._odiameter == approx(12.75 * sc.inch, abs=1.0E-3)
    assert p7._twall == approx((p7._odiameter - p7.idiameter) / 2.0)
    assert p7.flow_area == approx(6.5524E-2, abs=1.0E-6)
    assert p7.eroughness == approx(8.54E-4, abs=1.0E-6)

    p7 = Pipe(label='Pipe 7', length=10.0 * sc.foot, 
              schedule='80', nps=12, surface='cast iron')

    assert p7.eroughness == approx(2.59E-4)

    p7 = Pipe(label='Pipe 7', length=10.0 * sc.foot, 
              schedule='80', nps=12, surface='Steel tubes', is_clean=False)

    assert p7.eroughness == approx(1.0E-3)

    return
