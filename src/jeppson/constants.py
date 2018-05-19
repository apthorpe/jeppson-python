"""Hazen-Williams correlation constants needed by Jeppson applications.

The Hazen-Williams correlation constants are taken from Chapter 2 of *Steady
Flow Analysis of Pipe Networks: An Instructional Manual* and *Analysis of Flow
in Pipe Networks* by Roland W. Jeppson. Note that the Hazen-Williams
correlation is *dimensional*, not dimensionless. This requires particular care
when selecting leading coefficients for the flow and head loss
correrlations."""

from __future__ import absolute_import, division, print_function

import logging

from . import _logger, ureg, Q_

from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"

#: Exponent on *S* term of Hazen-Williams flow correlation, 0.54
eshw = 0.54

#: Exponent on *R* term of Hazen-Williams flow correlation, 0.63
erhw = 0.63

#: Exponent on *C* and *Q* term of Hazen-Williams head loss correlationm, about
#: 1.852
echw = 1.0 / eshw

#: Exponent on *D* term of Hazen-Williams head loss correlation, about 4.87
edhw = (2.0 + erhw) / eshw

#: Leading coefficient on Hazen-Williams head loss correlation, SI units, 10.67
ahws_si = 10.67

#: Leading coefficient on Hazen-Williams head loss correlation, US traditional
#: units (feet), 4.73
ahws_us = 4.73

#: Leading coefficient on Hazen-Williams flow correlation, SI units, 0.849
ahwq_si = 0.849

#: Leading coefficient on Hazen-Williams flow correlation, US traditional units
#: (feet), 1.318
ahwq_us = 1.318
