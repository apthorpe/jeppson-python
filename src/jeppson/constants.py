"""Constants needed by Jeppson applications such as Hazen-Williams correlation
contants"""

from __future__ import absolute_import, division, print_function

import logging

from . import _logger, ureg, Q_

from jeppson import __version__

__author__ = "Bob Apthorpe"
__copyright__ = "Bob Apthorpe"
__license__ = "mit"

# Exponent on \f[S\f] term of Hazen-Williams flow correlation, 0.54
eshw = 0.54

# Exponent on \f[R\f] term of Hazen-Williams flow correlation, 0.63
erhw = 0.63

# Exponent on \f[C\f] and \f[Q\f] term of Hazen-Williams head loss
# correlationm, about 1.852
echw = 1.0 / eshw

# Exponent on \f[D\f] term of Hazen-Williams head loss correlation, about
# 4.87
edhw = (2.0 + erhw) / eshw

# Leading coefficient on Hazen-Williams head loss correlation,
# US traditional units (feet), 4.73
ahws_us = 4.73

