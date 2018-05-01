# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound
from pint import UnitRegistry

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = 'jeppson-python'
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'

import logging
# Set default logging handler to avoid "No handler found" warnings.
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

_logger = logging.getLogger(__name__)
_logger.addHandler(NullHandler())

# Pint unit registry
# Share with all functions in an application  via 'from . import ureg, Q_
ureg = UnitRegistry()
Q_ = ureg.Quantity
