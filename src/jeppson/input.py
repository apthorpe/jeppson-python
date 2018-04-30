"""Input processing utilities for reading legacy Jeppson code input"""

from __future__ import absolute_import, division, print_function

import logging
# Set default logging handler to avoid "No handler found" warnings.
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

LOG = logging.getLogger(__name__)
LOG.addHandler(NullHandler())


class InputLine(dict):
    """Categorized and tokenized user input for Jeppson Ch. 2 friction factor
        and head loss calculator.

        All attributes are immutable except ``ipos``.

        Attributes:
            line (str): Original line of input text with newline(s) removed
            type (str): One of 'blank', 'comment', or 'data'
            typecode (str): One of 'B', 'C', or 'D', corresponding to type
            ipos (int): Line number in original file (count starts at 1)
            ntok (int): Number of tokens found (0 except for 'data' lines)
            token ([str]): Tokens parsed from line ('data' lines only,
              otherwise empty)

        Args:
            line (str): Original line of input text with newline(s) removed
            ipos (int): Line number in original file
            commentchar (str): leading character/string denoting that a line is
              a comment
    """

    def __init__(self, line, ipos=0, commentchar='#'):
        """Constructor """
        self.ipos = ipos

        self._line = line.rstrip('\r\n')
        self._type = 'unknown'
        self._ntok = 0
        self._token = ()

        tline = self._line.strip()
        ltline = len(tline)

        if ltline == 0:
            self._type = 'blank'
        else:
            if commentchar and self._line.startswith(commentchar):
                self._type = 'comment'
            else:
                self._type = 'data'
                self._token = tline.split()
                self._ntok = len(self.token)

    @property
    def line(self):
        """Line read accessor

            Returns:
                (str): Original input line stripped of line terminators
        """
        return self._line

#    @line.setter
#    def line(self, line):
#        """Line write accessor
#            Raises:
#                ValueError
#        """
#        return self._line

    @property
    def type(self):
        """Type read accessor

            Returns:
                (str): Type of input line; one of 'blank', 'comment', or 'data'
        """
        return self._type

    @property
    def typecode(self):
        """Type code read accessor

            Returns:
                (str): Type code of input line; one of 'B', 'C', or 'D',
                  corresponding to 'blank', 'comment', or 'data', respectively
        """
        return self._type[0].upper()

    @property
    def ntok(self):
        """Token count read accessor

            Returns:
                (int): Number of tokens found (0 except for 'data' lines)
        """
        return self._ntok

    @property
    def token(self):
        """Token list read accessor

            Returns:
                ([str]): Tokens parsed from line ('data' lines only, otherwise
                  empty)
        """
        return self._token

    def as_log(self, logfmt='{0:-6d} [{1:s}] {2:s}'):
        """Return line number (ipos), type code, and original line to assist
        in finding input errors. Type code is 'B', 'C', or 'D', corresponding
        to 'blank', 'comment', and 'data', respectively

            Args:
                logfmt (str): Format string for producing log output. Field 0
                  is the `ipos` attribute, field 1 is the type code, and field
                  2 is the `line` attribute

            Returns:
                (str): Formatted line with metadata
        """
        return logfmt.format(self.ipos, self.typecode, self._line)
