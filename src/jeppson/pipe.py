""" Piping classes """

from __future__ import absolute_import, division, print_function
# from future.utils import iteritems
# from pprint import pprint

# from collections import OrderedDict, defaultdict
import logging
# import pint

import fluids.vectorized as fv
import scipy.constants as sc
from tabulate import tabulate

LOG = logging.getLogger(__name__)


class Pipe(object):
    """Simple pipe segment class

    All parameters are optional except label and length.
    A sensible combination of NPS and schedule or diameter, wall
    thickness, and pipe schedule are needed to define radial pipe
    dimensions. Similarly, relative pipe roughness may be set directly or
    will be derived from surface and cleanliness or absolute roughness and
    inner diameter. Relative roughness will be taken as zero (smooth,
    clean pipe) unless otherwise specified.
    
    If insufficient parameters are specified in the constructor to completely
    define radial geometry, the object will be left in an inconsistent state
    which may cause the `as_table()` method to throw AttributeException errors.
    This may be resolved by defining additional radial geometry attributes (any
    two of `idiameter`, `odiameter`, and `twall`) or calling the
    `nearest_dimensions_from_schedule()` method.

    Attributes:
        label (str): text description
        length (float): length in meters
        idiameter (float): pipe inner diameter in meters
        odiameter (float): outer pipe diameter in meters
        twall (float): pipe wall thickness in meters
        nps (float): nominal pipe diameter in inches; empty if not scheduled
            pipe
        schedule (string): pipe schedule; empty if not scheduled pipe
        surface (str): inner pipe surface material. Read only.
        clean (str): pipe cleanliness. Read only.
        eroughness (float): pipe relative roughness

    Args:
        label (str): text description, Mandatory.
        length (float): length in meters. Mandatory.
        idiameter (float): pipe inner diameter in meters. Optional.
        odiameter (float): outer pipe diameter in meters. Optional.
        twall (float): pipe wall thickness in meters. Optional.
        nps (float): nominal pipe diameter, inches. Optional.
        schedule (string): pipe schedule. Optional.
        eroughness (float): relative pipe roughness. Optional.
        froughness (float): absolute pipe roughness, in meters. Optional.
        surface (str): inner pipe surface material. Optional.
        is_clean (bool): pipe cleanliness. Optional.
    """
# Not used yet; needed for facility modeling
#        zbottom (float): lower elevation (in meters) with respect to site
#                 zero, mandatory.
#        ztop (float): upper elevation (in meters) with respect to site zero,
#              mandatory.

    def __init__(self, label, length, idiameter=0.0, odiameter=0.0, twall=0.0,
                 nps=0.0, schedule='', eroughness=0.0, froughness=0.0,
                 surface='smooth', is_clean=True):
        """Simple pipe segment constructor

        All parameters are optional except label and length.
        A sensible combination of NPS and schedule or diameter, wall
        thickness, and pipe schedule are needed to define radial pipe
        dimensions. Similarly, relative pipe roughness may be set directly or
        will be derived from surface and cleanliness or absolute roughness and
        inner diameter. Relative roughness will be taken as zero (smooth,
        clean pipe) unless otherwise specified.

        Args:
            label (str): text description, Mandatory.
            length (float): length in meters. Mandatory.
            idiameter (float): pipe inner diameter in meters. Optional.
            odiameter (float): outer pipe diameter in meters. Optional.
            twall (float): pipe wall thickness in meters. Optional.
            nps (float): nominal pipe diameter, inches. Optional.
            schedule (string): pipe schedule. Optional.
            eroughness (float): relative pipe roughness. Optional.
            froughness (float): absolute pipe roughness, in meters. Optional.
            surface (str): inner pipe surface material. Optional.
            is_clean (bool): pipe cleanliness. Optional.

        Raises:
            ValueError: An error occurred setting an attribute.
        """

#            zbottom (float): lower elevation (in meters) with respect to site
#                     zero, mandatory.
#            ztop (float): upper elevation (in meters) with respect to site
#                zero, mandatory.

        self.label = label
        self.length = length
        self.schedule = schedule
        self.nps = nps
#        self._zbottom = zbottom
#        self._ztop = ztop

        # Set radial dimensions
        if nps > 0.0 and schedule:
            self.nearest_dimensions_from_schedule(schedule, nps)
        elif idiameter:
            self._idiameter = idiameter
            self._update_flow_area()
            if odiameter:
                self._odiameter = odiameter
                self._twall = (self._odiameter - self._idiameter) / 2.0
            elif twall:
                self._twall = twall
                self._odiameter = self.idiameter + 2.0 * twall
            elif schedule:
                self.nearest_dimensions_from_schedule(schedule)
                if (2.0 * abs(idiameter - self._idiameter)
                        / (idiameter + self._idiameter)) > 0.01:
                    LOG.info('Inside diameter of pipe "{0:s}" has '
                             'changed more than 1% from {1:0.4E} m to '
                             '{1:0.4E} m to match schedule {3:s}'
                             .format(label, idiameter, self._idiameter,
                                     schedule))

        # Set relative roughness
        if eroughness > 0.0:
            self.eroughness = eroughness
            self._surface = 'unknown'
            self._clean = 'unknown'
        elif froughness > 0.0:
            self.eroughness = froughness / self.idiameter
            self._surface = 'unknown'
            self._clean = 'unknown'
        else:
            if surface.strip().lower() == 'smooth':
                self._surface = 'smooth'
                self._clean = 'clean'
                self.eroughness = 0.0
            else:
                self._surface = surface.strip()
                self.nearest_material_roughness(surface, is_clean)

    def _update_flow_area(self):
        """Internal method to update flow area on a change to inner diameter"""
        self._flow_area = sc.pi * self._idiameter**2 / 4.0

    @property
    def idiameter(self):
        """Inner diameter read accessor

            Returns:
                float: inner diameter, meters"""
        return self._idiameter

    @idiameter.setter
    def idiameter(self, idiameter):
        """Inner diameter accessor - write

        Raises:
            ValueError: Unreasonable value for inner diameter."""
        if idiameter < 1.0E-3:
            raise ValueError('Inner diameter too small (<1mm)')
        elif idiameter > 10.0:
            raise ValueError('Inner diameter too large (>10m)')

        self._idiameter = idiameter
        self._update_flow_area()

        if hasattr(self, '_odiameter') and self._odiameter > 0.0:
            twall = (self._odiameter - idiameter) / 2.0
            if twall <= 0.0:
                raise ValueError('Inner diameter larger than outer diameter')
            elif twall < 1.0E-4:
                raise ValueError('Pipe wall too thin (<0.1mm)')
            else:
                self._twall = twall
        elif hasattr(self, '_twall') and self._twall > 0.0:
            self._odiameter = self._idiameter + 2.0 * self._twall

    @property
    def odiameter(self):
        """Outer diameter read accessor

            Returns:
                float: outer diameter, meters"""
        return self._odiameter

    @odiameter.setter
    def odiameter(self, odiameter):
        """Outer diameter accessor - write

        Raises:
            ValueError: Unreasonable value for inner diameter."""
        if odiameter < 1.0E-3:
            raise ValueError('Outer diameter too small (<1mm)')
        elif odiameter > 10.0:
            raise ValueError('Outer diameter too large (>10m)')

        self._odiameter = odiameter

        if hasattr(self, '_idiameter') and self._idiameter > 0.0:
            twall = (odiameter - self._idiameter) / 2.0
            if twall <= 0.0:
                raise ValueError('Outer diameter smaller than inner diameter')
            elif twall < 1.0E-4:
                raise ValueError('Pipe wall too thin (<0.1mm)')
            else:
                self._twall = twall
        elif hasattr(self, '_twall') and self._twall > 0.0:
            self._idiameter = self._odiameter - 2.0 * self._twall
            self._update_flow_area()

    @property
    def twall(self):
        """Wall thickness read accessor

            Returns:
                float: Wall thickness, meters"""
        return self._twall

    @twall.setter
    def twall(self, twall):
        """Wall thickness accessor - write

        Raises:
            ValueError: Unreasonable value for wall thickness."""
        if twall < 1.0E-4:
            raise ValueError('Wall thickness too small (<0.1mm)')
        elif twall > 1.0:
            raise ValueError('Wall thickness too large (>1m)')

        self._twall = twall
        if hasattr(self, '_idiameter') and self._idiameter > 0.0:
            self._odiameter = self._idiameter + 2.0 * self._twall
        elif hasattr(self, '_odiameter') and self._odiameter > 0.0:
            if 2.0 * self._twall - self._idiameter < 1.0E-4:
                raise ValueError('Outer diameter smaller than twice '
                                 'the wall thickness')
            else:
                self._idiameter = self._odiameter - 2.0 * self._twall
                self._update_flow_area()

    @property
    def length(self):
        """Pipe length read accessor

            Returns:
                float: pipe length, meters"""
        return self._length

    @length.setter
    def length(self, length):
        """Inner diameter accessor - write

        Raises:
            ValueError: Unreasonable value for pipe length."""
        if length < 1.0E-3:
            raise ValueError('Length too small (<1mm)')
        elif length > 1000.0:
            raise ValueError('Length too large (>1000m)')

        self._length = length

    @property
    def eroughness(self):
        """Pipe wall relative roughness read accessor

            Returns:
                float: pipe interior wall relative roughness, dimensionless"""
        return self._eroughness

    @eroughness.setter
    def eroughness(self, eroughness):
        """Pipe wall roughness accessor - write

        Raises:
            ValueError: Unreasonable value for relative roughness."""
        self._eroughness = eroughness

        if self._eroughness < 0.0:
            raise ValueError('Eroughness less than 0.0')
        elif self._eroughness > 0.1:
            raise ValueError('Eroughness larger than 0.1')

    @property
    def flow_area(self):
        """Pipe flow area read accessor

            Returns:
                float: pipe interior cross-sectional flow area, square
                    meters"""
        return self._flow_area

    @flow_area.setter
    def flow_area(self, flow_area):
        """Pipe flow area - write

        Raises:
            ValueError: Cannot set derived quantity. """
        raise ValueError('Cannot directly set flow area; '
                         'value is derived from inner diameter')

    @property
    def clean(self):
        """Pipe cleanliness read accessor

            Returns:
                str: Pipe cleanliness label ('clean', 'fouled', or 'unknown'
        """
        return self._clean

    @clean.setter
    def clean(self, clean):
        """Pipe cleanliness - write

        Raises:
            ValueError: Cannot set derived quantity. """
        raise ValueError('Cannot directly set pipe cleanliness; use '
                         'constructor or nearest_material_roughness() method')

    @property
    def surface(self):
        """Pipe interior surface material read accessor

            Returns:
                str: Pipe interior surface material"""
        return self._surface

    @surface.setter
    def surface(self, surface):
        """Pipe surface - write

        Raises:
            ValueError: Cannot set derived quantity. """
        raise ValueError('Cannot directly set pipe surface; use '
                         'constructor or nearest_material_roughness() method')

    def nearest_material_roughness(self, surface, is_clean):
        """Find nearest surface roughness by surface finish name
        and cleanliness

        Args:
            surface (str): text description of pipe surface
            is_clean (bool): pipe cleanliness

        Raises:
            ValueError: Failed to find surface roughness for material
                description."""

        surface_key = fv.nearest_material_roughness(surface, is_clean)

        lc_surface = str(surface).lower().strip()
        lc_surface_key = str(surface_key).lower().strip()
        if lc_surface not in lc_surface_key:
            self._surface = 'unknown'
            self._clean = 'unknown'
            raise ValueError('Surface specified "{0:s}" too different '
                             'from surface found "{1:s}"'
                             .format(surface, surface_key))

        eroughness = fv.material_roughness(surface_key)
        self.eroughness = eroughness

        self._surface = surface.strip()
        if is_clean:
            self._clean = 'clean'
        else:
            self._clean = 'fouled'

        msg = 'Note: Surface roughness of pipe "{0:s}" set to {1:0.4E}; ' \
              'used "{2:s}" based on specification "{3:s}, {4:s}"' \
              .format(self.label, self.eroughness, surface_key, self.surface,
                      self._clean)
        LOG.info(msg)

    def nearest_dimensions_from_schedule(self, schedule, dnominal=0.0):
        """Find dimensions closest to inner diameter or nominal diameter
        and given pipe schedule

        Args:
            schedule (str): Pipe schedule
            dnominal (float): Nominal pipe diameter, inches.

        Raises:
            ValueError: Failed to find dimensions for given size and pipe
                schedule."""
        if dnominal > 0.0:
            try:
                (NPS, Di, Do, t) = fv.nearest_pipe(NPS=dnominal,
                                                   schedule=schedule)
            except Exception:
                raise ValueError('Cannot find dimensions corresponding '
                                 'to {0:0.4f}-in nominal diameter and pipe '
                                 'schedule "{1:s}"'
                                 .format(dnominal, schedule))

        else:
            try:
                (NPS, Di, Do, t) = fv.nearest_pipe(Di=self.idiameter,
                                                   schedule=schedule)
            except ValueError:
                raise ValueError('Cannot find dimensions corresponding '
                                 'to {0:0.4e} m inner diameter and pipe '
                                 'schedule "{1:s}"'
                                 .format(self.idiameter, schedule))

        self.schedule = schedule
        self.nps = NPS
        self._twall = t
        self._odiameter = Do
        self._idiameter = Di
        self._update_flow_area()

    def as_table(self, tablefmt="psql",
                 headers=['length', 'idiameter', 'odiameter', 'twall',
                          'eroughness'],
                 extheaders=['nps', 'schedule', 'surface', 'clean']):
        """Display object as table

            Args:
                tablefmt (str): Table format; see `tabulate` documentation
                headers ([str]): Mandatory list of attributes to list in table
                extheaders ([str]): List of optional attributes to list if
                    they exist.

            Returns:
                str: Text table of pipe attributes."""
        tbl = "Pipe: {0:s}\n".format(self.label)
        cd_tuples = [(hdr, getattr(self, hdr)) for hdr in headers]
        for extkey in extheaders:
            if hasattr(self, extkey):
                # Filter empty 'noise' attributes from display
                if extkey == 'nps' and self.nps > 0.0:
                    cd_tuples.append(('nps', self.nps))
                elif extkey == 'schedule' and self.schedule != '':
                    cd_tuples.append(('schedule', self.schedule))
                else:
                    cd_tuples.append((extkey, getattr(self, extkey)))

        tbl += tabulate(cd_tuples, tablefmt=tablefmt)

        return tbl
