""" Piping classes """

from __future__ import absolute_import, division, print_function
# from future.utils import iteritems
# from pprint import pprint

# from collections import OrderedDict, defaultdict
# import logging
# import pint
import scipy.constants as sc
import fluids.vectorized as fv

from tabulate import tabulate


class Pipe(object):
    """Simple control volume class """
    def __init__(self, label, length, idiameter, eroughness,
                 zbottom, ztop, schedule='', twall=0.0, odiameter=0.0):
        """Simple pipe segment constructor

          All parameters are mandatory except schedule, twall, and odiameter.
          Only one of the three are mandatory to set wall thickness.

          :param label: Text description, string, Mandatory.
          :param length: Length in meters, real, Mandatory.
          :param idiameter: Pipe inner diameter in meters, real, Mandatory.
          :param eroughness: Pipe relative roughness, real, Mandatory.
          :param zbottom: Lower elevation (in meters) with respect to site
                          zero, real, Mandatory.
          :param ztop: Upper elevation (in meters) with respect to site zero,
                       real, Mandatory.
          :param schedule: Pipe schedule, string, Optional.
          :param twall: Pipe wall thickness in meters, real. Optional.
          :param odiameter: Outer pipe diameter in meters, real. Optional.
        """

        # Set atmospheric composition (pressure, T, mole fractions)
        self.label = label
        self.length = length
        self.idiameter = idiameter
        self.eroughness = eroughness
        self._zbottom = zbottom
        self._ztop = ztop

        if odiameter:
            self._odiameter = odiameter
            self._twall = (self._odiameter - self._idiameter) / 2.0
        elif twall:
            self._twall = twall
            self._odiameter = self.idiameter + 2.0 * twall
        elif schedule:
            self.nearest_dimensions_from_schedule(schedule)
# TODO Log that pipe ID has been adjusted beyond 1% of original ID
#                if abs(Di - self.idiameter) / self.idiameter > 0.01:
#                    log("")
#            except ValueError as err:
# TODO Handle nearest_pipe ValueError

    @property
    def idiameter(self):
        """ Inner diameter accessor - read """
        return self._idiameter

    @idiameter.setter
    def idiameter(self, idiameter):
        """ Inner diameter accessor - write """
        if idiameter < 1.0E-3:
            raise ValueError('Inner diameter too small (<1mm)')
        elif idiameter > 10.0:
            raise ValueError('Inner diameter too large (>10m)')

        self._idiameter = idiameter
        self._flow_area = sc.pi * self._idiameter**2 / 4.0

    @property
    def length(self):
        """ Inner diameter accessor - read """
        return self._length

    @length.setter
    def length(self, length):
        """ Inner diameter accessor - write """
        if length < 1.0E-3:
            raise ValueError('Length too small (<1mm)')
        elif length > 1000.0:
            raise ValueError('Length too large (>1000m)')

        self._length = length

    @property
    def eroughness(self):
        """ Pipe wall roughness accessor - read """
        return self._eroughness

    @eroughness.setter
    def eroughness(self, eroughness):
        """ Pipe wall roughness accessor - write """
        self._eroughness = eroughness

        if self._eroughness < 0.0:
            raise ValueError('Eroughness less than 0.0')
        elif self._eroughness > 0.1:
            raise ValueError('Eroughness larger than 0.1')

    @property
    def flow_area(self):
        """ Pipe flow area - read """
        return self._flow_area

    @flow_area.setter
    def flow_area(self, flow_area):
        """ Pipe flow area - write """
        raise ValueError('Cannot directly set flow area; '
                         'value is derived from inner diameter')

    def nearest_material_roughness(self, surface, clean):
        """Find nearest surface roughness by surface finish name
        and cleanliness"""
        surface_key = fv.nearest_material_roughness(surface, clean)
        if surface_key:
            eroughness = fv.material_roughness(surface_key)
#            if clean:
#                clean_tag = 'clean'
#            else:
#                clean_tag = 'fouled'
# TODO Log msg so assumption is clear
#            msg = "Note: Surface roughness set to {0:0.4E}; used '{1:s}' " \
#                  "based on specification '{2:s}, {3:s}'".format(
#                  eroughness, surface_key, surface, clean_tag)
            self.eroughness = eroughness
        else:
            raise ValueError('Cannot find surface finish corresponding '
                             'to "{0:s}"'.format(surface))

    def nearest_dimensions_from_schedule(self, schedule, dnominal=0.0):
        """Find dimensions closest to inner diameter or nominal diameter
        and given pipe schedule"""
        if dnominal > 0.0:
            try:
                (NPS, Di, Do, t) = fv.nearest_pipe(NPS=dnominal,
                                                   schedule=schedule)
                self._twall = t
                self._odiameter = Do
                self.idiameter = Di
            except:
                raise ValueError('Cannot find dimensions corresponding '
                                 'to {0:0.4f}-in nominal diameter and pipe '
                                 'schedule "{1:s}"'
                                 .format(dnominal, schedule))

        else:
            try:
                (NPS, Di, Do, t) = fv.nearest_pipe(Di=self.idiameter,
                                                   schedule=schedule)
                self._twall = t
                self._odiameter = Do
                self.idiameter = Di
            except:
                raise ValueError('Cannot find dimensions corresponding '
                                 'to {0:0.4e} m inner diameter and pipe '
                                 'schedule "{1:s}"'
                                 .format(self.idiameter, schedule))

#        raise ValueError('FAKE: NPS is {0:0.4f}'.format(NPS))

#    @property
#    def height(self):
#        """ Height accessor - read """
#        return self._height
#
#    @height.setter
#    def height(self, height):
#        """ Height accessor - write """
#        self._height = height
#        if hasattr(self, '_volume') and self._volume > 0.0:
#            self._cxarea = self._volume / self._height
#        elif hasattr(self, '_cxarea') and self._cxarea > 0.0:
#            self._volume = self._height * self._cxarea
#
#        if abs(self._volume - self._cxarea * self._height) > 1.0E-9:
#            raise ValueError('Setting height results in inconsistent volume')
#
#    @property
#    def pressure(self):
#        """ Pressure accessor - read """
#        return self._pressure
#
#    @pressure.setter
#    def pressure(self, pressure):
#        """ Pressure accessor - write """
#        self._pressure = pressure
#
#        if self._pressure <= 1.0E-9:
#            raise ValueError('Pressure too low')
#
#    @property
#    def temperature(self):
#        """ Temperature accessor - read """
#        return self._temperature
#
#    @temperature.setter
#    def temperature(self, temperature):
#        """ Temperature accessor - write """
#        self._temperature = temperature
#
#        if self._temperature <= 250.0:
#            raise ValueError('Temperature too low')
#
#    @property
#    def elevation(self):
#        """ Elevation accessor - read """
#        return self._elevation
#
#    @elevation.setter
#    def elevation(self, elevation):
#        """ Elevation accessor - write """
#        if abs(elevation) > 1.0E3:
#            raise ValueError('Elevation too high or too low: '
#                             '{0:0.4E} m'.format(elevation))
#        else:
#            self._elevation = elevation
#
#    def volume_from_floor(self, cxarea, height):
#        """ Set volume from rectangular dimensions """
#        if cxarea <= 0.001:
#            raise ValueError('Cross-sectional area must be '
#                             'greater than 0.001 m2')
#
#        if height <= 0.001:
#            raise ValueError('Height must be greater than 1mm')
#
#        self._height = height
#        self._cxarea = cxarea
#        self._volume = height * cxarea
#
#    def volume_as_box(self, length, width, height):
#        """ Set volume from rectangular dimensions """
#        if length <= 0.001:
#            raise ValueError('Length must be greater than 1mm')
#
#        if width <= 0.001:
#            raise ValueError('Width must be greater than 1mm')
#
#        if height <= 0.001:
#            raise ValueError('Height must be greater than 1mm')
#
#        self.volume_from_floor(cxarea=length*width, height=height)
#
#    def volume_as_cylinder(self, height, radius=None, diameter=None):
#        """ Set volume from cylindrical dimensions """
#        if radius and not diameter:
#            if radius > 0.001:
#                cxarea = sc.pi * radius**2
#            else:
#                raise ValueError('Radius must be greater than 1mm')
#        elif diameter and not radius:
#            if diameter > 0.002:
#                cxarea = sc.pi * diameter**2 / 4.0
#            else:
#                raise ValueError('Diameter must be greater than 2mm')
#        else:
#            raise ValueError('Cylindrical geometry must have only '
#                             'one of radius and diameter specified')
#
#        self.volume_from_floor(cxarea=cxarea, height=height)
#
#    def is_complete(self):
#        """ Perform basic sanity checking on the Volume object """
#        errfields = []
#        if self.pressure <= 1.0E-9 or self.pressure > 1.0E9:
#            errfields.append('pressure')
#
#        if self.temperature <= 250.0 or self.temperature > 5000.0:
#            errfields.append('temperature')
#
#        if self.height <= 1.0E-6 or self.height > 5000.0:
#            errfields.append('height')
#
#        if self.cxarea <= 1.0E-12 or self.cxarea > 2.5E7:
#            errfields.append('cxarea')
#
#        if self.volume <= 1.0E-6 or self.volume > 1.25E11:
#            errfields.append('volume')
#        else:
#            voldev = abs(self.volume - self.cxarea * self.height)
#            pcterr = 100.0 * voldev / self.volume
#            if pcterr > 0.001:
#                raise ValueError('Volume inconsistent with height and area,'
#                                 ' deviation of {0:0.3f}%'.format(pcterr))
#
#        if len(errfields) < 1:
#            errfields = None
#
#        return errfields
#
#    def v_at_z(self, z):
#        """ Partial volume from bottom to height z """
#        vz = 0.0
#        fz = z / self.height
#        if fz > 1.0:
#            if z - self.height < 0.001:
#                vz = self.volume
#            else:
#                raise ValueError('Elevation above volume height'
#                                 ' - {0:0.3f}% deviation'
#                                 .format(100.0 * (fz - 1.0)))
#        elif fz > 0.0:
#            vz = fz * self.volume
#        else:
#            if self.height - z < 0.001:
#                vz = 0.0
#            else:
#                raise ValueError('Elevation below volume base'
#                                 ' - {0:0.3f}% deviation'
#                                 .format(100.0 * fz))
#        return vz
#
#    def cxarea_at_z(self, z):
#        """ Cross-sectional area at a height z. Returns cxarea for
#        reasonable values of z """
#        fz = z / self.height
#        if fz > 1.0:
#            if z - self.height >= 0.001:
#                raise ValueError('Elevation above volume height'
#                                 ' - {0:0.3f}% deviation'
#                                 .format(100.0 * (fz - 1.0)))
#        elif fz < 0.0:
#            if self.height - z >= 0.001:
#                raise ValueError('Elevation below volume base'
#                                 ' - {0:0.3f}% deviation'
#                                 .format(100.0 * fz))
#        return self.cxarea
#
    def as_table(self, tablefmt="psql",
                 headers=['length', 'idiameter', 'odiameter',
                          'twall', 'schedule', 'eroughness',
                          'ztop', 'zbottom']):
        """ Display object as table """
        tbl = "Pipe: {0:s}\n".format(self.label)
        cd_tuples = [(hdr, getattr(self, hdr)) for hdr in headers]

        tbl += tabulate(cd_tuples, tablefmt=tablefmt)

        return tbl
