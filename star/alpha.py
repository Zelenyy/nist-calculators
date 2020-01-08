from typing import Union

import numpy as np
import tables

from star._data_converter import NIST_STAR_HDF5_PATH, ProtonNames, Names
from star.alpha_materials import AlphaMaterials
from xcom.interpolators import make_log_log_spline


class AlphaSTARCalculator:

    def __init__(self, material: AlphaMaterials):
        """

        """
        with tables.open_file(NIST_STAR_HDF5_PATH) as h5file:
            group = tables.Group(h5file.root, "helium_ions")
            data = h5file.get_node(group, material.name).read()
            energy = h5file.get_node(group, "energy").read()
            self.default_energy = energy
            self.splines = {}
            for name in data.dtype.names:
                self.splines[name] = make_log_log_spline(energy, data[name])

    def calculate_csda_ranges(self, energy: np.ndarray = None) -> np.ndarray:
        """
        CSDA ranges, in `g/cm2`
        Parameters
        ----------
        energy
                energies of helium ions in `MeV`, used `self.default_energy` by default
        Returns
        -------
        data : ndarray with  CSDA ranges, in `g/cm2`
        """
        if energy is None:
            energy = self.default_energy
        energy = np.asarray(energy)
        return self.splines[Names.CSDA_RANGE](energy)

    def calculate_detour_factors(self, energy: np.ndarray = None) -> np.ndarray:
        """
        Detour factors
        Parameters
        ----------
        energy
                energies of helium ions in `MeV`, used `self.default_energy` by default
        Returns
        -------
        data : ndarray with  detour factors
        """
        if energy is None:
            energy = self.default_energy
        energy = np.asarray(energy)
        return self.splines[ProtonNames.DETOUR_FACTOR](energy)

    def calculate_projected_range(self, energy: np.ndarray = None) -> np.ndarray:
        """
        Projected range (average penetration depths), `g/cm2`
        Parameters
        ----------
        energy
                energies of helium ions in `MeV`, used `self.default_energy` by default
        Returns
        -------
        data : ndarray with projected range, `g/cm2`
        """
        if energy is None:
            energy = self.default_energy
        energy = np.asarray(energy)
        return self.calculate_detour_factors(energy) * self.calculate_csda_ranges(energy)

    def calculate_electronic_stopping_powers(self, energy: Union[list, np.ndarray, float] = None) -> np.ndarray:
        """
        Electronic stopping powers, in `MeV cm2/g`
        Parameters
        ----------
        energy
                energies of helium ions in `MeV`, used `self.default_energy` by default
        Returns
        -------
        data : ndarray with  electronic stopping powers, in `MeV cm2/g`
        """
        if energy is None:
            energy = self.default_energy
        energy = np.asarray(energy)
        return self.splines[ProtonNames.ELECTRONIC_STOPPING_POWER](energy)

    def calculate_nuclear_stopping_powers(self, energy: Union[list, np.ndarray, float] = None) -> np.ndarray:
        """
        Nuclear stopping powers, in `MeV cm2/g`
        Parameters
        ----------
        energy
                energies of helium ions in `MeV`, used `self.default_energy` by default
        Returns
        -------
        data : ndarray with  nuclear stopping powers, in `MeV cm2/g`
        """
        if energy is None:
            energy = self.default_energy
        energy = np.asarray(energy)
        return self.splines[ProtonNames.NUCLEAR_STOPPING_POWER](energy)

    def calculate_total_stopping_powers(self, energy: Union[list, np.ndarray, float] = None) -> np.ndarray:
        """
        Total stopping powers, in `MeV cm2/g`
        Parameters
        ----------
        energy
                energies of helium ions in `MeV`, used `self.default_energy` by default
        Returns
        -------
        data : ndarray with  total stopping powers, in `MeV cm2/g`
        """
        if energy is None:
            energy = self.default_energy
        energy = np.asarray(energy)
        return self.calculate_electronic_stopping_powers(energy) + self.calculate_nuclear_stopping_powers(energy)

    def calculate_table(self, energy: Union[list, np.ndarray, float] = None):
        """
        Summary table:
            * Proton kinetic energy, `MeV`
            * Electronic stopping power, `MeV cm2/g`
            * Nuclear stopping power, `MeV cm2/g`
            * Total stopping power, `MeV cm2/g`
            * CSDA range, `g/cm2`
            * Projected range, `g/cm2`
            * Detour factor

                Parameters
        ----------
        energy
                energies of helium ions in `MeV`, used `self.default_energy` by default
        Returns
        -------
        data : record ndarray with  all data
        """
        if energy is None:
            energy = self.default_energy
        energy = np.asarray(energy)
        dtype = np.dtype(
            [
                (Names.ENERGY, "d"),
                (ProtonNames.ELECTRONIC_STOPPING_POWER, "d"),
                (ProtonNames.NUCLEAR_STOPPING_POWER, "d"),
                (ProtonNames.TOTAL_STOPPING_POWER, "d"),
                (Names.CSDA_RANGE, "d"),
                (ProtonNames.PROJECTED_RANGE, "d"),
                (ProtonNames.DETOUR_FACTOR, "d")
            ]
        )
        n = len(energy)
        data = np.zeros(n, dtype=dtype)
        data[Names.ENERGY] = energy
        data[ProtonNames.ELECTRONIC_STOPPING_POWER] = self.calculate_electronic_stopping_powers(energy)
        data[ProtonNames.NUCLEAR_STOPPING_POWER] = self.calculate_nuclear_stopping_powers(energy)
        data[ProtonNames.TOTAL_STOPPING_POWER] = data[ProtonNames.ELECTRONIC_STOPPING_POWER] + data[
            ProtonNames.NUCLEAR_STOPPING_POWER]
        data[Names.CSDA_RANGE] = self.calculate_csda_ranges(energy)
        data[ProtonNames.DETOUR_FACTOR] = self.calculate_detour_factors(energy)
        data[ProtonNames.PROJECTED_RANGE] = data[ProtonNames.DETOUR_FACTOR] * data[Names.CSDA_RANGE]
        return data
