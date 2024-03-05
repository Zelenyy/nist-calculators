"""
XCOM attenuation coefficients

Based on NIST XCOM data: https://www.nist.gov/pml/xcom-photon-cross-sections-database
"""

import warnings
from typing import Union

import numpy as np

from ._data_converter import NameProcess
from .interpolators import Interpolators, Material, MaterialFactory

_AVOGADRO = 0.60221367  # * 1e+24 mole^-1

# fmt: off
ENERGY_GRID_DEFAULT = np.array([1.0E+03, 1.5E+03, 2.0E+03, 3.0E+03, 4.0E+03, 5.0E+03,
                                6.0E+03, 8.0E+03, 1.0E+04, 1.5E+04, 2.0E+04, 3.0E+04, 4.0E+04,
                                5.0E+04, 6.0E+04, 8.0E+04, 1.0E+05, 1.5E+05, 2.0E+05, 3.0E+05,
                                4.0E+05, 5.0E+05, 6.0E+05, 8.0E+05, 1.0E+06, 1.022E+06, 1.25E+06,
                                1.5E+06, 2.0E+06, 2.044E+06, 3.0E+06, 4.0E+06, 5.0E+06, 6.0E+06,
                                7.0E+06, 8.0E+06, 9.0E+06, 1.0E+07, 1.1E+07, 1.2E+07, 1.3E+07,
                                1.4E+07, 1.5E+07, 1.6E+07, 1.8E+07, 2.0E+07, 2.2E+07, 2.4E+07,
                                2.6E+07, 2.8E+07, 3.0E+07, 4.0E+07, 5.0E+07, 6.0E+07, 8.0E+07,
                                1.0E+08, 1.5E+08, 2.0E+08, 3.0E+08, 4.0E+08, 5.0E+08, 6.0E+08,
                                8.0E+08, 1.0E+09, 1.5E+09, 2.0E+09, 3.0E+09, 4.0E+09, 5.0E+09,
                                6.0E+09, 8.0E+09, 1.0E+10, 1.5E+10, 2.0E+10, 3.0E+10, 4.0E+10,
                                5.0E+10, 6.0E+10, 8.0E+10, 1.0E+11], dtype='d')
# fmt: on

_INTERPOLATORS = Interpolators()


def calculate_attenuation(material: Material, energy: np.ndarray = None) -> np.ndarray:
    """
    Calculate attenuation (cm2/gram) for gamma-ray (at energies between 1 keV
    and 100 GeV) for next process:

        * Coherent scattering
        * Incoherent (Compton) scattering
        * Photoelectric absorption
        * Pair production in the field of the atomic nucleus and in the field of
          the atomic electrons

    Based on NIST XCOM data: https://www.nist.gov/pml/xcom-photon-cross-sections-database

    Parameters
    ----------
    material
            special class description simple material or compound

    energy
            energies of gamma-quanta in eV, used `ENERGY_GRID_DEFAULT` by default

    Returns
    -------
    data : ndarray with attenuation in cm2/gram
    """
    if not isinstance(material, Material):
        raise ValueError("Except material")

    if len(material) == 1:
        # Single constituent material
        element = material.elements_by_Z[0]
        data = calculate_cross_section(element, energy)
        # Attenuation coefficient = macro_cross_section/density = \
        # = micro_cross_section/atom_weight[gr]
        # atom_weight[gr] = atom_weight[amu] / AVOGADRO
        atom_weight = MaterialFactory.get_element_mass(element)
        for name in data.dtype.names:
            if name != "energy":
                data[name] *= _AVOGADRO / atom_weight
        return data
    elif len(material) > 1:
        # Mixture of materials
        atom_weights_amu = MaterialFactory.get_elements_mass_list(
            material.elements_by_Z
        )
        # Initialize data with data for first element
        data = calculate_cross_section(material.elements_by_Z[0], energy)
        for name in data.dtype.names:
            if name != "energy":
                data[name] *= material.weights[0] * _AVOGADRO / atom_weights_amu[0]
        # Add in other elements
        for atom_weight_amu, element, weight in zip(
            atom_weights_amu[1:], material.elements_by_Z[1:], material.weights[1:]
        ):
            temp = calculate_cross_section(element, energy)
            for name in data.dtype.names:
                if name != "energy":
                    data[name] += temp[name] * weight * _AVOGADRO / atom_weight_amu
        return data
    else:
        raise ValueError("Empty material")


def calculate_cross_section(
    element: Union[int, str], energy: np.ndarray = None
) -> np.ndarray:
    """
    Calculate cross-section (barn/atom) for gamma-ray (at energies between 1 keV
    and 100 GeV) for next process:

        * Coherent scattering
        * Incoherent (Compton) scattering
        * Photoelectric absorption
        * Pair production in the field of the atomic nucleus and in the field of
          the atomic electrons

    Based on NIST XCOM data: https://www.nist.gov/pml/xcom-photon-cross-sections-database

    Parameters
    ----------
    element
            atomic number or symbol of element

    energy
            energies of gamma-quanta in eV, used `ENERGY_GRID_DEFAULT` by default

    Returns
    -------
    data : ndarray with cross-section in barn/atom
    """
    if energy is None:
        energy = ENERGY_GRID_DEFAULT

    if (
        np.max(energy) > ENERGY_GRID_DEFAULT[-1]
        or np.min(energy) < ENERGY_GRID_DEFAULT[0]
    ):
        warnings.warn(
            "Energy or energies requested is outside of tabulated data, "
            + "using linear extrapolation",
            UserWarning,
        )

    if not isinstance(element, int):
        element = MaterialFactory.get_element_from_symbol(element)

    n = len(energy)
    dtype = np.dtype(
        [
            ("energy", "d"),
            (NameProcess.COHERENT, "d"),
            (NameProcess.INCOHERENT, "d"),
            (NameProcess.PHOTOELECTRIC, "d"),
            (NameProcess.PAIR_ATOM, "d"),
            (NameProcess.PAIR_ELECTRON, "d"),
            ("total_without_coherent", "d"),
            ("total", "d"),
        ]
    )

    data = np.zeros(n, dtype=dtype)
    data["energy"] = np.asarray(energy)

    for k, v in _INTERPOLATORS.get_interpolators(element).items():
        data[k] = v(data["energy"])
        data["total"] += data[k]
    data["total_without_coherent"] = data["total"] - data[NameProcess.COHERENT]
    return data
