"""Interpolators for tabulated data"""
import csv
import os
from typing import Callable, List, Union, Optional
import warnings

import numpy as np
import tables
from scipy.interpolate import CubicSpline, interp1d

from ._data_converter import NameProcess

ROOT_PATH = os.path.dirname(__file__)
DATA_PATH = os.path.join(ROOT_PATH, 'data')
NIST_XCOM_HDF5_PATH = os.path.join(DATA_PATH, 'NIST_XCOM.hdf5')
PERIODIC_TABLE_PATH = os.path.join(DATA_PATH, "PeriodicTableofElements.csv")

_THRESHOLD_PAIR_ELECTRON = 2.044014E+06  # eV
_THRESHOLD_PAIR_ATOM = 1.022007E+06  # eV


class MaterialFactory:
    """
    Class for creation compound by mass fraction
    """
    element_symbols = None

    def __init__(self):
        self.elements = []
        self.weights = []

    def add_element(self, element: Union[str, int], weight: float) -> 'MaterialFactory':
        """
        Add element and its mass fraction

        Parameters
        ----------
        element :  Union[str, int]
                atomic number or symbol of element
        weight : float
                 mass fraction, unnormed
        """
        self.elements.append(element)
        self.weights.append(weight)
        return self

    def add_material(self, material) -> 'MaterialFactory':
        """
        Add another material
        """
        self.elements += material.elements_by_Z
        self.weights += material.weights
        return self

    def build(self) -> "Material":
        """
        Build material from partition
        """
        elements = list(
            map(
                lambda x: self.get_element_from_symbol(x) if isinstance(x, str) else x,
                self.elements
            )
        )
        return Material(elements, weights=self.weights)

    @classmethod
    def from_formula(cls, formula) -> 'Material':
        """
        Create material from chemical formula

        Parameters
        ----------
        formula :
            Chemical formulas for compounds should be entered in standard chemical notation,
            with appropriate upper and lower case. However, because of hardware limitations,
            subscripts must be written on line. For example, the formula for calcium tungstate
            must be entered as CaWO4.
            Parentheses, spaces and dots may not be used.
            For example, the formula for calcium phosphate must be entered as Ca3P2O8 (and
            not as Ca3(PO4)2).

        Returns
        -------
        material : Material

        """
        name_list, value_list = [], []
        value = ''
        i = 0
        n = len(formula)
        while i < n:
            s = formula[i]
            # print(s)
            if s.isupper():
                i += 1
                if value != '':
                    value_list.append(int(value))
                if i == n:
                    name_list.append(s)
                    value_list.append(1)
                    break

                if formula[i].isdigit():
                    name_list.append(s)
                elif formula[i].isupper():
                    name_list.append(s)
                    value_list.append(1)
                else:
                    name_list.append(s + formula[i])
                    i += 1
                    if i == n:
                        value_list.append(1)
                        break

                    if formula[i].isupper():
                        value_list.append(1)
            elif s.isdigit():
                value += s
                i += 1
                if i == n:
                    value_list.append(int(value))
                    break

                if formula[i].isupper():
                    value_list.append(int(value))
                    value = ''
            else:
                break
        elements = []
        for name in name_list:
            elements.append(MaterialFactory.get_element_from_symbol(name))
        atomic_mass = MaterialFactory.get_elements_mass_list(elements)

        weights = []
        for mass, value in zip(atomic_mass, value_list):
            weights.append(mass * value)
        return Material(elements, weights)

    @classmethod
    def _prepare_element_symbol(cls):
        cls.element_symbols = {}
        with open(PERIODIC_TABLE_PATH, newline='', encoding='UTF-8') as csvfile:
            reader = csv.reader(csvfile, delimiter=",")
            next(reader)
            for row in reader:
                Z, _, symbol = row[:3]
                cls.element_symbols[symbol] = int(Z)

    @classmethod
    def get_element_from_symbol(cls, element: str) -> int:
        """
        Get atomic number of element based on symbol
        """
        if cls.element_symbols is None:
            cls._prepare_element_symbol()
        return cls.element_symbols[element]

    @staticmethod
    def get_element_mass(element: int) -> float:
        """
        Get element atomic mass in amu
        """
        if element <= 0 or element > 100:
            raise ValueError("Element must be from 1 ot 100")
        with tables.open_file(NIST_XCOM_HDF5_PATH) as h5file:
            group_name = f"/Z{str(element).rjust(3, '0')}"
            table = h5file.get_node(group_name, "data")
            return table.attrs['AtomicWeight']

    @staticmethod
    def get_elements_mass_list(elements: List[int]) -> np.ndarray:
        """
        Get list of elements atomic mass in amu
        """
        result = np.zeros(len(elements))
        with tables.open_file(NIST_XCOM_HDF5_PATH) as h5file:
            for indx, element in enumerate(elements):
                if element <= 0 or element > 100:
                    raise ValueError("Element must be from 1 ot 100")
                group_name = f"/Z{str(element).rjust(3, '0')}"
                table = h5file.get_node(group_name, "data")
                result[indx] = table.attrs['AtomicWeight']
            return result


class Material:
    """
    Define material for attenuation calculation
    """

    def __init__(self, elements: List[int], weights: Optional[List[float]] = None):
        """
        Parameters
        ----------

        elements : List[int]
                    List of atomic number of element
        weights : Optional[List[float]]
                    List of mass fraction of elements, not required for single element
        """
        self.elements_by_Z = elements
        if weights is not None:
            assert (len(self.elements_by_Z) == len(weights))
            sum_ = sum(weights)
            weights = list(map(lambda x: x / sum_, weights))
        self.weights = weights

    def __len__(self):
        return len(self.elements_by_Z)


def make_log_log_spline(x: np.ndarray, y: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """
    Create spline of log-log data
    """
    log_x = np.log(x)
    log_y = np.log(y)
    cs = CubicSpline(x=log_x, y=log_y, bc_type='natural')
    linear = interp1d(x=log_x, y=log_y,
                      kind='linear', fill_value="extrapolate")

    def spliner(x_sample: np.ndarray) -> np.ndarray:
        if np.min(x_sample) < np.min(x) or np.max(x_sample) < np.max(x):
            warnings.warn("Energy requested is outside of tabulated data, "+
                          "using linear extrapolation", UserWarning)
        return np.where(
            np.logical_and(x_sample > np.min(x), x_sample < np.max(x)),
            np.exp(cs(np.log(x_sample))),  # Use cubic within data range
            np.exp(linear(np.log(x_sample)))  # Extrapolate with linear
            )

    return spliner


def _interpolateAbsorptionEdge(data) -> Callable[[np.ndarray], np.ndarray]:
    data, h5file, group = data

    data_K = h5file.get_node(group, "K").read()
    cubicSplineThreshold = np.max(data_K[NameProcess.ENERGY]) * 1e6
    x = data[NameProcess.ENERGY]
    y = data[NameProcess.PHOTOELECTRIC]
    indx = x > cubicSplineThreshold

    cs = CubicSpline(np.log(x[indx]), np.log(y[indx]), bc_type='natural')
    linear = interp1d(np.log(x), np.log(y),
                      kind='linear', fill_value="extrapolate")

    def spliner(x_sample: np.ndarray) -> np.ndarray:
        if np.min(x_sample) < np.min(x) or np.max(x_sample) < np.max(x):
            warnings.warn("Energy requested is outside of tabulated data, "+
                          "using linear extrapolation", UserWarning)
        return np.where(
            np.logical_and(x_sample > cubicSplineThreshold, x_sample < np.max(x)),
            np.exp(cs(np.log(x_sample))),
            np.exp(linear(np.log(x_sample)))
        )

    return spliner


def make_pair_interpolator(x: np.ndarray, y: np.ndarray,
                           threshold: float) -> Callable[[np.ndarray], np.ndarray]:
    """
    Create spline of linearized log-log data
    """
    indx = x > threshold
    cs = CubicSpline(x=np.log(x[indx]),
                     y=np.log(y[indx] / (x[indx]*(x[indx] - threshold))**3),
                     bc_type='natural')
    linear = interp1d(x=np.log(x[indx]), y=np.log(y[indx]),
                      kind='linear', fill_value='extrapolate')

    def spliner(x_sample: np.ndarray) -> np.ndarray:
        if np.max(x_sample) < np.max(x):
            warnings.warn("Energy requested is outside of tabulated data, "+
                          "using linear extrapolation", UserWarning)
        indx = x_sample > threshold
        y = np.zeros(x_sample.shape[0])
        y[indx] = np.where(
            x_sample[indx] < np.max(x),
            np.exp(cs(np.log(x_sample[indx]))) * (x_sample[indx]*(x_sample[indx] - threshold))**3,
            np.exp(linear(np.log(x_sample[indx])))
        )
        return y

    return spliner


def create_coherent_interpolator(data: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """Create interpolator for coherent scattering from tabulated data"""
    return make_log_log_spline(data[NameProcess.ENERGY],
                               data[NameProcess.COHERENT])


def create_incoherent_interpolator(data: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """Create interpolator for incoherent scattering from tabulated data"""
    return make_log_log_spline(data[NameProcess.ENERGY],
                               data[NameProcess.INCOHERENT])


def create_pair_atom_interpolator(data: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """Create interpolator for pair production by an atomic nucleus from tabulated data"""
    return make_pair_interpolator(data[NameProcess.ENERGY],
                                  data[NameProcess.PAIR_ATOM],
                                  threshold=_THRESHOLD_PAIR_ATOM)


def create_pair_electron_interpolator(data: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """Create interpolator for pair production by an electron from tabulated data"""
    return make_pair_interpolator(data[NameProcess.ENERGY],
                                  data[NameProcess.PAIR_ELECTRON],
                                  threshold=_THRESHOLD_PAIR_ELECTRON)


def create_photoelectric_interpolator(data: np.ndarray,
                                      absorption_edge=False) -> Callable[[np.ndarray], np.ndarray]:
    """Create interpolator for photoelectric absorption from tabulated data"""
    if absorption_edge:
        return _interpolateAbsorptionEdge(data)

    return make_log_log_spline(data[NameProcess.ENERGY],
                               data[NameProcess.PHOTOELECTRIC])


class Interpolators:
    """Class for interpolating functions from tabulated cross-section data"""
    def __init__(self):
        self.cache = {}
        self.h5file = tables.open_file(NIST_XCOM_HDF5_PATH)

    def __del__(self):
        self.h5file.close()

    def get_interpolators(self, element: int):
        """Create interpolators for the provided element"""
        if element <= 0 or element > 100:
            raise ValueError("Element must be from 1 ot 100")
        try:
            return self.cache[element]
        except KeyError:
            group_name = f"/Z{str(element).rjust(3, '0')}"
            table = self.h5file.get_node(group_name, "data")
            data = table.read()
            if table.attrs['AbsorptionEdge']:
                group_name += "/AbsorptionEdge"
                data_phot = (data, self.h5file, group_name)
            else:
                data_phot = data
            temp = {
                NameProcess.COHERENT: create_coherent_interpolator(data),
                NameProcess.INCOHERENT: create_incoherent_interpolator(data),
                NameProcess.PAIR_ELECTRON: create_pair_electron_interpolator(data),
                NameProcess.PAIR_ATOM: create_pair_atom_interpolator(data),
                NameProcess.PHOTOELECTRIC:
                    create_photoelectric_interpolator(data_phot,
                                                      absorption_edge=table.attrs[
                                                          'AbsorptionEdge'])
            }
            self.cache[element] = temp
            return temp
