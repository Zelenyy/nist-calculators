import os
# from enum import Enum
# from typing import Union, List, Tuple, Callable
# import numpy as np
# from scipy.interpolate import CubicSpline, interp1d

# from xcom._dataconverter import convertMDATX3toHDF5
ROOT_PATH = os.path.dirname(__file__)
DATA_PATH = os.path.join(ROOT_PATH, 'data')
NIST_XCOM_HDF5_PATH = os.path.join(DATA_PATH, 'NIST_XCOM.hdf5')
# PERIODIC_TABLE_PATH = os.path.join(DATA_PATH, 'periodic_table.csv')


# if not os.path.exists(NIST_XCOM_HDF5_PATH):
#     convertMDATX3toHDF5
#
# class Element:
#     def __init__(self, element):
#         pass
#
# class Material:
#     pass
#
#
# class DataStorage:
#
#     _instance = None
#     def __new__(class_, *args, **kwargs):
#         if not isinstance(class_._instance, class_):
#             class_._instance = object.__new__(class_, *args, **kwargs)
#         return class_._instance
#
# class Unit(Enum):
#     """
#     Unit for data. Only alone elemental.
#
#     ALL_BARN_BY_ATOM: barn/atom for all,
#     OTHER_BARN_BY_ATOM - total data in cm^2/gr, other in barn/atom,
#     ALL_CM2_BY_GRAMM - cm^2/gr for all
#     """
#     ALL_BARN_BY_ATOM = 1
#     OTHER_BARN_BY_ATOM = 2
#     ALL_CM2_BY_GRAMM = 3
#
# class XCOM:
#     _TRESHOLD_PAIR_ELECTRON = 1.022007E+06 # eV
#     _TRESHOLD_PAIR_ATOM = 2.044014E+06 # eV
#     _AVOGADRO = 0.60221367 # * 1e+24 mole^-1
#
#     def __init__(self, material, path_to_data = None):
#         self.phys_process = ['Coherent', 'Incoherent', 'Photoelectric', 'PairAtom', 'PairElectron']
#
#     def calculate_coherent(self, energy: np.ndarray = None, unit : Unit = Unit.ALL_BARN_BY_ATOM) -> np.ndarray:
#         pass
#
#     def calculate_incoherent(self, energy: np.ndarray = None, unit : Unit = Unit.ALL_BARN_BY_ATOM) -> np.ndarray:
#         pass
#
#     def calculate_photoelectric(self, energy: np.ndarray = None, unit : Unit = Unit.ALL_BARN_BY_ATOM) -> np.ndarray:
#         pass
#
#     def calculate_pair_atom(self, energy: np.ndarray = None, unit : Unit = Unit.ALL_BARN_BY_ATOM) -> np.ndarray:
#         pass
#
#     def calculate_pair_electron(self, energy: Union[np.ndarray, list] = None, unit : Unit = Unit.ALL_BARN_BY_ATOM) -> np.ndarray:
#         pass
#
#     def calculate(self, energy: np.ndarray = None, unit : Unit = Unit.ALL_BARN_BY_ATOM) -> np.ndarray:
#         """
#         Расчитываем сечения и длинны пробегов для заданной энергии
#
#         Parameters
#         ----------
#         energy : energies of gamma-quant, used `ENERGY_GRID_DEFAULT` by default
#
#         Returns
#         -------
#         outData : data frame with cross-section and attenuatuon
#         """
#         if not self.cacheFreshness:
#             self._prepareInterpolation()
#
#         if energy is None:
#             energy = ENERGY_GRID_DEFAULT
#
#         n = len(energy)
#         energy = np.array(energy, dtype=[("energy", "d")])
#         total = np.zeros(n, dtype=[("total", "d")])
#         total_without_coherent = np.zeros(n, [("total_without_coherent", "d")])
#
#         if len(self.substanceSet) > 1:
#             for process in self.phys_process:
#                 outData[process] = np.zeros(n)
#                 for z in self.substanceSet:
#                     outData[process] += self.interpolation[z][process](energy) * self.weights[z] * self._AVOG / \
#                                         self.data[z].attrs['AtomicWeight']
#                 outData['Total'] += outData[process]
#             outData['TotalWithoutCoherent'] = outData['Total'] - outData['Coherent']
#         else:
#             z = list(self.substanceSet)[0]
#             for process in self.phys_process:
#                 outData[process] = self.interpolation[z][process](energy)
#                 outData['Total'] += outData[process]
#             if self.unit == 2:
#                 outData['Total'] = outData['Total'] * self._AVOG / self.data[z].attrs['AtomicWeight']
#                 outData['TotalWithoutCoherent'] = outData['TotalWithoutCoherent'] * self._AVOG / self.data[z].attrs[
#                     'AtomicWeight']
#             elif self.unit == 3:
#                 for name in self.phys_process + ['Total', 'TotalWithoutCoherent']:
#                     outData[name] = outData[name] * self._AVOG / self.data[z].attrs['AtomicWeight']
#
#         return outData
#
#
#
# class XCOM1:
#
#     def __init__(self):
#         self.data = h5py.File(NIST_XCOM_HDF5_PATH, 'r')
#         self.substanceList = []
#         self.weigthList = []
#         self.phys_process = ['Coherent', 'Incoherent', 'Photoelectric', 'PairAtom', 'PairElectron']
#         self.unit = 1
#
#     def addElement(self, Z: Union[int, str], weight: float = None):
#         """
#         Добавляет новый элемент к сущетсвующим
#
#         Parameters
#         ----------
#         Z : порядковый номер элемента или символ элемента
#         weight : указывает массовую долю элемента, по умолчанию использует массу элемента
#
#         Returns
#         -------
#         self
#         """
#         if type(Z) == str:
#             Z = self._findNumberOfElement(Z)
#         if (type(Z) == int) and (0 < Z < 101):
#             # TODO raise
#             if weight is None:
#                 weight = self.data[str(Z).rjust(3, '0')].attrs['AtomicWeight']
#             self._addSubstance([Z], [weight])
#         else:
#             pass  # TODO raise
#
#         return self
#
#     def addMaterial(self, formula: str, weight: float = None):
#         """
#         Добавляет новую молекулу для расчетов
#
#         Parameters
#         ----------
#         formula : Chemical formulas for compounds should be entered in standard chemical notation,
#             with appropriate upper and lower case. However, because of hardware limitations,
#             subscripts must be written on line. For example, the formula for calcium tungstate must be entered as CaWO4.
#             Parentheses, spaces and dots may not be used.
#             For example, the formula for calcium phosphate must be entered as Ca3P2O8 (and not as Ca3(PO4)2).
#         weight : указывает массовую долю элемента, по умолчанию использует массу элемента
#
#         Returns
#         -------
#
#         """
#         zList, weightFactor = self._parseFormula(formula)
#         # print(zList, weightFactor)
#         if weight is None:
#             weight = []
#             for indx, wf in enumerate(weightFactor):
#                 weight.append(wf * self.data[str(zList[indx]).rjust(3, '0')].attrs['AtomicWeight'])
#         else:
#             sum_ = sum(weightFactor)
#             weight = [weight * wf / sum_ for wf in weightFactor]
#         self._addSubstance(zList, weight)
#         return self
#
#     def addCompaund(self, compaund: dict):
#         for key, value in compaund.items():
#             self.addMaterial(key, value)
#         return self
#
#     def _addSubstance(self, Z: List[int], weight: List[float]) -> None:
#         self.substanceList += Z
#         self.weigthList += weight
#         self.cacheFreshness = False
#
#     def _interpolateLogLog(self, x: np.ndarray, y: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
#         """
#
#         Parameters
#         ----------
#         x
#         y
#
#         Returns
#         -------
#
#         """
#         x = np.log(x)
#         y = np.log(y)
#         cs = CubicSpline(x=x, y=y)
#
#         def makeSpliner(cs):
#             def spliner(x: np.ndarray) -> np.ndarray:
#                 return np.exp(cs(np.log(x)))
#
#             return spliner
#
#         return makeSpliner(cs)
#
#     def _interpolateAbsorptionEdge(self, data: h5py.Dataset) -> Callable[[np.ndarray], np.ndarray]:
#         cubicSplineTreshold = np.max(data['AbsorptionEdge']['K']['Energy']) * 1e6
#         x = np.log(data['Energy'])
#         y = np.log(data['Photoelectric'])
#         indx = x > np.log(cubicSplineTreshold)
#         cs = CubicSpline(x[indx], y[indx])
#
#         indx = np.logical_not(indx)
#         # print(x[indx], y[indx])
#         linear = interp1d(x[indx], y[indx], kind='linear')
#
#         def makeSpliner(cs, linear, treshold):
#             def spliner(x: np.ndarray) -> np.ndarray:
#                 x = np.log(x)
#                 indx = x > np.log(treshold)
#                 y = np.zeros(x.shape[0])
#                 y[indx] = cs(x[indx])
#                 indx = np.logical_not(indx)
#                 y[indx] = linear(x[indx])
#                 return np.exp(y)
#
#             return spliner
#
#         return makeSpliner(cs, linear, cubicSplineTreshold)
#
#     def _interpolatePair(self, x: np.ndarray, y: np.ndarray, treshold: float) -> Callable[[np.ndarray], np.ndarray]:
#         x = (1 - (x / treshold)) ** 3
#         indx = y > 0
#         x = x[indx][::-1]
#         y = np.log(y[indx])[::-1]
#
#         cs = CubicSpline(x=x, y=y)
#
#         def makeSpliner(cs, treshold):
#             def spliner(x: np.ndarray) -> np.ndarray:
#                 x = (1 - (x / treshold)) ** 3
#                 indx = x <= 0
#                 y = np.zeros(x.shape[0])
#                 y[indx] = np.exp(cs(x[indx]))
#                 return y
#
#             return spliner
#
#         return makeSpliner(cs, treshold)
#
#     def _prepareInterpolation(self):
#         if len(self.substanceList) == 0:
#             print('Add material')
#             # TODO raise
#         substanceList = [str(z).rjust(3, '0') for z in self.substanceList]
#
#         # Processing input substance and weigth
#         self.substanceSet = set(substanceList)
#         weights = {z: 0 for z in self.substanceSet}
#         sum_ = 0
#         for z, w in zip(substanceList, self.weigthList):
#             # massWeight = w * self.data[z].attrs['AtomicWeight']
#             weights[z] += w
#         sum_ += sum(self.weigthList)
#         self.weights = {key: value / sum_ for key, value in weights.items()}
#
#         # Create interpolation
#         self.interpolation = {}
#         for z in self.substanceSet:
#             if not (z in self.interpolation.keys()):
#                 temp = {}
#                 temp['Coherent'] = self._interpolateLogLog(self.data[z]['Energy'].value, self.data[z]['Coherent'].value)
#                 temp['Incoherent'] = self._interpolateLogLog(self.data[z]['Energy'].value,
#                                                              self.data[z]['Incoherent'].value)
#                 temp['PairAtom'] = self._interpolatePair(self.data[z]['Energy'].value, self.data[z]['PairAtom'].value,
#                                                          treshold=self._TRESHOLD_PAIR_ATOM)
#                 temp['PairElectron'] = self._interpolatePair(self.data[z]['Energy'].value,
#                                                              self.data[z]['PairElectron'].value,
#                                                              treshold=self._TRESHOLD_PAIR_ELECTRON)
#                 if self.data[z].attrs['AbsorptionEdge']:
#                     temp['Photoelectric'] = self._interpolateAbsorptionEdge(self.data[z])
#                 else:
#                     temp['Photoelectric'] = self._interpolateLogLog(self.data[z]['Energy'].value,
#                                                                     self.data[z]['Photoelectric'].value)
#                 self.interpolation[z] = temp
#                 # TODO global cache
#
#
#
#     def _findNumberOfElement(self, element: str) -> int:
#         data = pd.read_csv(PERIODIC_TABLE_PATH, sep=',', skipinitialspace=True)
#         indx = data.loc[lambda data: data['Символ'] == element]
#         if len(indx['Атомный номер'].values) != 0:
#             return int(indx['Атомный номер'].values[0])
#         else:
#             pass  # TODO raise
#
#     def _parseFormula(self, formula: str) -> Tuple[List[int], List[float]]:
#         data = pd.read_csv(PERIODIC_TABLE_PATH, sep=',', skipinitialspace=True)
#         nameList, valueList = [], []
#         value = ''
#         i = 0
#         n = len(formula)
#         while i < n:
#             s = formula[i]
#             # print(s)
#             if s.isupper():
#                 i += 1
#                 if value != '':
#                     valueList.append(int(value))
#                 if i == n:
#                     nameList.append(s)
#                     valueList.append(1)
#                     break
#                 else:
#                     if formula[i].isdigit():
#                         nameList.append(s)
#                     elif formula[i].isupper():
#                         nameList.append(s)
#                         valueList.append(1)
#                     else:
#                         nameList.append(s + formula[i])
#                         i += 1
#                         if i == n:
#                             valueList.append(1)
#                             break
#                         elif formula[i].isupper():
#                             valueList.append(1)
#             elif s.isdigit():
#                 value += s
#                 i += 1
#                 if i == n:
#                     valueList.append(int(value))
#                     break
#                 else:
#                     if formula[i].isupper():
#                         valueList.append(int(value))
#                         value = ''
#             else:
#                 break
#         # print(nameList, valueList)
#         zList = []
#         for name in nameList:
#             zList.append(self._findNumberOfElement(name))
#         return zList, valueList
#
# ENERGY_GRID_DEFAULT = np.array([1.0E+03, 1.5E+03, 2.0E+03, 3.0E+03, 4.0E+03, 5.0E+03,
#                                     6.0E+03, 8.0E+03, 1.0E+04, 1.5E+04, 2.0E+04, 3.0E+04, 4.0E+04,
#                                     5.0E+04, 6.0E+04, 8.0E+04, 1.0E+05, 1.5E+05, 2.0E+05, 3.0E+05,
#                                     4.0E+05, 5.0E+05, 6.0E+05, 8.0E+05, 1.0E+06, 1.022E+06, 1.25E+06,
#                                     1.5E+06, 2.0E+06, 2.044E+06, 3.0E+06, 4.0E+06, 5.0E+06, 6.0E+06,
#                                     7.0E+06, 8.0E+06, 9.0E+06, 1.0E+07, 1.1E+07, 1.2E+07, 1.3E+07,
#                                     1.4E+07, 1.5E+07, 1.6E+07, 1.8E+07, 2.0E+07, 2.2E+07, 2.4E+07,
#                                     2.6E+07, 2.8E+07, 3.0E+07, 4.0E+07, 5.0E+07, 6.0E+07, 8.0E+07,
#                                     1.0E+08, 1.5E+08, 2.0E+08, 3.0E+08, 4.0E+08, 5.0E+08, 6.0E+08,
#                                     8.0E+08, 1.0E+09, 1.5E+09, 2.0E+09, 3.0E+09, 4.0E+09, 5.0E+09,
#                                     6.0E+09, 8.0E+09, 1.0E+10, 1.5E+10, 2.0E+10, 3.0E+10, 4.0E+10,
#                                     5.0E+10, 6.0E+10, 8.0E+10, 1.0E+11], dtype='d')
#
#
# # import matplotlib.pyplot as plt
#
#
# # def main():
# #     # xcom = XCOMDataCalculator().addMaterial('H2O').setUnit(3).calculate()
# #     energy = np.linspace(1e3, 1e4, 1000)
# #     result = XCOMDataCalculator().addElement(50).setUnit(3).calculate(energy)
# #     # print(result.head(30))
# #     plt.plot(result['Energy'] / 1000, result['Total'])
# #     plt.show()
# #
# #     return 0
# #
# #
# # if __name__ == '__main__':
# #     main()
