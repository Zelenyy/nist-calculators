import unittest

import xcom
from xcom._data_converter import NameProcess
from xcom.interpolators import Interpolators


class XCOMCase(unittest.TestCase):
    def test_interpolators(self):
        inter = Interpolators().get_interpolators(1)
        value = inter[NameProcess.COHERENT](1e6)
        self.assertAlmostEquals(value, 4.625E-06)

    def test_interpolators(self):
        for i in range(1,100):
            Interpolators().get_interpolators(i)

    def test_attenuation(self):
        xcom.calculate_attenuation(xcom.Material([1]))

    def test_parser(self):
        material = xcom.MaterialFactory.from_formula("H2O")
        print(material.elements_by_Z, material.weights)


if __name__ == '__main__':
    unittest.main()
