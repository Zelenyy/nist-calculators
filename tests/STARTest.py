import logging
import unittest

import star
from star import ProtonPredefinedMaterials, ProtonSTARCalculator
from star.electron import PredefinedMaterials


class STARCase(unittest.TestCase):
    def test_electron(self):
        logging.root.setLevel(logging.DEBUG)
        data = star.electron.calculate_stopping_power(PredefinedMaterials.HYDROGEN)

    def test_electron_table(self):
        logging.root.setLevel(logging.DEBUG)
        data = star.electron.calculate_estar_table(PredefinedMaterials.HYDROGEN)

    def test_electron_all(self):
        logging.root.setLevel(logging.DEBUG)
        for material in PredefinedMaterials:
            star.electron.calculate_estar_table(material)

    def test_proton(self):
        material = ProtonPredefinedMaterials.BERYLLIUM
        calculator = ProtonSTARCalculator(material)

    def test_proton_all(self):
        for material in ProtonPredefinedMaterials:
            calculator = ProtonSTARCalculator(material)
            calculator.calculate_table()