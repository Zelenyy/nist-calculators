import unittest

import star
from star.electron import PredefinedMaterials


class STARCase(unittest.TestCase):
    def test_electron(self):
        data = star.electron.calculate_stopping_power(PredefinedMaterials.HYDROGEN)