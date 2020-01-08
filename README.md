# NIST calculators
My python implementaion of some NIST calculator.

Now implemented:

* [NIST XCOM](https://www.nist.gov/pml/xcom-photon-cross-sections-database) --- photon cross-sections database
* [NIST_STAR](https://www.nist.gov/pml/stopping-power-range-tables-electrons-protons-and-helium-ions) --- stopping-power & range tables for electrons, protons, and helium ions

## Installation

```shell script
pip install nist-calculators
```
## Usage

### XCOM

See [tutorial](./examples/XCOM.ipynb).

Example of usage:
```python
import xcom

data = xcom.calculate_cross_section("He", [1e6, 2e6, 5e6]) # Energy in eV
# or
import numpy as np
energy = np.array([1e6, 2e6, 5e6])
Z = 2
data = xcom.calculate_cross_section(Z, energy)
```

### ESTAR
See [tutorial](./examples/ESTAR.ipynb).

Example of usage:
```python
from star import electron
hydrogen = electron.PredefinedMaterials.HYDROGEN
data = electron.calculate_stopping_power(hydrogen, energy=[1e2,2e2,3e2]) # Energy in MeV
```
### PSTAR
See [tutorial](./examples/PSTAR.ipynb).

Example of usage:
```python
from star import ProtonSTARCalculator, ProtonMaterials

material = ProtonMaterials.BERYLLIUM
calculator = ProtonSTARCalculator(material)
total = calculator.calculate_total_stopping_powers( [10, 20, 50]) # Energy in MeV
```

### ASTAR
See [tutorial](./examples/ASTAR.ipynb).

Example of usage:
```python
from star import AlphaSTARCalculator, AlphaMaterials

material = AlphaMaterials.BERYLLIUM
calculator = AlphaSTARCalculator(material)
total = calculator.calculate_total_stopping_powers( [10, 20, 50]) # Energy in MeV
```