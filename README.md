# NIST calculators
My python implementaion of some NIST calculator.

Now implemented:

* [NIST XCOM](https://www.nist.gov/pml/xcom-photon-cross-sections-database) --- photon cross-sections database

## Installation

```shell script
pip install nist-calculators
```

## XCOM

See [example](./notebooks/XCOM.ipynb).

Usage:
```python
import xcom

data = xcom.calculate_cross_section("He", [1e6, 2e6, 5e6])
# or
import numpy as np
energy = np.array([1e6, 2e6, 5e6])
Z = 2
data = xcom.calculate_cross_section(Z, energy)
```

