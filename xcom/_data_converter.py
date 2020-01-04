import logging
import os
import sys
import numpy as np
import tables

ROOT_PATH = os.path.dirname(__file__)

class NameProcess:
    ENERGY = "energy"
    COHERENT = "coherent"
    INCOHERENT = "incoherent"
    PHOTOELECTRIC = "photoelectric"
    PAIR_ATOM = "pair_atom"
    PAIR_ELECTRON = "pair_electron"

# Atomic weight
ATWTS = [
    1.00794, 4.002602, 6.941, 9.012182,
    10.811, 12.011, 14.00674, 15.9994,
    18.9984032, 20.1797, 22.989768, 24.3050,
    26.981539, 28.0855, 30.973762, 32.066,
    35.4527, 39.948, 39.0983, 40.078,
    44.955910, 47.88, 50.9415, 51.9961,
    54.93805, 55.847, 58.93320, 58.69,
    63.546, 65.39, 69.723, 72.61,
    74.92159, 78.96, 79.904, 83.80,
    85.4678, 87.62, 88.90585, 91.224,
    92.90638, 95.94, 97.9072, 101.07,
    102.9055, 106.42, 107.8682, 112.411,
    114.82, 118.710, 121.75, 127.60,
    126.90447, 131.29, 132.90543, 137.327,
    138.9055, 140.115, 140.90765, 144.24,
    144.9127, 150.36, 151.965, 157.25,
    158.92534, 162.50, 164.93032, 167.26,
    168.93421, 173.04, 174.967, 178.49,
    180.9479, 183.85, 186.207, 190.2,
    192.22, 195.08, 196.96654, 200.59,
    204.3833, 207.2, 208.98037, 208.9824,
    209.9871, 222.0176, 223.0197, 226.0254,
    227.0278, 232.0381, 231.03588, 238.0289,
    237.0482, 239.0522, 243.0614, 247.0703,
    247.0703, 251.0796, 252.083, 257.0951]

SOURCE = 'MDATX3'


class XCOMData:
    IZ = None
    ATWT1 = None
    MAXEDG = None
    MAXE = None
    IDG = None
    ADG = None
    EDGEN = None
    E = None
    SCATCO = None
    SCATIN = None
    PHOT = None
    PAIRAT = None
    PAIREL = None
    ENG = None
    PHC = None
    MAXK = None


def read_data_for_element(path_to_xcom_data, Z) -> XCOMData:
    Z = int(Z)
    result = XCOMData()
    with open(os.path.join(path_to_xcom_data, SOURCE + '.' + str(Z).rjust(3, '0'))) as fin:
        file = []
        for line in fin.readlines():
            file += line.split()
    result.IZ = int(file[0])
    result.ATWT1 = float(file[1])
    result.MAXEDG, result.MAXE = tuple(map(int, file[2:4]))
    file_start = 4
    if result.MAXEDG > 0:
        result.IDG = np.fromiter(map(int, file[4:4 + result.MAXEDG]), dtype='i')[::-1]
        result.ADG = file[4 + result.MAXEDG:4 + result.MAXEDG + result.MAXEDG][::-1]
        result.EDGEN = np.fromiter(
            map(float, file[4 + result.MAXEDG + result.MAXEDG: 4 + result.MAXEDG + result.MAXEDG + result.MAXEDG]),
            dtype='d')[::-1]
        file_start = 4 + result.MAXEDG + result.MAXEDG + result.MAXEDG
    result.E = np.fromiter(map(float, file[file_start:file_start + result.MAXE]), dtype=[(NameProcess.ENERGY, 'd')])
    file_start += result.MAXE
    result.SCATCO = np.fromiter(map(float, file[file_start:file_start + result.MAXE]), dtype=[(NameProcess.COHERENT, 'd')])
    file_start += result.MAXE
    result.SCATIN = np.fromiter(map(float, file[file_start:file_start + result.MAXE]), dtype=[(NameProcess.INCOHERENT, 'd')])
    file_start += result.MAXE
    result.PHOT = np.fromiter(map(float, file[file_start:file_start + result.MAXE]), dtype=[(NameProcess.PHOTOELECTRIC, 'd')])
    file_start += result.MAXE
    result.PAIRAT = np.fromiter(map(float, file[file_start:file_start + result.MAXE]), dtype=[(NameProcess.PAIR_ATOM, 'd')])
    file_start += result.MAXE
    result.PAIREL = np.fromiter(map(float, file[file_start:file_start + result.MAXE]), dtype=[(NameProcess.PAIR_ELECTRON, 'd')])
    file_start += result.MAXE
    result.MAXK = 0

    if result.MAXEDG > 0:
        result.MAXK = result.MAXE - result.IDG[result.MAXEDG - 1] + 1
        LAX = int(file[file_start])
        file_start += 1
        KMX = np.fromiter(map(int, file[file_start:file_start + LAX]), dtype='i')
        file_start += LAX
        result.ENG = []
        for L in range(LAX):
            result.ENG.append(list(map(float, file[file_start: file_start + KMX[L]])))
            file_start += KMX[L]
        result.PHC = []
        for L in range(LAX):
            result.PHC.append(list(map(float, file[file_start: file_start + KMX[L]])))
            file_start += KMX[L]
        result.ENG = np.array(result.ENG)
        result.PHC = np.array(result.PHC)
    return result

def drop_doubling(data):
    _, indx = np.unique(data[NameProcess.ENERGY], return_index=True)
    return data[indx]

def convertMDATX3toHDF5(path_to_xcom_data):
    filename = os.path.join(ROOT_PATH, 'data', 'NIST_XCOM.hdf5')
    with tables.open_file(filename, "w") as h5file:

        filters = tables.Filters(complevel=3, fletcher32=True)

        for Z in range(1, 101):
            data = read_data_for_element(path_to_xcom_data, Z)
            MAXEDG = data.MAXEDG
            name_grp = "Z" + str(data.IZ).rjust(3, '0')
            group = h5file.create_group(h5file.root, name_grp)

            temp_tuple = (data.E, data.SCATCO, data.SCATIN, data.PHOT, data.PAIRAT, data.PAIREL)
            temp_names = [name for name in map(lambda x: x.dtype.names[0], temp_tuple)]
            temp_data = [it[name] for name, it in zip(temp_names, temp_tuple)]
            data_to_file = np.core.records.fromarrays(temp_data,names=",".join(temp_names))

            delta = np.diff(data_to_file[NameProcess.ENERGY])
            if np.any(delta <= 0):
                print("Drop doubling for {}".format(data.IZ))
                data_to_file = drop_doubling(data_to_file)


            table = h5file.create_table(group, "data", obj=data_to_file, filters=filters)
            table.attrs['energy_unit'] = 'eV'
            table.attrs["cross_section_unit"] = 'barn/atom'
            table.attrs['Z'] = data.IZ
            table.attrs['AtomicWeight'] = ATWTS[data.IZ - 1]
            table.attrs['AtomicWeightMDATX3'] = data.ATWT1
            if MAXEDG > 0:
                table.attrs['AbsorptionEdge'] = True
            else:
                table.attrs['AbsorptionEdge'] = False
            table.attrs['NumberAbsorptionEdgeEnergy'] = MAXEDG
            table.flush()
            if MAXEDG > 0:
                edge_group = h5file.create_group(group, 'AbsorptionEdge')
                name = ['index', 'name', 'EDGEN']
                data_to_file = np.core.records.fromarrays([data.IDG, np.array(data.ADG, dtype='S'), data.EDGEN], names=",".join(name))
                table = h5file.create_table(edge_group, "info", obj=data_to_file, filters=filters)
                table.attrs["energy_unit"] = 'eV'
                table.flush()
                for e, p, name in zip(data.ENG, data.PHC, data.ADG):
                    data_to_file = np.core.records.fromarrays([e,p], names="energy,photoelectric")
                    table = h5file.create_table(edge_group, name, obj=data_to_file, filters=filters)
                    table.attrs["energy_unit"] = 'MeV'
                    table.attrs["cross_section_unit"] = 'barn/atom'
                    table.flush()
    return filename


if __name__ == '__main__':
    path = sys.argv[1]
    print(path)
    convertMDATX3toHDF5(path)
