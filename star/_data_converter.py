import os
import sys

import numpy as np
import tables

from star.proton_materials import ProtonMaterials

ROOT_PATH = os.path.dirname(__file__)
DATA_PATH = os.path.join(ROOT_PATH, 'data')
NIST_STAR_HDF5_PATH = os.path.join(DATA_PATH, 'NIST_STAR.hdf5')


class Names:
    IONISATION_POTENTIAL = "ionisation_potential"
    DENSITY = "density"
    ZAG = "zag"
    ENERGY = "energy"
    ELEMENT = "element"
    FRACTION = "fraction"
    STOPPING_POWER_COLLISION_DELTA = "stopping_power_collision_delta"
    STOPPING_POWER_RADIATIVE = "stopping_power_radiative"
    STOPPING_POWER_TOTAL = "stopping_power_total"
    CSDA_RANGE = "csda_range"
    RADIATION_YIELD = "radiation_yield"
    DENSITY_EFFECT = "density_effect"


def convert_fedat_to_table(path):
    LKIN = 113

    dtype_radiation_loss = np.dtype(
        [
            (Names.ELEMENT, "i"),
            (Names.STOPPING_POWER_RADIATIVE, "d", (LKIN,))
        ]
    )
    data_radiation_loss = np.zeros(100, dtype=dtype_radiation_loss)
    data_radiation_loss[Names.ELEMENT] = range(1, 101)
    list_NC = []
    list_BD = []
    NMAX_array = np.zeros(100, "i")
    with open(path) as fin:
        for i in range(1, 101):
            NMAX, LKMAX = fin.readline().split()
            NMAX, LKMAX = int(NMAX), int(LKMAX)
            assert (LKMAX == LKIN)
            temp1 = []
            while True:
                temp1 += fin.readline().split()
                if len(temp1) == NMAX:
                    break
            temp2 = []
            while True:
                temp2 += fin.readline().split()
                if len(temp2) == NMAX:
                    break
            temp3 = []
            while True:
                temp3 += fin.readline().split()
                if len(temp3) == LKMAX:
                    break
            data_radiation_loss[i - 1][Names.STOPPING_POWER_RADIATIVE] = np.fromiter(map(lambda x: float(x), temp3),
                                                                                     "d")
            NC = np.fromiter(map(lambda x: float(x), temp1), "d")
            BD = np.fromiter(map(lambda x: float(x), temp2), "d")
            list_NC.append(NC)
            list_BD.append(BD)
            NMAX_array[i - 1] = NMAX
    return NMAX_array, list_NC, list_BD, data_radiation_loss


def convert_fcomp_to_table(path):
    KMAX = 279

    dtype_parametrs = np.dtype(
        [
            ("id", "i"),
            ("material", "S72"),
            ("number_of_components", "i"),
            ("zag", "d"),
            (Names.IONISATION_POTENTIAL, "d"),
            (Names.DENSITY, "d"),
        ]
    )
    dtype_composition = np.dtype(
        [
            (Names.ELEMENT, "i"),
            (Names.FRACTION, "d")
        ]
    )
    data_parameters = np.zeros(KMAX, dtype=dtype_parametrs)
    data_composition = []
    with open(path) as fin:
        count = 0
        while count < 279:
            material_id = count + 1
            material = fin.readline()[:-1].strip().replace(" ", "_")
            parameters = fin.readline()[:-1].split()
            number_of_components = int(parameters[0])
            compaund = []
            while True:
                compaund += fin.readline()[:-1].split()
                if len(compaund) == 2 * number_of_components:
                    break
            # print("Compaund: " + str(compaund[1::2]))
            data_parameters["id"][count] = material_id
            data_parameters["material"][count] = material
            data_parameters["number_of_components"][count] = number_of_components
            data_parameters[Names.ZAG][count] = float(parameters[1])
            data_parameters[Names.IONISATION_POTENTIAL][count] = float(parameters[2])
            data_parameters[Names.DENSITY][count] = float(parameters[3])
            data_composition.append(
                (
                    material_id,
                    np.array([(int(z), float(w)) for z, w in zip(compaund[::2], compaund[1::2])],
                             dtype=dtype_composition)
                )
            )
            count += 1
        return data_parameters, data_composition


class ProtonNames:
    ELECTRONIC_STOPPING_POWER = "electronic_stopping_power"
    NUCLEAR_STOPPING_POWER = "nuclear_stopping_power"
    TOTAL_STOPPING_POWER = "total_stopping_power"
    DETOUR_FACTOR = "detour_factor"
    PROJECTED_RANGE = "projected_range"


def convert_fprot(path, NMAX):
    result = []
    dtype = np.dtype(
        [
            (ProtonNames.ELECTRONIC_STOPPING_POWER, "d"),
            (ProtonNames.NUCLEAR_STOPPING_POWER, "d"),
            (Names.CSDA_RANGE, "d"),
            (ProtonNames.DETOUR_FACTOR, "d")
        ]
    )
    with open(path, newline="") as fin:
        data = list(map(lambda x: float(x), fin.read().split()))
    DATA_LENGTH = NMAX
    NUMBER_OF_COLS = 4
    start, stop = 0, DATA_LENGTH
    for material in ProtonMaterials:
        array_list = []
        for i in range(NUMBER_OF_COLS):
            array_list.append(data[start:stop])
            start = stop
            stop += DATA_LENGTH
        result.append(
            (
                material,
                np.core.records.fromarrays(array_list, dtype=dtype)
            )
        )
    return result


def get_mat_name(mat_id):
    return "M" + str(mat_id).rjust(3, "0")


def get_z_name(Z):
    return "Z" + str(Z).rjust(3, "0")


def convert_star_to_hdf5(path):
    with tables.open_file(NIST_STAR_HDF5_PATH, "w") as h5file:
        filters = tables.Filters(complevel=3, fletcher32=True)

        # FCOMP data (materials data)
        estar_path = os.path.join(path, "ESTAR")
        data_parameters, data_composition = convert_fcomp_to_table(os.path.join(estar_path, "FCOMP"))
        table = h5file.create_table(h5file.root, "material_parameters", obj=data_parameters, filters=filters)
        table.flush()
        group = h5file.create_group(h5file.root, "composition")
        for mat_id, data in data_composition:
            table = h5file.create_table(group, get_mat_name(mat_id), obj=data, filters=filters)

        # ESTAR data
        NMAX, list_NC, list_BD, data_radiation_loss = convert_fedat_to_table(
            os.path.join(estar_path, "FEDAT"))
        electron_group = h5file.create_group(h5file.root, "electrons", title="Data for ESTAR")
        table = h5file.create_table(electron_group, "radiation_loss", obj=data_radiation_loss, filters=filters)
        group_NC = h5file.create_group(electron_group, "NC")
        group_BD = h5file.create_group(electron_group, "BD")
        for i in range(0, 100):
            Z = i + 1
            name = get_z_name(Z)
            array = h5file.create_array(group_NC, name, obj=list_NC[i])
            array.flush()
            array = h5file.create_array(group_BD, name, obj=list_BD[i])
            array.flush()
        NMAX / 133
        # PSTAR snd ASTAR DATA
        parametrs = (
            ("PSTAR","FPROT","protons" , "ENG.PRO", 133),
            ("ASTAR","FALPH","helium_ions" , "ENG.ALF", 122))
        for param in parametrs:
            star_dir, data_file, particle, energy_file, NMAX = param
            pstar_path = os.path.join(path, star_dir)
            data_fprot = convert_fprot(os.path.join(pstar_path, data_file), NMAX)
            pstar_group = h5file.create_group(h5file.root, particle, title="Data for {}".format(star_dir))
            for material, data in data_fprot:
                table = h5file.create_table(pstar_group, material.name, filters=filters, obj = data, title="{} stopping-power table for {}".format(particle, material.name))
                table.attrs["cdsa_range_unit"] = "g/cm2"
                table.attrs[ProtonNames.ELECTRONIC_STOPPING_POWER + "_unit"] = "MeV cm2/g"
                table.attrs[ProtonNames.NUCLEAR_STOPPING_POWER + "_unit"] = "MeV cm2/g"
                table.flush()

            with open(os.path.join(pstar_path, energy_file), newline="") as fin:
                energy = map(lambda x: float(x), fin.read().split()[1:])
                energy = np.fromiter(energy, "d")
                array = h5file.create_array(pstar_group, "energy", obj=energy, title="Default energy of {} (corresponds to data in another table)".format(particle))
                array.attrs["unit"] = "MeV"
                array.flush()


if __name__ == '__main__':
    path = sys.argv[1]
    print(path)
    convert_star_to_hdf5(path)
    # print(convert_fprot(path))