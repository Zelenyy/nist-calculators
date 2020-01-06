import os
import sys

import numpy as np
import tables

ROOT_PATH = os.path.dirname(__file__)

class Names:
    IONISATION_POTENTIAL = "ionisation_potential"
    DENSITY = "density"
    ZAG = "zag"
    ENERGY = "energy"
    ELEMENT = "element"
    FRACTION = "fraction"
    STOPPING_POWER_COLLISION_DELTA =  "stopping_power_collision_delta"
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
    data_radiation_loss[Names.ELEMENT] = range(1,101)
    list_NC = []
    list_BD = []
    NMAX_array = np.zeros(100, "i")
    with open(path) as fin:
        for i in range(1,101):
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
            data_radiation_loss[i-1][Names.STOPPING_POWER_RADIATIVE] = np.fromiter(map(lambda x: float(x), temp3), "d")
            NC = np.fromiter(map(lambda x: float(x), temp1), "d")
            BD = np.fromiter(map(lambda x: float(x), temp2), "d")
            list_NC.append(NC)
            list_BD.append(BD)
            NMAX_array[i-1] = NMAX
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
                if len(compaund) == 2*number_of_components:
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
                    np.array([(int(z), float(w)) for z, w in zip(compaund[::2], compaund[1::2])], dtype=dtype_composition)
                )
            )
            count += 1
        return data_parameters, data_composition


def get_mat_name(mat_id):
    return  "M" + str(mat_id).rjust(3, "0")

def get_z_name(Z):
    return "Z" + str(Z).rjust(3, "0")

def convert_star_to_hdf5(path):
    data_parameters, data_composition = convert_fcomp_to_table(os.path.join(path,"FCOMP"))
    NMAX, list_NC, list_BD, data_radiation_loss = convert_fedat_to_table(
        os.path.join(path, "FEDAT"))


    with tables.open_file(os.path.join(ROOT_PATH, 'data', 'NIST_ESTAR.hdf5'), "w") as h5file:



        filters = tables.Filters(complevel=3, fletcher32=True)

        table = h5file.create_table(h5file.root, "material_parameters", obj=data_parameters, filters=filters)
        table.flush()
        group = h5file.create_group(h5file.root, "composition")
        for mat_id, data in data_composition:
            table = h5file.create_table(group, get_mat_name(mat_id), obj=data, filters=filters)

        table = h5file.create_table(h5file.root, "radiation_loss", obj=data_radiation_loss, filters=filters)
        group_NC = h5file.create_group(h5file.root, "NC")
        group_BD = h5file.create_group(h5file.root, "BD")
        for i in range(0,100):
            Z = i + 1
            name = get_z_name(Z)
            array = h5file.create_array(group_NC, name, obj=list_NC[i])
            array.flush()
            array = h5file.create_array(group_BD, name, obj=list_BD[i])
            array.flush()




if __name__ == '__main__':
    path = sys.argv[1]
    print(path)
    convert_star_to_hdf5(path)
