from star.electron import PredefinedMaterials

DATA_IND = [1,  2,  4,  6,  7,  8, 10, 13, 14, 18, 22, 26, 29, 32,
               36, 42, 47, 50, 54, 64, 74, 78, 79, 82, 92, 99,101,103,
              104,106,111,119,120,126,130,134,138,139,141,155,160,169,
              179,185,189,191,197,200,201,202,203,204,209,213,215,216,
 219,221,222,223,225,226,227,232,238,245,252,255,263,264,
266,276,277,279]


def generate_proton_materials():
    text = "class PredefinedMaterials(Enum):\n"
    for ind in DATA_IND:
        materials = PredefinedMaterials(ind)
        text += "    {} = {}\n".format(materials.name, "auto()")
    return text

if __name__ == '__main__':
    print(generate_proton_materials())
