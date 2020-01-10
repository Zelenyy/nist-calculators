import logging
import math
from dataclasses import dataclass
from enum import auto, Enum
from typing import List, Union, Tuple

import numpy as np
import tables
from scipy.interpolate import CubicSpline

from star._data_converter import Names, get_mat_name, get_z_name, NIST_STAR_HDF5_PATH

DATA_ATB = np.array([
    1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.011,
    14.00674, 15.9994, 18.9984032, 20.1797, 22.989768,
    24.3050, 26.981539, 28.0855, 30.973762, 32.066, 35.4527,
    39.948, 39.0983, 40.078, 44.955910, 47.88, 50.9415,
    51.9961, 54.93805, 55.847, 58.93320, 58.69, 63.546,
    65.39, 69.723, 72.61, 74.92159, 78.96, 79.904, 83.80,
    85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.94,
    97.9072, 101.07, 102.9055, 106.42, 107.8682, 112.411,
    114.82, 118.710, 121.75, 127.60, 126.90447, 131.29,
    132.90543, 137.327, 138.9055, 140.115, 140.90765, 144.24,
    144.9127, 150.36, 151.965, 157.25, 158.92534, 162.50,
    164.93032, 167.26, 168.93421, 173.04, 174.967, 178.49,
    180.9479, 183.85, 186.207, 190.2, 192.22, 195.08,
    196.96654, 200.59, 204.3833, 207.2, 208.98037, 208.9824,
    209.9871, 222.0176, 223.0197, 226.0254, 227.0278,
    232.0381, 231.03588, 238.0289, 237.0482, 239.0522,
    243.0614, 247.0703, 247.0703, 251.0796, 252.083,
    257.0951
])

DATA_ER = np.array([1.00E-03, 1.25E-03, 1.50E-03, 1.75E-03, 2.00E-03, 2.50E-03,
                    3.00E-03, 3.50E-03, 4.00E-03, 4.50E-03, 5.00E-03, 5.50E-03,
                    6.00E-03, 7.00E-03, 8.00E-03, 9.00E-03, 1.00E-02, 1.25E-02,
                    1.50E-02, 1.75E-02, 2.00E-02, 2.50E-02, 3.00E-02, 3.50E-02,
                    4.00E-02, 4.50E-02, 5.00E-02, 5.50E-02, 6.00E-02, 7.00E-02,
                    8.00E-02, 9.00E-02, 1.00E-01, 1.25E-01, 1.50E-01, 1.75E-01,
                    2.00E-01, 2.50E-01, 3.00E-01, 3.50E-01, 4.00E-01, 4.50E-01,
                    5.00E-01, 5.50E-01, 6.00E-01, 7.00E-01, 8.00E-01, 9.00E-01,
                    1.00E+00, 1.25E+00, 1.50E+00, 1.75E+00, 2.00E+00, 2.50E+00,
                    3.00E+00, 3.50E+00, 4.00E+00, 4.50E+00, 5.00E+00, 5.50E+00,
                    6.00E+00, 7.00E+00, 8.00E+00, 9.00E+00, 1.00E+01, 1.25E+01,
                    1.50E+01, 1.75E+01, 2.00E+01, 2.50E+01, 3.00E+01, 3.50E+01,
                    4.00E+01, 4.50E+01, 5.00E+01, 5.50E+01, 6.00E+01, 7.00E+01,
                    8.00E+01, 9.00E+01, 1.00E+02, 1.25E+02, 1.50E+02, 1.75E+02,
                    2.00E+02, 2.50E+02, 3.00E+02, 3.50E+02, 4.00E+02, 4.50E+02,
                    5.00E+02, 5.50E+02, 6.00E+02, 7.00E+02, 8.00E+02, 9.00E+02,
                    1.00E+03, 1.25E+03, 1.50E+03, 1.75E+03, 2.00E+03, 2.50E+03,
                    3.00E+03, 3.50E+03, 4.00E+03, 4.50E+03, 5.00E+03, 5.50E+03,
                    6.00E+03, 7.00E+03, 8.00E+03, 9.00E+03, 1.00E+04])


class PredefinedMaterials(Enum):
    HYDROGEN = auto()
    HELIUM = auto()
    LITHIUM = auto()
    BERYLLIUM = auto()
    BORON = auto()
    AMORPHOUS_CARBON = auto()  # (density 2.0 g / cm3)
    NITROGEN = auto()
    OXYGEN = auto()
    FLUORINE = auto()
    NEON = auto()
    SODIUM = auto()
    MAGNESIUM = auto()
    ALUMINUM = auto()
    SILICON = auto()
    PHOSPHORUS = auto()
    SULFUR = auto()
    CHLORINE = auto()
    ARGON = auto()
    POTASSIUM = auto()
    CALCIUM = auto()
    SCANDIUM = auto()
    TITANIUM = auto()
    VANADIUM = auto()
    CHROMIUM = auto()
    MANGANESE = auto()
    IRON = auto()
    COBALT = auto()
    NICKEL = auto()
    COPPER = auto()
    ZINC = auto()
    GALLIUM = auto()
    GERMANIUM = auto()
    ARSENIC = auto()
    SELENIUM = auto()
    BROMINE = auto()
    KRYPTON = auto()
    RUBIDIUM = auto()
    STRONTIUM = auto()
    YTTRIUM = auto()
    ZIRCONIUM = auto()
    NIOBIUM = auto()
    MOLYBDENUM = auto()
    TECHNETIUM = auto()
    RUTHENIUM = auto()
    RHODIUM = auto()
    PALLADIUM = auto()
    SILVER = auto()
    CADMIUM = auto()
    INDIUM = auto()
    TIN = auto()
    ANTIMONY = auto()
    TELLURIUM = auto()
    IODINE = auto()
    XENON = auto()
    CESIUM = auto()
    BARIUM = auto()
    LANTHANUM = auto()
    CERIUM = auto()
    PRASEODYMIUM = auto()
    NEODYMIUM = auto()
    PROMETHIUM = auto()
    SAMARIUM = auto()
    EUROPIUM = auto()
    GADOLINIUM = auto()
    TERBIUM = auto()
    DYSPROSIUM = auto()
    HOLMIUM = auto()
    ERBIUM = auto()
    THULIUM = auto()
    YTTERBIUM = auto()
    LUTETIUM = auto()
    HAFNIUM = auto()
    TANTALUM = auto()
    TUNGSTEN = auto()
    RHENIUM = auto()
    OSMIUM = auto()
    IRIDIUM = auto()
    PLATINUM = auto()
    GOLD = auto()
    MERCURY = auto()
    THALLIUM = auto()
    LEAD = auto()
    BISMUTH = auto()
    POLONIUM = auto()
    ASTATINE = auto()
    RADON = auto()
    FRANCIUM = auto()
    RADIUM = auto()
    ACTINIUM = auto()
    THORIUM = auto()
    PROTACTINIUM = auto()
    URANIUM = auto()
    NEPTUNIUM = auto()
    PLUTONIUM = auto()
    AMERICIUM = auto()
    CURIUM = auto()
    BERKELIUM = auto()
    CALIFORNIUM = auto()
    A_150_TISSUE_EQUIVALENT_PLASTIC = auto()
    ACETONE = auto()
    ACETYLENE = auto()
    ADENINE = auto()
    ADIPOSE_TISSUE_ICRP = auto()
    AIR_DRY_NEAR_SEA_LEVEL = auto()
    ALANINE = auto()
    ALUMINUM_OXIDE = auto()
    AMBER = auto()
    AMMONIA = auto()
    ANILINE = auto()
    ANTHRACENE = auto()
    B_100_BONE_EQUIVALENT_PLASTIC = auto()
    BAKELITE = auto()
    BARIUM_FLUORIDE = auto()
    BARIUM_SULFATE = auto()
    BENZENE = auto()
    BERYLLIUM_OXIDE = auto()
    BISMUTH_GERMANIUM_OXIDE = auto()
    BLOOD_ICRP = auto()
    BONE_COMPACT_ICRU = auto()
    BONE_CORTICAL_ICRP = auto()
    BORON_CARBIDE = auto()
    BORON_OXIDE = auto()
    BRAIN_ICRP = auto()
    BUTANE = auto()
    N_BUTYL_ALCOHOL = auto()
    C_552_AIR_EQUIVALENT_PLASTIC = auto()
    CADMIUM_TELLURIDE = auto()
    CADMIUM_TUNGSTATE = auto()
    CALCIUM_CARBONATE = auto()
    CALCIUM_FLUORIDE = auto()
    CALCIUM_OXIDE = auto()
    CALCIUM_SULFATE = auto()
    CALCIUM_TUNGSTATE = auto()
    CARBON_DIOXIDE = auto()
    CARBON_TETRACHLORIDE = auto()
    CELLULOSE_ACETATE_CELLOPHANE = auto()
    CELLULOSE_ACETATE_BUTYRATE = auto()
    CELLULOSE_NITRATE = auto()
    CERIC_SULFATE_DOSIMETER_SOLUTION = auto()
    CESIUM_FLUORIDE = auto()
    CESIUM_IODIDE = auto()
    CHLOROBENZENE = auto()
    CHLOROFORM = auto()
    CONCRETE_PORTLAND = auto()
    CYCLOHEXANE = auto()
    f1_2_DICHLOROBENZENE = auto()
    DICHLORODIETHYL_ETHER = auto()
    f1_2_DICHLOROETHANE = auto()
    DIETHYL_ETHER = auto()
    N_N_DIMETHYL_FORMAMIDE = auto()
    DIMETHYL_SULFOXIDE = auto()
    ETHANE = auto()
    ETHYL_ALCOHOL = auto()
    ETHYL_CELLULOSE = auto()
    ETHYLENE = auto()
    EYE_LENS_ICRP = auto()
    FERRIC_OXIDE = auto()
    FERROBORIDE = auto()
    FERROUS_OXIDE = auto()
    FERROUS_SULFATE_DOSIMETER_SOLUTION = auto()
    FREON_12 = auto()
    FREON_12_B2 = auto()
    FREON_13 = auto()
    FREON_13_B1 = auto()
    FREON_13_I1 = auto()
    GADOLINIUM_OXYSULFIDE = auto()
    GALLIUM_ARSENIDE = auto()
    GEL_IN_PHOTOGRAPHIC_EMULSION = auto()
    Pyrex_Glass = auto()
    GLASS_LEAD = auto()
    GLASS_PLATE = auto()
    GLUCOSE = auto()
    GLUTAMINE = auto()
    GLYCEROL = auto()
    GUANINE = auto()
    GYPSUM_PLASTER_OF_PARIS = auto()
    N_HEPTANE = auto()
    N_HEXANE = auto()
    KAPTON_POLYIMIDE_FILM = auto()
    LANTHANUM_OXYBROMIDE = auto()
    LANTHANUM_OXYSULFIDE = auto()
    LEAD_OXIDE = auto()
    LITHIUM_AMIDE = auto()
    LITHIUM_CARBONATE = auto()
    LITHIUM_FLUORIDE = auto()
    LITHIUM_HYDRIDE = auto()
    LITHIUM_IODIDE = auto()
    LITHIUM_OXIDE = auto()
    LITHIUM_TETRABORATE = auto()
    LUNG_ICRP = auto()
    M3_WAX = auto()
    MAGNESIUM_CARBONATE = auto()
    MAGNESIUM_FLUORIDE = auto()
    MAGNESIUM_OXIDE = auto()
    MAGNESIUM_TETRABORATE = auto()
    MERCURIC_IODIDE = auto()
    METHANE = auto()
    METHANOL = auto()
    MIX_D_WAX = auto()
    MS20_TISSUE_SUBSTITUTE = auto()
    MUSCLE_SKELETAL_ICRP = auto()
    MUSCLE_STRIATED_ICRU = auto()
    MUSCLE_EQUIVALENT_LIQUID_WITH_SUCROSE = auto()
    MUSCLE_EQUIVALENT_LIQUID_WITHOUT_SUCROSE = auto()
    NAPHTHALENE = auto()
    NITROBENZENE = auto()
    NITROUS_OXIDE = auto()
    NYLON_DU_PONT_ELVAMIDE_8062 = auto()
    NYLON_TYPE_6_AND_TYPE_6_6 = auto()
    NYLON_TYPE_6_10 = auto()
    NYLON_TYPE_11_RILSAN = auto()
    OCTANE_LIQUID = auto()
    PARAFFIN_WAX = auto()
    N_PENTANE = auto()
    PHOTOGRAPHIC_EMULSION = auto()
    PLASTIC_SCINTILLATOR_VINYLTOLUENE_BASED = auto()
    PLUTONIUM_DIOXIDE = auto()
    POLYACRYLONITRILE = auto()
    POLYCARBONATEMAKROLON_LEXAN = auto()
    POLYCHLOROSTYRENE = auto()
    POLYETHYLENE = auto()
    POLYETHYLENE_TEREPHTHALATE_MYLAR = auto()
    POLYMETHYL_METHACRALATE_LUCITE_PERSPEX = auto()
    POLYOXYMETHYLENE = auto()
    POLYPROPYLENE = auto()
    POLYSTYRENE = auto()
    POLYTETRAFLUOROETHYLENE_TEFLON = auto()
    POLYTRIFLUOROCHLOROETHYLENE = auto()
    POLYVINYL_ACETATE = auto()
    POLYVINYL_ALCOHOL = auto()
    POLYVINYL_BUTYRAL = auto()
    POLYVINYL_CHLORIDE = auto()
    POLYVINYLIDENE_CHLORIDE_SARAN = auto()
    POLYVINYLIDENE_FLUORIDE = auto()
    POLYVINYL_PYRROLIDONE = auto()
    POTASSIUM_IODIDE = auto()
    POTASSIUM_OXIDE = auto()
    PROPANE = auto()
    PROPANE_LIQUID = auto()
    N_PROPYL_ALCOHOL = auto()
    PYRIDINE = auto()
    RUBBER_BUTYL = auto()
    RUBBER_NATURAL = auto()
    RUBBER_NEOPRENE = auto()
    SILICON_DIOXIDE = auto()
    SILVER_BROMIDE = auto()
    SILVER_CHLORIDE = auto()
    SILVER_HALIDES_IN_PHOTOGRAPHIC_EMULSION = auto()
    SILVER_IODIDE = auto()
    SKIN_ICRP = auto()
    SODIUM_CARBONATE = auto()
    SODIUM_IODIDE = auto()
    SODIUM_MONOXIDE = auto()
    SODIUM_NITRATE = auto()
    STILBENE = auto()
    SUCROSE = auto()
    TERPHENYL = auto()
    TESTES_ICRP = auto()
    TETRACHLOROETHYLENE = auto()
    THALLIUM_CHLORIDE = auto()
    TISSUE_SOFT_ICRP = auto()
    TISSUE_SOFT_ICRU_FOUR_COMPONENT = auto()
    TISSUE_EQUIVALENT_GAS_METHANE_BASED = auto()
    TISSUE_EQUIVALENT_GAS_PROPANE_BASED = auto()
    TITANIUM_DIOXIDE = auto()
    TOLUENE = auto()
    TRICHLOROETHYLENE = auto()
    TRIETHYL_PHOSPHATE = auto()
    TUNGSTEN_HEXAFLUORIDE = auto()
    URANIUM_DICARBIDE = auto()
    URANIUM_MONOCARBIDE = auto()
    URANIUM_OXIDE = auto()
    UREA = auto()
    VALINE = auto()
    VITON_FLUOROELASTOMER = auto()
    WATER_LIQUID = auto()
    WATER_VAPOR = auto()
    XYLEN = auto()
    GRAPHITE = auto()  # (density 1.7 g / cm3)


@dataclass
class MaterialParameters:
    material_name: str
    number_of_components: int
    zag: float  # Z/A
    ionisation_potential: float  # eV
    density: float  # gram/cm3
    mz: List[int]
    wt: List[float]

    def potl(self):
        return math.log(self.ionisation_potential * 1e-6)


def load_radiation_loss(elements: Union[int, List[int]]) -> Tuple[list, list, np.ndarray]:
    if isinstance(elements, int):
        elements = [elements]
    with tables.open_file(NIST_STAR_HDF5_PATH) as h5file:
        electron_group = tables.Group(h5file.root, "electrons")
        table = h5file.get_node(electron_group, "radiation_loss")
        coord = np.asarray(elements) - 1
        data = table.read_coordinates(coord)
        list_NC = []
        list_BD = []
        group_NC = tables.Group(electron_group, "NC")
        group_BD = tables.Group(electron_group, "BD")
        for element in elements:
            name = get_z_name(element)
            array_NC = h5file.get_node(group_NC, name).read()
            array_BD = h5file.get_node(group_BD, name).read()
            list_NC.append(array_NC)
            list_BD.append(array_BD)
    return list_NC, list_BD, data


def load_material(material: PredefinedMaterials):
    with tables.open_file(NIST_STAR_HDF5_PATH) as h5file:
        table = h5file.get_node(h5file.root, "material_parameters")
        data = table[material.value - 1]
        name = get_mat_name(data["id"])
        composition = h5file.get_node("/composition", name).read()
    return MaterialParameters(
        material_name=data["material"],
        number_of_components=len(composition),
        zag=data[Names.ZAG],
        ionisation_potential=data[Names.IONISATION_POTENTIAL],
        density=data[Names.DENSITY],
        mz=composition[Names.ELEMENT],
        wt=composition[Names.FRACTION]
    )


def calculate_estar_table(material: Union[MaterialParameters, PredefinedMaterials]):
    """
     Stopping powers, ranges and radiation yields for standard energy grid. Return next columns:
        * Kinetic energy, `MeV`
        * Stopping power collision delta, `MeV cm2/g`
        * Stopping power radiative, `MeV cm2/g`
        * Stopping power total, `MeV cm2/g`
        * CSDA range, `g/cm2`
        * Radiation yield,
        * Density effect parameter delta

    Parameters
    ----------
    material
            special class description material or predefined material from `PredefinedMaterials` enumerations
    Returns
    -------
    data : ndarray with all data
    """
    if isinstance(material, PredefinedMaterials):
        material = load_material(material)
        logging.debug("Load predefiend material: {}".format(material))
    else:
        material.mz = np.asarray(material.mz)
        material.wt = np.asarray(material.wt)

    dtype = np.dtype([
        ("energy", "d"),
        ("stopping_power_collision_delta", "d"),
        ("stopping_power_radiative", "d"),
        ("stopping_power_total", "d"),
        ("csda_range", "d"),
        ("radiation_yield", "d"),
        ("density_effect", "d")
    ])

    data_stopping_power = calculate_stopping_power(material)
    data = np.zeros(len(DATA_ER), dtype=dtype)
    for name in data_stopping_power.dtype.names:
        data[name] = data_stopping_power[name]
    energy_log = np.log(DATA_ER)
    cs_rloss = CubicSpline(x=energy_log, y=np.log(data[Names.STOPPING_POWER_RADIATIVE]))
    cs_tloss = CubicSpline(x=energy_log, y=np.log(data[Names.STOPPING_POWER_TOTAL]))
    n = len(DATA_ER)
    RG = np.zeros(n, 'd')
    RAD = np.zeros(n, 'd')
    RG[0] = 0.5 * DATA_ER[0] / data[Names.STOPPING_POWER_TOTAL][0]
    RAD[0] = 0.5 * DATA_ER[0] * data[Names.STOPPING_POWER_RADIATIVE][0] / data[Names.STOPPING_POWER_TOTAL][0]
    MGRD = 21
    EDIFF = np.diff(DATA_ER) / (MGRD - 1)
    DET = EDIFF / 3.0
    for i in range(1, n):
        ETL = np.log(DATA_ER[i] - EDIFF[i - 1] * np.arange(MGRD))
        GRAND = np.exp(-cs_tloss(x=ETL))
        GRAND1 = np.exp(cs_rloss(x=ETL)) * GRAND
        STEP = GRAL(DET[i - 1], GRAND, MGRD)
        DRAD = GRAL(DET[i - 1], GRAND1, MGRD)
        RG[i] = RG[i - 1] + STEP
        RAD[i] = RAD[i - 1] + DRAD
    RAD /= DATA_ER
    data[Names.CSDA_RANGE] = RG
    data[Names.RADIATION_YIELD] = RAD
    # default grid cat head
    data = data[16:]
    return data


def calculate_stopping_power(material: Union[MaterialParameters, PredefinedMaterials], energy=None) -> np.ndarray:
    """
    Stopping powers only, for user-selected energy grid.
    Return next columns:
        * Kinetic energy, `MeV`
        * Stopping power collision delta, `MeV cm2/g`
        * Stopping power radiative, `MeV cm2/g`
        * Stopping power total, `MeV cm2/g`
        * Density effect parameter delta

    Parameters
    ----------
    material
            special class description material or predefined material from `PredefinedMaterials` enumerations
    Returns
    -------
    data : ndarray with all data
    """

    if isinstance(material, PredefinedMaterials):
        material = load_material(material)
        logging.debug("Load predefiend material: {}".format(material))
    else:
        material.mz = np.asarray(material.mz)
        material.wt = np.asarray(material.wt)

    dtype = np.dtype([
        ("energy", "d"),
        ("stopping_power_collision_delta", "d"),
        ("stopping_power_radiative", "d"),
        ("stopping_power_total", "d"),
        (Names.DENSITY_EFFECT, "d")
    ])

    # LKMAX = 113
    COFF = 0.307072
    RMASS = 0.510999906  # electron mass
    NUMQ = 50
    LMAX = 1101
    QBEG = 1.0e-4
    QFAC = 10 ** (1 / NUMQ)
    Q = np.zeros(LMAX, "d")
    Q[0] = QBEG
    for i in range(1, LMAX):
        Q[i] = Q[i - 1] * QFAC

    POTL = math.log(material.ionisation_potential * 1e-6)
    ERL = np.log(DATA_ER)
    if energy is None:
        energy = DATA_ER
    energy_log = np.log(energy)
    n = len(energy)
    data = np.zeros(n, dtype=dtype)
    data[Names.ENERGY] = np.asarray(energy)

    G = material.wt * material.mz / DATA_ATB[material.mz - 1]
    GTOT = np.sum(G)
    ZAV = GTOT
    logging.debug("ZAV = " + str(ZAV))
    HOM = 28.81593 * math.sqrt(material.density * ZAV)
    logging.debug("HOM = " + str(HOM))
    PHIL = 2.0 * math.log(material.ionisation_potential / HOM)
    CBAR = PHIL + 1.0
    G /= GTOT

    list_NC, list_BD, data_radiation_loss = load_radiation_loss(material.mz)
    RLOST = np.sum((material.wt * data_radiation_loss[Names.STOPPING_POWER_RADIATIVE].T), axis=1)

    for nc, bd in zip(list_NC, list_BD):
        if nc[-1] < 0:
            nc[-1] = -nc[-1]
            if material.number_of_components <= 1:
                bd[-1] = 0.0
    sum_ = np.array([np.sum(nc) for nc in list_NC])

    F, EN = [], []
    for i in range(material.number_of_components):
        nc = list_NC[i]
        bd = list_BD[i]
        for j in range(len(nc)):
            F.append(nc[j] * G[i] / sum_[i])
            EN.append(bd[j])
    F = np.asarray(F)
    EN = np.asarray(EN)
    logging.debug("EN = " + str(EN))
    RLOSTL = np.log(RLOST)
    # make spline with coeff ARL,BRL,CRL,DRL
    cs_rlostl = CubicSpline(x=ERL, y=RLOSTL)

    nmax = sum([len(x) for x in list_NC])
    ALF = np.repeat(2.0 / 3.0, nmax)
    if EN[-1] <= 0:
        ALF[-1] = 1.0
    EPS = (EN / HOM) ** 2
    ROOT = 1.0
    while True:
        FUN = -PHIL
        DER = 0.0
        TRM = ROOT * EPS + ALF * F
        FUN += np.sum(F * np.log(TRM))
        DER += np.sum(F * EPS / TRM)
        DROOT = FUN / DER
        ROOT = ROOT - DROOT
        if (abs(DROOT) < 0.00001):
            break
    # FACTOR = math.sqrt(ROOT) # debug variable in f77
    logging.debug("ROOT = " + str(ROOT))
    EPS = ROOT * EPS
    logging.debug("EPS = " + str(EPS))
    YCUT = 0.0

    if EN[-1] > 0:
        YCUT = 1 / np.sum(F / EPS)
    logging.debug("YCUT = " + str(YCUT))
    YQ = np.zeros(LMAX, "d")
    D = np.zeros(LMAX, "d")

    for i in range(0, LMAX):
        sum_ = np.sum(F / (EPS + Q[i]))
        YQ[i] = 1 / sum_
        sum_ = np.sum(F * np.log(1.0 + Q[i] / (EPS + ALF * F)))
        D[i] = sum_ - Q[i] / (YQ[i] + 1.0)
    YQL = np.log(YQ)
    # TCUT = RMASS*(math.sqrt(YCUT +1.0)-1) # debug value in f77
    cs_D = CubicSpline(x=YQL, y=D)

    TAU = energy / RMASS
    Y = TAU * (TAU + 2)
    BETQ = Y / ((TAU + 1.0) ** 2)
    indx = np.logical_and(Y >= YQ[0], Y > YCUT)  # check max energy Y < YQ[-1]
    DELTA = np.zeros(Y.size, "d")
    DELTA[indx] = cs_D(x=np.log(Y)[indx])
    SPART = energy_log - POTL + 0.5 * np.log(1 + 0.5 * TAU) - 0.5 * DELTA
    TERM = (1.0 - BETQ) * (1.0 + (TAU ** 2) / 8.0 - (2.0 * TAU + 1.0) * math.log(2.0))
    STNUM = SPART + 0.5 * TERM

    data[Names.STOPPING_POWER_COLLISION_DELTA] = COFF * ZAV * STNUM / BETQ
    data[Names.STOPPING_POWER_RADIATIVE] = np.exp(cs_rlostl(x=energy_log))
    data[Names.STOPPING_POWER_TOTAL] = data[Names.STOPPING_POWER_RADIATIVE] + data[Names.STOPPING_POWER_COLLISION_DELTA]
    data[Names.DENSITY_EFFECT] = DELTA
    return data


def GRAL(DELTA, G, N):
    NL1 = N - 1
    NL2 = N - 2
    if N % 2 == 1:
        if N - 1 <= 0:
            SIGMA = 0.0
        elif N - 3 <= 0:
            SIGMA = G[0] + 4.0 * G[1] + G[2]
        else:
            SUM4 = np.sum(G[1:NL1 - 1:2])
            SUM2 = np.sum(G[2:NL2 - 1:2])
            SIGMA = G[0] + 4.0 * SUM4 + 2.0 * SUM2 + G[N - 1]
    else:
        if N - 2 <= 0:
            SIGMA = 1.5 * (G[0] + G[1])
        elif N - 4 <= 0:
            SIGMA = 1.125 * (G[0] + 3.0 * G[1] + 3.0 * G[2] + G[3])
        elif N - 6 <= 0:
            SIGMA = G[0] + 3.875 * G[1] + 2.625 * G[2] + 2.625 * G[3] + 3.875 * G[4] + G[5]
        elif N - 8 <= 0:
            SIGMA = G[0] + 3.875 * G[1] + 2.625 * G[2] + 2.625 * G[3] + 3.875 * G[4] + 2 * G[5] + 4 * G[6] + G[7]
        else:
            SIG6 = G[0] + 3.875 * G[1] + 2.625 * G[2] + 2.625 * G[3] + 3.875 * G[4] + G[5]
            SUM4 = np.sum(G[6:NL1 - 1:2])
            SUM2 = np.sum(G[7:NL2 - 1:2])
            SIGMA = SIG6 + G[7] + 4.0 * SUM4 + 2.0 * SUM2 + G[N - 1]
    return DELTA * SIGMA
