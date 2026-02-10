import os
import random
import pickle
import pyard
import json
from runfile import run_impute
from reduce_loci import reduce_loci
from py_graph_imputation_mlo.run_grim import change_donor_file
from py_graph_imputation_mlo.filter_by_rest import change_output_by_extra_gl

# define the input and output directories as constants
INPUT_DIR = "./input_dir"
OUTPUT_DIR = "./output_dir"

# define the paths to the configuration file and the grim graph pickle file as constants
# Get from environment variable so you don't have to edit this file
PATH_TO_CONFIG = "./conf/conf.json"
print(f"Using configuration file: {PATH_TO_CONFIG}")

PATH_TO_GRIM_GRAPH = "./data/graph.pkl"
print(f"Using graph file: {PATH_TO_GRIM_GRAPH}")

# create the input and output directories if they don't already exist
os.makedirs(INPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# define lambdas to convert file numbers to their corresponding input/haplotype/genotype filenames
TO_INPUT_FILE = lambda num: f"input{num}.csv"
TO_HAPLOTYPE_FILE = lambda num: f"haplotype{num}.csv"
TO_GENOTYPE_FILE = lambda num: f"genotype{num}.csv"

with open(f"./data/freqs_dicts/all_freqs.pickle", "rb") as f:
    all_freqs = pickle.load(f)

# open the grim graph pickle file and load its contents into the grim_graph variable
with open(PATH_TO_GRIM_GRAPH, "rb") as f:
    grim_graph = pickle.load(f)

# initialize the pyard object
ard = pyard.init()

# define a mapping from allele names to their corresponding serology types
cast = {
    "A1": "A",
    "A2": "A",
    "B1": "B",
    "B2": "B",
    "C1": "C",
    "C2": "C",
    "DQB1_1": "DQB1",
    "DQB1_2": "DQB1",
    "DRB1_1": "DRB1",
    "DRB1_2": "DRB1",
    "DRB3_1": "DRB3",
    "DRB3_2": "DRB3",
    "DRB4_1": "DRB4",
    "DRB4_2": "DRB4",
    "DRB5_1": "DRB5",
    "DRB5_2": "DRB5",
    "DQA1_1": "DQA1",
    "DQA1_2": "DQA1",
    "DPA1_1": "DPA1",
    "DPA1_2": "DPA1",
    "DPB1_1": "DPB1",
    "DPB1_2": "DPB1",
    "DRBX_1": "DRBX",
    "DRBX_2": "DRBX",
}


def apply_ard_on_file(path):
    """
    Apply py-ard on the genos in the given file and write the output to the input file.

    Parameters:
    path (str): The path to the input file

    Raises:
    Exception: If the input file format is wrong, the function will raise an exception with the error message.

    Returns:
    None
    """
    new_lines = []
    with open(path) as f:
        lines = f.readlines()
        try:
            for line in lines:
                id, glstring, race1, race2 = line.split(",")
                ardstring = ard.redux(glstring, 'lgx')
                new_lines.append(f"{id},{ardstring},{race1},{race2}")
        except Exception as e:
            # raise Exception("Input file format wrong.", str(e))
            raise Exception(str(e))

    with open(path, "w") as f:
        f.writelines(new_lines)


def get_allele_type(allele: str):
    """
        Returns the type of the given allele string.

        Parameters:
        allele (str): The allele string.

        Returns:
        str: The type of the allele, or None if the type is not recognized.
        """
    if allele.startswith("A"):
        return "A"
    elif allele.startswith("B"):
        return "B"
    elif allele.startswith("C"):
        return "C"
    elif allele.startswith("DQB"):
        return "DQB"
    elif allele.startswith("DQA"):
        return "DQA"
    elif allele.startswith("DPA"):
        return "DPA"
    elif allele.startswith("DPB"):
        return "DPB"
    elif allele.startswith("DRX"):
        return "DRX"
    elif allele.startswith("DRB1"):
        return "DRB1"
    elif allele.startswith("DRB3"):
        return "DRB3"
    elif allele.startswith("DRB4"):
        return "DRB4"
    elif allele.startswith("DRB5"):
        return "DRB5"
    elif allele.startswith("DRBX"):
        return "DRBX"
    return None


def convert_allele_val(allele: str, val: str, is_genetic=True):
    """
        Converts the given allele and its value to the appropriate format for use in building a GL string.

        Parameters:
        allele (str): The allele string.
        val (str): The value string.
        is_genetic (bool): True if the value represents a genetic allele, False if it represents a serological allele.

        Returns:
        tuple[str, bool]: A tuple containing the converted string and a boolean indicating whether the conversion was successful.
        """
    if not val or not allele:
        return ""
    if ":" in val:
        return val
    if is_genetic:
        return f"{val}:XX"

    return f"{allele}{val}"  # serology


def build_glstring(alleles):
    """
        Builds a GL string from the given list of alleles.

        Parameters:
        alleles (List[str]): The list of alleles.

        Returns:
        str: The GL string.
        """
    alleles_by_type = {}
    for allele in alleles:
        atype = get_allele_type(allele)
        if atype in alleles_by_type:
            alleles_by_type[atype].append(allele)
        else:
            alleles_by_type[atype] = [allele]

    for atype, lst in alleles_by_type.items():
        # If there is only one allele from a locuse, duplicate it.
        if len(lst) == 1:
            lst = [lst[0], lst[0]]

        alleles_by_type[atype] = '+'.join(lst)

    return '^'.join(alleles_by_type.values())


def apply_ard(alleles: dict, is_genetic: bool):
    """
        Applies py-ard to the given alleles and returns the resulting GL string and ARD string.

        Parameters:
        alleles (dict): A dictionary of alleles and their values.
        is_genetic (bool): True if the values represent genetic alleles, False if they represent serological alleles.

        Returns:
        tuple[str, str]: A tuple containing the GL string and ARD string.
        """
    my_all_alleles = None
    if "glstring" in alleles and alleles["glstring"]:
        my_all_alleles = alleles["glstring"].replace("^", "+").split("+")
    else:
        my_all_alleles = [convert_allele_val(cast.get(allele, None), val, is_genetic) for allele, val in alleles.items()
                          if allele in cast and val]

    glstring = build_glstring(my_all_alleles)
    #print(ard.redux(glstring, 'lgx'))
    return glstring, ard.redux(glstring, 'lgx')


def get_random_number():
    """
        Generates a random number.

        Returns:
        int: The random number.
        """
    while True:
        num = random.randint(0, 10000)
        file_path = os.path.join(INPUT_DIR, TO_INPUT_FILE(num))
        if not os.path.exists(file_path):
            return num


def string_to_file(ard_string, race1="UNK", race2="UNK"):
    """
        Writes the given ARD string to a file with a random name and returns the name.

        Parameters:
        ard_string (str): The ARD string.
        race1 (str): The race of the patient.
        race2 (str): The race of the donor.

        Returns:
        int: The random number that was used as part of the filename.
        """
    line = f"D1,{ard_string},{race1},{race2}"

    while True:
        num = random.randint(0, 10000)
        file_path = os.path.join(INPUT_DIR, TO_INPUT_FILE(num))
        if not os.path.exists(file_path):
            with open(file_path, "w") as f:
                f.write(line)

            return num


def apply_grim_file(file):
    """
        Applies imputation to the given file using GRIM and returns the paths to the resulting genotype and haplotype files.

        Parameters:
        file (FileStorage): The file to apply imputation to.

        Returns:
        tuple[str, str]: A tuple containing the paths to the genotype and haplotype files.
        """
    num = get_random_number()
    input_path = os.path.join(INPUT_DIR, TO_INPUT_FILE(num))
    file.save(input_path)
    apply_ard_on_file(input_path)

    haplotype_path = os.path.join(OUTPUT_DIR, TO_HAPLOTYPE_FILE(num))
    genotype_path = os.path.join(OUTPUT_DIR, TO_GENOTYPE_FILE(num))

    run_impute(PATH_TO_CONFIG, grim_graph, input_path,
               output_haplotype_path=haplotype_path, output_genotype_path=genotype_path)
    os.remove(input_path)

    return genotype_path, haplotype_path


def apply_grim(alleles: dict, race, loci, is_genetic=True):
    """
        Applies py-ard and imputation to the given alleles and returns the resulting genotypes, haplotypes, GL string, and ARD string.

        Parameters:
        alleles (dict): A dictionary of alleles and their values.
        race (str): The race of the patient.
        is_genetic (bool): True if the values represent genetic alleles, False if they represent serological alleles.

        Returns:
        tuple[List[Tuple[str, str]], List[str], List[Tuple[str, str]], str, str]: A tuple containing the genotypes, haplotypes, GL string, and ARD string.
        """
    hap_pop_pair, dominant3 = True, True
    with open(PATH_TO_CONFIG, 'r') as f:
        config = json.load(f)
    glstring, ard_string = apply_ard(alleles, is_genetic)
    num = string_to_file(ard_string, race1=race, race2=race)

    input_path = os.path.join(INPUT_DIR, TO_INPUT_FILE(num))
    haplotype_path = os.path.join(OUTPUT_DIR, TO_HAPLOTYPE_FILE(num))
    genotype_path = os.path.join(OUTPUT_DIR, TO_GENOTYPE_FILE(num))

    gls, lines = change_donor_file(input_path)
    # Apply grim here
    run_impute(PATH_TO_CONFIG, grim_graph, input_path,
               output_haplotype_path=haplotype_path, output_genotype_path=genotype_path, hap_pop_pair=hap_pop_pair)


    path_pmug = haplotype_path
    path_umug = genotype_path
    path_umug_pops = os.path.join(
        config["imputation_out_path"], config["imputation_out_umug_pops_filename"]
    )
    path_pmug_pops = os.path.join(
        config["imputation_out_path"], config["imputation_out_hap_pops_filename"]
    )
    path_miss = os.path.join(
        config["imputation_out_path"], config["imputation_out_miss_filename"]
    )

    change_output_by_extra_gl(
        config, gls, path_pmug, path_umug, path_umug_pops, path_pmug_pops, path_miss
    )  # filter reasults in our origianl file, add miss to existing miss

    # changing to original donor file
    with open(input_path, "w") as file:
        for line in lines:
            file.write(line)
    file.close()

    os.remove(input_path) #?#

    output_file_hap, output_file_muug = reduce_loci(loci, genotype_path, haplotype_path)
    genotypes = read_genos(output_file_muug)
    haplotypes, haplotypes_pairs = read_haps(output_file_hap, race)
    return genotypes, haplotypes, haplotypes_pairs, glstring, ard_string


def read_genos(path):
    """
        Reads in HLA genotypes and probabilities from a file and returns them as a list of tuples.

        Parameters:
        path (str): Path to the input file.

        Returns:
        list: A list of tuples where each tuple contains the following four values:
              1. Index (str): The index of the HLA genotype.
              2. HLA genotype (str): The HLA genotype in the format 'A*01:01+B*08:01'.
              3. Probability (str): The probability of the HLA genotype in scientific notation with 2 decimal places.
              4. Index (str): The index of the HLA genotype (same as the first value).

        Example:
           >> read_genos('genos.txt')
        [('1', 'A*01:01+B*08:01', '1.23e-02', '1'), ('2', 'A*02:01+B*07:02', '2.34e-03', '2')]
        """
    hlas_and_probs = []

    # Of course use num here to get the results of this patient
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            _, hla, prob, index = line.split(",")
            hlas_and_probs.append([index, hla, f"{float(prob):0.2e}"])
            # hlas_and_probs[hla] = f"{float(prob):0.2e}"

    probs = [float(prob[2]) for prob in hlas_and_probs]
    #min_prob = min(probs)
    #max_prob = max(probs)
    sum_probs = sum(probs)

    for i, prob in enumerate(probs):
        normalized = prob / sum_probs
        #normalized = (prob - min_prob) / (max_prob - min_prob)
        hlas_and_probs[i].append(f"{normalized:.3f}")

    # os.remove(path)
    return hlas_and_probs


def read_haps(path, race_str):
    """
        Reads in HLA haplotypes and probabilities from a file and returns them as a tuple of two lists.

        Parameters:
        path (str): Path to the input file.

        Returns:
        tuple: A tuple of two lists where the first list contains dictionaries for each HLA haplotype and their probabilities for each race,
               and the second list contains tuples with the following four values:
                1. Index (str): The index of the HLA haplotype.
                2. First HLA allele (str): The first HLA allele in the haplotype, in the format 'A*01:01'.
                3. Second HLA allele (str): The second HLA allele in the haplotype, in the format 'B*08:01'.
                4. Probability (str): The probability of the HLA haplotype in scientific notation with 2 decimal places.

        Example:
            >> read_haps('haps.txt')
           ({'A*01:01': {'Asian': '1.23e-02', 'African': '3.45e-03'}, 'B*08:01': {'Asian': '5.67e-04', 'African': '7.89e-05'}},
            [('1', 'A*01:01', 'B*08:01', '1.23e-02'), ('2', 'A*01:01', 'B*08:02', '2.34e-03')])
        """
    print(grim_graph)
    hlas = set()
    hla_and_probs = {}
    hla_pairs = []
    race_list = race_str.split(";")

    # Of course use num here to get the results of this patient
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            _, hla, prob, index = line.split(",")
            hla1, hla2 = hla.split("+")
            hlas.add(hla1)
            hlas.add(hla2)
            hla_pairs.append([index, hla1, hla2, f"{float(prob):0.2e}"])

    haplos = [float(hap[3]) for hap in hla_pairs]
    sum_haplos = sum(haplos)

    for i, prob in enumerate(haplos):
        normalized = prob / sum_haplos
        hla_pairs[i].append(f"{normalized:.3f}")

    conf = dict()
    with open(PATH_TO_CONFIG, "r") as f:
        conf = json.load(f)

    populations = conf["populations"]

    """for hla in hlas:
        hla_and_probs[hla] = {}
        for race, dict_ in all_freqs.items():
            if hla in dict_:
                hla_and_probs[hla][race] = f"{dict_[hla]:0.2e}"
            else:
                freq = grim_graph.whole_graph.nodes[hla]["freq"][populations.index(race)]
                hla_and_probs[hla][race] = f"{freq:0.3e}" """
    for hla in hlas:
        hla_and_probs[hla] = {}
        for race, _ in all_freqs.items():
            if race in populations and hla in grim_graph.Vertices_attributes:
                freq = grim_graph.Vertices_attributes[hla][1][populations.index(race)]
                #freq = [str(round(f, 3)) for f in freq]
                hla_and_probs[hla][race] = f"{freq:0.8e}"
            elif race in populations and hla in grim_graph.Whole_Vertices_attributes:
                freq = grim_graph.Whole_Vertices_attributes[hla][1][populations.index(race)]
                hla_and_probs[hla][race] = f"{freq:0.8e}"
                #freq = [str(round(f,3)) for f in freq]
                #hla_and_probs[hla][race] = f"{freq:0.8e}"
                #hla_and_probs[hla][race] = populations.index(race)
            else:
                hla_and_probs[hla][race] = 0

                # os.remove(path)
    return hla_and_probs, hla_pairs

