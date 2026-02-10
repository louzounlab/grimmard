import json
import pathlib
import sys
import os
import pickle

sys.path.insert(0, os.path.join(".."))

from grim import grim
# grim.graph_freqs(configuration_file)

def run_impute(
        conf_file,
        graph,
        input_path,
        output_haplotype_path,
        output_genotype_path,
        project_dir_graph="",
        project_dir_in_file="",
        hap_pop_pair=False,
):
    configuration_file = conf_file

    # Read configuration file and load properties
    with open(configuration_file) as f:
        json_conf = json.load(f)

    graph_files_path = json_conf.get("graph_files_path")
    if graph_files_path[-1] != '/':
        graph_files_path += '/'
    output_dir = json_conf.get("imputation_out_path", "output")
    if output_dir[-1] != '/':
        output_dir += '/'
    config = {
        "planb": json_conf.get('planb', True),
        "pops": json_conf.get('populations'),
        "priority": json_conf.get('priority'),
        "epsilon": json_conf.get('epsilon', 1e-3),
        "number_of_results": json_conf.get('number_of_results', 1000),
        "number_of_pop_results": json_conf.get('number_of_pop_results', 100),
        "output_MUUG": json_conf.get("output_MUUG", True),
        "output_haplotypes": json_conf.get("output_haplotypes", False),
        "node_file": project_dir_graph + graph_files_path + json_conf.get("node_csv_file"),
        "top_links_file": project_dir_graph + graph_files_path + json_conf.get("top_links_csv_file"),
        "edges_file": project_dir_graph + graph_files_path +json_conf.get("edges_csv_file"),

        # We changed this line - input file path 
        "imputation_input_file": input_path,
        
        # We changed this line - output files paths
        "imputation_out_umug_freq_file": output_genotype_path,
        "imputation_out_hap_freq_file": output_haplotype_path,

        "imputation_out_umug_pops_file": output_dir + json_conf.get("imputation_out_umug_pops_filename"),
        "imputation_out_hap_pops_file": output_dir + json_conf.get("imputation_out_hap_pops_filename"),
        "imputation_out_miss_file": output_dir + json_conf.get("imputation_out_miss_filename"),
        "imputation_out_problem_file": output_dir + json_conf.get("imputation_out_problem_filename"),
        "factor_missing_data": json_conf.get("factor_missing_data", 0.01),
        "loci_map": json_conf.get("loci_map", {"A": 1, "B":3, "C": 2, "DQB1": 4, "DRB1": 5} ),
        "matrix_planb": json_conf.get("Plan_B_Matrix", [
                    [[1, 2, 3, 4, 5]],
                    [[1, 2, 3], [4, 5]],
                    [[1], [2, 3], [4, 5]],
                    [[1, 2, 3], [4], [5]],
                    [[1], [2, 3], [4], [5]],
                    [[1], [2], [3], [4], [5]]
                ]),
        "pops_count_file": project_dir_graph + json_conf.get("pops_count_file",'' ),
        "use_pops_count_file": json_conf.get("pops_count_file",False),
        "number_of_options_threshold": json_conf.get("number_of_options_threshold", 100000),
        "max_haplotypes_number_in_phase": json_conf.get("max_haplotypes_number_in_phase",100 ),
        "bin_imputation_input_file": project_dir_in_file + json_conf.get("bin_imputation_in_file", "None"),
        "nodes_for_plan_A": json_conf.get("Plan_A_Matrix", []),
        "save_mode": json_conf.get("save_space_mode", False),
        "UNK_priors" : json_conf.get("UNK_priors", "MR")

    }

    all_loci_set = set()
    for _, val in config["loci_map"].items():
        all_loci_set.add(str(val))

    config["full_loci"] = ''.join(sorted(all_loci_set))
    
    # Perform imputation
    config["pops_count_file"] = "../setup/" + config["pops_count_file"]
    imputation = grim.impute_instance(config, graph)

    # Create output directory if it doesn't exist
    pathlib.Path(output_dir).mkdir(parents=False, exist_ok=True)

    # Write out the results from imputation
    imputation.impute_file(config, em_mr=hap_pop_pair)
