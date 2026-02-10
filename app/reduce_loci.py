import argparse
import json


def write_best_prob_genotype(name_gl, res, fout, numOfResult=10):
    sorted_by_value = sorted(res.items(), key=lambda kv: kv[1], reverse=True)

    # write the output to file
    minBestResult = min(numOfResult, len(sorted_by_value))
    for k in range(minBestResult):
        fout.write(name_gl + ',' + str(sorted_by_value[k][0]) + ',' +
                   str(sorted_by_value[k][1]) + ',' + str(k) + '\n')

def reduce_muug_loci(file_in, file_out, loci, num_res):
    dict_res = {}
    with open(file_in) as six_file:
        for line in six_file:
            id, gl, prob, rank = line.strip().split(',')
            full_gl = gl.split('^')
            gl = []
            for locus in full_gl:
                if locus.split('*')[0] in loci:
                    gl.append(locus)
            gl = '^'.join(sorted(gl))
            if not id in dict_res:
                dict_res[id] = {}
            if gl in dict_res[id]:
                dict_res[id][gl]+= float(prob)
            else:
                dict_res[id][gl] = float(prob)

    six_file.close()

    f_out = open(file_out, 'w')
    for id in dict_res:
        write_best_prob_genotype(id, dict_res[id], f_out, num_res)
    f_out.close()

def reduce_haps_loci(file_in, file_out, loci, num_res):
    dict_res = {}
    with open(file_in) as six_file:
        for line in six_file:
            id, hap1, hap2, prob, rank = line.strip().split(',')
            haps = []
            for hap in [hap1, hap2]:
                hap = hap.split(';')[0]
                hap = hap.split("~")
                hap_tmp =[]
                for locus in hap:
                    if locus.split('*')[0] in loci:
                        hap_tmp.append(locus)
                haps.append(('~').join(sorted(hap_tmp)))


            haps = ('+').join(sorted(haps))
            if not id in dict_res:
                dict_res[id] = {}
            if haps in dict_res[id]:
                dict_res[id][haps]+= float(prob)
            else:
                dict_res[id][haps] = float(prob)

    six_file.close()

    f_out = open(file_out, 'w')
    for id in dict_res:
        write_best_prob_genotype(id, dict_res[id], f_out, num_res)
    f_out.close()


# parser = argparse.ArgumentParser()
# parser.add_argument("-c", "--config",
#                     required=False,
#                     default="../conf/minimal-configuration.json",
#                     help="Configuration JSON file",
#                     type=str)

def reduce_loci(loci,  imputation_out_muug_freqs, imputation_out_hap_freqs):
    num_res = 20

    reduce_muug_loci(imputation_out_muug_freqs, (".").join(imputation_out_muug_freqs.split('.')[:-1]) + "_reduced.txt", loci, num_res)
    path_mug = (".").join(imputation_out_muug_freqs.split('.')[:-1]) + "_reduced.txt"
    output_file =(".").join(imputation_out_hap_freqs.split('.')[:-1])+ "_reduced.txt"
    reduce_haps_loci(imputation_out_hap_freqs, (".").join(imputation_out_hap_freqs.split('.')[:-1])+ "_reduced.txt" , loci, num_res)
    return output_file, path_mug
