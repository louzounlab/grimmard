import os
import pickle
import shutil
import gzip
from graph_generation.generate_hpf import produce_hpf
from generate_config_dict import generate_dict_config
from grim import grim

# Step 1: Create HPF File
config_file = '../app/conf/conf.json'
produce_hpf(conf_file=config_file)
print("1. Produced: output_new/hpf.csv")

# Step 2: Create nodes/edges
grim.graph_freqs(config_file)
print("2. Produced nodes and edges: output_new/csv")

# Step 3: Create networkx Graph
dict_config = generate_dict_config(config_file)
# Create a graph instance based on your configuration
graph = grim.graph_instance(dict_config)

# Pickle it out
# Assuming 'graph' is your graph instance, save it as a pickle
os.makedirs("../app/data", exist_ok=True)
with open('../app/data/graph.pkl', 'wb') as pickle_file:
    pickle.dump(graph, pickle_file)
print("3. Produced Whole Graph: app/data/graph.pkl")

# Generate the freq dictionary file
os.makedirs("../app/data/freqs_dicts", exist_ok=True)

freqs_files = os.listdir("data/freq_9loci/")
all_freqs = {}

for file in freqs_files:
    print(f"Reading file: {file}")
    freqs = {}
    with gzip.open(f"data/freq_9loci/{file}", 'rt') as f:
        for i, line in enumerate(f):
            if i == 0:
              continue
            hap, _, prob = line.split(",")
            hap = "~".join(sorted(hap.split("~")))
            prob == 'Freq\n'
            freqs[hap] = float(prob)

    pop = file.split('.')[0]
    all_freqs[pop] = freqs

with open(f"../app/data/freqs_dicts/all_freqs.pickle", "wb") as f:
    pickle.dump(all_freqs, f)
print("4. Produced Pickled Freqs: app/data/freqs_dicts/all_freqs.pickle")

# Copy pop ratio file
shutil.copy('output_new/pop_counts_file.txt', '../app/data/pop_ratio.txt')
print("5. Produced app/data/pop_ratio.txt")