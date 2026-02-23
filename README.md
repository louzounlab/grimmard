# grimm-ard

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Example](#example)
3. [Updating the Website](#updating-the-website)
    - [Configuration File](#configuration-file)
    - [Generating HPF File](#hpf-file)
    - [Preparing the Incidence Graph](#preparing-the-incidence-graph)
    - [Configuring the application](#configure-app)
---

## Prerequisites<a name="prerequisites"></a>

- Access to the **grimm-ard** code repository on
  GitHub: [GitHub Repository](https://github.com/nmdp-bioinformatics/grimm-ard)
- `pip` install required libraries.
    ```
    pip install -r app/requirements.txt
    ```
- Necessary files and configurations as mentioned below.

---
### Example<a name="example"></a>

The following example shows a working application with minimal setup. See the next section for working with your own data set. 

#### Example of a minimal **grimm-ard** application working.

There are sample freq. files and configuration files available in `example` directory. To use those example files for **grimm-ard**, follow:

##### Produce the needed graph and freq. files

The script will create the needed intermediary files in `output` directory and generate a graph pickle file needed for
the web application.

Run the python script:

```
cd setup
python3 produce_example_graph_file.py
```

output as it produces the needed files:

```
****************************************************************************************************
Conversion to HPF file based on following configuration:
	Population: ['AFA', 'CAU']
	Frequency File Directory: example/data/freqs
	Output File: output/hpf.csv
****************************************************************************************************
Reading Frequency File:	 example/data/freqs/AFA.freqs.gz
Reading Frequency File:	 example/data/freqs/CAU.freqs.gz
Writing hpf File:	 output/hpf.csv
1. Produced: output/hpf.csv
****************************************************************************************************
Performing graph generation based on following configuration:
	Population: ['AFA', 'CAU']
	Freq File: output/hpf.csv
	Freq Trim Threshold: 1e-05
****************************************************************************************************
2. Produced nodes and edges: output/csv
****************************************************************************************************
Performing imputation based on:
	Population: ['AFA', 'CAU']
	Priority: {'alpha': 0.4999999, 'eta': 0, 'beta': 1e-07, 'gamma': 1e-07, 'delta': 0.4999999}
	UNK priority: SR
	Epsilon: 0.001
	Plan B: True
	Number of Results: 10
	Number of Population Results: 100
	Nodes File: output/csv/nodes.csv
	Top Links File: output/csv/edges.csv
	Input File: data/subjects/donor.csv
	Output UMUG Format: True
	Output UMUG Freq Filename: output/don.umug
	Output UMUG Pops Filename: output/don.umug.pops
	Output Haplotype Format: True
	Output HAP Freq Filename: output/don.pmug
	Output HAP Pops Filename: output/don.pmug.pops
	Output Miss Filename: output/don.miss
	Output Problem Filename: output/don.problem
	Factor Missing Data: 0.0001
	Loci Map: {'A': 1, 'B': 2, 'C': 3, 'DQB1': 4, 'DRB1': 5}
	Plan B Matrix: [[[1, 2, 3, 4, 5]], [[1, 2, 3], [4, 5]], [[1], [2, 3], [4, 5]], [[1, 2, 3], [4], [5]], [[1], [2, 3], [4], [5]], [[1], [2], [3], [4], [5]]]
	Pops Count File: output/pop_counts_file.txt
	Use Pops Count File: output/pop_counts_file.txt
	Number of Options Threshold: 100000
	Max Number of haplotypes in phase: 100
	Save space mode: False
****************************************************************************************************
3. Produced Whole Graph: app/graph.pickle
Reading file: CAU.freqs.gz
Reading file: AFA.freqs.gz
4. Produced Pickled Freqs: app/freqs_dicts/all_freqs.pickle
5. Produced app/pop_ratio.txt
```

##### Start the application

Change to `app` directory and start the web app.

```
cd ../app
python3 app.py
```

##### Use Application

Visit http://127.0.0.1:5000/ to test out the application locally.

---

## Updating the Website<a name="updating-the-website"></a>

To update the Grim_ard website, you will need to follow these steps:

### Configuration File<a name="configuration-file"></a>

Ensure you have a JSON configuration file prepared. This file is crucial for configuring various aspects of the website.

See example file: [Configuration File](example/conf/minimal-configuration.json)

### HPF File<a name="hpf-file"></a>

Generate a HPF file from the freq. file based on configuration.

```python
from graph_generation.generate_hpf import produce_hpf

config_file = 'app/conf/conf.json'
produce_hpf(conf_file=config_file)
```

This will produce a `setup/output/hpf.csv`. The output directory should look as:

```
output
├── csv
├── hpf.csv
└── pop_counts_file.txt
```

The `output/pop_counts_file.txt` needs to be copied to `app/pop_ratio.txt`:

### Create Nodes/Edges files

From HPF file, create nodes/edges CSV files.

```python
from grim import grim

config_file = 'example/conf/minimal-configuration.json'

grim.graph_freqs(config_file)
```

You should see the nodes/edges CSV files.
```
|-- csv
|   |-- edges.csv
|   |-- info_node.csv
|   |-- nodes.csv
|   `-- top_links.csv

```

### Preparing the Incidence Graph<a name="preparing-the-incidence-graph"></a>

Using nodes/edges CSV files, produce a graph

```python
from generate_config_dict import generate_dict_config
from grim import grim

config_file = 'example/conf/minimal-configuration.json'
dict_config = generate_dict_config(config_file)

graph = grim.graph_instance(dict_config)
```

Ensure you have the `graph.pickle` file containing your graph instance. This pickle file will be used during the website
update process.

Save the graph as a pickle file in the `data/` directory.

```python
import pickle

with open('data/graph.pickle', 'wb') as pickle_file:
    pickle.dump(graph, pickle_file)
```

The frequency files also need to be imported and pickled with haplotype->freq dictionary to create a pickle file
in `app/freqs_dicts/all_freqs.pickle`.

Adapt the [example Python script](#example-of-a-minimal-grimm-ard-application-working) to your dataset.

### Configure application<a name="configure-app"></a>

Create `WEBSITE_CONFIG_FILE` and `GRIM_GRAPH` environment variables to point to config and pickled graph object file
relative to `app/` directory.

Change `WEBSITE_CONFIG_FILE` to the path of your configuration file.

```
export WEBSITE_CONFIG_FILE=my-website-configuration.json
```

Change `GRIM_GRAPH` to the path of the frequency graph you created and saved e.g. 'my-graph.pickle'.

```
export GRIM_GRAPH=my-graph.pickle
```

Here’s the text you can add to your **README** file:

---

### **Note: Using an Existing Graph**
If you **already have a graph**, you only need to run the Python script **`produce_example_graph_file.py`** starting from **Step 3** (**skip Step 1 and Step 2**).

#### **Required File Locations**
Make sure to place the necessary files in the correct directories:
- **Graph-related files** should be inside the **`csv`** directory.
- The **path for `csv`** should be:
  ```
  setup/output_new/csv
  ```
- The **population counts file (`pop_counts_file.txt`)** should be located at:
  ```
  setup/output_new/pop_counts_file.txt
  ```
- The **configuration file (`conf.json`)** should be placed at:
  ```
  app/conf/conf.json
  ```

Ensure that all required files are correctly placed before running the script.

---

