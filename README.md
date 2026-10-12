# grimm-ard

A web application for HLA imputation. It normalizes a typing with
[py-ard](https://github.com/nmdp-bioinformatics/py-ard) and imputes genotypes and haplotypes with
[py-graph-imputation](https://github.com/nmdp-bioinformatics/py-graph-imputation) (GRIMM).

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Quick start](#quick-start)
3. [Using the application](#using-the-application)
4. [Using your own frequency data](#using-your-own-frequency-data)
    - [Configuration file](#configuration-file)
    - [Building the graph](#building-the-graph)
    - [Using an existing graph](#using-an-existing-graph)
5. [Project layout](#project-layout)

---

## Prerequisites<a name="prerequisites"></a>

- Python 3.10 or newer.
- A C compiler (e.g. `build-essential` on Debian/Ubuntu, Xcode Command Line Tools on macOS).
  `py-graph-imputation` may be built from source with Cython during installation.
- Internet access on the first run: py-ard downloads the IMGT/HLA database and caches it locally.

Clone the repository and install the dependencies, preferably in a virtual environment:

```
git clone https://github.com/louzounlab/grimmard.git
cd grimmard
python3 -m venv .venv
source .venv/bin/activate
pip install -r app/requirements.txt
```

`py-graph-imputation` 1.0.0 or newer is required.

---

## Quick start<a name="quick-start"></a>

The repository ships with 9-locus haplotype frequencies for two populations (`AFA`, `CAU`) in
`app/data/freq_9loci/`. To run the application on this example data:

```
cd app
python3 app.py
```

Then open http://127.0.0.1:5000/.

All paths are relative to `app/`, so **always start the application from the `app` directory**.

On the first start, `app.py` checks for `data/graph.pkl`, `data/pop_ratio.txt` and `data/freqs_dicts/`.
If any are missing it runs `produce_example_graph_file.py` to generate them, which takes a minute or so.
You can also run that step yourself beforehand:

```
cd app
python3 produce_example_graph_file.py
```

It prints the configuration it uses and finishes with:

```
1. Produced: output_new/hpf.csv
2. Produced nodes and edges: output_new/csv
3. Produced Whole Graph: data/graph.pkl
4. Produced Pickled Freqs: data/freqs_dicts/all_freqs.pickle
```

> **Note:** `app.py` runs Flask in debug mode on `0.0.0.0`. This is fine on your own machine, but do not
> expose it on a network: the debugger allows anyone who reaches the port to execute code.

> **Note:** The example frequency data is small, so many typings return few or no results. Use your own
> frequency data for meaningful imputation.

---

## Using the application<a name="using-the-application"></a>

There are two ways to submit typings:

- **Form:** Enter alleles per locus (or a GL string) and the population, and optionally choose which loci
  to show in the results. The results page lists the most likely genotypes and haplotype pairs.
- **File upload:** Upload a CSV file with one subject per line and no header:

  ```
  ID,GL_STRING,RACE1,RACE2
  ```

  For example:

  ```
  D1,A*01:01+A*02:01^B*08:01+B*07:02^DRB1*03:01+DRB1*15:01,CAU,CAU
  ```

  You get back a zip file with the genotype and haplotype results.

Populations must match the `populations` in the configuration file (`AFA` and `CAU` for the example data).

---

## Using your own frequency data<a name="using-your-own-frequency-data"></a>

### Configuration file<a name="configuration-file"></a>

The application reads its configuration from `app/conf/conf.json`. The parameters are described in
[app/conf/README.md](app/conf/README.md). The ones you most likely need to change are:

| Parameter | Description |
| --- | --- |
| `populations` | Population names. Each must have a frequency file `<POP>.freqs.gz` in `freq_data_dir`. |
| `freq_data_dir` | Directory with the haplotype frequency files (default `data/freq_9loci`). |
| `loci_map` | Maps each locus name to its index in the haplotype. |
| `Plan_A_Matrix` | Combinations of locus indexes that imputation can start from. |
| `Plan_B_Matrix` | Ways to split the full haplotype when no full-haplotype match is found. |
| `freq_trim_threshold` | Haplotypes below this frequency are left out of the graph. |
| `number_of_results` | Maximum number of genotypes/haplotype pairs returned per subject. |

Frequency files are gzipped CSVs with a header line and one haplotype per line, with alleles joined
by `~`:

```
Haplo,Count,Freq
A*33:03~C*02:10~B*15:03~DRB4*01:01~DRB1*04:04~DQA1*03:01~DQB1*03:02~DPA1*01:03~DPB1*04:01,1.629346,0.001172078432
```

### Building the graph<a name="building-the-graph"></a>

After placing your frequency files in `freq_data_dir` and updating `conf/conf.json`, delete any
previously generated files and rebuild:

```
cd app
rm -rf data/graph.pkl data/pop_ratio.txt data/freqs_dicts output_new
python3 produce_example_graph_file.py
```

`produce_example_graph_file.py` runs these steps:

1. Converts the frequency files into a haplotype-population-frequency file (`output_new/hpf.csv`) and
   writes population ratios to `data/pop_ratio.txt`.
2. Builds the node and edge CSV files in `output_new/csv/`.
3. Builds the imputation graph and saves it to `data/graph.pkl`.
4. Saves the haplotype frequencies per population to `data/freqs_dicts/all_freqs.pickle`.

### Using an existing graph<a name="using-an-existing-graph"></a>

If you already have the node and edge CSV files, you can skip steps 1 and 2. Place:

- the CSV files (`nodes.csv`, `edges.csv`, `top_links.csv`, `info_node.csv`) in `app/output_new/csv/`,
- the population ratio file in `app/data/pop_ratio.txt`,
- the configuration in `app/conf/conf.json`,

then comment out steps 1 and 2 in `produce_example_graph_file.py` and run it.

---

## Project layout<a name="project-layout"></a>

```
app/
├── app.py                         # Flask application
├── imputation.py                  # py-ard + GRIMM imputation used by the web routes
├── runfile.py                     # Builds the imputation config and runs GRIMM
├── reduce_loci.py                 # Restricts results to the loci selected by the user
├── produce_example_graph_file.py  # Generates the graph and frequency pickles
├── conf/conf.json                 # Configuration
├── data/freq_9loci/               # Example haplotype frequencies (AFA, CAU)
├── templates/, static/            # Web pages and assets
├── data/graph.pkl                 # Generated: imputation graph
├── data/pop_ratio.txt             # Generated: population ratios
├── data/freqs_dicts/              # Generated: pickled frequencies
└── output_new/                    # Generated: intermediate files and imputation output
```
