# grimm-ard

A web application for HLA imputation. It normalizes a typing with
[py-ard](https://github.com/nmdp-bioinformatics/py-ard) and imputes genotypes and haplotypes with
[py-graph-imputation](https://github.com/nmdp-bioinformatics/py-graph-imputation) (GRIMM).

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Quick start](#quick-start)
3. [Example](#example)
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

> **Note:** The example frequency data is small, so many typings return few or no results. Use your own
> frequency data for meaningful imputation.

---

## Example<a name="example"></a>

These examples use the bundled example data and were checked against it. Populations must match the
`populations` in the configuration file (`AFA` and `CAU` for the example data).

### Imputing one typing with the form

On the home page, enter this typing (as separate alleles, or as a GL string), select the `CAU`
population and the loci `A`, `B`, `C`, `DRB1` and `DQB1`:

```
A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*16:04^DQB1*02:01+DQB1*06:02^DRB1*03:01+DRB1*15:01
```

The results page shows the most likely genotype and haplotype pair, restricted to the selected loci:

| Genotype | Genotype Frequency | Normalized |
| --- | --- | --- |
| `A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*16:04^DQB1*02:01+DQB1*06:02^DRB1*03:01+DRB1*15:01` | 3.21e-06 | 1.000 |

| Haplotype 1 | Haplotype 2 | Genotype Frequency | Normalized |
| --- | --- | --- | --- |
| `A*01:01~B*08:01~C*07:01~DQB1*02:01~DRB1*03:01` | `A*02:01~B*07:02~C*16:04~DQB1*06:02~DRB1*15:01` | 3.21e-06 | 1.000 |

### Imputing a file

Create a CSV file with one subject per line and no header, in the format `ID,GL_STRING,RACE1,RACE2`.
For example, `subjects.csv`:

```
D1,A*01:01+A*02:01^B*08:01+B*07:02^DRB1*03:01+DRB1*15:01,CAU,CAU
```

Upload it on the home page. You get back `grimmard_results.zip`, which contains two CSV files with the
lines `ID,RESULT,PROBABILITY,RANK`:

- `genotype<N>.csv`, the imputed genotypes:

  ```
  D1,A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*16:04^DPA1*01:03+DPA1*01:03^DPB1*04:01+DPB1*294:01^DQA1*01:02+DQA1*05:01^DQB1*02:01+DQB1*06:02^DRB1*03:01+DRB1*15:01^DRB3*01:01+DRB5*01:01,3.210075733065604e-06,0
  ```

- `haplotype<N>.csv`, the imputed haplotype pairs:

  ```
  D1,A*01:01~B*08:01~C*07:01~DPA1*01:03~DPB1*294:01~DQA1*05:01~DQB1*02:01~DRB1*03:01~DRB3*01:01+A*02:01~B*07:02~C*16:04~DPA1*01:03~DPB1*04:01~DQA1*01:02~DQB1*06:02~DRB1*15:01~DRB5*01:01,3.210075733065604e-06,0
  ```

File results always cover all loci in the graph (9 loci for the example data).

> **Note:** With the default configuration, file upload only returns results for typings whose loci
> match an entry in `Plan_A_Matrix` (the full 9 loci, or any combination of up to 3 loci). Other
> subjects, for example a 5-locus A/B/C/DRB1/DQB1 typing, are skipped. The form does not have this
> limitation.

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
