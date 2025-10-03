<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-etd_logo_dark.png">
    <img alt="nf-core/etd" src="docs/images/nf-core-etd_logo_light.png">
  </picture>
</h1>

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

The **nf-core/etd** (Evolving Threat Detector) is a bioinformatics pipeline that analyzes antimicrobial resistome variation of bacterial genomes in comparison to closest relatives and provides contextual analysis of the results. The pipeline is implemented using Nextflow, a  workflow management system designed to execute complex analyses across diverse computing environments with high portability. It leverages Docker containers and conda environments to ensure strong reproducibility of results. As a pipeline, the etd provides:

- **Genome Sketching** 
	1. Genome sketching with Mash ([`Mash`](https://mash.readthedocs.io/en/latest/)) 

- **Antimicrobial Resistance Gene Annotation**

	1. AMR genes with AMRFinderPlus([`AMRFinderPlus`](https://github.com/ncbi/amr))

- **Moblie ELement Detection**
	1. Plasmids with MOB-Suite ([`MOB-Suite`](https://github.com/phac-nml/mob-suite))

	2. Prophages with PhiSpy ([`PhiSpy`](https://github.com/linsalrob/PhiSpy))

	3. Integrons with IntegronFinder ([`IntegronFinder`](https://github.com/gem-pasteur/Integron_Finder))

	4. ICEs with ICEberg ([`ICEberg`](https://ngdc.cncb.ac.cn/databasecommons/database/id/513)) using DIAMOND homology search ([`DIAMOND`](https://github.com/bbuchfink/diamond))

	5. Transposons with TnComp_finder ([`TnFinder`](https://github.com/danillo-alvarenga/tncomp_finder)) and Tn3+TA_finder ([`Tn3+TA_finder`](https://github.com/danillo-alvarenga/tn3-ta_finder)).

- **Closest genomic relative analysis** using Mash distance calculations
- **Resistome comparison** between query genomes and their closest genomic relatives
- **Genomic context mapping** to understand resistance gene mobility


The etd pipeline is organized into two main sub-workflows:

  1. MAKE_DB Sub-workflow

  - Processes reference genomes to build a reference resistance database
  - Annotates genomes using multiple tools for AMR genes and mobile elements
  - Creates Mash sketches for rapid similarity searching
  - Stores all annotations in a searchable SQLite database

  2. QUERY_DB Sub-workflow

  - Takes query genomes and compares them against the reference database
  - Identifies the closest phylogenetic relatives using Mash distances
  - Performs detailed resistome comparisons
  - Reports differences in resistance gene content and mobile element associations

## Dependency Management

- The ETD pipeline uses a hybrid profile called docker_conda to manage dependencies across its various components.

- Most nf-core modules bundled in this pipeline have pre-built Docker containers available via Biocontainers or Docker Hub, which allows for consistent and reproducible execution.

- Custom modules, including parse_gbk`, `Tn3+TA_finder`, `TnComp_finder`, `insert_db` and `query_resistome`, do not yet have Docker containers, and instead rely on Conda environments defined within the pipeline.

- The docker_conda profile allows both environments to coexist during execution. This means:

- Docker is used where stable containers already exist. 

- Conda is used for custom modules where containers have not yet been built.


## Installation

1. Install [`Nextflow`](https://nf-co.re/usage/installation). If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.
2. Install [`Docker`](https://www.docker.com/products/docker-desktop/) if you do not already have it installed.
3. Install [`Conda`](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) if you do not already have it installed.
4. Clone the repository

   ```bash
   git clone -b etd-v0.1 https://git.cs.dal.ca/research1/maguirelab/people/precious-adebola/etd.git

   cd etd
   ```

5. Test with a stub-run. The stub-run will ensure that the pipeline is able to download and use containers as well as execute in the proper logic.
 
   ``` bash
   nextflow run main.nf -profile test_stub,docker_conda -stub
   ```

## Usage

The etd pipeline accepts gbk files whose paths are specified in a samplesheet.csv file as input. 

To run your analysis, first, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,gbk,organism
SAMPLE_1,/path/to/gbk/file,Staphylococcus_aureus

```

Samplesheet.csv must be formatted as above, with the first column corresponding to sample names, the second column corresponding to the locations of the sample gbk file and the third column corresponding to an optional organism / species name if known. In the event that the organism name is unknown, your samplesheet would only have two - sample,gbk - columns.

> [!NOTE]
> Orgamism name when specified should correspond to the format for organism name specification indicated by AMRFinderPlus. (See [`Organism option`](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option))

- For analysis, users can either run individual sub-workflows - `make_db`, `query_db` or run the combined workflow - `all`. Each of these can be specified using the `--mode` parameter. 

1. **etd make_db** : This sub-workflow generates the etd reference database from a bunch of suppleid gennomes - gbk format.
   
```bash
nextflow run main.nf \
   --mode make_db \
   -profile docker_conda \
   --input_make reference_genomes.csv
```

2. **etd query_db**:  This sub-workflow analyzes query genome(s) against the closest reference sequences in the database. Must be run only after a reference database has been built using the make_db subworkflow.

```bash
nextflow run main.nf \
   --mode query_db \
   -profile docker_conda \
   --input_make reference_genomes.csv
```

3. **combined workflow**: This workflow runs both the make_db and query_db subworkflows sequentially.
   
```bash
nextflow run main.nf \
   --mode all \
   -profile docker_conda \
   --input_make reference_genomes.csv 
```

Parameters used:
-  `--mode` : **(Required)** specifies which subworkflow/workflow to run. 
-  `--input_make` : Path to the reference genomes input samplesheet.csv file. **Required** for the `make_db` and `all` modes.
-  `--input_query` : Path to the query genome(s) input samplesheet.csv file. **Required** for the `query_db` and `all` modes.

Optional parameters:

 - `--number` : Number of closest genomes to consider (default: 5)
 - `--output_format` {json, dataframe} : Output format (default:  json)
 - `--min_pident` : minimum percentage identity for diamond homology search (default: 60)
 - `--min_alignment_length` : minimum alignment length for diamond homology search (default: 60)

> [!NOTE]
> To override defaults for optional parameters, please provide pipeline parameters via the CLI. E.g. to override the `--number` and `--min_pident` default parameters:

```bash
nextflow run main.nf \
   --mode all \
   -profile docker_conda \
   --min_pident 80 \
   --number 20 \
   --input_make reference_genomes.csv
```



## Testing

Testing is a crutial part of the development process, ensuring that our tool behaves as expected and that new changes do not introduce bugs. Before running the tests, make sure you have the necessary dependencies installed. You can do this by following the instructions in the installation section.

To test the worklow on a minimal dataset you can use the test configuration (with -profile docker_conda) by executing the following command:

 ```bash
nextflow run main.nf -profile test,docker_conda
```


## Pipeline output

A successful run creates the parent output directory `etd_results` in which other sub annotation and analysis directories are stored. These other directories include:

- `parse/` : fasta, protein and .gff3 files derived from supplied input gbk files.
- `amrfinderplus/` : houses the updated amrfinderplus database and individual sample amr genes and mutations reports.
- `mash/` : combined sketch file for all reference database genomes, individual sample sketch files andmash distance report for query sample(s).
- `phispy/` : per sample .gbk and tsv report files of annotated prohpages, if present.
- `integronfinder/` : per sample report directories of annotated integrons, if present.
- `mobsuite/` : per sample report directories of predicted plamid results including chromosome and plasmid fasta files, if present, contig reports and mob_typer results. 
- `tncomp/` : per sample gbk and .txt reports of annotated composite transposons
- `tn3/` : per sample gnk and .txt reports of annotated tn3 transposons
- `get/` : ICEberg database proteins
- `diamond/` : reformatted iceberg database proteins, and per sample blastp .txt reports.
- `filter/`: report of predicted ICEs following homology search and threshold filtering.
- `insert/`: houses the created etd reference db and any log files.
- `resistome_analysis/` : per sample report of resistome differences between query and closest relatives.

## Credits

nf-core/etd is currently developed by Precious Osadebamwen.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).


## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/etd for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->


> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
