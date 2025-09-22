<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-etd_logo_dark.png">
    <img alt="nf-core/etd" src="docs/images/nf-core-etd_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/etd/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/etd/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23etd-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/etd)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

The **nf-core/etd** (Evolving Threat Detector) is a bioinformatics pipeline that analyzes antimicrobial resistome variation of bacterial genomes in comparison to closest relatives and provides contextual analysis of the results. The pipeline is implemented using Nextflow, a  workflow management system designed to execute complex analyses across diverse computing environments with high portability. It leverages Docker containers and conda environments to ensure strong reproducibility of results.

- Genome Sketching 
1. Genome sketching with Mash ([`Mash`](https://mash.readthedocs.io/en/latest/)) 

- Features Prediction 
2. AMR genes with AMRFinderPlus([`AMRFinderPlus`](https://github.com/ncbi/amr)) 
3. Plasmids with MOB-Suite ([`MOB-Suite`](https://github.com/phac-nml/mob-suite)) 
4. Prophages with PhiSpy ([`PhiSpy`](https://github.com/linsalrob/PhiSpy)) 
5. Integrons with IntegronFinder ([`IntegronFinder`](https://github.com/gem-pasteur/Integron_Finder)) 
6. ICEs with ICEberg ([`ICEberg`](https://ngdc.cncb.ac.cn/databasecommons/database/id/513)) using DIAMOND homology search ([`DIAMOND`](https://github.com/bbuchfink/diamond)) 
7. Transposons with TnComp_finder ([`TnFinder`](https://github.com/danillo-alvarenga/tncomp_finder)) and Tn3+TA_finder ([`Tn3+TA_finder`](https://github.com/danillo-alvarenga/tn3-ta_finder)).

The etd pipeline consists of two main sub-workflows:

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


## Installation

1. Install [`Nextflow`](https://nf-co.re/usage/installation). If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.
2. Install [`Docker`](https://www.docker.com/)
3. Install [`Conda`](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)
4. Clone the repository

   ```bash
   git clone

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

Orgamism name when specified should correspond to the format for organism name specification indicated by AMRFinderPlus. (See [`Organism option`](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option)

- For analysis, the etd pipeline is organized into two main sub-workflows - make_db, query_db and one combined workflow - all. Each of these can be specified using the --mode parameter. 

1. **etd make_db** : This sub-workflow generates the etd reference database from a bunch of suppleid gennomes - gbk format.
   
```bash
nextflow run main.nf \
   --mode make_db \
   - profile ddocker_conda \
   --input_make reference_genomes.csv
```

2. **etd query_db**:  This sub-workflow analyzes query genome(s) against the closest reference sequences in the database. Must be run only after a refernce database has been built using the make_db subworkflow.

```bash
nextflow run main.nf \
   --mode query_db \
   - profile ddocker_conda \
   --input_make reference_genomes.csv
```

3. **combined workflow**: This workflow runs both the make_db and query_db subworkflows sequentially.
   
```bash
nextflow run main.nf \
   --mode all \
   - profile ddocker_conda \
   --input_make reference_genomes.csv 
```

Parameters used:
-  `--mode` : **(Required)** specifies which subworkflow/workflow to run. 
-  `--input_make` : Path to the reference genomes input samplesheet.csv file. **Required** for the `make_db` and `all` modes.
-  `--input_query` : Path to the query genome(s) input samplesheet.csv file. **Required** for the `query_db` and `all` modes.

Optional parameters:

 - `--number` : Number of closest genomes to consider (default: 5)
 - `--output_forrmat` {json, dataframe} : Output format (default:  json)
 - `--min_pident` : minimum percentage identity for diamond homology search (default: 60)
 - `--min_alignment_length` : minimum alignment length for diamond homology search (default: 60)

> [!NOTE]
> To override defaults for optional parameters, please provide pipeline parameters via the CLI

## Testing

Testing is a crutial part of the development process, ensuring that our tool behaves as expected and that new changes do not introduce bugs. Before running the tests, make sure you have the necessary dependencies installed. You can do this by following the installation instructions in the usage section.

To test the worklow on a minimal dataset you can use the test configuration (with -profile docker_conda) by executing the following command:

 ```bash
nextflow run main.nf -profile test,docker_conda
```


## Pipeline output

A successful run creates the parent outpyt directory `etd_results` in which other sub annotation and analysis directories are stored. These other directories include:

- `parse/` : fasta, protein and .gff3 files derived from supplied input gbk files.
- `amrfinderplus/` ; houses the updated amrfinderplus database and individual sample amr genes and mutations reports.
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

nf-core/etd was originally written by Precious Osadebamwen.

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
