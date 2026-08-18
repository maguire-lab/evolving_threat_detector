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
        2. Batch sketch pasting for large reference databases

> [!NOTE]
> For reference databases exceeding 50,000 genomes, the pipeline automatically splits 
> sketches into batches of 50,000 before merging them into a final combined sketch. 
> This avoids hitting the Linux kernel's memory-mapped file limit 
> (`vm.max_map_count`) that causes single-pass `mash paste` to fail on very large 
> collections. The intermediate batch files are automatically removed after the 
> final sketch is produced, so no additional storage is consumed.

- **Antimicrobial Resistance Gene Annotation**

	1. AMR genes with AMRFinderPlus([`AMRFinderPlus`](https://github.com/ncbi/amr))

- **Moblie ELement Detection**
	1. Plasmids with MOB-Suite ([`MOB-Suite`](https://github.com/phac-nml/mob-suite))

	2. Prophages with PhiSpy ([`PhiSpy`](https://github.com/linsalrob/PhiSpy))

	3. Integrons with IntegronFinder ([`IntegronFinder`](https://github.com/gem-pasteur/Integron_Finder))

	4. ICEs with ICEberg ([`ICEberg`](https://ngdc.cncb.ac.cn/databasecommons/database/id/513)) using DIAMOND homology search ([`DIAMOND`](https://github.com/bbuchfink/diamond))

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

- The ETD pipeline uses docker to manage dependencies across its various components.

- nf-core modules bundled in this pipeline have pre-built Docker containers available via Biocontainers or Docker Hub, which allows for consistent and reproducible execution. This pipeline also supports Singularity/Apptainer for HPC environments where Docker is unavailable.

- Custom modules, including `parse_gbk`, `insert_db` and `query_resistome`, also use Docker containers.

## Databases and External Resources
Across its execution, the ETD interacts with a number of databases. These can be classified into three categories: databases downloaded at runtime, databases bundled within tool containers, and databases generated during the pipeline run. 

- Databases Downloaded at Runtime

- AMRFinderPlus Reference Database: at the start of each run, the ETD pipeline automatically downloads the latest version of the curated antimicrobial resistance gene database hosted by NCBI. This download requires internet access and is achieved via the `AMRFINDERPLUS_UPDATE` process.
- ICEberg Protein Database: The ICEberg datatbase houses experimentally verified integrative and conjugative element (ICE) protein sequences. Via the `GET_ICEBERG` process, the ETD pipeline downloads this database, reformats the FASTA headers for compatibility with DIAMOND, and builds a DIAMOND index from it. This download happens once per pipeline run and requires internet access.

- Other tools such as `MOB-suite`, `IntegronFinder` and `PhiSpy` have their databses and respective models bundled with the tool, and requires no user configuration.

- Databases Generated During Pipeline Execution
- ETD SQLite Database (etd.db): This relational database forms the core output of the `make_db` workflow and houses all annotated genomes with their AMR genes cross-referenced against corresponding mobile genetic element contexts. Location: `etd_results/insert`.
- Mash Sketch Reference (`all_genomes.msh`): A combined MinHash sketch of all reference genomes generated by `MASH_SKETCH`and `MASH_PASTE` processes.
- DIAMOND Index(`.dmnd`): A binary protein search index built from ICEberg FASTA by the `DAIMOND_MAKEDB` process. This file is not published to the output directory.

## Installation

1. Install [`Nextflow`](https://nf-co.re/usage/installation). If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.
2. Install [`Docker`](https://www.docker.com/products/docker-desktop/) if you do not already have it installed or Singularity/Apptainer if running on an HPC cluster where Docker is not available.
3. Install [`Conda`](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) if you do not already have it installed.
4. Clone the repository

   ```bash
   git clone -b etdv0.1 https://github.com/maguire-lab/evolving_threat_detector.git

   cd evolving_threat_detector
   ```

5. Test with a stub-run. The stub-run will ensure that the pipeline is able to download and use containers as well as execute in the proper logic.
 
   ``` bash
   nextflow run main.nf -profile test_stub,docker -stub
   ```

## Running on Shared HPC Clusters
On most shared HPC clusters, compute nodes do not have internet access. Pre-download all required databases and container images on a login node before submitting jobs.

1. Pre-download the AMRFinderPlus database:
   ``` bash
   amrfinder_update -d amrfinderdb
   ```
   Note the path to the resulting amrfinderdb directory for use with `--amrfinder_db`.

2. Pre-download and reformat the ICEberg database:
   ``` bash
   curl https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/ICE_aa_experimental.fas \
   --output ICE_aa_experimental.fas

   sed -E 's/>(ICEberg\|[0-9]+)\s+(gi.*\|)\s+(.*)\s+\[(.*)\]/>\1_\3_\2_[\4]/' \
   < ICE_aa_experimental.fas | tr ' ' '_' > ICE_aa_experimental_reformatted.fas
   ```
   Note the path to ICE_aa_experimental_reformatted.fas for use with `--iceberg_db`.

3. Pre-download container images:
   Set a persistent cache directory for Singularity images:
   ``` bash
   export NXF_SINGULARITY_CACHEDIR=/path/to/persistent/singularity_cache
   mkdir -p $NXF_SINGULARITY_CACHEDIR
   ```
   Then, run the test profile on a login node to pull all required container images:
   ``` bash
   nextflow run main.nf -profile test,singularity
   ```
   This executes the full pipeline on a minimal test dataset, forcing Nextflow to download and cache all Singularity images into NXF_SINGULARITY_CACHEDIR. On subsequent runs (including on compute nodes without internet), the cached images will be used automatically.

4. Set offline mode and submit:In your SLURM job script, set `NXF_OFFLINE` to `true` and `NXF_SINGULARITY_CACHEDIR` to `/path/to/persistent/singularity_cache`

## Usage

The etd pipeline accepts **annotated** GenBank files (`.gbk/.gbff/.gb`) whose paths are specified via a samplesheet.csv file as input. Unannotated GenBank files (i.e, those containing only raw nucleotide sequences without CDS features) will cause downstream tools to fail.

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
   -profile docker \
   --input_make reference_genomes.csv
```
**Multi-batch database building:**
To append to an existing database from a previous run, supply the path to
the previous database using `--db_path`:

```bash
nextflow run main.nf \
   --mode make_db \
   -profile docker \
   --input_make reference_genomes.csv \
   --db_path /path/to/previous_run/insert/etd.db
```

2. **etd query_db**:  This sub-workflow analyzes query genome(s) against the closest reference sequences in the database. Must be run only after a reference database has been built using the make_db subworkflow.

```bash
nextflow run main.nf \
   --mode query_db \
   -profile docker \
   --input_query query_genomes.csv \
   --outdir /path/to/make_db_results
```

> [!NOTE]
> On HPC clusters where compute nodes lack internet access, pass `--amrfinder_db` and `--iceberg_db` to provide pre-downloaded databases. This applies to both `make_db` and `query_db` modes, as the pipeline
> cannot download these databases without internet. On systems with internet, `make_db` downloads and stores these automatically, and `query_db` loads them from `--outdir`.

3. **combined workflow**: This workflow runs both the make_db and query_db subworkflows sequentially.
   
```bash
nextflow run main.nf \
   --mode all \
   -profile docker \
   --input_make reference_genomes.csv \
   --input_query query_genomes.csv \
   --keep_all
```

Parameters used:
-  `--mode` : **(Required)** specifies which subworkflow/workflow to run. 
-  `--input_make` : Path to the reference genomes input samplesheet.csv file. **Required** for the `make_db` and `all` modes.
-  `--input_query` : Path to the query genome(s) input samplesheet.csv file. **Required** for the `query_db` and `all` modes.
-  `--keep_all` : When specified, publishes all intermediae output files from individual process (default: `false`). By default, only essential files are published.

Optional parameters:

 - `--number` : Number of closest genomes to consider (default: 5)
 - `--output_format` {json, dataframe} : Output format (default:  json)
 - `--max_evalue` : maximum e_value cutoff for diamond homology search (default: 1e-5)
 - `--min_coverage` : minimum coverage of both queryy and subject length for diamond homology search (default: 0.5)
 - `--max_distance` : maximum co-location distance in base pairs between AMR and MGE (default: 5000) 
 - `--amrfinder_db` : path to a pre-downloaded AMRFinderPlus database directory. When provided, the pipeline skips the AMRFINDERPLUS_UPDATE step (useful for HPC compute nodes without internet).
 - `--iceberg_db` : path to a pre-downloaded ICEberg protein FASTA file. When provided, the pipeline skips the iceberg database download step (useful for HPC compute nodes without internet).

> [!NOTE]
> To override defaults for optional parameters, please provide pipeline parameters via the CLI. E.g. to override the `--number` and `--min_coverage` default parameters:

```bash
nextflow run main.nf \
   --mode all \
   -profile docker \
   --min_coverage 0.6 \
   --number 20 \
   --input_make reference_genomes.csv \
   --input_query query_genomes.csv
```

## Testing

Testing is a crutial part of the development process, ensuring that our tool behaves as expected and that new changes do not introduce bugs. Before running the tests, make sure you have the necessary dependencies installed. You can do this by following the instructions in the installation section.

To test the worklow on a minimal dataset you can use the test configuration (with -profile docker) by executing the following command:

 ```bash
nextflow run main.nf -profile test,docker
```

> On HPC systems using Singularity/Apptainer, replace `docker` with the appropriate Singularity-based profile. If on a compute node without internet, ensure databases and containers are pre-downloaded and pass `--amrfinder_db` and `--iceberg_db`.

``` bash
  nextflow run main.nf -profile test,singularity --amrfinder_db /path/to/amrfinderdb --iceberg_db /path/to/iceberg.fas
  ```


## Pipeline output

A successful run creates the parent output directory `etd_results` in which other sub annotation and analysis directories are stored. By default, only essential output files are published to minimize filesystem usage. When `--keep_all` is specified, all intermediate files are published.These other directories include:

- `parse/` : fasta, protein and .gff3 files derived from supplied input gbk files.
- `amrfinderplus/` : default: per-sample .tsv report only. With `--keep_all`: also includes mutation reports.
- `mash/` : default: combined sketch file for all reference database genomes, and mash distance report for query sample(s). With `--keep_all`: includes individual sketch files.
- `phispy/` : default: per-sample .tsv coordinates file only. With `--keep_all`: also includes annotated .gbk files.
- `integronfinder/` : default: summary and .integrons files only. With `--keep_all`: full report directories.
- `mobsuite/` : default: contig report only. With --keep_all: also includes mob_typer results, chromosome/plasmid FASTA files.
- `get/` : ICEberg database proteins
- `diamond/` : reformatted iceberg database proteins, and per sample blastp .txt reports.
- `filter/`: report of predicted ICEs following homology search and threshold filtering.
- `insert/`: houses the etd reference SQLite db and any log files. For multi-batch runs, use this database as input to subsequent batches via `--db_path`.
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
