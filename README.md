[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.1)

## Introduction

**nf-asseval** is a bioinformatics pipeline for evaluating genome assemblies. The pipeline is built using [Nextflow](https://www.nextflow.io/) and follows the [nf-core](https://nf-co.re/) community guidelines to ensure high-quality, reproducible, and portable analyses.

1. Run assembly evaluation with [`Quast`](https://github.com/ablab/quast)
2. Run k-mer counting with [`FastK`](https://github.com/thegenemyers/FASTK)
3. Run k-mer based evaluation of raw reads with GeneScopeFK [`GeneScopeFK`](https://github.com/thegenemyers/GENESCOPE.FK)
4. Run k-mer based evaluation of assemblies with [`MerquryFK`](https://github.com/thegenemyers/MERQURY.FK)
5. Generate mappability tracks with [`GenMap`](https://github.com/cpockrandt/genmap)
6. Map reads back to the assembly with [`minimap2`](https://github.com/lh3/minimap2)
7. Calculate coverage statistics based on mapped BAM files
8. Call variants with [`Clair3`](https://github.com/HKU-BAL/Clair3)
9. Generate window-wise plots of relative sequence depth and heterozygosity along chromosomes
10. Generate ERGA Assembly Reports (EARs)

## Usage

### Prerequisites

If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

All the tools used in the pipeline are containerized. Singularity/Apptainer needs to be installed on the machine running the Nextflow runner job. The following containers need to be build from definition files in `assets/containers`:

```
apptainer build pandepth_2.26.sif pandepth.def
```

### Samplesheet

First, prepare a CSV file with your assemblies that looks as follows:

`assemblies.csv`:
```csv
sample,type,status,fasta
genome1,primary,curated,/path/to/fasta
genome2,primary,contig,/path/to/fasta
```

Each row represents an assembly fasta file.

### Reads file

Next, prepare a CSV file with your raw sequence files:

`reads.csv`:
```csv
sample,type,runid,library,fastq_1
Pup3-NMR,hifi,run1,A,test/Pup3-NMR_hifi.fq.gz
Queen-NMR,hifi,run1,A,test/Queen-NMR_hifi.fq.gz
Dad-NMR,hifi,run1,A,test/Dad-NMR_hifi.fq.gz
```

Each row represents a fastq file. Files with the same sample id are merged during the analysis. The combination of sample, runid and library should be unique for each file.

### Parameter file

Lastly, provide a YAML file for the pipeline parameters, referring to the samplesheet with the `input` parameter.

`params.yaml`:
```yaml
# General settings:
#-----------------------
outdir: "results"

# Sample details:
#-----------------------
assemblies: "test/assemblies_curated.csv"
reads: "test/reads_hifi.csv,test/reads_ont_filtered.csv"

# Workflow settings:
#-----------------------
steps: "evaluation,mapping,stats,variant_calling,report"
kmer_size: 31
busco_lineage: "mammalia_odb10"
busco_lineages_path: "/path/to/busco_lineages"
mismatch_penalty: 4
mapq_filter: 0
min_coverage: 10
glnexus_config: "assets/clair3.yml"
publish_bam: false
```

### Running on the IBU cluster

To run the pipeline, start an interactive session on the IBU cluster:

```bash
srun --partition pibu_el8 --project pXXXX-YYYY --cpus-per-task=1 --mem=8000 --time=144:00:00 --pty bash
module load Java
export APPTAINER_CACHEDIR=$SCRATCH
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache
export NXF_TEMP=$SCRATCH
```

Now, you can run the pipeline using:

```bash
nextflow run main.nf \
  -profile unibe_ibu \
  -params-file test/params.yaml \
  --project pXXXX-YYYY
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Local execution

To run the pipeline locally on your computer, use the provided `local.config` file and adjust settings as needed.

```bash
nextflow run main.nf \
  -c local.config \
  -params-file test/params.yaml
```

Make sure to set the environment variable `NXF_SINGULARITY_CACHEDIR` to avoid having to download Singularity containers repeatedly:
```bash
echo $NXF_SINGULARITY_CACHEDIR
```

If not defined, set it before running the pipeline to a directory with sufficient free disk space:
```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache
```

### Clean-up

To completely clean up all previous pipeline runs, the following files and folders need to be deleted:

```bash
rm -r work            # The working cache of the pipeline
rm -r .nextflow       # The Nextflow database containing information for nextflow log
rm -r <outdir>        # The folder with published results as defined in --outdir
rm .nextflow.log*     # The log files of previous runs
```
After this, it is no longer possible to resume a previous run, so be careful what you delete. To just delete a specific pipeline run, use `nextflow log` to get the ids of all runs and remove the run with `nextflow clean <RUN_NAME> -f`.

## Pipeline output

All result files will be published to the directory defined by the `--outdir` parameter.

## Credits

nf-asseval was originally written by Alexander Nater.

We thank the following people for their extensive assistance in the development of this pipeline:
- 

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
