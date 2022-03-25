# ViReflow
ViReflow is a tool for constructing elastically-scaling parallelized automated AWS pipelines for viral consensus sequence generation. Given sequence data from a viral sample as well as information about the reference genome and primers, ViReflow generates a [Reflow](https://github.com/grailbio/reflow) file that contains all steps of the workflow, including AWS instance specifications. Because ViReflow is intended to be used with Reflow, the workflows that are developed by ViReflow automatically distribute independent tasks to be run in parallel as well as elastically scale AWS instances based on each individual step of the workflow. ViReflow makes use of compact minimal Docker images for each step of the viral analysis workflow, details about which can be found in the [Niema-Docker](https://github.com/Niema-Docker) GitHub organization.

## Workflow Summary
| Read Trimmers | Read Mappers | Variant Callers | Optional Analyses    |
| ------------- | ------------ | --------------- | -----------------    |
| fastp         | Bowtie2      | FreeBayes       | coronaSPAdes         |
| iVar Trim     | BWA-MEM      | iVar Variants   | MEGAHIT              |
| PRINSEQ       | Minimap2     | LoFreq          | metaviralSPAdes      |
| pTrimmer      |              |                 | minia                |
|               |              |                 | Pangolin (COVID-19)  |
|               |              |                 | rnaviralSPAdes       |
|               |              |                 | VirStrain            |
|               |              |                 | *π* Diversity Metric |

* **Trim the reads** using **[iVar](https://github.com/andersen-lab/ivar) (default)**, [fastp](https://github.com/OpenGene/fastp), [PRINSEQ](http://prinseq.sourceforge.net/), or [pTrimmer](https://github.com/DMU-lilab/pTrimmer)
* **Map the reads** using **[Minimap2](https://github.com/lh3/minimap2) (default)**, [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [BWA](http://bio-bwa.sourceforge.net/)
* **Generate a pile-up** from the trimmed mapped reads using **[samtools](http://www.htslib.org/)**
* **Call variants** using **[LoFreq](https://csb5.github.io/lofreq/) (default)**, [freebayes](https://github.com/freebayes/freebayes), [iVar](https://github.com/andersen-lab/ivar)
* **Calculate depth** from the trimmed mapped reads using **[samtools](http://www.htslib.org/)**
* **Call a consensus sequence** of high-depth regions from the variants using **[bcftools](http://samtools.github.io/bcftools/bcftools.html)**
* ***OPTIONAL:***
  * **Assign a COVID-19 lineage** using **[Pangolin](https://pangolin.cog-uk.io/)** (SARS-CoV-2 only) and/or **[VirStrain](https://github.com/liaoherui/VirStrain)**
  * **Assemble a *de novo* genome** using **[coronaSPAdes](https://cab.spbu.ru/software/coronaspades)**, **[metaviralSPAdes](https://doi.org/10.1093/bioinformatics/btaa490)**, **[rnaviralSPAdes](https://github.com/ablab/spades#supported-data-types)**, **[MEGAHIT](https://github.com/voutcn/megahit)**, and/or **[minia](https://github.com/GATB/minia)**

## Installation
ViReflow is written in Python 3. You can simply download [ViReflow.py](ViReflow.py) to your machine and make it executable:

```bash
wget "https://raw.githubusercontent.com/niemasd/ViReflow/master/ViReflow.py"
chmod a+x ViReflow.py
sudo mv ViReflow.py /usr/local/bin/ViReflow.py # optional step to install globally
```

While ViReflow itself only depends on Python 3, the pipelines it produces are [Reflow](https://github.com/grailbio/reflow) files that run on AWS. Thus, in order to run the pipelines ViReflow produces, one must first [install Reflow](../../wiki/Installing-Reflow).

## Usage
ViReflow can be used as follows:

```
usage: ViReflow.py [-o OUTPUT_RF] -d DESTINATION -rf REFERENCE_FASTA -rg REFERENCE_GFF -p PRIMER_BED [OPTIONAL ARGS] FASTQ1 [FASTQ2 ...]
```

For extensive details about each command line argument, see the [Command Line Argument Descriptions](../../wiki/Command-Line-Argument-Descriptions) section of the ViReflow wiki.

## Example Usage
We have provided [demo files](demo), and ViReflow can be executed as follows:

```bash
ViReflow.py -o demo.rf                                                                          `# output Reflow run file` \
            -d s3://my-s3-folder/vireflow-demo                                                  `# output S3 folder` \
            -rf https://github.com/niemasd/ViReflow/raw/main/demo/NC_045512.2.fas               `# reference genome (FASTA)` \
            -rg https://github.com/niemasd/ViReflow/raw/main/demo/NC_045512.2.gff3              `# reference genome annotation (GFF3)` \
            -p https://github.com/niemasd/ViReflow/raw/main/demo/sarscov2_v2_primers_swift.bed  `# primer coordinates file (BED)` \
            https://github.com/niemasd/ViReflow/raw/main/demo/test_R1.fastq                     `# FASTQ 1` \
            https://github.com/niemasd/ViReflow/raw/main/demo/test_R2.fastq                     `# FASTQ 2`
```

This will result in the creation of a file called `demo.rf`, which is the Reflow workflow file. Assuming Reflow is properly installed and configured, the workflow can now be run as follows:

```bash
reflow run demo.rf
```

## Batch Parallel Execution
In a given sequencing experiment, if you have multiple samples you want to run (e.g. `sample1`, `sample2`, ..., `sampleN`), you can use ViReflow to process all of them in parallel (assuming your AWS account has access to spin up sufficient EC2 instances). First, you need to use ViReflow to produce a Reflow run file (`.rf`) for each sample:

```bash
for s in sample1 sample2 [REST_OF_SAMPLES] sampleN ; do ViReflow.py -id $s -o $s.rf [REST_OF_VIREFLOW_ARGS] ; done
```

Alternatively, you can create a CSV file in the following format that, in which the first column contains the run ID, and all remaining columns denote the FASTQ files. You can then run ViReflow as follows to generate the Reflow files for all runs:

|         |                  |                                  |
| ------- | ---------------- | -------------------------------- |
| sample1 | sample1_R1.fastq | s3://my_samples/sample1_R2.fastq |
| sample2 | sample2_R1.fastq | s3://my_samples/sample2_R2.fastq |
| ...     | ...              | ...                              |

```
ViReflow.py [VIREFLOW_ARGS] my_samples.csv
```

Then, you can use the [`rf_batch.py`](rf_batch.py) script to create a batch Reflow run file that will execute all of the individual sample Reflow run files:

```bash
rf_batch.py -o batch_samples.rf sample1.rf sample2.rf [REST_OF_SAMPLES].rf sampleN.rf
```

Now, you can simply run Reflow on the newly-created `batch_samples.rf`, and it will automatically execute of the individual sample Reflow run files:

```bash
reflow run batch_samples.rf
```

# Citing ViReflow
If you use ViReflow in your work, please cite:

> Moshiri N, Fisch KM, Birmingham A, DeHoff P, Yeo GW, Jepsen K, Laurent LC, Knight R (2022). "The ViReflow pipeline enables user friendly large scale viral consensus genome reconstruction." *Scientific Reports*. 12:5077. [doi:10.1038/s41598-022-09035-w](https://doi.org/10.1038/s41598-022-09035-w)

Please also cite the mapper, trimmer, variant caller, and optional analysis tool(s) you used in your ViReflow run(s).
