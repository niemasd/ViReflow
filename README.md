# ViReflow
ViReflow is a tool for constructing elastically-scaling parallelized automated AWS pipelines for viral consensus sequence generation. Given sequence data from a viral sample as well as information about the reference genome and primers, ViReflow generates a [Reflow](https://github.com/grailbio/reflow) file that contains all steps of the workflow, including AWS instance specifications. Because ViReflow is intended to be used with Reflow, the workflows that are developed by ViReflow automatically distribute independent tasks to be run in parallel as well as elastically scale AWS instances based on each individual step of the workflow. ViReflow makes use of compact minimal Docker images for each step of the viral analysis workflow, details about which can be found in the [Niema-Docker](https://github.com/Niema-Docker) GitHub organization.

## Workflow Summary (TODO NEED TO UPDATE, OUT OF DATE)
The workflows produced by ViReflow have the following steps (**^** denotes default choice):
* **Trim the reads** using **^[fastp](https://github.com/OpenGene/fastp)**, [iVar](https://github.com/andersen-lab/ivar), [PRINSEQ](http://prinseq.sourceforge.net/), or [pTrimmer](https://github.com/DMU-lilab/pTrimmer)
* **Map the reads** using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [BWA](http://bio-bwa.sourceforge.net/), or **^[Minimap2](https://github.com/lh3/minimap2)**
* **Generate a pile-up** from the trimmed mapped reads using **[samtools](http://www.htslib.org/)**
* **Call variants** from the pile-up using [freebayes](https://github.com/freebayes/freebayes), [iVar](https://github.com/andersen-lab/ivar), or **^[LoFreq](https://csb5.github.io/lofreq/)*
* **Calculate depth** from the trimmed mapped reads using [samtools](http://www.htslib.org/)
* **Call a consensus sequence** of high-depth regions from the variants using **[bcftools](http://samtools.github.io/bcftools/bcftools.html)**

## Installation
ViReflow is written in Python 3. You can simply download [ViReflow.py](ViReflow.py) to your machine and make it executable:

```bash
wget "https://raw.githubusercontent.com/niemasd/ViReflow/master/ViReflow.py"
chmod a+x ViReflow.py
sudo mv ViReflow.py /usr/local/bin/ViReflow.py # optional step to install globally
```

While ViReflow itself only depends on Python 3, the pipelines it produces are [Reflow](https://github.com/grailbio/reflow) files that run on AWS. Thus, in order to run the pipelines ViReflow produces, one must first [install Reflow](../../wiki/Installing-Reflow).

## Usage (TODO NEED TO UPDATE, OUT OF DATE)
ViReflow can be used as follows:

```
usage: ViReflow.py [-h] -id RUN_ID -d DESTINATION -rf REFERENCE_FASTA -rg
                   REFERENCE_GFF -p PRIMER_BED [-o OUTPUT] [-t THREADS]
                   [-cl COMPRESSION_LEVEL] [--mapped_read_cap MAPPED_READ_CAP]
                   [--min_alt_freq MIN_ALT_FREQ] [--read_mapper READ_MAPPER]
                   [--read_trimmer READ_TRIMMER]
                   [--variant_caller VARIANT_CALLER] [--optional_pangolin]
                   [--optional_spades_coronaspades]
                   [--optional_spades_metaviralspades]
                   [--optional_spades_rnaviralspades] [-u]
                   FQ [FQ ...]

ViReflow: An elastically-scaling parallelized AWS pipeline for viral consensus
sequence generation

positional arguments:
  FQ                    Input FASTQ Files (s3/http/https/ftp; single
                        biological sample)

optional arguments:
  -h, --help            show this help message and exit
  -id RUN_ID, --run_id RUN_ID
                        Unique Run Identifier (for output file naming)
                        (default: None)
  -d DESTINATION, --destination DESTINATION
                        Destination for Results (s3 folder) (default: None)
  -rf REFERENCE_FASTA, --reference_fasta REFERENCE_FASTA
                        Reference Genome Sequence (s3/http/https/ftp to FASTA)
                        (default: None)
  -rg REFERENCE_GFF, --reference_gff REFERENCE_GFF
                        Reference Genome Annotation (s3/http/https/ftp to
                        GFF3) (default: None)
  -p PRIMER_BED, --primer_bed PRIMER_BED
                        Primer (s3/http/https/ftp to BED) (default: None)
  -o OUTPUT, --output OUTPUT
                        Output Reflow File (rf) (default: stdout)
  -t THREADS, --threads THREADS
                        Number of Threads (default: 1)
  -cl COMPRESSION_LEVEL, --compression_level COMPRESSION_LEVEL
                        Compression Level (1 = fastest, 9 = best) (default: 1)
  --mapped_read_cap MAPPED_READ_CAP
                        Successfully-Mapped Read Cap (default: None)
  --min_alt_freq MIN_ALT_FREQ
                        Minimum Alt Allele Frequency for consensus sequence
                        (default: 0.5)
  --read_mapper READ_MAPPER
                        Read Mapper (options: bowtie2, bwa, minimap2)
                        (default: minimap2)
  --read_trimmer READ_TRIMMER
                        Read Trimmer (options: fastp, ivar, prinseq, ptrimmer)
                        (default: ivar)
  --variant_caller VARIANT_CALLER
                        Variant Caller (options: freebayes, ivar, lofreq)
                        (default: lofreq)
  --optional_pangolin   Run Pangolin (optional) (default: False)
  --optional_spades_coronaspades
                        Run SPAdes in coronaSPAdes mode (optional) (default:
                        False)
  --optional_spades_metaviralspades
                        Run SPAdes in metaviralSPAdes mode (optional)
                        (default: False)
  --optional_spades_rnaviralspades
                        Run SPAdes in rnaviralSPAdes mode (optional) (default:
                        False)
  -u, --update          Update ViReflow (current version: 1.0.12) (default:
                        False)
```

For extensive details about each command line argument, see the [Command Line Argument Descriptions](../../wiki/Command-Line-Argument-Descriptions) section of the ViReflow wiki.

## Example Usage
We have provided [demo files](demo), and ViReflow can be executed as follows:

```bash
TODO UPDATE
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

Then, you can use the [`rf_batch.py`](rf_batch.py) script to create a batch Reflow run file that will execute all of the individual sample Reflow run files:

```bash
rf_batch.py -o batch_samples.rf sample1.rf sample2.rf [REST_OF_SAMPLES].rf sampleN.rf
```

Now, you can simply run Reflow on the newly-created `batch_samples.rf`, and it will automatically execute of the individual sample Reflow run files:

```bash
reflow run batch_samples.rf
```
