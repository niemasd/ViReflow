# ViReflow
ViReflow is a tool for constructing elastically-scaling parallelized automated AWS pipelines for viral consensus sequence generation. Given sequence data from a viral sample as well as information about the reference genome and primers, ViReflow generates a [Reflow](https://github.com/grailbio/reflow) file that contains all steps of the workflow, including AWS instance specifications. Because ViReflow is intended to be used with Reflow, the workflows that are developed by ViReflow automatically distribute independent tasks to be run in parallel as well as elastically scale AWS instances based on each individual step of the workflow. ViReflow makes use of compact minimal Docker images for each step of the viral analysis workflow, details about which can be found in the [Niema-Docker](https://github.com/Niema-Docker) GitHub organization.

## Workflow Summary
The workflows produced by ViReflow have the following steps:
* **Map the reads** using [Minimap2](https://github.com/lh3/minimap2)
    * Uses the [`niemasd/minimap2_samtools`](https://hub.docker.com/repository/docker/niemasd/minimap2_samtools) Docker image
* **Trim the mapped reads** using [iVar](https://github.com/andersen-lab/ivar)
    * Uses the [`niemasd/ivar`](https://hub.docker.com/repository/docker/niemasd/ivar) Docker image
* **Generate a pile-up** from the trimmed mapped reads using [samtools](http://www.htslib.org/)
    * Uses the [`niemasd/samtools`](https://hub.docker.com/repository/docker/niemasd/samtools) Docker image
* **Call variants** from the pile-up using [iVar](https://github.com/andersen-lab/ivar)
    * Uses the [`niemasd/ivar`](https://hub.docker.com/repository/docker/niemasd/ivar) Docker image
* **Call a consensus sequence** from the pile-up using [iVar](https://github.com/andersen-lab/ivar)
    * Uses the [`niemasd/ivar`](https://hub.docker.com/repository/docker/niemasd/ivar) Docker image
* **Calculate depth** from the trimmed mapped reads using [samtools](http://www.htslib.org/) (optional, recommended)
    * Uses the [`niemasd/samtools`](https://hub.docker.com/repository/docker/niemasd/samtools) Docker image

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
usage: ViReflow.py [-h] -id RUN_ID -d DESTINATION -rf REFERENCE_FASTA [-rm REFERENCE_MMI] -rg REFERENCE_GFF -p PRIMER_BED [-o OUTPUT] [-mt MAX_THREADS] [--include_depth] [-u] FQ [FQ ...]

flag arguments:
  -h, --help                                               show this help message and exit
  -id RUN_ID, --run_id RUN_ID                              Unique Run Identifier (for output file naming)
  -d DESTINATION, --destination DESTINATION                Destination for Results (s3 folder)
  -rf REFERENCE_FASTA, --reference_fasta REFERENCE_FASTA   Reference Genome Sequence (s3/http/https/ftp to FASTA)
  -rm REFERENCE_MMI, --reference_mmi REFERENCE_MMI         Reference Genome Minimap2 Index (s3 to MMI) (default: None)
  -rg REFERENCE_GFF, --reference_gff REFERENCE_GFF         Reference Genome Annotation (s3/http/https/ftp to GFF3)
  -p PRIMER_BED, --primer_bed PRIMER_BED                   Primer (s3/http/https/ftp to BED)
  -o OUTPUT, --output OUTPUT                               Output Reflow File (rf) (default: stdout)
  -mt MAX_THREADS, --max_threads MAX_THREADS               Max Threads (default: 32)
  --include_depth                                          Include Depth Calling (default: False)
  -u, --update                                             Update ViReflow (default: False)

positional arguments:
  FQ                    Input FASTQ Files (s3 paths; single biological sample)
```

For extensive details about each command line argument, see the [Command Line Argument Descriptions](../../wiki/Command-Line-Argument-Descriptions) section of the ViReflow wiki.

## Example Usage
We have provided [demo files](demo), and ViReflow can be executed as follows:

```bash
ViReflow.py -id demo                                            \ # specify run ID
            -rf s3://vireflow-demo/NC_045512.2.fas              \ # specify reference genome sequence
            -rm s3://vireflow-demo/NC_045512.2.fas.mmi          \ # specify reference genome Minimap2 index
            -rg s3://vireflow-demo/NC_045512.2.gff3             \ # specify reference genome annotation
            -p s3://vireflow-demo/sarscov2_v2_primers_swift.bed \ # specify primer BED
            -d s3://vireflow-demo/results_all_s3                \ # specify result destination directory
            --include_depth                                     \ # include depth calculation
            -o all_s3.rf                                        \ # specify output Reflow file
            -mt 2                                               \ # use at most 2 threads in any given step
            s3://vireflow-demo/test_R1.fastq                    \ # first FASTQ file
            s3://vireflow-demo/test_R2.fastq                      # second FASTQ file
```

This will result in the creation of a file called `all_s3.rf`, which is the Reflow workflow file. Assuming Reflow is properly installed and configured, the workflow can now be run as follows:

```bash
reflow run all_s3.rf
```

## Batch Parallel Execution
In a given sequencing experiment, if you have multiple samples you want to run (e.g. `sample1`, `sample2`, ..., `sampleN`), you can use ViReflow to process all of them in parallel (assuming your AWS account has access to spin up sufficient EC2 instances). First, you need to use ViReflow to produce a Reflow run file (`.rf`) for each sample:

```bash
for s in sample1 sample2 [REST_OF_SAMPLES] sampleN ; do ViReflow.py -id $s -o $s.rf [REST_OF_VIREFLOW_ARGS] ; done
```

Then, you can use the [`rf_batch.py`](rf_batch.py) script to create a batch Reflow run file (`.rf`) that will execute all of the individual sample Reflow run files:

```bash
rf_batch.py -o batch_samples.rf sample1.rf sample2.rf [REST_OF_SAMPLES].rf sampleN.rf
```

Now, you can simply run Reflow on the newly-created `batch_samples.rf`, and it will automatically execute of the individual sample Reflow run files:

```bash
reflow run batch_samples.rf
```
