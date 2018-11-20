# PAPER
**Python Aligns Paired-End Reads**

Python 3.5+

This is primarily for colleagues at Carnegie Institute for Science, Stanford

**This is meant to be run on a system managed with SLURM !!**

## Install
Currently no wheel, sorry! So just clone the repository and navigate to **src/** folder and run.
Essentially follow these steps:
1. Open a terminal and log into a server (i.e. Memex)
2. navigate to a directory that you are willing to store this tool in
3. Clone repo by typing this in terminal:

    ```git clone https://github.com/qks1lver/paper.git```
    
    This creates a folder called ```paper/``` with all the codes inside

**Now you have it!**

## Trinity
1. Navigate to ```paper/src/```
2. Run something like this:

    ```sbatch run.sh -a trinity -i pair_end_fastq_files_dir/ -o output_files_dir/ -c 12 -m 50```
    
    ```-i``` files in this directory can be .fastq or .fastq.gz
    
    ```-o``` directory of output files from Trinity, will creat folder if non-existent
    
    ```-c``` number of CPUs to allow Trinity to use
    
    ```-m``` amount of RAM in G to allow Trinity to use

## BWA
1. Navigate to ```paper/src/```
2. Run something like this:

    ```sbatch run.sh -a bwa -r reference_genome.fa -i pair_end_fastq_files_dir/ -o output_files_dir/ -c 12 -m 20```
    
    ```-r``` FASTA file reference genome
    
    ```-i``` files in this directory can be .fastq or .fastq.gz
    
    ```-o``` directory of output BAM files (fastq -> BWA -> SAMtools -> BAM), will creat folder if non-existent
    
    ```-c``` number of CPUs to allow BWA to use
    
    ```-m``` amount of RAM in G to allow BWA to use

## Extract unmapped reads
1. Navigate to ```paper/src/```
2. Run something like this:

    ```sbatch run.sh --unmap -i mapped_read_files_dir/ -o output_files_dir/ -c 4 -m 8```

    ```-i``` files in this directory can be BAM or SAM
    
    ```-o``` directory of output BAM files of unmapped reads, will creat folder if non-existent
    
    ```-c``` number of CPUs to allow SAMtools to use (default is the number of CPUs found.
    On Memex, this would go to 24, which is too much and can potentially make your job wait)
    
    ```-m``` amount of RAM in G to allow SAMtools to use (default is 4)

## Unpair paired BAMs
This could come after extracting the unmapped reads, which would generate one bam files with
both pairs of reads. Unpair the reads, so they can then go through Trinity.
1. Navigate to ```paper/src/```
2. Run something like this:

    ```sbatch run.sh --unpair -i unmapped_read_files_dir/ -o output_files_dir/ -c 4 -m 8```

    ```-i``` files in this directory can be BAM or SAM
    
    ```-o``` directory of output .fastq.gz files of unpaired unmapped reads, will creat folder if non-existent
    
    ```-c``` number of CPUs to allow SAMtools to use (default is the number of CPUs found.
    On Memex, this would go to 24, which is too much and can potentially make your job wait)
    
    ```-m``` amount of RAM in G to allow SAMtools to use (default is 4)

## Check jobs
To maximize parallelization, **Paper** submits one job for each file (or each set of paired-end fastq files).
Each job uses as many processors and RAM as specified by ```-c``` and ```-m```.

To check the progress of jobs, use SLURM's ```squeue```, for example, my account name on Memex is jyyen:

```squeue -u jyyen```

For more on SLURM, check out this great list of the most used commands:

https://www.rc.fas.harvard.edu/resources/documentation/convenient-slurm-commands/

## Contact
Jiun Yen (jyyen@carnegiescience.edu)
