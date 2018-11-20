# PAPER
**Python Align Pair-End Reads**

Python 3.5+

This is primarily for colleagues at Carnegie Institute for Science, Stanford

##Install
Currently no wheel, sorry! So just clone the repository and navigate to **src/** folder and run.
Essentially follow these steps:
1. Open a terminal and log into a SLURM server (i.e. Memex)
2. navigate to a directory that you are willing to store this tool in
3. Clone repo by typing this in terminal:

    ```git clone https://github.com/qks1lver/paper.git```
    
    This creates a folder called ```paper/``` with all the codes inside

**Now you have it!**

##Trinity
1. Navigate to ```paper/src/```
2. Run something like this:

    ```sbatch run.sh -a trinity -i pair_end_fastq_files_dir/ -o output_files_dir/ -c 12 -m 50```
    
    ```-c``` number of CPUs to allow Trinity to use
    
    ```-m``` amount of RAM in G to allow Trinity to use

## BWA
1. Navigate to ```paper/src/```
2. Run something like this:

    ```sbatch run.sh -a bwa -r reference_genome.fa -i pair_end_fastq_files_dir/ -o output_files_dir/ -c 12 -m 20```
    
    ```-c``` number of CPUs to allow BWA to use
    
    ```-m``` amount of RAM in G to allow BWA to use

## Extract unmapped reads
1. Navigate to ```paper/src/```
2. Run something like this:

    ```sbatch run.sh --unmap -i mapped_read_files_dir/ -o output_files_dir/ -c 4 -m 8```

    ```-c``` number of CPUs to allow SAMtools to use (default is the number of CPUs found.
    On Memex, this would go to 24, which is too much and can potentially make your job wait)
    
    ```-m``` amount of RAM in G to allow SAMtools to use (default is 4)

## Contact
Jiun Yen (jyyen@carnegiescience.edu)
