#!/usr/bin/env python3

# Imports
import subprocess
import os
from shutil import copyfile
import pdb

# Functions
def submit_slurm_bash(p_bash, i_try=1, max_tries=10):

    success = False

    if not p_bash.startswith('sbatch '):
        p_bash = 'sbatch ' + p_bash

    try:
        print('Submitting: %s ...' % p_bash)
        _ = subprocess.run(p_bash.split())
        print('\tSubmitted successfully!')
        success = True
    except:
        print('\tFailed to submit (%d / %d): %s' % (i_try, max_tries, p_bash))
        if i_try < max_tries:
            success = submit_slurm_bash(p_bash, i_try=i_try + 1, max_tries=max_tries)
        else:
            print('\tSubmission failed')

    return success

def write_slurm_bash(p_bash='', commands=None, modules='', p_out='', p_err='', jobname='slurm_job', partition='DPB', cpu=1, mem=1000, use_srun=True, verbose=True):

    if not p_bash:
        raise ValueError('Missing destination path (p_bash=)')

    if commands is None:
        raise ValueError('Missing command string/list (commands=)')
    elif not isinstance(commands, str) and not (isinstance(commands, list) and sum([isinstance(x, str) for x in commands]) == len(commands)):
        raise ValueError('Input for commands= must be either string or list of strings')

    if not p_out:
        p, f = os.path.split(p_bash)
        p += '/slurm_out/'
        if not os.path.isdir(p):
            os.makedirs(p)
        p_out = p + '/slurm_' + f + '.out'

    if not p_err:
        p, f = os.path.split(p_bash)
        p += '/slurm_err/'
        if not os.path.isdir(p):
            os.makedirs(p)
        p_err = p + '/slurm_' + f + '.err'

    if use_srun:
        srun = 'srun '
    else:
        srun = ''

    with open(p_bash, 'w+') as f:
        _ = f.write('#!/usr/bin/bash\n')
        _ = f.write('#SBATCH -J %s\n' % jobname)
        _ = f.write('#SBATCH -p %s\n' % partition)
        _ = f.write('#SBATCH -c %d\n' % cpu)
        _ = f.write('#SBATCH --mem=%d\n' % mem)
        _ = f.write('#SBATCH -o "%s"\n' % p_out)
        _ = f.write('#SBATCH -e "%s"\n' % p_err)
        _ = f.write('module load %s\n' % modules)

        if isinstance(commands, str):
            _ = f.write('%s%s\n' % (srun, commands))
        elif isinstance(commands, list):
            _ = f.write('\n'.join(['%s%s' % (srun, x) for x in commands]))

    if verbose:
        print('Created Slurm bash: %s' % p_bash)

    return p_out

# Classes
class Aligner:

    def __init__(self, target_dir='', out_dir='', reference='', mem=None, cpu=None, slurm_part='DPB'):

        self.target_dir = target_dir
        if not self.target_dir.endswith('/'):
            self.target_dir += '/'
        self.out_dir = out_dir
        if not self.out_dir.endswith('/'):
            self.out_dir += '/'
        self.reference = reference
        if mem is None:
            self.mem = 4
            print('Default with mem=%dG, probably not enough...' % self.mem)
        else:
            self.mem = int(mem)
        if cpu is None:
            self.cpu = os.cpu_count()
            print('Default with cpu=%d (number of CPUs found)' % self.cpu)
        else:
            self.cpu = int(cpu)
        self.slurm_part = slurm_part
        self.max_tries = 10

        self.key2fastqs = {}

    def trinity(self):

        if not self.out_dir:
            print('Need .out_dir - where the results will go')
            return False

        if not self.gen_key2fastqs():
            return False

        if not os.path.isdir(self.out_dir):
            os.makedirs(self.out_dir)
            print('Created folder: %s' % self.out_dir)

        for key, fq in self.key2fastqs.items():
            if 'left' in fq and 'right' in fq:
                # Paired reads
                p_bash = self.make_trinity_bash(key, fq)
                print('Paired-ends: %s ...' % key)
                success = submit_slurm_bash(p_bash, max_tries=self.max_tries)
                if success:
                    print('\tSubmitted %s' % key)

        return True

    def make_trinity_bash(self, key, fq):

        out_dir = self.out_dir + 'trinity_%s/' % key
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        bash_dir = self.out_dir + 'bash/'
        if not os.path.isdir(bash_dir):
            os.makedirs(bash_dir)
        p_bash = bash_dir + 'bash_%s.sh' % key

        _ = write_slurm_bash(
            p_bash=p_bash,
            commands='Trinity --seqType fq --max_memory %dG --left %s --right %s --CPU %d --output %s'
                     % (self.mem, fq['left'], fq['right'], self.cpu, out_dir),
            modules='Java/8 Trinity/2.3.2 Bowtie/2.2.9 Python/3.6.0',
            jobname='paper-trinity',
            partition=self.slurm_part,
            cpu=self.cpu,
            mem=self.mem*1000
        )

        return p_bash

    def bwa(self):

        if not self.out_dir:
            print('Need .out_dir - where the results will go')
            return False

        if not self.reference or not os.path.isfile(self.reference):
            print('Invalid reference genome: %s' % self.reference)
            return False

        if not self.gen_key2fastqs():
            return False

        if not os.path.isdir(self.out_dir):
            os.makedirs(self.out_dir)
            print('Created folder: %s' % self.out_dir)

        p_ref = self.bwa_index()

        for key, fq in self.key2fastqs.items():
            if 'left' in fq and 'right' in fq:
                # Paired reads
                p_bash = self.make_bwa_bash(key, fq, p_ref)
                print('Paired-ends: %s ...' % key)
                success = submit_slurm_bash(p_bash, max_tries=self.max_tries)
                if success:
                    print('\tSubmitted %s' % key)

        return True

    def bwa_index(self):

        print('BWA indexing reference genome ...')

        # Create a directory for reference genome index files
        ref_dir = self.out_dir + 'reference_genome/'
        if not os.path.isdir(ref_dir):
            os.makedirs(ref_dir)

        # Copy the reference genome fasta file to the index file directory
        p_ref = ref_dir + 'reference.fa'
        copyfile(self.reference, p_ref)

        # Run BWA index
        command = 'bwa index ' + p_ref
        subprocess.run(command.split())

        print('\tReference genome index files: %s' % ref_dir)

        return p_ref

    def make_bwa_bash(self, key, fq, p_ref):

        p_out = self.out_dir + 'bwa_%s.bam' % key

        bash_dir = self.out_dir + 'bash/'
        if not os.path.isdir(bash_dir):
            os.makedirs(bash_dir)
        p_bash = bash_dir + 'bash_%s.sh' % key

        _ = write_slurm_bash(
            p_bash=p_bash,
            commands='bwa mem -M -t %d %s %s %s | samtools sort -@%d -o %s -'
                     % (self.cpu, p_ref, fq['left'], fq['right'], self.cpu, p_out),
            modules='Python/3.6.0 Zlib/1.2.8 bwa/0.7.15 SAMtools/1.9',
            jobname='paper-bwa',
            partition=self.slurm_part,
            cpu=self.cpu,
            mem=self.mem * 1000,
            use_srun=False
        )

        return p_bash

    def gen_key2fastqs(self):

        if not self.target_dir:
            print('Need .target_dir - where all the .fastq.gz or .fastq files are stored')
            return False

        print('Looking for FASTQ files in %s ...' % self.target_dir)

        self.key2fastqs = {}
        p_files = os.listdir(self.target_dir)
        for p_file in p_files:
            if self.check_fastq(p_file):
                f_dir, f_name, key, direction = self.parse_fastq(p_file)
                if not f_dir:
                    f_dir = self.target_dir
                if not f_dir.endswith('/'):
                    f_dir += '/'
                f_full_name = f_dir + f_name
                if key not in self.key2fastqs:
                    self.key2fastqs[key] = {direction:f_full_name}
                else:
                    self.key2fastqs[key][direction] = f_full_name

                print('Added %s to dictionary' % f_name)

        return True

    @staticmethod
    def parse_fastq(p_file, left_tag='_1_pf', right_tag='_2_pf'):

        f_dir, f_name = os.path.split(p_file)

        key = f_name.replace('.gz', '')
        key = key.replace('.fastq', '')

        if left_tag in key:
            key = key.replace(left_tag, '')
            return f_dir, f_name, key, 'left'

        if right_tag in key:
            key = key.replace(right_tag, '')
            return f_dir, f_name, key, 'right'

    @staticmethod
    def check_fastq(p_file=''):

        if p_file.endswith('.fastq.gz') or p_file.endswith('.fastq'):
            return True
        else:
            return False

class Analyzer:

    def __init__(self, target_dir='', out_dir='', cpu=None, mem=None, slurm_part='DPB'):

        self.target_dir = target_dir
        if not self.target_dir.endswith('/'):
            self.target_dir += '/'
        self.out_dir = out_dir
        if not self.out_dir.endswith('/'):
            self.out_dir += '/'
        if cpu is None:
            self.cpu = os.cpu_count()
            print('Default with cpu=%d (number of CPUs found)' % self.cpu)
        else:
            self.cpu = int(cpu)
        if mem is None:
            self.mem = 4
            print('Default with mem=%dG' % self.mem)
        else:
            self.mem = int(mem)
        self.slurm_part = slurm_part

        self.max_tries = 10

    def get_unmapped(self):

        p_files = os.listdir(self.target_dir)
        if not p_files:
            print('Source directory is empty: %s' % self.target_dir)
            return False

        if not os.path.isdir(self.out_dir):
            os.makedirs(self.out_dir)

        for p_file in p_files:
            if (p_file.endswith('.bam') or p_file.endswith('.sam')) and os.path.isfile(self.target_dir + p_file):
                p_bash = self.make_sam_bash(p_bam=p_file)
                submit_slurm_bash(p_bash=p_bash, max_tries=self.max_tries)

        return True

    def make_sam_bash(self, p_bam=''):

        if not p_bam:
            raise ValueError('Missing BAM/SAM file path (p_bam=)')

        opt_sam = ''
        if not p_bam.endswith('.bam'):
            if p_bam.endswith('.sam'):
                opt_sam = ' -S'
            else:
                raise ValueError('p_bam must be either .bam or .sam')

        p_in = self.target_dir + p_bam
        if not os.path.isfile(p_in):
            raise ValueError('File does not exist: %s' % p_in)

        bash_dir = self.out_dir + 'bash/'
        if not os.path.isdir(bash_dir):
            os.makedirs(bash_dir)
        if not opt_sam:
            p_out = self.out_dir + p_bam.replace('.bam', '_unmapped.sam')
            p_bash = bash_dir + 'bash_' + p_bam.replace('.bam', '.sh')
        else:
            p_out = self.out_dir + p_bam.replace('.sam', '_unmapped.sam')
            p_bash = bash_dir + 'bash_' + p_bam.replace('.sam', '.sh')

        _ = write_slurm_bash(
            p_bash=p_bash,
            commands='samtools view%s -f12 %s > %s' % (opt_sam, p_in, p_out),
            modules='Python/3.6.0 SAMtools/1.9',
            jobname='paper-sam',
            partition=self.slurm_part,
            cpu=self.cpu,
            mem=self.mem * 1000,
            use_srun=False
        )

        return p_bash
    