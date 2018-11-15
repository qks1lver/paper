#!/usr/bin/env python3

import subprocess
import os
import pdb

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
                print('Submitting Trinity job - paired-ends: %s ...' % key)
                command = 'sbatch %s' % p_bash
                _ = subprocess.run(command.split())
                print('\tSubmitted %s' % key)

        return True

    def make_trinity_bash(self, key, fq):

        out_dir = self.out_dir + 'trinity_%s/' % key
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        p = self.out_dir + 'bash_%s.sh' % key
        with open(p, 'w+') as f:
            _ = f.write('#!/usr/bin/bash\n')
            _ = f.write('#SBATCH -J paper\n')
            _ = f.write('#SBATCH -p %s\n' % self.slurm_part)
            _ = f.write('#SBATCH -c %d\n' % self.cpu)
            _ = f.write('#SBATCH --mem=%d\n' % self.mem*1000)
            _ = f.write('module load Java/8 Trinity/2.3.2 Bowtie/2.2.9 Python/3.6.0\n')
            _ = f.write('srun Trinity --seqType fq --max_memory %dG --left %s --right %s --CPU %d --output %s 2>&1\n' % (self.mem, fq['left'], fq['right'], self.cpu, out_dir))

        print('Created bash: %s' % p)
        return p

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
