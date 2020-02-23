import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-cfp', help = 'config file path', type=str)
parser.add_argument('-icp', help='inish condish path', type=str)
parser.add_argument('-noj', help='number of jobs', type=int)
parser.add_argument('-start', help='start number', type=int)
parser.add_argument('-end', help='end number', type=int)
parser.add_argument('-q', help='which queue', type=str)
opts = parser.parse_args()
for jobnum in range(opts.start, opts.end):
    preamble = f"""#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=1
#PBS -N geod{jobnum}
#PBS -q {opts.q}

cd $PBS_O_WORKDIR

module load gsl/2.5
"""

    cmd = f"""python jobmanager.py -cfp {opts.cfp} -icp {opts.icp} -pr {jobnum}"""

    jobfile = open(f'geod{jobnum}.pbs', 'w')
    jobfile.write(preamble)
    jobfile.write(cmd)
    jobfile.close()
