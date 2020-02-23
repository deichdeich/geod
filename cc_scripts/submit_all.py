import os

for i in range(65):
    os.system(f'qsub geod{i}.pbs')
