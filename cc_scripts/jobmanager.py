import numpy
import pygeod
import pickle
import argparse
import time

parser = argparse.ArgumentParser()

parser.add_argument('-cfp', help='the path to the config file', type=str)
parser.add_argument('-icp', help='the path the the i.c. file', type=str)
parser.add_argument('-pr', help='the processor number', type=int)
args = parser.parse_args()
print(args.icp)
pkl_file = open(args.icp,'rb')
ics = pickle.load(pkl_file)
pkl_file.close()
y0 = ics[args.pr][:5]

t0 = time.time()
pygeod.get_geodesic(f'out{args.pr}.bin', init_state=y0, overwrite_header = True, config_file = args.cfp, statelist = args.icp)
print(f"{(time.time()-t0)/60} minutes")
