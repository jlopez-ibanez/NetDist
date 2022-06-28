#!/usr/bin/env python3
import os,sys
import argparse
import numpy as np
import networkx as nx
from os import path
from itertools import repeat
import concurrent.futures as cf
#from netdist_functions import *
from netdist_functions import get_lcc,sif2nx,labels2idx,load_nodes,get_degree_binning,nodesfromh5, random_nodes_from,get_sp,calc_max_ccv,calculate_distance

methods_list = ['closest','separation','shortest','kernel','iclosest','ikernel','center']
avm="".join(map("'{}',".format,methods_list)).strip(',')
pydir = path.curdir if path.dirname(sys.argv[0]) == '' else path.dirname(sys.argv[0])

parser=argparse.ArgumentParser(description='Calculates a series of network distances among one or two sets of nodes.')
parser.add_argument('network_file', help='A network file (interactome) in ".sif" format')
parser.add_argument('origin_nodes', help='Nodes FROM which calculate the distance.')
parser.add_argument('destiny_nodes', help='Nodes TO which calculate the distance.')
parser.add_argument('--method','-m', choices=methods_list, help=f'One of the available methods for calculating distances: {avm}')
parser.add_argument('--n_random','-nr', type=int, default=1000, help='Perform this number of randomizations when calculating the null distribution.')
parser.add_argument('--seed','-s', type=int, default=3285, help='A seed to generate reproducible results.')
parser.add_argument('--no_lcc', action='store_true', default=False, help='Consider the whole network and not only the largest connected component.')
parser.add_argument('--bin_size', type=int, default=2, help='Size of the bins when doing degree-aware randomization. Default is twice the maximum number of nodes in input sets.')
parser.add_argument('--n_cpus', type=int, default=4, help='To speed-up some calculations, the script tries to parallelize critical steps. Set to zero to disable this behaviour')

args = parser.parse_args()

file_name,nodes_fromf,nodes_tof = args.network_file,args.origin_nodes,args.destiny_nodes
method = args.method
only_lcc,n_random,seed = not args.no_lcc,args.n_random,args.seed
ncpus = args.n_cpus

network = sif2nx(file_name)
if only_lcc and not file_name.endswith('.lcc'):
 file_name = get_lcc(network,file_name)
 network = sif2nx(file_name)

all_spf,indexf = file_name+"_sp.npy",file_name+"_idx.npy"
n = len(network)
N = n*(n-1)//2

nfrom,nodes_from = load_nodes(nodes_fromf,network)
if nfrom==0: sys.exit(nodes_from)
nto,nodes_to = load_nodes(nodes_tof,network)
if nto == 0: sys.exit(nodes_to)
ncommon = len(set(nodes_from).intersection(nodes_to))
nf_all,nt_all = load_nodes(nodes_fromf)[0],load_nodes(nodes_tof)[0]

bin_size = 2*max((nfrom,nto)) if args.bin_size == 2 else args.bin_size
bins = get_degree_binning(network, bin_size, None)

if not path.isfile(all_spf):
 print(f"WARNING: File with shortest_paths ('{all_spf}') not found. Depending on the size of the network, calculations may take some time...\n",file=sys.stderr)
 index,sp = None,network
 nf_random = random_nodes_from(network, nodes_from, n_random, bins, index, seed)
 nt_random = random_nodes_from(network, nodes_to, n_random, bins, index, seed)
 
else:
 print(f"Loading distances and indexes of nodes from: '{indexf}' '{all_spf}'...", file=sys.stderr)
 index,sp = np.load(indexf),np.load(all_spf)
 index = dict((label,i) for i,label in enumerate(index))
 nf_random = random_nodes_from(network, nodes_from, n_random, bins, index, seed)
 nt_random = random_nodes_from(network, nodes_to, n_random, bins, index, seed)
 nodes_to = labels2idx(nodes_to,index)
 nodes_from = labels2idx(nodes_from,index)

setids = '_'.join([path.basename(nodes_fromf),path.basename(nodes_tof)])
source,target = [path.basename(nodes_fromf),path.basename(nodes_tof)]
header = ["setid","distance","zdistance","mean",'std','p-value','method','nfrom','nto','common']

if method == 'center':
 center_nodes = calc_max_ccv(nodes_to,sp,n)
 d = calculate_distance(nodes_from,center_nodes,sp,n,method)
 with cf.ProcessPoolExecutor(ncpus) as executor:
  futures = executor.map(calc_max_ccv, nt_random,repeat(sp),repeat(n), chunksize=n_random//ncpus)
  nt_random = [future for future in futures]
 print("DONE",file=sys.stderr)
else:
 d = calculate_distance(nodes_from,nodes_to,sp,n,method)

values=np.zeros(n_random,dtype=float)
if n_random*nfrom > 1e5 or n_random > 1000:
 optchunk = n_random//ncpus
 N_sp,N_n,N_method = repeat(sp),repeat(n),repeat(method)
 with cf.ProcessPoolExecutor(ncpus) as executor:
  results = executor.map(calculate_distance,nf_random,nt_random,N_sp,N_n,N_method,chunksize=optchunk)
  for i,result in enumerate(results):
   values[i] = result
else:
 for i,(nf_random,nt_random) in enumerate(zip(nf_random,nt_random)):
  values[i] = calculate_distance(nf_random, nt_random,sp,n,method)

pval = (values <= d).mean() # needs high number of n_random
m, s = np.mean(values), np.std(values)
z = 0.0 if s == 0 else (d - m) / s
print(*header,sep='\t')
print(setids,d, z, m, s,pval,method,nfrom,nto,ncommon,sep='\t')

