#!/usr/bin/env python3
import os,sys
import argparse
import numpy as np
import networkx as nx
from os import path
from itertools import repeat
import concurrent.futures as cf
#from netdist_functions import *
from netdist_functions import get_lcc,sif2nx,labels2idx,load_nodes,get_degree_binning,nodestoh5, random_nodes_from,get_sp,calc_max_ccv,calculate_distance,sp_ij

parser=argparse.ArgumentParser(description='Calculates the shortest paths among all nodes of a network and saves them in sparse format along with their "k" index in a file.')
parser.add_argument('network_file', help='A network file (interactome) in ".sif" format')
parser.add_argument('--no_lcc', default=False, action='store_true',  help='Consider the whole network and not only the largest connected component.')
parser.add_argument('--n_cpus', type=int, default=4, help='To speed-up some calculations, the script tries to parallelize critical steps. Set to zero to disable this behaviour')
parser.add_argument('--force', '-f', default=False, action='store_true', help='Force execution ignoring requirements')

args = parser.parse_args()

file_name = args.network_file
only_lcc = not args.no_lcc
ncpus,maxcpus = args.n_cpus,len(os.sched_getaffinity(0))

network = sif2nx(file_name)
if only_lcc and not file_name.endswith('.lcc'):
 file_name = get_lcc(network,file_name)
 #volverla a cargar es la forma mas rapida de que siempre esten en el mismo orden los nodos...
 network = sif2nx(file_name)

all_spf, indexf = file_name + "_sp.npy",file_name + "_idx.npy"
n = len(network)
N = n*(n-1) // 2 #sum(range(n)) #((n*n)-n)//2

if not path.isfile(all_spf):
 print(f"File with shortest_paths ('{all_spf}') not found. Calculating...", end=' ',flush=True)
 index = list(network.nodes)
 N_n, N_network, N_index = repeat(n), repeat(network), repeat(index)
 if N > 1e7 and ncpus < 12:
  if not args.force: sys.exit(f"W: A total of '{N}' distances are to be calculated.\nAdjust value of 'n_cpus' (up to: {maxcpus}) or use '--force' argument.")
 sp = np.zeros(N,dtype=float)
 #Nncpus = 40
 if (chunks:=int(N/(ncpus-1))) > 1e6: chunks = int(1e6)
 with cf.ProcessPoolExecutor(ncpus) as executor:
  futures = executor.map(sp_ij,range(N),N_n,N_network,N_index,chunksize=chunks)
  #[future for future in futures]
  for i,future in enumerate(futures): sp[i] = future
 print('DONE')
 np.save(indexf,index)
 np.save(all_spf,sp)
else:
 print(f"File with shortest paths already exists: {all_spf}")
