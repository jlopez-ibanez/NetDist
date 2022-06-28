if __name__ == "__main__":
 print("Functions used with the NetDist scripts")
import sys,os
import gzip
import h5py
import random
from os import path
from math import sqrt
import numpy as np
import networkx as nx
#from itertools import repeat
#import concurrent.futures as cf


##### Random Node Choice #####
def get_degree_binning(g, bin_size, lengths=None):
    degree_to_nodes = {}
    for node, degree in g.degree(): #.items(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = list(degree_to_nodes.keys())
    values.sort()
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        #Here BOTH LISTS ARE NOW LINKED. The val.extend() below WILL ALSO UPDATE degree_to_nodes[values[i]] (python will keep this reference despite modifying 'i' in the next lines
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            #!!Using val.extend() updates also the list of nodes of degree 'i' declared in the outer loop. DOESN'T AFFECT THE FINAL RESULT as I won't take that degree again.
            val.extend(degree_to_nodes[values[i]])
            #val=val+degree_to_nodes[values[i]]
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1 
        #print(i, low, high, len(val))
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins

def get_degree_equivalents(seeds, bins, g):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        #Dif not seed in g:sys.exit(f"{seed}")
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes

def random_nodes_from(network,nodes_selected,n_random,bins, node2idx=None, seed=None, degree_aware=True, connected=False):
 """
    Modified from pick_random_nodes_matching_selected() for returning idx of labels instead of labels. Removed also conditions for "degree_aware = False" and for "degree_aware == True and connected"
 """
 if seed is not None: random.seed(seed)
 values = []
 nodes = network.nodes()
 for i in range(n_random):
  if degree_aware and connected: raise ValueError("Not implemented!")
  elif degree_aware:
   nodes_random = set()
   node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
   for node, equivalent_nodes in node_to_equivalent_nodes.items():
    #nodes_random.append(random.choice(equivalent_nodes))
    chosen = random.choice(equivalent_nodes)
    for k in range(20): # Try to find a distinct node (at most 20 times)
     if chosen in nodes_random:  chosen = random.choice(equivalent_nodes)
    nodes_random.add(chosen)
  #si no es degree_aware Y connected ni degree_aware sera normal, connected o no
  elif connected:
   nodes_random = [ random.choice(nodes) ]
   k = 1
   while True:
    if k == len(nodes_selected): break
    node_random = random.choice(nodes_random)
    node_selected = random.choice(network.neighbors(node_random))
    if node_selected in nodes_random: continue
    nodes_random.append(node_selected)
    k += 1
  else:
   nodes_random = random.sample(nodes, len(nodes_selected))
  nodes_random = list(nodes_random)
  if node2idx is None:  values.append(nodes_random)
  else: values.append([node2idx[label] for label in nodes_random])
 return values

##### END Random Node Choice #####

##### Network related #####

def labels2idx(nodes,index,invert=False):
 if invert: index=dict((v,k) for k,v in index.items())
 return [index[node] for node in nodes]

def get_lcc(network,file_name,out_format='sif'):
 lcc=max(nx.connected_components(network),key=len)
 network = network.subgraph(lcc)
 print(file_name,file=sys.stderr)
 if file_name.endswith('.lcc'):
  print('USING PATH',file=sys.stderr)
  basename = path.splitext(file_name)[0]
 if out_format == 'sif':
  file_name = path.splitext(file_name)[0]+'.sif.lcc'
  rows = [(u,'1',v) for u,v in network.edges()]
 elif not out_format is None:
  #this requires to check the format for loading the file again...(?)
  sys.exit("ERROR: Saving to another format not implemented...")
  file_name = path.splitext(file_name)[0]+f'_edges.lcc'
  rows = [(u,v) for u,v in network.edges()]
 with open(file_name,'w') as outf:
  print(*map(' '.join,rows),sep='\n', file=outf)
 print(f"Saved '{file_name}'",file=sys.stderr)
 return file_name

##### END Network related #####

#### Load Stuff ####
def sif2nx(siffile):
  """
  Load sif file as networkx network
  default is to (try to) convert values to float
  use TAB as delim in sif file if you are expecting whitespaces in names
  by definition, the middle column should always exist BUT I may not want those values...
  nodes dictionary will contain all unconnected nodes...
  """
  nodes,edges = {},{}
  with open(siffile) as inpf:
    for line in inpf:
      if line[0] == '#': continue
      l=line.strip().split()
      if len(l)>=3:
        n1,value,n2=l[:3]
        if not value is None:
          try: value=float(value)
          #print("Value for node/edge '{}' can't be interpreted as a number")
          except ValueError: pass
        edges[(n1,n2)]=value
      elif len(l)>0:
        n1,value=(l[0],None) if len(l)==1 else l[:2]
        if not value is None:
          try: value=float(value)
          except ValueError: pass
        nodes[n1]=value
  network=nx.Graph()
  if len(edges)>0:
    if all((v is None for v in edges.values())): network.add_edges_from(edges.keys())
    else: network.add_weighted_edges_from((nodes+(v,) for nodes,v in edges.items()),weight='weight')
      #for edge,w in edges.items(): network.add_edge(*edge),weight=w)
  if len(nodes)>0:
    if all((v is None for v in nodes.values())):network.add_nodes_from(nodes.keys())
    else:network.add_nodes_from(((node,{'value':v}) for node,v in nodes.items()))
      #for node,v in nodes.items(): network.add_node(node,value=v)
  return network

def edges2nx(edgef, w=False, attr=None,attr_label='label'):
 #to be done: checking if weights are correctly loaded and/or need conversion to number 
 #rewritten from the one used in heatdiffusion to make it simpler
 if isinstance(attr,list): sys.exit("ERROR! List of attributes not supported")
 filelist,weights=[],[]
 try:
  with open(edgef) as inpf:
   for line in inpf:
    if line[0] == '#': continue
    filelist.append(line.strip().split())
 except IOError as e:
  if not path.isfile(edgef): sys.exit("ERROR! Couldn't locate '{}'.".format(edgef))
  else: sys.exit("ERROR! Couldn't open '{}'".format(edgef))
 ncols = set(map(len,filelist))
 if min(ncols)<2: sys.exit("ERROR! File format error. Expected at least two columns in every row")
 if w and len(ncols)>1: sys.exit("ERROR! Found rows with different number of columns. (Not supported when loading weights)")
 if w and not (mncols:=max(ncols))>2 : sys.exit(f"ERROR! Expected at least three columns to load weights, found: {mncols}")
 filecols = list(zip(*filelist))
 G = nx.Graph()
 if w:
  s,t,weights = filecols[:3]
  G.add_weighted_edges_from(edges)
 else:
  s,t = filecols[:2]
  G.add_edges_from(zip(s,t))
 if not attr is None:
  G.add_nodes_from(set((node for edge in edges for node in edge)))
  for node in G.nodes(): G.nodes[node][attr_label]=attr
 return G

def load_nodes(nodes_listf ,network_nodes=[]):
 if not path.isfile(nodes_listf): return(0, f"ERROR! Couldn't locate '{nodes_listf}'")
 nodes_list = [__.strip() for __ in open(nodes_listf).readlines()]
 if not network_nodes==[]:
  nodes_list = sorted(list(set(nodes_list).intersection(network_nodes)))
 else:
  nodes_list.sort()
 n = len(nodes_list)
 if n == 0 and network_nodes == []: return(0,f"ERROR! Not found any node in '{nodes_listf}'")
 elif n == 0: return(0, f"ERROR! Not found any node from '{nodes_listf}' in the network")
 return (n,nodes_list)

def nodesfromh5(nodes_fromf,nodes_tof,network_file,seed=None,idx=False):
 nameset = '_'.join([path.basename(nf) for nf in (nodes_fromf,nodes_tof)])
 if network_file.endswith('.lcc'): network_file = path.splitext(network_file)[0]
 h5file = f"{path.splitext(network_file)[0]}-{nameset}"
 if not seed is None: h5file += "-"+str(seed)
 h5file = path.join(path.dirname(nodes_fromf),h5file)+'.h5.gz'
 if path.isfile(h5file):
  try:
   with gzip.open(h5file,'rb') as inpf:
    with h5py.File(inpf,'r') as hf:
     nodes = [np.array(hf[grpname+'/idx']) if idx else np.array(hf[grpname+'/lbl']).astype(str) for grpname in ['nodes_from','nodes_to']]
  except Exception as e: sys.exit(f"ERROR while opening '{h5file}':\n{e}")
 else:
  h5file = path.splitext(h5file)[0]
  try:
   with h5py.File(h5file,'r') as hf:
    nodes = [np.array(hf[grpname+'/idx']) if idx else np.array(hf[grpname+'/lbl']).astype(str) for grpname in ['nodes_from','nodes_to']]
  except Exception as e: sys.exit(f"ERRROR while opening '{h5file}':\n {e}")
 return nodes

def nodestoh5(nodesets,setsfiles,network_file,seed=None,compressed=True):
 if not isinstance(nodesets,tuple) or not isinstance(setsfiles,tuple): sys.exit("E! Arguments not passed as tuple!")
 nameset='_'.join([path.basename(nf) for nf in setsfiles])
 if network_file.endswith('.lcc'): network_file = path.splitext(network_file)[0]
 h5file = f"{path.splitext(network_file)[0]}-{nameset}"
 if not seed is None: h5file+="-"+str(seed)
 h5file =path.join(path.dirname(setsfiles[0]),h5file)+'.h5'
 hf = h5py.File(h5file,'w')
 #hf.attrs.create('nname',network_file, dtype=h5py.string_dtype())
 #hf.attrs.create('seed',str(seed), dtype=h5py.string_dtype())
 try:
  for grpname,nodeset in zip(['nodes_from','nodes_to'],nodesets):
   g1 = hf.create_group(grpname)
   for i,(dsetname,dsettype) in enumerate(zip(['lbl','idx'], [h5py.string_dtype(),int])):
    datos = np.stack(nodeset[i]).astype(dsettype)
    g1.create_dataset(dsetname, data=datos, dtype=dsettype)
    #g1 = hf.create_group('labels')
    #g1.create_dataset(f"r{i+1}",data=nodes1[::-1], dtype=h5py.string_dtype())
 except Exception as e:
  #hf.close()
  sys.exit(f"E! Creation of {h5file} failed: {e}")
 hf.close()
 if compressed:
  x = run(['gzip','-f',h5file])
  if x.returncode == 0: h5file += '.gz'
 return h5file

#### END Load Stuff ####

#### Indexes manipulation ####
def sp_ij(k,n,G,nodes):
 i = int(n - 2 - int(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5))
 j = int(k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2)
 return nx.shortest_path_length(G,nodes[i],nodes[j])

def ij2k(i,j,n):
 #np.unravel_index(k,(n,n))
 if i > j: i,j = j,i
 return int((n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1)

def k2ij(k,n):
 #https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
 #https://stackoverflow.com/questions/19143657/linear-indexing-in-symmetric-matrices
 #np.unravel_index(k,(n,n))
 i = int(n - 2 - int(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5))
 j = int(k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2)
 return(i,j)
#### END Indexes manipulation ####

#### DistanceCalculations ####
def get_sp(nodes_from,nodes_to,network,network_size,not_self=False):
 self_value = np.nan if not_self else 0
 sp = np.zeros((len(nodes_from),len(nodes_to)),dtype=float)
 if len(network) == network_size:
  for i,nf in enumerate(nodes_from):
   for j,nt in enumerate(nodes_to):
    sp[i,j] = self_value if nf == nt else nx.shortest_path_length(network,nf,nt)
 else: #network === sp
  for i,nf in enumerate(nodes_from):
   for j,nt in enumerate(nodes_to):
    sp[i,j] = self_value if nf == nt else network[ij2k(nf,nt,network_size)]
 return sp


def calc_max_ccv(nodes,network,n, return_values=False):
 #If length of second argument (network) is the same as third (n) this calculates shortest paths. If not, assumes second argument is the list of precalculated shortest paths.
 #NOTE that in either case 'nodes' argument would differ (labels for calculating SP and index of labels for retrieving the precalculated SP) 
 nsp = len(network)
 center,center_d,center_values = None,None,None
 for nf in nodes:
  values = np.zeros(len(nodes),dtype=float)
  for i,nt in enumerate(nodes):
   if nsp == n:
    values[i] = 0 if nf == nt else nx.shortest_path_length(network,nf,nt)
   else:
    values[i] = 0 if nf == nt else network[ij2k(nf,nt,n)]
  d = values.sum()
  if center_d is None or d < center_d:
   center = [ nf ]
   center_d,center_values = d,values
   #center_values = [ values ]
  elif d == center_d:
   center.append(nf)
   center_values += values
   #center_values.append(values)
 center_values /= len(center)
 return (center,center_values) if return_values else center

def max_ccv(i_nodes,lengths,n,return_values=False):
 #Lo mismo que el codigo de Guney pero cambiando algunas cosas para hacerlo más compacto 
 #En este "i_nodes" son LOS INDICES DE LOS NODOS que quiero calcular su centralidad en la matriz de distancias
 center,center_d,center_values = None,None,None
 for nf in i_nodes:
  values = np.zeros(len(i_nodes),dtype=float)
  for i,nt in enumerate(i_nodes):
   values[i] = 0 if nf == nt else lengths[ij2k(nf,nt,n)]
  d = values.sum()
  if center_d is None or d < center_d:
   center = [ nf ]
   center_d,center_values = d,values
   #center_values = [ values ]
  elif d == center_d:
   center.append(nf)
   center_values += values
   #center_values.append(values)
 center_values /= len(center)
 if return_values: return (center,center_values)
 else: return center

def calculate_distance(nodes_from,nodes_to,network,network_size,method='closest'):
 if method[0]=='i':
  method = method[1:]
  nodes_from,nodes_to = nodes_to,nodes_from
 lengths_to = get_sp(nodes_from,nodes_to,network,network_size)
 if method == 'separation':
  '''
  Separation metric proposed by Menche et al. 2015 => DOI: 10.1126/science.1257601
  '''
  #dAB => dispersion
  dAB = np.concatenate((lengths_to.min(1),lengths_to.min(0))).mean()
  inner_ds = []
  for nodes in (nodes_from,nodes_to):
   lengths_to = get_sp(nodes,nodes,network,network_size,not_self=True)
   masked_lt = np.ma.masked_array(lengths_to,np.isnan(lengths_to))
   mm = np.min(masked_lt,axis=1)
   inner_ds.append(mm.filled(np.nan).mean())
  dAA,dBB = inner_ds
  d = dAB - (dAA + dBB) / 2.0
  return d
 #Dlengths_to = get_sp(nodes_from,nodes_to,network,network_size)
 if method == 'closest' or method == 'closest-min':
  return lengths_to.min(1).mean()
 elif method == 'shortest' or method == 'shortest-min':
  return lengths_to.mean()
 elif method == 'shortest2':
  return np.sqrt((lengths_to*lengths_to)).mean()
 elif method == 'kernel' or method == 'kernel-min':
  #val = -np.log(np.mean([np.exp(-value-1) for value in values]))
  return -np.log(np.exp(-lengths_to-1).mean(1)).mean()
 elif method == 'kernel2':
  return np.exp(lengths_to.mean(1)).mean()
 elif method == 'center':
  return lengths_to.mean()
#### END DistanceCalculations ####

