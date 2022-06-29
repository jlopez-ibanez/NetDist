# NetDist
Calculate distances among sets of nodes in a network. 
Methods used here are based in [Guney *et al.*](https://doi.org/10.1038/ncomms10331) and [Menche *et al.*](https://doi.org/10.1126/science.1257601) papers. 
Also selected functions (e.g *degree-aware* randomization of nodes) have been adapted to python3 from [toolbox](https://github.com/emreg00/toolbox) and [proximity](https://github.com/emreg00/proximity) repositories.

## Requirements
This code has been tested using:
 - python 3.8
   - networkx (v2.6.2)
   - numpy (v1.21)

## Usage
***NetDist*** consists of complementary scripts that can be used independently. Run them with argument *--help* to get detailed information about available options.
 
### netdist.py
This is the main script. It takes a network file in **SIF** format (check files in ***example*** folder) and two lists of nodes. Nodes should exist in the network. Default is to extract the largest connected component (LCC) from the network and save it in the same folder as the original file with the suffix *.lcc*.

Use option *--method* to indicate one of the available distance methods (closest, separation, shortest, kernel, center).

To calculate a *p-value* of the resulting measure this generates a number of randomized sets having the same degree as the nodes in the original sets and calculate the selected distance to each of them. By default the number of randomized sets to generate is **one thousand**.

For large networks and/or large sets of nodes, it is recommended to pre-calculate the matrix of all shortest paths among the nodes (see **calculate_sp.py**)

An example for calculating the **closest** measure among nodes from **set_A** and **set_B** in **network_small.sif**:

		python3 netdist.py example/network_small.sif example/set_A example/set_B -m closest
		
### calculate_sp.py

Having pre-calculated the shortest paths among all nodes in the network may speed up the calculations, specially if you want to compute different distance methods within the same network.

		python3 calculate_sp.py example/network_small.sif

## References

- Menche J, Sharma A, Kitsak M, Ghiassian SD, Vidal M, Loscalzo J, Barabási AL. **Disease networks. Uncovering disease-disease relationships through the incomplete interactome**. *Science*. 2015 Feb 20;347(6224):1257601. DOI: [10.1126/science.1257601](https://doi.org/10.1126/science.1257601) 
- Guney E, Menche J, Vidal M, Barábasi AL. **Network-based in silico drug efficacy screening**. *Nat Commun*. 2016 Feb 1;7:10331. DOI: [10.1038/ncomms10331](https://doi.org/10.1038/ncomms10331)
- [toolbox](https://github.com/emreg00/toolbox) repository
- [proximity](https://github.com/emreg00/proximity) repository
