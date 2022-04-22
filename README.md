# TIGER2hs
TIGER2 hybrid solvent replica exchange implementations for NAMD

## Usage

Setup a NAMD simulation with TIGER2hs<sup>PE</sup> for one of the example system (AAQAA, HP7, TrpCage):

    #copy files
    cp -r <TIGER2 Repository Path>/systems/<system> simulation & cd simulation
    cp -r <TIGER2 Repository Path>/src/TIGER2hs_PE remd
    cp -r <TIGER2 Repository Path>/tools .
    cp remd/template.sh simulation.sh
    
Edit the `simulation.sh` and add values according to system's README. Then run the simulation:
    
    ./simulation.sh
    
You can monitor the sampling and exchange process by `get_samping.tcl` and `exmat.tcl`:

    ./tools/exmat.tcl
    ./remd/get_sampling.tcl <timestep>
    
Further information see: [tutorial/tutorial.md](tutorial/tutorial.md)

## Required software

We normally run our simulations with an MPI built of NAMD-2.13 or later. Our implementations work fine with CUDA enabled NAMD builts. Launch options in `template.sh` may need updates for `charmrun` and specialized compute clusters. 
The hybrid solvent energy is now generated by OpenMM and MDAnalysis and hence requires a python 3.7 environment with these packages installed.

## TODO

- [ ] add parameter to choose between wrapped or unwrapped output 
- [ ] enable option write the trajectory for all replicas

## Cite this work

We are currently preparing the manuscript for the new updates.

## See also

* [Replica-Based Protein Structure Sampling Methods II: Advanced Hybrid Solvent TIGER2hs](https://doi.org/10.1021/acs.jpcb.9b03134)
* [Replica-Based Protein Structure Sampling Methods: Compromising between Explicit and Implicit Solvents](https://doi.org/10.1021/acs.jpcb.8b05178)
* [TIGER2 with solvent energy averaging (TIGER2A): An accelerated sampling method for large molecular systems with explicit representation of solvent](https://doi.org/10.1063/1.4932341)
* [TIGER2: An improved algorithm for temperature intervals with global exchange of replicas](https://doi.org/10.1063/1.3129342)
* [An improved replica-exchange sampling method: Temperature intervals with global energy reassignment](https://doi.org/10.1063/1.2780152)

## Applications

* [A Hypothesized Mechanism for Chronic Pancreatitis Caused by the N34S Mutation of Serine Protease Inhibitor Kazal-Type 1 Based on Conformational Studies](https://doi.org/10.2147/JIR.S304787)
* [Nuclear import of BCL11B is mediated by a classical nuclear localization signal and not the Krüppel-like zinc fingers](https://doi.org/10.1242/jcs.258655)
* [Birefringent Silk Fibroin Hydrogel Constructed via Binary Solvent-Exchange-Induced Self-Assembly](https://doi.org/10.1021/acs.biomac.1c00065)
* [Unveiling the N-Terminal Homodimerization of BCL11B by Hybrid Solvent Replica-Exchange Simulations](https://doi.org/10.3390/ijms22073650)
* [Phosphorylation of fibronectin influences the structural stability of the predicted interchain domain](https://doi.org/10.1021/acs.jcim.9b00555)
* [Effect of the silica nanoparticle size on the osteoinduction of biomineralized silk-silica nanocomposites](https://doi.org/10.1016/j.actbio.2020.10.043)
