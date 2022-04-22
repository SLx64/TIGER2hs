#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import logging as log

from openmm.app import *
from openmm.vec3 import Vec3
from openmm import *
import simtk.unit as u

import MDAnalysis

import warnings

warnings.filterwarnings("ignore")
AMBER_GB = {"HCT": HCT, "OBC1": OBC1, "OBC2": OBC2,
            "GBn": GBn, "GBn2": GBn2}

LOGLEVEL = log.INFO
log.basicConfig(level=LOGLEVEL,
                format='OpenMM: %(levelname)-5s %(message)s')


def timer(f):
    def timed_function(*args, **kwargs):
        from timeit import default_timer as timer

        start = timer()
        result = f(*args, **kwargs)
        end = timer()

        log.debug("Total time spend: {:.3f} seconds.".format(end - start))
        return result

    if LOGLEVEL <= log.DEBUG:
        return timed_function
    else:
        return f


def parse_arguments():
    parser = argparse.ArgumentParser(description="Hybrid Solvent Potential Energy")
    
    parser.add_argument("-t", "--topology",
                        action="store", dest="topology")
    
    parser.add_argument("-pf", "--parameter", nargs="+",
                        action="store", dest="parameters")
    
    parser.add_argument("-c", "--coordinates", 
                        action="store", dest="coordinates")
    
    parser.add_argument("-o", "--output", action="store", 
                        dest="output", default="potE.mdout")
    
    parser.add_argument("-b", "--box", action="store", 
                        dest="box", default=None)
    
    parser.add_argument("-p", "--platform", action="store", 
                        dest="platform", default="CPU")
    
    parser.add_argument("-gb", "--gbmodel", action="store", 
                        dest="gb", default="OBC2")
    
    parser.add_argument("-sc", "--saltcon", action="store",
                        type=float, dest="saltcon", default=0)
    
    parser.add_argument("--nogbsa", action="store_true",
                        dest="nogbsa")
    
    return parser.parse_args()


def read_xsc(f):
    box = np.genfromtxt(f, skip_header=2, delimiter=" ")
    if len(box) < 13:
        return [None] * 3
    else:
        return box[[1, 5, 9]]


def cell_vectors(a, b, c):
    basis_vec1 = Vec3(a, 0, 0)
    basis_vec2 = Vec3(0, b, 0)
    basis_vec3 = Vec3(0, 0, c)
    return Vec3(basis_vec1, basis_vec2, basis_vec3) * u.angstrom


@timer
def main(argv):
    if argv.gb not in AMBER_GB:
        log.warning("GB model not found... using default.")

    GBMODEL = AMBER_GB.get(argv.gb, OBC2)
    GBSA = None if argv.nogbsa else "ACE"

    topology = CharmmPsfFile(argv.topology)
    parameters = CharmmParameterSet(*argv.parameters)
    
    universe = MDAnalysis.Universe(argv.topology, argv.coordinates, format="pdb")
    
    platform = Platform.getPlatformByName(argv.platform)

    if argv.box:
        a, b, c = read_xsc(argv.box)
        universe.atoms.wrap(compound="fragments", box=(a, b, c, 90, 90, 90))
        topology.setBox(a, b, c)
        
        if a and b and c:
            cutoff = np.min([a, b, c]) / 2 - 0.0001
            cutoff *= u.angstrom
            
            system = topology.createSystem(parameters,
                                           nonbondedCutoff=cutoff,
                                           nonbondedMethod=CutoffPeriodic,
                                           implicitSolvent=GBMODEL,
                                           implicitSolventSaltConc=argv.saltcon * (u.moles/u.liter),
                                           gbsaModel=GBSA)
        else:
            log.error("Failed to create periodic system!")
            raise
    else:
        log.warning("Using a non-periodic system...")
        
        system = topology.createSystem(parameters,
                                       nonbondedMethod=NoCutoff,
                                       implicitSolvent=GBMODEL,
                                       implicitSolventSaltConc=argv.saltcon * (u.moles/u.liter),
                                       gbsaModel=GBSA)

    integrator = CustomIntegrator(0) # placeholder since we don't run a simulation
    coordinates = universe.atoms.positions * u.angstrom

    simulation = Simulation(topology.topology, system, integrator, platform)
    simulation.context.setPositions(coordinates)
    
    log.info("System setup complete.")

    state = simulation.context.getState(getEnergy=True, enforcePeriodicBox=True)
    potE = state.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
    
    log.info("Potential energy: {:.4f} kcal/mol".format(potE))
    
    if argv.output:
        with open(argv.output, "w") as output:
            output.write("{}\n".format(potE))


if __name__ == '__main__':
    argv = parse_arguments()
    main(argv)


