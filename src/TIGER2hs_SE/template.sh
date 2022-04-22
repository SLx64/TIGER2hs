#!/bin/bash -l
# AUTO RESTARTING TEMPLATE
#################################################################
jobname="<INPUT>"					    #jobname
numreplicas="8"							#number of replicas
mintemp="300"							#lower bound for temperature
maxtemp="450"							#upper bound for temperature
tigerheat="0"							#tiger2 heat steps
tigersample="5000"						#tiger2 sample steps
tigerquench="2500"						#tiger2 quenching steps
tigersolute="<INPUT>"					#resid selection for solute
tigersolvent="<INPUT>"					#resid selection for solvent
tigerignore=""							#resid selection to ignore part of the solute during shellsearch
tigeromm="remd/omm_impl_spe.py"  #OpenMM implicit solvent energy evaluation
tigerimplplatform="Reference"           #OpenMM platform used for energy calculation [CPU, CUDA, OpenCL]
tigerimplgb="OBC2"                      #Implicit solvent model [HCT, OBC1, OBC2, GBn, GBn2] 
tigerimplsaltcon="0"                    #implicit solvent salt concentration
tigerimplpbc="1"                        #activate pbc for implicit solvent calculations
tigerimpltop="${jobname}_gb.top"        #topology for the implicit solvent system
tigershell="<INPUT>"                    #number of shell solvents
tigerspace="3"                          #Spacing for domain decomposition in shell search
tigerconheat="1"                        #ConHeat mod
ommpre=""                               #command prefix for OpenMM call
ommsuff=""                              #command suffix for OpenMM call
minruns="1"							    #perform a minimization of this number of runs before REMD
numruns="500000"						#number of runs
runsperrestart="10"						#number of runs between writing restart files
CellX="<INPUT>"							#Cell dimensions
CellY="<INPUT>"
CellZ="<INPUT>"
timestep="4"
stepspercycle="20"
langevinpiston="1"                      #enable NPT
remdpressuregen="1"                     #generate temperature dependent target pressures for langevinpiston
langevinpistontarget="1.01325"
langevinpistondecay="100"
langevinpistonperiod="200"
######################################################################

######################### BINARY ENVIRONMENT #########################
# module load openmpi-4.0.4
# module load namd-2.14-cuda
# conda activate t2hs_openmm
RUNOPT="-n 8" #mpirun flags
NAMDOPT=""    #example namd2 flags
######################################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
######################################################################
######## DO NOT CHANGE BELOW UNLESS YOU KNOW WHAT U R DOING ##########
######################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 
############################# Prepearation ##########################
remd/make_output_dirs.sh output $numreplicas
rm output/*.mdout; #remove old OpenMM files
######################################################################

############################# BASE CONFIG ############################
cat > ${jobname}_base.namd <<END
# input options

amber on
parmfile ${jobname}_HMassR.top

# pbc 
cellBasisVector1     ${CellX} 0      0
cellBasisVector2     0      ${CellY} 0
cellBasisVector3     0      0      ${CellZ}
cellOrigin           0      0      0 

# simulation options

# Basic dynamics
exclude                 scaled1-4
1-4scaling              0.8333
dielectric              1

# Simulation space partitioning
switching 	   on
switchdist		9
cutoff         10
pairlistdist   12
rigidbonds	  all

# Multiple timestepping
timestep                ${timestep}
stepspercycle           ${stepspercycle}
fullElectFrequency      1

# PME
pme yes
pmegridspacing 1

# temperature control
langevin on
langevinDamping 1.0

# Constant Pressure Control (variable volume)
if {${langevinpiston}} {
  LangevinPiston on
  if {!${remdpressuregen}} {
    LangevinPistonTarget ${langevinpistontarget}
  }
  LangevinPistonPeriod ${langevinpistonperiod}
  LangevinPistonDecay  ${langevinpistondecay}
  useGroupPressure yes
}

wrapWater on
margin 5
END
######################################################################

######################## PARALLEL TEMPERING CONF ####################
cat > ${jobname}_remd.conf <<END
set num_replicas ${numreplicas}
set min_temp ${mintemp}
set max_temp ${maxtemp}
set tigerheat ${tigerheat}
set tigersample ${tigersample}
set tigerquench ${tigerquench}
set tigersolute ${tigersolute}
set tigersolvent ${tigersolvent}
set tigerignore "${tigerignore}"
set tigerconheat ${tigerconheat}
set tigershell ${tigershell}
set tigerspace ${tigerspace}
set tigeromm ${tigeromm}
set tigerimplgb ${tigerimplgb}
set tigerimplsaltcon ${tigerimplsaltcon}
set tigerimplpbc ${tigerimplpbc}
set tigerimpltop ${tigerimpltop}
set tigerimplplatform ${tigerimplplatform}
set ommpre "${ommpre}"
set ommsuff "${ommsuff}"
set num_runs ${numruns}
set runs_per_restart ${runsperrestart}
set namd_config_file "${jobname}_base.namd"
set output_root "output/%s/${jobname}"
set remdpressuregen ${remdpressuregen}
set minruns ${minruns}
set initial_pdb    "${jobname}.pdb"
END
######################################################################

############################# JOBS ###################################
jobnum=0
while ls output/${jobname}.job${jobnum}.restart*.tcl > /dev/null 2>&1
do
  jobnum=$(($jobnum+1))
  restartfile=$(ls -tr output/${jobname}.job$(($jobnum-1)).restart*.tcl | tail -n 1)
done

cat > job${jobnum}.conf <<END
source ${jobname}_remd.conf

if {$jobnum > 0} {
    source [format ${restartfile} ""]
}

if { ! [catch numPes] } { source remd/replica.namd }
END
mpirun $RUNOPT namd2 $NAMDOPT +idlepoll +replicas $numreplicas job${jobnum}.conf +stdout output/%d/job${jobnum}.%d.log >> output/job${jobnum}.log 2> output/job${jobnum}.out
######################################################################