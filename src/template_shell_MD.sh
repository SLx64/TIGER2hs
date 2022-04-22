#!/bin/bash -l
# AUTO RESTARTING MD TEMPLATE

#################################################################
jobname="<INPUT>"
dcdfreq="5000"
restartfreq="50000"
CellX="<INPUT>"
CellY="<INPUT>"
CellZ="<INPUT>"
timestep="4"
stepspercycle="20"
langevindamping="1.0"
langevinpistontarget="1.01325"
langevinpistondecay="100"
langevinpistonperiod="200"
#--------------------------------------------------------------------
ensamblenames="min     h300     prod"
ensembletemps="0        300      300"
ensamblesteps="50000 250000   500000"
ensambletypes="MIN      NVT      NPT"
ensamblerestarts="none    o        o"; #o for finished pre-stage; r for unfinished pre-stage; a for auto repitition
#################################################################

################### CUSTOM ENSAMBLE TYPE INCLUDES ###############
cat > "h300.namd" <<END
    #stage based stuff
END

cat > "NPT.namd" <<END
    #ensemble based stuff
    margin 5
END
################################################################

######################### BINARY ENVIRONMENT ###################
# module load openmpi-4.0.4
# module load namd-2.14-cuda
RUNOPT="" #mpirun flags
NAMDOPT=""    #example namd2 flags
################################################################
# 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
################################################################
##### DO NOT CHANGE BELOW UNLESS YOU KNOW WHAT U R DOING #######
################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

############################# Prepearation hacks ###############
Aensamblenames=($ensamblenames)
Aensembletemps=($ensembletemps)
Aensamblesteps=($ensamblesteps)
Aensambletypes=($ensambletypes)
Aensamblerestarts=($ensamblerestarts)
################################################################

######################### SCRIPT GEN & EXEC ####################
enameold="false"
estepsold="false"
etypeold="false"
erestartold="false"
for ((i=0; i<${#Aensamblenames[*]}; i++))
do
  ename=${Aensamblenames[$i]}
  etemp=${Aensembletemps[$i]}
  esteps=${Aensamblesteps[$i]}
  etype=${Aensambletypes[$i]}
  erestart=${Aensamblerestarts[$i]}
  nexterestart=${Aensamblerestarts[$i+1]}
  
  #automatic restarting
  if [[ $erestart == "a" ]]; then
    round=1
    while [ -a "${jobname}_${ename}_r$round.restart.coor" ]
    do
        round=$((round+1))
    done

    if [ $round -gt 1 ]; then
        enameold="${ename}_r$((round-1))"
    fi
    
    ename="${ename}_r$round"
  fi

  #not yet done?
  if [[ ! (-a "${jobname}_$ename.coor") && ! (($nexterestart == "r" || $nexterestart == "a") && -a "${jobname}_$ename.restart.coor") ]]; then

cat > ${jobname}_${ename}.namd <<END
####### basic input #######
amber on
parmfile ${jobname}_HMassR.top
coordinates ${jobname}.pdb
###########################

####### restart input #######
if {"${erestart}" == "o"} {
  bincoordinates ${jobname}_${enameold}.coor
  if {"${etypeold}" != "MIN"} {
    binvelocities ${jobname}_${enameold}.vel
  }
  extendedsystem ${jobname}_${enameold}.xsc
} elseif {"${erestart}" == "r" || "${erestart}" == "a"} {
  bincoordinates ${jobname}_${enameold}.restart.coor
  if {"${etypeold}" != "MIN"} {
    binvelocities ${jobname}_${enameold}.restart.vel
  }
  extendedsystem ${jobname}_${enameold}.restart.xsc
} else {
  cellBasisVector1     ${CellX} 0        0
  cellBasisVector2     0        ${CellY} 0
  cellBasisVector3     0        0        ${CellZ}
  cellOrigin           0        0        0 
}
############################

########### output options ################
outputname              ${jobname}_${ename}
binaryoutput            yes
binaryrestart           yes
outputEnergies          1000
outputTiming            1000
restartfreq             ${restartfreq}
dcdfreq                 ${dcdfreq}
###########################################

########### simulation options ############
# Basic dynamics
exclude                 scaled1-4
1-4scaling              0.8333
dielectric              1

# Simulation space partitioning
switching      on
switchdist      9
cutoff         10
pairlistdist   12
rigidbonds    all

# Multiple timestepping
timestep                ${timestep}
stepspercycle           ${stepspercycle}
fullElectFrequency      1

# PME
pme yes
pmegridspacing 1
###########################################

############# Thermostat ##################
if {"${etype}" == "NVT" || "${etype}" == "NPT"} {
  langevin                on
  langevinTemp            $etemp
  langevinHydrogen        on
  langevinDamping         ${langevindamping}
}
if {"${etypeold}" == "MIN" || "${erestart}" == "none"} {
  # initial temperature
  temperature $etemp
}
###########################################

######## Constant Pressure Control ########
if {"${etype}" == "NPT"} {
  LangevinPiston                     on
  LangevinPistonTarget               ${langevinpistontarget}
  LangevinPistonPeriod               ${langevinpistonperiod}
  LangevinPistonDecay                ${langevinpistondecay}
  LangevinPistonTemp                 $etemp
}
###########################################

########## custom includes ###############
if {[file exists ${etype}.namd]} {
  source ${etype}.namd
}
if {[file exists ${ename}.namd]} {
  source ${ename}.namd
}
###########################################

################# RUN #####################
if {"${etype}" == "MIN"} {
  minimize ${esteps}
} else {
  run ${esteps}
}
###########################################
END
    mpirun $RUNOPT namd2 $NAMDOPT  +idlepoll ${jobname}_${ename}.namd > ${jobname}_${ename}.out 2> ${jobname}_${ename}.e
fi
  enameold=$ename
  estepsold=$steps
  etypeold=$etype
  erestartold=$erestart
done
#################################################################
