#!/usr/bin/tclsh
puts "Usage: inpdb outpdb solute solvent ignoresel dspace xsc maxsol"

proc timer_read {} {
  set clocksource "/proc/uptime"
  if {[file exists $clocksource]} {
    set fp [open $clocksource r]
    set raw [read $fp]
    close $fp
    
    return [lindex $raw 0]
  } else {
    puts "Warn: TIMING INFORMATION ONLY AVAILABLE TO UNIX BASED SYSTEMS ($clocksource)"
    return 0
  }
}

proc timer_init {} {
  global TIMERSTAMPS
  array set TIMERSTAMPS {}
}

proc timer_start {name} {
  global TIMERSTAMPS
  set TIMERSTAMPS($name) [timer_read]
}

proc timer_now {name} {
  global TIMERSTAMPS
  if {[info exists TIMERSTAMPS($name)]} {
    set lstamp [expr [timer_read]-$TIMERSTAMPS($name)]
    set TIMERSTAMPS($name.val) $lstamp
    return $lstamp
  } else {
    puts "Warn: TIMER $name NOT RUNNING"
    return 0
  }
}

proc timer_stop {name} {
  global TIMERSTAMPS
  if {[info exists TIMERSTAMPS($name)]} {
    set lstamp [expr [timer_read]-$TIMERSTAMPS($name)]
    set TIMERSTAMPS($name.val) $lstamp
    set TIMERSTAMPS($name) 0
    return $lstamp
  } else {
    puts "Warn: TIMER $name NOT STARTED"
    return 1
  }
}

proc timer_get {name} {
  global TIMERSTAMPS
  if {[info exists TIMERSTAMPS($name.val)]} {
    return $TIMERSTAMPS($name.val)
  } else {
    puts "Warn: TIMER $name NOT STOPPED"
    return 1
  }
}

proc s2human {secs {prec 2}} {
  set dur ""
  set factors {1. 60. 60. 24.}
  set captions {"s" "m" "h" "d"}
  set values 0
  lset values 0 $secs
  
  for {set i 1} {$i<[llength $factors]} {incr i} {
    lappend values [expr int([lindex $values [expr $i-1]]/[lindex $factors $i])]
    if {[lindex $values $i] > 0} {
      lset values [expr $i-1] [expr [lindex $values [expr $i-1]] - [lindex $values $i]*[lindex $factors $i]]
    } else {
      set i [llength $factors]
    }
  }
  for {set i 0} {$i<[llength $values]} {incr i} {
    if {[lindex $values $i] > 0} {
      set tprec [expr $i > 0 ? 0 : $prec]
      set dur [concat [format %.${tprec}f [lindex $values $i]][lindex $captions $i] $dur]
    }
  }
  return $dur
}

proc sel2map {sel} {
  array set map {}
  set ranges [split $sel ","]
  foreach range $ranges {
    set items [split $range "-"]
    if {[llength $items] > 1} {
      for {set i [lindex $items 0]} {$i <= [lindex $items 1]} {incr i} {
	set map($i) 1
      }
    } else {
      set map($items) 1
    }
  }
  return [array get map]
}

proc getBoundaries {atom_array natom} {
  array set atoms  $atom_array
  array set minmax {x.center 0 y.center 0 z.center 0}
  for {set i 1} {$i <= $natom} {incr i} {
    set atom $atoms($i)
    set x [lindex $atom 3]
    set y [lindex $atom 4]
    set z [lindex $atom 5]
    
    set minmax(x.center) [expr {$minmax(x.center) + $x}]
    set minmax(y.center) [expr {$minmax(y.center) + $y}]
    set minmax(z.center) [expr {$minmax(z.center) + $z}]
    
    if {![info exists minmax(x.min)] || $x < $minmax(x.min)} { set minmax(x.min) $x }
    if {![info exists minmax(y.min)] || $y < $minmax(y.min)} { set minmax(y.min) $y }
    if {![info exists minmax(z.min)] || $z < $minmax(z.min)} { set minmax(z.min) $z }
    if {![info exists minmax(x.max)] || $x > $minmax(x.max)} { set minmax(x.max) $x }
    if {![info exists minmax(y.max)] || $y > $minmax(y.max)} { set minmax(y.max) $y }
    if {![info exists minmax(z.max)] || $z > $minmax(z.max)} { set minmax(z.max) $z }
  }
  set minmax(x.center) [expr $minmax(x.center) / $natom]
  set minmax(y.center) [expr $minmax(y.center) / $natom]
  set minmax(z.center) [expr $minmax(z.center) / $natom]
  
  set minmax(x.dim) [expr $minmax(x.max) - $minmax(x.min)]
  set minmax(y.dim) [expr $minmax(y.max) - $minmax(y.min)]
  set minmax(z.dim) [expr $minmax(z.max) - $minmax(z.min)]
  
  return [array get minmax]
}

proc move {atom_array natom vec} {
  array set atoms $atom_array
  for {set i 1} {$i <= $natom} {incr i} {
    set atom $atoms($i)
    set x [lindex $atom 3]
    set y [lindex $atom 4]
    set z [lindex $atom 5]
    
    set atom [lset atom 3 [expr {$x + [lindex $vec 0]}]]
    set atom [lset atom 4 [expr {$y + [lindex $vec 1]}]]
    set atom [lset atom 5 [expr {$z + [lindex $vec 2]}]]
    set atoms($i) $atom
  }

  return [array get atoms]
}

proc tcl::mathfunc::pmodulo {n m} {
  return [expr {$n-$m * floor(1.* $n/$m)}]
}

proc wrap {atom_array natom boxvec} {
  array set atoms $atom_array
  set boxx [lindex $boxvec 0]
  set boxy [lindex $boxvec 1]
  set boxz [lindex $boxvec 2]
  
  for {set i 1} {$i <= $natom} {incr i} {
    set atom $atoms($i)
    set x [lindex $atom 3]
    set y [lindex $atom 4]
    set z [lindex $atom 5]
    
    #wrapping on-demand saves time alot
    if {$x < 0 || $x >= $boxx} { set atom [lset atom 3 [expr {pmodulo($x,$boxx)}]] }
    if {$y < 0 || $y >= $boxy} { set atom [lset atom 4 [expr {pmodulo($y,$boxy)}]] }
    if {$z < 0 || $z >= $boxz} { set atom [lset atom 5 [expr {pmodulo($z,$boxz)}]] }
    
    set atoms($i) $atom
  }
  
  return [array get atoms]
}

#perform domain decomposition with given pbc box 
#using the same box for two molecules and
#the generated domain neighbor lists, fast atomistic nearest  
#neighbors searches are easily realized
#Hint: If box is periodic, it's safe to have multiple solutes
proc decompose {atom_array natom dspace boxvec {neighborlists 0}} {
  array set atoms  $atom_array
   
  #get grid sizes
  set xdim   [lindex $boxvec 0]
  set ydim   [lindex $boxvec 1]
  set zdim   [lindex $boxvec 2]
  set xs     [expr int(ceil($xdim / $dspace))]
  set ys     [expr int(ceil($ydim / $dspace))]
  set zs     [expr int(ceil($zdim / $dspace))]
  puts "Info: TIGER2hs - DOMAIN GRID     [format %8s $xs] [format %8s $ys] [format %8s $zs]"
  
  #build domains 
  array set domains {}
  array set xocc    {}
  array set yocc    {}
  array set zocc    {}
  for {set i 1} {$i <= $natom} {incr i} {
    set atom $atoms($i)
    set x [lindex $atom 3]
    set y [lindex $atom 4]
    set z [lindex $atom 5]
    set ix [expr {int($x / $dspace)}]
    set iy [expr {int($y / $dspace)}]
    set iz [expr {int($z / $dspace)}]
    
    incr    domains($ix.$iy.$iz.natom)
    lappend domains($ix.$iy.$iz) $i
    
    set xocc($ix) 1
    set yocc($iy) 1
    set zocc($iz) 1
  }
  set xsc [llength [array names xocc]]
  set ysc [llength [array names yocc]]
  set zsc [llength [array names zocc]]
  puts "Info: TIGER2hs - DOMAIN OCCUPIED [format %8s $xsc] [format %8s $ysc] [format %8s $zsc]"
  
  #generate neighbor lists
  if {$neighborlists} {
    for {set ix 0} {$ix <= $xs} {incr ix} {
      for {set iy 0} {$iy <= $ys} {incr iy} {
	for {set iz 0} {$iz <= $zs} {incr iz} {
	  if {[info exists domains($ix.$iy.$iz)] && $domains($ix.$iy.$iz.natom) > 0} {
	    #add respective lower and higher neighbors in 
	    #all dimensions also considering periodicity
	    array set temp {}
	    set temp(ix) $ix
	    set temp(iy) $iy
	    set temp(iz) $iz
	    set temp(lix) [expr {int(pmodulo($ix-1,$xs))}]
	    set temp(hix) [expr {int(pmodulo($ix+1,$xs))}]
	    set temp(liy) [expr {int(pmodulo($iy-1,$ys))}]
	    set temp(hiy) [expr {int(pmodulo($iy+1,$ys))}]
	    set temp(liz) [expr {int(pmodulo($iz-1,$zs))}]
	    set temp(hiz) [expr {int(pmodulo($iz+1,$zs))}]
	    
	    foreach xitem {ix lix hix} {
	      foreach yitem {iy liy hiy} {
		foreach zitem {iz liz hiz} {
		  lappend domains($ix.$iy.$iz.neighbors) $temp($xitem).$temp($yitem).$temp($zitem)
		}
	      }
	    }
	    array unset temp
	  }
	}
      }
    }
  }
  
  #return data 
  set domains(x.max) $xs
  set domains(y.max) $ys
  set domains(z.max) $zs
  return [array get domains]
}

proc readxsc {xsc} {
  #read input xsc
  set fp [open "$xsc" r]
  set raw [read $fp]
  close $fp
  
  #read just the 1st non-commend line
  set lines [split $raw "\n"]
  foreach line $lines {
    if {[regexp {^[0-9]} $line]} {
      set x [lindex $line 1]
      set y [lindex $line 5]
      set z [lindex $line 9]
      return [list $x $y $z]
    }
  }
  
  #if we get here, the xsc is bad
  return 1
}

proc writepdb {atom_array natom outpdb {boxvec {0 0 0}}} {
  array set atoms $atom_array
  set fp [open "$outpdb" w]
  
  #write cryst1 record
  puts -nonewline $fp "CRYST1"
  puts -nonewline  $fp [format %9.3f [lindex $boxvec 0]]
  puts -nonewline  $fp [format %9.3f [lindex $boxvec 1]]
  puts -nonewline  $fp [format %9.3f [lindex $boxvec 2]]
  puts -nonewline  $fp [format %7.2f 90]
  puts -nonewline  $fp [format %7.2f 90]
  puts -nonewline  $fp [format %7.2f 90]
  puts $fp ""
  
  for {set i 1} {$i <= $natom} {incr i} {
    set atom $atoms($i)
    set name  [lindex $atom 0]
    set res   [lindex $atom 1]
    set resid [lindex $atom 2]
    set x     [lindex $atom 3]
    set y     [lindex $atom 4]
    set z     [lindex $atom 5]
    
    puts -nonewline $fp "ATOM  "
    puts -nonewline $fp [format %5u [expr $i % 100000]]
    puts -nonewline $fp " "
    puts -nonewline $fp [format %4s $name]
    puts -nonewline $fp " "
    puts -nonewline $fp [format %4s $res]
    puts -nonewline $fp " "
    puts -nonewline $fp [format %4u [expr $resid % 10000]]
    puts -nonewline $fp " "
    puts -nonewline $fp "   "
    puts -nonewline  $fp [format %8.3f $x]
    puts -nonewline  $fp [format %8.3f $y]
    puts -nonewline  $fp [format %8.3f $z]
    puts -nonewline  $fp [format %6.2f 0]
    puts -nonewline  $fp [format %6.2f 0]
    puts $fp ""
  }
  close $fp
}

proc filterPDB {inpdb outpdb solutesel solventsel ignoresel dspace boxvec maxsol} {
  puts "Info: TIGER2hs - PROCESSING $inpdb"
  puts "Info: TIGER2hs - SOLUTE  RESIDUES  $solutesel"
  puts "Info: TIGER2hs - SOLVENT RESIDUES  $solventsel"
  puts "Info: TIGER2hs - IGNORE FOR SHELL  $ignoresel"
  
  #read input pdb
  set fp [open "$inpdb" r]
  set raw [read $fp]
  close $fp
  
  #setup selection maps
  array set solute  [sel2map $solutesel]
  array set ignore  [sel2map $ignoresel]
  array set solvent [sel2map $solventsel]
  
  #read selected atoms from inpdb --------------------------------------------------------------------------------
  set csolute  0
  set csolute2 0
  set csolvent 0
  set cres     0
  set l  ""
  array set solvent_atoms {}
  array set solute_atoms  {}
  array set solute_atoms2 {}
  set lines [split $raw "\n"]
  foreach line $lines {
    if {[regexp {^ATOM} $line]} {
      #get columns from pdb line
      set name  [string trim [string range $line 12 15]]
      set res   [string trim [string range $line 17 20]]
      set resid [string trim [string range $line 22 26]]
      set x     [string trim [string range $line 30 37]]
      set y     [string trim [string range $line 38 45]]
      set z     [string trim [string range $line 46 53]]
      
      #renum resids in case of hex numbers etc.
      #so solutesel and solventsel are reliable
      if {$resid != $l} {
	set l $resid
	incr cres
      }
      set resid $cres

      #store atom if selected
      if {[info exists solute($resid)] && ![info exists ignore($resid)]} {
	incr csolute
	set solute_atoms($csolute) [list $name $res $resid $x $y $z]
      } 
      if {[info exists solute($resid)]} {
	incr csolute2
	set solute_atoms2($csolute2) [list $name $res $resid $x $y $z]
      }
      if {[info exists solvent($resid)]} {
      	incr csolvent
	set solvent_atoms($csolvent) [list $name $res $resid $x $y $z]
      }
    }
  }
  set tignore [expr $csolute2-$csolute]

  puts "Info: TIGER2hs - SOLUTE  ATOMS $csolute ($tignore ignored during shellsearch)"
  puts "Info: TIGER2hs - SOLVENT ATOMS $csolvent"
  #---------------------------------------------------------------------------------------------------------------
  
  #center solute and wrap coordinates ---------------------------------------------------------------------------
  array set solute_minmax [getBoundaries [array get solute_atoms] $csolute]
  set boxx   [lindex $boxvec 0]
  set boxy   [lindex $boxvec 1]
  set boxz   [lindex $boxvec 2]
  set shiftX [expr ($boxx/2) - $solute_minmax(x.center)]
  set shiftY [expr ($boxy/2) - $solute_minmax(y.center)]
  set shiftZ [expr ($boxz/2) - $solute_minmax(z.center)]
  puts "Info: TIGER2hs - CENTERING SOLUTE $shiftX $shiftY $shiftZ"
  array set solute_atoms  [move [array get solute_atoms ] $csolute  [list $shiftX $shiftY $shiftZ]]
  array set solute_atoms2 [move [array get solute_atoms2] $csolute2 [list $shiftX $shiftY $shiftZ]]
  array set solvent_atoms [move [array get solvent_atoms] $csolvent [list $shiftX $shiftY $shiftZ]]
  puts "Info: TIGER2hs - APPLYING PBC BOX $boxvec"
  array set solute_atoms  [wrap [array get solute_atoms ] $csolute  $boxvec]
  array set solute_atoms2 [wrap [array get solute_atoms2] $csolute2 $boxvec]
  array set solvent_atoms [wrap [array get solvent_atoms] $csolvent $boxvec]
  #---------------------------------------------------------------------------------------------------------------

  #perform domain decomposition ----------------------------------------------------------------------------------
  puts "Info: TIGER2hs - DOMAIN DECOMPOSITION AND NEIGHBOR LIST SETUP SPACING=$dspace..."
  array set solute_domains  [decompose [array get solute_atoms ] $csolute  $dspace $boxvec 1]
  array set solvent_domains [decompose [array get solvent_atoms] $csolvent $dspace $boxvec 0]
  #---------------------------------------------------------------------------------------------------------------
  
  #perform nearest neighbor search based on domains --------------------------------------------------------------
  puts "Info: TIGER2hs - ENTERING FAST NEAREST NEIGHBOR SEARCH..."
  set exhaust [expr $csolute * $csolvent]
  set pairs   0
  array set dists {}
  #loop solute domains
  for {set ix 0} {$ix <= $solute_domains(x.max)} {incr ix} {
    for {set iy 0} {$iy <= $solute_domains(y.max)} {incr iy} {
      for {set iz 0} {$iz <= $solute_domains(z.max)} {incr iz} {
	#domain exists and not empty?
	if {[info exists solute_domains($ix.$iy.$iz)] && $solute_domains($ix.$iy.$iz.natom) > 0} {
	  #loop possible neighbors from solvent
	  foreach neighbor $solute_domains($ix.$iy.$iz.neighbors) {
	    #domain exists and not empty?
	    if {[info exists solvent_domains($neighbor)] && $solvent_domains($neighbor.natom) > 0} {
	      #perform all-atom measures on occupied neighboring domain
	      foreach item1 $solute_domains($ix.$iy.$iz) {
		set atom1 $solute_atoms($item1)
		set x1 [lindex $atom1 3]
		set y1 [lindex $atom1 4]
		set z1 [lindex $atom1 5]
		foreach item2 $solvent_domains($neighbor) {
		  set atom2 $solvent_atoms($item2)
		  set resid [lindex $atom2 2]
		  set x2    [lindex $atom2 3]
		  set y2    [lindex $atom2 4]
		  set z2    [lindex $atom2 5]
		  #shorter periodic distance
		  set dx  [expr {pmodulo($x1-$x2+($boxx/2.),$boxx) - ($boxx/2.)}]
		  set dy  [expr {pmodulo($y1-$y2+($boxy/2.),$boxy) - ($boxy/2.)}]
		  set dz  [expr {pmodulo($z1-$z2+($boxz/2.),$boxz) - ($boxz/2.)}]
		  set dsq [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
		  lappend dists($dsq) $resid
		  incr pairs
		}
	      }
	    }
	  }
	}
      }
    }
  }
  puts "Info: TIGER2hs - COMPLETED SHELL SEARCH - DD COMPUTED ONLY $pairs PAIRS OUT OF $exhaust"
  #-----------------------------------------------------------------------------------------------------------------
  
  #identify solvent shell by distance ----------------------------------------------------------------------------
  set sorted [lsort -real -increasing [array names dists]]
  array set shell  {}
  array set scount {}
  foreach item $sorted {
    foreach resid $dists($item) {
      if {![info exists shell($resid)] && [llength [array names shell]] < $maxsol} {
	set shell($resid) 1
	puts "Info: TIGER2hs - ASSIGNED $resid TO SHELL WITH DISTANCE [expr sqrt($item)]"
      }
      set scount($resid) 1
    }
  }
  set sc [llength [array names scount]]
  puts "Info: TIGER2hs - TOTAL NUMBER OF SHELL CANDIDATES $sc"
  
  #fast neighbor search tweaks
  if {$sc < $maxsol} {
    puts stderr "Error: TIGER2hs - NEAREST NEIGHBOR SEARCH RETURNED LESS SOLVENTS THAN REQUESTED $sc < $maxsol"
    puts stderr "Error: TIGER2hs - INCREASE dpace"
    exit 1
  } elseif {[expr ($sc > ($maxsol * 1.66))]} {
    puts "Warn: TIGER2hs - NEAREST NEIGHBOR SEARCH RETURNED MUCH MORE SOLVENTS THAN REQUESTED"
    puts "Warn: TIGER2hs - DECREASING dpace MIGHT IMPROVE PERFORMANCE"
  }
  #---------------------------------------------------------------------------------------------------------------
  
  #join solute2 with selected solvents and write pdb --------------------------------------------------------------
  set ac $csolute2
  for {set i 1} {$i <= $csolvent} {incr i} {
    set resid [lindex $solvent_atoms($i) 2]
    if {[info exists shell($resid)]} {
      incr ac
      set solute_atoms2($ac) $solvent_atoms($i)
    }
  }
  
  writepdb [array get solute_atoms2] $ac $outpdb
  puts "Info: TIGER2hs - DONE CREATING SHELL PDB $outpdb WITH $maxsol SOLVENTS"
  #---------------------------------------------------------------------------------------------------------------
}

timer_init

timer_start "shellsearch"
#         inpdb            outpdb           solute           solvent          solute_ignoresel dspace           xsc                        maxsol          "
filterPDB [lindex $argv 0] [lindex $argv 1] [lindex $argv 2] [lindex $argv 3] [lindex $argv 4] [lindex $argv 5] [readxsc [lindex $argv 6]] [lindex $argv 7]
timer_stop "shellsearch"

puts [s2human [timer_get "shellsearch"]]