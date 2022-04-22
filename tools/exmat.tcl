#!/usr/bin/tclsh
# INFO: Show global exchange matrix, replica sampling across temperatures and K/ex monitoring (Norman Geist 2021)

set B [exec tput bold]
set N [exec tput sgr0]
proc tcl::mathfunc::round2 {value digits} {expr {round(10**$digits*$value)/10.0**$digits}};

#setup vars
set i_job       0
set runs        0
array set t2r   {}
array set r2t   {}
array set temps {}
array set history_files {}

#read histories
while {[file exists job${i_job}.conf]} {
    set go 1
    
    #probe restart job
    catch {source job${i_job}.conf}
    set outputbase "$output_root.job$i_job"
    set history_filename_format "${outputbase}.%d.history"
    
    puts "Info) Browsing history for job $i_job"
    #open file pointer for each replica
    for {set r 0} {$r < $num_replicas} {incr r} {
        set filename [format $history_filename_format $r $r]
        if {[file exists $filename]} {
            set history_files($r) [open $filename "r"]
        } else {
            puts "Error) Could not find history file $filename!"
            set go 0
            break
        }
    }
    #read runs line by line for each replica
    while {$go} {
        for {set r 0} {$r < $num_replicas} {incr r} {
            gets $history_files($r) vals
            
            if {[llength $vals] >= 5} {
                set t2r([concat $runs " " [lindex $vals 1]]) $r
                set r2t([concat $runs " " $r]) [lindex $vals 1]
                set temps([lindex $vals 1]) [lindex $vals 2]
            } else {
                puts "Warn) Incomplete run $runs in job $i_job!"
                set go 0
                break
            }
        }
        if {$go} { incr runs }
    }
    #close file pointers
    foreach fp [array names history_files] {
        close $history_files($fp)
    }
    incr i_job
}
incr runs -1

# Output only makes sense after at least two runs
if {$runs < 2} {
    puts "Info) Need at least two runs, exiting!"
    exit
}
puts "Info) Read $runs runs! ($num_replicas replicas)"

#run through data and construct exmat, locmat and diffs
array set exmat {}
array set diffs {}
array set locmat {}
for {set i 1} {$i < $runs} {incr i} {

    #look foreach temperature if and were it has moved
    for {set t 0} {$t < $num_replicas} {incr t} {
        set rold $t2r([concat [expr $i-1] " " $t])
        set rnew $t2r([concat       $i    " " $t])
        #-
        set told $r2t([concat [expr $i-1] " " $rold])
        set tnew $r2t([concat [expr $i-1] " " $rnew])
        
        #register temperature distribution across physical replicas
        incr locmat([concat $rold $told])
        
        if {$rold != $rnew} {
            #exchange detected
            #puts "EXCHANGE ACCEPT RUN [expr $i-1] $told $tnew"
            incr exmat([concat $told $tnew])
            
            #Remember T change
            if {![info exists diffs($rold)]} { set diffs($rold) 0 }
            set diffs($rold) [expr $diffs($rold)+abs($temps($told)-$temps($tnew))]
        }
    }
}

#output exchange matrix
set avg 0.
puts ""
puts "-= GLOBAL EXCHANGE MATRIX =-"
for {set r1 0} {$r1 < $num_replicas} {incr r1} {

    #Column label
    if {$r1 == 0} {
        for {set r 0} {$r < $num_replicas} {incr r} { puts -nonewline "\t$r" }
    }
    
    puts ""
    for {set r2 0} {$r2 < $num_replicas} {incr r2} {
        #Row label
        if {$r2 == 0} {
            puts -nonewline "$r1\t"
        }
        
        #diagonal
        if {$r1 == $r2} {
            set csum 0
            for {set r3 0} {$r3 < $num_replicas} {incr r3} {
                if { [info exists exmat([concat $r1 $r3])] } { set temp $exmat([concat $r1 $r3]) } else { set temp 0}
                set csum [expr $csum + $temp]
            }
            set csum [expr 1.*($runs-$csum)/$runs]
            set avg [expr $avg + (1.-$csum)]
            set csum [expr round2($csum,2)]
            puts -nonewline "$B$csum$N\t"
        } else {
            if { [info exists exmat([concat $r1 $r2])] } { set temp $exmat([concat $r1 $r2]) } else { set temp 0}
            set temp [expr round2(1.*$temp/$runs,2)]
            puts -nonewline "$temp\t"
        }
    }
}
puts ""
set avg [expr $avg/$num_replicas]
puts "AVG EX. RATE ${B}$avg$N"

#output location matrix
puts ""
puts "-= DISTRIBUTION MATRIX =- (T=Temperature;R=Rank)"
for {set r1 0} {$r1 < $num_replicas} {incr r1} {

    #Column label
    if {$r1 == 0} {
        for {set r 0} {$r < $num_replicas} {incr r} { puts -nonewline "\tR$r" }
    }
    
    puts ""
    for {set r2 0} {$r2 < $num_replicas} {incr r2} {
        #Row label
        if {$r2 == 0} {
            puts -nonewline "T$r1\t"
        }
        
        if { [info exists locmat([concat $r1 $r2])] } { set temp $locmat([concat $r1 $r2]) } else { set temp 0}
        set temp [expr round2(1.*$temp/($runs-1),2)]
        puts -nonewline "$temp\t"
    }
}
puts ""

#output K/ex
puts ""
puts "-= RANDOM WALK SPEED =- (dT/X)"
set diffavg 0
for {set r 0} {$r < $num_replicas} {incr r} {
    if { [info exists diffs($r)] } { set temp [expr $diffs($r) / ($runs-1)] } else { set temp 0 }
    set diffavg [expr $diffavg+$temp]
    puts "R$r\t$temp"
}
puts "-----------------------"
set diffavg [expr $diffavg / $num_replicas]
puts "AVG\t$diffavg"

puts ""
puts "Exiting!"