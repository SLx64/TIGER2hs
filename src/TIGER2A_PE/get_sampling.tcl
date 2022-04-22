#!/usr/bin/tclsh
if {$argc < 1} {
  puts "args: <timestep_fs>"
  exit;
}

set tcl_precision 4

#setup vars
set timestep [lindex $argv 0]
set i_job       0
set runs_total  0
set time_total  0
set time_quench 0
set time_avg    0
set samples     0

#read history
while {[file exists job${i_job}.conf]} {
  
  catch {source job${i_job}.conf}

  set runs 0
  set outputbase "$output_root.job$i_job"
  set history_filename_format "${outputbase}.%d.history"

  set dt [expr ((1.+$tigerheat+$tigersample) * $timestep) / 1000]
  puts -nonewline "Info) Browsing history for job $i_job (Δt = $dt ps w/o quenching & averaging)"
  for {set replica_id 0} {$replica_id < $num_replicas} {incr replica_id} {
    set history_file [open [format $history_filename_format $replica_id $replica_id] "r"]
    while {[gets $history_file vals] >= 0} {
      incr runs
    }
    close $history_file
  }
  puts ""
  
  #add current sampling time with current steps_per_run (in case it changes between runs)
  set time_total  [expr 1.* $time_total   + ($runs * ($tigerheat+$tigersample+$tigerquench+$tigeravg) * $timestep)]
  set time_quench [expr 1.* $time_quench  + ($runs *                          $tigerquench            * $timestep)]
  set time_avg    [expr 1.* $time_avg     + ($runs *                                       $tigeravg  * $timestep)]
  set samples     [expr 1. *  $samples    +  $runs]
  
  incr i_job
  incr runs_total [expr $runs / $num_replicas]
}
puts "Info) Read $runs_total runs! ($num_replicas replicas)"

set time_total         [expr $time_total  / 1000000.]
set time_quench        [expr $time_quench / 1000000.]
set time_avg           [expr $time_avg    / 1000000.]
set time_diff          [expr $time_total - $time_quench - $time_avg]
set time_base          [expr $time_total  / $num_replicas]
set time_base_quench   [expr $time_quench / $num_replicas]
set time_base_avg      [expr $time_avg / $num_replicas]
set time_base_diff     [expr $time_base - $time_base_quench - $time_base_avg]
set samples_base       [expr $samples    / $num_replicas]
set base_samples_ratio [expr 1 / ((1. * $samples_base) / ($time_total * 1000))]

puts "Info) Total sampling 	$time_diff ns	+ $time_quench ns quenching 	+ $time_avg ns averaging 	= $time_total ns"
puts "Info) Total base 	$time_base_diff ns	+ $time_base_quench ns quenching 	+ $time_base_avg ns averaging 	= $time_base ns ($samples_base structures)"
puts "Info) Ratio $base_samples_ratio ps*replicas/sample"
