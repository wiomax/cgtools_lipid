# ::cgtools::analysis::analyze_stress_tensor --
#
#  Calculate the pressure tensor of the system.  
#  by Zun-Jing Wang Jan. 2009

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::stress_tensor {
    variable all_particles
    variable f_tvsstress_tensor
    variable verbose
    variable av_stress_tensor { 0 0 0 0 0 0 0 0 0 }
    variable av_stress_tensor_i 0
    namespace export setup_stress_tensor
    namespace export analyze_stress_tensor
    namespace export printav_stress_tensor
    namespace export resetav_stress_tensor
}

proc ::cgtools::analysis::stress_tensor::resetav_stress_tensor { } {
    variable av_stress_tensor
    variable av_stress_tensor_i
    set av_stress_tensor { 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 }
    set av_stress_tensor_i 0
}

proc ::cgtools::analysis::stress_tensor::printav_stress_tensor { } {
    global ::cgtools::analysis::time
    variable av_stress_tensor
    variable f_tvsstress_tensor
    variable av_stress_tensor_i
    variable verbose

    if { $av_stress_tensor_i > 0 } {
	puts -nonewline $f_tvsstress_tensor "$time "
	for { set v 0 } { $v < [llength $av_stress_tensor] } {incr v} {
	    puts -nonewline $f_tvsstress_tensor "[expr [lindex $av_stress_tensor $v]/($av_stress_tensor_i*1.0)] "
	    
	}
	puts $f_tvsstress_tensor " "
    } else {
	mmsg::warn [namespace current] "can't print average stress_tensor"
	flush stdout		       
    }
    flush $f_tvsstress_tensor
}

proc ::cgtools::analysis::stress_tensor::setup_stress_tensor { args } {
    # Use the global command to access variables in the std_analysis namespace
    global ::cgtools::analysis::n_particles
    global ::cgtools::analysis::suffix
    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::iotype

    variable verbose
    variable all_particles
    variable f_tvsstress_tensor

    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_stress_tensor verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    for { set j 0 } { $j < $n_particles } { incr j } {
	lappend all_particles $j
    }
    mmsg::debug [namespace current] "opening $outputdir/time_vs_stress_tensor$suffix "
    
    if { [file exists "$outputdir/time_vs_stress_tensor$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvsstress_tensor [open "$outputdir/time_vs_stress_tensor$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsstress_tensor "\# Components of the total pressure tensor in row major order"
	puts $f_tvsstress_tensor "\# pxx pxy pxz pyx pyy pyz pzx pzy pzz"
    }
}

proc ::cgtools::analysis::stress_tensor::analyze_stress_tensor {  } {
    variable all_particles
    variable av_stress_tensor
    variable av_stress_tensor_i
    variable verbose
    mmsg::send [namespace current] "analyzing stress_tensor"
    
    set blen [setmd box_l]
    set tot_volume [expr [lindex $blen 0]*[lindex $blen 1]*[lindex $blen 2]]
    
    set ik1_tmp [analyze stress_tensor ]
#    puts "k1_tmp is $ik1_tmp"
    
    
    
    set tp 0
    for { set i 0 } { $i < 9 } { incr i } {
	set component  "[lindex [lindex $ik1_tmp $tp] [expr $i+1] ]"
	lset av_stress_tensor  $i [expr [lindex $av_stress_tensor $i] + $component]
    }
    
    incr av_stress_tensor_i
    
    if { $verbose } {
	set avpx [expr [lindex $av_stress_tensor 0]/($av_stress_tensor_i*1.0)]
	set avpy [expr [lindex $av_stress_tensor 4]/($av_stress_tensor_i*1.0)]
	set avpz [expr [lindex $av_stress_tensor 8]/($av_stress_tensor_i*1.0)]
	mmsg::send [namespace current] "<p> : xx:$avpx yy:$avpy zz:$avpz"
	
	set matrixout [format "\n%.5f %.5f %.5f \n%.5f %.5f %.5f \n%.5f %.5f %.5f \n" [lindex $av_stress_tensor 0] [lindex $av_stress_tensor 1] [lindex $av_stress_tensor 2] [lindex $av_stress_tensor 3] [lindex $av_stress_tensor 4] [lindex $av_stress_tensor 5]  [lindex $av_stress_tensor 6] [lindex $av_stress_tensor 7] [lindex $av_stress_tensor 8]]
	mmsg::send [namespace current] $matrixout
    }
    mmsg::debug [namespace current] "done"
    
}

