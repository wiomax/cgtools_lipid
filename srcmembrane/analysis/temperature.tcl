# ::cgtools::analysis::analyze_temperature --
#
# Extract the temperature from espresso
# Author: Zun-Jing Wang 
# Dec. 2008 

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::temperature {
    variable av_temperature 0.0
    variable av_temperature_i 0
    variable f_tvstemp

    namespace export setup_temperature
    namespace export analyze_temperature
    namespace export printav_temperature
    namespace export resetav_temperature
}

proc ::cgtools::analysis::temperature::resetav_temperature { } {
    variable av_temperature 
    variable av_temperature_i 
    set av_temperature 0.0 
    set av_temperature_i 0
}

proc ::cgtools::analysis::temperature::printav_temperature { } {
    variable av_temperature_i
    variable av_temperature
    global  ::cgtools::analysis::outputdir
    global ::cgtools::analysis::suffix
    variable f_tvstemp
    global ::cgtools::analysis::time
    
    if { $av_temperature_i > 0 } {
	set avtemp [expr $av_temperature/$av_temperature_i]
    	set f_tvstemp [open "$outputdir/time_vs_temperature$suffix" a]
   	puts $f_tvstemp "$time $avtemp"
   	close $f_tvstemp
    } else {
	::mmsg::warn [namespace current] "can't print temperature"
	flush stdout
    }
    #flush $f_tvstemp
}

proc ::cgtools::analysis::temperature::setup_temperature { } {
    global  ::cgtools::analysis::outputdir
    variable f_tvstemp
    global ::cgtools::analysis::iotype
    global ::cgtools::analysis::suffix
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_temperature$suffix "
    
    if { [file exists "$outputdir/time_vs_temperature$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvstemp [open "$outputdir/time_vs_temperature$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvstemp "\# Time Temperature"
    }
    close $f_tvstemp
}

proc ::cgtools::analysis::temperature::analyze_temperature {  } {
    #::mmsg::send [namespace current] "analyzing temperature"
    variable av_temperature_i
    variable av_temperature

    set deg_free 3 
    if { [regexp "ROTATION" [code_info]] } {
	set deg_free 6
    } 

    set npart [setmd n_part]
    #puts "deg_free is $deg_free, npart is $npart"

    set inst_temperature [expr [analyze energy kinetic]/(($deg_free/2.0)*$npart)]
    #puts "inst_temperature is $inst_temperature"

    set av_temperature [expr $av_temperature + $inst_temperature]

    incr av_temperature_i

    ::mmsg::debug [namespace current] "done"
}
