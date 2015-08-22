# ::cgtools::analysis::analyze_box_len --
#
# Extract the box dimensions from espresso
# Author: Zun-Jing Wang
# 2008 Nov.

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::boxl {
    variable av_boxl { 0 0 0 }
    variable av_boxl_i 0
    variable f_tvsbl
    variable verbose

    namespace export setup_boxl
    namespace export analyze_boxl
    namespace export printav_boxl
    namespace export resetav_boxl
}

proc ::cgtools::analysis::boxl::resetav_boxl { } {
    variable av_boxl 
    variable av_boxl_i 
    set av_boxl { 0.0 0.0 0.0 }
    set av_boxl_i 0
}

proc ::cgtools::analysis::boxl::printav_boxl { } {
    variable av_boxl_i
    variable av_boxl
    global  ::cgtools::analysis::outputdir
    global ::cgtools::analysis::suffix
    variable f_tvsbl
    global ::cgtools::analysis::time
    
    if { $av_boxl_i > 0 } {
	set avblx [expr [lindex $av_boxl 0]/($av_boxl_i*1.0)]
	set avbly [expr [lindex $av_boxl 1]/($av_boxl_i*1.0)]
	set avblz [expr [lindex $av_boxl 2]/($av_boxl_i*1.0)]
    	set f_tvsbl [open "$outputdir/time_vs_boxl$suffix" a]
	puts $f_tvsbl "$time $avblx $avbly $avblz"
    	close $f_tvsbl
    } else {
	::mmsg::warn [namespace current] "can't print average box length"
	flush stdout
    }
    #flush $f_tvsbl
}

proc ::cgtools::analysis::boxl::setup_boxl { args } {
    global  ::cgtools::analysis::outputdir
    variable f_tvsbl
    global ::cgtools::analysis::iotype
    global ::cgtools::analysis::suffix
    variable verbose
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_boxl$suffix "
    
     set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_boxl verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_boxl$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvsbl [open "$outputdir/time_vs_boxl$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsbl "\# Time boxx boxy boxz"
    }
    close $f_tvsbl

}

proc ::cgtools::analysis::boxl::analyze_boxl {  } {
    #::mmsg::send [namespace current] "analyzing box l"
    variable av_boxl_i
    variable av_boxl
    variable verbose

    set inst_boxl [setmd box_l]
    lset av_boxl 0 [expr [lindex $av_boxl 0] + [lindex $inst_boxl 0] ]
    lset av_boxl 1 [expr [lindex $av_boxl 1] + [lindex $inst_boxl 1] ]
    lset av_boxl 2 [expr [lindex $av_boxl 2] + [lindex $inst_boxl 2] ]

    incr av_boxl_i

    if { $verbose } {
	set avblx [expr [lindex $av_boxl 0]/($av_boxl_i*1.0)]
	set avbly [expr [lindex $av_boxl 1]/($av_boxl_i*1.0)]
	set avblz [expr [lindex $av_boxl 2]/($av_boxl_i*1.0)]
	::mmsg::send [namespace current]  "L: [lindex $inst_boxl 0] [lindex $inst_boxl 1] [lindex $inst_boxl 2] :: <L> $avblx $avbly $avblz"
	flush stdout
    }
    ::mmsg::debug [namespace current] "done"
}

