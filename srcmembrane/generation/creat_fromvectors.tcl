# 
# Routines for generating system fromvectors
# Author: Zun-Jing Wang
# June 20-25 2010

namespace eval cgtools::generation {}

namespace eval ::cgtools::generation::fromvectors {
    namespace export create_fromvectors
}

#::cgtools::generation::fromvectors::create_fromvectors
#
# Place the initial configuration of particles from vectors
#
# Arguments: 
# topo: The topology in espresso format. It should be
# sorted in ascending order according to molecule type (not bead id).
#
proc ::cgtools::generation::fromvectors::create_fromvectors { args } {
    variable ::cgtools::generation::topology

    ::mmsg::send [namespace current] "placing initial positions from vectorlist"
    set options {
	{readfile.arg  ""    "vectorlist containing initial particle positions"}
    }
    set usage "Usage: create_fromvectors \[readfile] "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]

    # set all the particle positions from $params(readfile)
    set nvect [llength $params(readfile)]
    for {set ivect 0 } { $ivect <  $nvect} { incr ivect} {
    	lappend vectors [lindex $params(readfile) $ivect]
    }
    lappend vectorlist $vectors
    #puts "::cgtools::generation::fromvectors::create_fromvectors::nvect =  $nvect"
    #puts "::cgtools::generation::fromvectors::create_fromvectors::vectorlist =  $vectorlist"
    set mollists [::cgtools::generation::placeparticles_vectors $vectorlist]
    #puts "Zunjing Zunjing Zunjing Zunjing Zunjing Zunjing Zunjing"

    # number of molecules
    set imol 0
    foreach mol $topology {
	set partlist [lindex $mollists $imol]
        incr imol 
    	#puts "mol= $mol"
	#puts "partlist= $partlist"
	::cgtools::generation::placemol $mol $partlist -changepos 0 
    }
    
    # Check particle consistency
    if { [setmd n_part] != [expr [::cgtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::cgtools::utils::maxpartid $topology] +1] were specified in topology "
    }
}
