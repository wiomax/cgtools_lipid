# 
#  Routines for generating system fromfile
#
# Author: Zun-Jing Wang
# Sep.26 - Oct 1 2008 Done 


namespace eval cgtools::generation {}

namespace eval ::cgtools::generation::fromfile {
    namespace export create_fromfile
}

#::cgtools::generation::fromfile::create_fromfile
#
# Place the initial configuration of lipids from a pdb file 
# Notice the pdb file must be consistent with the pregenerated topology.
#
# Arguments: 
# topo: The topology in espresso format. It should be
# sorted in ascending order according to molecule type (not bead id).
#
proc ::cgtools::generation::fromfile::create_fromfile { args } {
    variable ::cgtools::generation::topology

    ::mmsg::send [namespace current] "placing initial positions of lipids from a pdb file"
    set options {
	{readfile.arg  ""    "readpdb file name containing initial particle positions"}
    }
    set usage "Usage: create_fromfile \[readfile] "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]


    # Read all the particle positions from pdb file $params(readfile)
    
    set mollists [::cgtools::generation::placeparticles_all $params(readfile)]

    # number of molecules
    set imol 0
    foreach mol $topology {
	set partlist [lindex $mollists $imol]
        incr imol 
	::cgtools::generation::placemol $mol $partlist -changepos 0 
    }
    
    # Check particle consistency
    if { [setmd n_part] != [expr [::cgtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::cgtools::utils::maxpartid $topology] +1] were specified in topology "
    }
}
