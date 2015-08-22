# 
# Routines for generating random system
# Author: Zun-Jing Wang
# Sep.26 - Oct 1 2008 Done 
#

namespace eval cgtools::generation {}

namespace eval ::cgtools::generation::random {
    namespace export create_random
}

#::cgtools::generation::random::create_random
#
# Now the code only works for single component system 
#
# Arguments: 
# topo: The topology in espresso format. It should be
# sorted in ascending order according to molecule type (not bead id).
#
proc ::cgtools::generation::random::create_random { args } {

    variable ::cgtools::generation::boxl
    variable ::cgtools::generation::topology
    
    ::mmsg::send [namespace current] "randomly placing initial positions of lipids"
    set options {
	{readfile.arg  ""    "readpdb file name containing initial particle positions"}
    }
    set usage "Usage: create_random \[readfile] "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]

    # Read one mol from pdb file $params(readfile)
    set mollists_template [::cgtools::generation::placeparticles_template $params(readfile)]
    #::mmsg::send [namespace current] "$mollists_template"
    set partlist_template [lindex $mollists_template 0]
    set moltype_template [lindex [lindex $partlist_template 0] 4]
    #::mmsg::send [namespace current] "$moltype_template"

    set bx [lindex $boxl 0]
    set by [lindex $boxl 1]
    set bz [lindex $boxl 2]

    # number of molecules
    set imol 0
    #::mmsg::send [namespace current] "$topology"
    foreach mol $topology {
    	set moltype [lindex $mol 0]
	if {$moltype != $moltype_template} {
		mmsg::err [namespace current] "molecular type read from pdb file $$params(readfile) is not right "
	}
	if {$moltype == $moltype_template} {
		#mass of center of the molecule_template
		set pos_sum { 0. 0. 0. }
	    	foreach curpart $partlist_template {
    		   set curpos [lindex $curpart 1]
		   set pos_add [::cgtools::utils::add_vecs $pos_sum $curpos]
   	       	   set pos_sum $pos_add
		}
		set nbead [llength $partlist_template]
		set svalue [expr -1.0/($nbead * 1.0)] 
		set pos_old [::cgtools::utils::scalevec $pos_sum $svalue] 
    	#::mmsg::send [namespace current] "$pos_old"

		#mass of center of the molecule setted up
		set pos_new 0
		unset pos_new
    		set curx [expr $bx*([t_random]-0.5)]
    		set cury [expr $by*([t_random]-0.5)]
    		set curz [expr $bz*([t_random]-0.5)]
		lappend pos_new $curx
		lappend pos_new $cury
		lappend pos_new $curz
		#vector in moving center of mass
		set mvector [::cgtools::utils::add_vecs $pos_new $pos_old]
    		#::mmsg::send [namespace current] "$mvector"

		#vector in rotating the molecule
		set rmatrix 0
		unset rmatrix
		set theta [expr 3.14*[t_random]]
		set phi [expr 3.14*2.0*[t_random]]
		set psi [expr 3.14*2.0*[t_random]]
		set rmatrix [::cgtools::utils::get_matrix_rotation $theta $phi $psi]
    		#::mmsg::send [namespace current] "$rmatrix"
		
		#set one molecule
		::cgtools::generation::placemol $mol $partlist_template -changepos 1 -move $mvector -rotate $rmatrix
    		#::mmsg::send [namespace current] "hello, seting one mol is successful"
	}
    }
    
    # Check particle consistency
    if { [setmd n_part] != [expr [::cgtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::cgtools::utils::maxpartid $topology] +1] were specified in topology "
    }
}
				      


