# cgtools::generation --
#
# Generate molecular topology, parameter and positions
# 
# Author: Zun-Jing Wang
# Sep. 25 2008

package require ::mmsg 
package require ::cgtools::utils

package provide ::cgtools::generation 1.0.0

# Create the empty namespace to which we shall add our routines
namespace eval ::cgtools::generation {
    # Global routine for setting up the system
    namespace export generate_system

    # Global variables for setup routines
    variable topology
    variable boxl
    variable notopo

    variable trappedmols
    variable trappedmolsupdated
    variable userfixedparts
    variable interbead


    # Read in all the routines
    # Files with separate namespaces corresponding to geometries
    source [file join [file dirname [info script]] creat_fromfile.tcl]
    source [file join [file dirname [info script]] creat_random.tcl]
    source [file join [file dirname [info script]] creat_fromvectors.tcl]


    # Helper function files
    source [file join [file dirname [info script]] topologies.tcl]
    source [file join [file dirname [info script]] placemol.tcl]
    source [file join [file dirname [info script]] placeparticles.tcl]
}

# ::cgtools::generation::generate_system -- 
#
# A large scale wrapper routine for all of the system setup commands
# in this module.
# 
# This routine should allow flexible use of multiple geometrical
# structures and topologies in a simulation.  A complete topology and
# geometry specification is provided for each object (eg. singlemol,
# bilayer etc) in the system and these are grouped into the list
# <system_specs>.  generate_system loops through each of these objects
# setting each up independently and then combining topologies at the
# end.
#
# Arguments:
# 
# system_specs: A list of specified system objects.  Each object is
# itself a list of the form < geometry n_lipidslist beads_per_mol >
# where geometry may be any allowable geometry (eg singlemol, bilayer 
# etc). <n_lipidslist> is a list containing the number of molecules of
# each of the molecule types in this object and beads_per_mol
# specifies the number of particles in each of the molecule types.
#
# 
proc ::cgtools::generation::generate_system { system_specs iboxl } {
    ::mmsg::send [namespace current] "setting up system "

    # The molecule types spec should be globally accessible
    variable boxl
    variable topology 
    variable notopo

    set boxl $iboxl

    set topolist 0
    unset topolist

    set topologieslist 0
    unset topologieslist

    # Starting value for particle ids
    set startp 0
    # Current particle id
    set currpid $startp


    foreach spec $system_specs {
	#puts "::cgtools::generation::generate_system::spec = $spec"
	# Flags for tracking input settings
	set geometryset 0
	set n_molslistset 0
	set n_molslist 0
	set notopo 0
	foreach item $spec {
	    switch [lindex $item 0] {
		"geometry" {
		    set geometry [lindex $item 1]
		    set geometryreadfile [lindex $item 2]
		    set geometryset 1
		}
		"n_molslist" {
		    set n_molslist [lindex $item 1]
		    set n_molslistset 1
		}
		"default" {
		    ::mmsg::warn [namespace current] "unknown item [lindex $item 0] in system spec. allowed values are: \n geometry \n n_molslist  "
		}
	    }
	}


	if { !$geometryset } {
	    mmsg::err [namespace current] "geometry not specified"
	}
	if { !$n_molslistset } {
	    mmsg::err [namespace current] "n_molslist not specified for [lindex $geometry 0]"
	}

	# Generate a topology from a list of the number and size of
	# each molecule
	foreach mtp $n_molslist {
	    set thismoltypeid [lindex $mtp 0]
	    set nmols [lindex $mtp 1]
	    set tpspec [::cgtools::utils::matchtype [lindex $mtp 0]]
	    set nbeads_mol [llength [lindex [lindex $tpspec 2] 0]]
	    # Create the topology for this lipid type
	    set topo [create_simple_topo $nmols $nbeads_mol -moltype  $thismoltypeid -startpart $currpid ]		
		
	    # Just in case zero molecules were specified we need
	    # to check if topo was actually created at all by the
	    # last command
	    if { $topo == 0 } {
		::mmsg::err [namespace current] "no topo created for molecule type $thismoltypeid"
	    } else {
		lappend topolist $topo
		set currpid [expr [::cgtools::utils::maxpartid $topo ] + 1]
	    }

	}

	# Join all of the previously made topologies
	set first 1
	foreach topo $topolist {
	    if { $first } { 
		set topology $topo 
		set first 0
	    } else {
		set topology [join_topos $topology $topo]
	    }
	    
	}
	unset topolist


	# Now wrap the topology onto a specified geometry and perform any
	# other geometry specific tasks

	# Shuffle the topology
	#set topology [shuffle_topo $topology ]
	# Now run the creation command for the specified geometry
	set createprefix "create_"
	set namespaceprefix "::cgtools::generation::"
	# Construct the name of the create command
	set command $geometry
	set geometry [lindex [split $geometry " "] 0]
	set createcommand "$namespaceprefix$geometry\:\:$createprefix$command "
	::mmsg::debug [namespace current] "executing $command"
	#puts "::cgtools::generation::generate_system::geometryreadfile = $geometryreadfile"
    	#exit
	#eval $createcommand -readfile $geometryreadfile
	if { [catch  {eval $createcommand -readfile $geometryreadfile} errm ] } {
	    mmsg::err [namespace current] "couldn't execute creation command for $command \n $errm"
	}

	if {!$notopo} {
	    lappend topologieslist $topology
	}
	#puts "topology: $topology"
	
    }

    # Join all of the previously made topologies
    set first 1
    foreach topo $topologieslist {
	if { $first } { 
	    set topology $topo 
	    set first 0
	} else {
	    set topology [join_topos $topology $topo]
	}
	
    }

    #puts "topology= $topology"
    #set topology [sort_topo $topology]

    #puts "topology= $topology"
    return $topology

}

proc ::cgtools::generation::get_trappedmols {  } {
    variable trappedmols
    variable topology

    variable trappedmolsupdated

    if { [catch { set dum $trappedmols } ] } {
	::mmsg::warn [namespace current] "no trappedmols defined"
	return -1
    } else {
	if { !$trappedmolsupdated } {
	    set didntfindmol 1
	    for { set j 0 } { $j < [llength $trappedmols] } { incr j } {
		# Update trappedmols
		set fmol [lindex $trappedmols $j]
		for { set i 0 } { $i < [llength $topology] } { incr i } {
		    set mol [lindex $topology $i]
		    if { [lindex $fmol 0] == [lindex $mol 1] } {
			lset trappedmols $j 0 $i
			set didntfindmol 0
			break			
		    }
		}
		if { $didntfindmol } {
		    ::mmsg::err [namespace current] "could not get_trappedmols unable to find the corresponding particles"
		}

	    }
	}
	set trappedmolsupdated 1
	return $trappedmols
    }
}

proc ::cgtools::generation::get_userfixedparts {  } {
    variable userfixedparts

    if { [catch { set dum $userfixedparts } ] } {
	::mmsg::warn [namespace current] "no user fixed particles defined"
	return -1
    } else {
	return $userfixedparts
    }
}

proc ::cgtools::generation::get_interbead {  } {
    variable interbead
    return $interbead
}
