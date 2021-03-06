# 
#  Proceedures for creating various types of topologies.  In this case
#  topology refers to the espresso feature with the same name and
#  describes the groupings of atoms into molecules.
#
# Obtained from MBTools, modified by Zun-Jing Wang
# Sep. 2008 Done 


namespace eval cgtools::generation {}


#----------------------------------------------------------#
# ::cgtools::generation::create_simple_topo --
#
# Creates a basic topology from a single molecule type
#
proc ::cgtools::generation::create_simple_topo {n_mols beads_per_mol  args } {

    set topo 0 
    unset topo

    set options {
	{moltype.arg 0 "type of molecule to create" }
	{startpart.arg 0 "The starting id for particles" }
    }
    set usage "Usage: ceate_simple_topo n_mols beads_per_mol \[startmol:startpart:]"
    array set params [::cmdline::getoptions args $options $usage]

    ::mmsg::send [namespace current] "creating a topology for moltype $params(moltype)" 

    # First we create a very simple completely ordered and sorted
    # topology 
    set mol_type $params(moltype)
    set partid $params(startpart)

    if { $n_mols <= 0 } {
	return 0
    }

    for { set i 0 } { $i < $n_mols } { incr i } {
	# construct the topo_entry for this lipid

	set topo_entry $mol_type
	for { set j 0 } { $j < $beads_per_mol } { incr j } {
	    lappend topo_entry $partid
	    incr partid
	}
	lappend topo $topo_entry
    }
    return $topo

}


# 
# ::cgtools::generation::shuffle_topo--
#  
# Performs a shuffle on the topology so that molecules are no longer
# in blocks of the same type
#
#
# Arguments:
#
# topo: The topology to be shuffled.  
# 

proc ::cgtools::generation::shuffle_topo {topo} {
    mmsg::send [namespace current] "shuffling topology" 

    set typetemplate 0
    unset typetemplate

    set shuffledtopo 0
    unset shuffledtopo

    set startmol [::cgtools::utils::minmoltype $topo ]
    set startpart [::cgtools::utils::minpartid $topo ]

    set full_n_molslist [::cgtools::utils::listnmols $topo]
    set n_molstotal [llength $topo]

    set remaining_molslist $full_n_molslist

    # Check that we have there are no discontinuities in the particle ids

    # To do the shuffle properly we really need to completely regenerate the topology


    # The first step is to construct a list with randomly allocated
    # molecule types to serve as a template
    set n_remaining $n_molstotal
    mmsg::send [namespace current] "constructing topo template " nonewline 
    flush stdout
    set dotfreq [expr int(floor($n_molstotal/10.0))]
    if { $dotfreq < 1 } { set dotfreq 1}
    for { set i 0 } { $i < $n_molstotal } { incr i } {
	
	# According to the molecule proportions determine the type of
	# the next molecule
	set lpick [expr $n_remaining*[t_random]]
	for {set lnum 0 } { $lnum <  [llength $remaining_molslist] } { incr lnum } {
	    set lpick [expr $lpick - [lindex $remaining_molslist $lnum 1]]

		if  { $lpick <= 0 } {
		    set thislipidtype [lindex $remaining_molslist $lnum 0]
		    # Subtract the picked lipid from our molslist
		    set n_before [lindex  $remaining_molslist $lnum 1]
		    lset remaining_molslist $lnum 1 [expr  $n_before -1]
		    set n_remaining [expr $n_remaining -1 ]
		    break
		}
	}

	lappend typetemplate $thislipidtype

        if { $i%$dotfreq ==0 } {
	    mmsg::send [namespace current] "." nonewline
	    flush stdout
	}
    }
    mmsg::send [namespace current] "done" 

    # We need a list telling us how many atoms are in each molecule
    set molsizes [::cgtools::utils::listmollengths $topo ]
    set pnum $startpart
    foreach type $typetemplate {
	foreach t $molsizes {
	    if { [lindex $t 0] == $type  } {
		set thislen [lindex $t 1]
	    }
	}
	set thismol $type
	for { set p 0 } { $p < $thislen } { incr p } {
	    lappend thismol $pnum
	    incr pnum
	}
	lappend shuffledtopo $thismol
    }
 

    flush stdout
    return $shuffledtopo
}

#
# ::cgtools::generation::sort_topo -- 
#
# Sort a topology into ascending molecule type order
#
#
proc ::cgtools::generation::sort_topo { topo } {
    set sortedtopo 0
    unset sortedtopo 
   
    set n_molstotal [llength $topo]
    set maxtp [::cgtools::utils::maxmoltypeid $topo]

    for { set l 0 } { $l <= $maxtp } { incr l } {
	for { set i 0 } { $i < $n_molstotal } { incr i } {
	    if { [lindex [lindex $topo $i] 0 ] == $l } {
		lappend sortedtopo [lindex $topo $i]
	    }
	}
    }
    return $sortedtopo
}


#
# ::cgtools::generation::join_topos--
#
# Join two topologies
#
proc ::cgtools::generation::join_topos { topo1 topo2 args } {
    mmsg::send [namespace current] "joining topologies"

    set joinedtopo 0
    unset joinedtopo

    set options {
	{sortmolid  "sort by molid" }
    }
    set usage "Usage: join_topos \[sortmolid]"
    array set params [::cmdline::getoptions args $options $usage]
    foreach mol $topo1 {
	lappend joinedtopo $mol
    }
    foreach mol $topo2 {
	lappend joinedtopo $mol
    }

    return $joinedtopo
}

