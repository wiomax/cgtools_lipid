#
# Routines for placing Positions, Bonds, Angles, Dihedrals of a molecule 
#
# Author: Zun-Jing Wang
# Sep.26 - Oct 1 2008 Done 


namespace eval cgtools::generation {}

# ::cgtools::generation::placemol-- 
# general routine for placing molecules
# Arguments:
# mol : list of atoms in the molecule
# partlist : list of particle information of the atoms 
proc ::cgtools::generation::placemol { mol partlist args } { 

    set options {
	{changepos.arg "0" "if change positions the particles"  }
	{move.arg   "{ 0. 0. 0. }" "movement vector of the center of mass of mol"  }
	{rotate.arg  "{ { 1. 0. 0.} { 0. 1. 0.} { 0. 0. 1.} }" "rotation matrix for the mol " }
    }
    set usage "Usage: create_bilayer topo boxl \[changepos:move:rotate]"
    array set params [::cmdline::getoptions args $options $usage]
    

    # Retrieve the molecule information for this molecule type	
    set typekey [::cgtools::utils::matchtype [lindex $mol 0]]

    # Place the molecule depending on type
    switch [lindex $typekey 1] {
	"POPC" {
		if {$params(changepos)==0} {
	    		place_lipid $mol $partlist
		}
		if {$params(changepos)!=0} {
	        	set partlist_new [place_posmoverotation $partlist $params(move) $params(rotate)]
	    		place_lipid $mol $partlist_new 
		}
	}
	"DOPC" {
		if {$params(changepos)==0} {
	    		place_lipid $mol $partlist
		}
		if {$params(changepos)!=0} {
	        	set partlist_new [place_posmoverotation $partlist $params(move) $params(rotate)]
	    		place_lipid $mol $partlist_new 
		}
	}
	"DPPC" {
		if {$params(changepos)==0} {
	    		place_lipid $mol $partlist
		}
		if {$params(changepos)!=0} {
	        	set partlist_new [place_posmoverotation $partlist $params(move) $params(rotate)]
	    		place_lipid $mol $partlist_new 
		}
	}
	"PART" {
		if {$params(changepos)==0} {
	    		place_part $mol $partlist
		}
		if {$params(changepos)!=0} {
	        	set partlist_new [place_posmoverotation $partlist $params(move) $params(rotate)]
	    		place_part $mol $partlist_new 
		}
	}
	"default" {
	    ::mmsg::err [namespace current] "couldn't place molecule of type [lindex $typekey 1], possibilities are: \n lipid \n hollowsphere \n sphericalconstraint"
	}
    }

    return
}

# ::cgtools::generation::place_lipid-- 
# Place positions of particles and bonds between particles of a lipid molecule 
# Arguments:
# mol : list of atoms in the molecule
# partlist : list of particle information of the atoms 
# Note that both linkbond and bendbond will have previously been set by set_bonded_interactions
proc ::cgtools::generation::place_lipid { mol partlist } {

    set moltype [lindex $mol 0]
    set typeinfo [::cgtools::utils::matchtype $moltype]
   
    set partbondlists [lindex $typeinfo 2]
    set partbondtypelists [lindex $typeinfo 3]

    #particle positions
    set beadlists [lindex $partbondlists 0]
    set nbeads [llength $beadlists]

    set beadtypelists [lindex $partbondtypelists 0]

    for { set b 0 } { $b < $nbeads } {incr b } {

	#current positions of particles
	set partnum [lindex $mol [ expr $b + 1] ]
	set parttype [lindex $beadlists $b]
	set parttypeinfo [lindex $beadtypelists $parttype]
	set partmass [lindex $parttypeinfo 2]

	set curpart [lindex $partlist $b]
	set curpos [lindex $curpart 1]
	set posx [lindex $curpos 0]
	set posy [lindex $curpos 1]
	set posz [lindex $curpos 2]
	
	part $partnum pos $posx $posy $posz type $parttype mass $partmass
	#puts "part $partnum pos $posx $posy $posz type $parttype mass $partmass"
    }

    #bonds
    set bondlists [lindex $partbondlists 1]
    set nbonds [llength $bondlists]
    for { set b 0 } { $b < $nbonds } {incr b } {
	set curbond [lindex $bondlists $b ]
    	# bond type i.e. bondid 
	set btype [lindex $curbond 0]
    	# index of the particles inside mol 
	set partlists_inmol [lindex $curbond 1]
	set part1_inmol [lindex $partlists_inmol 0]
	set part2_inmol [lindex $partlists_inmol 1]
	# pid [0:npart-1] to link the bond 
	set partnum1 [lindex $mol [ expr $part1_inmol + 1] ]
	set partnum2 [lindex $mol [ expr $part2_inmol + 1] ]
	 
	part $partnum2 bond $btype $partnum1
        #puts "part $partnum2 bond $btype $partnum1"
    }
	
    #angl 
    set angllists [lindex $partbondlists 2]
    set nangls [llength $angllists]
    for { set b 0 } { $b < $nangls } {incr b } {
	set curbond [lindex $angllists $b ]
    	# bond type i.e. bondid 
	set btype [lindex $curbond 0]
    	# index of the particles inside mol 
	set partlists_inmol [lindex $curbond 1]
	set part1_inmol [lindex $partlists_inmol 0]
	set part2_inmol [lindex $partlists_inmol 1]
	set part3_inmol [lindex $partlists_inmol 2]
	# pid [0:npart-1] to link the angle 
	set partnum1 [lindex $mol [ expr $part1_inmol + 1] ]
	set partnum2 [lindex $mol [ expr $part2_inmol + 1] ]
	set partnum3 [lindex $mol [ expr $part3_inmol + 1] ]
	 
	part $partnum2 bond $btype $partnum1 $partnum3
        #puts "$btype $partnum1 $partnum2 $partnum3"
    }
    
    #dihe 
    set dihelists [lindex $partbondlists 3]
    set ndihes [llength $dihelists]
    for { set b 0 } { $b < $ndihes } {incr b } {
	set curbond [lindex $dihelists $b ]
    	# bond type i.e. bondid 
	set btype [lindex $curbond 0]
    	# index of the particles inside mol 
	set partlists_inmol [lindex $curbond 1]
	set part1_inmol [lindex $partlists_inmol 0]
	set part2_inmol [lindex $partlists_inmol 1]
	set part3_inmol [lindex $partlists_inmol 2]
	set part4_inmol [lindex $partlists_inmol 3]
	# pid [0:npart-1] to link the dihedral
	set partnum1 [lindex $mol [ expr $part1_inmol + 1] ]
	set partnum2 [lindex $mol [ expr $part2_inmol + 1] ]
	set partnum3 [lindex $mol [ expr $part3_inmol + 1] ]
	set partnum4 [lindex $mol [ expr $part4_inmol + 1] ]
	 
	part $partnum2 bond $btype $partnum1 $partnum3 $partnum4
        #puts "$btype $partnum1 $partnum2 $partnum3 $partnum4"
    }
}

proc ::cgtools::generation::place_posmoverotation {partlist move rotate} {
    	#::mmsg::send [namespace current] "entering place_posmoverotation"
	
	set npart [llength $partlist]

	set partlist_new 0	
	unset partlist_new
    	for { set i 0 } { $i < $npart } {incr i } {
		set pos_old [lindex [lindex $partlist $i] 1]
		set pos_move [::cgtools::utils::add_vecs $pos_old $move]
 		set pos_rotate [::cgtools::utils::matrix_vec_multiply $rotate $pos_move]

	        set ninfo [llength [lindex $partlist $i]]
		set poslist_new 0	
		unset poslist_new
    		for { set j 0 } { $j < $ninfo } {incr j } {
    			if { $j!=1 } {
				lappend poslist_new [lindex [lindex $partlist $i] $j]
			}
    			if { $j==1 } {
				lappend poslist_new $pos_rotate
			}
		}
		lappend partlist_new $poslist_new
	}

	return $partlist_new
}

# ::cgtools::generation::place_part-- 
# Place positions of particles for sigle-part mol
# Arguments:
# mol : list of atoms in the molecule
# partlist : list of particle information of the atoms 
# Note that both linkbond and bendbond will have previously been set by set_bonded_interactions
proc ::cgtools::generation::place_part { mol partlist } {

    set moltype [lindex $mol 0]
    set typeinfo [::cgtools::utils::matchtype $moltype]
    #puts "typeinfo= $typeinfo"
   
    set partbondlists [lindex $typeinfo 2]
    set partbondtypelists [lindex $typeinfo 3]

    #particle positions
    set beadlists [lindex $partbondlists 0]
    set nbeads [llength $beadlists]
    set beadtypelists [lindex $partbondtypelists 0]
    set itype_begin [lindex [lindex $beadtypelists 0] 0]
    #puts "itype_begin= $itype_begin"
    #puts "beadlists= $beadlists"
    #puts "nbeads= $nbeads"
    #puts "beadtypelists= $beadtypelists"

    for { set b 0 } { $b < $nbeads } {incr b } {

	#current positions of particles
	set partnum [lindex $mol [ expr $b + 1] ]
	set parttype [lindex $beadlists $b]
	#puts "parttype= $parttype"
	set parttypeinfo [lindex $beadtypelists [expr $parttype - $itype_begin]]
	set partmass [lindex $parttypeinfo 2]

	set curpart [lindex $partlist $b]
	set curpos [lindex $curpart 1]
	set posx [lindex $curpos 0]
	set posy [lindex $curpos 1]
	set posz [lindex $curpos 2]
	
	part $partnum pos $posx $posy $posz type $parttype mass $partmass
	#puts "part $partnum pos $posx $posy $posz type $parttype mass $partmass"
	#exit
    }
}
