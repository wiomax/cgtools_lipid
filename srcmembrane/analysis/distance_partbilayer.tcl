#
# Proceedures for compute order parmeter: distance between residue-particle and bilayer-membrane
# Author: Zun-Jing Wang
# 2010 June 
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::distance_partbilayer {
    variable dis_partbilayer
    variable f_dispb
    namespace export printav_distance_partbilayer
    namespace export setup_distance_partbilayer
    namespace export analyze_distance_partbilayer
    namespace export resetav_distance_partbilayer
}

proc ::cgtools::analysis::distance_partbilayer::resetav_distance_partbilayer { } {
    # Do nothing 
}

proc ::cgtools::analysis::distance_partbilayer::printav_distance_partbilayer { } {
    # Do nothing
}

proc ::cgtools::analysis::distance_partbilayer::setup_distance_partbilayer { args } {
    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::iotype
    variable f_dispb

    if { [file exists "$outputdir/distance_partbilayer.dat"] } {
        set newfile 0
    } else {
        set newfile 1
    }
    set f_dispb [open "$outputdir/distance_partbilayer.dat" $iotype]
    if { $newfile || $iotype == "w"} {
        puts $f_dispb "\#distance_part1bilayer distance_part2bilayer"
    }
    close $f_dispb
 
}

proc ::cgtools::analysis::distance_partbilayer::analyze_distance_partbilayer { } {
    variable ::cgtools::analysis::topology
    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::iotype
    variable f_dispb

    set memcomz [::cgtools::utils::compute_membrane_comz $topology]
    set ibead 0
    foreach mol $topology {
      set moltype [lindex $mol 0]
      set typeinfo [::cgtools::utils::matchtype $moltype]
      set molname [lindex $typeinfo 1]
      if {$molname == "PART"} { 
    	set partbondlists [lindex $typeinfo 2]
    	set partbondtypelists [lindex $typeinfo 3]
    	set partbondcharmmlists [lindex $typeinfo 4]
    
    	#set and write particle list
    	set beadlists [lindex $partbondlists 0]
    	set nbeads [llength $beadlists]
    	set beadcharmmlists [lindex $partbondcharmmlists 0]
	for { set b 0 } { $b < $nbeads } {incr b } {
        	set partnum [lindex $mol [ expr $b + 1] ]
	        #puts "partnum= $partnum"
		set posvec [part $partnum print pos]
		set posz [lindex $posvec 2]
		set distance [expr abs($posz - $memcomz) ]
		if {$ibead == 0} {
			set dis_partbilayer $distance
		} else {
			lappend dis_partbilayer $distance
		}
		incr ibead		
	}
      } 
    } 
    set f_dispb [open "$outputdir/distance_partbilayer.dat" a]
    puts $f_dispb "$dis_partbilayer"
    close $f_dispb

    return $dis_partbilayer
}

