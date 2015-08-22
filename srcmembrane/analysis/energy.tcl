# ::cgtools::analysis::analyze_energy --
#
# Calculate the total energy of the system and break it into components
# Author: Zun-Jing Wang 
# Jan 2009   

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::energy {
    variable av_components_en 0
    variable av_total_en 0
    variable av_kin_en 0
    variable av_overlap_en 0
    variable av_nb_en 0
    variable av_en_i 0

    variable f_tvsen
    variable verbose

    namespace export printav_energy
    namespace export setup_energy
    namespace export analyze_energy
    namespace export resetav_energy
}

proc ::cgtools::analysis::energy::resetav_energy { } {
    variable av_components_en 
    variable av_total_en 
    variable av_kin_en 
    variable av_overlap_en 
    variable av_nb_en 
    variable av_en_i 
    for {set i 0 } {$i < [llength $av_components_en] } { incr i } {
	lset av_components_en $i 0.0
    }
    set av_total_en 0
    set av_kin_en 0
    set av_overlap_en 0
    set av_nb_en 0
    set av_en_i 0
}

proc ::cgtools::analysis::energy::printav_energy { } {
    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::suffix
    variable av_components_en 
    variable av_total_en 
    variable av_kin_en 
    variable av_overlap_en 
    variable av_nb_en 
    variable av_en_i     
    variable f_tvsen
    global ::cgtools::analysis::time

    if { $av_en_i > 0 } {
        set f_tvsen [open "$outputdir/time_vs_energy$suffix" a ]
	puts -nonewline $f_tvsen "$time [expr $av_total_en/(1.0*$av_en_i)] [expr $av_kin_en/(1.0*$av_en_i)] [expr $av_overlap_en/(1.0*$av_en_i)] [expr $av_nb_en/(1.0*$av_en_i)]"
	foreach comp $av_components_en {
	    puts -nonewline $f_tvsen " [expr $comp/(1.0*$av_en_i)]"
	}
	puts $f_tvsen ""
        close $f_tvsen
    } else {
        ::mmsg::warn [namespace current] "can't print energy"
        flush stdout
    }
}

proc ::cgtools::analysis::energy::setup_energy { args } {
    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::suffix
    global ::cgtools::analysis::iotype
    variable f_tvsen
    variable verbose
    variable av_components_en

    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_energy verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_energy$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    ::mmsg::debug  [namespace current]  "opening $outputdir/time_vs_energy$suffix "
    set f_tvsen [open "$outputdir/time_vs_energy$suffix" $iotype ]
    # First work out the names of components
    set raw [analyze energy]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsen "\# Components of the energy "
	puts -nonewline $f_tvsen "\# Time total_en kinetic_en overlap_en nonbonded_en"
    }
    unset av_components_en
    for { set k 0 } { $k < [llength $raw ] } { incr k } {
	set tmp [lindex $raw $k]
	set ntmp [llength $tmp]	
	if { [ regexp "nonbonded" $tmp ] } {
	    puts -nonewline $f_tvsen " [lrange $tmp 0 end-1]"
	    lappend av_components_en 0.0
	}
    }
    close $f_tvsen
}

proc ::cgtools::analysis::energy::analyze_energy {  } {
    #::mmsg::send [namespace current] "analyzing energy"
    variable av_components_en
    variable av_total_en
    variable av_kin_en
    variable av_overlap_en
    variable av_nb_en
    variable av_en_i
    variable verbose
    set energy_all [analyze energy]
#    puts $energy_all
    set nb_en 0
    set nbcount 0
    for { set i 0 } { $i < [llength $energy_all ] } { incr i } {
	set tmp [lindex $energy_all $i]
	set ntmp [llength $tmp]	
	if { [ regexp "energy" $tmp ] } {
	    set total_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "kinetic" $tmp ] } {
	    set kin_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "overlapped" $tmp ] } {
	    set overlap_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "nonbonded" $tmp ] } {
	    set nb_en [expr $nb_en + [lindex $tmp [expr $ntmp -1]]]
#	    puts "$nbcount [llength $av_components_en]"
	    lset av_components_en $nbcount [expr [lindex $av_components_en $nbcount] + [lindex $tmp [expr $ntmp -1]] ]

	    incr nbcount
	}
	
    }
    incr av_en_i
    set av_total_en  [expr $av_total_en + $total_en/(1.0)]
    set av_kin_en  [expr $av_kin_en + $kin_en/(1.0)]
    set av_overlap_en  [expr $av_overlap_en + $overlap_en/(1.0)]
    set av_nb_en  [expr $av_nb_en + $nb_en/(1.0)]
    if { $verbose } {
	::mmsg::send [namespace current] "energy: [expr $av_total_en/($av_en_i*1.0)] kinetic: [expr $av_kin_en/($av_en_i*1.0)] OVERLAPPED: [expr $av_overlap_en/($av_en_i*1.0)] nonbonded: [expr $av_nb_en/($av_en_i*1.0)]"
    }
    
    ::mmsg::debug [namespace current] "done"

}


