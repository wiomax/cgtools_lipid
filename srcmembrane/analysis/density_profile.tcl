#
# Proceedures for determining the vertical density profile of a bilayer
# Author: Zun-Jing Wang
# 2008 Nov.
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::density_profile {
    variable av_densities
    variable av_densities_i 0
    variable nogrid
    variable verbose
    variable nbins
    variable hrange
    variable beadtypes
    variable colloidmoltypes
    variable colloidradii
    variable density_cntfile
    namespace export printav_density_profile
    namespace export setup_density_profile
    namespace export analyze_density_profile
    namespace export resetav_density_profile
}

proc ::cgtools::analysis::density_profile::resetav_density_profile { } {
    # Do nothing because we want to average this over the entire simulation
    
}

proc ::cgtools::analysis::density_profile::printav_density_profile { } {
    variable f_densityprof
    variable av_densities_i
    variable beadtypes
    variable av_densities
    variable hrange
    variable nbins
    variable density_cntfile
    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::suffix

    set binwidth [expr $hrange*2.0/(1.0*$nbins)]
    if {  $av_densities_i > 0 } {
	
	set f_densityprof [open "$outputdir/av_zprof$suffix" "w" ]
	# Write a header to the file
	puts -nonewline $f_densityprof "\# zheight "
	foreach bt $beadtypes {
	    puts -nonewline $f_densityprof "|| $bt up "
	}
	for { set bt [expr [llength $beadtypes] -1 ] } { $bt >= 0 } { incr bt -1 } {
	    puts -nonewline $f_densityprof "|| $bt down "
	}
	puts $f_densityprof ""
	
	# Write out the density profiles to a file
	for { set bin 0 } { $bin < [llength $av_densities] } { incr bin } {
	    set currbin [lindex $av_densities $bin]
	    puts -nonewline $f_densityprof "[expr $bin*$binwidth+($binwidth/2.0)] "
	    for { set bt 0 } { $bt < [llength $currbin] } { incr bt } {
		puts -nonewline $f_densityprof "[expr [lindex $currbin $bt]/(1.0*$av_densities_i)] "
	    }
	    puts $f_densityprof ""
	}		    
	close $f_densityprof

         # backup <$outputdir/av_zprof$suffix> 
        incr density_cntfile 
        set newfile 0
        while { $newfile == 0 } {
                if { [file exists "$outputdir/av_zprof$suffix.$density_cntfile"] } {
                        incr density_cntfile
                } else {
                        set newfile 1
                }
        }
        catch { exec cp $outputdir/av_zprof$suffix  $outputdir/av_zprof$suffix.$density_cntfile }

    }
}

proc ::cgtools::analysis::density_profile::setup_density_profile { args } {

    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::suffix
    global ::cgtools::analysis::iotype
    variable nogrid
    variable nbins
    variable hrange
    variable beadtypes
    variable verbose
    variable colloidradii
    variable colloidmoltypes

    variable av_densities
    variable av_densities_i
    variable density_cntfile

    set options {
	{nogrid  "dont use height grid for density profile analysis" }
	{verbose "print out lots of stuff" }
	{nbins.arg "100" "Number of bins used for density profile analysis" }
	{hrange.arg "50.0" "Range over which to calculate density profile"}
	{beadtypes.arg "0" "Identity of beads to use for density profile analysis" }
	{colloidmoltypes.arg "0" "Molecule types assigned to colloids"}
	{r.arg "0" "radii of colloids" }
    }
    set usage "Usage: setup_density_profile nogrid:verbose:nbins:hrange:beadtypes:coloidmoltypes:colloidradii "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)
    set beadtypes $params(beadtypes)
    set nbins $params(nbins)
    set hrange $params(hrange)
    set nogrid $params(nogrid)
    set colloidmoltypes $params(colloidmoltypes)
    set colloidradii $params(r)

    set density_cntfile 0
    #Initialize av_densities
    set thisbinlist 0.0
    unset thisbinlist
    for { set bn 0 } { $bn < $nbins } { incr bn } {
	for { set bt 0 } { $bt < [expr 2*[llength $beadtypes]] } { incr bt } {
	    lappend thisbinlist 0.0
	}
	lappend av_densities $thisbinlist
	unset thisbinlist
	
    }
    mmsg::send [namespace current] "setup density profile with for beads $beadtypes with $nbins bins and $hrange range"

}

# ::cgtools::analysis::analyze_density_profile --
# 
# Calculates the density of each lipid bead type as a function of
# vertical distance relative to the bilayer midplane
#
proc ::cgtools::analysis::density_profile::analyze_density_profile { } {
    variable ::cgtools::analysis::topology

    variable av_densities
    variable av_densities_i
    variable hrange
    variable beadtypes
    variable nbins
    variable nogrid
    variable colloidmoltypes
    variable colloidradii
    variable verbose

    ::mmsg::send [namespace current] "analyzing density profile"

    # First check to see if we have a spherical situation
    if { $colloidradii != 0 } {

	if { [llength $colloidmoltypes] == 0 } {
	    # Special case for floppy vesicles
	    set colloidcenter [analyze centermass 0 ]
	} else {
	    set colloidcentersmass [::cgtools::utils::calc_centersofmass_bymoltype $colloidmoltypes $topology]

	    # At the moment we support just one colloid but we could in principle support many
	    if { [llength colloidcentersmass] > 1 || [llength $colloidradii] > 1 } {
		::mmsg::err [namespace current] "no support for multiple colloids yet in analyze density profile"
	    } else {
		set colloidcenter [lindex $colloidcentersmass 0]
	    }

	    set colloidcenter [list [lindex $colloidcenter 0 0]  [lindex $colloidcenter 0 1] [lindex $colloidcenter 0 2] ]
	}

	set densities [analyze bilayer_density_profile $hrange $nbins $beadtypes withsphere $colloidradii [lindex $colloidcenter 0] [lindex $colloidcenter 1] [lindex $colloidcenter 2] ]
    } else {
	# Normal situation No COlloid

	if { $nogrid } {
	    puts "hrange = $hrange,   nbins = $nbins"
	    set densities [ analyze bilayer_density_profile $hrange $nbins $beadtypes nogrid ]
	} else {
	    set densities [ analyze bilayer_density_profile $hrange $nbins $beadtypes  ]
	}
    }

    # Select the actual data and skip the label
    set densities [lindex $densities 1]

    # Bin up the data
    for { set bn 0 } { $bn < $nbins } { incr bn } {
	for { set bt 0 } { $bt < [expr 2*[llength $beadtypes]] } { incr bt } { 
	    lset av_densities $bn $bt [ expr [lindex $densities $bt $bn] + [lindex $av_densities $bn $bt] ]
	}
    }
    incr av_densities_i

}


