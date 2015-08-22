# ::cgtools::analysis::analyze_size_gelclusters_largest --
#
# Extract the box dimensions from espresso
# Author: Zun-Jing Wang
# 2009 August

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::gelclusters_largest {
    variable size_gelclusters_largest 
    variable f_cluster

    namespace export setup_gelclusters_largest
    namespace export analyze_gelclusters_largest
    namespace export printav_gelclusters_largest
    namespace export resetav_gelclusters_largest
}

proc ::cgtools::analysis::gelclusters_largest::resetav_gelclusters_largest { } {
    variable size_gelclusters_largest 
}

proc ::cgtools::analysis::gelclusters_largest::printav_gelclusters_largest { } {
    variable size_gelclusters_largest 
    global  ::cgtools::analysis::outputdir
    global ::cgtools::analysis::suffix
    variable f_cluster
    global ::cgtools::analysis::time
    
    set f_cluster [open "$outputdir/time_vs_gelclusters_largest$suffix" a]
    puts $f_cluster "$time $gelclusters_largest" 
    close $f_cluster
    #flush $f_cluster
}

proc ::cgtools::analysis::gelclusters_largest::setup_gelclusters_largest { args } {
    global  ::cgtools::analysis::outputdir
    variable f_cluster
    global ::cgtools::analysis::iotype
    global ::cgtools::analysis::suffix
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_gelclusters_largest$suffix "
    
    set usage "Usage: setup_gelclusters_largest "
    array set params [::cmdline::getoptions args $options $usage]

    if { [file exists "$outputdir/time_vs_gelclusters_largest$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_cluster [open "$outputdir/time_vs_gelclusters_largest$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
    	puts $f_cluster "\# Time size_cluster_largest"
    }
    close $f_cluster

}

proc ::cgtools::analysis::gelclusters_largest::analyze_gelclusters_largest {  } {
    #::mmsg::send [namespace current] "analyzing gelclusters_largest"
    variable ::cgtools::analysis::topology
    variable ::cgtools::system_specs
    variable size_gelclusters_largest 

    #========= get the number of the lipids of all the gel clusters
    set gelclusters_all [::cgtools::utils::compute_gelcluster_list $system_specs $topology]
    set gelclusters_largest [::cgtools::utils::find_the_largest_cluster $gelclusters_all] 
    set size_gelclusters_largest [llength $gelclusters_largest]

    #========= get the number of the lipids of all the gel clusters
    # set gelclusters_largest []

    ::mmsg::debug [namespace current] "done"
}

