# ::cgtools::analysis::analyze_gelclusters_all --
#
# Extract the box dimensions from espresso
# Author: Zun-Jing Wang
# 2009 August

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::gelclusters_all {
    variable gelclusters_all 
    variable f_cluster

    namespace export setup_gelclusters_all
    namespace export analyze_gelclusters_all
    namespace export printav_gelclusters_all
    namespace export resetav_gelclusters_all
}

proc ::cgtools::analysis::gelclusters_all::resetav_gelclusters_all { } {
    variable gelclusters_all 
}

proc ::cgtools::analysis::gelclusters_all::printav_gelclusters_all { } {
    variable gelclusters_all
    global  ::cgtools::analysis::outputdir
    global ::cgtools::analysis::suffix
    variable f_cluster
    global ::cgtools::analysis::time
    
    set f_cluster [open "$outputdir/time_vs_gelclusters_all$suffix" a]
    set ncluster [llength $gelclusters_all]
    puts $f_cluster "$time $ncluster" nonewline
    for { set b 0 } { $b < $ncluster } {incr b } {
   	puts $f_cluster "[lindex $size_gelclusters $b] " nonewline 
    }
    puts $f_cluster ""
    close $f_cluster
    #flush $f_cluster
}

proc ::cgtools::analysis::gelclusters_all::setup_gelclusters_all { args } {
    global  ::cgtools::analysis::outputdir
    variable f_cluster
    global ::cgtools::analysis::iotype
    global ::cgtools::analysis::suffix
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_gelclusters_all$suffix "
    
    set usage "Usage: setup_gelclusters_all "
    array set params [::cmdline::getoptions args $options $usage]

    if { [file exists "$outputdir/time_vs_gelclusters_all$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_cluster [open "$outputdir/time_vs_gelclusters_all$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
    	puts $f_cluster "\# Time number_of_gel_clusters size_cluster_1 size_cluster_2 size_cluster_3 ...."
    }
    close $f_cluster

}

proc ::cgtools::analysis::gelclusters_all::analyze_gelclusters_all {  } {
    #::mmsg::send [namespace current] "analyzing gelclusters_all"
    variable ::cgtools::analysis::topology
    variable ::cgtools::system_specs
    variable gelclusters_all

    #========= get the number of the lipids of all the gel clusters
    set gelclusters_all [::cgtools::utils::compute_gelcluster_list $system_specs $topology]

    ::mmsg::debug [namespace current] "done"
}

