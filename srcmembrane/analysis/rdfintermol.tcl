# ::cgtools::analysis::analyze_rdfintermol --
#
# Calculate the rdfintermol of the system.  
# Author: Zun-Jing Wang 
# Dec. 2008

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::rdfintermol {
    variable rdfintermol_infolist
    variable rlist 
    variable avg_rdfintermollist_list
    variable rdfintermol_cnt
    variable rdfintermol_cntfile
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdfintermol
#newnew
    variable rdfintermol_output_infolist 
    variable rdfintermol_merge_coeflist 

    namespace export setup_rdfintermol
    namespace export analyze_rdfintermol
    namespace export printav_rdfintermol
    namespace export resetav_rdfintermol
}

proc ::cgtools::analysis::rdfintermol::resetav_rdfintermol { } {
    variable avg_rdfintermollist_list
    variable rdfintermol_cnt
    variable n_bin 
    variable n_rdfintermol

    set rdfintermol_cnt 0
    set avg_rdfintermollist_list [::cgtools::utils::init_matrix $n_rdfintermol $n_bin] 
}

proc ::cgtools::analysis::rdfintermol::printav_rdfintermol { } {
   global ::cgtools::analysis::outputdir
   variable rdfintermol_infolist
   variable rlist 
   variable avg_rdfintermollist_list
   variable rdfintermol_cnt
   variable rdfintermol_cntfile
   variable range_min 
   variable range_max 
   variable n_bin 
   variable n_rdfintermol
#newnew
   variable rdfintermol_output_infolist 
   variable rdfintermol_merge_coeflist 

   mmsg::send [namespace current] "printing rdfintermol at rdfintermol_cnt = $rdfintermol_cnt"

   if {  $rdfintermol_cnt > 0 } {

     # If <$outputdir/rdfcgintermol> doesn't exist then create it
     catch { exec mkdir $outputdir/rdfcgintermol }

#newnew
     set n_rdfintermoloutput [llength $rdfintermol_output_infolist]
     mmsg::send [namespace current] "n_rdfintermoloutput = $n_rdfintermoloutput"
     for { set i_rdfintermoloutput 0 } { $i_rdfintermoloutput  < $n_rdfintermoloutput } { incr i_rdfintermoloutput } {

        set current_rdfintermollist [lindex $rdfintermol_output_infolist $i_rdfintermoloutput]
        set nrdfintermol_plus [llength $current_rdfintermollist]

        set i_rdfintermol [lindex $current_rdfintermollist 0]

        set rdfintermol_info_now [lindex $rdfintermol_infolist $i_rdfintermol]
        set typei [lindex $rdfintermol_info_now 0]
        set typej [lindex $rdfintermol_info_now 1]
        set rdfintermolfile [lindex $rdfintermol_info_now 2]

        set nbin_rlist [llength $rlist]

        set avg_rdfintermollist [lindex $avg_rdfintermollist_list $i_rdfintermol]
        set rscale_rdfintermol [expr 1.0/($rdfintermol_cnt*1.0)]
        set avg_rdfintermollist [::cgtools::utils::scalevec $avg_rdfintermollist $rscale_rdfintermol]

        if { $nrdfintermol_plus > 1 } {
           #mmsg::send [namespace current] "merge rdfintermol $rdfintermolfile, # of files is  $nrdfintermol_plus"
           if { $nrdfintermol_plus == 2 } {
                set coeffnow [lindex $rdfintermol_merge_coeflist 0]
                if { [llength $coeffnow] != $nrdfintermol_plus } {
                       mmsg::err [namespace current] "coeffnow is wrong: $coeffnow "
                }
           }
           if { $nrdfintermol_plus == 3 } {
                set coeffnow [lindex $rdfintermol_merge_coeflist 1]
                if { [llength $coeffnow] != $nrdfintermol_plus } {
                       mmsg::err [namespace current] "coeffnow is wrong: $coeffnow "
                } 
           }
           if { $nrdfintermol_plus != 3 && $nrdfintermol_plus != 2 } {
                mmsg::err [namespace current] "nrdfintermol_plus is wrong: nrdfintermol_plus = $nrdfintermol_plus "
           }

           set f_rdfintermol [open "$rdfintermolfile.0" w]
           puts $f_rdfintermol "\# $nbin_rlist     [lindex $rlist 0]       [lindex $rlist [expr $nbin_rlist - 1] ]"
           foreach r $rlist rdfintermol $avg_rdfintermollist {
                puts $f_rdfintermol "$r $rdfintermol"
           }
           close $f_rdfintermol
            
           set avg_rdfintermollist [::cgtools::utils::scalevec $avg_rdfintermollist [lindex $coeffnow 0] ]


           for { set irdfintermol_index 1 } { $irdfintermol_index  < $nrdfintermol_plus } { incr irdfintermol_index} {

                set i_rdfintermol [lindex $current_rdfintermollist $irdfintermol_index]
                set avg_rdfintermolpluslist [lindex $avg_rdfintermollist_list $i_rdfintermol]

                set rscale_rdfintermol [expr 1.0/($rdfintermol_cnt*1.0)]
                set avg_rdfintermolpluslist [::cgtools::utils::scalevec $avg_rdfintermolpluslist $rscale_rdfintermol]

           	set f_rdfintermol [open "$rdfintermolfile.$irdfintermol_index" w]
           	puts $f_rdfintermol "\# $nbin_rlist     [lindex $rlist 0]       [lindex $rlist [expr $nbin_rlist - 1] ]"
           	foreach r $rlist rdfintermol $avg_rdfintermolpluslist {
                	puts $f_rdfintermol "$r $rdfintermol"
           	}
           	close $f_rdfintermol

                set coeffvalue [lindex $coeffnow $irdfintermol_index]
                set avg_rdfintermolpluslist [::cgtools::utils::scalevec $avg_rdfintermolpluslist $coeffvalue ]

                set avg_rdfintermollist [::cgtools::utils::add2vec $avg_rdfintermollist $avg_rdfintermolpluslist]
           }
        }

        set f_rdfintermol [open "$rdfintermolfile" w]
        puts $f_rdfintermol "\# $nbin_rlist     [lindex $rlist 0]       [lindex $rlist [expr $nbin_rlist - 1] ]"
        foreach r $rlist rdfintermol $avg_rdfintermollist {
                puts $f_rdfintermol "$r $rdfintermol"
        }
        close $f_rdfintermol
        #after 100
     } 
     #end for { set i_rdfintermoloutput 0 }

     # backup <$outputdir/rdfcgintermol> 
     incr rdfintermol_cntfile

     # checkwith if the filedir exists already 
     set newdir 0
     while { $newdir == 0 } { 
    	if { [file isdirectory  $outputdir/rdfcgintermol.$rdfintermol_cntfile ] } {
		incr rdfintermol_cntfile
    	} else {
		set newdir 1
	}
     }

     puts "rdfintermol_cntfile = $rdfintermol_cntfile"
     catch { exec cp -r $outputdir/rdfcgintermol  $outputdir/rdfcgintermol.$rdfintermol_cntfile }

   }
}

#newnew
proc ::cgtools::analysis::rdfintermol::setup_rdfintermol { rdfintermolcglist rdfintermolcgoutputlist args } {
    variable rdfintermol_infolist
    variable avg_rdfintermollist_list
    variable rdfintermol_cnt
    variable rdfintermol_cntfile
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdfintermol
#newnew
    variable rdfintermol_output_infolist 
    variable rdfintermol_merge_coeflist 
    variable ::cgtools::analysis::topology

    set options {
	{rmin.arg "0."  "minimum distange in rdfintermol computation"}
	{rmax.arg "15.01" "maximum distange in rdfintermol computation"}
	{nbin.arg "1501" "number of bins in rdfintermol computation "}
    }
    set usage "Usage: ::cgtools::analysis::rdfintermol::setup_rdfintermol rdfintermolcglist rdfintermolcgoutputlist \[rmin:rmax:nbin]"
    array set params [::cmdline::getoptions args $options $usage]

    set rdfintermol_infolist $rdfintermolcglist
    set n_rdfintermol [llength $rdfintermol_infolist]
#newnew
    set rdfintermol_output_infolist $rdfintermolcgoutputlist

    #set rdfintermol_merge_coeflist which is used in printav_rdf
    set nmols [llength $topology]
    mmsg::send [namespace current] "nmols = $nmols"

    set coeflist_now 0.5
    lappend coeflist_now 0.5
    lappend rdfintermol_merge_coeflist $coeflist_now

    set coeflist_now [expr 0.5*($nmols*1.0-1.0)/(2.*$nmols-1.) ]
    lappend coeflist_now [expr ($nmols * 1.0)/(2.*$nmols-1.) ]
    lappend coeflist_now [expr 0.5*($nmols*1.0-1.0)/(2.*$nmols-1.) ]
    lappend rdfintermol_merge_coeflist $coeflist_now
    mmsg::send [namespace current] "$rdfintermol_merge_coeflist"


    set rdfintermol_cnt 0
    set rdfintermol_cntfile 0
    set range_min $params(rmin)
    set range_max $params(rmax)
    set n_bin $params(nbin)

    set avg_rdfintermollist_list [::cgtools::utils::init_matrix $n_rdfintermol $n_bin] 
}

proc ::cgtools::analysis::rdfintermol::analyze_rdfintermol { } {
    variable rdfintermol_infolist
    variable rlist 
    variable avg_rdfintermollist_list
    variable rdfintermol_cnt
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdfintermol

    mmsg::send [namespace current] "analyzing rdfintermol"

    set rdfintermollist_list_now "" 
    for { set i_rdfintermol 0 } { $i_rdfintermol < $n_rdfintermol } { incr i_rdfintermol } {
	set rdfintermol_info_now [lindex $rdfintermol_infolist $i_rdfintermol]
	set typei [lindex $rdfintermol_info_now 0]
	set typej [lindex $rdfintermol_info_now 1]
	set rdfintermolfile [lindex $rdfintermol_info_now 2]

    	set rdfintermol [analyze rdf-intermol $typei $typej $range_min $range_max $n_bin]
    	set rlist ""
    	set rdfintermollist ""
    	foreach value [lindex $rdfintermol 1] {
   		lappend rlist [lindex $value 0]
# time factor 0.5 to the current rdfintermol, since the boxsize is 2 times big of the AA simulations
    		lappend rdfintermollist [expr [lindex $value 1]*0.5] 
	}
	lappend rdfintermollist_list_now $rdfintermollist
    }
    #mmsg::send [namespace current] "$rdfintermollist_list_now"
    set avg_rdfintermollist_list [::cgtools::utils::add_matrixs $avg_rdfintermollist_list $rdfintermollist_list_now] 

    incr rdfintermol_cnt

    mmsg::send [namespace current] "analyzing rdfintermol done"
}
