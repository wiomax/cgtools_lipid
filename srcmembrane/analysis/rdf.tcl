# ::cgtools::analysis::analyze_rdf --
#
# Calculate the rdf of the system.  
# Author: Zun-Jing Wang 
# Dec. 2008

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::rdf {
    variable rdf_infolist
    variable rlist 
    variable avg_rdflist_list
    variable rdf_cnt
    variable rdf_cntfile
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdf
    variable rdfoutput_infolist 
    variable rdfmerge_coeflist 

    namespace export setup_rdf
    namespace export analyze_rdf
    namespace export printav_rdf
    namespace export resetav_rdf
}

proc ::cgtools::analysis::rdf::resetav_rdf { } {
    variable avg_rdflist_list
    variable rdf_cnt
    variable n_bin 
    variable n_rdf

    set rdf_cnt 0
    set avg_rdflist_list [::cgtools::utils::init_matrix $n_rdf $n_bin] 
}

proc ::cgtools::analysis::rdf::printav_rdf { } {
   global ::cgtools::analysis::outputdir
   variable rdf_infolist
   variable rlist 
   variable avg_rdflist_list
   variable rdf_cnt
   variable rdf_cntfile
   variable range_min 
   variable range_max 
   variable n_bin 
   variable n_rdf
   variable rdfoutput_infolist 
   variable rdfmerge_coeflist 

   mmsg::send [namespace current] "printing rdf at rdf_cnt = $rdf_cnt"

   if {  $rdf_cnt > 0 } {

     # If <$outputdir/rdfcg> doesn't exist then create it
     catch { exec mkdir $outputdir/rdfcg }

     set n_rdfoutput [llength $rdfoutput_infolist]
     #mmsg::send [namespace current] "n_rdfoutput = $n_rdfoutput"
     for { set i_rdfoutput 0 } { $i_rdfoutput  < $n_rdfoutput } { incr i_rdfoutput } {
        
        set current_rdflist [lindex $rdfoutput_infolist $i_rdfoutput] 
        set nrdf_plus [llength $current_rdflist]     
        
	set i_rdf [lindex $current_rdflist 0]     

        set rdf_info_now [lindex $rdf_infolist $i_rdf]
        set typei [lindex $rdf_info_now 0]
        set typej [lindex $rdf_info_now 1]
        set rdffile [lindex $rdf_info_now 2]

	set nbin_rlist [llength $rlist]

	set avg_rdflist [lindex $avg_rdflist_list $i_rdf]
	set rscale_rdf [expr 1.0/($rdf_cnt*1.0)]
	set avg_rdflist [::cgtools::utils::scalevec $avg_rdflist $rscale_rdf]

        if { $nrdf_plus > 1 } {
           #mmsg::send [namespace current] "merge rdf $rdffile, # of files is  $nrdf_plus"
           if { $nrdf_plus == 2 } {
                set coeffnow [lindex $rdfmerge_coeflist 0]
                if { [llength $coeffnow] != $nrdf_plus } {
		       mmsg::err [namespace current] "coeffnow is wrong: $coeffnow "
		} 
           }
           if { $nrdf_plus == 3 } {
                set coeffnow [lindex $rdfmerge_coeflist 1]
                if { [llength $coeffnow] != $nrdf_plus } {
		       mmsg::err [namespace current] "coeffnow is wrong: $coeffnow "
		} 
           }
           if { $nrdf_plus != 3 && $nrdf_plus != 2 } {
		mmsg::err [namespace current] "nrdf_plus is wrong: nrdf_plus = $nrdf_plus "
           }
           
           set f_rdf [open "$rdffile.0" w]
           puts $f_rdf "\#	$nbin_rlist	[lindex $rlist 0]	[lindex $rlist [expr $nbin_rlist - 1] ]"
           foreach r $rlist rdf $avg_rdflist {
		puts $f_rdf "$r $rdf" 
	   }
	   close $f_rdf

	   set avg_rdflist [::cgtools::utils::scalevec $avg_rdflist [lindex $coeffnow 0] ]
           #mmsg::send [namespace current] "avg_rdflist: $avg_rdflist"

           for { set irdf_index 1 } { $irdf_index  < $nrdf_plus } { incr irdf_index} {
	        
		set i_rdf [lindex $current_rdflist $irdf_index]
		set avg_rdfpluslist [lindex $avg_rdflist_list $i_rdf]
                #mmsg::send [namespace current] "avg_rdfpluslist: $avg_rdfpluslist"

		set rscale_rdf [expr 1.0/($rdf_cnt*1.0)]
		set avg_rdfpluslist [::cgtools::utils::scalevec $avg_rdfpluslist $rscale_rdf]

                set f_rdf [open "$rdffile.$irdf_index" w]
                puts $f_rdf "\#	$nbin_rlist	[lindex $rlist 0]	[lindex $rlist [expr $nbin_rlist - 1] ]"
                foreach r $rlist rdf $avg_rdfpluslist {
		     puts $f_rdf "$r $rdf" 
	        }
	        close $f_rdf


                set coeffvalue [lindex $coeffnow  $irdf_index]
		set avg_rdfpluslist [::cgtools::utils::scalevec $avg_rdfpluslist $coeffvalue ]
                #mmsg::send [namespace current] "avg_rdfpluslist: $avg_rdfpluslist"

                #mmsg::send [namespace current] "[llength $avg_rdflist]"
                #mmsg::send [namespace current] "[llength $avg_rdfpluslist]"
		set avg_rdflist [::cgtools::utils::add2vec $avg_rdflist $avg_rdfpluslist]
                #mmsg::send [namespace current] "avg_rdflist: $avg_rdflist"
	   }
        }
	
        set f_rdf [open "$rdffile" w]
        puts $f_rdf "\#	$nbin_rlist	[lindex $rlist 0]	[lindex $rlist [expr $nbin_rlist - 1] ]"
        foreach r $rlist rdf $avg_rdflist {
		puts $f_rdf "$r $rdf" 
     	}

        close $f_rdf
        #after 100
     }

     # backup <$outputdir/rdfcg> 
     incr rdf_cntfile

     # checkwith if the filedir exists already 
     set newdir 0
     while { $newdir == 0} { 
    	if { [ file isdirectory $outputdir/rdfcg.$rdf_cntfile ] } {
		incr rdf_cntfile
    	} else {
		set newdir 1
	}
     }

     puts "rdf_cntfile = $rdf_cntfile"
     catch { exec cp -r $outputdir/rdfcg  $outputdir/rdfcg.$rdf_cntfile }

   }
}

proc ::cgtools::analysis::rdf::setup_rdf { rdfcglist rdfcgoutputlist args } {
    variable rdf_infolist
    variable avg_rdflist_list
    variable rdf_cnt
    variable rdf_cntfile
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdf
    variable rdfoutput_infolist 
    variable rdfmerge_coeflist 
    variable ::cgtools::analysis::topology

    set options {
	{rmin.arg "0."  "minimum distange in rdf computation"}
	{rmax.arg "15.01" "maximum distange in rdf computation"}
	{nbin.arg "1501" "number of bins in rdf computation "}
    }
    set usage "Usage: ::cgtools::analysis::rdf::setup_rdf rdfcglist rdfcgoutputlist \[rmin:rmax:nbin]"
    array set params [::cmdline::getoptions args $options $usage]

    set rdf_infolist $rdfcglist
    set n_rdf [llength $rdf_infolist]

    set rdfoutput_infolist $rdfcgoutputlist
    
    #set rdfmerge_coeflist which is used in printav_rdf
    set nmols [llength $topology]
    #mmsg::send [namespace current] "nmols = $nmols"

    set coeflist_now 0.5
    lappend coeflist_now 0.5
    lappend rdfmerge_coeflist $coeflist_now 

    set coeflist_now [expr 0.5*($nmols*1.0-1.0)/(2.*$nmols-1.) ]
    lappend coeflist_now [expr ($nmols * 1.0)/(2.*$nmols-1.) ]
    lappend coeflist_now [expr 0.5*($nmols*1.0-1.0)/(2.*$nmols-1.) ]
    lappend rdfmerge_coeflist $coeflist_now 
    #mmsg::send [namespace current] "$rdfmerge_coeflist"

    set rdf_cnt 0
    set rdf_cntfile 0
    set range_min $params(rmin)
    set range_max $params(rmax)
    set n_bin $params(nbin)

    set avg_rdflist_list [::cgtools::utils::init_matrix $n_rdf $n_bin] 
}

proc ::cgtools::analysis::rdf::analyze_rdf { } {
    variable rdf_infolist
    variable rlist 
    variable avg_rdflist_list
    variable rdf_cnt
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdf

    mmsg::send [namespace current] "analyzing rdf"

    set rdflist_list_now "" 
    for { set i_rdf 0 } { $i_rdf < $n_rdf } { incr i_rdf } {
	set rdf_info_now [lindex $rdf_infolist $i_rdf]
	set typei [lindex $rdf_info_now 0]
	set typej [lindex $rdf_info_now 1]
	set rdffile [lindex $rdf_info_now 2]

    	set rdf [analyze rdf $typei $typej $range_min $range_max $n_bin]
    	set rlist ""
    	set rdflist ""
    	foreach value [lindex $rdf 1] {
   		lappend rlist [lindex $value 0]
# time factor 0.5 to the current rdf, since the boxsize is 2 times big of the AA simulations
    		lappend rdflist [expr [lindex $value 1]*0.5] 
	}
	lappend rdflist_list_now $rdflist
    }
    set avg_rdflist_list [::cgtools::utils::add_matrixs $avg_rdflist_list $rdflist_list_now] 

    incr rdf_cnt

    mmsg::send [namespace current] "analyzing rdf done"
}
