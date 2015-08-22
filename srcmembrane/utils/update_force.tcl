# Update forcetables.  
# Author: Zun-Jing Wang 
# Nov. 17 2008
#

namespace eval ::cgtools::utils {
    namespace export update_force 
    namespace export readtablines
}

# ::cgtools::utils::update_force --
#
# This routine is designed to update all the force tables for the optimization of the force field
# by matching the rdf between AA and CG simulations. 
#
proc ::cgtools::utils::update_force { rdfcglist rdfaalist forcetabdir forcetabnamelist args } {
    set options {
	{rmin.arg "0."  "minimum distange in rdf computation"}
	{rmax.arg "15.01" "maximum distange in rdf computation"}
	{nbin.arg "1501" "number of bins in rdf computation "}
    }
    set usage "Usage: ::cgtools::utils::update_force rdfcglist rdfaalist forcetabdir forcetabnamelist \[rmin:rmax:nbin]"
    array set params [::cmdline::getoptions args $options $usage]

    mmsg::send [namespace current] "updating force tables"

    #set range_min $params(rmin)
    #set range_max $params(rmax)
    #set n_bin $params(nbin)

    set errcode [ catch { exec mkdir $forcetabdir/old } ]
    foreach filename $forcetabnamelist {
    	set errcode [ catch { exec cp $forcetabdir/$filename  $forcetabdir/old/ } ]
    	if { $errcode } {
        	::mmsg::warn [namespace current]  "couldn't cp $filename from $forcetabdir to $forcetabdir/old/"
    	} 
    }
	

    set n_rdf [llength $rdfcglist]
    for { set i_rdf 0 } { $i_rdf < $n_rdf } { incr i_rdf } {
        #read the rdf_aa 
	set rdfaa_info_now [lindex $rdfaalist $i_rdf]
	set rdffile_aa [lindex $rdfaa_info_now 2]
	set rdfaa_lines_list [::cgtools::utils::readtablines $rdffile_aa]
	set nbin_rdfaa [llength $rdfaa_lines_list]
	set rdfaa_list ""
	foreach iline $rdfaa_lines_list {
		lappend rdfaa_list [lindex $iline 1]
	}

        #read the rdf_cg 
	set rdfcg_info_now [lindex $rdfcglist $i_rdf]
	set rdffile_cg [lindex $rdfcg_info_now 2]
	set rdfcg_lines_list [::cgtools::utils::readtablines $rdffile_cg]
	set nbin_rdfcg [llength $rdfcg_lines_list]
    	if { $nbin_rdfcg !=  $nbin_rdfaa } {
       		::mmsg::warn [namespace current]  "the rdfcg doesn't match rdfaa"
    	} 
	set rdfcg_list ""
	foreach iline $rdfcg_lines_list {
		lappend rdfcg_list [lindex $iline 1]
	}

        #read the force_old 
	set tabname [lindex $forcetabnamelist $i_rdf]
	set tabfile_old "$forcetabdir/old/$tabname"
	set forceold_lines_list [::cgtools::utils::readtablines $tabfile_old]
	set nbin_forceold [llength $forceold_lines_list]
    	if { $nbin_forceold !=  [expr $nbin_rdfaa -1] } {
       		::mmsg::warn [namespace current]  "the foretable doesn't match rdfaa"
    	} 
	set potold_list ""
	set forceold_list ""
	set r_list ""
	foreach iline $forceold_lines_list {
		lappend r_list [lindex $iline 0]
		lappend forceold_list [lindex $iline 1]
		lappend potold_list [lindex $iline 2]
	}
	lappend potold_list 0.

        #compute the pot_new
    	set EACC 1.0E-4 
        set alpha 0.15
	set potnew_templist ""
	foreach rdfcg $rdfcg_list rdfaa $rdfaa_list potold $potold_list {
		if { $rdfcg < $EACC } { set rdfcg $EACC}
		if { $rdfaa < $EACC } { set rdfaa $EACC}
		set potnew_temp [expr $potold + $alpha * log( $rdfcg / $rdfaa )]
		lappend potnew_templist $potnew_temp
	}
	
	set potnew_temp_cutoff [lindex $potnew_templist [expr $nbin_rdfaa - 1]]
	set potnew_list ""
	foreach potnew_temp $potnew_templist {
		set potnew [expr $potnew_temp - $potnew_temp_cutoff]
		lappend potnew_list $potnew
	}

        #compute the force_new
	set dr [expr [lindex $r_list 1] - [lindex $r_list 0]]
	set forcenew_list ""
        for { set ibin 0 } { $ibin < $nbin_forceold } { incr ibin } {
		set ibin_idecrhalf $ibin 
		set ibin_iplushalf [expr $ibin + 1]
		set potnew_idecrhalf [lindex $potnew_list $ibin_idecrhalf]
		set potnew_iplushalf [lindex $potnew_list $ibin_iplushalf] 
		set forcenew [expr ($potnew_idecrhalf - $potnew_iplushalf)/$dr]
		lappend forcenew_list $forcenew
	}

        #delete the last element in potnew_list
	set ibin_end [expr [llength $potnew_list] - 1]
	set ibin_endminus1 [expr [llength $potnew_list] - 2]
	set potnew_list [lreplace $potnew_list $ibin_endminus1 $ibin_end [lindex $potnew_list $ibin_endminus1]] 

        #print the force_new 
	set tabfile_new "$forcetabdir/$tabname"
    	set nl_rlist [llength $r_list]

        set f_tabnew [open "$tabfile_new" w]
        puts $f_tabnew "\#	$nl_rlist	[lindex $r_list 0]	[lindex $r_list [expr $nl_rlist-1] ]"
	foreach r $r_list force $forcenew_list pot $potnew_list {
		puts $f_tabnew [format "%10.5f %20.16e %20.16e" $r $force $pot] 
	}
	close $f_tabnew
    }

    mmsg::send [namespace current] "done"
}

proc ::cgtools::utils::readtablines { readfile } {

    set readlist 0
    unset readlist

    # open readfile
    set errcode [ catch { set readchannel [open $readfile "r"] } ]
    if { $errcode } {
        ::mmsg::err [namespace current]  "couldn't open file $readfile"
    }
 
    ::mmsg::send [namespace current] "reading tabs from the file $readfile "

    # read lines in pdb file 
    while {[gets $readchannel fileline] >= 0} {
        if { [lindex $fileline 0] != "#" } {
             lappend readlist $fileline
        }
    }
    close $readchannel

    return $readlist
}
