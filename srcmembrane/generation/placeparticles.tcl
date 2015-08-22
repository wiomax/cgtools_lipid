# Routines for reading positions of particles 
#
# Author: Zun-Jing Wang
# Sep. 2008 Done 

namespace eval cgtools::generation {}

#::cgtools::generation::placeparticles_all
# Reading positions of all particles from a pdb file 
proc ::cgtools::generation::placeparticles_all { readfile } {

   variable topology
 
   # readpdbfile 
   #set linelist [::cgtools::utils::readpdb $readfile]
   set linelist [::cgtools::utils::readcrd $readfile]

   set atomNumber 0
   set molNumber 0
   set mollists 0
   unset mollists 

   foreach mol $topology {
   	set moltype [lindex $mol 0]
    	set typeinfo [::cgtools::utils::matchtype $moltype ]
  
 	set partbondlists [lindex $typeinfo 2]
    	set beadlists [lindex $partbondlists 0]
	set nbeads_mol [llength $beadlists]
    	
	set partbondtypelists [lindex $typeinfo 3]
    	set beadtypelists [lindex $partbondtypelists 0]

	set partlists 0
	unset partlists
    	for { set i 0 } { $i < $nbeads_mol } {incr i } {
		set curline [lindex $linelist $atomNumber] 
		#puts "curline= $curline"
	        set posx [lindex $curline 5] 
	        set posy [lindex $curline 6] 
	        set posz [lindex $curline 7]

	        set beadname_read [lindex $curline 2]
	        regsub ES $beadname_read E beadname_read 
	        set beadname [string range $beadname_read 0 1]

		set beadtype [lindex $beadlists $i]
		set beadname_right [lindex [lindex $beadtypelists $beadtype] 1]
		#puts "beadtype= $beadtype"
		#puts "beadtypelists= $beadtypelists"
		if { $beadname != $beadname_right } {
       			::mmsg::warn [namespace current]  "sequence of the particles is wrong in $readfile: line $atomNumber"
			#puts "beadname= $beadname"
			#puts "beadname_right= $beadname_right"
			exit

		}
		set posvec [list $posx $posy $posz]
		set curpart [list $atomNumber $posvec $beadtype $beadname $moltype $molNumber] 
		lappend partlists $curpart 
	
		incr atomNumber
		
	}
	lappend mollists $partlists 

	incr molNumber 
   }

   # Check particle consistency
   if { $atomNumber != [expr [::cgtools::utils::maxpartid $topology] + 1] } {
       	::mmsg::err [namespace current] "$readfile has $atomNumber particles but [expr [::cgtools::utils::maxpartid $topology] +1] were specified in topology "
   }

   return $mollists
}

#::cgtools::generation::placeparticles_template
# Reading positions of all particles from a pdb file
# Now it only works for one type of mol 
proc ::cgtools::generation::placeparticles_template { readfile } {

   variable topology
 
   # readpdbfile 
   #set linelist [::cgtools::utils::readpdb $readfile]
   set linelist [::cgtools::utils::readcrd $readfile]

   set atomNumber 0
   set molNumber 0
   set mollists 0
   unset mollists 

   set mol [lindex $topology 0]
   set moltype [lindex $mol 0]
   set typeinfo [::cgtools::utils::matchtype $moltype]
   
   set partbondlists [lindex $typeinfo 2]
   set beadlists [lindex $partbondlists 0]
   set nbeads_mol [llength $beadlists]
    	
   set partbondtypelists [lindex $typeinfo 3]
   set beadtypelists [lindex $partbondtypelists 0]

   set partlists 0
   unset partlists
   for { set i 0 } { $i < $nbeads_mol } {incr i } {
	set curline [lindex $linelist $atomNumber] 
        set posx [lindex $curline 5] 
        set posy [lindex $curline 6] 
        set posz [lindex $curline 7]

        set beadname_read [lindex $curline 2]
        regsub ES $beadname_read E beadname_read 
        set beadname [string range $beadname_read 0 1]

	set beadtype [lindex $beadlists $i]
	set beadname_right [lindex [lindex $beadtypelists $beadtype] 1]

	if { $beadname != $beadname_right } {
       		::mmsg::warn [namespace current]  "sequence of the particles is wrong in $readfile: line $atomNumber"
	}
	set posvec [list $posx $posy $posz]
	set curpart [list $atomNumber $posvec $beadtype $beadname $moltype $molNumber] 
	lappend partlists $curpart 
	
	incr atomNumber
   }
   lappend mollists $partlists 
   incr molNumber 

   # Check particle consistency
   if { $atomNumber != [expr [llength $mol] -1] } {
       	mmsg::warn [namespace current] "$readfile has $atomNumber particles but [expr [llength $mol] -1] were specified in the temple mol "
   }

   return $mollists
}


#::cgtools::generation::placeparticles_vectors
# place positions of all particles from a vectorlist
proc ::cgtools::generation::placeparticles_vectors { vectorlist } {

   variable topology
 
   set atomNumber [setmd n_part]
   set molNumber 0
   set mollists 0
   unset mollists 

   #puts "topology= $topology"

   set allvectors [lindex $vectorlist 0]
   foreach mol $topology {
   	set moltype [lindex $mol 0]
    	set typeinfo [::cgtools::utils::matchtype $moltype ]
	#puts "::cgtools::generation::placeparticles_vectors:moltype = $moltype"
	#puts "::cgtools::generation::placeparticles_vectors:typeinfo = $typeinfo"
  
 	set partbondlists [lindex $typeinfo 2]
    	set beadlists [lindex $partbondlists 0]
	set nbeads_mol [llength $beadlists]
	#puts "partbondlists= $partbondlists"
	#puts "beadlists= $beadlists"
	#puts "nbeads_mol= $nbeads_mol"
    	
	set partbondtypelists [lindex $typeinfo 3]
    	set beadtypelists [lindex $partbondtypelists 0]
	#puts "partbondtypelists = $partbondtypelists"
	#puts "beadtypelists= $beadtypelists"

	set posmollist [lindex $allvectors $molNumber]
	#puts "posmollist= $posmollist"
	set partlists 0
	unset partlists
    	for { set i 0 } { $i < $nbeads_mol } {incr i } {
	        set posvec [lindex $posmollist $i]
	        set beadname [lindex [lindex $beadtypelists $i] 1]
	        set beadtype [lindex [lindex $beadtypelists $i] 0]
		set curpart [list $atomNumber $posvec $beadtype $beadname $moltype $molNumber] 
		#puts "curpart= $curpart"
		lappend partlists $curpart 
	
		incr atomNumber
		
	}
	lappend mollists $partlists 

	incr molNumber 
   }
   #puts "mollists= $mollists"

   # Check particle consistency
   if { $atomNumber != [expr [::cgtools::utils::maxpartid $topology] + 1] } {
       	::mmsg::err [namespace current] "$vectorlist has $atomNumber particles but [expr [::cgtools::utils::maxpartid $topology] +1] were specified in topology "
   }

   return $mollists
}
