# ::cgtools::espresso
#
#------------------------------------------------------------#
# Routines for ffs process (Forward Flux Sampling)
# Author: Zun-Jing Wang 
# Aug 9 2010
# Has not been done yet
# -----------------------------------------------------------#
namespace eval cgtools {
    namespace eval espresso {

        #################################################################################
        # proc ffs_init
        proc ffs_init { } {
	    #############################################################
	    # local variable: 	
	    #		  	f_topo
	    #                 	starttime
	    #                 	i 
    	    #			startj
    	    #			startk
	    #                 	cutoff 
	    #                 	timingstart 
	    #                 	timingcurr
	    #                 	elapsedtime
	    #############################################################
    		#nnnnnn : number of the FFS trial runs 
    		global nnnnnn 
    		#ssssss : successful number of the FFS trial runs 
    		global ssssss 
    		#jjjjjj : number of integrated \tau  per FFS run
    		global jjjjjj
    		#kkkkkk : number of timesteps per \tau within a FFS run
    		global kkkkkk
		global topology 
    		global checkpointexists
		global errorInfo errorCode 

	     	set this [namespace current]
                ::mmsg::send $this "Feeding lipid parameters into ESPResSo running..."

    		# Attempt to read a checkpoint file
    		set checkpointexists [ ::cgtools::utils::readcheckpoint $cgtools::outputdir ]
 
		# ----------- Initialization ------------------ -----------#
		# Start a new computation, if no checkpoint
		if { !$checkpointexists } {
			# No checkpoint exists so we need to setup everything from scratch
    		        set startj 0
    			set startk 0

    			# Setup the output directory by creating it and copying forcetables and overlapped potcoffs to it
    			::cgtools::utils::setup_outputdir  $cgtools::outputdir -paramsfile $cgtools::paramsfile \
				-tabdir $cgtools::tabledir -tabnames $cgtools::tablenames -coffdir $cgtools::overlapdir \
				-coffnames $cgtools::overlapnames -readpdbdir $cgtools::readpdbdir \
				-readpdbname $cgtools::readpdbname

    			# Set the box dimensions : no need becuase the box_l is now set in ::cgtools::utils::readcrd
    			setmd box_l [lindex $cgtools::setbox_l 0] [lindex $cgtools::setbox_l 1] [lindex $cgtools::setbox_l 2]
                	setmd periodic 1 1 1
			if {$cgtools::linetension } {setmd periodic 1 0 1}
			puts "period is [setmd periodic]"
	
    			# Specify the bonded interactions
    			#puts "bonded_parms $cgtools::bonded_parms"
    			::cgtools::utils::set_bonded_interactions $cgtools::bonded_parms

    			# Specify any other non-bonded interactions
    			if { [ catch { ::cgtools::utils::set_nb_interactions $cgtools::nb_interactions } ] } {
				::mmsg::send $this "no non-bonded interactions used"
    			}

    			set cutoff [setmd max_cut] 
    			puts "max_cut is $cutoff"

    			# Initialize variable moltypelists in the namespace ::cgtools::utils
    			::cgtools::utils::initmoltypeskey $cgtools::moltypelists 

    			# Initialize topology and read crd files
    			set topology [::cgtools::generation::generate_system $cgtools::system_specs $cgtools::setbox_l]
    			# Set the generated topology into the internals of espresso.
    			::cgtools::utils::set_topology $topology
				
	    		
			# set exclustions for the bonded particles 1: 2-bodyinteraction  2: 3-body interaction
    			part auto_exclusions 1

    			#Initialise Random Number Generator
    			::cgtools::utils::init_random $cgtools::nprocessors

    			#Set velect list shell-size
			setmd skin      $cgtools::verlet_skin
    
    			#write topology file
    			set f_topo [open "$cgtools::outputdir/$cgtools::ident.top" w]
    			blockfile_write_topology $f_topo write topology   
    			close $f_topo

    			# Check if there are any extra vmdcommands and if not initialize a default
    			::cgtools::utils::initialize_vmd $cgtools::use_vmd $cgtools::outputdir $cgtools::ident \
				$topology -extracommands $cgtools::vmdcommands

                        # Setup analysis
                        ::cgtools::analysis::setup_analysis $cgtools::analysis_flags -outputdir  $cgtools::outputdir \
                                 -g $cgtools::mgrid -str $cgtools::stray_cut_off
			
    			# Make sure that we start exactly from where the checkpoint was written
    			set startn 0
    			set starts 0
    			set startj 0
    			set startk 0

                	# Warm up to the systemtemp
	                setmd time_step $cgtools::warm_time_step
        	        thermostat langevin  $cgtools::systemtemp $cgtools::langevin_gamma
  	                ::mmsg::send $this "warming up without fixed particles at  [setmd temp]"
 			::cgtools::utils::warmup  $cgtools::warmsteps $cgtools::warmtimes $topology -cfgs 10 \
                                -outputdir $cgtools::outputdir

			# Reset the time to a starttime (usually zero) after warmup
                        setmd time $cgtools::startmdtime

		}
 
		# Resume a computation, if exists checkpoint
		if { $checkpointexists } {
    			# A checkpoint exists so all we need to do is reset the moltypelists, topology and setup analysis again
    
    			::cgtools::utils::initmoltypeskey $cgtools::moltypelists 
    			::cgtools::utils::read_topology "$cgtools::outputdir/$cgtools::topofile"
			
			# set exclustions for the bonded particles 1: 2-bodyinteraction  2: 3-body interaction
			set topology [analyze set] 
    			::cgtools::utils::initialize_vmd $cgtools::use_vmd $cgtools::outputdir $cgtools::ident \
				$topology
			#puts "$topology"

    			part auto_exclusions 1
                        # Setup analysis
                        ::cgtools::analysis::setup_analysis $cgtools::analysis_flags -outputdir  $cgtools::outputdir \
                                 -g $cgtools::mgrid -str $cgtools::stray_cut_off
			

    			# Make sure that we start exactly from where the checkpoint was written
    			set startn $nnnnnn
    			set starts [expr $ssssss + 1]
    			set startj [expr $jjjjjj + 1]
    			set startk [expr $kkkkkk + 1]

		}


                # Start ffs parameters after warming
                setmd time_step $cgtools::main_time_step
		
                set nnnnnn $startn
                set ssssss $starts
                set kkkkkk $startk
               	set jjjjjj $startj
	        puts "nnnnnn = $nnnnnn"
	        puts "ssssss = $ssssss"
	        puts "kkkkkk = $kkkkkk"
	        puts "jjjjjj = $jjjjjj"

		# checkwith if the $cgtools::ffs_nextintersavefolder exists in $cgtools::outputdir already 
               	if { ![ file isdirectory $cgtools::outputdir/$cgtools::ffs_nextintersavefolder] } {
     			catch { exec mkdir $cgtools::outputdir/$cgtools::ffs_nextintersavefolder }
		}
		# checkwith if the $ssssss right
	        #set newfile 0
		#set interfacesave_cntfile 0
                #while { $newfile == 0} { 
                #	if { [ file exists $cgtools::outputdir/$cgtools::ffs_nextintersavefolder/checkpoint.$ssssss] } {
	        #      	 	incr interfacesave_cntfile
               	#	} else {
        	#        	set newfile 1
               	#	}
	        #}

		while { $nnnnnn < $cgtools::ffs_numbertrialloops_thisinterface} {

		    set go_out 0 

		    # Set Thermostat and Barastat 
		    set temp_current $cgtools::systemtemp
		    thermostat off
    		    thermostat langevin $temp_current $cgtools::langevin_gamma
               	    puts "ffs  to  [setmd temp]"
       	            if { $cgtools::npt == "on" } {
               	        integrate set npt_isotropic $cgtools::p_ext $cgtools::piston_mass 1 1 0
	             	puts "npt integrator has been set"
        	        flush stdout
       	       		#-cubic_box
        		thermostat set npt_isotropic $temp_current  $cgtools::gamma_0  $cgtools::gamma_v
            	    }

	       	    #puts "kkkkkk= $kkkkkk"
                    #puts "cgtools::int_n_times=$cgtools::int_n_times"

               	    #Perform one FFS run 
		    #------------------------------------------------------------------------------
		    for {set kkkkkk $startk } { $kkkkkk < $cgtools::int_n_times } { incr kkkkkk} {
		    	while { $go_out != 1 } {
                        	#puts "run $kkkkkk at time=[setmd time]"

	                        # Do the real work of integrating equations of motion
   	                        integrate $cgtools::int_steps

                     		# Call all of the analyze routines that we specified when setting up our analysis
	                        ::cgtools::analysis::do_analysis

       		                # If kkkkkk is a multiple of analysis_write_frequency then write the analysis results to file
             		        if { [expr $kkkkkk + 1] % $cgtools::analysis_write_frequency ==0 } {
                        	        ::cgtools::analysis::print_averages
					### ======== compute order parameter: orderparameter_now
					## get the list of all clusters
					set gelclusterlist_all [ ::cgtools::utils::compute_gelcluster_list $cgtools::system_specs $topology ]
					puts "gelclusterlist_all = $gelclusterlist_all"

					## get the laregest cluster : local order parameter
					set gelcluster_largest [ ::cgtools::utils::find_the_largest_cluster $gelclusterlist_all ]
					## get all the gel molecules : global order parameter
					set gelmolall_list [ ::cgtools::utils::find_all_gel_molecules $gelclusterlist_all ]
					puts "gelmolall_list : $gelmolall_list"
					set composition_gel [ ::cgtools::utils::compute_composition $cgtools::system_specs $gelmolall_list]
					puts "composition_gel: $composition_gel"

				        set liqmolall_list [ ::cgtools::utils::find_all_liquid_molecules $gelclusterlist_all $topology]
					puts "liqmolall_list : $liqmolall_list"

					set composition_liq [ ::cgtools::utils::compute_composition $cgtools::system_specs $liqmolall_list]
					puts "composition_liq : $composition_liq"

				        set orderparameter_now [llength $gelcluster_largest]
					puts "orderparameter_now = $orderparameter_now"
    					::cgtools::utils::writepdb_cluster_charmm "$cgtools::outputdir/$cgtools::ident.cluster.vmd[format %04d $jjjjjj].pdb" $topology $gelclusterlist_all $gelcluster_largest -periodbox 1 -computecomz 1
                                        if { $orderparameter_now >= $cgtools::ffs_orderparameter_thisinterface } {	
					    ### ======== save the current coordinates and velocities to the inter-int-nextint#
					    # set cgtools::ffs_nextintersavefolder
					    # count nsavenextinter in the folder at the beginning
					    checkpoint_set "$cgtools::outputdir/$cgtools::ffs_nextintersavefolder/checkpoint.$ssssss.OP_$orderparameter_now"

					    incr ssssss
	
					    set go_out 1
					}
                    		} 
				#end of if { [expr $kkkkkk + 1] % $cgtools::analysis_write_frequency ==0 }

               		        # If kkkkkk is a multiple of write_frequency then write out a full particle configuration
	                        if { [expr $kkkkkk + 1] % $cgtools::write_frequency ==0 } {

                                	if { $cgtools::use_vmd == "offline" } {
                                    		::cgtools::utils::writecrd_charmm \
						"$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].crd" $topology \
						-periodbox 1 -computecomz 1

                                   		  ::cgtools::utils::writepdb_charmm \
						"$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].pdb" $topology \
						-periodbox 1 -computecomz 1
                                	}

                                	# Write a checkpoint to allow restarting.  Overwrites previous checkpoint
                                	::mmsg::send $this "setting checkpoint $kkkkkk [setmd time] $jjjjjj"
        	                        #checkpoint_set "$cgtools::outputdir/checkpoint.latest.out.$jjjjjj"

			                if { [file exists "$cgtools::outputdir/checkpoint.latest.chk"] } {
						catch { exec rm $cgtools::outputdir/checkpoint.latest.chk }	
					}
			                if { [file exists "$cgtools::outputdir/checkpoint.latest.out"] } {
				                 catch { exec mv $cgtools::outputdir/checkpoint.latest.out \ 
        	     			            $cgtools::outputdir/checkpoint.latest.out.old }
					}
					checkpoint_set "$cgtools::outputdir/checkpoint.latest.out"

                                	incr jjjjjj
	                        }
       		                #end of if { [expr $kkkkkk + 0] % $cgtools::write_frequency ==0 }
   		   	}
		    	#while { $go_out != 1 }
             	    }
            	    #end of MD integration for {set kkkkkk $startk } { $kkkkkk <  $cgtools::int_n_times } { incr k}
		    set startk 0
	            incr nnnnnn 

    		    # ======= Initialize a new FFS run 
	    	    # READ INITIAL COORDINATES 
	    	    # RESET the INITIAL box dimensions : no need becuase the box_l is now set in ::cgtools::utils::readcrd
	    	    # setmd box_l [lindex $cgtools::setbox_l 0] [lindex $cgtools::setbox_l 1] [lindex $cgtools::setbox_l 2]
		    # Initialise Random Number Generator
		    ::cgtools::utils::init_random $cgtools::nprocessors
		}
		#end of nnnnnn

		return
    	} 
	#end of proc espresso_init
     }
     #end of namespace eval espresso
}
#end of namespace eval cgtools 
