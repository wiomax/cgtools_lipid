# ::cgtools::espresso
#
#------------------------------------------------------------#
# Routines for annealfast process (for removing metastable status, and speed up equilibrium)
# Author: Zun-Jing Wang 
# Apr 27 2010
#------------------------------------------------------------#
#
namespace eval cgtools {
    namespace eval espresso {

        #################################################################################
        # proc annealfast_init
        proc annealfast_init { } {
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
		global topology 
    		global jjjjjj
    		global kkkkkk
    		global tttttt 
    		global nnnnnn 
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

    			# Set the box dimensions
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

    			# Initialize topology
    			set topology [::cgtools::generation::generate_system $cgtools::system_specs $cgtools::setbox_l]
    			# Set the generated topology into the internals of espresso.
    			::cgtools::utils::set_topology $topology
	    		
			# set exclustions for the bonded particles 1: 2-bodyinteraction  2: 3-body interaction
    			part auto_exclusions 1

    			#Initialise Random Number Generator
    			::cgtools::utils::init_random $cgtools::nprocessors

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
    			set startt 0
    			set startj 0
    			set startk 0

                	# Warm up to the systemtemp
	                setmd time_step $cgtools::warm_time_step
        	        thermostat langevin  $cgtools::systemtemp $cgtools::langevin_gamma
  	                ::mmsg::send $this "warming up without fixed particles at  [setmd temp]"
 			::cgtools::utils::warmup  $cgtools::warmsteps $cgtools::warmtimes $topology -cfgs 10 \
                                -outputdir $cgtools::outputdir
			#equlibrate after warming NVT at $cgtools::systemtemp
			integrate $cgtools::annealfast_int_n_times

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
    			set startt $tttttt 
    			set startj [expr $jjjjjj + 1]
    			set startk [expr $kkkkkk + 1]

		}

                # Start annealfast after warming
                setmd time_step $cgtools::main_time_step
		
                set nnnnnn $startn
                set tttttt $startt
                set kkkkkk $startk
               	set jjjjjj $startj
	        puts "nnnnnn = $nnnnnn"
	        puts "tttttt = $tttttt"
	        puts "kkkkkk = $kkkkkk"
	        puts "jjjjjj = $jjjjjj"


		set ntotal_annealfast 3
		while { $nnnnnn <= $cgtools::annealfast_loops } {
		    while { $tttttt < $ntotal_annealfast } {

                	puts "tttttt = $tttttt"
			puts "box_current = [setmd box_l]"
                	puts "annealfast  from  [setmd temp]"
			if { $tttttt == 0 } {
			     set temp_current $cgtools::annealfast_temph
               		     integrate set nvt 
	               	     puts "nvt integrator has been set"
        	             flush stdout
			     thermostat off
    			     thermostat langevin $temp_current $cgtools::langevin_gamma

		    	} else {
			   if { $tttttt == 1 } {
			     set temp_current $cgtools::systemtemp
                             integrate set nvt
                             puts "nvt integrator has been set"
                             flush stdout
                             thermostat off
                             thermostat langevin $temp_current $cgtools::langevin_gamma
			   } else {
			     set temp_current $cgtools::systemtemp
       	        	     if { $cgtools::npt == "on" } {
               		        integrate set npt_isotropic $cgtools::p_ext $cgtools::piston_mass 1 1 0
	               		puts "npt integrator has been set"
        	                flush stdout
       	       		        #-cubic_box
        		        thermostat set npt_isotropic $temp_current  $cgtools::gamma_0  $cgtools::gamma_v
            	             }
			  }
			}
	        	puts "temp_current = $temp_current"
                	puts "annealfast  to  [setmd temp]"

                	#Perform Annealing 
			#------------------------------------------------------------------------------
			for {set kkkkkk $startk } { $kkkkkk < $cgtools::int_n_times } { incr kkkkkk} {
                        	#puts "run $kkkkkk at time=[setmd time]"

	                        # Do the real work of integrating equations of motion
   	                        if { $tttttt == 0 } { 
					integrate $cgtools::annealfast_int_n_times
				} else {
					if { $tttttt == 1 } {
						integrate $cgtools::annealfast_int_n_times
					} else {
						integrate $cgtools::int_steps
					}
				}

                     		# Call all of the analyze routines that we specified when setting up our analysis
	                        ::cgtools::analysis::do_analysis

       		                # If kkkkkk is a multiple of analysis_write_frequency then write the analysis results to file
             		        if { [expr $kkkkkk + 1] % $cgtools::analysis_write_frequency ==0 } {
                        	        ::cgtools::analysis::print_averages
                    		}


               		        # If kkkkkk is a multiple of write_frequency then write out a full particle configuration
	                        if { [expr $kkkkkk + 1] % $cgtools::write_frequency ==0 } {

                                	if { $cgtools::use_vmd == "offline" } {
                                    		::cgtools::utils::writecrd_charmm \
                                                "$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].crd" $topology

                                   		  ::cgtools::utils::writepdb_charmm \
                                                "$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].pdb" $topology
                                	}

                                	# Write a checkpoint to allow restarting.  Overwrites previous checkpoint
                                	::mmsg::send $this "setting checkpoint $kkkkkk [setmd time] $jjjjjj"


        	                        checkpoint_set "$cgtools::outputdir/checkpoint.latest.out.$jjjjjj"
                                	incr jjjjjj
	                        }
       		                #end of if { [expr $kkkkkk + 0] % $cgtools::write_frequency ==0 }
             		}
            		#end of MD integration for {set kkkkkk $startk } { $kkkkkk <  $cgtools::int_n_times } { incr k}
			set startk 0
		        incr tttttt
   		   }
		   #end of tttttt 
		   set tttttt 0
		   incr nnnnnn 
		}
		#end of nnnnnn

		return
    	} 
	#end of proc espresso_init
     }
     #end of namespace eval espresso
}
#end of namespace eval cgtools 
