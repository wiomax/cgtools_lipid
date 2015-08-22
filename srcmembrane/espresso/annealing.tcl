# ::cgtools::espresso
#
#------------------------------------------------------------#
# Routines for annealing process (Slow rate, for phase transformation)
# Author: Zun-Jing Wang 
# Apr 27 2010
#
# -----------------------------------------------------------#
namespace eval cgtools {
    namespace eval espresso {

        #################################################################################
        # proc annealing_init
        proc annealing_init { } {
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


                # Start annealing after warming
                setmd time_step $cgtools::main_time_step
		
                set nnnnnn $startn
                set tttttt $startt
                set kkkkkk $startk
               	set jjjjjj $startj
		#if {$kkkkkk == $cgtools::int_n_times} {
	 	#	set startk 0
		#	incr tttttt
		#}
	        puts "nnnnnn = $nnnnnn"
	        puts "tttttt = $tttttt"
	        puts "kkkkkk = $kkkkkk"
	        puts "jjjjjj = $jjjjjj"

		set n_annealingstart [expr int(($cgtools::systemtemp-$cgtools::annealing_templ)/$cgtools::annealing_deltaT)+1]
		set deltaT_annealingstart [expr ($cgtools::systemtemp-$cgtools::annealing_templ)/($n_annealingstart*1.0)]

		set n_annealing [expr int(($cgtools::annealing_temph-$cgtools::annealing_templ)/$cgtools::annealing_deltaT)+1]
		set deltaT_annealing [expr ($cgtools::annealing_temph-$cgtools::annealing_templ)/($n_annealing*1.0)]
	
		puts "$n_annealingstart, $deltaT_annealingstart, $n_annealing, $deltaT_annealing"
		while { $nnnnnn <= $cgtools::annealing_loops } {
		    if { $nnnnnn == 0 } { 
		    # if $nnnnnn==0,  cooling from $systemp to $annealingtempl 
			set ntotal_annealing [expr $n_annealingstart+1] 
		    }
		    if { $nnnnnn != 0 } { 
		    # if $nnnnnn!=0, first heating from $annealingtempl to $annealingtemph; 
		    #		 then cooling from $annealingtemph to $annealingtempl
			#set ntotal_annealing [expr 2*$n_annealing]
		    #### only heating, no cooling
			set ntotal_annealing [expr $n_annealing + 1]
		    }

		    puts "ntotal_annealing= $ntotal_annealing"
		    puts "tttttt = $tttttt"

		    while { $tttttt < $ntotal_annealing } {

			if { $nnnnnn == 0 } {
			    set temp_current [expr $cgtools::systemtemp-$deltaT_annealingstart*$tttttt]
			}

			if { $nnnnnn != 0 } {
		    	# first heating $annealingtempl-->$annealingtemph; then cooling $annealingtemph --> $annealingtempl
			    #if { $tttttt <  $n_annealing } {
				#set temp_current [expr $cgtools::annealing_templ+$deltaT_annealing*($tttttt+1)]
			    #}
			    #if { $tttttt >= $n_annealing } {
				#set temp_current [expr $cgtools::annealing_temph-$deltaT_annealing*($tttttt-$n_annealing+1)]
			    #}
			    
		       #### simulate at $cgtools::annealing_templ for a while, then heating $annealingtempl-->$annealingtemph 
			    set temp_current [expr $cgtools::annealing_templ+$deltaT_annealing*$tttttt]
			}

                	puts "tttttt = $tttttt"
	        	puts "temp_current = $temp_current"
                	puts "annealing  from  [setmd temp]"
			thermostat off
    			thermostat langevin $temp_current $cgtools::langevin_gamma
                	puts "annealing  to  [setmd temp]"
       	        	if { $cgtools::npt == "on" } {
               		        integrate set npt_isotropic $cgtools::p_ext $cgtools::piston_mass 1 1 0
	               		puts "npt integrator has been set"
        	                flush stdout
       	       		        #-cubic_box
        		        thermostat set npt_isotropic $temp_current  $cgtools::gamma_0  $cgtools::gamma_v
            	        }

	        	#puts "kkkkkk= $kkkkkk"
                	#puts "cgtools::int_n_times=$cgtools::int_n_times"

                	#Perform Annealing 
			#------------------------------------------------------------------------------
			for {set kkkkkk $startk } { $kkkkkk < $cgtools::int_n_times } { incr kkkkkk} {
                        	#puts "run $kkkkkk at time=[setmd time]"

	                        # Do the real work of integrating equations of motion
   	                        integrate $cgtools::int_steps

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
						"$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].crd" $topology \
						-periodbox 1 -computecomz 1

                                   		  ::cgtools::utils::writepdb_charmm \
						"$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].pdb" $topology \
						-periodbox 1 -computecomz 1
                                	}

                                	# Write a checkpoint to allow restarting.  Overwrites previous checkpoint
                                	::mmsg::send $this "setting checkpoint $kkkkkk [setmd time] $jjjjjj"

			                #if { [file exists "$cgtools::outputdir/checkpoint.latest.chk"] } {
					#	catch { exec rm $cgtools::outputdir/checkpoint.latest.chk }	
					#}

			                #if { [file exists "$cgtools::outputdir/checkpoint.latest.out"] } {
				        #         catch { exec mv $cgtools::outputdir/checkpoint.latest.out \ 
        	     			#            $cgtools::outputdir/checkpoint.latest.out.old }
					#}

        	                        checkpoint_set "$cgtools::outputdir/checkpoint.latest.out.$jjjjjj"

        	                        #checkpoint_set "$cgtools::outputdir/checkpoint.latest.out"

                                	#catch { exec cp $cgtools::outputdir/checkpoint.latest.out \ 
	                                #        $cgtools::outputdir/checkpoint.latest.out.$jjjjjj }

                	                # Try to copy a checkpoint to the backup checkpoint folder
                 	                # Usefull if the program crashes while writing a checkpoint
           	                        #if { [ catch { exec cp $cgtools::outputdir/checkpoint.latest.out \
	                                #        $cgtools::outputdir/checkpoint_bak/ } ] } {
        	                        #        ::mmsg::warn $this "warning: couldn't copy backup checkpoint"
                	                #}

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
