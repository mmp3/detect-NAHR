#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <algorithm>

#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include <utility>
#include <string.h>
#include <vector>
#include <map>
#include <list>
#include <set>

#include <omp.h>


#include <boost/numeric/interval.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>


// #define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>


#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>



#include <mpreal.h>
#include <mpfr.h>



#include <boost/random.hpp>

#include <boost/mpi/collectives.hpp>


#include <api/BamAlignment.h>
#include <api/BamReader.h>



#include <general_typedefs.h>
#include <SD_Entry_general.h>
#include <Interlocking_entry.h>
#include <templates.h>
#include <Coords_and_homologous_region.h>
#include <translations_through_alignment.h>
#include <Sparse_map.h>
#include <Readgroup_statistics.h>



#include "Conn_Comp.h"
#include "globals.h"
#include "Event.h"
#include "Partial_sum_function.h"
#include "other_functions.h"
#include "io_functions.h"
#include "Visitation_object.h"
#include "Paired_end_read.h"

#include "MiniRegion.h"
#include "Star_alignment.h"



// #include <unordered_map>;
// #include <boost/tr1/unordered_map.hpp>
#include <boost/graph/graph_concepts.hpp>














//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////











int Conn_Comp::perform_complete_variable_elimination_schedule
			(BamTools::BamReader &my_BAM_reader,
			 const bool &consider_GeneConversion_breakpoints_in_Heuristic,
			 const bool &consider_GeneConversion_in_Full_compute)
{      
    print_line_of_markers("#");
    
    std::cerr << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nrank " << my_MPI_rank
		<< ":  perform_complete_variable_elimination_schedule for CC_ID = " << Conn_Comp_ID 
		<< " !!!\n\t\t\t\t size = " << events.size() << "\n\n\n";
    
		
		
    //Heuristic:
    double time_begin = omp_get_wtime();
                
    type_set_uint excluded_Events;
    const bool success_heuristic = HEURISTIC__test_each_specific_variational_position_possibility(
						excluded_Events, consider_GeneConversion_breakpoints_in_Heuristic, my_BAM_reader);
    
    std::cerr<< "time Heuristic total = "  << omp_get_wtime() - time_begin << "\n";
    
    if (!success_heuristic)
    {
        if (my_MPI_rank == 0)
	{  append_to_skipped_connected_component_file(*this, "Heuristic on variational positions failed - during uploading PERs for some event");  }
        return -2;
    }
    
    
    
    
    
    if (my_MPI_rank == 0)
    {
	save_excluded_Events_to_file(excluded_Events);            
	save_remaining_original_breakpoints_to_file(events);
    }


    
    
    
    //Perform schedule...        
    for (type_list_Visitation_object::iterator it_Visit = variable_elimination_schedule.begin();
         it_Visit != variable_elimination_schedule.end();
         ++it_Visit)
    {                                                             
        const type_map_uint_to_Event::iterator event_being_summed_out = events.find( it_Visit->event_for_elimination->UID );  // readability.         
        
        if (my_MPI_rank == 0)                
	{
            std::fprintf(stderr, "\n\n\n\nevent_being_summed_out UID = %u\n\tprofile_length = %u\n\tremaining_coordinates_for_consideration.size() = %u\n\n",
                                event_being_summed_out->second.UID,
                                event_being_summed_out->second.profile_length,
                                event_being_summed_out->second.remaining_coordinates_for_consideration.size()  );  
	}
        
    
        if (my_MPI_rank == 0)
        {            
            char dest_dirname[option_size];                                       
            std::sprintf(dest_dirname, "%s/Collect_Likelihoods_for_Event_%u",
                                        scratch_job_output_dir.c_str(),
                                        event_being_summed_out->second.UID   );          
        
            boost::filesystem::create_directory(boost::filesystem::path(dest_dirname));        
        }
        world_MPI_boost_ptr->barrier();  // so that the directory is created in time!        
        
        
	
	
	
	
	


           
        if (my_MPI_rank == 0)                         
	{  event_being_summed_out->second.print_this_entry();  }
                                
        
        // just in case we need them, we define these here to have them in scope...
	type_map_uint_to_Breakpoint_complex  relevant_Breakpoint_complex__neighbors_and_self;
	type_map_uint_to_set_uint relevant_variational_positions__neighbors_and_self;
	
	
        type_set_uint haploid_UID_arguments_for_new_PSF;                                                
        
        
        

        
        
        if (!event_being_summed_out->second.remaining_coordinates_for_consideration.empty())
        {                        
            event_being_summed_out->second.pre_upload_homologous_regions_from_the_Reference_genome(); 
                                   
            time_begin = omp_get_wtime();
            if (my_MPI_rank == 0)
	    {  std::fprintf(stderr, "rank %d:  partition_DAR_by_homology    and    create_natural_poisson_intervals.\n", my_MPI_rank);  }
            
            event_being_summed_out->second.upload_or_create_natural_poisson_intervals___gender_INsensitive();

            if (my_MPI_rank == 0)
	    {  std::cerr<< "\n\nrank " << my_MPI_rank <<  "   time    =     "  << omp_get_wtime() - time_begin << "\n\n";  }
		
	
            check_that_relevant_regions_are_well_defined__and__hide_undefined_regions_if_necessary(event_being_summed_out, events);
	}//remaining coords
	    
	    
	    
	
	
	std::cerr << "\tcalculate_total_GC_rate_for_every_base__by_summing_over_all_Readgroup_GC_rates_at_each_position.\n";		
	event_being_summed_out->second.calculate_total_HAPLOID_GC_rate_for_every_base_by_summing_over_all_Readgroup_GC_rates_at_each_position(); 
	
	
	std::cerr << "\tcalculate_base_diploid_fragmentation_rate___haploid_sensitive___gender_sensitive__and_base_cumulative_haploid_fragmentation_rate.\n";		
	event_being_summed_out->second.calculate_base_diploid_fragmentation_rate___haploid_sensitive___gender_sensitive__and_base_cumulative_haploid_fragmentation_rate();       
	
	

	print_map_keys_and_values<uint, longdouble>(event_being_summed_out->second.base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive,
						    "base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive", NULL, true);
	
								
							
	time_begin = omp_get_wtime();
	std::cerr << "\tcalculate_observed_read_depth_for_natural_poisson_intervals\n";
	
	event_being_summed_out->second.calculate_observed_read_depth_for_natural_poisson_intervals(my_BAM_reader);
	
	std::cerr<< "time to calculate observed read depth  =   "  << omp_get_wtime() - time_begin << " s\n"; 
				

	
	
	std::cerr<< "\tset_up_all_PERs:\n";
	if (!event_being_summed_out->second.remaining_coordinates_for_consideration.empty())
	{//set up
	    const int success_setup = (int)event_being_summed_out->second.set_up_all_PERs(consider_GeneConversion_in_Full_compute, my_BAM_reader);
	    
	    if (!success_setup)
	    {
		if (my_MPI_rank == 0)
		{
		    std::stringstream reason_strm;
		    reason_strm << "failed to set up all PERs for Event " << event_being_summed_out->second.UID;                 
		    append_to_skipped_connected_component_file( *this, reason_strm.str()  );
		}
		return -1;
	    }
	}//set up
    
	
					    
	
	std::cerr << "\tdetermine_homologous_variational_positions_on_immediate_neighbors_and_self...\n";
	relevant_variational_positions__neighbors_and_self = event_being_summed_out->second.determine_homologous_variational_positions_on_immediate_neighbors_and_self();
	
	for (type_map_uint_to_set_uint::const_iterator it_ev_rel = relevant_variational_positions__neighbors_and_self.begin();
		it_ev_rel != relevant_variational_positions__neighbors_and_self.end();
		++it_ev_rel)
	{
	    relevant_Breakpoint_complex__neighbors_and_self[it_ev_rel->first].variational_positions__ALL = it_ev_rel->second;
	}
	
	
	std::cerr << "\tdetermine_relevant_NAHR_and_GeneConversion_breakpoints_for_each_interlocking_event_and_self_according_to_relevant_variational_positions...\n";
	
	determine_relevant_NAHR_and_GeneConversion_breakpoints_for_each_interlocking_event_and_self_according_to_relevant_variational_positions(
									*this,
									relevant_variational_positions__neighbors_and_self,
									relevant_Breakpoint_complex__neighbors_and_self);
	
	
	{//print
	    std::stringstream print_strm;
	    for (type_map_uint_to_Breakpoint_complex::const_iterator it_pev = relevant_Breakpoint_complex__neighbors_and_self.begin();
		    it_pev != relevant_Breakpoint_complex__neighbors_and_self.end();
		    ++it_pev)
	    {
		print_strm << "\tEvent:  " << it_pev->first << "\n";
		it_pev->second.print_this_Breakpoint_complex(&print_strm);
	    }
	    
	    std::cerr << "\n" << print_strm.str() << "\n";		
	}//print		
		
				
							    
	haploid_UID_arguments_for_new_PSF = extract_keys_of_map_and_return_as_set<uint, Breakpoint_complex>
						(relevant_Breakpoint_complex__neighbors_and_self); // INCLUDES the eliminating Event UID !!!                                                                             
    
						

	print_set<uint>(haploid_UID_arguments_for_new_PSF, "haploid_UID_arguments_for_new_PSF");
                

 
        
        
  
        
        
        
        
        event_being_summed_out->second.uploaded_pieces_of_chromos.clear();
        
        
                
        
        
        
        
        
        
        
        
        
        //Partial Sum preparation:
        
        //Identify relevant partial sum functions
        type_list_Partial_sum_fcn_iter relevant_partial_sum_functions;
        for (type_list_Partial_sum_fcn::iterator it_all_psf = partial_sum_functions.begin();
                it_all_psf != partial_sum_functions.end();
                ++it_all_psf)        
	{
            if (it_all_psf->check_if_event_is_an_argument_to_this_function(event_being_summed_out->second.UID))
	    {
                relevant_partial_sum_functions.push_back(it_all_psf);
	    }
	}
            
        if (my_MPI_rank == 0)
	{  fprintf(stderr, "debug:  relevant_partial_sum_functions identified.  size = %u\n", relevant_partial_sum_functions.size() );  }                  
        
                                
        
        it_Visit->my_Forward_partial_sum_function.set_eliminating_Event__and__haploid_argument_UIDs_according_to_relevant_PSFs(
                                                        event_being_summed_out->second.UID,
                                                         relevant_partial_sum_functions,
                                                         haploid_UID_arguments_for_new_PSF);
                                                         

         
                                                         
  
                        
        
        //load permissibile state vectors:
        if (my_MPI_rank == 0)
	{  std::fprintf(stderr, "\n\ndebug:  read_permissible_state_vectors_from_file_for_summing_out_a_given_UID = %u\n", event_being_summed_out->second.UID);  }
        
        {
            const bool success_read_state_vectors = read_permissibile_state_vectors_from_file_for_summing_out_a_given_UID_and_tack_on_empty_state_vector((*it_Visit));                
            if (!success_read_state_vectors)
            {
                if (my_MPI_rank == 0)
                {
                    std::stringstream  reason_strm;
                    reason_strm << "failed to read permissible state vectors from file for Event " << event_being_summed_out->second.UID;
                    append_to_skipped_connected_component_file(  *this,  reason_strm.str() );
                }
                return -3;            
            }
        }
        
        
        
        
        
        {//sanity check     
            if (it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs != it_Visit->neighbors)
            {
                std::stringstream error_strm;
                print_line_of_markers("ERROR! ", &error_strm);
                print_line_of_markers("(", &error_strm);
                error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                
                error_strm << "\n\nERROR!   recorded neighbors from Visitation schedule do not match with neighbors during \"perform complete variable elimination schedule\" "  << "\n\n";
                
                print_set<uint>(it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs, "it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs", &error_strm);
                print_set<uint>(it_Visit->neighbors, "it_Visit->neighbors", &error_strm);
                                
                
                print_line_of_markers(")", &error_strm);
                std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                        
            }
        
        }//sanity check
                
                
                
                
                
                
                
                
                
        
        
        it_Visit->refine_state_vectors_according_to_excluded_Events(excluded_Events);	
                    
            
	
        
        
        
        std::vector<type_vector_Sparse_map::const_iterator>   iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors;
        std::vector<type_vector_Sparse_map::const_iterator>   iterators_to_haploid_1_X_or_Y_sensitive_state_vectors; 		
	
	make_gender_sensitive_iterators_for_state_vectors(
					    iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors,
					    iterators_to_haploid_1_X_or_Y_sensitive_state_vectors,
					    events,
					    it_Visit);
					    
	if (!consider_GeneConversion_in_Full_compute)
	{
	    eliminate_GeneConversion_event_outcomes(iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors);
	    eliminate_GeneConversion_event_outcomes(iterators_to_haploid_1_X_or_Y_sensitive_state_vectors);	    	
	}// !consider_GeneConversion_in_Full_compute

                         
        
        const int size_state_vec_for_haploid_0 = iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.size();
        const int size_state_vec_for_haploid_1 = iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.size();
        
        
        
       
        const int haploid_state_compute_size = it_Visit->permissible_haploid_state_vectors_sparse.size();
        
        if (my_MPI_rank == 0)
	{
            std::cerr << "\nthaploid_state_compute_size = " << haploid_state_compute_size
			<< "\n\tOLD  it_Visit->diploid_compute_size  = " << it_Visit->diploid_compute_size
			<< "\n\tNEW gender sensitive compute size =  " << size_state_vec_for_haploid_0 
			<< "  *  " << size_state_vec_for_haploid_1
			<< "  =  " << size_state_vec_for_haploid_0  *  size_state_vec_for_haploid_1  << "\n\n\n";
	}
        
                        
        const bool  compute_size_small_enough_to_display_diagnostic_information  =    (size_state_vec_for_haploid_0  *  size_state_vec_for_haploid_1   <  150);                                
        
        if (my_MPI_rank == 0)
	{
            if (!compute_size_small_enough_to_display_diagnostic_information)
                std::cerr << "\n\tdiploid compute size for Event " << event_being_summed_out->second.UID
                        << "  is " << size_state_vec_for_haploid_0  *  size_state_vec_for_haploid_1
                        << ",  which is too big to display diagnostics (will slow down the compute dramatically)\n\n";   
	}
                        
			

	for (type_map_string_to_PER::iterator it_PER = event_being_summed_out->second.PERs_on_this_profile.begin();
		it_PER != event_being_summed_out->second.PERs_on_this_profile.end();
		++it_PER)
	{
	    bool has_empty = false;
	    
	    for (type_map_uint_to_BI::const_iterator it_mr_no = it_PER->second.map_MiniID__to__nohybrid_almgt_interval_absolute.begin();
		    it_mr_no != it_PER->second.map_MiniID__to__nohybrid_almgt_interval_absolute.end();
		    ++it_mr_no)
	    {
		if (BOOST_empty(it_mr_no->second))
		{  has_empty = true;  break;  }
			
	    }//no
	    
		if (has_empty)
		{
		    std::stringstream error_strm; 
		    error_strm << "\nit_PER = " << it_PER->first << "\n";
		    print_map_keys_and_values<uint, BOOST_Interval>(it_PER->second.map_MiniID__to__nohybrid_almgt_interval_absolute, "it_PER->second.map_MiniID__to__nohybrid_almgt_interval_absolute", &error_strm, false);
		}		    
	}
			
			
		
			
			
			
                        
        std::cerr.flush();std::fflush(stderr);                  
                         
        Partial_sum_function rank_specific___piece_of_Visitation_PSF(it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs,
                                                                     it_Visit->event_for_elimination->UID);  
		// This is the partial sum term:  Sum_{states of Eliminating Event} (  all likelihood terms Phi that involve the eliminiating Event , times priors and other PSFs involving eliminating Event    )
		// It is rank specific because each MPI rank computes different parts of the state space (states = outcome of Events).
        
	Partial_sum_function rank_specific___piece_of_new_PSF(haploid_UID_arguments_for_new_PSF,
								relevant_variational_positions__neighbors_and_self);
		// This is the product of likelihood terms likelihood involving eliminating_Event:
		//    Prod_{i such that genome (i.e. data) partition G_i (i.e. D_i) involves eliminating Event}  (    P(D_i | E_1, ..., En)    )              
    
    
    
	const double time_states_loop_begin = omp_get_wtime();
	
	for (int state_vec_hap_0_ctr = my_MPI_rank;
		state_vec_hap_0_ctr < size_state_vec_for_haploid_0;
		state_vec_hap_0_ctr += size_MPI_WORLD )
	for (int state_vec_hap_1_ctr = 0;
		state_vec_hap_1_ctr < size_state_vec_for_haploid_1;
		++state_vec_hap_1_ctr)
	{	    
	    const type_vector_Sparse_map::const_iterator it_states_0 = iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.at(state_vec_hap_0_ctr);
	    const type_vector_Sparse_map::const_iterator it_states_1 = iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.at(state_vec_hap_1_ctr);
	    
	    if (test_if_diploid_state_vectors_contain_non_identifiable_diploid_outcomes(*it_states_0, *it_states_1))
	    {  continue;  }
	    	    
	    if (!event_being_summed_out->second.remaining_coordinates_for_consideration.empty())
	    {                        
		
		//sum over breakpoints:  
		type_map_uint_to_vector_BI  breakpoint_space__hap_0;
		type_map_uint_to_vector_BI  breakpoint_space__hap_1;

		const uint number_of_haploid_breakpoints_0
			= create_haploid_breakpoint_space
				    (breakpoint_space__hap_0,
				    *it_states_0,
				    relevant_Breakpoint_complex__neighbors_and_self);
		    
		const uint number_of_haploid_breakpoints_1					
			= create_haploid_breakpoint_space
				    (breakpoint_space__hap_1,
				    *it_states_1,
				    relevant_Breakpoint_complex__neighbors_and_self);
				    
				    
		rank_specific___piece_of_new_PSF.allocate_diploid_breakpoints_for_given_diploid_states(
								*it_states_0, *it_states_1,
								breakpoint_space__hap_0, breakpoint_space__hap_1);
		
						    
				    
		{//scope
		    type_map_uint_to_BI breakpoints__haploid_0;
		    initialize_map_with_keys_from_singlecontainer<uint, BOOST_Interval, type_set_uint>(
					breakpoints__haploid_0,
					extract_keys_of_map_and_return_as_set<uint,uint>(it_states_0->sparse_map_UID_to_value), empty_BI);
		    
		    type_map_uint_to_BI breakpoints__haploid_1;
		    initialize_map_with_keys_from_singlecontainer<uint, BOOST_Interval, type_set_uint>(
					breakpoints__haploid_1,
					extract_keys_of_map_and_return_as_set<uint,uint>(it_states_1->sparse_map_UID_to_value), empty_BI);
				    
		    set_breakpoints_using_breakpoint_space(breakpoints__haploid_0, breakpoint_space__hap_0, 0);
		    set_breakpoints_using_breakpoint_space(breakpoints__haploid_1, breakpoint_space__hap_1, 0);                            
											
		    //read depth:                                        
		    if (compute_size_small_enough_to_display_diagnostic_information)
		    {
			std::stringstream debug_strm;
			debug_strm << "rank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
			print_map_keys_and_values<uint, uint>(it_states_0->sparse_map_UID_to_value, "it_states_0", &debug_strm);                     
			print_map_keys_and_values<uint, uint>(it_states_1->sparse_map_UID_to_value, "it_states_1", &debug_strm);			
			
			event_being_summed_out->second.calculate_or_get_Read_Depth_for_these_state_vectors_and_breakpoints
									(*it_states_0, *it_states_1,
									breakpoints__haploid_0, breakpoints__haploid_1,
									&debug_strm);                                                                                                                             
		    
			std::fprintf(stderr, "\n\n%s\n\n", debug_strm.str().c_str() );
		    }
		}//scope       					
		    
		    
		#pragma omp parallel
		{//parallel			
		    type_map_uint_to_BI breakpoints__haploid_0;
		    initialize_map_with_keys_from_singlecontainer<uint, BOOST_Interval, type_set_uint>(
					breakpoints__haploid_0,
					extract_keys_of_map_and_return_as_set<uint,uint>(it_states_0->sparse_map_UID_to_value), empty_BI);
		    
		    type_map_uint_to_BI breakpoints__haploid_1;
		    initialize_map_with_keys_from_singlecontainer<uint, BOOST_Interval, type_set_uint>(
					breakpoints__haploid_1,
					extract_keys_of_map_and_return_as_set<uint,uint>(it_states_1->sparse_map_UID_to_value), empty_BI);
		
		    #pragma omp for collapse(2) schedule(dynamic,10)
		    for (uint brkpt_ctr_0 = 0; brkpt_ctr_0 < number_of_haploid_breakpoints_0; ++brkpt_ctr_0)
		    {
			for (uint brkpt_ctr_1 = 0; brkpt_ctr_1 < number_of_haploid_breakpoints_1; ++brkpt_ctr_1)
			{  
			    set_breakpoints_using_breakpoint_space(breakpoints__haploid_0, breakpoint_space__hap_0, brkpt_ctr_0);
			    set_breakpoints_using_breakpoint_space(breakpoints__haploid_1, breakpoint_space__hap_1, brkpt_ctr_1);
			    
			    //read depth:               
// 			    double time_begin_check = omp_get_wtime();
			    
			    real prod_PER__diploid_sum_m__mPER(
						event_being_summed_out->second.calculate_or_get_Read_Depth_for_these_state_vectors_and_breakpoints
										(*it_states_0, *it_states_1,
										breakpoints__haploid_0, breakpoints__haploid_1));
			    
// 			    std::cerr << "\ttime RD calc:  " << omp_get_wtime() - time_begin_check << "\n";			    
// 			    time_begin_check = omp_get_wtime();
			    
			    const double PER_loop_time_begin = omp_get_wtime();
			    for (type_map_string_to_PER::iterator it_PER = event_being_summed_out->second.PERs_on_this_profile.begin();
				    it_PER != event_being_summed_out->second.PERs_on_this_profile.end();
				    ++it_PER)
			    {
				const type_uint__real haploid_0__sum_P_mPER(
						    it_PER->second.calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints
											(false,
											*it_states_0,
											breakpoints__haploid_0));
											
				const type_uint__real haploid_1__sum_P_mPER(
						    it_PER->second.calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints
											(true,
											*it_states_1,
											breakpoints__haploid_1));
				
				const uint total_diploid_contributions =  haploid_0__sum_P_mPER.first + haploid_1__sum_P_mPER.first; 
					    //i.e. this is a prior on the number of mapping locations (~ # LCRs) that a read could have come from.
				
				switch (total_diploid_contributions)
				{
				    case 0: 
					break;
				    case 1:
					prod_PER__diploid_sum_m__mPER *= (haploid_0__sum_P_mPER.second + haploid_1__sum_P_mPER.second);
					break;                                            
				    default:
					prod_PER__diploid_sum_m__mPER *= (haploid_0__sum_P_mPER.second + haploid_1__sum_P_mPER.second);
					prod_PER__diploid_sum_m__mPER /=  total_diploid_contributions; 
					break;                                                                        
				}   
			    }  // end for-loop   PER
			    
// 			    std::cerr << "\ttime PER-loop:  " << omp_get_wtime() - time_begin_check << "\n";
			    
														
                            //this is omp safe:
			    rank_specific___piece_of_new_PSF.save_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints
										    (*it_states_0, *it_states_1,
										    breakpoints__haploid_0, breakpoints__haploid_1,
										    prod_PER__diploid_sum_m__mPER);                                                                
							
			}  // end  while   breakpoint_iterator_1                            
		    }  // end  while   breakpoint_iterator_0
		    
		    mpfr_free_cache();
		}//parallel
	    } //remaining_coordinates  is not empty

	    const double time_cart_begin = omp_get_wtime();
	    rank_specific___piece_of_Visitation_PSF.
		    save_Cartesian_product_of_related_PSFs_wrt_eliminating_Event_and_include_probabilities_of_Breakpoint_and_Event
									(event_being_summed_out->second,
									*it_states_0,   *it_states_1,                                           
									relevant_partial_sum_functions,
									rank_specific___piece_of_new_PSF);
									
	    if (compute_size_small_enough_to_display_diagnostic_information)
	    {  std::cerr << "time Cartesian-product:  " << omp_get_wtime() - time_cart_begin << "  s.\n";  }
									
	}     // states
	
	
	std::cerr << "time states loop = " << omp_get_wtime() - time_states_loop_begin << "  s.\n";
	
	
	const double time_housekeeping_begin = omp_get_wtime();
	
			    
			    
	rank_specific___piece_of_new_PSF.copy_individual_likelihood_files_to_common_likelihood_directory_for_this_Event(
											    event_being_summed_out->second.UID);
	rank_specific___piece_of_new_PSF.clear_data();
	

	mpfr_free_cache();    






        
        iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.clear();
        iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.clear(); 


        haploid_UID_arguments_for_new_PSF.clear();

        
        
        event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.clear();
        event_being_summed_out->second.base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive.clear();
	event_being_summed_out->second.base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.clear();
        event_being_summed_out->second.observed_read_depth_for_natural_poisson_intervals.clear();         

        event_being_summed_out->second.PERs_on_this_profile.clear();
        event_being_summed_out->second.filter_for_NON_large_gap_of_neighbors_and_self__that_affect_Read_depth.clear();
        event_being_summed_out->second.RD_affecting_neighbors_self.clear();
        
        event_being_summed_out->second.large_gaps_of_profile.clear();
        
                        
        
        it_Visit->neighbors.clear();
        
        it_Visit->permissible_haploid_state_vectors_sparse.clear();

            

        
        
        world_MPI_boost_ptr->barrier(); 
        if (my_MPI_rank == 0)
        {
            Partial_sum_function::combine_individual_likelihood_files_into_common_likelihood_file_for_this_Event(event_being_summed_out->second.UID);
        }
        
        
        
        
        world_MPI_boost_ptr->barrier();  //because  before we erase an old PSF, all ranks must be done using it (i.e. done with above loop).
        
        std::cerr <<"\n\terasing old PSFs...\n";
        
        
        
        
        //Partial sum function:                        
        for (type_list_Partial_sum_fcn_iter::iterator it_old_psf = relevant_partial_sum_functions.begin();
                it_old_psf != relevant_partial_sum_functions.end();
                ++it_old_psf)          
        {
            partial_sum_functions.erase(*it_old_psf);        
        }
        
        
        
        
        
        


        if (size_MPI_WORLD > 1)
	{  std::cerr << "\nrank  " << my_MPI_rank << " is waiting at ::gather...\n";  }
        
        
        
	{//all_gather
	    std::vector<  Partial_sum_function >  gathered_PSFs_of_Visitation;	
	    boost::mpi::all_gather<Partial_sum_function>( *world_MPI_boost_ptr, 
							rank_specific___piece_of_Visitation_PSF,
							gathered_PSFs_of_Visitation );
							
	    rank_specific___piece_of_Visitation_PSF.clear_data();
			
	    for (int rank=0; rank < size_MPI_WORLD; ++rank)
	    {
		it_Visit->my_Forward_partial_sum_function.absorb_piece_of_Partial_Sum_Function_into_one_Partial_Sum_Function(
							gathered_PSFs_of_Visitation.at(rank));
		
		gathered_PSFs_of_Visitation.at(rank).clear_data();
	    }	
	    	    	    	 						  
	}// gather
						  

              
	
	
	
        if (my_MPI_rank == 0)	
	{  std::cerr << "\nit_Visit->my_Forward_partial_sum_function.marginalize_over_eliminating_Event_breakpoints()\n";  }
	
	it_Visit->my_Forward_partial_sum_function.marginalize_over_eliminating_Event_breakpoints_and_outcomes();  	
            
            
	
	
	
                
        if (my_MPI_rank == 0) 
        {
            std::cerr << "\n\nrank " << my_MPI_rank << "  it_Visit->save_PSF_to_archive()... ";
            it_Visit->save_PSF_to_archive(); 
            std::cerr << "it_Visit->save_PSF_to_archive()-----done!\n\n";
        }                
                

        
        
        
        
        
        
        //housekeeping:
        if (my_MPI_rank == 0)
	{  print_line_of_markers("housekeeping ");  }
                    

        
        
        
                                                                 
        
    

	
        
     
                    
        
        
        { // copy of PSF
            partial_sum_functions.push_back(it_Visit->my_Forward_partial_sum_function);
            const type_list_Partial_sum_fcn::iterator  propagating_PSF  = --partial_sum_functions.end();
	    	    	    
            
            std::cerr << "\n\nrank " << my_MPI_rank  << " remove_eliminated_Event_from_storage___MPI \n";            
            propagating_PSF->remove_eliminated_Event_from_storage(event_being_summed_out->second.UID);
	    

	    std::cerr << "\n\nrank " << my_MPI_rank  
			<< " make_numerically_stable_by_multiplying_every_partial_sum_value_by_the_inverse_of_the_smallest_partial_sum_value\n ";	    
	    propagating_PSF->make_numerically_stable_by_multiplying_every_partial_sum_value_by_the_inverse_of_the_smallest_partial_sum_value();     
	    // We do NOT need to keep track of this normalizing constant!!!!!
	    
// 	    propagating_PSF->print_this_PSF__states_only();                    		    
        } // copy of PSF
        
        
        
        
        
        
        
        it_Visit->my_Forward_partial_sum_function.clear_data();


      
        if (my_MPI_rank == 0)        
	{  std::fprintf(stderr, "debug:  erasing positions from all related events...\n");  }
        //FINALLY,   ERASE GENOME POSITIONS FROM ALL RELATED EVENTS (SO THE SAME GENOME POSITIONS ARENT CONSIDERED MORE THAN ONCE)
        //What about genome positions of the current event that weren't on anybody elses homology lists?  (i.e. unique to this LCR-pair)???????                        
                        
        type_map_uint_to_Interlock* local_and_global_list[2];
        local_and_global_list[0] = &(event_being_summed_out->second.local_interlocking_Events);
        local_and_global_list[1] = &(event_being_summed_out->second.global_interlocking_Events);

        for (uint lg = 0; lg<2; ++lg)                           
            for (type_map_uint_to_Interlock::iterator ev_it = local_and_global_list[lg]->begin();
                    ev_it != local_and_global_list[lg]->end();
                    ++ev_it)   // for every interlocking event...
            {
                
                const type_map_uint_to_Event::iterator it_interlocking_Event =  events.find(ev_it->first);
//                 Event *const interlocking_Event_ptr = &(events.at(ev_it->first));
                                                                
                //erase yourself from this interlocking event's lists:
                it_interlocking_Event->second.local_interlocking_Events.erase(event_being_summed_out->second.UID);
                it_interlocking_Event->second.global_interlocking_Events.erase(event_being_summed_out->second.UID);                
                
            
                    
                //and erase your very coordinates.  (necessary, for say, Inbetween region).
                if (event_being_summed_out->second.chromos[0] ==  it_interlocking_Event->second.chromos[0])  
                {
                    for (uint pos =  event_being_summed_out->second.region_between_and_including_the_LCRs_themselves.lower();
                              pos <= event_being_summed_out->second.region_between_and_including_the_LCRs_themselves.upper();
                              ++pos)                   
                        it_interlocking_Event->second.remaining_coordinates_for_consideration.erase(pos);               
                }
                
                
                    
                    
//                 for (type_set_uint::const_iterator it_remaining = event_being_summed_out->second.remaining_coordinates_for_consideration.begin();
//                         it_remaining != event_being_summed_out->second.remaining_coordinates_for_consideration.end();
//                         ++it_remaining)
//                     it_interlocking_Event->second.remaining_coordinates_for_consideration.erase(*it_remaining);
                
                
                
                
                for (type_list_Coords_and_homologous_region::const_iterator it_homol
                                = event_being_summed_out->second.regions_homologous_to_directly_affected_region.begin();
                        it_homol != event_being_summed_out->second.regions_homologous_to_directly_affected_region.end();
                        ++it_homol)              
                    if (it_homol->chromosome_of_homologous_region == it_interlocking_Event->second.chromos[0]
                        and  BOOST_overlap( get_endpoints_of_keys_of_compressed_map(it_homol->compressed_map_homologous_region_to_profile),
                                            it_interlocking_Event->second.region_between_and_including_the_LCRs_themselves )      )
                    {
                        const type_map_uint_to_uint full_map( convert_compressed_map_to_full_map(it_homol->compressed_map_homologous_region_to_profile) );                        
                        for (type_map_uint_to_uint::const_iterator map_it = full_map.begin();
                                map_it != full_map.end();
                                ++map_it)                                                
                            it_interlocking_Event->second.remaining_coordinates_for_consideration.erase( map_it->first );                                                                      
                    }  
            
            
            
                           
                           
                           
                {//homol reg restriction
                    type_list_Coords_and_homologous_region::iterator it_homol
                                = it_interlocking_Event->second.regions_homologous_to_directly_affected_region.begin();                                
                    while(  it_homol != it_interlocking_Event->second.regions_homologous_to_directly_affected_region.end()  )
                    {
                        it_homol->compressed_map_homologous_region_to_profile
                                =  get_compressed_map_inverse_of_compressed_map(
                                        convert_full_map_to_compressed_map(
                                                get_map_intersection_of_compressed_map_keys_and_set(
                                                        get_compressed_map_inverse_of_compressed_map(   it_homol->compressed_map_homologous_region_to_profile    ),
                                                        it_interlocking_Event->second.remaining_coordinates_for_consideration )  )  ); 
                                            
                        if ( it_homol->compressed_map_homologous_region_to_profile.at(0).empty() )
                            it_homol = it_interlocking_Event->second.regions_homologous_to_directly_affected_region.erase( it_homol );
                        else
                            ++it_homol;                                            
                    }                                                        
                }//homol reg restriction                           
             
            
            
            
            
            
                //and erase variational positions from its profile.
                if (!event_being_summed_out->second.remaining_coordinates_for_consideration.empty())
                {
                    if (relevant_variational_positions__neighbors_and_self.count(it_interlocking_Event->second.UID) > 0)
                    {
                        for (type_set_uint::const_iterator it_rel_varpos = relevant_variational_positions__neighbors_and_self.at(it_interlocking_Event->second.UID).begin();
                                it_rel_varpos != relevant_variational_positions__neighbors_and_self.at(it_interlocking_Event->second.UID).end();
                                ++it_rel_varpos)
			{  it_interlocking_Event->second.variational_positions_of_profile.erase(*it_rel_varpos);  }
                    }
                    else
                    {
                        std::stringstream warning_strm;
                        print_line_of_markers("WARNING! ", &warning_strm);
                        print_line_of_markers("(", &warning_strm);
                        warning_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                        
                        warning_strm << "\n\nWARNING!    in   \" erase variational positions from its profile.\"  area,    relevant_variational_positions__neighbors_and_self.count(it_interlocking_Event->second.UID)  ==  0      for     it_interlocking_Event->second.UID  ="  <<  it_interlocking_Event->second.UID  << "\n\n";
                        
                        print_line_of_markers(")", &warning_strm);
                        std::fprintf(stderr, "\n\n%s\n\n", warning_strm.str().c_str() );                    
                    }                        
                }
                
                    
                for (uint qs=0; qs<2; ++qs)     
		{
                    it_interlocking_Event->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)                        
                        =  get_compressed_map_inverse_of_compressed_map(
                                convert_full_map_to_compressed_map(
                                    get_map_intersection_of_compressed_map_keys_and_set(
                                            get_compressed_map_inverse_of_compressed_map(
                                                    it_interlocking_Event->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)),
                                            it_interlocking_Event->second.variational_positions_of_profile)));   
		}
                                                        
                            
                                        
                if (my_MPI_rank == 0)  
                {                		    
                    std::fprintf(stderr, "UID %u", it_interlocking_Event->second.UID);                        
                    print_set<uint>(it_interlocking_Event->second.variational_positions_of_profile, "it_interlocking_Event->second.variational_positions_of_profile");
                }
                
//                 //and erase your homologous coordinates:
//                 for (type_list_vector_vector_int::const_iterator compress_it = ev_it->second.compressed_maps_homol_regions_to_profile.begin();
//                         compress_it != ev_it->second.compressed_maps_homol_regions_to_profile.end();  //maps from interlocking event ev_it to event_being_summed_out
//                         ++compress_it)
//                 {                    
//                     type_map_uint_to_uint full_map_homol_to_profile(  convert_compressed_map_to_full_map(*compress_it)  ); 
//                      
//                     for (type_map_uint_to_uint::const_iterator map_it = full_map_homol_to_profile.begin();
//                             map_it != full_map_homol_to_profile.end();
//                             ++map_it)                                                
//                         it_interlocking_Event->second.remaining_coordinates_for_consideration.erase( map_it->first  );                    
//                 } // end for-loop   compress_it                                                                
                        
            }  // end for-loop   ev_it
                      
                      
                      
                      
                      
        
        // clean up:
        if (my_MPI_rank == 0)
	{  std::fprintf(stderr, "debug:  deleting  rest of  vectors...\n");  }    
	
        //This variable has been summed out and no longer contains any information.  delete all of its vectors.
        event_being_summed_out->second.global_interlocking_Events.clear();
        event_being_summed_out->second.local_interlocking_Events.clear();
        event_being_summed_out->second.remaining_coordinates_for_consideration.clear();
        event_being_summed_out->second.my_original_Breakpoint_complex.clear_all_data();
        
        event_being_summed_out->second.regions_homologous_to_LCR_profile.clear();
        event_being_summed_out->second.regions_homologous_to_directly_affected_region.clear();
        event_being_summed_out->second.variational_positions_of_profile.clear();  
        event_being_summed_out->second.compressed_map_LCR_to_profile__for_each_LCR.clear();
        event_being_summed_out->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.clear();    
	event_being_summed_out->second.absolute_coordinates_inside_of_large_gaps_of_profile.first.clear();
	event_being_summed_out->second.absolute_coordinates_inside_of_large_gaps_of_profile.second.clear();
        
	
	
	
	
        
        
        std::cerr << "time housekeeping = " << omp_get_wtime() - time_housekeeping_begin << "  s.\n";
                        
        
        if (my_MPI_rank == 0)
        {
            print_line_of_markers("*-+-*");
            print_line_of_markers("Visitation Schedule step complete!");  
            print_line_of_markers("*-+-*");
        }
        
    }//variable_elimination_schedule



    return 0;  //success

}  //  end of   perform_complete_variable_elimination_schedule

























//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
















bool Conn_Comp::HEURISTIC__test_each_specific_variational_position_possibility
                                                    (type_set_uint &excluded_Events,
						     const bool &consider_GeneConversion_breakpoints,
						    BamTools::BamReader &my_BAM_reader)
{   
    
    std::fprintf(stderr, "\n\n\nrank %d:   begin \"HEURISTIC__test_each_specific_variational_position_possibility\"...\n\n",
                 my_MPI_rank);                 
                     
    type_set_uint::iterator insert_excluded_Events_it = excluded_Events.begin();
            

    //set up PERs
    for (type_map_uint_to_Event::iterator it_ev = events.begin();
            it_ev != events.end();
            ++it_ev)
    {                   
        it_ev->second.PERs_on_this_profile.clear();
                
        {//find all PERs        
            {//search intervals and upload
		const double time_begin_upload = omp_get_wtime();
                                               
                const type_map_uint_to_list_BI search_intervals(
                                it_ev->second.create_search_intervals_from_DAR_and_homologous_regs_for_uploading_PERs
                                                (maximum_average_fragment_length, search_block_size)   );
                                                                                                                                                                                 
                                std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() 
                                            << "   for event " << it_ev->second.UID  
                                            << "   upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals...\n ";   
                            
                                            
                const bool success_in_upload = it_ev->second.upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals_____MPI_and_OMP
                                                        (search_intervals,  20, my_BAM_reader);
                                                        
                if (!success_in_upload)
		{  return false;  }
                                                                                         
		std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() 
			    << "   DONE WITH upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals.     time = "
			    << omp_get_wtime() -  time_begin_upload << "\n ";                                                                                        
            }//upload   
                                        

		        
	    std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() 
			<< "\n number_of_PERs = " << it_ev->second.PERs_on_this_profile.size()  << " \n";
    
	    const double time_orient = omp_get_wtime();                                                                    
    
	    std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "   for_all_PERs____set_region_sequences___and___orient...\n ";               
			                    
                    
            it_ev->second.pre_upload_homologous_regions_from_the_Reference_genome();                
             
            it_ev->second.for_all_PERs____set_region_sequences___and___orient___MPI_and_OMP();
            
            it_ev->second.uploaded_pieces_of_chromos.clear();
        
            
	    std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num()
			<<  "  done orienting.  for EV    " <<  it_ev->second.UID
			<< "  , time  =     "  << omp_get_wtime() - time_orient << "  seconds.\n";
        }//find all PERs
        
             
        

        
        
	
	set_up_and_align_PERs_using_GPU____OMP
			(it_ev->second.PERs_on_this_profile,
			 &it_ev->second,
			 true,//create_and_consider_MiniEvents,
			 consider_GeneConversion_breakpoints,  
			 false,//sum_over_all_unaffected_LCRs,  //"false" for heuristic
			 true,//special_save_breakpoints_actually_captured_in_alignment); //"true" for heuristic    		                      
			false);//do_PER_orient_formatting
                                                      
			
			
                                
        type_map_uint_to_uint map_NAHR_varpos_to_number_of_decent_reads__AB;
	type_map_uint_to_uint map_NAHR_varpos_to_number_of_decent_reads__BA;
		             
	{//get supported breakpoints
	    const std::vector<type_map_string_to_PER::iterator>  PERs___loop = get_loopable_iterators<type_map_string_to_PER>(it_ev->second.PERs_on_this_profile);
        
	    double time_begin = omp_get_wtime();                     
	    std::cerr << "\nheuristic on varpos...\n";
	    
	    mpfr_free_cache(); 
	    #pragma omp parallel
	    {//parallel     
	    
		type_map_uint_to_uint map_NAHR_varpos_to_number_of_decent_reads__AB____thread;
		type_map_uint_to_uint map_NAHR_varpos_to_number_of_decent_reads__BA____thread;
	    
		#pragma omp for  schedule(dynamic,20)  nowait
		for (uint j = 0; j < it_ev->second.PERs_on_this_profile.size(); ++j)
		{              
		    PERs___loop.at(j)->second.HEURISTIC__determine_which_breakpoints_are_strongly_supported_by_this_PER_for_a_given_Event____using_pseudo_posterior(
							it_ev->second.UID,
							heuristic_posterior_PROBABILITY_threshold,
							    map_NAHR_varpos_to_number_of_decent_reads__AB____thread,
							    map_NAHR_varpos_to_number_of_decent_reads__BA____thread);
		}//for PERs___loop                            
				
		
		#pragma omp critical (aggregate_possibly_supported_varpos_per_thread_into__main_sets__01)
		for (type_map_uint_to_uint::const_iterator it_poss_thread = map_NAHR_varpos_to_number_of_decent_reads__AB____thread.begin();
			it_poss_thread != map_NAHR_varpos_to_number_of_decent_reads__AB____thread.end();
			++it_poss_thread)
		{  map_NAHR_varpos_to_number_of_decent_reads__AB[it_poss_thread->first] += it_poss_thread->second;  }
		    
		    
		#pragma omp critical (aggregate_possibly_supported_varpos_per_thread_into__main_sets__10)
		for (type_map_uint_to_uint::const_iterator it_poss_thread = map_NAHR_varpos_to_number_of_decent_reads__BA____thread.begin();
			it_poss_thread != map_NAHR_varpos_to_number_of_decent_reads__BA____thread.end();
			++it_poss_thread)
		{  map_NAHR_varpos_to_number_of_decent_reads__BA[it_poss_thread->first] += it_poss_thread->second;  }
		
		mpfr_free_cache();
	    }//parallel
	    mpfr_free_cache(); 
	    
	    std::cerr  <<  "\ntime for heuristic posterior for Event  " <<  it_ev->second.UID << "  =   "  << omp_get_wtime() - time_begin << "  s.\n";  
	                                
	}//get supported breakpoints
            
            
        it_ev->second.PERs_on_this_profile.clear();
        
        
        
        
      
            
        
	print_map_keys_and_values<uint,uint>(map_NAHR_varpos_to_number_of_decent_reads__AB, "map_NAHR_varpos_to_number_of_decent_reads__AB", NULL, true);
	print_map_keys_and_values<uint,uint>(map_NAHR_varpos_to_number_of_decent_reads__BA, "map_NAHR_varpos_to_number_of_decent_reads__BA", NULL, true);
	
	
	
	
	
	
	//Now we have support varpos.  For Gene Conversion consistency, we REQUIRE >= 2 possible breakpoints per Event.	                                    
            
        //now find sufficiently good varpos  
	
    
	
    
	std::cerr <<  "  refining supported breakpoints...\n\n";   
	{//refine

	    if (it_ev->second.recomb_type == recomb_class__DupDel)
	    {		
	
		const type_set_uint supported_vp__AB(
					refine_breakpoints_using_heuristic(
								    it_ev->second, 
								    maximum_number_of_breakpoints_to_consider_after_heuristic,
								    map_NAHR_varpos_to_number_of_decent_reads__AB));
		
		const type_set_uint supported_vp__BA(
					refine_breakpoints_using_heuristic(
								    it_ev->second, 
								    maximum_number_of_breakpoints_to_consider_after_heuristic,
								    map_NAHR_varpos_to_number_of_decent_reads__BA));	
		

		it_ev->second.my_original_Breakpoint_complex.breakpoints_NAHR__AB = convert_set_of_uint_into_set_of_BI(supported_vp__AB);
		it_ev->second.my_original_Breakpoint_complex.breakpoints_NAHR__BA = convert_set_of_uint_into_set_of_BI(supported_vp__BA);
		
		
		
		type_map_BI_to_uint map_GeneConv_interval_to_decent_read_count__ABA;
		type_map_BI_to_uint map_GeneConv_interval_to_decent_read_count__BAB;
		
		for (type_set_uint::const_iterator it_AB = supported_vp__AB.begin();
			it_AB != supported_vp__AB.end();
			++it_AB)
		{
		    for (type_set_uint::const_iterator it_BA = supported_vp__BA.begin();
			    it_BA != supported_vp__BA.end();
			    ++it_BA)
		    {
			if (*it_AB < *it_BA)
			{
			    const BOOST_Interval new_GeneConv_interval(*it_AB, *it_BA);
			    it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA.insert(new_GeneConv_interval);
			    
			    if (map_NAHR_varpos_to_number_of_decent_reads__AB.count(*it_AB) > 0)
			    {  map_GeneConv_interval_to_decent_read_count__ABA[new_GeneConv_interval] += map_NAHR_varpos_to_number_of_decent_reads__AB.at(*it_AB);  }
			    
			    if (map_NAHR_varpos_to_number_of_decent_reads__BA.count(*it_BA) > 0)
			    {  map_GeneConv_interval_to_decent_read_count__ABA[new_GeneConv_interval] += map_NAHR_varpos_to_number_of_decent_reads__BA.at(*it_BA);  }	    
			    
			}
			else if (*it_BA < *it_AB)
			{
			    const BOOST_Interval new_GeneConv_interval(*it_BA, *it_AB);
			    it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB.insert(new_GeneConv_interval);
			    
			    if (map_NAHR_varpos_to_number_of_decent_reads__AB.count(*it_AB) > 0)
			    {  map_GeneConv_interval_to_decent_read_count__BAB[new_GeneConv_interval] += map_NAHR_varpos_to_number_of_decent_reads__AB.at(*it_AB);  }
			    
			    if (map_NAHR_varpos_to_number_of_decent_reads__BA.count(*it_BA) > 0)
			    {  map_GeneConv_interval_to_decent_read_count__BAB[new_GeneConv_interval] += map_NAHR_varpos_to_number_of_decent_reads__BA.at(*it_BA);  }	    
			}
		    }//10  
		}//01
		
		
		prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained<BOOST_Interval,uint, compare_BI>(
					it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA,
					invert_map_into_multimap<BOOST_Interval,uint, compare_BI>(map_GeneConv_interval_to_decent_read_count__ABA),
					maximum_number_of_breakpoints_to_consider_after_heuristic);
			
		prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained<BOOST_Interval,uint, compare_BI>(
					it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB,
					invert_map_into_multimap<BOOST_Interval,uint, compare_BI>(map_GeneConv_interval_to_decent_read_count__BAB),
					maximum_number_of_breakpoints_to_consider_after_heuristic);
			
		
		//if empty, just insert anything.
		if (it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA.empty())
		{
		    add_sparse_vp_to_form_Geneconversion_intervals
				(it_ev->second,
				maximum_number_of_breakpoints_to_consider_after_heuristic,
				map_NAHR_varpos_to_number_of_decent_reads__AB,
				map_NAHR_varpos_to_number_of_decent_reads__BA,
				it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA);
		}
		
		if (it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB.empty())
		{
		    add_sparse_vp_to_form_Geneconversion_intervals
				(it_ev->second,
				maximum_number_of_breakpoints_to_consider_after_heuristic,
				map_NAHR_varpos_to_number_of_decent_reads__BA,
				map_NAHR_varpos_to_number_of_decent_reads__AB,
				it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB);
		}		
	    }//dupdel
	    else
	    {//inv
	    
		//for inversion breakpoints:
		type_set_uint  refined_varpos__AB;
		type_map_uint_to_uint  refined_map_varpos_to_count_AB;
	    
		preliminary_refine_breakpoints_for_inversion_heuristic(
				    it_ev->second,
				    maximum_number_of_breakpoints_to_consider_after_heuristic,
				    map_NAHR_varpos_to_number_of_decent_reads__AB,
				    refined_varpos__AB,
				    refined_map_varpos_to_count_AB);
		
		
		type_set_uint  refined_varpos__BA;
		type_map_uint_to_uint  refined_map_varpos_to_count_BA;
	    
		preliminary_refine_breakpoints_for_inversion_heuristic(
				    it_ev->second,
				    maximum_number_of_breakpoints_to_consider_after_heuristic,
				    map_NAHR_varpos_to_number_of_decent_reads__BA,
				    refined_varpos__BA,
				    refined_map_varpos_to_count_BA);
		
		type_set_uint  double_supported_vp(get_set_intersection_of_two_sets<uint>(refined_varpos__AB, refined_varpos__BA));
		
// 		print_set<uint>(double_supported_vp, "double_supported_vp", NULL, true);
// 		print_set<uint>(refined_varpos__AB, "refined_varpos__AB", NULL, true);
// 		print_set<uint>(refined_varpos__BA, "refined_varpos__BA", NULL, true);
// 		print_map_keys_and_values<uint,uint>(refined_map_varpos_to_count_AB, "refined_map_varpos_to_count_AB", NULL, true);
// 		print_map_keys_and_values<uint,uint>(refined_map_varpos_to_count_BA, "refined_map_varpos_to_count_BA", NULL, true);
		
		type_map_uint_to_uint map_varpos_to_total_count;
		for (type_set_uint::const_iterator it_vp = double_supported_vp.begin();
			it_vp != double_supported_vp.end();
			++it_vp)
		{
		    map_varpos_to_total_count[*it_vp] = refined_map_varpos_to_count_AB.at(*it_vp) + refined_map_varpos_to_count_BA.at(*it_vp);
		}
		
		
		const type_set_uint supported_vp__INV(
					refine_breakpoints_using_heuristic(
								    it_ev->second, 
								    maximum_number_of_breakpoints_to_consider_after_heuristic,
								    map_varpos_to_total_count));										
		
		it_ev->second.my_original_Breakpoint_complex.breakpoints_NAHR__AB = convert_set_of_uint_into_set_of_BI(supported_vp__INV);
		it_ev->second.my_original_Breakpoint_complex.breakpoints_NAHR__BA = convert_set_of_uint_into_set_of_BI(supported_vp__INV);
		
		
		
		type_map_BI_to_uint map_GeneConv_interval_to_decent_read_count__ABA;
		type_map_BI_to_uint map_GeneConv_interval_to_decent_read_count__BAB;
		
		for (type_set_uint::const_iterator it_AB = refined_varpos__AB.begin();
			it_AB != refined_varpos__AB.end();
			++it_AB)
		{
		    for (type_set_uint::const_iterator it_BA = refined_varpos__BA.begin();
			    it_BA != refined_varpos__BA.end();
			    ++it_BA)
		    {
			if (*it_AB < *it_BA)
			{
			    const BOOST_Interval new_GeneConv_interval(*it_AB, *it_BA);
			    it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA.insert(new_GeneConv_interval);
			    
			    if (map_NAHR_varpos_to_number_of_decent_reads__AB.count(*it_AB) > 0)
			    {  map_GeneConv_interval_to_decent_read_count__ABA[new_GeneConv_interval] += map_NAHR_varpos_to_number_of_decent_reads__AB.at(*it_AB);  }
			    
			    if (map_NAHR_varpos_to_number_of_decent_reads__BA.count(*it_BA) > 0)
			    {  map_GeneConv_interval_to_decent_read_count__ABA[new_GeneConv_interval] += map_NAHR_varpos_to_number_of_decent_reads__BA.at(*it_BA);  }	    
			    
			}
			else if (*it_BA < *it_AB)
			{
			    const BOOST_Interval new_GeneConv_interval(*it_BA, *it_AB);
			    it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB.insert(new_GeneConv_interval);
			    
			    if (map_NAHR_varpos_to_number_of_decent_reads__AB.count(*it_AB) > 0)
			    {  map_GeneConv_interval_to_decent_read_count__BAB[new_GeneConv_interval] += map_NAHR_varpos_to_number_of_decent_reads__AB.at(*it_AB);  }
			    
			    if (map_NAHR_varpos_to_number_of_decent_reads__BA.count(*it_BA) > 0)
			    {  map_GeneConv_interval_to_decent_read_count__BAB[new_GeneConv_interval] += map_NAHR_varpos_to_number_of_decent_reads__BA.at(*it_BA);  }	    
			}
		    }//total
		}//total
		
		prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained<BOOST_Interval,uint, compare_BI>(
					it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA,
					invert_map_into_multimap<BOOST_Interval,uint, compare_BI>(map_GeneConv_interval_to_decent_read_count__ABA),
					maximum_number_of_breakpoints_to_consider_after_heuristic);
		
		prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained<BOOST_Interval,uint, compare_BI>(
					it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB,
					invert_map_into_multimap<BOOST_Interval,uint, compare_BI>(map_GeneConv_interval_to_decent_read_count__BAB),
					maximum_number_of_breakpoints_to_consider_after_heuristic);		
					
		//if empty, just insert anything.
		if (it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA.empty())
		{
		    add_sparse_vp_to_form_Geneconversion_intervals
				(it_ev->second,
				maximum_number_of_breakpoints_to_consider_after_heuristic,
				map_NAHR_varpos_to_number_of_decent_reads__AB,
				map_NAHR_varpos_to_number_of_decent_reads__BA,
				it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA);
		}
		
		if (it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB.empty())
		{
		    add_sparse_vp_to_form_Geneconversion_intervals
				(it_ev->second,
				maximum_number_of_breakpoints_to_consider_after_heuristic,
				map_NAHR_varpos_to_number_of_decent_reads__BA,
				map_NAHR_varpos_to_number_of_decent_reads__AB,
				it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB);
		}				
		
	    }//inv
	    
	    
	    //record vp
	    {
		const type_set_uint  vp_from_ab(
				extract_all_BOOST_Interval_endpoints_to_set(it_ev->second.my_original_Breakpoint_complex.breakpoints_NAHR__AB));
		
		const type_set_uint  vp_from_ba(
				extract_all_BOOST_Interval_endpoints_to_set(it_ev->second.my_original_Breakpoint_complex.breakpoints_NAHR__BA));		    		    
		
		const type_set_uint  vp_from_aba(
				extract_all_BOOST_Interval_endpoints_to_set(it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA));
		
		const type_set_uint  vp_from_bab(
				extract_all_BOOST_Interval_endpoints_to_set(it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB));

		it_ev->second.variational_positions_of_profile.clear();
		it_ev->second.variational_positions_of_profile.insert(vp_from_aba.begin(), vp_from_aba.end());
		it_ev->second.variational_positions_of_profile.insert(vp_from_bab.begin(), vp_from_bab.end());
		it_ev->second.variational_positions_of_profile.insert(vp_from_ab.begin(), vp_from_ab.end());
		it_ev->second.variational_positions_of_profile.insert(vp_from_ba.begin(), vp_from_ba.end());
		
		it_ev->second.my_original_Breakpoint_complex.variational_positions__ALL = it_ev->second.variational_positions_of_profile;        
	    }	    

	}//refine
	

	
// 	it_ev->second.variational_positions_of_profile = it_ev->second.my_original_Breakpoint_complex.variational_positions__ALL;
	
	

	
    
	std::cerr << "\n\t\t\tEvent " << it_ev->second.UID << ":\n";
	it_ev->second.my_original_Breakpoint_complex.print_this_Breakpoint_complex();
	
	
        
        for (uint qs=0; qs<2; ++qs)    
	{
            it_ev->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)                        
                =  get_compressed_map_inverse_of_compressed_map(  //LCR to profile varpos
                        convert_full_map_to_compressed_map(  
                            get_map_intersection_of_compressed_map_keys_and_set(
                                    get_compressed_map_inverse_of_compressed_map(  //profile varpos to LCR
                                            it_ev->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)  ),
                                    it_ev->second.variational_positions_of_profile)));            
	}
	
            
        if (my_MPI_rank == 0)
	{  std::cerr  << "DONE \"HEURISTIC__test_each_specific_variational_position_possibility\"   for   Event   " << it_ev->first << "\n\n";  }
        
    } // it_ev



    
    if (my_MPI_rank == 0)
    {  std::cerr << "\n\n\nDONE \"HEURISTIC__test_each_specific_variational_position_possibility\"!!!   For all Events!\n\n";  }
    
    
    return true;
    
}  //  HEURISTIC__test_each_specific_variational_position_possibility

































void Conn_Comp::output_reads_covering_certain_breakpoint_regions                             
				(const Sampled_diploid_Event_data &some_result,
				 const bool &do_likelihood_ratio_test_on_alignments,
				BamTools::BamReader &my_BAM_reader,
				const bool &display_ALL_reads_not_just_good_ones)				
{
    
    read_Connected_Component_ALL_Event_data_from_file(*this);
   
    const type_map_uint_to_Event::iterator Ev_it = events.find(some_result.event_UID);    
    
    const bool is_a_homozygous_Event_outcome = (some_result.the_diploid_Event_outcome.first == some_result.the_diploid_Event_outcome.second);
    const bool is_a_homozygous_breakpoint = is_a_homozygous_Event_outcome  and  BOOST_equal(some_result.the_diploid_profile_brkpts.first, some_result.the_diploid_profile_brkpts.second);
                  
    const uint RD_breakpoint_zygous_multiplier =    is_a_homozygous_breakpoint ?   2 : 1;
    const uint haploid_upper_bound  =    is_a_homozygous_breakpoint   ?   1 : 2;
    
//     std::cerr << "\n\n\ninside \"output_reads_covering_certain_breakpoint_regions\"\n\tis_a_homozygous_Event_outcome = " << is_a_homozygous_Event_outcome
// 	    << "\n\tis_a_homozygous_breakpoint = " << is_a_homozygous_breakpoint
// 	    << "\n\tRD_zygous_multiplier = " << RD_breakpoint_zygous_multiplier
// 	    << "\n\thaploid_upper_bound = " << haploid_upper_bound
// 	    << "\n\thaploid_coverage = " << haploid_coverage 
// 	    << "\n\n";
	    
    for (uint hap=0; hap < haploid_upper_bound; ++hap)
    {        	
	
	const type_haploid_outcome the_hap_outcome = pair_at<type_haploid_outcome>(some_result.the_diploid_Event_outcome, hap); // readability
	
	if (the_hap_outcome == hap_outcome__None)
	{  continue;  }
	
	if (my_MPI_rank == 0)
	{
	    print_line_of_markers("==");
	    std::cerr << " \n\n\n\n\n\noutput_reads_covering_certain_breakpoint_regions:  CC " << some_result.cc_of_event << ",   Event   " << some_result.event_UID 
		    << ",  outcome = " << convert_haploid_outcome_to_string(the_hap_outcome) << ",  profile breakpoint  " << pair_at<BOOST_Interval>(some_result.the_diploid_profile_brkpts,hap) << "\n\n";
	}
	
	
	
	
	
	read_Event_ALL_data_from_file(*this,  Ev_it->second);
	Ev_it->second.remaining_coordinates_for_consideration.clear();
	
	
	
	const bool is_NAHR = test_if_haploid_outcome_is_NAHR(the_hap_outcome);
	
	const BOOST_Interval  desired_breakpoint(pair_at<BOOST_Interval>(some_result.the_diploid_profile_brkpts,hap));
	
	if( !((is_NAHR and  BOOST_is_a_point(desired_breakpoint))  or  (!is_NAHR  and !BOOST_is_a_point(desired_breakpoint))) )
	{
	    std::stringstream error_strm;
	    error_strm << "ERROR!  in   \"output_reads_covering_certain_breakpoint_regions\"!  Illegal outcome/breakpoint combination:\n"
			<< "is_NAHR = " << is_NAHR << "desired_breakpoint = " << desired_breakpoint << "\n\n";
	    error_message(error_strm,false);
	    return;
	}
	
	
	
	BOOST_Interval padded_breakpoint_region__by_LCR[2];
	BOOST_Interval  breakpoint_region__by_LCR[2];
	{
	    const type_BI__BI  absolute_breakpoints_across_both_LCRs(
					convert_each_profile_breakpoint_to_absolute_breakpoints(desired_breakpoint, 
												Ev_it->second.compressed_map_LCR_to_profile__for_each_LCR));	    
													    
	    BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(absolute_breakpoints_across_both_LCRs, breakpoint_region__by_LCR);
	}
	
	std::cerr << "\n\n\t\tbreakpoint_region__by_LCR = " << breakpoint_region__by_LCR[0] << "   and   " << breakpoint_region__by_LCR[1] << "\n\n";
	    
	
	
	type_vector_int  padded_LCR_around_brkpts__vector[2];        
	for (uint qs = 0; qs < 2; ++qs)
	{                                                                                                                 
	    padded_breakpoint_region__by_LCR[qs].set( safe_subtract_base_1(breakpoint_region__by_LCR[qs].lower(),
									maximum_average_fragment_length*2 + maximum_frag_length_absolute_dev*7),
						    safe_chromo_add_base_1(breakpoint_region__by_LCR[qs].upper(),
									    maximum_average_fragment_length*2 + maximum_frag_length_absolute_dev*7,
									    Ev_it->second.chromos[0]));
				
	    padded_LCR_around_brkpts__vector[qs].clear();
	    padded_LCR_around_brkpts__vector[qs].reserve(BOOST_width_inclusive(padded_breakpoint_region__by_LCR[qs]));         
	    for (int pos_ctr = (int)padded_breakpoint_region__by_LCR[qs].lower();  pos_ctr <= (int)padded_breakpoint_region__by_LCR[qs].upper(); ++pos_ctr)
	    {  padded_LCR_around_brkpts__vector[qs].push_back(pos_ctr);  }
	}//qs    
	
	
	
	
	    
	
	for (uint lu=0; lu<2; ++lu)
	{		                                          									    
	    const uint lower_bound_LCR =  safe_subtract_base_1(padded_breakpoint_region__by_LCR[lu].lower(),
							       maximum_average_fragment_length + maximum_frag_length_absolute_dev*3);
	    const uint upper_bound_LCR =  safe_chromo_add_base_1(padded_breakpoint_region__by_LCR[lu].upper(),
								  maximum_average_fragment_length + maximum_frag_length_absolute_dev*3,
														Ev_it->second.chromos[0]);
					
	    const BOOST_Interval  selected_LCR_coordinates_for_consideration(lower_bound_LCR, upper_bound_LCR);
				    
	    type_set_uint::iterator insert_it_rem = Ev_it->second.remaining_coordinates_for_consideration.begin();        
	    for (uint LCR_pos = selected_LCR_coordinates_for_consideration.lower(); LCR_pos <= selected_LCR_coordinates_for_consideration.upper();  ++LCR_pos)
	    {  insert_it_rem = Ev_it->second.remaining_coordinates_for_consideration.insert(insert_it_rem,  LCR_pos);  }        
	}
	

	{// homol regs restriction
	    type_list_Coords_and_homologous_region::iterator it_homol =  Ev_it->second.regions_homologous_to_directly_affected_region.begin();
	    while (it_homol != Ev_it->second.regions_homologous_to_directly_affected_region.end())
	    {        
		it_homol->compressed_map_homologous_region_to_profile =   
		    get_compressed_map_inverse_of_compressed_map(
			convert_full_map_to_compressed_map(
			    get_map_intersection_of_map_keys_and_set<uint,uint>(
				    convert_compressed_map_to_full_map( 
						get_compressed_map_inverse_of_compressed_map(it_homol->compressed_map_homologous_region_to_profile)),
				    Ev_it->second.remaining_coordinates_for_consideration)));
				    
		if (it_homol->compressed_map_homologous_region_to_profile.at(0).empty())
		    it_homol = Ev_it->second.regions_homologous_to_directly_affected_region.erase(it_homol);  
		else            
		    ++it_homol;               		    
	    }// homol regs restriction
	}// homol regs restriction
	
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	std::cerr<< "\n\n\nsearching   selected_coordinates_for_consideration....\n\n";
	    
	
	
	
	{//search
	    const type_map_uint_to_list_BI search_intervals(Ev_it->second.create_search_intervals_from_DAR_and_homologous_regs_for_uploading_PERs(
								maximum_average_fragment_length + padding_for_region_homologous_to_a_var_pos, search_block_size));    
		    
							    std::cerr << "\n upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals...\n";
	    const bool success_upload = Ev_it->second.upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals_____MPI_and_OMP
							(search_intervals, 20, my_BAM_reader);
							
	    if (!success_upload)                                            
	    {  error_message("failed to upload correctly in \"output reads\"", true);  }
	}//search
				    
				
				
				    
	Ev_it->second.pre_upload_homologous_regions_from_the_Reference_genome();		    
	
			std::cerr << "\n pseudo-validation:  number PERs = " <<  Ev_it->second.PERs_on_this_profile.size()
				<< "    for_all_PERs____set_region_sequences___and___orient...\n";
	Ev_it->second.for_all_PERs____set_region_sequences___and___orient___MPI_and_OMP();
	
	
	
	
	
	
	
	
	
	const uint number_of_PERs = Ev_it->second.PERs_on_this_profile.size();
	const std::vector<type_map_string_to_PER::iterator>  PERs___loop = get_loopable_iterators<type_map_string_to_PER>(Ev_it->second.PERs_on_this_profile);
	    	    
	
	
	
	
	if (my_MPI_rank > 0)
	{  continue;  }
	
	
	
	
	std::cerr << "\n\n\n preparing Star alignments...\n\n\n";
	
	
	
	
	std::cerr << "\n\npreparing Star MRs:    LCRs 0, 1...\n\n";
	
	
	MiniRegion mini_0( Ev_it->second.chromos[0],
				    padded_breakpoint_region__by_LCR[0],
				    true,
				    false);
	    
	MiniRegion mini_1(Ev_it->second.chromos[1],
				    padded_breakpoint_region__by_LCR[1],
				    (Ev_it->second.recomb_type == recomb_class__DupDel),
				    !(Ev_it->second.recomb_type == recomb_class__DupDel));     
	    
				    
	mini_0.set_region_sequence_from_preloaded_sequence(Ev_it->second.uploaded_pieces_of_chromos);                         
	mini_1.set_region_sequence_from_preloaded_sequence(Ev_it->second.uploaded_pieces_of_chromos);    
			
	
	Ev_it->second.uploaded_pieces_of_chromos.clear();	
	
	
	
	
	
	
	
	
	
	
	
	Star_alignment  star_alignment__LCR_0;
	{//LCR 0
	    std::stringstream LCR_name;
	    LCR_name << "LCR A:  " << "chr " << Ev_it->second.chromos[0] << ":  [" << Ev_it->second.LCRs[0].lower() << " - " << Ev_it->second.LCRs[0].upper() << "]";	    
	    
	    star_alignment__LCR_0.init_Star_alignment
					    (padded_LCR_around_brkpts__vector[0],
					    mini_0.region_sequence,
					    LCR_name.str(),
					    BOOST_make_point(BOOST_width_inclusive(padded_breakpoint_region__by_LCR[0])),
					    true,   //beyond end
					    0); 
	}//LCR 0
					
	
					
	
	
	
	
	
	
	Star_alignment  star_alignment__LCR_1;
	{//LCR 1
	    std::stringstream LCR_name;
	    LCR_name << "LCR B:  " << "chr " << Ev_it->second.chromos[1] << ":  [" << Ev_it->second.LCRs[1].lower() << " - " << Ev_it->second.LCRs[1].upper() << "]"; 
	
	    std::string prepared_LCR_1_seq(mini_1.region_sequence);
	    if (Ev_it->second.recomb_type == recomb_class__Inv)
	    {
		padded_LCR_around_brkpts__vector[1] = get_reverse_of_vector<int>( padded_LCR_around_brkpts__vector[1] );        
		prepared_LCR_1_seq = get_reverse_complement_of_sequence___ie_inversion(prepared_LCR_1_seq);
	    }	
	
	    star_alignment__LCR_1.init_Star_alignment
					    (padded_LCR_around_brkpts__vector[1],
					    prepared_LCR_1_seq,
					    LCR_name.str(),
					    BOOST_make_point(BOOST_width_inclusive(padded_breakpoint_region__by_LCR[1])),
					    true,//beyond end  
					    1);
	}//LCR 1
	
                            
	
	
	
	
	
	
	
	
// 	//hybrids:
// 	Star_alignment  star_alignment__hybrid_AB_ABA;  // or "ABA"
// 	Star_alignment  star_alignment__hybrid_BA_BAB;  // or "BAB"
	
	
	
	
	const std::string  hybrid_suffix__AB__ABA =   is_NAHR  ?   "AB" : "ABA";
	const std::string  hybrid_suffix__BA__BAB =   is_NAHR  ?   "BA" : "BAB";
	
	Star_alignment  star_alignment__hybrid_AB_ABA;
	Star_alignment  star_alignment__hybrid_BA_BAB;	
	{//hybrids
			    
	    const int nahr_vs_Geneconv__adder =   is_NAHR   ?   0  :  2;
	    
	    std::cerr << "\n\npreparing Star MRs:    hybrid AB...\n\n";
	    
	    if (test_if_outcome_should_display__AB_ABA(the_hap_outcome))
	    {// AB/ABA
		const type_string__BI hybrid_AB__ABA(  
				get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
						(mini_0,
						breakpoint_region__by_LCR[0], //lower_brkpt_on_LCR[0],  //not to be included!!
						false,
						false,   //  0 = lower() is boundary,  1= upper() is boundary
						mini_1,
						breakpoint_region__by_LCR[1], // included_on_LCR_1_for_AB, //this one is included!!
						!(Ev_it->second.recomb_type == recomb_class__DupDel),                                     
						!(Ev_it->second.recomb_type == recomb_class__DupDel),
						is_NAHR));  //  0 = lower() is boundary,  1= upper() is boundary
		
	
		
		const type_vector_int hybrid_int_vector_AB__ABA(	
					    label_absolute_coordinates_on_hybrid
							(mini_0,
							breakpoint_region__by_LCR[0],
							mini_1,
							breakpoint_region__by_LCR[1],
							is_NAHR,
							false,
							Ev_it->second.recomb_type));
		
		std::stringstream hybrid_name;
		hybrid_name << "Hybrid " << hybrid_suffix__AB__ABA << ":  brkpt  " << desired_breakpoint;
		    
		star_alignment__hybrid_AB_ABA.init_Star_alignment
						(hybrid_int_vector_AB__ABA,
						hybrid_AB__ABA.first,
						hybrid_name.str(),
						hybrid_AB__ABA.second,
						is_NAHR,
						2 + nahr_vs_Geneconv__adder);
	    }// AB/ABA	
	
	
	
	
	    std::cerr << "\n\npreparing Star MRs:    hybrid BA...\n\n";
	    if(test_if_outcome_should_display__BA_BAB(the_hap_outcome))
	    {// BA/BAB  
		const type_string__BI hybrid_BA__BAB(  
				get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
						(mini_1,
						breakpoint_region__by_LCR[1] , //included_on_LCR_1___BA,    //not to be included!!
						!(Ev_it->second.recomb_type == recomb_class__DupDel),                                     
						!(Ev_it->second.recomb_type == recomb_class__DupDel),         //  0 = lower() is boundary,  1= upper() is boundary    
						mini_0,
						breakpoint_region__by_LCR[0],//lower_brkpt_on_LCR[0],  
						false,             //always false in this case since this is LCR  0  !!!!!!
						false,
						is_NAHR));   //  0 = lower() is boundary,  1= upper() is boundary  		
		
		const type_vector_int hybrid_int_vector_BA__BAB(
					    label_absolute_coordinates_on_hybrid
							(mini_0,
							breakpoint_region__by_LCR[0],
							mini_1,
							breakpoint_region__by_LCR[1],
							is_NAHR,  //else, GeneConv
							true,
							Ev_it->second.recomb_type));
						    
		std::stringstream hybrid_name;
		hybrid_name << "Hybrid " << hybrid_suffix__BA__BAB << ":  brkpt  " << desired_breakpoint;         
		
		star_alignment__hybrid_BA_BAB.init_Star_alignment
						    (hybrid_int_vector_BA__BAB,
						    hybrid_BA__BAB.first,
						    hybrid_name.str(),
						    hybrid_BA__BAB.second,
						    is_NAHR,
						    3 + nahr_vs_Geneconv__adder);
	    }// BA/BAB
	}//hybrids
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	std::cerr << "\nabout to do alignments...\n";	    	
	//NO GPU!!!!
	#pragma omp parallel
	{//parallel
	
	    #pragma omp for schedule(dynamic,3)
	    for (uint j=0; j<number_of_PERs; ++j)    
	    {      
		{
		    const type_map_uint_to_uint_to_set_uint MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY(
							    PERs___loop[j]->second.create_MiniEvents_and_identify_gender_of_MRs(&Ev_it->second));

		    PERs___loop[j]->second.make_partners_out_of_MiniRegions_for_each_MiniEvent
								    (MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY);
		}   
		if (!display_ALL_reads_not_just_good_ones)
		{
			const bool return_alignment_val = 
			PERs___loop[j]->second.align_to_every_possible_hybrid_outcome_AND_DISPLAY(
								some_result.event_UID, desired_breakpoint, the_hap_outcome,
								star_alignment__LCR_0, star_alignment__LCR_1,
								star_alignment__hybrid_AB_ABA,  star_alignment__hybrid_BA_BAB);	
		}
		else
		{
			const bool return_alignment_val = 
			PERs___loop[j]->second.align_to_every_possible_hybrid_outcome_AND_DISPLAY_REGARDLESS_OF_NOHYBRIDS(
								some_result.event_UID, desired_breakpoint, the_hap_outcome,
								star_alignment__LCR_0, star_alignment__LCR_1,
								star_alignment__hybrid_AB_ABA,  star_alignment__hybrid_BA_BAB);				
			
		}
	    }//j
	    
	    mpfr_free_cache();
	}//parallel







	
	real focused_likelihood___null_hypothesis(1);
	real focused_likelihood___alternative_hypothesis(1);
	
	if (do_likelihood_ratio_test_on_alignments)
	{//likelihood ratio test
	    std::cerr << "\tdoing likelihood ratio test...\n";
	
	    
	    std::cerr << "\tnull...\n";
	    {//null hypothesis
		for (uint j=0; j<number_of_PERs; ++j)    
		{
		    const type_uint__real haploid_0__sum_P_mPER(
					PERs___loop[j]->second.calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints
									    (false,
									    Sparse_map(),
									    type_map_uint_to_BI()));
									    
		    const type_uint__real haploid_1__sum_P_mPER(
					PERs___loop[j]->second.calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints
									    (true,
									    Sparse_map(),
									    type_map_uint_to_BI()));
		    
		    const uint total_diploid_contributions = haploid_0__sum_P_mPER.first + haploid_1__sum_P_mPER.first;
		    
		    switch (total_diploid_contributions)
		    {
			case 0: 
			    break;
			case 1:
			    focused_likelihood___null_hypothesis *= (haploid_0__sum_P_mPER.second + haploid_1__sum_P_mPER.second);
			    break;                                            
			default:
			    focused_likelihood___null_hypothesis *= (haploid_0__sum_P_mPER.second + haploid_1__sum_P_mPER.second);
			    focused_likelihood___null_hypothesis /=  total_diploid_contributions; 
			    break;                                                                        
		    }   
		}//j	    	    	    
	    }//null hypothesis
	    
	    std::cerr << "\talternative...\n";
	    {//alternative
		Sparse_map dummy_states[2];
		type_map_uint_to_BI dummy_breakpoints[2];
		
		for (uint state_hap = 0; state_hap < 2; ++state_hap)
		{
		    if (pair_at<BOOST_Interval>(some_result.the_diploid_Event_outcome,state_hap) != hap_outcome__None)
		    {
			dummy_states[state_hap].sparse_map_UID_to_value[some_result.event_UID] = pair_at<type_haploid_outcome>(some_result.the_diploid_Event_outcome,state_hap);
			dummy_breakpoints[state_hap][some_result.event_UID] = pair_at<BOOST_Interval>(some_result.the_diploid_profile_brkpts,state_hap);
		    }
		}    	    

								
		for (uint j=0; j<number_of_PERs; ++j)    
		{
		    const type_uint__real haploid_0__sum_P_mPER(
					PERs___loop[j]->second.calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints
									    (false,
									    dummy_states[0],
									    dummy_breakpoints[0]));
									    
		    const type_uint__real haploid_1__sum_P_mPER(
					PERs___loop[j]->second.calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints
									    (true,
									    dummy_states[1],
									    dummy_breakpoints[1]));
		    
		    const uint total_diploid_contributions =  haploid_0__sum_P_mPER.first + haploid_1__sum_P_mPER.first;
		    
		    switch (total_diploid_contributions)
		    {
			case 0: 
			    break;
			case 1:
			    focused_likelihood___alternative_hypothesis *= (haploid_0__sum_P_mPER.second + haploid_1__sum_P_mPER.second);
			    break;                                            
			default:
			    focused_likelihood___alternative_hypothesis *= (haploid_0__sum_P_mPER.second + haploid_1__sum_P_mPER.second);
			    focused_likelihood___alternative_hypothesis /=  total_diploid_contributions; 
			    break;                                                                        
		    }   
		}//j
	    }//alternative
	    
	}//likelihood ratio test
	
	









	
	
	Ev_it->second.PERs_on_this_profile.clear();                        
			
	Ev_it->second.uploaded_pieces_of_chromos.clear();
	
	
	
	
	
	
	
	
	
	
	
	std::cerr << "\nabout to do save star alignments...\n";		
	{//save star alignments   
			
	    char star_algmts_filename[option_size];        
	    std::sprintf(star_algmts_filename, "%s/Star_alignments/star_alignments___cc_%u__ev_%u__%s__breakpoint_%u_%u___%s__%s.html",
					    output_dir.c_str(),				
					    Conn_Comp_ID, some_result.event_UID,
					    convert_haploid_outcome_to_string(the_hap_outcome).c_str(),
					    desired_breakpoint.lower(), desired_breakpoint.upper(),
					    genome_name.c_str(),
					    population_abbreviation.c_str());       
				
	    std::ofstream  save_star_algmts_filestream(star_algmts_filename);        
					if (  !save_star_algmts_filestream.is_open())
					{
					    std::fprintf(stderr, "\n\n\n\n\nERROR: Command [%s] could not be opened for \"save group of star_alignments_to_file_as_html\".\n\n\n\n\n\n", 
								star_algmts_filename);
					    continue;
					}       
					
					
	    save_star_algmts_filestream << "\n<pre>\n<font size=1>\n\n";
	    save_star_algmts_filestream << "\n<span style=\"font-size: 18pt\">Genome: " <<  genome_name 
					<< "\nConnected component: " <<  Conn_Comp_ID
					<< "\nEvent: " << some_result.event_UID 
					<< "\nbreakpoint: " << desired_breakpoint
					<< "\ncall: " << convert_haploid_outcome_to_string(the_hap_outcome)
					<< "</span>\n\n";                                    
					
					
	    
	    {// LCR A	 
		std::cerr << "saving LCR A...\n";
		save_star_algmts_filestream << "\n\n<span style=\"font-size: 16pt\">LCR  A</span>\n\n";
		    
		type_vector_string  html_output;
		const bool success_create_html =  star_alignment__LCR_0.create_HTML_output_for_star_alignment(   
									    (int)Ev_it->second.chromos[0], -1,
									    true,
									    html_output);            
										    
		for (uint row=0; row < html_output.size(); ++row)
		{  save_star_algmts_filestream << html_output.at(row) << "\n";  }                            
	    
		
		star_alignment__LCR_0.clear_all_data(); 	    
	    }// LCR A
	    
	    
	    
	    
	   
	    {// LCR B	    
		std::cerr << "saving LCR B...\n";
		save_star_algmts_filestream << "\n\n<span style=\"font-size: 16pt\">LCR  B</span>\n\n";
		    
		type_vector_string  html_output;
		const bool success_create_html =  star_alignment__LCR_1.create_HTML_output_for_star_alignment(   
									    (int)Ev_it->second.chromos[1], -1,
									    true,
									    html_output);            
										    
		for (uint row=0; row < html_output.size(); ++row)
		    save_star_algmts_filestream << html_output.at(row) << "\n";                                                           
	    
		
		star_alignment__LCR_1.clear_all_data();	    
	    }// LCR B        
					
					
					
					
			
	    
	    if (test_if_outcome_should_display__AB_ABA(the_hap_outcome))
	    {// Hybrid AB	    
		std::cerr << "saving Hybrid AB/ABA...\n";
	    
		const type_longdouble__longdouble expected_read_count____informative___GConly(
				    Ev_it->second.calculate_expected_number_of_PERs_at_breakpoint(
							desired_breakpoint,
							mini_0,
							mini_1,
							false));
	    
		save_star_algmts_filestream << "\n\n<span style=\"font-size: 16pt\">Hybrid  " << hybrid_suffix__AB__ABA  << "</span>\n";
		save_star_algmts_filestream
			<< "<span style=\"font-size: 12pt\">naive expected number of observed reads overlapping breakpoint:  " << haploid_coverage * RD_breakpoint_zygous_multiplier
			<< "\nGC-bias adjusted expected number of observed reads overlapping breakpoint:  " << expected_read_count____informative___GConly.second * RD_breakpoint_zygous_multiplier
			<< "\nGC-bias adjusted expected number of INFORMATIVE reads at breakpoint:  " << expected_read_count____informative___GConly.first * RD_breakpoint_zygous_multiplier
			<<  "</span>\n\n";
		    
		type_vector_string  html_output;
		const bool success_create_html =  star_alignment__hybrid_AB_ABA.create_HTML_output_for_star_alignment(   
									    (int)Ev_it->second.chromos[0],
									    (int)Ev_it->second.chromos[1],
									    true,
									    html_output    );
		
										    
		for (uint row=0; row < html_output.size(); ++row)
		{  save_star_algmts_filestream << html_output.at(row) << "\n";  }
	    
		
		star_alignment__hybrid_AB_ABA.clear_all_data();	    
	    }// Hybrid AB                                      
			    
			    
			    
			    
	    if (test_if_outcome_should_display__BA_BAB(the_hap_outcome))
	    {// Hybrid BA	
		std::cerr << "saving Hybrid BA/BAB...\n";
	    
		const type_longdouble__longdouble expected_read_count____informative___GConly(
				    Ev_it->second.calculate_expected_number_of_PERs_at_breakpoint(desired_breakpoint,
									      mini_0,
									      mini_1,
									      true));
	    
		save_star_algmts_filestream << "\n\n<span style=\"font-size: 16pt\">Hybrid  " << hybrid_suffix__BA__BAB  << "</span>\n";
		save_star_algmts_filestream
			<< "<span style=\"font-size: 10pt\">naive expected number of observed reads overlapping breakpoint:  " << haploid_coverage * RD_breakpoint_zygous_multiplier
			<< "\nGC-bias adjusted expected number of observed reads overlapping breakpoint:  " << expected_read_count____informative___GConly.second * RD_breakpoint_zygous_multiplier
			<< "\nGC-bias adjusted expected number of INFORMATIVE reads at breakpoint:  " << expected_read_count____informative___GConly.first * RD_breakpoint_zygous_multiplier
			<<  "</span>\n\n";	    
	    	    
		    
		type_vector_string  html_output;
		const bool success_create_html =  star_alignment__hybrid_BA_BAB.create_HTML_output_for_star_alignment(   
									    (int)Ev_it->second.chromos[1],
									    (int)Ev_it->second.chromos[0],
									    true,
									    html_output);
		
										    
		for (uint row=0; row < html_output.size(); ++row)
		{  save_star_algmts_filestream << html_output.at(row) << "\n";   }                                                        
	    
		
		star_alignment__hybrid_BA_BAB.clear_all_data();
	    
	    }// Hybrid BA                           
		
		
		
		
		
		
	    if (do_likelihood_ratio_test_on_alignments)
	    {//likelihood ratio test
		save_star_algmts_filestream
		    << "\n\n<span style=\"font-size: 8pt\">"
		    << "P(data | null hypothesis)  =  " <<  focused_likelihood___null_hypothesis.toString(15)
		    << "\nP(data | alternative hypothesis)  =  " <<  focused_likelihood___alternative_hypothesis.toString(15)
		    << "\nodds ratio:   P(data | null) / P(data | alternative)  =  "
		    << (focused_likelihood___null_hypothesis / focused_likelihood___alternative_hypothesis).toString(15)
		    << "\nlog-odds:  " 
		    << (mpfr::log10(focused_likelihood___null_hypothesis) - mpfr::log10(focused_likelihood___alternative_hypothesis)).toString(15)
		    << "\n\n</span>";    
	    }//likelihood ratio test



	    save_star_algmts_filestream << "\n\n</font>\n</pre>\n";
	    save_star_algmts_filestream.close();               

	}//save star alignments    
	
	
	
	
	if (  (test_if_outcome_should_display__AB_ABA(the_hap_outcome) and star_alignment__hybrid_AB_ABA.failed_alignment)
		or   (test_if_outcome_should_display__BA_BAB(the_hap_outcome)  and star_alignment__hybrid_BA_BAB.failed_alignment)
		or  star_alignment__LCR_0.failed_alignment  or star_alignment__LCR_1.failed_alignment)
	{
	    std::cerr << "\n\n\n\nWARNING!!! THERE WAS A FAILED STAR ALIGNMENT IN \"output_reads_covering_certain_breakpoint_regions\" !!!!\n\n"
		    << "star_alignment__hybrid_AB_ABA.failed_alignment = " << star_alignment__hybrid_AB_ABA.failed_alignment << "\n"
		    << "star_alignment__hybrid_BA_BAB.failed_alignment = " << star_alignment__hybrid_AB_ABA.failed_alignment << "\n"
		    << "star_alignment__LCR_0.failed_alignment = " << star_alignment__hybrid_AB_ABA.failed_alignment << "\n"
		    << "star_alignment__LCR_1.failed_alignment = " << star_alignment__hybrid_AB_ABA.failed_alignment << "\n";		    
		    
	    Ev_it->second.print_this_entry();	
	}
	
	
				
	
	if (my_MPI_rank == 0)
	{  print_line_of_markers("==");  }	
		        
    }//hap
    
} //  output_reads_covering_certain_breakpoint_regions


































void  Conn_Comp::create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights
				(BamTools::BamReader &my_BAM_reader,
				 const Sampled_diploid_Event_data &some_result,
				 const std::string &outdir)
{
    
    print_line_of_markers("==");
    std::cerr << " \n\n\n\n\n\ncreate_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights:   Event   " << some_result.event_UID 
	<< "\n\n";
			
	
    some_result.print_this_Sampled_diploid_Event_data();
    
    
    const type_map_uint_to_Event::iterator result_Ev_it = events.find(some_result.event_UID); 
    
    
    if (result_Ev_it->second.chromos[0] == 24  and   gender_of_individual_is_female)
    {
	std::cerr << "Event is on Y chromosome, but individual is female.  Skipping read-depth graph in in \"create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights\".\n";
	return;
    }
    
    const BOOST_Interval expanded_region_of_interest( 
			    safe_subtract_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.lower(), 100000),
			    safe_chromo_add_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.upper(), 100000, result_Ev_it->second.chromos[0]));        
			    
    //identify all events intersecting this region:
    type_map_uint_to_Event_ptr  all_events_near_expanded_region;
    for (type_map_uint_to_Event::iterator it_ev = events.begin();
	    it_ev != events.end();
	    ++it_ev)
    {	    
	if (it_ev->second.chromos[0] == result_Ev_it->second.chromos[0]
	    and  (BOOST_overlap(it_ev->second.LCRs[0],expanded_region_of_interest)
		    or  BOOST_overlap(it_ev->second.LCRs[1],expanded_region_of_interest)))
	{
	    all_events_near_expanded_region[it_ev->first] = &it_ev->second;
	}
    }//ev
    
    print_map_keys<uint, Event*>(all_events_near_expanded_region, "all_events_near_expanded_region", NULL, true);
    
    
    
    
    
    
    
      
    
    type_map_string_to_PER  grand_collection_of_PERs;
    type_map_uint_to_list_BI  all_search_intervals;    
    
    for (type_map_uint_to_Event_ptr::iterator it_ev = all_events_near_expanded_region.begin();
	    it_ev != all_events_near_expanded_region.end();
	    ++it_ev)
    {
	//we are not going to test for any hybrid read business.
	it_ev->second->variational_positions_of_profile.clear();
	it_ev->second->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.clear();
	it_ev->second->my_original_Breakpoint_complex.clear_all_data();
		
	const type_map_uint_to_list_BI search_intervals(it_ev->second->create_search_intervals_from_DAR_and_homologous_regs_for_uploading_PERs
							(maximum_average_fragment_length + padding_for_region_homologous_to_a_var_pos, search_block_size) );    
		
							std::cerr << "\n upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals...\n";
	const bool success_upload = it_ev->second->upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals_____MPI_and_OMP
						    (search_intervals, 20, my_BAM_reader);						      
						    
	for (type_map_uint_to_list_BI::const_iterator it_sub = search_intervals.begin();
		it_sub != search_intervals.end();
		++it_sub)
	{
	    for (type_list_BI::const_iterator it_bi = it_sub->second.begin();
		    it_bi != it_sub->second.end();
		    ++it_bi)
	    {  all_search_intervals[it_sub->first].push_back(*it_bi);  } 
	}
						    
	if (!success_upload)                                            
	{  error_message("failed to upload correctly in \"inferred RD\"", true);  }
	
	it_ev->second->pre_upload_homologous_regions_from_the_Reference_genome();
	
			std::cerr << "\n pseudo-RD:  number PERs = " <<  it_ev->second->PERs_on_this_profile.size()
				<< "    for_all_PERs____set_region_sequences___and___orient...\n";
	it_ev->second->for_all_PERs____set_region_sequences___and___orient___MPI_and_OMP();  	
	
	grand_collection_of_PERs.insert(it_ev->second->PERs_on_this_profile.begin(), it_ev->second->PERs_on_this_profile.end());
	it_ev->second->PERs_on_this_profile.clear();
    }//ev

		

    
	

    std::cerr << "\nabout to do alignments...\n";	 
    
    set_up_and_align_PERs_using_GPU____OMP
			(grand_collection_of_PERs,
			 &result_Ev_it->second,
			 false,//create_and_consider_MiniEvents,
			 false,//consider_GeneConversion_breakpoints,  
			 false,//sum_over_all_unaffected_LCRs,  //"false" for heuristic
			 false,//special_save_breakpoints_actually_captured_in_alignment); //"true" for heuristic    
			false);//do_PER_orient_formatting
    
    
    
    std::cerr<< "\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "  done uploading and preparing.  Beginning in \"pseudovalidation\" broadcasting...\n";
    
    
    for (type_map_uint_to_Event_ptr::iterator it_ev = all_events_near_expanded_region.begin();
	    it_ev != all_events_near_expanded_region.end();
	    ++it_ev)
    {
	it_ev->second->uploaded_pieces_of_chromos.clear();    
    }
    
        
    
    
    
    const uint number_of_PERs = grand_collection_of_PERs.size();
    std::vector< type_map_string_to_PER::iterator >  PERs___loop = get_loopable_iterators<type_map_string_to_PER>(grand_collection_of_PERs);                    		
    
    
    {//broadcasting    
	world_MPI_boost_ptr->barrier();  //NECESSARY!
	
	//MPI
	for (int rank=0; rank < size_MPI_WORLD; ++rank) 
	{
	    if (my_MPI_rank == rank)
	    {
		type_map_string_to_PER sending_PERs;
		type_map_string_to_PER::iterator insert_it = sending_PERs.begin();
		for (int p = my_MPI_rank; p < number_of_PERs; p += size_MPI_WORLD)
		{            
		    insert_it = sending_PERs.insert(insert_it, *(PERs___loop[p])  );
		    insert_it->second.deep_copy_from_serialize( PERs___loop[p]->second, this  );
		}
		
		boost::mpi::broadcast<type_map_string_to_PER>(*world_MPI_boost_ptr, sending_PERs, rank);                                            
	    } // sending rank
	    else
	    {
		type_map_string_to_PER received_PERs;
		boost::mpi::broadcast<type_map_string_to_PER>(*world_MPI_boost_ptr, received_PERs, rank);
		
		type_map_string_to_PER::iterator it_recv = received_PERs.begin();
		while (it_recv   !=   received_PERs.end())
		{
		    const type_map_string_to_PER::iterator found_it = grand_collection_of_PERs.find( it_recv->first );                
		    found_it->second.deep_copy_from_serialize(it_recv->second, this);
		    
		    received_PERs.erase(it_recv++);   //post-increment necessary.
		}                            
	    } // receiving rank
	}//rank
	
	
	world_MPI_boost_ptr->barrier();  //NECESSARY!    
    }//broadcasting
    
	
	
	
	std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "    counting reads... ";
	
    type_map_uint_to_uint counts_per_position_of_expanded_region__NOT_from_uploaded_PER_regions;    
    {//counts
	initialize_map_across_interval__via_insert<uint,uint>(counts_per_position_of_expanded_region__NOT_from_uploaded_PER_regions,
							      expanded_region_of_interest, 0);  
							      
			
	my_BAM_reader.SetRegion(create_BAM_region(result_Ev_it->second.chromos[0], expanded_region_of_interest, 3));	
	count_decent_reads_with_left_endpoints_in_set_region(
						    my_BAM_reader,
						    expanded_region_of_interest,
						    counts_per_position_of_expanded_region__NOT_from_uploaded_PER_regions);	
    
	//screen    
	const type_map_uint_to_list_BI::const_iterator it_found_chr = all_search_intervals.find(result_Ev_it->second.chromos[0]);    
	if (it_found_chr != all_search_intervals.end())
	{
	    for (type_list_BI::const_iterator it_bi = it_found_chr->second.begin();
		    it_bi != it_found_chr->second.end();
		    ++it_bi)
	    {
		const BOOST_Interval searched_part_of_expanded_reg(BOOST_intersect(expanded_region_of_interest, *it_bi));
		for (uint j=searched_part_of_expanded_reg.lower(); j <= searched_part_of_expanded_reg.upper(); ++j)
		{  counts_per_position_of_expanded_region__NOT_from_uploaded_PER_regions.at(j) = 0;  }	    
	    }//bi    
	}//chr
    }//counts
    
    std::cerr << "\n\n\t\tTotal number of counts in flanking non-repeat region =  "
		<< sum_over_map_values<uint,uint>(counts_per_position_of_expanded_region__NOT_from_uploaded_PER_regions) << "\n\n";
		
		
	
    
    
    
    
    
    if (my_MPI_rank > 0)
    {  return;  }

	


	
	

    
    //Doing an RD hypothesis test on the region of the called Event, with or without inffered Read Depth, is easy.  The above code sets you up for that.
    //But we also like to display a large amount of flanking on either side of the events for effect/context.  This is trickier, because 100,000 of flanking is very likely to include parts of other Connected Components/Events.  Some of these other CCs or EVs may even be huge/uncomputable!!!
    //WE don't care too much about the flanking - it is just a visual.  The test is in the called region.  Thus, for the called region, ww use the RD inferred from the uploaded reads and alignments.  But for the flanking, we don't upload or any any reads.  Instead, we just get counts-per-position.  So as not to interfere with the Rd inferred rigorously from the uploaded reads, we mask this "counting approximation" by the intervals from which we uploadeed the reads (set the counts at these positions to be 0.  Thus, we synthesize the counts-per-position and the infferred-from-read-alignment read-depths to get a big picture.
    

    
	

    
    
    //get "inferred" observed values.
    type_map_uint_to_longdouble map_Ref_Genome_position_to_per_read_posterior_sum;
    
    #pragma omp parallel
    {//parallel	
	type_map_uint_to_longdouble map_Ref_Genome_position_to_per_read_posterior_sum__THREAD_specific;
    
	#pragma omp for schedule(dynamic,100)
	for (uint j=0; j<number_of_PERs; ++j)    
	{
	    type_map_uint_to_uint_to_real  map_chromo_to_absolute_position_to_likelihood;
	    
	    for (type_map_uint_to_MiniRegion::const_iterator it_MR = PERs___loop[j]->second.my_MiniRegions.begin();
		    it_MR != PERs___loop[j]->second.my_MiniRegions.end();
		    ++it_MR)
	    {
		map_chromo_to_absolute_position_to_likelihood
			[it_MR->second.chromosome_of_region]
			    [PERs___loop[j]->second.map_MiniID__to__nohybrid_almgt_interval_absolute.at(it_MR->first).lower()]
				= PERs___loop[j]->second.map_affected_MiniID__to__P_nohybrid.at(it_MR->first);
	    }//MR
	    
	    //prepare for sum
	    type_list_real  likelihoods_in_bayes_denominator;
	    for (type_map_uint_to_uint_to_real::const_iterator it_chr = map_chromo_to_absolute_position_to_likelihood.begin();
		    it_chr != map_chromo_to_absolute_position_to_likelihood.end();
		    ++it_chr)
	    {
		for (type_map_uint_to_real::const_iterator it_pos = it_chr->second.begin();
			it_pos != it_chr->second.end();
			++it_pos)
		{
		    likelihoods_in_bayes_denominator.push_back(it_pos->second);
		}//pos
	    }//chr	    	    
	    
	    likelihoods_in_bayes_denominator.sort();
	    const real bayes_denominator(sum_over_container<type_list_real, real>(likelihoods_in_bayes_denominator));
		//no need to multiply by the "prior" = 1/likelihoods_in_bayes_denominator.size()    because  it will cancel with numerator.;
	    
	    //distribute posterior to genome positions:
	    const type_map_uint_to_uint_to_real::const_iterator it_find_chr = map_chromo_to_absolute_position_to_likelihood.find(result_Ev_it->second.chromos[0]);
	    if (it_find_chr != map_chromo_to_absolute_position_to_likelihood.end())
	    {
		for (type_map_uint_to_real::const_iterator it_pos = it_find_chr->second.begin();
			it_pos != it_find_chr->second.end();
			++it_pos)
		{
		    if (BOOST_in(it_pos->first, expanded_region_of_interest))
		    {
			map_Ref_Genome_position_to_per_read_posterior_sum__THREAD_specific[it_pos->first] += (it_pos->second/bayes_denominator).toLDouble();
		    }
		}//pos
	    }//chr
	}//PERs
	
	
	#pragma omp critical (add_inferred_bayes_posterior_read_depth_counts_back_into_master_counts)
	for (type_map_uint_to_longdouble::const_iterator it_pos = map_Ref_Genome_position_to_per_read_posterior_sum__THREAD_specific.begin();
		it_pos != map_Ref_Genome_position_to_per_read_posterior_sum__THREAD_specific.end();
		++it_pos)
	{
	    map_Ref_Genome_position_to_per_read_posterior_sum[it_pos->first] += it_pos->second;	    
	}//pos	
	
	
	mpfr_free_cache();
    }//parallel

    std::cerr << "\n\n\t\tmap_Ref_Genome_position_to_per_read_posterior_sum.size() = " << map_Ref_Genome_position_to_per_read_posterior_sum.size() << "\n\n";

    
    grand_collection_of_PERs.clear();
    
    
    
    
    //combine
    for (type_map_uint_to_uint::const_iterator it_count = counts_per_position_of_expanded_region__NOT_from_uploaded_PER_regions.begin();
	    it_count != counts_per_position_of_expanded_region__NOT_from_uploaded_PER_regions.end();
	    ++it_count)
    {  map_Ref_Genome_position_to_per_read_posterior_sum[it_count->first] += it_count->second;  }//count
    
    
    
    
    
    
    std::cerr << "map_Ref_Genome_position_to_per_read_posterior_sum.size() = " << map_Ref_Genome_position_to_per_read_posterior_sum.size() << "\n\n";
    
    //get expected values for different hypotheses...
    const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble map_diploid_outcome_to_absolute_position_to_expected_RD(
				construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region(
							    result_Ev_it->second.chromos[0],
							     expanded_region_of_interest,
							    some_result.the_diploid_Event_outcome,
							    some_result.the_diploid_profile_brkpts,
							    result_Ev_it->second));    
    
    
    
    
    std::string checked_outdir(outdir);
    checked_outdir.append("/Read_Depth");
    
    
    char out_RD_core_fname[option_size];	
    std::sprintf(out_RD_core_fname, "%s/RD_hypothesis_test__cc_%u__ev_%u___%s_%s___%s", 
				checked_outdir.c_str(),
				Conn_Comp_ID, result_Ev_it->second.UID,
				convert_haploid_outcome_to_string(some_result.the_diploid_Event_outcome.first).c_str(),
				convert_haploid_outcome_to_string(some_result.the_diploid_Event_outcome.second).c_str(),
				genome_name.c_str());
				
				
    {//save
	char out_data_fname[option_size];
	std::sprintf(out_data_fname, "%s.data.observed", out_RD_core_fname);	
	write_map_to_file<uint, longdouble>(out_data_fname, map_Ref_Genome_position_to_per_read_posterior_sum);
	
	
	std::sprintf(out_data_fname, "%s.data.expected", out_RD_core_fname);
	std::ofstream out_expected_fs(out_data_fname);
	if (check_fstream_for_goodness(out_expected_fs, out_data_fname, "create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights", false))
	{
	    const type_haploid_outcome__haploid_outcome alternative_dip_state(some_result.the_diploid_Event_outcome);
	    std::cerr << "alternative_dip_state = " << convert_haploid_outcome_to_string(alternative_dip_state.first) << ", " << convert_haploid_outcome_to_string(alternative_dip_state.second) << "\n";
	    
	    const type_haploid_outcome__haploid_outcome null_dip_state(type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None));
	    
	    
	    const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_null_dip_state_expected_map 
				    = map_diploid_outcome_to_absolute_position_to_expected_RD.find(null_dip_state);
	    const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_alternative_dip_state_expected_map 
				    = map_diploid_outcome_to_absolute_position_to_expected_RD.find(alternative_dip_state);
	    

	    if (it_null_dip_state_expected_map != map_diploid_outcome_to_absolute_position_to_expected_RD.end()
		and  it_alternative_dip_state_expected_map != map_diploid_outcome_to_absolute_position_to_expected_RD.end())
	    {
		for (type_map_uint_to_longdouble::const_iterator it_null = it_null_dip_state_expected_map->second.begin(),
								it_expect = it_alternative_dip_state_expected_map->second.begin();
			it_null != map_diploid_outcome_to_absolute_position_to_expected_RD.at(null_dip_state).end();
			++it_null, ++it_expect)
		{  out_expected_fs << it_null->first << "\t" << it_null->second << "\t" << it_expect->second << "\n";  }
	    }
	    else
	    {  std::cerr << "\n\tunrecognized alternative dip sates!!!\n\n";  }
	    
	    out_expected_fs.close();	
	}
    }//save

	
	
    std::cerr << "\t\t\tidentifying nonuniqueness in genome...\n\n";
    const type_set_uint  nonunique_positions_in_region_of_interest(
					    get_nonunique_positions_in_interval(
							load_Non_unique_regions_of_genome(),
							result_Ev_it->second.chromos[0], expanded_region_of_interest));
    
    std::cerr << "nonunique_positions_in_region_of_interest.size = " << nonunique_positions_in_region_of_interest.size() << "\n";
    
    {//save
	char out_data_fname[option_size];
	std::sprintf(out_data_fname, "%s.data.nonunique", out_RD_core_fname);
	write_set_to_file<uint>(out_data_fname, nonunique_positions_in_region_of_interest);
    }//save    
    
    
    
	
	
	

	
	

    {//RD hypothesis test

	BOOST_Interval affected_region(empty_BI);
	
	if (some_result.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None))
	{
	    affected_region = result_Ev_it->second.region_between_and_including_the_LCRs_themselves;	    
	}
	else
	{
	    for (uint hap=0; hap<2; ++hap)
	    {
		if (pair_at<type_haploid_outcome>(some_result.the_diploid_Event_outcome,hap) != hap_outcome__None)
		{
		    const type_BI__BI affected_regions_each_LCR(
					    BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(
						    convert_each_profile_breakpoint_to_absolute_breakpoints(
							    pair_at<BOOST_Interval>(some_result.the_diploid_profile_brkpts,hap),
							    result_Ev_it->second.compressed_map_LCR_to_profile__for_each_LCR)));
					    
		    const BOOST_Interval affected_hull(affected_regions_each_LCR.first.lower(), affected_regions_each_LCR.second.upper());
																	    
		    if (BOOST_empty(affected_region))
		    {  affected_region = affected_hull;  }
		    else
		    {  affected_region = BOOST_hull(affected_region, affected_hull);  }
		}//occurs
	    }//hap
	}
	
	std::cerr << "\n\naffected_region = " << affected_region << "\n\n";
	
	//restrict to affected region:
	const type_map_uint_to_longdouble  map_affected_Ref_Genome_position_to_per_read_posterior_sum(
		    get_map_intersection_of_map_keys_and_interval<uint,longdouble>(map_Ref_Genome_position_to_per_read_posterior_sum, affected_region));
		    
	type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble  map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD;		    
	for (type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_full_out 
			    = map_diploid_outcome_to_absolute_position_to_expected_RD.begin();
		it_full_out != map_diploid_outcome_to_absolute_position_to_expected_RD.end();
		++it_full_out)
	{
	    map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD[it_full_out->first]
		    = get_map_intersection_of_map_keys_and_interval<uint, longdouble>(it_full_out->second, affected_region);
	}//out
	

	
	std::cerr << "\n\n\tmaking expected_RD_of_region__NULL.  map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.size() = "
						<< map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.size() << "\n\n";
	
	

						
						
	const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_expected_RD_of_region__NULL
		= map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.find(type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None));
		
	const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_expected_RD_of_region__ALTERNATIVE
		= map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.find(some_result.the_diploid_Event_outcome);
		
	const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_expected_RD_of_region__def_alt_del
		= map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.find(type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__Del));
		
	const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_expected_RD_of_region__def_alt_dup
		= map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.find(type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__Dup));			
		
		
		
	if (it_expected_RD_of_region__NULL == map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.end()
		or  it_expected_RD_of_region__ALTERNATIVE == map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.end())	
	{
	    std::stringstream errorstrm;
	    errorstrm << "ERROR!  in \"create_inferred...\"\n"
			<< "it_expected_RD_of_region__NULL == map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.end()     =    "  
			<< (it_expected_RD_of_region__NULL == map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.end())
			<< "\nit_expected_RD_of_region__ALTERNATIVE == map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.end()     =    "  
			<< (it_expected_RD_of_region__ALTERNATIVE == map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.end())
			<< "\n\n";
		
	    errorstrm << "keys of \"map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD\":\n";
	    for (type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_full_out 
				= map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.begin();
		    it_full_out != map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.end();
		    ++it_full_out)
	    {
		errorstrm << convert_haploid_outcome_to_string(it_full_out->first.first) << ",   " << it_full_out->first.second << "\n";	    	    
	    }
	    
	    errorstrm << "\n\nkeys of \"map_diploid_outcome_to_absolute_position_to_expected_RD\":\n";
	    for (type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_full_out 
				= map_diploid_outcome_to_absolute_position_to_expected_RD.begin();
		    it_full_out != map_diploid_outcome_to_absolute_position_to_expected_RD.end();
		    ++it_full_out)
	    {
		errorstrm << convert_haploid_outcome_to_string(it_full_out->first.first) << ",   " << it_full_out->first.second << "\n";	    	    
	    }	   
	    
	    error_message(errorstrm, false);
	    return;
	}
		
		
			
	std::cerr << "\n\n creating params\n\n";		
	
	const longdouble observed_counts__inferred = sum_over_map_values<uint,longdouble>(map_affected_Ref_Genome_position_to_per_read_posterior_sum);
	const longdouble observed_counts__UNIQUE_ONLY
		= sum_over_map_values<uint,longdouble>(
		    erase_set_from_map_keys<uint,longdouble>(map_affected_Ref_Genome_position_to_per_read_posterior_sum, nonunique_positions_in_region_of_interest));

	const longdouble null_hypothesis__mean___inferred = sum_over_map_values<uint,longdouble>(it_expected_RD_of_region__NULL->second);
	const longdouble alternative_hypothesis__mean___inferred = sum_over_map_values<uint,longdouble>(it_expected_RD_of_region__ALTERNATIVE->second);
	const longdouble def_alt_del_hypothesis__mean___inferred = sum_over_map_values<uint,longdouble>(it_expected_RD_of_region__def_alt_del->second);
	const longdouble def_alt_dup_hypothesis__mean___inferred = sum_over_map_values<uint,longdouble>(it_expected_RD_of_region__def_alt_dup->second);
	
	const longdouble null_hypothesis__mean___unique
		    = sum_over_map_values<uint, longdouble>(
			erase_set_from_map_keys<uint,longdouble>(it_expected_RD_of_region__NULL->second, nonunique_positions_in_region_of_interest));
	const longdouble alternative_hypothesis__mean___unique
		    = sum_over_map_values<uint, longdouble>(
			erase_set_from_map_keys<uint,longdouble>(it_expected_RD_of_region__ALTERNATIVE->second, nonunique_positions_in_region_of_interest));
	const longdouble def_alt_del_hypothesis__mean___unique
		    = sum_over_map_values<uint, longdouble>(
			erase_set_from_map_keys<uint,longdouble>(it_expected_RD_of_region__def_alt_del->second, nonunique_positions_in_region_of_interest));    
	const longdouble def_alt_dup_hypothesis__mean___unique
		    = sum_over_map_values<uint, longdouble>(
			erase_set_from_map_keys<uint,longdouble>(it_expected_RD_of_region__def_alt_dup->second, nonunique_positions_in_region_of_interest));			    

		    
		    
	const longdouble sigma_inferred = get_standard_deviation_from_mean_for_RD(null_hypothesis__mean___inferred);
	const longdouble sigma_UNIQUE = get_standard_deviation_from_mean_for_RD(null_hypothesis__mean___unique);
	
	const longdouble standardized_Z_score___inferred = (observed_counts__inferred - null_hypothesis__mean___inferred)/sigma_inferred;
	const longdouble standardized_Z_score___UNIQUE = (observed_counts__UNIQUE_ONLY - null_hypothesis__mean___unique)/sigma_UNIQUE;
	
	const real standard_normal_tail___inferred(calculate_standard_normal_tail(standardized_Z_score___inferred));
	const real standard_normal_tail___UNIQUE(calculate_standard_normal_tail(standardized_Z_score___UNIQUE));
	
	
	const real null_likelihood__inferred(calculate_normal_pdf_or_cdf(
							observed_counts__inferred,
							null_hypothesis__mean___inferred,
							get_standard_deviation_from_mean_for_RD(null_hypothesis__mean___inferred),
							false));
	
	const real null_likelihood__unique(calculate_normal_pdf_or_cdf(
							observed_counts__UNIQUE_ONLY,
							null_hypothesis__mean___unique,
							get_standard_deviation_from_mean_for_RD(null_hypothesis__mean___unique),
							false));
	
	const real alternative_likelihood__inferred(calculate_normal_pdf_or_cdf(
							observed_counts__inferred,
							alternative_hypothesis__mean___inferred,
							get_standard_deviation_from_mean_for_RD(alternative_hypothesis__mean___inferred),
							false));		
	
	const real alternative_likelihood__unique(calculate_normal_pdf_or_cdf(
							observed_counts__UNIQUE_ONLY,
							alternative_hypothesis__mean___unique,
							get_standard_deviation_from_mean_for_RD(alternative_hypothesis__mean___unique),
							false));
	
	const real def_alt_del_likelihood__inferred(calculate_normal_pdf_or_cdf(
							observed_counts__inferred,
							def_alt_del_hypothesis__mean___inferred,
							get_standard_deviation_from_mean_for_RD(def_alt_del_hypothesis__mean___inferred),
							false));
	
	const real def_alt_del_likelihood__unique(calculate_normal_pdf_or_cdf(
							observed_counts__UNIQUE_ONLY,
							def_alt_del_hypothesis__mean___unique,
							get_standard_deviation_from_mean_for_RD(def_alt_del_hypothesis__mean___unique),
							false));	
	

	const real def_alt_dup_likelihood__inferred(calculate_normal_pdf_or_cdf(
							observed_counts__inferred,
							def_alt_dup_hypothesis__mean___inferred,
							get_standard_deviation_from_mean_for_RD(def_alt_dup_hypothesis__mean___inferred),
							false));	
	
	const real def_alt_dup_likelihood__unique(calculate_normal_pdf_or_cdf(
							observed_counts__UNIQUE_ONLY,
							def_alt_dup_hypothesis__mean___unique,
							get_standard_deviation_from_mean_for_RD(def_alt_dup_hypothesis__mean___unique),
							false));
	
	
	const real  odds_ratio__alt___inferred(null_likelihood__inferred/alternative_likelihood__inferred);
	const real  odds_ratio__alt___unique(null_likelihood__unique/alternative_likelihood__unique);
	
	const real  odds_ratio__def_del___inferred(null_likelihood__inferred/def_alt_del_likelihood__inferred);
	const real  odds_ratio__def_del___unique(null_likelihood__unique/def_alt_del_likelihood__unique);
	
	const real  odds_ratio__def_dup___inferred(null_likelihood__inferred/def_alt_dup_likelihood__inferred);
	const real  odds_ratio__def_dup___unique(null_likelihood__unique/def_alt_dup_likelihood__unique);
	
	const real odds_ratio__best_default__inferred(mpfr::min(odds_ratio__def_del___inferred, odds_ratio__def_dup___inferred));
	const real odds_ratio__best_default__unique(mpfr::min(odds_ratio__def_del___unique, odds_ratio__def_dup___unique));
	
	
	
	std::cerr << "\n\nwhat will be written to file:\n"
		    << observed_counts__inferred
		    << "\n" << observed_counts__UNIQUE_ONLY
		    
		    << "\n" << null_hypothesis__mean___inferred
		    << "\n" << alternative_hypothesis__mean___inferred
		    << "\n" << def_alt_del_hypothesis__mean___inferred
		    << "\n" << def_alt_dup_hypothesis__mean___inferred			    

		    << "\n" << null_hypothesis__mean___unique
		    << "\n" << alternative_hypothesis__mean___unique
		    << "\n" << def_alt_del_hypothesis__mean___unique
		    << "\n" << def_alt_dup_hypothesis__mean___unique //10
					    
		    << "\n" << sigma_inferred
		    << "\n" << standardized_Z_score___inferred			    
		    << "\n" << sigma_UNIQUE
		    << "\n" << standardized_Z_score___UNIQUE	

		    << "\n" << standard_normal_tail___inferred.toString(10)
		    << "\n" << standard_normal_tail___UNIQUE.toString(10)
		    
		    << "\n" << null_likelihood__inferred.toString(10)
		    << "\n" << null_likelihood__unique.toString(10)
		    << "\n" << alternative_likelihood__inferred.toString(10)
		    << "\n" << alternative_likelihood__unique.toString(10)  //20
		    
		    << "\n" << def_alt_del_likelihood__inferred.toString(10)
		    << "\n" << def_alt_del_likelihood__unique.toString(10)
		    << "\n" << def_alt_dup_likelihood__inferred.toString(10)
		    << "\n" << def_alt_dup_likelihood__unique.toString(10)			    
		    
		    << "\n" << odds_ratio__alt___inferred.toString(10)
		    << "\n" << odds_ratio__alt___unique.toString(10) //26
		    
		    << "\n" << odds_ratio__def_del___inferred.toString(10)
		    << "\n" << odds_ratio__def_del___unique.toString(10)
		    
		    << "\n" << odds_ratio__def_dup___inferred.toString(10)
		    << "\n" << odds_ratio__def_dup___unique.toString(10)
		    
		    << "\n" << odds_ratio__best_default__inferred.toString(10)
		    << "\n" << odds_ratio__best_default__unique.toString(10);		
	
	
	
	{//save to file			    		
	    char out_RD_hypothesis_tests__fname[option_size];	
	    std::sprintf(out_RD_hypothesis_tests__fname, "%s.tests", out_RD_core_fname);
					
	    std::ofstream  outfs(out_RD_hypothesis_tests__fname);
	    if ( check_fstream_for_goodness(outfs, out_RD_hypothesis_tests__fname,
					"create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights", false))
	    {				
		outfs 
		    << observed_counts__inferred
		    << "\n" << observed_counts__UNIQUE_ONLY
		    
		    << "\n" << null_hypothesis__mean___inferred
		    << "\n" << alternative_hypothesis__mean___inferred
		    << "\n" << def_alt_del_hypothesis__mean___inferred
		    << "\n" << def_alt_dup_hypothesis__mean___inferred			    

		    << "\n" << null_hypothesis__mean___unique
		    << "\n" << alternative_hypothesis__mean___unique
		    << "\n" << def_alt_del_hypothesis__mean___unique
		    << "\n" << def_alt_dup_hypothesis__mean___unique //10
					    
		    << "\n" << sigma_inferred
		    << "\n" << standardized_Z_score___inferred			    
		    << "\n" << sigma_UNIQUE
		    << "\n" << standardized_Z_score___UNIQUE	

		    << "\n" << standard_normal_tail___inferred.toString(10)
		    << "\n" << standard_normal_tail___UNIQUE.toString(10)
		    
		    << "\n" << null_likelihood__inferred.toString(10)
		    << "\n" << null_likelihood__unique.toString(10)
		    << "\n" << alternative_likelihood__inferred.toString(10)
		    << "\n" << alternative_likelihood__unique.toString(10)  //20
		    
		    << "\n" << def_alt_del_likelihood__inferred.toString(10)
		    << "\n" << def_alt_del_likelihood__unique.toString(10)
		    << "\n" << def_alt_dup_likelihood__inferred.toString(10)
		    << "\n" << def_alt_dup_likelihood__unique.toString(10)			    
		    
		    << "\n" << odds_ratio__alt___inferred.toString(10)
		    << "\n" << odds_ratio__alt___unique.toString(10) //26
		    
		    << "\n" << odds_ratio__def_del___inferred.toString(10)
		    << "\n" << odds_ratio__def_del___unique.toString(10)
		    
		    << "\n" << odds_ratio__def_dup___inferred.toString(10)
		    << "\n" << odds_ratio__def_dup___unique.toString(10)
		    
		    << "\n" << odds_ratio__best_default__inferred.toString(10)
		    << "\n" << odds_ratio__best_default__unique.toString(10);	
		    
		outfs.close();
	    }//outfs
	}//save to file    
	
	    
	
	

	
	{//save to file
	    char out_RD_info__fname[option_size];	
	    std::sprintf(out_RD_info__fname, "%s.info", out_RD_core_fname);
					
	    std::ofstream  outfs(out_RD_info__fname);
	    if (check_fstream_for_goodness(outfs, out_RD_info__fname,
					"create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights", false))
	    {
		const type_BI__BI approp_brkpts_absolute_first(
				    convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
							    some_result.the_diploid_profile_brkpts.first,
							    result_Ev_it->second.compressed_map_LCR_to_profile__for_each_LCR,
							    some_result.the_diploid_Event_outcome.first));
					    
		const type_BI__BI approp_brkpts_absolute_second(
				    convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
							    some_result.the_diploid_profile_brkpts.second,
							    result_Ev_it->second.compressed_map_LCR_to_profile__for_each_LCR,
							    some_result.the_diploid_Event_outcome.second));

		outfs 
		    << genome_name  << "\n"
		    << population_abbreviation  << "\n"
		    << some_result.cc_of_event  << "\n"
		    << some_result.event_UID  << "\n"
		    << convert_haploid_outcome_to_string(some_result.the_diploid_Event_outcome.first)  << "\n"
		    << convert_haploid_outcome_to_string(some_result.the_diploid_Event_outcome.second)  << "\n"
		    
		    << result_Ev_it->second.chromos[0]  << "\n"
		    
		    << approp_brkpts_absolute_first.first.lower()  << "\n"
		    << approp_brkpts_absolute_first.first.upper()  << "\n"
		    << approp_brkpts_absolute_first.second.lower()  << "\n"
		    << approp_brkpts_absolute_first.second.upper()  << "\n"
		    
		    << approp_brkpts_absolute_second.first.lower()  << "\n"
		    << approp_brkpts_absolute_second.first.upper()  << "\n"
		    << approp_brkpts_absolute_second.second.lower()  << "\n"
		    << approp_brkpts_absolute_second.second.upper()  << "\n"

		    << some_result.the_diploid_profile_brkpts.first.lower()  << "\n"
		    << some_result.the_diploid_profile_brkpts.first.upper()  << "\n"
		    << some_result.the_diploid_profile_brkpts.second.lower()  << "\n"
		    << some_result.the_diploid_profile_brkpts.second.upper();
		    
		outfs.close();
	    }//outfs
	}//save to file  	
						
    
    }//RD hypothesis test
    
    
	
	
	
	
				
	
    if (my_MPI_rank == 0)
    {  print_line_of_markers("==");  }	

}//create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights






































































type_map_uint_to_Sampled_diploid_Event_data  Conn_Comp::get_and_save_MAP_estimate()
{
              
    type_map_uint_to_Sampled_diploid_Event_data map_estimate_for_this_CC;
        
    
    std::stringstream output_MAP_verbose_strm;        
    
    output_MAP_verbose_strm << "\n\n\n\n\n\nMAP posterior:\n\tconnected component    "  <<  Conn_Comp_ID  << "\n";
    output_MAP_verbose_strm.precision(10);
    
    real P_MAP(1.00L);    
    
    for (type_list_Visitation_object::reverse_iterator it_Visit = variable_elimination_schedule.rbegin();
            it_Visit != variable_elimination_schedule.rend();
            ++it_Visit)
    {        
        
        
        output_MAP_verbose_strm << "\n\t\t\tEvent   " << it_Visit->event_for_elimination->UID << "\n";
                               
	it_Visit->load_PSF_from_archive();

        
        //First randomly draw a state for the eliminated Event:
        {//outcome
            const type_multimap_real_to_uint__uint EE_outcome_CDF( it_Visit->my_Forward_partial_sum_function.calculate_and_return_eliminated_Event_diploid_outcome_CDF
                                                                                        (map_estimate_for_this_CC));              
            
            output_MAP_verbose_strm << "\t\t\t\t\t\t\tEE_outcome_CDF:\n";
            for (type_multimap_real_to_uint__uint::const_iterator it_EE_outcome_CDF = EE_outcome_CDF.begin();
                    it_EE_outcome_CDF != EE_outcome_CDF.end();
                    ++it_EE_outcome_CDF)   
	    {
                output_MAP_verbose_strm << "\t\t\t\t\t\t\t\t\t\t\t\t\t"
                            << it_EE_outcome_CDF->second.first 
                            << " &&  " << it_EE_outcome_CDF->second.second 
                            <<"    =   " <<  it_EE_outcome_CDF->first.toString(10,10) << "\n";   
	    }
            
            const type_uint__uint drawn_EE_diploid_states(  (--EE_outcome_CDF.end())->second   );  //since is sorted.    
            const type_haploid_outcome__haploid_outcome drawn_EE_diploid_states__cast(
					    static_cast<type_haploid_outcome>(drawn_EE_diploid_states.first),
					    static_cast<type_haploid_outcome>(drawn_EE_diploid_states.second)  );
            
                    
            //tack on to running sample:  
	    map_estimate_for_this_CC[it_Visit->event_for_elimination->UID]
				= Sampled_diploid_Event_data(Conn_Comp_ID, it_Visit->event_for_elimination->UID,
							     drawn_EE_diploid_states__cast,
							     empty_BI__BI);	    
            
            if ( (--EE_outcome_CDF.end()) != EE_outcome_CDF.begin() )
	    {  P_MAP *=  (   (--EE_outcome_CDF.end())->first  -  (--(--EE_outcome_CDF.end()))->first    );  }
            //otherwise, there is only one possibility, so it has probability 1.00, and there is no need to multiply.
             
        }//outcome
        
        
        
                
        
        //Next, draw breakpoints for these EE states.

        {//breakpoint
            const type_multimap_real_to__BI__BI EE_breakpoint_CDF( it_Visit->my_Forward_partial_sum_function.calculate_and_return_eliminated_Event_diploid_breakpoint_CDF
                                                                                    (map_estimate_for_this_CC)     );
                        
            output_MAP_verbose_strm << "\n\t\t\t\t\t\t\t\t\t\t\tEE_breakpoint_CDF (few with highest probability):\n";
            
            type_multimap_real_to__BI__BI::const_iterator best_few_breakpoints_in_CDF = EE_breakpoint_CDF.end();
            for (uint j=0; j<10; ++j)
	    {
                if (best_few_breakpoints_in_CDF != EE_breakpoint_CDF.begin())
                    --best_few_breakpoints_in_CDF;
                else
                    break;
	    }
            
            
            for (type_multimap_real_to__BI__BI::const_iterator it_EE_brkpt_CDF = best_few_breakpoints_in_CDF;
                    it_EE_brkpt_CDF != EE_breakpoint_CDF.end();
                    ++it_EE_brkpt_CDF) 
	    {
                output_MAP_verbose_strm << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t" <<  it_EE_brkpt_CDF->second.first 
                            << "  &&  " << it_EE_brkpt_CDF->second.second << "    =   "<< it_EE_brkpt_CDF->first <<"\n";
	    }
                            
            output_MAP_verbose_strm << "\n\n";
                                                            

            
            if ( (--EE_breakpoint_CDF.end()) != EE_breakpoint_CDF.begin() )
	    {  P_MAP *=  (   (--EE_breakpoint_CDF.end())->first  -  (--(--EE_breakpoint_CDF.end()))->first    ); }
            //otherwise, there is only one possibility, so it has probability 1.00, and there is no need to multiply.            
                                                                     
            
            const type_BI__BI drawn_EE_diploid_breakpoints( (--EE_breakpoint_CDF.end())->second  );      	    	    
	    
            //tack on to running sample:  
	    map_estimate_for_this_CC[it_Visit->event_for_elimination->UID].the_diploid_profile_brkpts = drawn_EE_diploid_breakpoints;	    	                                                                                                                       
        }//breakpoint
                                        
                                 
        it_Visit->my_Forward_partial_sum_function.clear_data();
    
    } // end of    for-loop  it_Visit
        
        
        
        
        
    if (my_MPI_rank == 0)
    {                
        output_MAP_verbose_strm << "\n\t\t\t\t\t\t\t\t\t\t\t\t\t\tP(MAP)   =   " << P_MAP.toString(10) << "\n\n";
            
        std::fprintf(stderr, "\n\n%s\n\n", output_MAP_verbose_strm.str().c_str() );         
            

        {//posterior
            char posterior_filename[option_size];
            std::sprintf(posterior_filename, "%s/posterior_conditional_CDFs_for_MAP_estimates",
                                            output_dir.c_str());
            std::ofstream posterior_filestrm(  posterior_filename,  std::ios::app );
            posterior_filestrm << output_MAP_verbose_strm.str();
            posterior_filestrm.close();
        }//posterior
        
        
        //write to file
        char save_filename[option_size];
        std::sprintf(save_filename, "%s/MAP_estimates",
                                    output_dir.c_str());
                                    
        FILE *save_file = std::fopen(save_filename, "a");
	if (save_file == NULL)
	{
	    std::stringstream errorstrm;
	    errorstrm << "\n\n\n\n\nERROR: unable to open file [" << save_filename  << "] in \"save MAP estimates\".\n\n";
	    error_message(errorstrm, false);
	}        
	else
	{                
	    std::fprintf(save_file, "\n\nConn_Comp:   %u\n",
					Conn_Comp_ID);
				    
	    save_sampled_Connected_Component_to_file( save_file, map_estimate_for_this_CC );         
	    
	    std::fprintf(save_file, "\t\tP_MAP  =  %s\n\n", 
				    P_MAP.toString().c_str() );
	    
	    
	    std::fclose(save_file);
	}
    } // 0
    
    
    
    
    
    
    return   map_estimate_for_this_CC;
      
}  // end  of   get_MAP_estimate

























type_map_uint_to_Sampled_diploid_Event_data Conn_Comp::get_a_single_sample
							    ()  const
{
    
    static boost::mt19937 mersenne_engine(71);
    static boost::uniform_01<boost::mt19937> my_unif_gen(mersenne_engine);      
    
    type_map_uint_to_Sampled_diploid_Event_data  a_sample_from_the_posterior;
    
    for (type_list_Visitation_object::const_reverse_iterator it_Visit = variable_elimination_schedule.rbegin();
            it_Visit != variable_elimination_schedule.rend();
            ++it_Visit)
    {        
        type_uint__uint drawn_EE_diploid_states;
	
        //First randomly draw a state for the eliminated Event:
        {//outcome
            const type_multimap_real_to_uint__uint EE_outcome_CDF(
		    it_Visit->my_Forward_partial_sum_function.calculate_and_return_eliminated_Event_diploid_outcome_CDF
                                                                                        (a_sample_from_the_posterior)      );
                                    
            //draw a uniform random number and see where it falls:
            const longdouble random_draw = my_unif_gen();
            
                
            
            for (type_multimap_real_to_uint__uint::const_iterator it_CDF = EE_outcome_CDF.begin();
                    it_CDF != EE_outcome_CDF.end();
                    ++it_CDF)   
	    {
                if ( it_CDF->first  >=  random_draw   )
                {
                    drawn_EE_diploid_states = it_CDF->second;                
                    break;
                }  
	    }
                                                                                                                                                                                                                                        
        }//outcome
        
        
        a_sample_from_the_posterior[it_Visit->event_for_elimination->UID]
		= Sampled_diploid_Event_data( it_Visit->event_for_elimination->my_Conn_Comp->Conn_Comp_ID,
						it_Visit->event_for_elimination->UID,
					      type_haploid_outcome__haploid_outcome(static_cast<type_haploid_outcome>(drawn_EE_diploid_states.first),
										    static_cast<type_haploid_outcome>(drawn_EE_diploid_states.second)),
					       empty_BI__BI);
						      						     						      
        
	type_BI__BI drawn_EE_diploid_breakpoints; 
        
        //Next, draw breakpoints for these EE states.
        {//breakpoint
            const type_multimap_real_to__BI__BI EE_breakpoint_CDF( 
			it_Visit->my_Forward_partial_sum_function.calculate_and_return_eliminated_Event_diploid_breakpoint_CDF
                                                                                    (a_sample_from_the_posterior)     );
                                                                                                                    
            //draw a uniform random number and see where it falls:
            const longdouble random_draw = my_unif_gen();        
        
                                
        
            for (type_multimap_real_to__BI__BI::const_iterator it_CDF = EE_breakpoint_CDF.begin();
                    it_CDF != EE_breakpoint_CDF.end();
                    ++it_CDF)   
	    {
                if (it_CDF->first >= random_draw )
                {
                    drawn_EE_diploid_breakpoints = it_CDF->second;
                    break;
                }                            
	    }
                                                                                                                     
        }//breakpoint  
        
		
		
	a_sample_from_the_posterior[it_Visit->event_for_elimination->UID].the_diploid_profile_brkpts = drawn_EE_diploid_breakpoints;
    
    } // end of    for-loop  it_Visit
    
    
    
    return a_sample_from_the_posterior;
    
}  // end  of   get_a_single_sample







































// typedef  std::pair<type_map_uint_to_uint, type_map_uint_to_uint>  type_haploid_sample; // event outcomes , breakpoints.
// typedef  std::pair<type_haploid_sample, type_haploid_sample>  type_diploid_sample;
// typedef  std::vector< type_diploid_sample >  type_vector_diploid_sample;




// typedef  std::vector< type_map_uint_to__uint__uint__2 >  type_vector__map_uint_to__uint__uint__2;
// // first  = Event outcomes.  second = breakpoints for those outcomes.
// 
// 
type_map_uint_to_Marginal_Event_posterior_distribution  Conn_Comp::compute_marginal_probabilities_from_posterior_sample
                                                (const type_vector_map_uint_to_Sampled_diploid_Event_data &posterior_samples)  const
{
        
    type_map_uint_to_Marginal_Event_posterior_distribution   marginal_probability_distribution_per_Event;              
                for (type_map_uint_to_Event::const_iterator it_ev = events.begin();
                        it_ev != events.end();
                        ++it_ev)
		{  marginal_probability_distribution_per_Event[it_ev->first] = Marginal_Event_posterior_distribution(Conn_Comp_ID, it_ev->first);  }		
         

    
    const uint number_of_samples = posterior_samples.size();
    
    for (uint j=0; j<number_of_samples; ++j)    
    {
	for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_sampled_ev = posterior_samples.at(j).begin();
		it_sampled_ev != posterior_samples.at(j).end();
		++it_sampled_ev)
	{	
	    type_haploid_outcome__haploid_outcome haploid_neutral__outcomes(it_sampled_ev->second.the_diploid_Event_outcome);	    
	    type_BI__BI  haploid_neutral__brkpts(it_sampled_ev->second.the_diploid_profile_brkpts);
						   
	    if ((int)haploid_neutral__outcomes.first > (int)haploid_neutral__outcomes.second)
	    {
		std::swap<type_haploid_outcome>(haploid_neutral__outcomes.first, haploid_neutral__outcomes.second);
		std::swap<BOOST_Interval>(haploid_neutral__brkpts.first, haploid_neutral__brkpts.second);
	    }
						   
	    
	    ++marginal_probability_distribution_per_Event.at(it_sampled_ev->second.event_UID).map_diploid_outcomes_to_posterior_probability[haploid_neutral__outcomes];
	    ++marginal_probability_distribution_per_Event.at(it_sampled_ev->second.event_UID).map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability
			[haploid_neutral__outcomes][haploid_neutral__brkpts];
	}//sampled_ev	    	    
    }//sample_ctr
    
    
        
        
    std::cerr << "normalize by number of samples:\n";
    //normalize by number of Events:
    for (type_map_uint_to_Marginal_Event_posterior_distribution::iterator it_ev = marginal_probability_distribution_per_Event.begin();
            it_ev != marginal_probability_distribution_per_Event.end();
            ++it_ev)  
    {
        for (type_map__hap_outcome__hap_outcome__to__longdouble::iterator it_outcome = it_ev->second.map_diploid_outcomes_to_posterior_probability.begin();
                it_outcome != it_ev->second.map_diploid_outcomes_to_posterior_probability.end();
                ++it_outcome)
        {
            it_outcome->second /= number_of_samples;
            
            for (type_map__BI__BI__to__longdouble::iterator it_brkpts
				= it_ev->second.map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability.at(it_outcome->first).begin();
                    it_brkpts != it_ev->second.map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability.at(it_outcome->first).end();
                    ++it_brkpts)
	    {  it_brkpts->second /= number_of_samples;  }
        }//outcome
    }//ev
       
             
    return marginal_probability_distribution_per_Event;


}  //  compute_marginal_probabilities_from_posterior_sample





























void Conn_Comp::sample_from_posterior__and__save____marginal_distributions__and__outcome_Centroid__and_Shannon__Entropies
                                                (const uint &number_of_samples)
{
    
    
    std::fprintf(stderr, "\n\n\nsampling from posterior...\n");
    
    type_vector_map_uint_to_Sampled_diploid_Event_data samples_from_posterior(number_of_samples);    
    {
	//load all PSFs	
	for (type_list_Visitation_object::iterator it_Visit = variable_elimination_schedule.begin();
		it_Visit != variable_elimination_schedule.end();
		++it_Visit)                                  
	{  it_Visit->load_PSF_from_archive();  }
    
	#pragma omp parallel
	{//parallel          
	    
	    #pragma omp for  schedule(dynamic,20)
	    for (uint s=0; s < number_of_samples; ++s)
	    {  samples_from_posterior.at(s) = get_a_single_sample();  }
	    
	    mpfr_free_cache();
	}//parallel
	
	
	//clear all PSFs	
	for (type_list_Visitation_object::iterator it_Visit = variable_elimination_schedule.begin();
		it_Visit != variable_elimination_schedule.end();
		++it_Visit)                                  
	{  it_Visit->my_Forward_partial_sum_function.clear_data();  }	    
    }

    save_posterior_samples_to_file(samples_from_posterior, Conn_Comp_ID, number_of_samples, output_dir);
    

    
    

   
    double time_marginals_begin = omp_get_wtime();   
    std::fprintf(stderr, "\n\n\nmarginal_distributions...\n"); 
    const type_map_uint_to_Marginal_Event_posterior_distribution marginal_distributions
                                                (compute_marginal_probabilities_from_posterior_sample(samples_from_posterior));       
    std::cerr<< "\n\nrank " << my_MPI_rank <<  "   time    =     "  << omp_get_wtime() - time_marginals_begin << "\n\n";       
    
    
    save_marginal_posterior_distributions_to_file(marginal_distributions, output_dir);
     
    
    


    std::cerr << "\nsave table\n";
        
        
    std::stringstream output_marginal_Positive_calls_strm;
    std::stringstream output_marginal_Negative_calls_strm;
    std::stringstream output_marginal_all_calls_strm;
    
    output_marginal_Positive_calls_strm.precision(10);
    output_marginal_Negative_calls_strm.precision(10);
    output_marginal_all_calls_strm.precision(10);
    
    
    
    make_MAP_calls_from_marginal_posterior_distributions_per_Event
                                (marginal_distributions,
                                 *this,
                                 output_marginal_Positive_calls_strm,
                                output_marginal_Negative_calls_strm,
                                 output_marginal_all_calls_strm);    
        
        
        
                                 
//     save_calls_tables_to_file
//         (output_marginal_Positive_calls_strm,
//             output_marginal_Negative_calls_strm,
//             output_marginal_all_calls_strm,
//          output_dir);  
                                    

    std::cerr << "\n\n\ndone!\n\n\n\n\n\n";    

}  // sample_from_posterior__and__save____marginal_distributions__and__outcome_Centroid__and_Shannon__Entropies















































































void Conn_Comp::simulate_elimination_schedule_and_save_natural_poisson_intervals()
{
    
    
    std::cerr << "\n\ncreating natural poisson subdirectory of Visitation Schedules subdir [" << visitation_schedule_subdirname << "]...\n\n";
    {    
        char sub_out_dirs[option_size];         
        std::sprintf(sub_out_dirs, "%s/%s/natural_poisson_intervals",
                                data_directory.c_str(), visitation_schedule_subdirname.c_str() );  
                                
        boost::filesystem::create_directory(sub_out_dirs);                              
    }

    
    
    
    
    
    for (type_list_Visitation_object::iterator it_Visit = variable_elimination_schedule.begin();
         it_Visit != variable_elimination_schedule.end();
         ++it_Visit)
    {      
                        
        
//         Event *const event_being_summed_out = it_Visit->event_for_elimination;   // readability.
        const type_map_uint_to_Event::iterator  event_being_summed_out = events.find(  it_Visit->event_for_elimination->UID  );   // readability.
                   
        std::fprintf(stderr, "\n\n\n\nnatural poisson intervals for:     event_being_summed_out UID = %u\n\tprofile_length = %u\n\tremaining_coordinates_for_consideration.size() = %u\n\n",
                            event_being_summed_out->second.UID,
                            event_being_summed_out->second.profile_length,
                            event_being_summed_out->second.remaining_coordinates_for_consideration.size()  );
        
                    
        event_being_summed_out->second.print_this_entry();
                                
        
        // just in case we need them, we define these here to have them in scope...
        type_map_uint_to_set_uint variational_positions_on_neighbors_and_self;  
        
        
        if ( !event_being_summed_out->second.remaining_coordinates_for_consideration.empty() )
        {
            
                      

            //Read-depth preparation:
            
            //set base multiplicity:   

       
            double time_begin = omp_get_wtime();
            if (my_MPI_rank == 0)
                std::fprintf(stderr, "rank %d:  partition_DAR_by_homology    and    create_natural_poisson_intervals.\n", my_MPI_rank);
            
            event_being_summed_out->second.partition_DAR_by_homology___and___create_natural_poisson_intervals___NOT__gender_sensitive();
            
            std::cerr<< "\n\nrank " << my_MPI_rank <<  "   time    =     "  << omp_get_wtime() - time_begin << "\n\n";               
            
            
            
            
            {//save natural poisson intervals
                   
                char nat_poisson_filename[option_size];
                
                std::sprintf(nat_poisson_filename, "%s/%s/natural_poisson_intervals/natural_poisson_intervals_for_Event_%u",
                                                    data_directory.c_str(),
                                                    visitation_schedule_subdirname.c_str(),
                                                    event_being_summed_out->second.UID );
                                                                                                        
                std::ofstream out_fs( nat_poisson_filename,   std::ios::out  );  //  std::ios_base::binary     
                                        if ( !out_fs.good()  )
                                        {
                                            std::stringstream error_strm;
                                            print_line_of_markers("ERROR! ", &error_strm);
                                            print_line_of_markers("(", &error_strm);
                                            error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                                            
                                            error_strm << "\n\nERROR!   not good ofstream   in  \"save natural poisson interval\",\n\t\t nat_poisson_filename =    "  
                                            <<  nat_poisson_filename  << "\n\n";
                                            
                                            print_line_of_markers(")", &error_strm);
                                            std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
                                            return;
                                        }    

        
//                 boost::archive::text_oarchive out_archive( out_fs );
//     
//                 out_archive << event_being_summed_out->second.natural_poisson_intervals;

                out_fs << "number_of_colors = " << event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.size();

                for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_nat = event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.begin();
                        it_nat != event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.end();
                        ++it_nat)
                {
                    out_fs << "\ncolor:  " << it_nat->first << "      number_chromos =  " << it_nat->second.size();
                    for (type_map_uint_to_uint_to_longdouble::const_iterator it_chr = it_nat->second.begin();
                            it_chr != it_nat->second.end();
                            ++it_chr)
                    {                                                
                        const type_set_uint individ_coords__set(  extract_keys_of_map_and_return_as_set<uint,longdouble>(it_chr->second)   );
                        
                        //for ease of conversion
                        type_map_uint_to_uint storage_coords;
                                                {//fill
                                                    type_map_uint_to_uint::iterator insert_it = storage_coords.begin();
                                                    for (type_set_uint::const_iterator it_coord = individ_coords__set.begin();
                                                            it_coord != individ_coords__set.end();
                                                            ++it_coord)
                                                        insert_it = storage_coords.insert(   insert_it,   type_uint__uint(*it_coord, *it_coord)  );
                                                }//fill
                                                
                        const type_vector_vector_int compressed_coords(    convert_full_map_to_compressed_map(  storage_coords  )      );    
                        
                        out_fs << "\n\t\tchromo  " << it_chr->first << "       compressed_size = " << compressed_coords.at(0).size();
                        
                        for (int comp = 0; comp < compressed_coords.at(0).size(); ++comp)
                            out_fs << "\n\t\t\t" << compressed_coords.at(0).at(comp)  << "      " << compressed_coords.at(1).at(comp);                   
                    } // chr                                
                }//nat
                                
                
                out_fs.close();                  
                         
            }//save natural poisson intervals
            
            
            
            
            
            variational_positions_on_neighbors_and_self = event_being_summed_out->second.determine_homologous_variational_positions_on_immediate_neighbors_and_self();   
            
            if (my_MPI_rank == 0)
                print_map_to_set<uint, uint>(variational_positions_on_neighbors_and_self, "variational_positions_on_neighbors_and_self");
            

            
        } //remaining coordinates not empty.
        
            
            
            
            
        
        

        
        event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.clear();
        event_being_summed_out->second.base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive.clear();
	event_being_summed_out->second.base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.clear();
        event_being_summed_out->second.observed_read_depth_for_natural_poisson_intervals.clear();         

        event_being_summed_out->second.PERs_on_this_profile.clear();
        event_being_summed_out->second.filter_for_NON_large_gap_of_neighbors_and_self__that_affect_Read_depth.clear();
        event_being_summed_out->second.RD_affecting_neighbors_self.clear();
        
        event_being_summed_out->second.large_gaps_of_profile.clear();
        
                        

        

 
        
        
        
        //housekeeping:
        if (my_MPI_rank == 0)
            print_line_of_markers("housekeeping ");                        
                    

        
        
        
     
        
     
        
        
        

      
        if (my_MPI_rank == 0)        
            std::fprintf(stderr, "debug:  erasing positions from all related events...\n");                           
        //FINALLY,   ERASE GENOME POSITIONS FROM ALL RELATED EVENTS (SO THE SAME GENOME POSITIONS ARENT CONSIDERED MORE THAN ONCE)
        //What about genome positions of the current event that weren't on anybody elses homology lists?  (i.e. unique to this LCR-pair)???????                        
                        



        type_map_uint_to_Interlock* local_and_global_list[2];
        local_and_global_list[0] = &(event_being_summed_out->second.local_interlocking_Events);
        local_and_global_list[1] = &(event_being_summed_out->second.global_interlocking_Events);

        for (uint lg = 0; lg<2; ++lg)                           
            for (type_map_uint_to_Interlock::iterator ev_it = local_and_global_list[lg]->begin();
                    ev_it != local_and_global_list[lg]->end();
                    ++ev_it)   // for every interlocking event...
            {
                		
                const  type_map_uint_to_Event::iterator  it_interlocking_Event =   events.find(ev_it->first);
                                                                
                //erase yourself from this interlocking event's lists:
                it_interlocking_Event->second.local_interlocking_Events.erase(event_being_summed_out->second.UID);
                it_interlocking_Event->second.global_interlocking_Events.erase(event_being_summed_out->second.UID);                
                
            
                    
                //and erase your very coordinates.  (necessary, for say, Inbetween region).
                if (event_being_summed_out->second.chromos[0] ==  it_interlocking_Event->second.chromos[0])  
                {
                    for (uint pos =  event_being_summed_out->second.region_between_and_including_the_LCRs_themselves.lower();
                              pos <= event_being_summed_out->second.region_between_and_including_the_LCRs_themselves.upper();
                              ++pos)                   
                        it_interlocking_Event->second.remaining_coordinates_for_consideration.erase(pos);               
                }
                
                
                    
                    
//                 for (type_set_uint::const_iterator it_remaining = event_being_summed_out->second.remaining_coordinates_for_consideration.begin();
//                         it_remaining != event_being_summed_out->second.remaining_coordinates_for_consideration.end();
//                         ++it_remaining)
//                     it_interlocking_Event->second.remaining_coordinates_for_consideration.erase(*it_remaining);
                
                
                
                
                for (type_list_Coords_and_homologous_region::const_iterator it_homol
                                = event_being_summed_out->second.regions_homologous_to_directly_affected_region.begin();
                        it_homol != event_being_summed_out->second.regions_homologous_to_directly_affected_region.end();
                        ++it_homol)              
                    if (it_homol->chromosome_of_homologous_region == it_interlocking_Event->second.chromos[0]
                        and  BOOST_overlap( get_endpoints_of_keys_of_compressed_map(it_homol->compressed_map_homologous_region_to_profile),
                                            it_interlocking_Event->second.region_between_and_including_the_LCRs_themselves )      )
                    {
                        const type_map_uint_to_uint full_map( convert_compressed_map_to_full_map(it_homol->compressed_map_homologous_region_to_profile) );                        
                        for (type_map_uint_to_uint::const_iterator map_it = full_map.begin();
                                map_it != full_map.end();
                                ++map_it)                                                
                            it_interlocking_Event->second.remaining_coordinates_for_consideration.erase( map_it->first );                                                                      
                    }  
            
            
            
                           
                           
                           
                {//homol reg restriction
                    type_list_Coords_and_homologous_region::iterator it_homol
                                = it_interlocking_Event->second.regions_homologous_to_directly_affected_region.begin();                                
                    while(  it_homol != it_interlocking_Event->second.regions_homologous_to_directly_affected_region.end()  )
                    {
                        it_homol->compressed_map_homologous_region_to_profile
                                =  get_compressed_map_inverse_of_compressed_map(
                                        convert_full_map_to_compressed_map(
                                                get_map_intersection_of_compressed_map_keys_and_set(
                                                        get_compressed_map_inverse_of_compressed_map(   it_homol->compressed_map_homologous_region_to_profile    ),
                                                        it_interlocking_Event->second.remaining_coordinates_for_consideration )  )  ); 
                                            
                        if ( it_homol->compressed_map_homologous_region_to_profile.at(0).empty() )
                            it_homol = it_interlocking_Event->second.regions_homologous_to_directly_affected_region.erase( it_homol );
                        else
                            ++it_homol;                                            
                    }                                                        
                }//homol reg restriction                           
             
            
            
            
            
            
                //and erase variational positions from its profile.
                if ( !event_being_summed_out->second.remaining_coordinates_for_consideration.empty() )                              
                    for (type_set_uint::const_iterator it_rel_varpos = variational_positions_on_neighbors_and_self.at(it_interlocking_Event->second.UID).begin();
                            it_rel_varpos != variational_positions_on_neighbors_and_self.at(it_interlocking_Event->second.UID).end();
                            ++it_rel_varpos)
                        it_interlocking_Event->second.variational_positions_of_profile.erase( *it_rel_varpos );
                
                    
                for (uint qs=0; qs<2; ++qs)                
                    it_interlocking_Event->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)                        
                        =  get_compressed_map_inverse_of_compressed_map(
                                convert_full_map_to_compressed_map(
                                    get_map_intersection_of_compressed_map_keys_and_set(
                                            get_compressed_map_inverse_of_compressed_map(
                                                    it_interlocking_Event->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)  ),
                                            it_interlocking_Event->second.variational_positions_of_profile)  ));   
                                                        
                            
                                        
                if (my_MPI_rank == 0)  
                {                
                    std::fprintf(stderr, "UID %u", it_interlocking_Event->second.UID);                        
                    print_set<uint>(it_interlocking_Event->second.variational_positions_of_profile, "it_interlocking_Event->second.variational_positions_of_profile");
                }
                
//                 //and erase your homologous coordinates:
//                 for (type_list_vector_vector_int::const_iterator compress_it = ev_it->second.compressed_maps_homol_regions_to_profile.begin();
//                         compress_it != ev_it->second.compressed_maps_homol_regions_to_profile.end();  //maps from interlocking event ev_it to event_being_summed_out
//                         ++compress_it)
//                 {                    
//                     type_map_uint_to_uint full_map_homol_to_profile(  convert_compressed_map_to_full_map(*compress_it)  ); 
//                      
//                     for (type_map_uint_to_uint::const_iterator map_it = full_map_homol_to_profile.begin();
//                             map_it != full_map_homol_to_profile.end();
//                             ++map_it)                                                
//                         it_interlocking_Event->second.remaining_coordinates_for_consideration.erase( map_it->first  );                    
//                 } // end for-loop   compress_it                                                                
                        
            }  // end for-loop   ev_it
                      
                      
                      
        
        // clean up:
        if (my_MPI_rank == 0)
            std::fprintf(stderr, "debug:  deleting  rest of  vectors...\n");                           
        //This variable has been summed out and no longer contains any information.  delete all of its vectors.
        event_being_summed_out->second.global_interlocking_Events.clear();
        event_being_summed_out->second.local_interlocking_Events.clear();
        event_being_summed_out->second.remaining_coordinates_for_consideration.clear();
        event_being_summed_out->second.my_original_Breakpoint_complex.clear_all_data();
        
        event_being_summed_out->second.regions_homologous_to_LCR_profile.clear();
        event_being_summed_out->second.regions_homologous_to_directly_affected_region.clear();
        event_being_summed_out->second.variational_positions_of_profile.clear();  
        event_being_summed_out->second.compressed_map_LCR_to_profile__for_each_LCR.clear();
        event_being_summed_out->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.clear();        
        
        
    }      // end    for-loop    variable_elimination_schedule






    return;  
    
    
                           

}  //  simulate_elimination_schedule_and_save_natural_poisson_intervals















































int Conn_Comp::use_saved_likelihoods_to_perform_complete_variable_elimination_schedule_with_current_prior
				(const std::string &old_run_data_dir,
				 const type_set_uint &excluded_Events_from_old_run)
{		
    return 0;
    
//     std::fprintf(stderr, "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nrank %d:  use_saved_likelihoods_to_perform_complete_variable_elimination_schedule_with_current_prior  for CC_ID = %u !!!\n\t\t\t\t size = %u\n\n\n",
//                  my_MPI_rank, Conn_Comp_ID, events.size());  
// 		 
// 
// 
//     load_original_breakpoints_after_heuristic_for_all_Events_in_Connected_Component_from_file( old_run_data_dir, events );
// 	
//     
// //     if (my_MPI_rank == 0)
// //     {
// // 	print_set<uint>(excluded_Events_from_old_run, "excluded_Events_from_old_run", NULL, true);
// // 	
// // 	for (type_map_uint_to_Event::const_iterator it_ev = events.begin();
// // 		it_ev != events.end();
// // 		++it_ev)
// // 	{
// // 	    std:: cerr << "\n\tEvent " << it_ev->second.UID << ":\n";
// // 	    print_set<uint>(it_ev->second.my_original_Breakpoint_complex.variational_positions__ALL, "uploaded original_variational_positions");
// // 	}
// //     }
//                  
//          
//          
//     
//     //Perform schedule...
//     
//     
//     for (type_list_Visitation_object::iterator it_Visit = variable_elimination_schedule.begin();
//          it_Visit != variable_elimination_schedule.end();
//          ++it_Visit)
//     {                                                      
// 	
//         const type_map_uint_to_Event::iterator   event_being_summed_out = events.find( it_Visit->event_for_elimination->UID );  // readability. 
//                         
//         
//         if (my_MPI_rank == 0)                
//             std::fprintf(stderr, "\n\n\n\nevent_being_summed_out UID = %u\n\tprofile_length = %u\n\tremaining_coordinates_for_consideration.size() = %u\n\n",
//                                 event_being_summed_out->second.UID,
//                                 event_being_summed_out->second.profile_length,
//                                 event_being_summed_out->second.remaining_coordinates_for_consideration.size()  );
//            
//         if (my_MPI_rank == 0)                         
//             event_being_summed_out->second.print_this_entry();
//                                 
//         
// 	
//         // just in case we need them, we define these here to have them in scope...
//         type_set_uint haploid_UID_arguments_for_new_PSF;                                                        
// 		{//scope	    
// 		    haploid_UID_arguments_for_new_PSF.insert( event_being_summed_out->second.UID );
// 		    
// 		    for (type_map_uint_to_Interlock::const_iterator it_neighb = event_being_summed_out->second.local_interlocking_Events.begin();
// 			    it_neighb != event_being_summed_out->second.local_interlocking_Events.end();
// 			    ++it_neighb)
// 			haploid_UID_arguments_for_new_PSF.insert(   it_neighb->first   );  
// 		    
// 		    for (type_map_uint_to_Interlock::const_iterator it_neighb = event_being_summed_out->second.global_interlocking_Events.begin();
// 			    it_neighb != event_being_summed_out->second.global_interlocking_Events.end();
// 			    ++it_neighb)
// 			haploid_UID_arguments_for_new_PSF.insert(   it_neighb->first   );  	    
// 		}//scope	
//         
// 
// 
//             
//         
//         
//         
//   
//         
//         
//         
//         
//         //Partial Sum preparation:
//         
//         //Identify relevant partial sum functions
//         type_list_Partial_sum_fcn_iter relevant_partial_sum_functions;
//         for (type_list_Partial_sum_fcn::iterator it_all_psf = partial_sum_functions.begin();
//                 it_all_psf != partial_sum_functions.end();
//                 ++it_all_psf)        
//             if (it_all_psf->check_if_event_is_an_argument_to_this_function(event_being_summed_out->second.UID))
//                 relevant_partial_sum_functions.push_back( it_all_psf );
//             
//         if (my_MPI_rank == 0)
//             fprintf(stderr, "debug:  relevant_partial_sum_functions identified.  size = %u\n", relevant_partial_sum_functions.size() );                  
//         
//                                 
//         
//         it_Visit->my_Forward_partial_sum_function.set_eliminating_Event__and__haploid_argument_UIDs_according_to_relevant_PSFs
//                                                         (event_being_summed_out->second.UID,
//                                                          relevant_partial_sum_functions,
//                                                          haploid_UID_arguments_for_new_PSF);
//                                                          
// 
// 
//                         
//         
//         //load permissibile state vectors:
//         if (my_MPI_rank == 0)
//             std::fprintf(stderr, "\n\ndebug:  read_permissible_state_vectors_from_file_for_summing_out_a_given_UID = %u\n", event_being_summed_out->second.UID);
//         
//         {
//             const bool success_read_state_vectors = read_permissibile_state_vectors_from_file_for_summing_out_a_given_UID_and_tack_on_empty_state_vector( (*it_Visit) );                
//             if (!success_read_state_vectors)
//             {
//                 if (my_MPI_rank == 0)
//                 {
//                     std::stringstream  reason_strm;
//                     reason_strm << "failed to read permissible state vectors from file for Event " << event_being_summed_out->second.UID;
//                     append_to_skipped_connected_component_file(  *this,  reason_strm.str() );
//                 }
//                 return -3;            
//             }
//         }
//         
//         
//         
//         
//         
//         {//sanity check     
// //             type_set_uint PSF_neighbors( it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs );            
// //             for (type_set_uint::const_iterator it_saved = it_Visit->neighbors.begin();
// //                     it_saved != it_Visit->neighbors.end();
// //                     ++it_saved)
// //                 PSF_neighbors.erase( *it_saved );
// //             
// //             if  ( ! PSF_neighbors.empty() ) //( ! (PSF_neighbors.size() == 1  and  *(PSF_neighbors.begin()) == event_being_summed_out->second.UID)  )
//             if (  it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs  !=  it_Visit->neighbors   )
//             {
//                 std::stringstream error_strm;
//                 print_line_of_markers("ERROR! ", &error_strm);
//                 print_line_of_markers("(", &error_strm);
//                 error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
//                 
//                 error_strm << "\n\nERROR!   recorded neighbors from Visitation schedule do not match with neighbors during \"perform complete variable elimination schedule\" "  << "\n\n";
//                 
//                 print_set<uint>(it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs, "it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs", &error_strm);
//                 print_set<uint>(it_Visit->neighbors, "it_Visit->neighbors", &error_strm);
//                                 
//                 
//                 print_line_of_markers(")", &error_strm);
//                 std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                        
//             }
//         
//         }//sanity check
//                 
//                 
//                 
//         
//         
//         
//         
//         it_Visit->refine_state_vectors_according_to_excluded_Events( excluded_Events_from_old_run );
// 
//             
//                 
//         std::vector< type_vector_Sparse_map::const_iterator >   iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors;
//         std::vector< type_vector_Sparse_map::const_iterator >   iterators_to_haploid_1_X_or_Y_sensitive_state_vectors; 
// 	
// 	
// 	make_gender_sensitive_iterators_for_state_vectors(
// 					    iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors,
// 					    iterators_to_haploid_1_X_or_Y_sensitive_state_vectors,
// 					    events,
// 					    it_Visit );
//                          
//         
//         const int size_state_vec_for_haploid_0 = iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.size();
//         const int size_state_vec_for_haploid_1 = iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.size();
//         
//         
//         
//       
//         
//         
//         
//        
//         const int haploid_state_compute_size = it_Visit->permissible_haploid_state_vectors_sparse.size();
//         
//         if (my_MPI_rank == 0)
//             std::cerr << "\nthaploid_state_compute_size = " << haploid_state_compute_size
// 			<< "\n\tOLD  it_Visit->diploid_compute_size  = " << it_Visit->diploid_compute_size
// 			<< "\n\tNEW gender sensitive compute size =  " << size_state_vec_for_haploid_0 
// 			<< "  *  " << size_state_vec_for_haploid_1
// 			<< "  =  " << size_state_vec_for_haploid_0  *  size_state_vec_for_haploid_1  << "\n\n\n";
//         
//                         
//         const bool  compute_size_small_enough_to_display_diagnostic_information  =    (size_state_vec_for_haploid_0  *  size_state_vec_for_haploid_1   <  150);                                
//         
//         if (my_MPI_rank == 0)
// 	{
//             if (!compute_size_small_enough_to_display_diagnostic_information)
//                 std::cerr << "\n\tdiploid compute size for Event " << event_being_summed_out->second.UID
//                         << "  is " << size_state_vec_for_haploid_0  *  size_state_vec_for_haploid_1
//                         << ",  which is too big to display diagnostics (will slow down the compute dramatically)\n\n";                        
// 	}
//                         
//                         
//         std::cerr.flush();std::fflush(stderr);         
// 	
// 	
// 	
// 	
// 	
// 	
// 	
// 	
// 	
//         Partial_sum_function rank_specific___piece_of_Visitation_PSF(it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs,
//                                                                      it_Visit->event_for_elimination->UID ); 
// 								     
//         #pragma omp parallel
//         { // parallel
//         
//             const int my_thread_ID = omp_get_thread_num();
//         
//             Partial_sum_function thread_specific___piece_of_Visitation_PSF(it_Visit->my_Forward_partial_sum_function.haploid_argument_UIDs,
//                                                                            it_Visit->event_for_elimination->UID ); 
// 									   
// 						    
//             
// 									   
// 	    Partial_sum_function thread_specific___piece_of_new_PSF(  haploid_UID_arguments_for_new_PSF  );
// 	    
// 	    thread_specific___piece_of_new_PSF.upload_all_likelihood_factors_from_file_for_a_given_Event
// 				    (old_run_data_dir,
// 				    event_being_summed_out->second.UID);                                   
// 	    
// 	    
//         
//             #pragma omp for  collapse(2)  schedule(dynamic,1) nowait
//             for (int state_vec_hap_0_ctr = my_MPI_rank;
//                     state_vec_hap_0_ctr < size_state_vec_for_haploid_0;
//                     state_vec_hap_0_ctr += size_MPI_WORLD )
//             for (int state_vec_hap_1_ctr = 0;
//                     state_vec_hap_1_ctr < size_state_vec_for_haploid_1;
//                     ++state_vec_hap_1_ctr)                    
//             {
//                 
//                 const type_vector_Sparse_map::const_iterator it_states_0 = iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.at(state_vec_hap_0_ctr);
//                 const type_vector_Sparse_map::const_iterator it_states_1 = iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.at(state_vec_hap_1_ctr);
//                 
//                 const int total_diploid_states_ctr =   (  ((state_vec_hap_0_ctr+1) - 1)  *  size_state_vec_for_haploid_1  )  +  (state_vec_hap_1_ctr+1);
//                                                             // +1 because we adjust from base 0 to base 1 for ease of computation.        
// 							    
// 		   
// 		thread_specific___piece_of_Visitation_PSF.
// 		    save_Cartesian_product_of_related_PSFs_wrt_eliminating_Event_and_include_probabilities_of_Breakpoint_and_Event
// 									(event_being_summed_out->second,
// 									 *it_states_0,   *it_states_1,                                           
// 									relevant_partial_sum_functions,
// 									thread_specific___piece_of_new_PSF);  
//                                                                         
//             }     // states
//                              
//             
// 
//             #pragma omp critical (absorbing_into_rank_visitation_PSF)
//             rank_specific___piece_of_Visitation_PSF.absorb_piece_of_Partial_Sum_Function_into_one_Partial_Sum_Function
//                                                                         (thread_specific___piece_of_Visitation_PSF);
// 
// 
//             mpfr_free_cache();                
//         } // omp parallel
// 
// 
// 
// 
//     
// 
// 
// 
// 
// 
// 
//         
//         iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.clear();
//         iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.clear(); 
// 
// 
//         haploid_UID_arguments_for_new_PSF.clear();        
//         
//         event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.clear();
//         event_being_summed_out->second.base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive.clear();
//         event_being_summed_out->second.observed_read_depth_for_natural_poisson_intervals.clear();         
// 
//         event_being_summed_out->second.PERs_on_this_profile.clear();
//         event_being_summed_out->second.filter_for_NON_large_gap_of_neighbors_and_self__that_affect_Read_depth.clear();
//         event_being_summed_out->second.RD_affecting_neighbors_self.clear();

//         
//         event_being_summed_out->second.large_gaps_of_profile.clear();
//         
//                         
//         
//         it_Visit->neighbors.clear();
//         
//         it_Visit->permissible_haploid_state_vectors_sparse.clear();
// 
//             
// 
//  
//         
//         
//         world_MPI_boost_ptr->barrier();  //because  before we erase an old PSF, all ranks must be done using it (i.e. done with above loop).
//         
//         std::fflush(stderr);std::cerr.flush();
//         std::fprintf(stderr, "\nrank  %d   erasing old PSFs...\n\n", my_MPI_rank);
//         std::fflush(stderr);std::cerr.flush();
//         
//         
//         
//         
//         //Partial sum function:                        
//         for (type_list_Partial_sum_fcn_iter::iterator it_old_psf = relevant_partial_sum_functions.begin();
//                 it_old_psf != relevant_partial_sum_functions.end();
//                 ++it_old_psf)          
//         {
//             partial_sum_functions.erase( *it_old_psf );        
//         }
//         
//                         
//                      
//                      
//         std::fprintf(stderr, "\nrank  %d   is waiting at ::broadcast...\n\n", my_MPI_rank);                    
//                 
// 	
//         for (int rank_ctr = 0; rank_ctr < size_MPI_WORLD; ++rank_ctr)
//         {
// 	    
//             Partial_sum_function received_PSF;
// 	    	    
//             if (my_MPI_rank == rank_ctr) 
// 	    {		 
//                 received_PSF = rank_specific___piece_of_Visitation_PSF;
// 	    }	    
//                             
//             boost::mpi::broadcast< Partial_sum_function >
//                                             (*world_MPI_boost_ptr,
//                                             received_PSF,
//                                             rank_ctr); 
// 					          	    
// 	    it_Visit->my_Forward_partial_sum_function
// 			.absorb_piece_of_Partial_Sum_Function_into_one_Partial_Sum_Function(received_PSF);					
//                             
//             world_MPI_boost_ptr->barrier();
// 	    
//             if (my_MPI_rank == rank_ctr) 
// 	    {		 
//                 rank_specific___piece_of_Visitation_PSF.clear_data(); 
// 	    }	                
//         } 	                              
//                   
//                   
//                   
//                                                                  
//         
//     
//         if (my_MPI_rank == 0)	
//             std::fprintf(stderr, "\nit_Visit->my_Forward_partial_sum_function.marginalize_over_eliminating_Event_breakpoints()\n"); 
// 	
// 	it_Visit->my_Forward_partial_sum_function.marginalize_over_eliminating_Event_breakpoints_and_outcomes();  	
//         
//      	
//       
//             
//             
//             
//                 
//         if (my_MPI_rank == 0)  //so multiple ranks don't simultaneously write!!!
//         {
//             std::cerr << "\n\nrank " << my_MPI_rank << "  it_Visit->save_PSF_to_archive()... ";
//             it_Visit->save_PSF_to_archive();  //this makes a lot of room.
//             std::cerr << "it_Visit->save_PSF_to_archive()-----done!\n\n";
//         }                
//                 
// 
//         
//         
//         //housekeeping:
//         if (my_MPI_rank == 0)
//             print_line_of_markers("housekeeping ");                        
//                     
// 
//                     
//         
//         
//         { // copy of PSF
//             partial_sum_functions.push_back(  it_Visit->my_Forward_partial_sum_function  );
//             const type_list_Partial_sum_fcn::iterator  propagating_PSF  = --partial_sum_functions.end();
// 	    	    	    
//             
//             std::cerr << "\n\nrank " << my_MPI_rank  << " remove_eliminated_Event_from_storage___MPI \n";            
//             propagating_PSF->remove_eliminated_Event_from_storage( event_being_summed_out->second.UID );
// 	    
// 
// 	    std::cerr << "\n\nrank " << my_MPI_rank  
// 			<< " make_numerically_stable_by_multiplying_every_partial_sum_value_by_the_inverse_of_the_smallest_partial_sum_value\n ";	    
// 	    propagating_PSF->make_numerically_stable_by_multiplying_every_partial_sum_value_by_the_inverse_of_the_smallest_partial_sum_value();                
// 	    
//                     // We do NOT need to keep track of this normalizing constant!!!!!	            
// 		    
//         } // copy of PSF
//         
//         
//         
//         
//         
//         
//         
//         it_Visit->my_Forward_partial_sum_function.clear_data();
// 
// 
// 
// 
// 
// 
//       
//         if (my_MPI_rank == 0)        
//             std::fprintf(stderr, "debug:  erasing positions from all related events...\n");                           
//         //FINALLY,   ERASE GENOME POSITIONS FROM ALL RELATED EVENTS (SO THE SAME GENOME POSITIONS ARENT CONSIDERED MORE THAN ONCE)
//         //What about genome positions of the current event that weren't on anybody elses homology lists?  (i.e. unique to this LCR-pair)???????                        
//                         
//         type_map_uint_to_Interlock* local_and_global_list[2];
//         local_and_global_list[0] = &(event_being_summed_out->second.local_interlocking_Events);
//         local_and_global_list[1] = &(event_being_summed_out->second.global_interlocking_Events);
// 
//         for (uint lg = 0; lg<2; ++lg)                           
//             for (type_map_uint_to_Interlock::iterator ev_it = local_and_global_list[lg]->begin();
//                     ev_it != local_and_global_list[lg]->end();
//                     ++ev_it)   // for every interlocking event...
//             {
//                 
//                 const type_map_uint_to_Event::iterator it_interlocking_Event =  events.find(ev_it->first);
//                                                                 
//                 //erase yourself from this interlocking event's lists:
//                 it_interlocking_Event->second.local_interlocking_Events.erase(event_being_summed_out->second.UID);
//                 it_interlocking_Event->second.global_interlocking_Events.erase(event_being_summed_out->second.UID);                
//                 
//             
//                     
//                 //and erase your very coordinates.  (necessary, for say, Inbetween region).
//                 if (event_being_summed_out->second.chromos[0] ==  it_interlocking_Event->second.chromos[0])  
//                 {
//                     for (uint pos =  event_being_summed_out->second.region_between_and_including_the_LCRs_themselves.lower();
//                               pos <= event_being_summed_out->second.region_between_and_including_the_LCRs_themselves.upper();
//                               ++pos)                   
//                         it_interlocking_Event->second.remaining_coordinates_for_consideration.erase(pos);               
//                 }
//                 
//                 
//                     
//                 
//                 for (type_list_Coords_and_homologous_region::const_iterator it_homol
//                                 = event_being_summed_out->second.regions_homologous_to_directly_affected_region.begin();
//                         it_homol != event_being_summed_out->second.regions_homologous_to_directly_affected_region.end();
//                         ++it_homol)              
//                     if (it_homol->chromosome_of_homologous_region == it_interlocking_Event->second.chromos[0]
//                         and  BOOST_overlap( get_endpoints_of_keys_of_compressed_map(it_homol->compressed_map_homologous_region_to_profile),
//                                             it_interlocking_Event->second.region_between_and_including_the_LCRs_themselves )      )
//                     {
//                         const type_map_uint_to_uint full_map( convert_compressed_map_to_full_map(it_homol->compressed_map_homologous_region_to_profile) );                        
//                         for (type_map_uint_to_uint::const_iterator map_it = full_map.begin();
//                                 map_it != full_map.end();
//                                 ++map_it)                                                
//                             it_interlocking_Event->second.remaining_coordinates_for_consideration.erase( map_it->first );                                                                      
//                     }  
//             
//             
//             
//                            
//                            
//                 {//homol reg restriction
//                     type_list_Coords_and_homologous_region::iterator it_homol
//                                 = it_interlocking_Event->second.regions_homologous_to_directly_affected_region.begin();                                
//                     while(  it_homol != it_interlocking_Event->second.regions_homologous_to_directly_affected_region.end()  )
//                     {
//                         it_homol->compressed_map_homologous_region_to_profile
//                                 =  get_compressed_map_inverse_of_compressed_map(
//                                         convert_full_map_to_compressed_map(
//                                                 get_map_intersection_of_compressed_map_keys_and_set(
//                                                         get_compressed_map_inverse_of_compressed_map(   it_homol->compressed_map_homologous_region_to_profile    ),
//                                                         it_interlocking_Event->second.remaining_coordinates_for_consideration )  )  ); 
//                                             
//                         if ( it_homol->compressed_map_homologous_region_to_profile.at(0).empty() )
//                             it_homol = it_interlocking_Event->second.regions_homologous_to_directly_affected_region.erase( it_homol );
//                         else
//                             ++it_homol;                                            
//                     }                                                        
//                 }//homol reg restriction                           
//              
// 
//                 
//                     
//                 for (uint qs=0; qs<2; ++qs)                
//                     it_interlocking_Event->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)                        
//                         =  get_compressed_map_inverse_of_compressed_map(
//                                 convert_full_map_to_compressed_map(
//                                     get_map_intersection_of_compressed_map_keys_and_set(
//                                             get_compressed_map_inverse_of_compressed_map(
//                                                     it_interlocking_Event->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)  ),
//                                             it_interlocking_Event->second.variational_positions_of_profile)  ));   
//                                                         
//                             
//                                         
//                 if (my_MPI_rank == 0)  
//                 {                
//                     std::fprintf(stderr, "UID %u", it_interlocking_Event->second.UID);                        
//                     print_set<uint>(it_interlocking_Event->second.variational_positions_of_profile, "it_interlocking_Event->second.variational_positions_of_profile");
//                 }
//                 
//                                                             
//                         
//             }  // end for-loop   ev_it
//                       
//                       
//                       
//                       
//                       
//         
//         // clean up:
//         if (my_MPI_rank == 0)
//             std::fprintf(stderr, "debug:  deleting  rest of  vectors...\n");                           
//         //This variable has been summed out and no longer contains any information.  delete all of its vectors.
//         event_being_summed_out->second.global_interlocking_Events.clear();
//         event_being_summed_out->second.local_interlocking_Events.clear();
//         event_being_summed_out->second.remaining_coordinates_for_consideration.clear();
//         event_being_summed_out->second.original_variational_positions.clear();
//         
//         event_being_summed_out->second.regions_homologous_to_LCR_profile.clear();
//         event_being_summed_out->second.regions_homologous_to_directly_affected_region.clear();
//         event_being_summed_out->second.variational_positions_of_profile.clear();  
//         event_being_summed_out->second.compressed_map_LCR_to_profile__for_each_LCR.clear();
//         event_being_summed_out->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.clear();        
//         
//         
//         
//         
//         
//         
//         
//         
//         if (my_MPI_rank == 0)
//         {
//             print_line_of_markers("*-+-*");
//             print_line_of_markers("Visitation Schedule step complete!");  
//             print_line_of_markers("*-+-*");
//         }
//         
//     }      // end    for-loop    variable_elimination_schedule
// 
// 
// 
// 
//     return 0;  //success


}  //  end of   use_saved_likelihoods_to_perform_complete_variable_elimination_schedule_with_current_prior










