#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <algorithm>

#include <limits>

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


#include <mpreal.h>
#include <mpfr.h>

// #include <gmp-i386.h>


#include <general_typedefs.h>
#include <Sparse_map.h>
#include <templates.h>
#include <translations_through_alignment.h>

#include "globals.h"
#include "Event.h"
#include "other_functions.h"
#include "Partial_sum_function.h"
#include "Conn_Comp.h"





void Partial_sum_function::print_this_PSF__states_only
			    (std::stringstream *const &some_ss) const
{
    std::stringstream output_str;
    
    const bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str; 
    
    
    (*ss_ptr) << "\n\n\n\n\n";
    print_line_of_markers("[[[", ss_ptr);    
    (*ss_ptr) << "printing this PSF - states only!:\n\n";
    
    (*ss_ptr) << "\t\t\teliminating event:  "  <<  eliminating_Event_UID  << "\n";        
    print_set<uint>(haploid_argument_UIDs, "haploid_argument_UIDs", ss_ptr);    
    
    print_map_to_set<uint,uint>(current_relevant_variational_positions_per_Event_in_PSF, "current_relevant_variational_positions_per_Event_in_PSF", ss_ptr, false);
    
    for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_dip_neighb_states
                                =  map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.begin();
            it_dip_neighb_states != map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end();
            ++it_dip_neighb_states)
    {
        (*ss_ptr) << "\n\n\n\t\t\t\tdiploid neighbor states:\n";
        print_map_keys_and_values<uint, uint>(it_dip_neighb_states->first.first, "neighbor states hap 0:", ss_ptr);
        print_map_keys_and_values<uint, uint>(it_dip_neighb_states->first.second, "neighbor states hap 1:", ss_ptr);
	
	(*ss_ptr) << "\n\t\t\t\t\t\tnumber diploid breakpoints for this diploid state:  " << it_dip_neighb_states->second.size() << "\n\n";
    }    
        
    
    
    print_line_of_markers("]]]", ss_ptr);
    (*ss_ptr) << "\n\n\n\n\n";    

    if (!provided_with_valid_sstream)
    {  std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );  }    

}//print_this_PSF__states_only





void Partial_sum_function::print_this_PSF(std::stringstream *const &some_ss) const
{
    
    std::stringstream output_str;
    
    const bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str; 
    
    
    (*ss_ptr) << "\n\n\n\n\n";
    print_line_of_markers("[[[", ss_ptr);
    
    (*ss_ptr) << "printing this PSF:\n\n";
    

    (*ss_ptr) << "\t\t\teliminating event:  "  <<  eliminating_Event_UID  << "\n";
    
    
    print_set<uint>(haploid_argument_UIDs, "haploid_argument_UIDs", ss_ptr);
    
    print_map_to_set<uint,uint>(current_relevant_variational_positions_per_Event_in_PSF, "current_relevant_variational_positions_per_Event_in_PSF", ss_ptr, false);
    
    
    for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_dip_neighb_states
                                =  map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.begin();
            it_dip_neighb_states != map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end();
            ++it_dip_neighb_states)
    {
        (*ss_ptr) << "\n\n\n\t\t\t\tdiploid neighbor states:\n";
        print_map_keys_and_values<uint, uint>(it_dip_neighb_states->first.first, "neighbor states hap 0:", ss_ptr);
        print_map_keys_and_values<uint, uint>(it_dip_neighb_states->first.second, "neighbor states hap 1:", ss_ptr);
    }
    
    
    
    for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_dip_neighb_states
                                =  map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.begin();
            it_dip_neighb_states != map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end();
            ++it_dip_neighb_states)
    {    
        (*ss_ptr) << "\n\n\n\t\t\t\tdiploid neighbor states:\n";
        print_map_keys_and_values<uint, uint>(it_dip_neighb_states->first.first, "neighbor states hap 0:", ss_ptr);
        print_map_keys_and_values<uint, uint>(it_dip_neighb_states->first.second, "neighbor states hap 1:", ss_ptr);
        
        for (type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_dip_neighb_brkpts = it_dip_neighb_states->second.begin();
                it_dip_neighb_brkpts != it_dip_neighb_states->second.end();
                ++it_dip_neighb_brkpts)
        {
            (*ss_ptr) << "\n\t\t\t\t\t\t\t\tdiploid neighbor breakpoints:\n";
            print_map_keys_and_values<uint, BOOST_Interval>(it_dip_neighb_brkpts->first.first, "neighbor breakpoints hap 0:", ss_ptr);
            print_map_keys_and_values<uint, BOOST_Interval>(it_dip_neighb_brkpts->first.second, "neighbor breakpoints hap 1:", ss_ptr); 
            
            
            (*ss_ptr) << "\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t PARTIAL SUM VALUE  =   " << it_dip_neighb_brkpts->second.partial_sum_value_S  << "\n\n\n\n";            
                        
//             for (traceback_data_storage::type_map_eliminating_Event_diploid_outcomes_to_hybrid_factor::const_iterator it_EE_outcome
//                                 = it_dip_neighb_brkpts->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.begin();
//                     it_EE_outcome != it_dip_neighb_brkpts->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.end();
//                     ++it_EE_outcome)
//             {
//                 std::fprintf(stderr, "\n\t\t\t\t\t\t\t\t\t\t\t\tEE outcomes: { %u  &&   %u}\n",
//                              it_EE_outcome->first.first, it_EE_outcome->first.second);
//                              
//                 for (type_map_BI__BI___to__real::const_iterator it_EE_brkpts = )             
//             
//             
//             }                            
        }
    
    
    
    
    }
    

    print_line_of_markers("]]]", ss_ptr);
    (*ss_ptr) << "\n\n\n\n\n";
    

    if (!provided_with_valid_sstream)
    {  std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );  }
    
    
}   //  print_this_PSF














void Partial_sum_function::set_eliminating_Event__and__haploid_argument_UIDs_according_to_relevant_PSFs
                                                ( const uint &in_eliminating_Event_UID,
                                                  const type_list_Partial_sum_fcn_iter &list_of_related_PSF,
                                                  const type_set_uint &haploid_UID_arguments_for_new_PSF)
{

    eliminating_Event_UID = (int)in_eliminating_Event_UID;
    
    haploid_argument_UIDs = haploid_UID_arguments_for_new_PSF;
    
    for (type_list_Partial_sum_fcn_iter::const_iterator it_rel_PSF = list_of_related_PSF.begin();
            it_rel_PSF != list_of_related_PSF.end();
            ++it_rel_PSF)
    {        
        type_set_uint rel_PSF_arguments(    (*it_rel_PSF)->retrieve_haploid_argument_UIDs()    );          
        haploid_argument_UIDs.insert( rel_PSF_arguments.begin(), rel_PSF_arguments.end() );                
    }
    
    haploid_argument_UIDs.erase(eliminating_Event_UID);


}  //  set_eliminating_Event_and_haploid_argument_UIDs

















Partial_sum_function::type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator  
	Partial_sum_function::get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints
				(const Sparse_map &oversized_haploid_sparse_state_vector_0,
				const Sparse_map &oversized_haploid_sparse_state_vector_1,
				const type_map_uint_to_BI &oversized_haploid_breakpoints_0,
				const type_map_uint_to_BI &oversized_haploid_breakpoints_1)  const
{
    //the breakpoints requested should be subsets of the breakpoints in this PSF
    
    const type_map_uint_to_uint argument_state_vector_0(
                            get_map_intersection_of_map_keys_and_set<uint, uint>(oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value,
                                                                                             haploid_argument_UIDs));
    
    const type_map_uint_to_uint argument_state_vector_1(
                            get_map_intersection_of_map_keys_and_set<uint, uint>(oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value,
                                                                                             haploid_argument_UIDs));   

                
    
    const type_map_uint_to_BI argument_haploid_breakpoints_0(
                                    get_map_intersection_of_map_keys_and_set<uint, BOOST_Interval>(oversized_haploid_breakpoints_0,
                                                                                                     haploid_argument_UIDs));
    const type_map_uint_to_BI argument_haploid_breakpoints_1(
                                    get_map_intersection_of_map_keys_and_set<uint, BOOST_Interval>(oversized_haploid_breakpoints_1,
                                                                                                     haploid_argument_UIDs));        

    const type_map_uint_to_BI* argument_haploid_breakpoints[2];
    argument_haploid_breakpoints[0] = &argument_haploid_breakpoints_0;
    argument_haploid_breakpoints[1] = &argument_haploid_breakpoints_1;
    
    
    const Sparse_map* oversized_state_vectors[2];
    oversized_state_vectors[0] = &oversized_haploid_sparse_state_vector_0;
    oversized_state_vectors[1] = &oversized_haploid_sparse_state_vector_1;    
                                                                                                     
                                                                                                     
    //This is an important point.  Since we are not summing our graphical model "directly,, but rather summing it "indirectly."  That is, we are summing out by Events, when we should really be summing out by "duplicons".
    //As such, while two non-interlocking Events may become dependent on each through through a partial sum (during variable eliminiation), their dependence is only on certain of their coordinates.  Thus, when querying a Partial Sum Function for some value, that PSF may not actually contain the COORDINATES we seek.
    
    
    const type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_found_states
		    = map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage
			    .find(type_map_uint_to_uint__2(argument_state_vector_0, argument_state_vector_1));
                      
		if (it_found_states == map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end())  
		{                    
		    std::stringstream error_strm;        
		    error_strm << "unable to find state vector in saved PSF in \"get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints\"!!!\n\n";
		    
		    print_map_keys_and_values<uint, uint>(argument_state_vector_0, "argument_state_vector_0", &error_strm);
		    print_map_keys_and_values<uint, uint>(argument_state_vector_1, "argument_state_vector_1", &error_strm);
		    
		    print_map_keys_and_values<uint, uint>(oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value,
							"oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value", &error_strm);
		    print_map_keys_and_values<uint, uint>(oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value,
							"oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value", &error_strm);                                                              
		    
		    print_this_PSF(&error_strm);
		    error_message(error_strm, true);
		}
                    
       


       

    type_map_uint_to_set_BI__2 all_haploid_brkpts_per_Event;
    initialize_map_with_keys_from_singlecontainer<uint, type_set_BI, type_set_uint>(
		all_haploid_brkpts_per_Event.first, extract_keys_of_map_and_return_as_set<uint, uint>(it_found_states->first.first));
    initialize_map_with_keys_from_singlecontainer<uint, type_set_BI, type_set_uint>(
		all_haploid_brkpts_per_Event.second, extract_keys_of_map_and_return_as_set<uint, uint>(it_found_states->first.second));
    
    for (type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_stored_dip_brkpts = it_found_states->second.begin();
            it_stored_dip_brkpts != it_found_states->second.end();
            ++it_stored_dip_brkpts)    
    {
        for (uint hap=0; hap<2; ++hap)
        {                                                           
            for (type_map_uint_to_BI::const_iterator it_stored_hap_brkpts = pair_at<type_map_uint_to_BI>(it_stored_dip_brkpts->first, hap).begin();
                    it_stored_hap_brkpts != pair_at<type_map_uint_to_BI>(it_stored_dip_brkpts->first, hap).end();  
                    ++it_stored_hap_brkpts)   
	    {  pair_at<type_map_uint_to_set_BI>(all_haploid_brkpts_per_Event,hap).at(it_stored_hap_brkpts->first).insert(it_stored_hap_brkpts->second);  }
	}//hap
    }//dip_brkpts
    
    
    
    //Recall that to take the "Cartesian product" of some PSFs, we first aggregate all of the variational positions contained in each of them, and then we go about making NAHR/GeneConversion breakpoints.
    //This implies that the NAHR breakpoints being search for during the Cartesian product are, by definition, finer (in a topological sense) than ALL of the breakpoints in previous PSFs.
    //Therefore, when we search for sufficient NAHR positions, we should always look for an interval which contains the breakpoint being sought out.  Since it is NAHR, we don't have to worry about a "best/tightest fitting" such interval, since no two NAHR intervals overlap.
    
    //Gene Conversion is much trickier - by having more variational positions in the current PSF compared to all previous ones, then some of the Gene Conversion breakpoint intervals that we searcfh for are finer than in any previous PSF, and some of the intervals are coarser. 
    
    
    type_map_uint_to_BI__2 best_haploid_brkpts_per_Event;
    initialize_map_with_keys_from_singlecontainer<uint, BOOST_Interval, type_set_uint>(
		best_haploid_brkpts_per_Event.first, extract_keys_of_map_and_return_as_set<uint, uint>(it_found_states->first.first), empty_BI);
    initialize_map_with_keys_from_singlecontainer<uint, BOOST_Interval, type_set_uint>(
		best_haploid_brkpts_per_Event.second, extract_keys_of_map_and_return_as_set<uint, uint>(it_found_states->first.second), empty_BI);
		
    for (uint hap=0; hap<2; ++hap)
    {     
	for (type_map_uint_to_BI::iterator it_best = pair_at<type_map_uint_to_BI>(best_haploid_brkpts_per_Event, hap).begin();
		it_best != pair_at<type_map_uint_to_BI>(best_haploid_brkpts_per_Event, hap).end();
		++it_best)
	{		    
	    const bool is_NAHR = test_if_haploid_outcome_is_NAHR(static_cast<type_haploid_outcome>(oversized_state_vectors[hap]->get_state_for_UID(it_best->first)));
				    
	    if (is_NAHR)
	    {//NAHR	 
		//This is the finest-grained search.  We want anything which contains this point
		for (type_set_BI::const_iterator it_all_bi = pair_at<type_map_uint_to_set_BI>(all_haploid_brkpts_per_Event,hap).at(it_best->first).begin();
			it_all_bi != pair_at<type_map_uint_to_set_BI>(all_haploid_brkpts_per_Event,hap).at(it_best->first).end();
			++it_all_bi)
		{
		    if (BOOST_subset(argument_haploid_breakpoints[hap]->at(it_best->first), *it_all_bi))
		    {
			it_best->second = *it_all_bi;
			break;
		    }//subset			   
		}//all bi
	    }//NAHR
	    else
	    {//GeneConv
		uint best__minimum_distance_between_superset_brkpts = std::numeric_limits<int>::max(); 
		for (type_set_BI::const_iterator it_all_bi = pair_at<type_map_uint_to_set_BI>(all_haploid_brkpts_per_Event,hap).at(it_best->first).begin();
			it_all_bi != pair_at<type_map_uint_to_set_BI>(all_haploid_brkpts_per_Event,hap).at(it_best->first).end();
			++it_all_bi)
		{
		    if (BOOST_subset(argument_haploid_breakpoints[hap]->at(it_best->first), *it_all_bi))
		    {
			const uint set_difference_size = BOOST_width_inclusive(*it_all_bi) - BOOST_width_inclusive(argument_haploid_breakpoints[hap]->at(it_best->first));
			if (set_difference_size < best__minimum_distance_between_superset_brkpts)
			{
			    best__minimum_distance_between_superset_brkpts = set_difference_size;
			    it_best->second = *it_all_bi;	
			    
			    if (best__minimum_distance_between_superset_brkpts == 0)
			    {  break;  }
			}//diff
		    }//subset			   
		}//all bi		
	    }//GeneConv
	    
	    		
	    
	    if (BOOST_empty(it_best->second)  and  !pair_at<type_map_uint_to_set_BI>(all_haploid_brkpts_per_Event,hap).at(it_best->first).empty())		
	    {
		std::stringstream error_strm;
		error_strm << "\nin \"get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints\", unable to find acceptable diploid brkpts\n"
			<< "\t\tfor Event = " << it_best->first
			<< "\n\thap = " << hap 
			<< "\n\tsearched-interval = argument_haploid_breakpoints[hap]->at(it_best->first) = " << argument_haploid_breakpoints[hap]->at(it_best->first)
			<<  "\n\n";
		
		print_map_keys_and_values<uint, BOOST_Interval>(*argument_haploid_breakpoints[hap], "*argument_haploid_breakpoints[hap]", &error_strm, true);
								
		print_set<BOOST_Interval>(pair_at<type_map_uint_to_set_BI>(all_haploid_brkpts_per_Event, hap).at(it_best->first), 
						"pair_at<type_map_uint_to_set_BI>(all_haploid_brkpts_per_Event, hap).at(it_best->first)",
						&error_strm, true);
			    
		print_map_keys_and_values<uint,uint>(oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value, "oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value", &error_strm);
		print_map_keys_and_values<uint,uint>(oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value, "oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value", &error_strm);
		print_map_keys_and_values<uint, BOOST_Interval>( oversized_haploid_breakpoints_0, "oversized_haploid_breakpoints_0", &error_strm);
		print_map_keys_and_values<uint, BOOST_Interval>( oversized_haploid_breakpoints_1, "oversized_haploid_breakpoints_1", &error_strm);
		
		print_this_PSF(&error_strm);     
		
		error_message(error_strm, true);
	    }//error-check	    
	}//ev			
    }//hap
    
    
    
    
    
    
    
    const type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_found_best_brkpts = it_found_states->second.find(best_haploid_brkpts_per_Event);
    
	    if (it_found_best_brkpts == it_found_states->second.end())
	    {
		//if we made it this far, then everything is messed up: 
		std::stringstream error_strm;		
		error_strm << "\nin \"get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints\", determine best set of superset breakpoints, but unable to find such a combination in the existing map.\n\n"
		<< "\nbest_haploid_brkpts_per_Event, by haploid:\n\n"; 
		
		print_map_keys_and_values<uint, BOOST_Interval>(best_haploid_brkpts_per_Event.first, "best_haploid_brkpts_per_Event.first", &error_strm, true);
		print_map_keys_and_values<uint, BOOST_Interval>(best_haploid_brkpts_per_Event.second, "best_haploid_brkpts_per_Event.second", &error_strm, true);	
			
		print_map_keys_and_values<uint,uint>(oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value, "oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value", &error_strm);
		print_map_keys_and_values<uint,uint>(oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value, "oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value", &error_strm);
		print_map_keys_and_values<uint, BOOST_Interval>( oversized_haploid_breakpoints_0, "oversized_haploid_breakpoints_0", &error_strm);
		print_map_keys_and_values<uint, BOOST_Interval>( oversized_haploid_breakpoints_1, "oversized_haploid_breakpoints_1", &error_strm);
		
		print_this_PSF(&error_strm);   
		
		error_message(error_strm, true);    
	    }    
    
    
    return  it_found_best_brkpts;

}//get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints





real Partial_sum_function::get_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints
                                        ( const Sparse_map &oversized_haploid_sparse_state_vector_0,
                                          const Sparse_map &oversized_haploid_sparse_state_vector_1,
                                          const type_map_uint_to_BI &oversized_haploid_breakpoints_0,
                                          const type_map_uint_to_BI &oversized_haploid_breakpoints_1 ) const
{    
    return  get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints
		    (oversized_haploid_sparse_state_vector_0,
		     oversized_haploid_sparse_state_vector_1,
		     oversized_haploid_breakpoints_0,
		     oversized_haploid_breakpoints_1)->second.partial_sum_value_S;
					  
}//get_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints








void  Partial_sum_function::allocate_diploid_breakpoints_for_given_diploid_states
		    ( const Sparse_map &oversized_haploid_sparse_state_vector_0,
		      const Sparse_map &oversized_haploid_sparse_state_vector_1,
		      const type_map_uint_to_vector_BI  &breakpoint_space__hap_0,
		      const type_map_uint_to_vector_BI  &breakpoint_space__hap_1)		      
{
    
    const type_map_uint_to_uint argument_state_vector_0(
                            get_map_intersection_of_map_keys_and_set<uint, uint>(oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value,
                                                                                             haploid_argument_UIDs));
    
    const type_map_uint_to_uint argument_state_vector_1(
                            get_map_intersection_of_map_keys_and_set<uint, uint>(oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value,
                                                                                             haploid_argument_UIDs));              
					
				
    const type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::iterator  it_insert_diploid_states
		    =  map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage
			    .insert(std::pair<type_map_uint_to_uint__2, type_map_diploid_neighbors_brkpts_to_traceback_data_storage>
				    (type_map_uint_to_uint__2(argument_state_vector_0,argument_state_vector_1),
				     type_map_diploid_neighbors_brkpts_to_traceback_data_storage() )  ).first;
				
				     
    uint number_of_breakpoints_0 = 1;
    if (!argument_state_vector_0.empty())
    {  number_of_breakpoints_0 = breakpoint_space__hap_0.begin()->second.size();  }
    
    uint number_of_breakpoints_1 = 1;
    if (!argument_state_vector_1.empty())
    {  number_of_breakpoints_1 = breakpoint_space__hap_1.begin()->second.size();  }
    
    
    
    type_map_uint_to_BI  hap_brkpts_0;
    type_map_uint_to_BI  hap_brkpts_1;    
    initialize_map_with_keys_from_another_map<uint, BOOST_Interval, type_vector_BI>(hap_brkpts_0, breakpoint_space__hap_0, empty_BI);
    initialize_map_with_keys_from_another_map<uint, BOOST_Interval, type_vector_BI>(hap_brkpts_1, breakpoint_space__hap_1, empty_BI);    
    
    type_map_diploid_neighbors_brkpts_to_traceback_data_storage::iterator insert_brkpts_it = it_insert_diploid_states->second.begin();
    for (uint brkpt_ctr_0 = 0; brkpt_ctr_0 < number_of_breakpoints_0; ++brkpt_ctr_0)
    {	
	set_breakpoints_using_breakpoint_space(hap_brkpts_0, breakpoint_space__hap_0, brkpt_ctr_0);	
	
	for (uint brkpt_ctr_1 = 0; brkpt_ctr_1 < number_of_breakpoints_1; ++brkpt_ctr_1)
	{	    
	    set_breakpoints_using_breakpoint_space(hap_brkpts_1, breakpoint_space__hap_1, brkpt_ctr_1);	    
	    
	    insert_brkpts_it = it_insert_diploid_states->second.insert(
						insert_brkpts_it,
						std::pair<type_map_uint_to_BI__2, traceback_data_storage>
						    (type_map_uint_to_BI__2(hap_brkpts_0, hap_brkpts_1), traceback_data_storage()));	    	    	    
	}//brkpt_ctr_1
    }//brkpt_ctr_0        
    
        
}//allocate_diploid_breakpoints_for_given_diploid_states












void Partial_sum_function::save_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints
                                        ( const Sparse_map &oversized_haploid_sparse_state_vector_0,
                                          const Sparse_map &oversized_haploid_sparse_state_vector_1,
                                          const type_map_uint_to_BI &oversized_haploid_breakpoints_0,
                                          const type_map_uint_to_BI &oversized_haploid_breakpoints_1,
                                          const real &contribution_to_partial_sum_value_S)
{ 
    
    const type_map_uint_to_uint argument_state_vector_0(
                            get_map_intersection_of_map_keys_and_set<uint, uint>(oversized_haploid_sparse_state_vector_0.sparse_map_UID_to_value,
                                                                                             haploid_argument_UIDs));
    
    const type_map_uint_to_uint argument_state_vector_1(
                            get_map_intersection_of_map_keys_and_set<uint, uint>(oversized_haploid_sparse_state_vector_1.sparse_map_UID_to_value,
                                                                                             haploid_argument_UIDs));                   
    
                                                                                             
                                                                                             
    const type_map_uint_to_BI argument_haploid_breakpoints_0(
                                    get_map_intersection_of_map_keys_and_set<uint, BOOST_Interval>(oversized_haploid_breakpoints_0,
                                                                                                     haploid_argument_UIDs));
    const type_map_uint_to_BI argument_haploid_breakpoints_1(
                                    get_map_intersection_of_map_keys_and_set<uint, BOOST_Interval>(oversized_haploid_breakpoints_1,
                                                                                                     haploid_argument_UIDs));        
                                             					     
												     
//     if (map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage
// 			    .count(type_map_uint_to_uint__2(argument_state_vector_0, argument_state_vector_1)) == 0
// 	or   map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage
// 			    .at(type_map_uint_to_uint__2(argument_state_vector_0, argument_state_vector_1))
// 			    .count(type_map_uint_to_BI__2(argument_haploid_breakpoints_0, argument_haploid_breakpoints_1))  == 0)						
//     {
// 	std::cerr << "\n\nsave_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints DNE  !!! ";
// 	print_map_keys_and_values<uint,uint>(argument_state_vector_0,"argument_state_vector_0");
// 	print_map_keys_and_values<uint,uint>(argument_state_vector_1,"argument_state_vector_0");
// 	print_map_keys_and_values<uint,BOOST_Interval>(argument_haploid_breakpoints_0,"argument_haploid_breakpoints_0");
// 	print_map_keys_and_values<uint,BOOST_Interval>(argument_haploid_breakpoints_1,"argument_haploid_breakpoints_1");
// 	print_this_PSF();
//     }
												     
    map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage
		    .at(type_map_uint_to_uint__2(argument_state_vector_0, argument_state_vector_1))
			    .at(type_map_uint_to_BI__2(argument_haploid_breakpoints_0, argument_haploid_breakpoints_1))
				    .partial_sum_value_S = contribution_to_partial_sum_value_S;
                         
} // end of   save_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints
































void Partial_sum_function::marginalize_over_eliminating_Event_breakpoints_and_outcomes()
{                
    
    for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::iterator it_diploid_neighbor_states 
					= map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.begin();
	    it_diploid_neighbor_states != map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end();
	    ++it_diploid_neighbor_states)
    {
	const std::vector<type_map_diploid_neighbors_brkpts_to_traceback_data_storage::iterator>
		diploid_brkpt___loop(get_loopable_iterators<type_map_diploid_neighbors_brkpts_to_traceback_data_storage>(
						it_diploid_neighbor_states->second));
						
	#pragma omp parallel
	{//parallel
	
	    #pragma omp for schedule(static)
	    for (uint brkpt_ctr = 0; brkpt_ctr < diploid_brkpt___loop.size(); ++brkpt_ctr)
	    {				
                type_list_real EE_outcome_probabilities;            
                            
                for (traceback_data_storage::type_map_eliminating_Event_diploid_outcomes_to_hybrid_factor::iterator it_EE_outcome
				    = diploid_brkpt___loop[brkpt_ctr]->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.begin();
                        it_EE_outcome != diploid_brkpt___loop[brkpt_ctr]->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.end();
                        ++it_EE_outcome)   
                {
                                                    
                    type_list_real E_brkpts_probabilities(
                                    extract_values_of_map_and_return_as_list<type_BI__BI, real>(
                                            it_EE_outcome->second.
                                map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts)  );
                                                                                                    
                                
                    E_brkpts_probabilities.sort();  //necessary because of underflow problems.
                                    
                    const real EE_brkpt_sum(sum_over_list<real>(E_brkpts_probabilities));
                                    
                    
                    it_EE_outcome->second.P_D__cond__diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states = EE_brkpt_sum;

                    EE_outcome_probabilities.push_back(EE_brkpt_sum);
                }//EE_outcome     
                
                                        
                EE_outcome_probabilities.sort();  
                
                const real EE_outcome_sum(sum_over_list<real>(EE_outcome_probabilities));
                
                diploid_brkpt___loop[brkpt_ctr]->second.partial_sum_value_S = EE_outcome_sum;             
            }//neighbor brkpts                 
            
            mpfr_free_cache();
        }//parallel                
    }//neighbor states
    
    mpfr_free_cache();

}  //   marginalize_over_eliminating_Event_breakpoints_and_outcomes






















bool Partial_sum_function::check_if_event_is_an_argument_to_this_function(const uint &test_UID) const
{
  
    return (bool)haploid_argument_UIDs.count(test_UID);  

}   //  end of   check_if_event_is_an_argument_to_this_function






type_set_uint Partial_sum_function::retrieve_haploid_argument_UIDs() const
{
    
    return haploid_argument_UIDs;  

}  //  end of   retrieve_haploid_argument_UIDs




                                              
            
            
            
            
            
            
            
            
            

            
// typedef   std::pair<type_haploid_outcome, BOOST_Interval>  type_hapoutcome__BI;
// typedef   std::map<uint, type_hapoutcome__BI>   type_map_uint_to__hapoutcome__BI;
// typedef   std::pair<type_map_uint_to__hapoutcome__BI,type_map_uint_to__hapoutcome__BI>   type_map_uint_to__hapoutcome__BI__2;
                                                
type_multimap_real_to_uint__uint Partial_sum_function::calculate_and_return_eliminated_Event_diploid_outcome_CDF
                                      ( const type_map_uint_to_Sampled_diploid_Event_data &map_estimate_for_this_CC)  const
{
            
    type_map_uint_to_uint general_state_vectors[2];         
    type_map_uint_to_BI general_breakpoint_vectors[2];
    type_map_uint_to_uint argument_state_vector[2];
    type_map_uint_to_BI argument_haploid_breakpoints[2];

    for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_ev = map_estimate_for_this_CC.begin();
		it_ev != map_estimate_for_this_CC.end();
		++it_ev)                       
    {
	for (uint hap=0; hap<2; ++hap)
	{
	    if (pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap) != hap_outcome__None)
	    {	   
		general_state_vectors[hap].insert(type_uint__uint(it_ev->first, (uint)pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap)));
		general_breakpoint_vectors[hap].insert(type_uint__BI(it_ev->first, pair_at<BOOST_Interval>(it_ev->second.the_diploid_profile_brkpts, hap)));
	    }
	}    
    }//ev
    
    

    for (uint hap=0; hap<2; ++hap)
    {
	argument_state_vector[hap] = get_map_intersection_of_map_keys_and_set<uint, uint>(general_state_vectors[hap],
											    haploid_argument_UIDs);	
	
	argument_haploid_breakpoints[hap] = get_map_intersection_of_map_keys_and_set<uint, BOOST_Interval>(general_breakpoint_vectors[hap],
													    haploid_argument_UIDs);
													
    }//hap
  
                    
    
    const type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator  it_requested_diploid_states
		    =  map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.find(  
						    type_map_uint_to_uint__2(argument_state_vector[0], argument_state_vector[1])  );
												  		    
    if (it_requested_diploid_states == map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end())
    {
	std::stringstream error_strm;
	
	error_strm << "\n\nERROR!   in  \"calculate_and_return_eliminated_Event_diploid_outcome_CDF\" " 
		    << "    unable to find  argument state vectors in \"map_diploid_neighbor_states_to_filename_UUID\"\n\n";
		    
	print_map_keys_and_values<uint,uint>(argument_state_vector[0], "argument_state_vector[0]", &error_strm, true);
	print_map_keys_and_values<uint,uint>(argument_state_vector[1], "argument_state_vector[1]", &error_strm, true);
	
	error_message(error_strm, true);
    }

               

    
   
   
    const type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_requested_diploid_brkpts
		= get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints(
						    general_state_vectors[0],
						    general_state_vectors[1],
						    general_breakpoint_vectors[0],
						    general_breakpoint_vectors[1]);
    
    if (it_requested_diploid_brkpts == it_requested_diploid_states->second.end())
    {
        std::stringstream error_strm;
        print_line_of_markers("ERROR! ", &error_strm);
        print_line_of_markers("(", &error_strm);
        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
        
        error_strm << "\n\nERROR!    requested_traceback_storage  == NULL.   failed to find sufficient breakpoints  in  \"calculate_and_return_eliminated_Event_diploid_outcome_CDF\"!!!\n\n";
        
        print_this_PSF( &error_strm );
        
        print_line_of_markers(")", &error_strm);
        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );        
    }
            
            
            
            
        
    
    
                  
    type_multimap_real_to_uint__uint EE_state_probabilities;            

    
    for (traceback_data_storage::type_map_eliminating_Event_diploid_outcomes_to_hybrid_factor::const_iterator it_EE_outcomes
                                            = it_requested_diploid_brkpts->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.begin();
            it_EE_outcomes != it_requested_diploid_brkpts->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.end();
            ++it_EE_outcomes)
    {
        EE_state_probabilities.insert( std::pair<real, type_uint__uint> 
                                        (it_EE_outcomes->second.P_D__cond__diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states,
                                        it_EE_outcomes->first));
    }
                                        



    type_multimap_real_to_uint__uint  EE_state_CDF;
    type_multimap_real_to_uint__uint::iterator insert_it_CDF = EE_state_CDF.begin();
    real running_sum( 0.00 );  
    
    for (type_multimap_real_to_uint__uint::const_iterator it_probs = EE_state_probabilities.begin();
            it_probs != EE_state_probabilities.end();
            ++it_probs)
    {
        running_sum += (it_probs->first / it_requested_diploid_brkpts->second.partial_sum_value_S);
        insert_it_CDF = EE_state_CDF.insert(insert_it_CDF, std::pair<real, type_uint__uint>
                                                            (running_sum, it_probs->second));                                                 
    }                                            



    if ( std::abs( (running_sum.toLDouble() - 1.00L) ) >  0.00001L) 
    {
        std::stringstream error_strm;
        print_line_of_markers("ERROR! ", &error_strm);
        print_line_of_markers("(", &error_strm);        
        error_strm  << "\n\nWARNING: outcome_CDF running_sum = " << running_sum <<"  is far from 1.00  in \"calculate_and_return_eliminated_Event_diploid_outcome_CDF\" !!!\n\n\n";
	
	print_this_PSF( &error_strm );
	
        print_line_of_markers(")", &error_strm);
        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
    }
    
     
    mpfr_free_cache(); 
            
    return EE_state_CDF;

} // end of calculate_and_return_eliminated_Event_diploid_outcome_CDF
            
            
            
            
            
            
            
    
    
    
    
    
    
    
    
    
    
    
    
    
type_multimap_real_to__BI__BI Partial_sum_function::calculate_and_return_eliminated_Event_diploid_breakpoint_CDF
                                      ( const type_map_uint_to_Sampled_diploid_Event_data &map_estimate_for_this_CC )  const
{
    
    type_map_uint_to_uint general_state_vectors[2];         
    type_map_uint_to_BI general_breakpoint_vectors[2];
    type_map_uint_to_uint argument_state_vector[2];
    type_map_uint_to_BI argument_haploid_breakpoints[2];

    for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_ev = map_estimate_for_this_CC.begin();
		it_ev != map_estimate_for_this_CC.end();
		++it_ev)                       
    {
	for (uint hap=0; hap<2; ++hap)
	{
	    if (pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap) != hap_outcome__None)
	    {	   
		general_state_vectors[hap].insert(type_uint__uint(it_ev->first, (uint)pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap)));
		general_breakpoint_vectors[hap].insert(type_uint__BI(it_ev->first, pair_at<BOOST_Interval>(it_ev->second.the_diploid_profile_brkpts, hap)));
	    }
	}    
    }//ev
    
    

    for (uint hap=0; hap<2; ++hap)
    {
	argument_state_vector[hap] = get_map_intersection_of_map_keys_and_set<uint, uint>( general_state_vectors[hap],
											    haploid_argument_UIDs );	
	
	argument_haploid_breakpoints[hap] = get_map_intersection_of_map_keys_and_set<uint, BOOST_Interval>( general_breakpoint_vectors[hap],
													    haploid_argument_UIDs);
													
    }//hap
                    
    
    const type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator  it_requested_diploid_states
		    =  map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.find(  
						    type_map_uint_to_uint__2(argument_state_vector[0], argument_state_vector[1])  );
												  		    
    if (it_requested_diploid_states == map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end())
    {
	std::stringstream error_strm;
	
	error_strm << "\n\nERROR!   in  \"calculate_and_return_eliminated_Event_diploid_outcome_CDF\" " 
		    << "    unable to find  argument state vectors in \"map_diploid_neighbor_states_to_filename_UUID\"\n\n";
		    
	print_map_keys_and_values<uint,uint>(argument_state_vector[0], "argument_state_vector[0]", &error_strm, true);
	print_map_keys_and_values<uint,uint>(argument_state_vector[1], "argument_state_vector[1]", &error_strm, true);
	
	error_message(error_strm, true);
    }

               

    
   
   
    const type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_requested_diploid_brkpts
		= get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints(
						    general_state_vectors[0],
						    general_state_vectors[1],
						    general_breakpoint_vectors[0],
						    general_breakpoint_vectors[1]);									          
    
    if (it_requested_diploid_brkpts == it_requested_diploid_states->second.end())
    {
        std::stringstream error_strm;
        print_line_of_markers("ERROR! ", &error_strm);
        print_line_of_markers("(", &error_strm);
        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
        
        error_strm << "\n\nERROR!    requested_traceback_storage  == NULL.   failed to find sufficient breakpoints  in  \"calculate_and_return_eliminated_Event_diploid_outcome_CDF\"!!!\n\n";
        
        print_this_PSF( &error_strm );
        
        print_line_of_markers(")", &error_strm);
        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );        
    }    
    
    
        
  
         
    

    
    
 
 
    
    
    const type_uint__uint eliminated_Event_diploid_outcome( map_estimate_for_this_CC.at(eliminating_Event_UID).the_diploid_Event_outcome.first,
							    map_estimate_for_this_CC.at(eliminating_Event_UID).the_diploid_Event_outcome.second);
    
    const Partial_sum_function::traceback_data_storage::type_map_eliminating_Event_diploid_outcomes_to_hybrid_factor::const_iterator  it_requested_hybrid_factor
                =  it_requested_diploid_brkpts->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.find(  eliminated_Event_diploid_outcome  );
    

    
  
    real running_sum( 0.00L );
        

                 
    //first convert to multi-map for summing (underflow problems).
    type_multimap_real_to__BI__BI  sorted_EE_brkpts;
    
    
    
    for (type_map_BI__BI___to__real::const_iterator it_EE_brkpts  = it_requested_hybrid_factor->second.
                    map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts.begin();
            it_EE_brkpts  != it_requested_hybrid_factor->second.
                    map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts.end();
            ++it_EE_brkpts)  
    {
        sorted_EE_brkpts.insert( std::pair<real, type_BI__BI>
                                  (it_EE_brkpts->second,   type_BI__BI( it_EE_brkpts->first.first,
                                                                            it_EE_brkpts->first.second   )     ) );
    }
    
               
    
    //now sum over multimap
    type_multimap_real_to__BI__BI EE_breakpoint_CDF;
    type_multimap_real_to__BI__BI::iterator it_insert_EE_brkpt_CDF = EE_breakpoint_CDF.begin();        
    
    for (type_multimap_real_to__BI__BI::const_iterator it_sorted_EE_brkpts = sorted_EE_brkpts.begin();
            it_sorted_EE_brkpts != sorted_EE_brkpts.end();
            ++it_sorted_EE_brkpts)
    {
        running_sum += (   it_sorted_EE_brkpts->first
                           / it_requested_hybrid_factor->second.P_D__cond__diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states   );
        it_insert_EE_brkpt_CDF = EE_breakpoint_CDF.insert( it_insert_EE_brkpt_CDF,   std::pair<real, type_BI__BI>
                                                                                   ( running_sum, it_sorted_EE_brkpts->second )    );    
    }                                
                 



          
                
                
    if ( std::abs( (running_sum.toLDouble() - 1.00L) ) >  0.00001L ) 
    {
        std::stringstream error_strm;
        print_line_of_markers("ERROR! ", &error_strm);
        print_line_of_markers("(", &error_strm);        
        error_strm  << "\n\nWARNING: outcome_CDF running_sum = " << running_sum <<"  is far from 1.00   in \"calculate_and_return_eliminated_Event_diploid_breakpoint_CDF\" !!!\n\n\n";
	
	print_this_PSF( &error_strm );
	
        print_line_of_markers(")", &error_strm);
        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
    }
                
                     
    mpfr_free_cache();                   
                     
    return EE_breakpoint_CDF;

}      // end of   calculate_and_return_eliminated_Event_diploid_breakpoint_CDF
    
    
    
    
    
    
    
            































void Partial_sum_function::save_Cartesian_product_of_related_PSFs_wrt_eliminating_Event_and_include_probabilities_of_Breakpoint_and_Event
                            (const Event  &eliminating_Event,
			     const Sparse_map &full_haploid_sparse_state_vector_0,
                            const Sparse_map &full_haploid_sparse_state_vector_1,                                           
                            const type_list_Partial_sum_fcn_iter &list_of_related_PSF,
                            const Partial_sum_function &thread_specific____new_PSF)                            
{
    
    //the following conditionals in setting up are for the case that there are no remaining positions for consideration, in whcih case "thread_specific____new_PSF" will be empty and thus should be ignored.       
    
    const uint total_number_of_related_PSF
	=  thread_specific____new_PSF.map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage
					.count(type_map_uint_to_uint__2(full_haploid_sparse_state_vector_0.sparse_map_UID_to_value,
									full_haploid_sparse_state_vector_1.sparse_map_UID_to_value))  == 0  ?
			    list_of_related_PSF.size()    :   (list_of_related_PSF.size() + 1);

    std::vector<const Partial_sum_function*>  list_of_related_PSF__and__thread_specific____new_PSF;
    list_of_related_PSF__and__thread_specific____new_PSF.reserve(total_number_of_related_PSF);    
    
    {//scope
        for (type_list_Partial_sum_fcn_iter::const_iterator it_PSF = list_of_related_PSF.begin();
                it_PSF != list_of_related_PSF.end();
                ++it_PSF)
	{  list_of_related_PSF__and__thread_specific____new_PSF.push_back(&(**it_PSF));  }
        
        if (total_number_of_related_PSF == (list_of_related_PSF.size() + 1))
	{  list_of_related_PSF__and__thread_specific____new_PSF.push_back(&thread_specific____new_PSF);  }
	
	assert(list_of_related_PSF__and__thread_specific____new_PSF.size() == total_number_of_related_PSF);
    }//scope
            
                          
    
    
    
    
    type_map_uint_to_Breakpoint_complex total_Breakpoint_complex;    
    {//total varpos    
    
	type_map_uint_to_set_uint total_map_Event_to_relevant_varpos;
    
	for (uint l_PSF = 0; l_PSF < total_number_of_related_PSF; ++l_PSF)        
	{
	    total_map_Event_to_relevant_varpos
		    = combine_two__maps_to_set<uint,uint>(total_map_Event_to_relevant_varpos,
							    list_of_related_PSF__and__thread_specific____new_PSF.at(l_PSF)->current_relevant_variational_positions_per_Event_in_PSF);
	}//other PSFs
	
	
	for (type_map_uint_to_set_uint::const_iterator it_ev_rel = total_map_Event_to_relevant_varpos.begin();
		it_ev_rel != total_map_Event_to_relevant_varpos.end();
		++it_ev_rel)
	{  total_Breakpoint_complex[it_ev_rel->first].variational_positions__ALL = it_ev_rel->second;  }
		
	
	determine_relevant_NAHR_and_GeneConversion_breakpoints_for_each_interlocking_event_and_self_according_to_relevant_variational_positions(
									*eliminating_Event.my_Conn_Comp,
									total_map_Event_to_relevant_varpos,
									total_Breakpoint_complex);
	
	//record it for yourself too:
	record_additional_relevant_variational_positions_per_Event(total_map_Event_to_relevant_varpos);
	
    }//total varpos
  
    
                         
                    
                    
     //NOTE FOR ANY SET OF STATE VECTORS, THERE WILL ALWAYS BE ~FOUR STATE VECTORS   WHO HAVE ALL DIPLOID STATES EQUAL EXCEPT FOR THE (DIPLOID) VALUE OF THE ELIMINATING EVENT, WHICH WILL CAN RANGE FROM (0,0)...(1,1).  BUT ALL OF THESE STATE VECTORS WILL BE GROUPED UNDER THE SAME NEIGHBOR DIPLOIID STATE, I.E. HAVE THE SAME FILE UUID.
                    
                    
                    
            
                    
    const uint EE_UID = eliminating_Event_UID;          
        
    type_map_diploid_neighbors_brkpts_to_traceback_data_storage diploid_breakpoints_to_data;
    type_map_uint_to_uint__2 the_diploid_states;        
    {//set states    
        type_map_uint_to_uint non_EE_full_haploid_state_vector_0(full_haploid_sparse_state_vector_0.sparse_map_UID_to_value);
        type_map_uint_to_uint non_EE_full_haploid_state_vector_1(full_haploid_sparse_state_vector_1.sparse_map_UID_to_value);
        
        non_EE_full_haploid_state_vector_0.erase(EE_UID);    
        non_EE_full_haploid_state_vector_1.erase(EE_UID);  
        
        the_diploid_states = type_map_uint_to_uint__2(non_EE_full_haploid_state_vector_0, non_EE_full_haploid_state_vector_1);
        
        const type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_find
				    = map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.find(the_diploid_states);
        if (it_find != map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end())
        {
	    diploid_breakpoints_to_data = it_find->second;
        }                                 
    }//set states
    
    
    
    
                
    const type_uint__uint diploid_EE_state(full_haploid_sparse_state_vector_0.get_state_for_UID(EE_UID),
                                            full_haploid_sparse_state_vector_1.get_state_for_UID(EE_UID));                                       
    
                                            
                                            
                                            
    //Here we  incorporate prior probabiltiies   P( E | theta)  and   P( s | E,theta )  , but only for the eliminating event.   
    
    
    real P__E_diploid__cond_theta(1.00L); 
    
    for (uint hap=0; hap < 2; ++hap)
    {
	switch( static_cast<type_haploid_outcome>(pair_at<uint>(diploid_EE_state, hap)))
	{
	    case hap_outcome__None:
		P__E_diploid__cond_theta *= P_E_absent_cond_theta;
		break;
	    case hap_outcome__Del:
	    case hap_outcome__Dup:
		P__E_diploid__cond_theta *= P_E_dup_or_del_cond_theta;
		break;
	    case hap_outcome__Inv:
		P__E_diploid__cond_theta *= P_E_inv_cond_theta;
		break;
	    case hap_outcome__GeneConv_ABA:
	    case hap_outcome__GeneConv_BAB:
		P__E_diploid__cond_theta *= P_E_GeneConv_onetype_cond_theta;
		break;
	    default:
		std::cerr << "\n\n\n\n\n\tERROR!  this hap outcome = " << pair_at<uint>(diploid_EE_state, hap) << "  shouldn't be considered yet!!!\n\n";
		exit(1);
		break;
	};//switch
    }//hap
    


    
    
    
    type_real__real  P__s_haploid__cond__E_haploid_and_theta(1,1);
    for (uint hap=0; hap<2; ++hap)
    {	
	switch(static_cast<type_haploid_outcome>(pair_at<uint>(diploid_EE_state,hap)))
	{
	    case hap_outcome__None:
		//do nothing.  already = 1.
		break;
	    case hap_outcome__Del:
		pair_at<real>(P__s_haploid__cond__E_haploid_and_theta,hap) /= eliminating_Event.my_original_Breakpoint_complex.breakpoints_NAHR__AB.size();
		break;
	    case hap_outcome__Dup:
		pair_at<real>(P__s_haploid__cond__E_haploid_and_theta,hap) /= eliminating_Event.my_original_Breakpoint_complex.breakpoints_NAHR__BA.size();
		break;
	    case hap_outcome__Inv:
		pair_at<real>(P__s_haploid__cond__E_haploid_and_theta,hap) /= eliminating_Event.my_original_Breakpoint_complex.breakpoints_NAHR__AB.size(); // either one
		break;
	    case hap_outcome__GeneConv_ABA:
		pair_at<real>(P__s_haploid__cond__E_haploid_and_theta,hap) /= eliminating_Event.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA.size();
		break;
	    case hap_outcome__GeneConv_BAB:
		pair_at<real>(P__s_haploid__cond__E_haploid_and_theta,hap) /= eliminating_Event.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB.size();
		break;
	    default:
		std::stringstream error_strm;
		error_strm << "ERROR!   unrecognized hapoutcome = " << convert_haploid_outcome_to_string(static_cast<type_haploid_outcome>(pair_at<uint>(diploid_EE_state,hap)))
			    << "  in  \"type_real__real  P__s_haploid__cond__E_haploid_and_theta(1,1);\" part of PSFCartesian product.\n";
		error_message(error_strm, true);
		break;
	}	
    }//hap
    
    

    
    //diploid___grouped__NAHR__GeneConv;//[haploid][NAHR or GeneConv]             
    //sum over breakpoints:  
    type_map_uint_to_vector_BI  breakpoint_space__hap_0;
    type_map_uint_to_vector_BI  breakpoint_space__hap_1;

    const uint number_of_haploid_breakpoints_0
	    = create_haploid_breakpoint_space
			(breakpoint_space__hap_0,
			full_haploid_sparse_state_vector_0,
			total_Breakpoint_complex);
	
    const uint number_of_haploid_breakpoints_1
	    = create_haploid_breakpoint_space
			(breakpoint_space__hap_1,
			full_haploid_sparse_state_vector_1,
			total_Breakpoint_complex);
			
			
    allocate_diploid_breakpoints_for_given_diploid_states(
					full_haploid_sparse_state_vector_0, full_haploid_sparse_state_vector_1,
					breakpoint_space__hap_0, breakpoint_space__hap_1);


			
    #pragma omp parallel	
    {//parallel
    
	type_map_diploid_neighbors_brkpts_to_traceback_data_storage thread_specific__diploid_breakpoints_to_data(diploid_breakpoints_to_data);    
	{
	    #pragma omp barrier  
	    // sync threads.
	}

	type_map_uint_to_BI breakpoints__haploid_0;
	initialize_map_with_keys_from_singlecontainer<uint, BOOST_Interval, type_set_uint>(
			    breakpoints__haploid_0,
			    extract_keys_of_map_and_return_as_set<uint,uint>(full_haploid_sparse_state_vector_0.sparse_map_UID_to_value), empty_BI);
	
	type_map_uint_to_BI breakpoints__haploid_1;
	initialize_map_with_keys_from_singlecontainer<uint, BOOST_Interval, type_set_uint>(
			    breakpoints__haploid_1,
			    extract_keys_of_map_and_return_as_set<uint,uint>(full_haploid_sparse_state_vector_1.sparse_map_UID_to_value), empty_BI);    
    
	#pragma omp for collapse(2) schedule(dynamic,50)		
	for (uint brkpt_ctr_0 = 0; brkpt_ctr_0 < number_of_haploid_breakpoints_0; ++brkpt_ctr_0)
	{
	    for (uint brkpt_ctr_1 = 0; brkpt_ctr_1 < number_of_haploid_breakpoints_1; ++brkpt_ctr_1)
	    {			    
		set_breakpoints_using_breakpoint_space(breakpoints__haploid_0, breakpoint_space__hap_0, brkpt_ctr_0);
		set_breakpoints_using_breakpoint_space(breakpoints__haploid_1, breakpoint_space__hap_1, brkpt_ctr_1);    
				
		
		real prod_PSF(P__E_diploid__cond_theta);

		for (uint l_PSF = 0; l_PSF < total_number_of_related_PSF; ++l_PSF )  
		{
		    prod_PSF *= list_of_related_PSF__and__thread_specific____new_PSF[l_PSF]
					->get_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints(
							    full_haploid_sparse_state_vector_0,
							    full_haploid_sparse_state_vector_1,
							    breakpoints__haploid_0,
							    breakpoints__haploid_1);
		}

		    
		type_BI__BI diploid_EE_brkpt(empty_BI__BI);                
		{
		    const type_map_uint_to_BI::const_iterator it_EE_brk = breakpoints__haploid_0.find(EE_UID);
		    if (it_EE_brk != breakpoints__haploid_0.end()) 
		    {
			diploid_EE_brkpt.first = it_EE_brk->second;
			prod_PSF *= P__s_haploid__cond__E_haploid_and_theta.first;
		    }
		}                
		{
		    const type_map_uint_to_BI::const_iterator it_EE_brk = breakpoints__haploid_1.find(EE_UID);
		    if (it_EE_brk != breakpoints__haploid_1.end())  
		    {
			diploid_EE_brkpt.second = it_EE_brk->second;
			prod_PSF *= P__s_haploid__cond__E_haploid_and_theta.second;
		    }
		}
		
		
		
		type_map_uint_to_BI__2 non_EE_haploid_breakpoints(breakpoints__haploid_0, breakpoints__haploid_1);
		non_EE_haploid_breakpoints.first.erase(EE_UID);
		non_EE_haploid_breakpoints.second.erase(EE_UID);                                
		    
		
		
		thread_specific__diploid_breakpoints_to_data[non_EE_haploid_breakpoints]
		    .map_eliminating_Event_diploid_outcomes_to_hybrid_factor
			[diploid_EE_state]
			    .map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts
				[diploid_EE_brkpt] = prod_PSF;
				
	    }//breakpoint_iterator  1        
	}//breakpoint_iterator  0		
	
	#pragma omp critical (absorb_diploid_neighbor_brkpts_per_thread_inside_Cartesian_product_of_PSFs)
	absorb_piece_of_diploid_neighbor_brkpts(thread_specific__diploid_breakpoints_to_data, diploid_breakpoints_to_data);
	
	mpfr_free_cache();
    }//parallel
    
	
    map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage[the_diploid_states] = diploid_breakpoints_to_data;

} // save_Cartesian_product_of_related_PSFs_wrt_eliminating_Event_and_include_probabilities_of_Breakpoint_and_Event___THREAD_SPECIFIC




























void Partial_sum_function::remove_eliminated_Event_from_storage
                                                    (const uint &eliminated_Event)
{
    
    haploid_argument_UIDs.erase(eliminated_Event);  
    eliminating_Event_UID = -1;
    
    for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::iterator it_dip_neighb_states 
				= map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.begin();
	    it_dip_neighb_states != map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end();
	    ++it_dip_neighb_states)
    {
	//remove eliminiating Event data:
	for (type_map_diploid_neighbors_brkpts_to_traceback_data_storage::iterator it_diploid_neighbor_brkpts = it_dip_neighb_states->second.begin();
		it_diploid_neighbor_brkpts != it_dip_neighb_states->second.end();
		++it_diploid_neighbor_brkpts)                                            
	{ 
	    it_diploid_neighbor_brkpts->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.clear();  
	}//brkpts
    }//states
                        
                        
    current_relevant_variational_positions_per_Event_in_PSF.erase(eliminated_Event);
            
    mpfr_free_cache();
       
}  // remove_eliminated_Event_from_storage













void Partial_sum_function::make_numerically_stable_by_multiplying_every_partial_sum_value_by_the_inverse_of_the_smallest_partial_sum_value()
{
    
    
    long int power_of_2_for_exponent;    
    
    //find minimal value.
    
    if (my_MPI_rank == 0)
    {
    
	real minimum_value_in_PSF(0.00L);

	for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::iterator it_diploid_states
					= map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.begin();
		    it_diploid_states != map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end();
		    ++it_diploid_states)
	{
	    for (type_map_diploid_neighbors_brkpts_to_traceback_data_storage::iterator it_diploid_neighbor_brkpts = it_diploid_states->second.begin();
		    it_diploid_neighbor_brkpts != it_diploid_states->second.end();
		    ++it_diploid_neighbor_brkpts) 
	    {
		if ( !mpfr::iszero(it_diploid_neighbor_brkpts->second.partial_sum_value_S) )                 
		{
		    if ( mpfr::iszero(minimum_value_in_PSF)
			    or   it_diploid_neighbor_brkpts->second.partial_sum_value_S < minimum_value_in_PSF  )
		    {  minimum_value_in_PSF = it_diploid_neighbor_brkpts->second.partial_sum_value_S;  }
		}
	    }
	} // states
         
	    
	if (mpfr::iszero(minimum_value_in_PSF)  or    minimum_value_in_PSF <= 0.00L)     
	{	
	    std::stringstream error_strm;	    
	    error_strm << "\n\nWARNING!    minimum_value_in_PSF  is  zero  or negative !!!\n\n"
			<< "\t\tminimum_value_in_PSF = " << minimum_value_in_PSF
			<< "\t\tsetting to   std::numeric_limits<unsigned long int>::max() ^ 2   =      ";	    
	    minimum_value_in_PSF =    (  std::numeric_limits<long int>::max()   );
	    minimum_value_in_PSF = mpfr::pow(minimum_value_in_PSF, 2);
	    minimum_value_in_PSF = mpfr::pow(minimum_value_in_PSF, -1);
	    
	    error_strm << minimum_value_in_PSF << "\n\n\n";	    
	    warning_message(error_strm, false);     
	}
         
            
            
            
	std::cerr <<  "\n\n\nmaking numerically stable:\n\t\tminimum_value_in_PSF  =    " << minimum_value_in_PSF;			
	{
	    const longdouble base_d = mpfr_get_ld_2exp( &power_of_2_for_exponent,  minimum_value_in_PSF.mpfr_ptr(), MPFR_RNDN );
	    	    	    
	    if (base_d == 0.00L  or  power_of_2_for_exponent == 0L  or  power_of_2_for_exponent == NAN  or  base_d == NAN)
	    {	    	    
		power_of_2_for_exponent = 100000L;  // 100,000
		if ( minimum_value_in_PSF >= 1.00L)
		    power_of_2_for_exponent *= -1;
	    }
	    else if ( minimum_value_in_PSF >= 1.00L)
	    {
		power_of_2_for_exponent = std::max<long int>(-power_of_2_for_exponent, -100000L);	    
	    }
	    else
	    {
		power_of_2_for_exponent = std::min<long int>(-power_of_2_for_exponent, 100000L);	    
	    }	    	    	    
	    
	    std::cerr << "\n\t\t\tbase_d = " << base_d  <<  "   and  power_of_2_for_exponent = " << power_of_2_for_exponent << "\n\n";        
	}
    
    }// 0
    
    
    boost::mpi::broadcast<long int>( *world_MPI_boost_ptr, power_of_2_for_exponent, 0 );
        

    


    for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::iterator it_diploid_states
				    = map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.begin();
                it_diploid_states != map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end();
                ++it_diploid_states)
    {
        //remove eliminiating Event data:
        for (type_map_diploid_neighbors_brkpts_to_traceback_data_storage::iterator it_diploid_neighbor_brkpts = it_diploid_states->second.begin();
                it_diploid_neighbor_brkpts != it_diploid_states->second.end();
                ++it_diploid_neighbor_brkpts)  
	{  it_diploid_neighbor_brkpts->second.partial_sum_value_S <<= power_of_2_for_exponent;  }  //(  <<=   means    fast mult by power of 2)                                                                    
    } // states
                                           
         
         
    mpfr_free_cache();
             
}   // make_numerically_stable_by_multiplying_every_partial_sum_value_by_the_inverse_of_the_smallest_such_value































void  Partial_sum_function::absorb_piece_of_diploid_neighbor_brkpts
				(const type_map_diploid_neighbors_brkpts_to_traceback_data_storage &diploid_neighbor_brkpts___source,
				 type_map_diploid_neighbors_brkpts_to_traceback_data_storage &diploid_neighbor_brkpts___destination)
{
    
    for (type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_diploid_neighbor_brkpts = diploid_neighbor_brkpts___source.begin();
	    it_diploid_neighbor_brkpts != diploid_neighbor_brkpts___source.end();
	    ++it_diploid_neighbor_brkpts)
    {	
	const type_map_diploid_neighbors_brkpts_to_traceback_data_storage::iterator  the_traceback_data_storage_for_saving
			= diploid_neighbor_brkpts___destination.insert(
				    std::pair<type_map_uint_to_BI__2, Partial_sum_function::traceback_data_storage>
					(it_diploid_neighbor_brkpts->first,  Partial_sum_function::traceback_data_storage())).first;
				    
	for (traceback_data_storage::type_map_eliminating_Event_diploid_outcomes_to_hybrid_factor::const_iterator it_EE_outcome
				= it_diploid_neighbor_brkpts->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.begin();
		it_EE_outcome != it_diploid_neighbor_brkpts->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor.end();
		++it_EE_outcome)
	{
	    const traceback_data_storage::type_map_eliminating_Event_diploid_outcomes_to_hybrid_factor::iterator   updating_hybrid_factor
			=  the_traceback_data_storage_for_saving->second.map_eliminating_Event_diploid_outcomes_to_hybrid_factor
				.insert(std::pair<type_uint__uint, traceback_data_storage::hybrid_factor>
					    (it_EE_outcome->first,  traceback_data_storage::hybrid_factor())).first;
						

	    type_map_BI__BI___to__real::iterator insert_it_for_saving =  updating_hybrid_factor->second.
			map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts.begin();
						    
	    for (type_map_BI__BI___to__real::const_iterator it_EE_brkpts = it_EE_outcome->second.
			map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts.begin();
		    it_EE_brkpts != it_EE_outcome->second.
			map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts.end();
		    ++it_EE_brkpts)
	    {
		insert_it_for_saving = updating_hybrid_factor->second.
					    map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts
						    .insert(insert_it_for_saving, *it_EE_brkpts); 
	    }
	}//EE_outcome
    }  //brkpts         

}//absorb_piece_of_diploid_neighbor_brkpts





void Partial_sum_function::absorb_piece_of_Partial_Sum_Function_into_one_Partial_Sum_Function
						(const Partial_sum_function &piece_of_new_Visitation_PSF)
{

            
    for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_diploid_neighbor_states
                                        =  piece_of_new_Visitation_PSF.map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.begin();
            it_diploid_neighbor_states != piece_of_new_Visitation_PSF.map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end();
            ++it_diploid_neighbor_states)
    {                
        //check if already exists
        const type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::iterator it_found
                                = map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.find(it_diploid_neighbor_states->first);                             
                                
        if ( it_found == map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.end() )  
        {
	    map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.insert(*it_diploid_neighbor_states);
        }
        else
        {                	    
	    absorb_piece_of_diploid_neighbor_brkpts(it_diploid_neighbor_states->second, it_found->second);	    
        } // else already exists, so update.                
    } // states       
    
    current_relevant_variational_positions_per_Event_in_PSF
		= combine_two__maps_to_set<uint,uint>(current_relevant_variational_positions_per_Event_in_PSF, 
						      piece_of_new_Visitation_PSF.current_relevant_variational_positions_per_Event_in_PSF);

}   //  absorb_pieces_of_Partial_Sum_Function_into_one_Partial_Sum_Function















void Partial_sum_function::clear_data()
{
    map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage.clear();    
    haploid_argument_UIDs.clear();       
    current_relevant_variational_positions_per_Event_in_PSF.clear();
    
}//clear_data
           
           
           
           
           
           
           
           
           
           
           
           
           























void  Partial_sum_function::copy_individual_likelihood_files_to_common_likelihood_directory_for_this_Event
                                            (const uint &given_likelihood_eliminating_Event_UID )  const
{        	

    char dest_filename[option_size];
    std::sprintf(dest_filename, "%s/Collect_Likelihoods_for_Event_%u/all_diploid_states_and_brkpts_from_rank_%d__thread_%d",
				scratch_job_output_dir.c_str(),
				given_likelihood_eliminating_Event_UID,
				my_MPI_rank, omp_get_thread_num() );
    
    std::ofstream out_fs( dest_filename,   std::ios::binary | std::ios::out  );  //  std::ios_base::binary     
			if ( !out_fs.good()  )
			{
			    std::stringstream error_strm;
			    print_line_of_markers("ERROR! ", &error_strm);
			    print_line_of_markers("(", &error_strm);
			    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
			    
			    error_strm << "\n\nERROR!   not good ofstream   in  \"copy_individual_likelihood_files_to_common_likelihood_directory_for_this_Event\" -  no IO,\n\t\t dest_filename =    "                              
					<<  dest_filename  << "\ngiven_likelihood_eliminating_Event_UID = " << given_likelihood_eliminating_Event_UID << "\n\n";
					
					
// 			    print_map_keys<type_map_uint_to_uint__2, type_map_diploid_neighbors_brkpts_to_traceback_data_storage, compare_pairs_of_maps<uint, uint, uint, uint> >( map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage, "map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage", &error_strm, true);
					
			    print_line_of_markers(")", &error_strm);
			    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
			}    
    
	
	
    boost::archive::binary_oarchive out_archive( out_fs );
    
    
    out_archive << map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage;
	    
    out_fs.close();	

} // copy_individual_likelihood_files_to_common_likelihood_directory_for_this_Event














void  Partial_sum_function::combine_individual_likelihood_files_into_common_likelihood_file_for_this_Event
                                           (const uint &given_likelihood_eliminating_Event_UID) 
{

    char collection_dirname[option_size];        
    std::sprintf(collection_dirname, "%s/Collect_Likelihoods_for_Event_%u",
				scratch_job_output_dir.c_str(),
				given_likelihood_eliminating_Event_UID   );                                
    const boost::filesystem::path  collection_dir(collection_dirname);
				
    
    char likelihoods_filename[option_size];        
    std::sprintf(likelihoods_filename, "%s/Likelihoods/Likelihoods_for_Event_%u",
				output_dir.c_str(), 
				given_likelihood_eliminating_Event_UID   );                   
    
	
	    
    std::ofstream out_fs( likelihoods_filename,   std::ios::binary | std::ios::out  );
			if ( !out_fs.good()  )
			{
			    std::stringstream error_strm;
			    print_line_of_markers("ERROR! ", &error_strm);
			    print_line_of_markers("(", &error_strm);
			    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
			    
			    error_strm << "\n\nERROR!   not good ofstream   in  \"combine_individual_likelihood_files_into_common_likelihood_file_for_this_Event\",\n\t\t likelihoods_filename =    ["                              
					<<  likelihoods_filename  << "]\n\n";                                                                                             
					
			    print_line_of_markers(")", &error_strm);
			    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
			}                                    
				
    boost::archive::binary_oarchive out_archive( out_fs );
    
    

    type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage  all_diploid_states_to_brkpts;
				
    const boost::filesystem::directory_iterator end_dir_it;                                                            
    for (boost::filesystem::directory_iterator it_dir( collection_dir  );
	    it_dir != end_dir_it;
	    ++it_dir)
    {                  	
	// open the archive
	std::ifstream in_fs(  (*it_dir).path().c_str(),   std::ios::binary | std::ios::in );  //  std::ios_base::binary    
				if ( !in_fs.good()  )
				{
				    std::stringstream error_strm;
				    print_line_of_markers("ERROR! ", &error_strm);
				    print_line_of_markers("(", &error_strm);
				    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
				    
				    error_strm << "\n\nERROR!   not good ifstream   in  \"combine_individual_likelihood_files_into_common_likelihood_file_for_this_Event\",\n\t\t (*it_dir).path().c_str() =    "  
				    <<  (*it_dir).path().c_str()  << "\n\n";
				    
				    print_line_of_markers(")", &error_strm);
				    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
				}	    
	
	
	boost::archive::binary_iarchive in_archive( in_fs );
			
	{
	    type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage in_map_diploid_states_to_brkpts;	    
	    in_archive >> in_map_diploid_states_to_brkpts;
	    	    
	    deep_insert(all_diploid_states_to_brkpts, in_map_diploid_states_to_brkpts);
	}
		
	in_fs.close();	
	    
    } // it_dir
    	
		
		
		
		
 
    out_archive <<  all_diploid_states_to_brkpts;
	    	    
    out_fs.close(); 
            
    
    boost::filesystem::remove_all(  collection_dirname  );
       

} // combine_individual_likelihood_files_into_common_likelihood_file_for_this_Event








































void Partial_sum_function::upload_all_likelihood_factors_from_file_for_a_given_Event
				    (const std::string &other_job_run_data_dir,
				    const uint &some_eliminating_Event_UID)
{ 
    
    char likelihoods_filename[option_size];        
    std::sprintf(likelihoods_filename, "%s/Likelihoods/Likelihoods_for_Event_%u",
                                other_job_run_data_dir.c_str(), 
                                some_eliminating_Event_UID   );                   
    

    // open the archive
    std::ifstream in_fs(  likelihoods_filename,   std::ios::binary | std::ios::in );  //  std::ios_base::binary    
			    if ( !in_fs.good()  )
			    {
				//This could mean there is an error, but it might also mean that there is no likelihood for this Event alone - i.e. all of its coordinates were considered in previous steps of the variable eleiminiation schedule.
				if (my_MPI_rank == 0)
				{
				    std::stringstream warning_strm;
				    print_line_of_markers("ERROR! ", &warning_strm);
				    print_line_of_markers("(", &warning_strm);
				    warning_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
				    
				    warning_strm << "\n\nWARNING!   not good ifstream   in  \"upload_all_likelihood_factors_from_file_for_a_given_Event\",\n\t\t likelihoods_filename =    "  
				    <<  likelihoods_filename  << "\n\n";
				    
				    print_line_of_markers(")", &warning_strm);
				    std::fprintf(stderr, "\n\n%s\n\n", warning_strm.str().c_str() ); 	
				}
				
				return;			   
			    }

    boost::archive::binary_iarchive in_archive( in_fs );
    
    

    type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage  in_map_states_to_brkpts;    
    in_archive  >>  in_map_states_to_brkpts;
    
    map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage = in_map_states_to_brkpts;

    
    in_fs.close();         
                  

} // upload_likelihood_factor_from_file_for_a_given_Event































void Partial_sum_function::deep_insert
				(type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage &receiver,
				 const type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage &donor)
{
    for (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_dip_states = donor.begin();
	    it_dip_states != donor.end();
	    ++it_dip_states)
    {
	const bool did_something_new = receiver.insert(*it_dip_states).second;
	if (!did_something_new)
	{
	    for (type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator it_dip_brkpts = it_dip_states->second.begin();
		    it_dip_brkpts != it_dip_states->second.end();
		    ++it_dip_brkpts)
	    {
		receiver.at(it_dip_states->first).insert(*it_dip_brkpts);	    
	    }//brkpts	
	}//must go deeper    
    }//states

}//deep_insert






























void Partial_sum_function::record_additional_relevant_variational_positions_per_Event
					(const type_map_uint_to_set_uint &map_more_relevant_varpos_per_Event)
{
    current_relevant_variational_positions_per_Event_in_PSF = combine_two__maps_to_set<uint,uint>(map_more_relevant_varpos_per_Event, current_relevant_variational_positions_per_Event_in_PSF);    

}//record_additional_relevant_variational_positions_per_Event





