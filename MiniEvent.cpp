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


#include <mpreal.h>
#include <mpfr.h>



#include <general_typedefs.h>
#include <SD_Entry_general.h>
#include <templates.h>
#include <translations_through_alignment.h>

#include "globals.h"
#include "Paired_end_read.h"
#include "other_functions.h"
#include "Event.h"
#include "MiniRegion.h"








void MiniEvent::print_this_MiniEvent(std::stringstream *const &some_ss)  const
{
    
    std::stringstream output_str;
    
    const bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str; 
    
    print_line_of_markers("[", ss_ptr);    
    
    (*ss_ptr) <<  "\n\n\tprinting MiniEvent:\n\tfull_Event->UID =   " <<  full_Event->UID  <<  "\n";
    print_set<uint>(affected_MiniRegions, "\t\taffected_MiniRegions", ss_ptr);
    print_map_keys_and_values<uint, uint>(map_MiniRegions_ON_LCRs__to__intersecting_LCR, "\t\tmap_MiniRegions_ON_LCRs__to__intersecting_LCR", ss_ptr);
    
//     print_map_keys_and_pair_values<BOOST_Interval, real, real>(map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA,
//                                                      "map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA",
//                                                      ss_ptr);
                                                     
            
    (*ss_ptr) << "\n\t\tmap_affected_MiniID_to_relevant_breakpoints_for_that_MiniRegion:\n";
    for (type_map_uint_to_set_uint::const_iterator it = map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion.begin();
        it != map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion.end();
        ++it)            
    {
        (*ss_ptr) << "\t\t\t" << it->first << ":  ";
        print_set<uint>(it->second, NULL, ss_ptr);                                
    }
    
    
    print_set<uint>(MiniRegions_in_Inbetween_of_this_MiniEvent, "MiniRegions_in_Inbetween_of_this_MiniEvent", ss_ptr);
    print_set<uint>(MiniRegions_in_insertion_of_an_LCR, "MiniRegions_in_insertion_of_an_LCR", ss_ptr);
    
    
    print_map_keys_and_pair_values<BOOST_Interval, real, real, compare_BI>(
	    map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA, "map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA", ss_ptr, true);
    
    
    
    print_line_of_markers("]", ss_ptr);  
    
    if (!provided_with_valid_sstream)
        std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );
          
}   // end of   print_this_MiniEvent






























type_uint__real MiniEvent::
            calculate__haploid_sum___P_mPER__cond_this_Event_occurred_with_these_relevant_profile_breakpoints____ONLY_SUM_OVER_MINIREGIONS_ON_THE_LCRS
                                                          ( const type_haploid_outcome &miniEvent_outcome,
                                                            const type_set_uint &set_MiniRegions_on_LCRs_of_this_MiniEvent,
							    const BOOST_Interval &miniEvent_breakpoint)  const
{
//Given a SPARSE state vector and breakpoints, we will call this function for each of the (occurring) Events in the state vector.  After calling this function for every occurring event, we will be left with a set of MiniIDs ("already_considered_MiniIDs") whose mappings have already contributed to the sum.  Thus, after calling this function for every Event on the state vector, we know which MiniRegions remain for their non-hybrid probabilities to be added to the running sum.

//CAREFUL:  WE NEED TO KNOW WHICH OF THE POTENTIALLY-AFFECTED MINIREGIONS HAVE ALREADY HAD THEIR NON-HYBRID PROBABILITIES CONTRIBUTE TO THE SUM!!! 
  
//Recall that we only call this function for MiniEvents  that OCCURRED.  
    
    real sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent(0.00L);
    uint number_of_contributing_regions = 0;            

    switch (miniEvent_outcome)
    {
	case hap_outcome__Del:  // deletion    -->    LCR_0 and LCR_1 hybridize to become  LCR_01. 
	{
	    bool already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint = false;
	    bool dummy_second_monitor = false;
	    
	    for (type_set_uint::const_iterator it_MR_on_LCR = set_MiniRegions_on_LCRs_of_this_MiniEvent.begin();
		    it_MR_on_LCR != set_MiniRegions_on_LCRs_of_this_MiniEvent.end();
		    ++it_MR_on_LCR)
	    {
		const real returned_value( 
				get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion
						    (true,   //is NAHR
						     false,    //  01
						    *it_MR_on_LCR,
						    miniEvent_breakpoint,
						     already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint,
						     dummy_second_monitor)    );  
		if (returned_value > 0.00L)                                    
		{
		    sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent +=  returned_value;
		    ++number_of_contributing_regions;
		}
	    }//MR     
	    
	    break;
	}//del
                
	case hap_outcome__Dup:  //duplication
	{
	    //for a duplication - both original LCRs (LCR_0 and LCR_1) still exist, but now there is an additional LCR_10.
	    //For clarity, first we add up all the non-hybrid mappings (since the original LCRs still exist).
	    //then we add up the probabilities LCR_10.  this follows the deletion case above.
	    
	    bool already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint = false;
	    bool dummy_second_monitor = false;
	    
	    for (type_set_uint::const_iterator it_MR_on_LCR = set_MiniRegions_on_LCRs_of_this_MiniEvent.begin();
		    it_MR_on_LCR != set_MiniRegions_on_LCRs_of_this_MiniEvent.end();
		    ++it_MR_on_LCR)      
	    {                      
		
		//normal copies of LCRs A and B:
		sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent
						+= parent_PER->map_affected_MiniID__to__P_nohybrid.at(*it_MR_on_LCR); 
		++number_of_contributing_regions;
		
		//hybrid BA:
		const real returned_value(
			    get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion
						(true,   //is NAHR
						 true,    //  10
						*it_MR_on_LCR,
						miniEvent_breakpoint,
						 already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint,
						 dummy_second_monitor)    );
						
						
		if (returned_value > 0.00L)                                    
		{
		    sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent +=  returned_value;
		    ++number_of_contributing_regions;
		}  
	    }//MR
	    break;
	}//dup
	
	case hap_outcome__GeneConv_ABA:
	case hap_outcome__GeneConv_BAB:
	{
	    //Similar idea to duplication:
	    bool already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__LOWER = false;
	    bool already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__UPPER = false;
	    
	    for (type_set_uint::const_iterator it_MR_on_LCR = set_MiniRegions_on_LCRs_of_this_MiniEvent.begin();
		    it_MR_on_LCR != set_MiniRegions_on_LCRs_of_this_MiniEvent.end();
		    ++it_MR_on_LCR)      
	    {                      
		//add normal copy of donor LCR.  if "ABA", then donor = LCR "B", so add normal for LCR "B".  if "BAB"  then add normal for "A".
		const uint intersecting_LCR_qs = map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(*it_MR_on_LCR);
		
		if (    (miniEvent_outcome == hap_outcome__GeneConv_ABA  and  intersecting_LCR_qs == 1)
		     or (miniEvent_outcome == hap_outcome__GeneConv_BAB  and  intersecting_LCR_qs == 0)  )
		{
		    sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent
						+= parent_PER->map_affected_MiniID__to__P_nohybrid.at(*it_MR_on_LCR);     
		    ++number_of_contributing_regions;
		}
		
		//hybrid:					    		
		const real returned_value(
			    get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion
						(false,   //is NAHR
						 miniEvent_outcome == hap_outcome__GeneConv_BAB,  //false "ABA", true "BAB"
						*it_MR_on_LCR,
						miniEvent_breakpoint,
						 already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__LOWER,
						 already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__UPPER)    );
						
		if (returned_value > 0.00L)                                    
		{
		    sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent +=  returned_value;
		    ++number_of_contributing_regions;
		}
	    }//MR
	    
	    break;
	}//GeneConv
		          
		          
	case hap_outcome__Inv:
	{
	    //in this case, we will have TWO hybrids - BOTH hybrid_AB and hybrid_BA !!!!        
	    bool already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_a_breakpoint_01 = false;
	    bool already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_a_breakpoint_10 = false;	
	    bool dummy_second_monitor = false;
	    
	    for (type_set_uint::const_iterator it_MR_on_LCR = set_MiniRegions_on_LCRs_of_this_MiniEvent.begin();
		    it_MR_on_LCR != set_MiniRegions_on_LCRs_of_this_MiniEvent.end();
		    ++it_MR_on_LCR)    
	    {                     
		const real returned_value_0(
			    get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion
							(true,   //is NAHR
							false,    //  01
							*it_MR_on_LCR,
							miniEvent_breakpoint,
							already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_a_breakpoint_01,
							dummy_second_monitor)      );    
								
		const real returned_value_1(          
			    get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion
							(true,   //is NAHR
							true,    //  10    
							*it_MR_on_LCR,
							miniEvent_breakpoint,
							already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_a_breakpoint_10,
							dummy_second_monitor)        );  
							
		if (returned_value_0 > 0.00L)                                        
		{
		    sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent += returned_value_0;
		    ++number_of_contributing_regions;           
		}
		
		if (returned_value_1 > 0.00L)                                        
		{
		    sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent += returned_value_1;
		    ++number_of_contributing_regions;           
		}                                                                                                                                
	    }   // end     for-loop        
		
	    break;
	}//inv 
    
    
	default:
	{
	    std::stringstream error_strm;
	    error_strm << "\n\nin \"calculate__haploid_sum___P_mPER__cond_this_Event_occurred_with_these_relevant_profile_breakpoints____ONLY_SUM_OVER_MINIREGIONS_ON_THE_LCRS\".\n\n unrecognized recomb type = " << full_Event->recomb_type << "  !!!\n\n\n";
				    
	    error_strm << "\n\n\nMiniEvent generating this error:\n\n";
	    print_this_MiniEvent(&error_strm);
	    
	    error_strm << "\n\n\n\nPER generating this error:\n\n";
	    parent_PER->print_this_PER(&error_strm);  
	    
	    error_message(error_strm, true);
	    break;                           
	}
    
    
    }; //switch - outcome
        
    return type_uint__real(number_of_contributing_regions, sum_P_mappings_for_those_mappings_to_MiniRegions_affected_by_this_MiniEvent);

}  // end of     calculate__haploid_sum___P_mPER__cond_this_Event_occurred_with_these_proposed_profile_breakpoints




















real MiniEvent::get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion
                                                (const bool &is_NAHR,
						 const bool &hybrid_type,      // false  <==>   "AB" (NAHR)  /  "ABA" (GeneConv)
                                                 const uint &affected_MiniID,  
                                                 const BOOST_Interval &miniEvent_profile_breakpoint,
						 bool &already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__LOWER,
						 bool &already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__UPPER)  const						 
{
//hybrid_type = false   <==>   "AB" (NAHR)  /  "ABA" (GeneConv)
//	      = true    <==>   "BA" (NAHR)  /  "BAB" (Gene Conv) 
    
//should only be called for MiniRegions on LCR of MiniEvent!!!    
    
//this is only for MiniRegions that intersect the profile.  (i.e., if MiniEvent is a dupdel and the MiniRegion lies in the DAR, then returns 0.00, i.e. that MiniRegion is ignored here).

//"already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint"     this is necessary because there will be either exactly 2 regions or exactly 0 regions overlapping a breakpoint.  If 2 regions intersect the breakpoint, then in the case of the hybrid, we only want to consider one of the MiniRegions that overlap the breakpoint.
                     
    
	    
    //right off the bat, check if the desired breakpoint exists....
    const type_map_BI_to__real__real::const_iterator it_found = map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.find(miniEvent_profile_breakpoint);    
    if (it_found != map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.end())//breakpoint is in MiniRegion
    {
	if (already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__LOWER
	    and already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__UPPER)
	{  return 0;  }
	else
	{                
	    already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__LOWER = true;                                
	    already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__UPPER = true;
	    //Note that   "map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA"  was constructed so that AB ALWAYS means from LCR_0 to LCR_1.  Thus we don't actually have to check
	    //obviously      0 --> 01 --> AB,     and     1 --> 10 --> BA.     
	    
	    //hybrid_type = false   <==>   "AB" (NAHR)  /  "ABA" (GeneConv)
	    //	      = true    <==>   "BA" (NAHR)  /  "BAB" (Gene Conv)	    
	    if (!hybrid_type)
		return it_found->second.first;   // "AB"  /  "ABA"
	    else
		return it_found->second.second;    // "BA" / "BAB"
	}   		    
    }//found            
    else
    {//not found    
	
	const type_map_BI_to__real__real::const_iterator it_found_lower 
		    = map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.find(BOOST_make_point(miniEvent_profile_breakpoint.lower()));

	const type_map_BI_to__real__real::const_iterator it_found_upper 
		    = map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.find(BOOST_make_point(miniEvent_profile_breakpoint.upper()));
		
	//check if is Geneconversion and we see one of the sides.
	if (!is_NAHR  and  (it_found_lower != map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.end()
			    or  it_found_upper != map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.end()))
	{//GeneConv and we see one of the sides
	    if (it_found_lower != map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.end())
	    {
		if (already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__LOWER)
		{  return 0;  }
		else
		{
		    already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__LOWER = true;
		    if (!hybrid_type) //"ABA", we see "AB"
			return  it_found_lower->second.first;
		    else  //"BAB", we see "BA"
			return it_found_lower->second.second;
		}	    
	    }
	    else //found upper
	    {
		if (already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__UPPER)
		{  return 0;  }
		else
		{
		    already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__UPPER = true;
		    if (!hybrid_type) //"ABA", we see "BA"
			return  it_found_upper->second.second;
		    else  //"BAB", we see "AB"
			return it_found_upper->second.first;
		}	    
	    }	    	
	}//GeneConv and we see one of the sides
	else
	{//completely not found
	    //Thus, we do not observe the breakpoint, and it is not GeneConversion and we only see one of the sides.
	    //So the breakpoint must be "away" from the MiniRegion.
	    //Since the breakpoint is away form this MiniRegion, then we might as well just look at the "nohybrid" alignment for this MR. (since with or without the breakpoint won't have any impact on the alignment of this read.)
	
// 	    assert(map_MiniRegions_ON_LCRs__to__intersecting_LCR.count(affected_MiniID) > 0);
	    const uint intersecting_LCR_of_MR__qs = map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(affected_MiniID);
	
	    BOOST_Interval MR_nohybrid_algmt_on_Event_profile(empty_BI);
	    {//set MR_nohybrid_algmt_on_Event_profile
		const type_map_uint_to_BI::const_iterator it_absolute_nohybrid_algmt = parent_PER->map_MiniID__to__nohybrid_almgt_interval_absolute.find(affected_MiniID);
			    if (it_absolute_nohybrid_algmt == parent_PER->map_MiniID__to__nohybrid_almgt_interval_absolute.end())
			    {
				std::stringstream error_strm;
				error_strm << "\nERROR    affected_MiniID  = " << affected_MiniID 
					    << "\nit_absolute_nohybrid_algmt == parent_PER->map_MiniID__to__nohybrid_almgt_interval_absolute.end()  !!!\n\n";
					print_map_keys_and_values<uint, BOOST_Interval>(parent_PER->map_MiniID__to__nohybrid_almgt_interval_absolute,
											"parent_PER->map_MiniID__to__nohybrid_almgt_interval_absolute",
											&error_strm, true);
				print_this_MiniEvent(&error_strm);
				parent_PER->print_this_PER(&error_strm);
				error_message(error_strm, true);
			    }
		
		if (BOOST_subset(it_absolute_nohybrid_algmt->second, full_Event->LCRs[intersecting_LCR_of_MR__qs]))
		{		    
// 		    MR_nohybrid_algmt_on_Event_profile = get_compressed_map_value_endpoints_of_intersection_of_compressed_map_keys_and_interval(
// 									full_Event->compressed_map_LCR_to_profile__for_each_LCR.at(intersecting_LCR_of_MR__qs),
// 									it_absolute_nohybrid_algmt->second);
		    
		    MR_nohybrid_algmt_on_Event_profile = get_endpoints_of_values_of_compressed_map(
							    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
								    full_Event->compressed_map_LCR_to_profile__for_each_LCR.at(intersecting_LCR_of_MR__qs),
								    it_absolute_nohybrid_algmt->second));
		}
		else
		{
		    const bool MR_orientation_agrees_with_profile = (full_Event->recomb_type != recomb_class__Inv  or  intersecting_LCR_of_MR__qs == 0);
				    
		    uint profile_val_from_upper; 
		    uint profile_val_from_lower;
	    
// 		    if (it_absolute_nohybrid_algmt->second.upper() > full_Event->LCRs[intersecting_LCR_of_MR__qs].upper())
		    {
			if (it_absolute_nohybrid_algmt->second.upper() > full_Event->LCRs[intersecting_LCR_of_MR__qs].upper())
			    profile_val_from_upper = MR_orientation_agrees_with_profile  ? full_Event->profile_length  :    0;
			else if (it_absolute_nohybrid_algmt->second.upper() < full_Event->LCRs[intersecting_LCR_of_MR__qs].lower())
			    profile_val_from_upper = MR_orientation_agrees_with_profile  ?  0 : full_Event->profile_length;
			else
			{
// 			    profile_val_from_upper = get_compressed_map_value_endpoints_of_intersection_of_compressed_map_keys_and_interval(	
// 								full_Event->compressed_map_LCR_to_profile__for_each_LCR.at(intersecting_LCR_of_MR__qs),
// 								BOOST_make_point(it_absolute_nohybrid_algmt->second.upper())).upper();
			    
			    profile_val_from_upper = 
				    get_endpoints_of_values_of_compressed_map(
					    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
						    full_Event->compressed_map_LCR_to_profile__for_each_LCR.at(intersecting_LCR_of_MR__qs),
						    BOOST_make_point(it_absolute_nohybrid_algmt->second.upper()))).upper();
			} 			    
			
			if (it_absolute_nohybrid_algmt->second.lower() > full_Event->LCRs[intersecting_LCR_of_MR__qs].upper())
			    profile_val_from_lower = MR_orientation_agrees_with_profile  ? full_Event->profile_length  :    0;
			else if (it_absolute_nohybrid_algmt->second.lower() < full_Event->LCRs[intersecting_LCR_of_MR__qs].lower() )
			    profile_val_from_lower = MR_orientation_agrees_with_profile  ?   0  :  full_Event->profile_length;
			else
			{
// 			    profile_val_from_lower = get_compressed_map_value_endpoints_of_intersection_of_compressed_map_keys_and_interval(
// 								full_Event->compressed_map_LCR_to_profile__for_each_LCR.at(intersecting_LCR_of_MR__qs),
// 								BOOST_make_point(it_absolute_nohybrid_algmt->second.lower())).lower();
			    profile_val_from_lower = 
				    get_endpoints_of_values_of_compressed_map(
					    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
						    full_Event->compressed_map_LCR_to_profile__for_each_LCR.at(intersecting_LCR_of_MR__qs),
						    BOOST_make_point(it_absolute_nohybrid_algmt->second.lower()))).lower();
			}   
		    }
		    
		    MR_nohybrid_algmt_on_Event_profile.set( std::min<uint>(profile_val_from_lower, profile_val_from_upper),
							    std::max<uint>(profile_val_from_lower, profile_val_from_upper));
		}//not subset
		
			    if (BOOST_empty(MR_nohybrid_algmt_on_Event_profile))
			    {//error-check
				std::stringstream  error_strm;
				error_strm << "\n\nERROR:   MR_nohybrid_algmt_on_Event_profile = " << MR_nohybrid_algmt_on_Event_profile 
					    << "\n\tit_absolute_nohybrid_algmt = " << it_absolute_nohybrid_algmt->second
					    << "\n\tBOOST_subset(it_absolute_nohybrid_algmt->second, full_Event->LCRs[intersecting_LCR_of_MR__qs]) = " 
					    << BOOST_subset(it_absolute_nohybrid_algmt->second, full_Event->LCRs[intersecting_LCR_of_MR__qs])
					    << "\n\tintersecting_LCR_of_MR__qs = " << intersecting_LCR_of_MR__qs 
					    << "\n\n";

				error_message(error_strm, true);
			    }//error-check
	    }//set MR_nohybrid_algmt_on_Event_profile
	    
		    
	    
	    
	    
	    if (is_NAHR)
	    {//NAHR	
		if (MR_nohybrid_algmt_on_Event_profile.upper() < miniEvent_profile_breakpoint.lower())
		{
		    if (!hybrid_type) //AB		
			return   (intersecting_LCR_of_MR__qs == 0)   ?
				    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)   :   0;		
		    else //BA		
			return   (intersecting_LCR_of_MR__qs == 1)   ?
				    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)   :   0;					
		}//PER algmt   <  breakpoint
		else if (MR_nohybrid_algmt_on_Event_profile.lower() > miniEvent_profile_breakpoint.upper())
		{
		    if (!hybrid_type) //AB		
			return   (intersecting_LCR_of_MR__qs == 1)   ?
				    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)   :   0;
		    else //BA		
			return   (intersecting_LCR_of_MR__qs == 0)   ?
				    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)   :   0;					    
		}
		else 
		{
		    //can happen for the "beyond last varpos" scenario, for example.
// 		    std::stringstream error_strm;
// 		    error_strm << "error in  \"get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion\"  for PER name = [" 
// 				<< parent_PER->name << "].  NAHR,  miniEvent_profile_breakpoint = " << miniEvent_profile_breakpoint 
// 				<< ", not supposed to make it this far!!!  overlaps breakpoint, but somehow breakpoint was not recorded.\nMR_nohybrid_on_Event_profile  = " << MR_nohybrid_algmt_on_Event_profile << "\n\n";
// 		    error_message(error_strm, false);		    
		    return   parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID);	    
		}		
	    }//NAHR
	    else
	    {//GeneConv
		if (MR_nohybrid_algmt_on_Event_profile.upper() < miniEvent_profile_breakpoint.lower()
		    or   MR_nohybrid_algmt_on_Event_profile.lower() > miniEvent_profile_breakpoint.upper())
		{
		    if (!hybrid_type) // ABA		
			return   (intersecting_LCR_of_MR__qs == 0)   ?
				    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)   :   0;
		    else  //BAB		
			return   (intersecting_LCR_of_MR__qs == 1)   ?
				    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)   :   0;
						
		} // PER algmt < brkpt  or  > brkpt
		else if (BOOST_subset(MR_nohybrid_algmt_on_Event_profile, miniEvent_profile_breakpoint))
		{
		    if (!hybrid_type) // ABA		
			return   (intersecting_LCR_of_MR__qs == 1)   ?
				    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)   :   0;
		    else  //BAB		
			return   (intersecting_LCR_of_MR__qs == 0)   ?
				    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)   :   0; 
		}//subset
		else //if (BOOST_subset(miniEvent_profile_breakpoint, MR_nohybrid_algmt_on_Event_profile))
		{	
		    const uint intersect_size = BOOST_width_inclusive(BOOST_intersect(miniEvent_profile_breakpoint, MR_nohybrid_algmt_on_Event_profile));
		    
		    if ((double)BOOST_width_inclusive(MR_nohybrid_algmt_on_Event_profile)*0.05 > intersect_size) //barely affects MR
		    {
			if (!hybrid_type) // ABA
			    return   (intersecting_LCR_of_MR__qs == 0)   ?
					    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)  :  0;
			else //BAB
			    return  (intersecting_LCR_of_MR__qs == 1)   ?
					    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)  :  0;
		    }
		    else  //substantially affects MR
		    {
			if (!hybrid_type) // ABA
			    return   (intersecting_LCR_of_MR__qs == 1)   ?
					    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)  :  0;
			else //BAB
			    return  (intersecting_LCR_of_MR__qs == 0)   ?
					    parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID)  :  0;
		    }
// 		    std::stringstream error_strm;
// 		    error_strm << "error in  \"get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion\"  for PER name = [" 
// 				<< parent_PER->name << "].  GeneConv,  miniEvent_profile_breakpoint = " << miniEvent_profile_breakpoint 
// 				<< ", not supposed to make it this far!!! \nMR_nohybrid_on_Event_profile  = " << MR_nohybrid_algmt_on_Event_profile  <<  "\n\n";
// 				
// 		    print_set<BOOST_Interval, compare_BI>(	
// 				extract_keys_of_map_and_return_as_set<BOOST_Interval, type_real__real, compare_BI>(
// 					    map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA),
// 				"map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA-keys", &error_strm, true);
// 				
// 		    error_message(error_strm, false);
// 		    
// 		    return parent_PER->map_affected_MiniID__to__P_nohybrid.at(affected_MiniID);
		}			    
	    }//GeneConv			
	}//completely not found
    }//not found
    
    
   
}  // end of     get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion

