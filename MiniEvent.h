#ifndef MINIEVENT_H_INCLUDED
#define MINIEVENT_H_INCLUDED





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



#include <api/BamAlignment.h>



#include <general_typedefs.h>
#include <Sparse_map.h>
#include <Readgroup_statistics.h>



#include "globals.h"

#include <CUPA.h>





class Event;
class Paired_end_read;


class MiniEvent
{
    public:
	//constructor:
	MiniEvent () : full_Event(NULL), parent_PER(NULL)
	{ }
	
	MiniEvent ( const Event *const &in_full_Event,
		    const Paired_end_read *const &in_parent_PER) : full_Event(in_full_Event),
								    parent_PER(in_parent_PER)
	{ }
	
// 	//copy constructor
// 	MiniEvent(const MiniEvent &another_ME)
// 			:  full_Event(another_ME.full_Event),
// 			    parent_PER(another_ME.parent_PER),
// 			    affected_MiniRegions(another_ME.affected_MiniRegions),
// 			    map_MiniRegions_ON_LCRs__to__intersecting_LCR(another_ME.map_MiniRegions_ON_LCRs__to__intersecting_LCR),
// 			    MiniRegions_in_Inbetween_of_this_MiniEvent(another_ME.MiniRegions_in_Inbetween_of_this_MiniEvent),
// 			    partners_for_MiniRegions_on_LCRs(another_ME.partners_for_MiniRegions_on_LCRs),
// 			    MiniRegions_in_insertion_of_an_LCR(another_ME.MiniRegions_in_insertion_of_an_LCR),
// 			    map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
// 						(another_ME.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA),
// 			    map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion
// 					(another_ME.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion)
// 	{ }
	
	//destructor:
	~MiniEvent()
	{ }
				    
	
	//variables:          
	const Event *full_Event;
	const Paired_end_read *parent_PER;
	
	type_set_uint  affected_MiniRegions;
			//includes MiniRegions that intersect either LCR of this MiniEvent, AND includes MiniRegions that intersect the region "Inbetween" the two LCRs of the MiniEvent (if is intra-chromosomal Event, of course).
	
	
	type_map_uint_to_uint map_MiniRegions_ON_LCRs__to__intersecting_LCR;
	//  ===>  <====   ====>   <====   ====>   <====   ====>   <====   ====>   <====   ====>   <====   ====>   <====   ====>   <====   ====>   <====  //
	type_set_uint MiniRegions_in_Inbetween_of_this_MiniEvent;
	
	
	
	//Note:  "(keys of)map_MiniRegions_ON_LCRs__to__intersecting_LCR"  DISJOINT UNION  "MiniRegions_in_Inbetween_of_this_MiniEvent"  =  affected_MiniRegions.
	
	
	type_map_uint_to_set_uint partners_for_MiniRegions_on_LCRs;   // for MiniRegions on LCRs only, of course.
	type_set_uint MiniRegions_in_insertion_of_an_LCR;  //that are in an "insertion" of one of the LCRs wrt the other LCR.
	    // Note:  this is different from "MiniRegions_in_Inbetween_of_this_MiniEvent".  Both are "partnerless", but are pertinent to different regions of this MiniEvent.
	    // Note:  this is a subset of   "affected_MiniRegions"  and   "map_MiniRegions_ON_LCRs__to__intersecting_LCR".  It is disjoint from "MiniRegions_in_Inbetween_of_this_MiniEvent".
						
	
	type_map_BI_to__real__real map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA;                                                
	//New strategy:
	//	We always want to test Gene Conversion because it helps avoid making bogus NAHR calls.
	//	Thus, our set of possible "breakpolints" should be singleton coordinates.  Then, whenever we loop through breakpoints, we always do an ordered-double-for-loop to get all possible APPROPRIATE pairs.
	//	But then for paired-end reads, we need to store probabilities related to PAIRS of breakpoints, where "degenerate" or "singleton pairs, i.e. closed interval whose bounds are equal, are recognized as NAHR breakpoints (just a single breakpoint), where as non-degenerate intervals are recogznied as Gene Conversion breakpoints.              
	
	type_map_uint_to_set_uint map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion;
												// (brkpts wrt this MiniEvent, of course).
			//Not that by construction (in "PER.align_to_every_possible_..."), MiniRegions in the INBEWTEEN are NOT contained in this map.
	
	
	type_set_BI  relevant_breakpoints_for_which_this_PER_actually_shows_a_hybrid_pattern;//only for the heuristic. 
	
	
	//functions:
	void print_this_MiniEvent(std::stringstream *const &some_ss = NULL)  const;
	
	type_uint__real  calculate__haploid_sum___P_mPER__cond_this_Event_occurred_with_these_relevant_profile_breakpoints____ONLY_SUM_OVER_MINIREGIONS_ON_THE_LCRS
						    ( const type_haploid_outcome &miniEvent_outcome,
						    const type_set_uint &set_MiniRegions_on_LCRs_of_this_MiniEvent,
						    const BOOST_Interval &miniEvent_breakpoint)  const;                                           
	    // depends on outcome (.e.g. dup vs del)!!!!!!
	    //this function should be called for all occurring events.    Probability for mappings from non-occurring events should only be called after!! (so that we have the list of "already_considered_MiniIDs".
	    
		


	

    private:             
	real get____P__m_PER_____where_m_maps_to_the_hybrid_LCR_created_from_given_MiniRegion
					(const bool &is_NAHR,
					    const bool &hybrid_type,      // false  <==>   "AB" (NAHR)  /  "ABA" (GeneConv)
					    const uint &affected_MiniID,  
					    const BOOST_Interval &miniEvent_profile_breakpoint,
					    bool &already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__LOWER,
					    bool &already_added_probability_of_mapping_for_one_of_the_MiniRegions_that_contains_the_breakpoint__UPPER)  const;
						
				
				
    public:                             
	
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{                
	    
	    ar  &  affected_MiniRegions;
	    
	    ar  &  map_MiniRegions_ON_LCRs__to__intersecting_LCR;
	    ar  &  MiniRegions_in_Inbetween_of_this_MiniEvent;
	    
	    ar  &  partners_for_MiniRegions_on_LCRs;
	    ar  &  MiniRegions_in_insertion_of_an_LCR;
	    
	    ar  &  map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA;
	    ar  &  map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion;
	};
	

					    
};//MiniEvent




      
typedef std::map<uint, MiniEvent>   type_map_uint_to_MiniEvent;
typedef std::pair<uint, MiniEvent>  type_uint__MiniEvent;









#endif