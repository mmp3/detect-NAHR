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
#include <boost/filesystem.hpp>


#include <gmpfrxx.h>
#include <boost/math/bindings/mpfr.hpp>


#include <mpreal.h>
#include <mpfr.h>


#include <api/BamAlignment.h>

#include <general_typedefs.h>



#include "Event.h"
#include "other_functions.h"
#include "globals.h"
#include "MiniRegion.h"
#include "Paired_end_read.h"
#include "MiniEvent.h"


#include <CUPA.h>

#include "PER_bundle_wrapper.h"











void PER_bundle_wrapper__space::PER_bundle_wrapper::print_this_PER_bundle_wrapper
					(std::stringstream *const &some_ss)  const
{
    
    std::stringstream output_str;
    
    const bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str; 

    (*ss_ptr) << "\n\nprinting PER_bundle_wrapper...\n\n\n";
    
    (*ss_ptr) << "\nname = " << ptr__PER->name << "\n";
    (*ss_ptr) << "\nis_hybrid = " << is_hybrid << "\n";
    
//     (*ss_ptr) << "\nmini_A_it:\n";
//     mini_A_it->second.print_this_MiniRegion(ss_ptr);
    
    it_gpu_result__AB->print_Paired_read_alignment_bundle(ss_ptr);
    
    if (is_hybrid)    
    {
// 	(*ss_ptr) << "\nmini_B_it:\n";
// 	mini_B_it->second.print_this_MiniRegion(ss_ptr);  
	
// 	(*ss_ptr) << "\naffecting ME:\n";
// 	MiniEvent_intersecting_minis_AB->second.print_this_MiniEvent(ss_ptr);
	
	(*ss_ptr) << "\nprofile_brkpt_interval = " << profile_brkpt_interval
		<< "\nhybrid_AB__breakpoint = " << hybrid_AB__breakpoint
		<< "\nhybrid_BA__breakpoint = " << hybrid_BA__breakpoint
		<< "\nabsolute_brkpt_interval__on_A = " << absolute_brkpt_interval__on_A
		<< "\nabsolute_brkpt_interval__on_B = " << absolute_brkpt_interval__on_B
		<< "\nis_NAHR = " << is_NAHR
		<< "\nmini_A_orientation_agrees_with_MiniEvent_profile = " << mini_A_orientation_agrees_with_MiniEvent_profile
		<< "\nmini_B_orientation_agrees_with_MiniEvent_profile = " << mini_B_orientation_agrees_with_MiniEvent_profile
		<< "\n\n";
		
	it_gpu_result__BA->print_Paired_read_alignment_bundle(ss_ptr);
    }
    
    (*ss_ptr) << "\n\nDONE printing PER_bundle_wrapper\n\n";

    if (!provided_with_valid_sstream)
    {
        std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );       
        std::cerr.flush();std::fflush(stderr);
    }
    
           
}//print_this_PER_bundle_wrapper











void PER_bundle_wrapper__space::PER_bundle_wrapper::process_results
			    (const bool &special_save_breakpoints_actually_captured_in_alignment)  const
{        
    
    if (!is_hybrid)
    {
	//assume only "AB" and "A" variables have been defined.
	
	const mpfr::mpreal P_read_given_location(get_likelihood_of_Viterbi_alignment(it_gpu_result__AB));
	
	const BOOST_Interval algmt_endpoints_hull(make_BI(it_gpu_result__AB->get_hull_of_mate_algmt_endpoints_on_relative_Reference()));
	
	ptr__PER->map_affected_MiniID__to__P_nohybrid.insert(type_uint__real(mini_A_it->first, P_read_given_location));
	
	ptr__PER->map_MiniID__to__nohybrid_almgt_interval_absolute.insert(
		    type_uint__BI(mini_A_it->first, BOOST_Interval(algmt_endpoints_hull.lower() + mini_A_it->second.region_interval.lower(),
								    algmt_endpoints_hull.upper() + mini_A_it->second.region_interval.lower()) ));
		
    
    }//nohybrid
    else
    {//hybrid

	//Here, we determine if a breakpoint is relevant or not.  If the breakpoint did not oiccur within the alignment of the PER to the hybrid sequence, then the alignment is (in theory) equivalent to aligning to the "no-hybrid" version of the appropriate seuqnece, which was done at the beginning of this function.  Thus, there is no need to save this probability again, and we shouldn't mark such a breakpoint as relevant.  Otherwise, we save it    
    
	const BOOST_Interval algmt_endpoints_hull__AB(make_BI(it_gpu_result__AB->get_hull_of_mate_algmt_endpoints_on_relative_Reference()));
	const BOOST_Interval algmt_endpoints_hull__BA(make_BI(it_gpu_result__BA->get_hull_of_mate_algmt_endpoints_on_relative_Reference()));
		
	const bool breakpoint_is_relevant_for_AB = BOOST_overlap(hybrid_AB__breakpoint, algmt_endpoints_hull__AB);                                                    
	const bool breakpoint_is_relevant_for_BA = BOOST_overlap(hybrid_BA__breakpoint, algmt_endpoints_hull__BA);
	
	

	if (breakpoint_is_relevant_for_AB   or   breakpoint_is_relevant_for_BA)
	{
	    const mpfr::mpreal P_hybrid_AB(get_likelihood_of_Viterbi_alignment(it_gpu_result__AB));
	    const mpfr::mpreal P_hybrid_BA(get_likelihood_of_Viterbi_alignment(it_gpu_result__BA));
	    
	    for (uint lu=0; lu<2; ++lu)
	    {
		MiniEvent_intersecting_minis_AB->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[mini_A_it->second.MiniID]
							.insert(pair_at(profile_brkpt_interval,lu));    
							
		MiniEvent_intersecting_minis_AB->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[mini_B_it->second.MiniID]
							.insert(pair_at(profile_brkpt_interval,lu));     
	    }
	    
	    
	    if (MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID) == 0)
	    {
		MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
				    .insert(type_BI__real__real(profile_brkpt_interval, type_real__real(P_hybrid_AB, P_hybrid_BA)));
	    }
	    else
	    {
		MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
				    .insert(type_BI__real__real(profile_brkpt_interval, type_real__real(P_hybrid_BA, P_hybrid_AB)));                                   
	    }
	
	
			    
	    //For heuristic:
	    if (special_save_breakpoints_actually_captured_in_alignment)
	    {		
		const Event *const Event_intersecting_minis_AB = MiniEvent_intersecting_minis_AB->second.full_Event;   //readability
		
		const uint mini_A_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID);
                const uint mini_B_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_B_it->second.MiniID);
			
		
						
		type_set_uint absolute_variational_positions__on_A;
		{//varpos on A		
		    const type_map_uint_to_uint map_MiniEvent_profile_varpos_only_to_mini_A(
                            convert_compressed_map_to_full_map(
                                    get_compressed_map_inverse_of_compressed_map(
                                        get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                Event_intersecting_minis_AB->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(mini_A_qs),
                                                mini_A_it->second.region_interval))));
				
		    absolute_variational_positions__on_A
				= extract_values_of_map_and_return_as_set<uint,uint>(map_MiniEvent_profile_varpos_only_to_mini_A);
		
		    absolute_variational_positions__on_A.insert(
		pair_at<type_set_uint>(Event_intersecting_minis_AB->absolute_coordinates_inside_of_large_gaps_of_profile, mini_A_qs).lower_bound(mini_A_it->second.region_interval.lower()),
		pair_at<type_set_uint>(Event_intersecting_minis_AB->absolute_coordinates_inside_of_large_gaps_of_profile, mini_A_qs).upper_bound(mini_A_it->second.region_interval.upper()));	
		}//varpos on A
		
				
		
		type_set_uint absolute_variational_positions__on_B;
		{//varpos on B		
		    const type_map_uint_to_uint map_MiniEvent_profile_varpos_only_to_mini_B(
			    convert_compressed_map_to_full_map(
				    get_compressed_map_inverse_of_compressed_map(
					get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
						Event_intersecting_minis_AB->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(mini_B_qs),
						mini_B_it->second.region_interval   )  ) ));
						
		   absolute_variational_positions__on_B
				= extract_values_of_map_and_return_as_set<uint,uint>(map_MiniEvent_profile_varpos_only_to_mini_B);
				
		    absolute_variational_positions__on_B.insert(
		pair_at<type_set_uint>(Event_intersecting_minis_AB->absolute_coordinates_inside_of_large_gaps_of_profile, mini_B_qs).lower_bound(mini_B_it->second.region_interval.lower()),
		pair_at<type_set_uint>(Event_intersecting_minis_AB->absolute_coordinates_inside_of_large_gaps_of_profile, mini_B_qs).upper_bound(mini_B_it->second.region_interval.upper()));						
		
		}//varpos on B
			    
			    

			    

    
		
		type_set_uint hybrid_varpos_indeces__AB[2];
		{//varpos on AB		    		    
		    get_variational_indeces_mapped_onto_hybrid_sequence(
			    mini_A_it->second,
			    absolute_brkpt_interval__on_A,    //not to be included!!
			    !mini_A_orientation_agrees_with_MiniEvent_profile,   //always want beginning of A
			    mini_B_it->second,
			    absolute_brkpt_interval__on_B,
			    !mini_B_orientation_agrees_with_MiniEvent_profile,
			    absolute_variational_positions__on_A,
			    absolute_variational_positions__on_B,
			    is_NAHR,
			    hybrid_varpos_indeces__AB[0],
			    hybrid_varpos_indeces__AB[1]);
		}//varpos on AB
			
			
		type_set_uint hybrid_varpos_indeces__BA[2];
		{//varpos on BA
		    get_variational_indeces_mapped_onto_hybrid_sequence
			    (mini_B_it->second,
			    absolute_brkpt_interval__on_B,  //profile to mini_b.   //not to be included!!
			    !mini_B_orientation_agrees_with_MiniEvent_profile,    // tricky!
			    mini_A_it->second,
			    absolute_brkpt_interval__on_A,  //inclusive
			    !mini_A_orientation_agrees_with_MiniEvent_profile,
			    absolute_variational_positions__on_B,
			    absolute_variational_positions__on_A,
			    is_NAHR,
			    hybrid_varpos_indeces__BA[1],
			    hybrid_varpos_indeces__BA[0]);
		}//varpos on BA
	
		
	
	    
		type_map_uint_to_uint absolute_positions_to_hybrid_indeces__AB;
		{//AB
		    type_vector_int absolute_labels_of_hybrid_indeces__AB;
		    
		    if (mini_A_qs == 0)
		    {	    	    
			absolute_labels_of_hybrid_indeces__AB
						= label_absolute_coordinates_on_hybrid
							(mini_A_it->second,
							absolute_brkpt_interval__on_A,
							mini_B_it->second,
							absolute_brkpt_interval__on_B,
							is_NAHR,  //else, GeneConv
							false, // AB
							Event_intersecting_minis_AB->recomb_type);
		    }
		    else
		    {
			absolute_labels_of_hybrid_indeces__AB
						= label_absolute_coordinates_on_hybrid
							(mini_B_it->second,
							absolute_brkpt_interval__on_B,
							mini_A_it->second,
							absolute_brkpt_interval__on_A,
							is_NAHR,  //else, GeneConv
							false, // AB
							Event_intersecting_minis_AB->recomb_type);
		    }
							
		    type_map_uint_to_uint::iterator insert_it_abs = absolute_positions_to_hybrid_indeces__AB.begin();
		    for (uint j=0; j < absolute_labels_of_hybrid_indeces__AB.size(); ++j)
		    {
			insert_it_abs = absolute_positions_to_hybrid_indeces__AB.insert(insert_it_abs, 
											type_uint__uint((uint)absolute_labels_of_hybrid_indeces__AB.at(j), j));
		    }
		}//AB
	
	
		type_map_uint_to_uint absolute_positions_to_hybrid_indeces__BA;	
		{//BA	
		    type_vector_int absolute_labels_of_hybrid_indeces__BA;
		    
		    if (mini_A_qs == 0)
		    {	    	    	
			absolute_labels_of_hybrid_indeces__BA
						= label_absolute_coordinates_on_hybrid
							    (mini_A_it->second,
							    absolute_brkpt_interval__on_A,
							    mini_B_it->second,
							    absolute_brkpt_interval__on_B,
							    is_NAHR,  //else, GeneConv
							    true,  //BA
							    Event_intersecting_minis_AB->recomb_type);
		    }
		    else
		    {
			absolute_labels_of_hybrid_indeces__BA
						    = label_absolute_coordinates_on_hybrid
							    (mini_B_it->second,
							    absolute_brkpt_interval__on_B,
							    mini_A_it->second,
							    absolute_brkpt_interval__on_A,
							    is_NAHR,  //else, GeneConv
							    true,  //BA
							    Event_intersecting_minis_AB->recomb_type);
		    }
							    
		    type_map_uint_to_uint::iterator insert_it_abs = absolute_positions_to_hybrid_indeces__BA.begin();
		    for (uint j=0; j < absolute_labels_of_hybrid_indeces__BA.size(); ++j)
		    {
			insert_it_abs = absolute_positions_to_hybrid_indeces__BA.insert(insert_it_abs, 
											type_uint__uint((uint)absolute_labels_of_hybrid_indeces__BA.at(j), j));
		    }							
		}//BA		    





		const BOOST_Interval mate_endpts__A(make_BI(it_gpu_result__AB->get_mate_alignment_endpoints__A()));
		const BOOST_Interval mate_endpts__B(make_BI(it_gpu_result__AB->get_mate_alignment_endpoints__B()));
		
	
		bool hybrid_pattern_is_actually_displayed_in_alignment___AB;
		{//AB
		    const bool  informs_LCR_0 = ( !get_intersection_of_set_and_interval(hybrid_varpos_indeces__AB[0], mate_endpts__A).empty()
						    or   !get_intersection_of_set_and_interval(hybrid_varpos_indeces__AB[0], mate_endpts__B).empty() );
		
		    const bool  informs_LCR_1 = ( !get_intersection_of_set_and_interval(hybrid_varpos_indeces__AB[1], mate_endpts__A).empty()
						    or   !get_intersection_of_set_and_interval(hybrid_varpos_indeces__AB[1], mate_endpts__B).empty() );
			
		    const BOOST_Interval hybrid_indeces_below_LCR_0(
						get_positions_of_hybrid_lying_BELOW_LCR(absolute_positions_to_hybrid_indeces__AB, Event_intersecting_minis_AB->LCRs[0]));
					
		    const BOOST_Interval hybrid_indeces_inbetween_LCRs(
					get_positions_of_hybrid_lying_INBETWEEN_LCRs(absolute_positions_to_hybrid_indeces__AB, 
										    Event_intersecting_minis_AB->region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs));

		    const BOOST_Interval hybrid_indeces_above_LCR_1(
					get_positions_of_hybrid_lying_ABOVE_LCR(absolute_positions_to_hybrid_indeces__AB, Event_intersecting_minis_AB->LCRs[1], 
										Event_intersecting_minis_AB->chromos[1]));

		    const bool  intersects_below_LCR_0(BOOST_overlap(mate_endpts__A,hybrid_indeces_below_LCR_0) 
							or  BOOST_overlap(mate_endpts__B,hybrid_indeces_below_LCR_0));
		    const bool  intersects_between_LCRs(BOOST_overlap(mate_endpts__A,hybrid_indeces_inbetween_LCRs)  
							or  BOOST_overlap(mate_endpts__B,hybrid_indeces_inbetween_LCRs));
		    const bool  intersects_above_LCR_1(BOOST_overlap(mate_endpts__A,hybrid_indeces_above_LCR_1)  
							or  BOOST_overlap(mate_endpts__B,hybrid_indeces_above_LCR_1));
							
		    hybrid_pattern_is_actually_displayed_in_alignment___AB
			    =   (    (informs_LCR_0  and  informs_LCR_1)
				or   (informs_LCR_0  and  (intersects_between_LCRs or intersects_above_LCR_1))
				or   (informs_LCR_1  and  (intersects_between_LCRs or intersects_below_LCR_0)) );  //then is informative of this breakpoint.	
		}//AB
		
		
		
		
		bool hybrid_pattern_is_actually_displayed_in_alignment___BA;
		{//BA
		    const bool  informs_LCR_0 = ( !get_intersection_of_set_and_interval(hybrid_varpos_indeces__BA[0], mate_endpts__A).empty()
						    or   !get_intersection_of_set_and_interval(hybrid_varpos_indeces__BA[0], mate_endpts__B).empty() );
		
		    const bool  informs_LCR_1 = ( !get_intersection_of_set_and_interval(hybrid_varpos_indeces__BA[1], mate_endpts__A).empty()
						    or   !get_intersection_of_set_and_interval(hybrid_varpos_indeces__BA[1], mate_endpts__B).empty() );

			
		    const BOOST_Interval hybrid_indeces_below_LCR_0(
						get_positions_of_hybrid_lying_BELOW_LCR(absolute_positions_to_hybrid_indeces__BA, Event_intersecting_minis_AB->LCRs[0]));
					
		    const BOOST_Interval hybrid_indeces_inbetween_LCRs(
					get_positions_of_hybrid_lying_INBETWEEN_LCRs(absolute_positions_to_hybrid_indeces__BA, 
										    Event_intersecting_minis_AB->region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs));

		    const BOOST_Interval hybrid_indeces_above_LCR_1(
					get_positions_of_hybrid_lying_ABOVE_LCR(absolute_positions_to_hybrid_indeces__BA, Event_intersecting_minis_AB->LCRs[1], 
										Event_intersecting_minis_AB->chromos[1]));

		    const bool  intersects_below_LCR_0(BOOST_overlap(mate_endpts__A,hybrid_indeces_below_LCR_0) 
							or  BOOST_overlap(mate_endpts__B,hybrid_indeces_below_LCR_0));
		    
		    const bool  intersects_between_LCRs(BOOST_overlap(mate_endpts__A,hybrid_indeces_inbetween_LCRs)  
							or  BOOST_overlap(mate_endpts__B,hybrid_indeces_inbetween_LCRs));
		    
		    const bool  intersects_above_LCR_1(BOOST_overlap(mate_endpts__A,hybrid_indeces_above_LCR_1)  
							or  BOOST_overlap(mate_endpts__B,hybrid_indeces_above_LCR_1));
							
		    hybrid_pattern_is_actually_displayed_in_alignment___BA
			    =   (    (informs_LCR_0  and  informs_LCR_1)
				or   (informs_LCR_0  and  (intersects_between_LCRs or intersects_above_LCR_1))
				or   (informs_LCR_1  and  (intersects_between_LCRs or intersects_below_LCR_0)) );  //then is informative of this breakpoint.	
		}//BA	
		
		
		if (hybrid_pattern_is_actually_displayed_in_alignment___AB  or  hybrid_pattern_is_actually_displayed_in_alignment___BA)
		{
		    MiniEvent_intersecting_minis_AB->second.relevant_breakpoints_for_which_this_PER_actually_shows_a_hybrid_pattern.insert(profile_brkpt_interval);
		}
	
	    }//special_save_breakpoints_actually_captured_in_alignment		
	}//save it            
    }//hybrid

}//process_results


































mpfr::mpreal  PER_bundle_wrapper__space::get_likelihood_of_Viterbi_alignment
			   (const CUPA::type_list_Paired_read_alignment_bundle::const_iterator &it_some_gpu_result)
{
    return   mpfr::exp10(mpfr::mpreal(it_some_gpu_result->log10Viterbi_likelihood__mate_A) + mpfr::mpreal(it_some_gpu_result->log10Viterbi_likelihood__mate_B));
}












void  PER_bundle_wrapper__space::prepare_PER_mates_for_GPU_alignment
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
			     const std::string &reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
			     std::string &prepared_read_A,
			     std::string &prepared_base_Qualities__A,
			     std::string &prepared_read_B,			     
			     std::string &prepared_base_Qualities__B,
			     std::string &prepared_reference)
{
    //NOTE:  THIS IS DIFFERENT FROM THE REGULAR PER ALIGNER!!!!!
    //The convention for the GPU is:   mate A and REF agree in complementariy and orientation, and mate B runs in opposite orientation of REF and has opposite cfomplementarity.
    //Thus, we MAY need to COMPLEMENT the REF to meet this convention.  However, in the regular (cpu verison) of PER_aligner, the REF always remains unchanged and solely thre mates are adjusted.
    //Note that our adjustment is not bad since we are ONLY potentially complementing thje REF.  Thus, we can "naively" determine the alignment endpoints on the REF, regardless of which mate was reversed and/or comeplemented and whether the ref was complemented or not.


    //prepare the mates as if the REF is fixed:
    
    bool prepared_read_A__is__FORWARD_strand =  !mate_A.IsReverseStrand();    
        
    
    prepared_read_A = mate_A.QueryBases;
    prepared_read_B = mate_B.QueryBases;
            
    //seqs:
    if (!orientation_of_PER_and_MiniRegion_agree)
    {
        prepared_read_A = get_reverse_of_sequence(prepared_read_A);
        prepared_read_B = get_reverse_of_sequence(prepared_read_B);
        std::swap<std::string>(prepared_read_A, prepared_read_B);
        
        prepared_read_A__is__FORWARD_strand = !prepared_read_A__is__FORWARD_strand;  //flip
    }

    if (complement_PER_to_get_to_MiniRegion)
    {             
        prepared_read_A = get_complement_of_sequence(prepared_read_A);
        prepared_read_B = get_complement_of_sequence(prepared_read_B);
        prepared_read_A__is__FORWARD_strand = !prepared_read_A__is__FORWARD_strand;   //flip
    }   
    
    
    //Quals:
    prepared_base_Qualities__A = mate_A.Qualities;
    prepared_base_Qualities__B = mate_B.Qualities;
    
    if (!orientation_of_PER_and_MiniRegion_agree)
    {
	prepared_base_Qualities__A = get_reverse_of_sequence(prepared_base_Qualities__A);
	prepared_base_Qualities__B = get_reverse_of_sequence(prepared_base_Qualities__B);
	std::swap<std::string>(prepared_base_Qualities__A, prepared_base_Qualities__B);                    
    }   
    
    
    //Right now, the mates have been properly ordered and complemented so that they both match the given Reference sequence.  In other words, "prepared_read_A" aligns to REF as is, and "prepared_read_B" also aligns to REF as is  (i.e. just take the strings "prepared_read_A/B" and "ref" AS IS and feed them into any old pairwise alignment program and you would get good alignments).        
    
    //Now adjust for GPU convention:  mate A and REF always agree in orientation and complementarity, and mate B is always opposite orientation and oppposite complementairty of Ref.
    //But also, since we use context-specific alignment, then it is important that both MATES "look" like they did when they were generated (recall that ".QueryBases" has always reversed and/or complemented the mate sequences as necessary so that if you take the Referecne and either mate sequence's ".QueryBases" and just naively feed them into any old pairwise aligner, then they will align well.
    
    prepared_reference = reference;
    
    if (prepared_read_A__is__FORWARD_strand)
    {
	//make A "look like" it did when it was generated:
	//do nothing.
	
	//Recall the GPU convention that mate A and REF should agree in complementairty.
	//Maintain GPU convention:
	//do nothing.
	    
	//make B "look like" it did when it was generated:
	prepared_read_B = get_reverse_complement_of_sequence___ie_inversion(prepared_read_B); // safe?
	prepared_base_Qualities__B = get_reverse_of_sequence(prepared_base_Qualities__B);
    }
    else
    {
	//make A "look like" it did when it was generated:
	prepared_read_A = get_complement_of_sequence(prepared_read_A);
	
	
	//Recall the GPU convention that mate A and REF should agree in complementairty.
	//Maintain GPU convention:
	prepared_reference = get_complement_of_sequence(prepared_reference);
		
	//make B "look like" it did when it was generated:
	prepared_read_B = get_reverse_of_sequence(prepared_read_B);
	prepared_base_Qualities__B = get_reverse_of_sequence(prepared_base_Qualities__B);
	// (Note that by complementing the REF above, we now have that mate B and REF have oppossite  complementairt - as required by GPU convention.  Thus we really do only need to reverse, NOT complement, mate B).
    }    
    
    remove_N_tails_from_read_and_quality(prepared_read_A, prepared_base_Qualities__A);
    remove_N_tails_from_read_and_quality(prepared_read_B, prepared_base_Qualities__B);
                                                                                      

}//prepare_PER_mates_for_GPU_alignment










void PER_bundle_wrapper__space::remove_N_tails_from_read_and_quality
				    (std::string &prepared_read_bases,
				     std::string &prepared_read_qualities)
{
    int last_good_base_index = prepared_read_qualities.size()-1;
    while (last_good_base_index >= 0  and  prepared_read_qualities.at(last_good_base_index)  <= BASE_QUAL_CUTOFF)
    {  --last_good_base_index;  }
    
    
    int first_good_base_index = 0;
    while (first_good_base_index < prepared_read_qualities.size()  and  prepared_read_qualities.at(first_good_base_index) <= BASE_QUAL_CUTOFF)
    {  ++first_good_base_index;  }
    
    
    const int size_of_good_subsring = last_good_base_index - first_good_base_index + 1;
    if (size_of_good_subsring > 0)
    {
	prepared_read_bases = prepared_read_bases.substr(first_good_base_index, size_of_good_subsring);
	prepared_read_qualities = prepared_read_qualities.substr(first_good_base_index, size_of_good_subsring);    
    }    

}//remove_N_tails_from_read_and_quality




void PER_bundle_wrapper__space::remove_N_tails_from_read_and_quality
				    (std::string &prepared_read_bases,
				     type_vector_int &prepared_read_qualities)
{
    int last_good_base_index = prepared_read_qualities.size()-1;
    while (last_good_base_index >= 0  and  prepared_read_qualities.at(last_good_base_index) <= BASE_QUAL_CUTOFF)
    {  --last_good_base_index;  }
    
    
    int first_good_base_index = 0;
    while (first_good_base_index < prepared_read_qualities.size()  and  prepared_read_qualities.at(first_good_base_index) <= BASE_QUAL_CUTOFF)
    {  ++first_good_base_index;  }
    
    
    const int size_of_good_subsring = last_good_base_index - first_good_base_index + 1;
    if (size_of_good_subsring > 0)
    {
	prepared_read_bases = prepared_read_bases.substr(first_good_base_index, size_of_good_subsring);
	
	type_vector_int new_quals(prepared_read_bases.size(), 0);
	uint ctr=0;
	for (int j = first_good_base_index; j <= last_good_base_index; ++j)
	{
	    new_quals.at(ctr) = prepared_read_qualities.at(j); 
	    ++ctr;
	}
	
	prepared_read_qualities = new_quals; 
    }    

}//remove_N_tails_from_read_and_quality










CUPA::type_list_Paired_read_alignment_bundle::const_iterator   PER_bundle_wrapper__space::format_and_submit_alignment_work_to_GPU_workload
                            (CUPA::GPU_alignment_workload &GPU_work,			     
			     const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
			     const std::string &reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion)
{
    std::string prepared_read_A;
    std::string prepared_base_Qualities__A;
    std::string prepared_read_B;
    std::string prepared_base_Qualities__B;       
    std::string prepared_reference;
			 
    prepare_PER_mates_for_GPU_alignment
                            (mate_A, mate_B,
			     reference,
                            orientation_of_PER_and_MiniRegion_agree,
                            complement_PER_to_get_to_MiniRegion,
			     prepared_read_A, prepared_base_Qualities__A,
			     prepared_read_B, prepared_base_Qualities__B,
			    prepared_reference);
			     
    return  GPU_work.add_alignment_work(prepared_reference,				  
					prepared_read_A, prepared_base_Qualities__A,
					prepared_read_B, prepared_base_Qualities__B,
					ASCII_quality_offset);
			     
}//format_and_submit_alignment_work_to_GPU_workload
































void  PER_bundle_wrapper__space::process_orienting_of_an_illdefined_PER
		    (PER_bundle_wrapper__space::type_map_string_to_list_PER_bundle_wrapper::const_iterator it_some_PER)
{
    assert(it_some_PER->second.size() == 4);
    
    //determine best configuration
    
    type_list_PER_bundle_wrapper::const_iterator it_algmt_result = it_some_PER->second.begin();
    
    
    
    const longdouble  normal__log10prob = (longdouble)it_algmt_result->it_gpu_result__AB->log10Viterbi_likelihood__mate_A 
					+ (longdouble)it_algmt_result->it_gpu_result__AB->log10Viterbi_likelihood__mate_B;
    ++it_algmt_result;
    
    
    const longdouble  revcomp__log10prob = (longdouble)it_algmt_result->it_gpu_result__AB->log10Viterbi_likelihood__mate_A 
					 + (longdouble)it_algmt_result->it_gpu_result__AB->log10Viterbi_likelihood__mate_B;    
    ++it_algmt_result;
    
    
    const longdouble  oppnormal__log10prob = (longdouble)it_algmt_result->it_gpu_result__AB->log10Viterbi_likelihood__mate_A 
					   + (longdouble)it_algmt_result->it_gpu_result__AB->log10Viterbi_likelihood__mate_B;    
    ++it_algmt_result;
    
    
    const longdouble  opprevcomp__log10prob = (longdouble)it_algmt_result->it_gpu_result__AB->log10Viterbi_likelihood__mate_A 
					    + (longdouble)it_algmt_result->it_gpu_result__AB->log10Viterbi_likelihood__mate_B;   
//     ++it_algmt_result;    
    
					    
					    
					    
    const longdouble*  best_config_prob = &normal__log10prob;     
    int best_config = 0;  // 0 -normal, 1- revcomp second,   2- opp normal,  3 - opp revcomp   
    
    if ( *best_config_prob < revcomp__log10prob)
    {
        best_config_prob = &revcomp__log10prob;
        best_config = 1;        
    }        
    
    if ( *best_config_prob < oppnormal__log10prob)
    {
        best_config_prob = &oppnormal__log10prob;
        best_config = 2;        
    }       
    
    if ( *best_config_prob < opprevcomp__log10prob)
    {
        best_config_prob = &opprevcomp__log10prob;
        best_config = 3;        
    }        
            
            
            
    Paired_end_read *const the_PER_ptr = it_some_PER->second.begin()->ptr__PER;    
    
    //assign
    switch (best_config)  // 0 normal, 1 revcomp second,   2 opp normal,  3 opp revcomp  
    {
	case 0:                
	    break;
	case 1:
//                 mates.second = revcomp_second;
	    the_PER_ptr->mates.second.QueryBases = get_reverse_complement_of_sequence___ie_inversion(the_PER_ptr->mates.second.QueryBases);
	    the_PER_ptr->mates.second.Qualities = get_reverse_of_sequence(the_PER_ptr->mates.second.Qualities);                
	    break;
	case 2:  
	    std::swap<BamTools::BamAlignment>(the_PER_ptr->mates.first, the_PER_ptr->mates.second); 
	    break;
	case 3:                
	    std::swap<BamTools::BamAlignment>(the_PER_ptr->mates.first, the_PER_ptr->mates.second);    
//                 mates.first = revcomp_second;
	    the_PER_ptr->mates.first.QueryBases = get_reverse_complement_of_sequence___ie_inversion(the_PER_ptr->mates.first.QueryBases);
	    the_PER_ptr->mates.first.Qualities = get_reverse_of_sequence(the_PER_ptr->mates.first.Qualities);
	    break;      
	default:
	{
	    std::stringstream error_strm;       
	    error_strm << "\n\nERROR!    best_algmt_config_first_MR   =  "  <<   best_config   << "!!!\n\n";                
	    the_PER_ptr->print_this_PER(&error_strm);                
	    error_message(error_strm, false);
	    break;            
	}
    }//switch    
    
    
}//process_orienting_of_an_illdefined_PER



