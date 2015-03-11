#include <stdlib.h>
#include <math.h>
#include <ctype.h>/*
#include <assert.h>*/
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

#include <mpreal.h>
#include <mpfr.h>



#include <api/BamAlignment.h>




#include <general_typedefs.h>
#include <SD_Entry_general.h>
#include <Interlocking_entry.h>
#include <templates.h>
#include <translations_through_alignment.h>
#include <Sparse_map.h>
#include <Readgroup_statistics.h>


#include "globals.h"
#include "Paired_end_read.h"
#include "other_functions.h"
#include "Event.h"
#include "MiniRegion.h"
#include "PER_aligner.h"
#include "Conn_Comp.h"


#include "io_functions.h"

#include "Star_alignment.h"

#include <CUPA.h>
#include "PER_bundle_wrapper.h"

#include "MiniEvent.h"








void Paired_end_read::print_this_PER(std::stringstream *const &some_ss)  const
{
    
    std::stringstream output_str;
    
    const bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str; 

    print_line_of_markers("[[", ss_ptr);

    (*ss_ptr) << "\nprinting PER:\n\tname =  " <<  name << "\n";
    print_Balgmt(mates.first, ss_ptr);
    print_Balgmt(mates.second, ss_ptr);
//     (*ss_ptr) << "mate_sequences.first = " << mate_sequences.first << "\nmate_sequences.second = " << mate_sequences.second << "\n";

    print_map_keys<uint, MiniRegion>(my_MiniRegions, "my_MiniRegions", ss_ptr);

    for (type_map_uint_to_MiniRegion::const_iterator it = my_MiniRegions.begin();
            it != my_MiniRegions.end();
            ++it)
    {  it->second.print_this_MiniRegion(ss_ptr);  }

    (*ss_ptr) << "\n\n\n\n";
    print_map_keys_and_values<uint, real>(map_affected_MiniID__to__P_nohybrid, "map_affected_MiniID__to__P_nohybrid", ss_ptr);
    
    (*ss_ptr) << "\t\t\t\thaploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal =   " 
            <<  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal.first  <<  ",  " 
            << haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal.second  << "\n\n";
    
    (*ss_ptr) << "\t\t\t\thaploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome =   " 
            <<  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome.first  <<  ",  " 
            << haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome.second  << "\n\n";
            
    (*ss_ptr) << "\t\t\t\thaploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome =   " 
            <<  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome.first  <<  ",  " 
            << haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome.second  << "\n\n";            
          
    print_map_keys<uint, MiniEvent>(my_MiniEvents, "my_MiniEvents", ss_ptr);

    for (type_map_uint_to_MiniEvent::const_iterator it = my_MiniEvents.begin();
            it != my_MiniEvents.end();
            ++it)
    {  it->second.print_this_MiniEvent(ss_ptr);  }
    
    
    
//     for (type_map_uint_to_MiniEvent::const_iterator it = my_MiniEvents.begin();
//             it != my_MiniEvents.end();
//             ++it)   
//     {
//         std::fprintf(stderr, "\n\nMiniEvent %u:\n", it->first );
//         print_map_keys_and_pair_values<uint, real, real>(it->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA,
//                                                          "it->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA");
//     }        


    print_line_of_markers("]]", ss_ptr);

    if (!provided_with_valid_sstream)
    {
        std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );       
        std::cerr.flush();std::fflush(stderr);
    }
    
    

}  // end of   print_this_PER

















void Paired_end_read::init_PER
		     (const BamTools::BamAlignment &some_Balgmt,
                     const BamTools::BamAlignment &another_Balgmt)
{
    mates.first = some_Balgmt;
    mates.second = another_Balgmt;    
              
    const type_map_string_to_Readgroup_stats::const_iterator it_find_rg = identify_Readgroup_from_name_of_PER(some_Balgmt.Name);
        
    
    if (it_find_rg == Readgroup_stats.end())
    {
		failed_PER_flag = true;
		my_Readgroup_ptr = NULL;
		return;
    }
    else
    {
		failed_PER_flag = false;
		my_Readgroup_ptr = &(it_find_rg->second);
    }
    
    
    for (uint m=0;m<2;++m)
    {
	if (pair_at<BamTools::BamAlignment>(mates,m).QueryBases.empty()
	    or  pair_at<BamTools::BamAlignment>(mates,m).Qualities.empty())
	{
	    failed_PER_flag = true;  
	    return;	    
	}
	else
	{
	    const int qualcutoff_as_appropriate_char = (int)(ASCII_quality_offset + BASE_QUAL_CUTOFF);
	    uint number_of_decent_bases = 0;
	    for (uint j=0; j < pair_at<BamTools::BamAlignment>(mates,m).Qualities.size(); ++j)
	    {
		if ((int)pair_at<BamTools::BamAlignment>(mates,m).Qualities[j] > qualcutoff_as_appropriate_char)
		{
		    ++number_of_decent_bases;
		    if (number_of_decent_bases >= 20)
		    {  break;  }
		}
	    }
	    
	    if (number_of_decent_bases < 20)
	    {
		failed_PER_flag = true;
		return;
	    }
	}
    }
    
}//init_PER














uint Paired_end_read::determine_observed_frag_length()  const
{    
    return   (std::max<int>(mates.first.Position + mates.first.Length - 1, mates.second.Position + mates.second.Length - 1)
	     - std::min<int>(mates.first.Position, mates.second.Position));
    
} // determine_observed_frag_length





















void Paired_end_read::orient_and_complement_each_mate_according_to_its_spawning_Event_profile()
{
//Thus, after calling this function, it no longer matters if a read was "unmapped" in the SAM/BAMfile.  In this function, we performed both possible aliognments to the profile and determined which one was better and called that better one the real alignment.
    
    
    //First, orient the PER to its spawning MiniRegion.
                //After, (at the end of this function), we will then orient the PER according to the spawning Event profile.
    
    if (mates.first.IsProperPair())  // then only has one MR, which is the spawning MR.
    {                
        //for consistency:
        if (mates.first.Position > mates.second.Position)
	{  std::swap<BamTools::BamAlignment>(mates.first, mates.second);  }
        
        
        if (mates.first.IsReverseStrand()  ==  mates.second.IsReverseStrand())
        {
            std::stringstream warning_strm;
            print_line_of_markers("WARNING! ", &warning_strm);
            print_line_of_markers("(", &warning_strm);
            warning_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
            
            warning_strm << "\n\nWARNING!   properly paired read   " << name << "  has same strandedness for mates!!!\n\n";                        
            
            print_line_of_markers(")", &warning_strm);
            std::fprintf(stderr, "\n\n%s\n\n", warning_strm.str().c_str() );
            
            
            mates.second.SetIsReverseStrand( !mates.first.IsReverseStrand() );                    
        }           
    }
    else  //.IsPaired, but NOT .IsProperPair as far as I am concerned
    {           
	//do nothing.  already did the work on  the GPU.             
    }
                                
    
    adjust_orientation_and_complementairy_of_PER_mates_according_its_spawning_Event_profile();


} // end of  orient_and_complement_each_mate_according_to_its_spawning_MiniRegion










void Paired_end_read::adjust_orientation_and_complementairy_of_PER_mates_according_its_spawning_Event_profile()
{
    //now fix orientation/complementarity, so that the PER is aligned to the spawning Event profile.
    if (!my_MiniRegions.begin()->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree)
    {
        mates.first.QueryBases = get_reverse_of_sequence(mates.first.QueryBases);
        mates.first.Qualities = get_reverse_of_sequence(mates.first.Qualities);
        
        mates.second.QueryBases = get_reverse_of_sequence(mates.second.QueryBases);
        mates.second.Qualities = get_reverse_of_sequence(mates.second.Qualities);
        
        std::swap<BamTools::BamAlignment>(mates.first, mates.second);        
    }
    
    if (my_MiniRegions.begin()->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile)
    {        
        mates.first.SetIsReverseStrand(!mates.first.IsReverseStrand());
        mates.second.SetIsReverseStrand(!mates.second.IsReverseStrand());
        
        mates.first.QueryBases =  get_complement_of_sequence(mates.first.QueryBases);
        mates.second.QueryBases = get_complement_of_sequence(mates.second.QueryBases);     
    }                  
}//adjust_orientation_and_complementairy_of_PER_mates_according_its_spawning_Event_profile
















// int Paired_end_read::determine_best_alignment_configuration_for_mates_according_to_this_MR
//                                     (MiniRegion &some_MR)
// {
//     
//     if (  !( mates.first.RefID   ==   map_chromosome_value_to_BAM_Ref_IDs.at(some_MR.chromosome_of_region)
//             and  BOOST_in(mates.first.Position+1, some_MR.region_interval)     )        )
//         return return_best_alignment_configuration_for_mates
//                                 (mates.first,
//                                  mates.second,
//                                  some_MR);
//     else
//         return return_best_alignment_configuration_for_mates
//                                 (mates.second,
//                                  mates.first,
//                                  some_MR);        
// 
// }




void Paired_end_read::submit_PER_orienting_work_to_GPU
			    (const MiniRegion &some_MR,
			     CUPA::GPU_alignment_workload &GPU_work,
			    std::map<std::string, std::list<PER_bundle_wrapper__space::PER_bundle_wrapper> >  &record_of_PER_algmt_jobs)
			    
{            
    
    if (!mates.first.IsProperPair())
    {    
	mates.second.SetIsReverseStrand(!mates.first.IsReverseStrand());  //for illumina, always must be opposite.   first mate is always properly placed                 

	
	BamTools::BamAlignment revcomp_second(mates.second);
	revcomp_second.QueryBases = get_reverse_complement_of_sequence___ie_inversion(revcomp_second.QueryBases);
	revcomp_second.Qualities = get_reverse_of_sequence(revcomp_second.Qualities);
// 	revcomp_second.SetIsReverseStrand(!mates.second.IsReverseStrand());    
	
	
	
	//submit work:
	CUPA::type_list_Paired_read_alignment_bundle::const_iterator it_submitted_work;
	
	
	//normal
	it_submitted_work = PER_bundle_wrapper__space::format_and_submit_alignment_work_to_GPU_workload
									(GPU_work,
									mates.first, mates.second,
									some_MR.region_sequence,
									true, false);
									
	record_of_PER_algmt_jobs[name].push_back(PER_bundle_wrapper__space::PER_bundle_wrapper(
							it_submitted_work,  this));
	

	//revcomp
	it_submitted_work = PER_bundle_wrapper__space::format_and_submit_alignment_work_to_GPU_workload
									(GPU_work,
									mates.first, revcomp_second,
									some_MR.region_sequence,
									true, false);
									
	record_of_PER_algmt_jobs[name].push_back(PER_bundle_wrapper__space::PER_bundle_wrapper(
							it_submitted_work,  this));
	
	
	//opp normal
	it_submitted_work = PER_bundle_wrapper__space::format_and_submit_alignment_work_to_GPU_workload
									(GPU_work,
									mates.second, mates.first,
									some_MR.region_sequence,
									true, false);
									
	record_of_PER_algmt_jobs[name].push_back(PER_bundle_wrapper__space::PER_bundle_wrapper(
							it_submitted_work,  this));
	
	
	//opp revcomp
	it_submitted_work = PER_bundle_wrapper__space::format_and_submit_alignment_work_to_GPU_workload
									(GPU_work,
									revcomp_second, mates.first,
									some_MR.region_sequence,
									true, false);
									
	//record wrapper:
	record_of_PER_algmt_jobs[name].push_back(PER_bundle_wrapper__space::PER_bundle_wrapper(
							it_submitted_work,  this));      
    }
    
}//submit_PER_orienting_work_to_GPU












































type_map_uint_to_uint_to_set_uint Paired_end_read::create_MiniEvents_and_identify_gender_of_MRs
                                                                        (const Event *const &spawning_Event)
{

//     std::fprintf(stderr, "\n%s.create_MiniEvents_and_populate_variational_positions_for_checking", name.c_str() );

    
//     std::fprintf(stderr, "\t\t\t\tpopulating variational positioins according to universal variational positions.");
//     for (type_map_uint_to_MiniRegion::iterator it_mini = my_MiniRegions.begin();
//             it_mini != my_MiniRegions.end();
//             ++it_mini)            
//         it_mini->second.my_variational_positions_for_checking
//                     = get_intersection_of_set_and_interval( universal_variational_positions.at(it_mini->second.chromosome_of_region),
//                                                             it_mini->second.region_interval  );
    
    {//gender
        type_set_uint::iterator insert_it_XXXX = set_of_MiniRegions_on__XXXX_chromosome.begin();
        type_set_uint::iterator insert_it_YYYY = set_of_MiniRegions_on__YYYY_chromosome.begin();
                                                        
        for (type_map_uint_to_MiniRegion::const_iterator it_MR = my_MiniRegions.begin();
                it_MR != my_MiniRegions.end();
                ++it_MR)         
            switch (it_MR->second.chromosome_of_region)
            {
                case 23:
                    insert_it_XXXX = set_of_MiniRegions_on__XXXX_chromosome.insert(  insert_it_XXXX,  it_MR->first );
                    break;
                case 24:
                    insert_it_YYYY = set_of_MiniRegions_on__YYYY_chromosome.insert(  insert_it_YYYY,  it_MR->first );
                    break;     
                default:
                    break;                        
            };
    }//gender
    
    
    
    
    
    //convert all Interlocking_entries to pointers - we don't use anything else from them anyway.
    const uint relevant_Events_length = spawning_Event->local_interlocking_Events.size() + spawning_Event->global_interlocking_Events.size() + 1;    
    const Event* relevant_Events[relevant_Events_length];
    {//scope
                    uint j_insert = 0;
                    const type_map_uint_to_Interlock* local_or_global[2];
                    local_or_global[0] = &(spawning_Event->local_interlocking_Events);
                    local_or_global[1] = &(spawning_Event->global_interlocking_Events);

                    for (uint lg=0; lg<2; ++lg)
                        for (type_map_uint_to_Interlock::const_iterator it_interlock = local_or_global[lg]->begin();
                                it_interlock != local_or_global[lg]->end();
                                ++it_interlock)
                            relevant_Events[j_insert++] = it_interlock->second.entry;

                    relevant_Events[j_insert++] = spawning_Event;

                    if ( j_insert != relevant_Events_length)
                    {
                        std::fprintf(stderr, "\n\n\n\nERROR:  you can't count!   j_insert = %u != %u =  relevant_Events_length     in  \"identify_all_relevant_interlocking_Events\"\n\n\n\n",
                                    j_insert,
                                    relevant_Events_length);
                        exit(1);
                    }
    }//scope






//For each interlocking event, scroll through a PER's MiniRegions and decide if either of its LCRs intersect any of these MiniRegions.  Record information when it does.

    type_map_uint_to_uint_to_set_uint MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY;    
    type_map_uint_to_MiniEvent::iterator insert_ME_it = my_MiniEvents.begin();    


    for (uint j=0; j<relevant_Events_length; ++j)
    {

        const Event *const event_being_considered = relevant_Events[j];

        
        MiniEvent new_ME( event_being_considered, this);

        type_set_uint::iterator affected_MiniID_insert_it = new_ME.affected_MiniRegions.begin();
        
        type_map_uint_to_uint::iterator MiniID_to_LCR_insert_it = new_ME.map_MiniRegions_ON_LCRs__to__intersecting_LCR.begin();
        type_set_uint::iterator MiniID_in_DAR_insert_it = new_ME.MiniRegions_in_Inbetween_of_this_MiniEvent.begin();

        
        
        for (type_map_uint_to_MiniRegion::iterator it_mini_A = my_MiniRegions.begin();
             it_mini_A != my_MiniRegions.end();
             ++it_mini_A)
        {
                                                                        
            //else, check if this MiniRegion mostly intersects one of the LCRs of the MiniEvent.
            for (uint ev_qs = 0; ev_qs < 2; ++ev_qs)
            {

                if ( it_mini_A->second.chromosome_of_region != event_being_considered->chromos[ev_qs] )
                    continue;

                
                const BOOST_Interval Mini_ev_qs_absolute_intersection(
                                        BOOST_intersect(it_mini_A->second.region_interval, event_being_considered->LCRs[ev_qs]) );
                
                //and actually intersecting:
                if ( BOOST_empty(Mini_ev_qs_absolute_intersection) )
                    continue;



                
                type_map_uint_to_uint MR_absolute_coords_to_relevant_profile_positions(
                                    convert_compressed_map_to_full_map(
                                        get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                event_being_considered->compressed_map_LCR_to_profile__for_each_LCR.at(ev_qs),
                                                Mini_ev_qs_absolute_intersection   )  ));
                                            
                if ( MR_absolute_coords_to_relevant_profile_positions.empty() )                            
                    continue;
                    
                
                const type_set_uint relevant_positions_on_profile(
                                            extract_values_of_map_and_return_as_set<uint, uint>(MR_absolute_coords_to_relevant_profile_positions) );        
                MR_absolute_coords_to_relevant_profile_positions.clear();
                
                
                
                MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY[event_being_considered->UID]
                                                    .insert(  type_uint__set_uint(it_mini_A->first, relevant_positions_on_profile)  );                                                    
                                                            
                affected_MiniID_insert_it = new_ME.affected_MiniRegions.insert( affected_MiniID_insert_it,  it_mini_A->first );                

                MiniID_to_LCR_insert_it = new_ME.map_MiniRegions_ON_LCRs__to__intersecting_LCR
                                                            .insert( MiniID_to_LCR_insert_it,  type_uint__uint(it_mini_A->first, ev_qs)  );                
                                                                               
          
                                        
                //We DO NOT worry about relevant breakpoints yet, since "relevant" is determined by the mapping of the PER, and the PER isn't mapped until the 'next' step!  But we get a preliminary idea of the breakpoints for testing...
                

//                 std::fprintf(stderr, "MR_coords_to_relevant_breakpoints.size() = %u...", MR_absolute_coords_to_relevant_profile_positions.size() );

                //new_ME.relevant_profile_breakpoints.insert( relevant_breakpoints.begin(), relevant_breakpoints.end() );
        //         new_ME.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[it_mini_A->first]
        //                                                                     .insert( relevant_breakpoints.begin(), relevant_breakpoints.end() );
                                                   
                //We just recorded, in a MiniEvent, the breakpoints of that Event on  a MiniRegion.  That will be used for creating Hybrid LCRs.

//                 std::fprintf(stderr, "done.\n");
            }  //end for-loop  ev_qs < 2 
            
            
            //check intervening sequence...(if approrpriate)
            if (    (event_being_considered->recomb_type == recomb_class__Inv  or   event_being_considered->recomb_type == recomb_class__DupDel)
                    and   new_ME.affected_MiniRegions.count(it_mini_A->first) == 0   //i.e. doesn't intersect one of the LCRs
                    and   it_mini_A->second.chromosome_of_region == event_being_considered->chromos[0]                              
                    and   BOOST_overlap(it_mini_A->second.region_interval,
                                       event_being_considered->region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs)    )
                {                       
                    affected_MiniID_insert_it = new_ME.affected_MiniRegions.insert( affected_MiniID_insert_it, it_mini_A->first);
                    
                    MiniID_in_DAR_insert_it = new_ME.MiniRegions_in_Inbetween_of_this_MiniEvent.insert( MiniID_in_DAR_insert_it, it_mini_A->first);                    
                }                        
                        

        } // end for-loop   it_mini_A




        
//         if ( !new_ME.affected_MiniRegions.empty() )
//         {            
//             my_MiniEvents.insert( type_uint__MiniEvent(event_being_considered->UID, new_ME) );
//             std::fprintf(stderr, "\n\t\t\t\tMiniEvent %u intersects Miniregions:\n", new_ME.full_Event->UID);
//             new_ME.print_this_MiniEvent();
//         }


//         // error-checking:
//         if (   ! new_ME.affected_MiniRegions.empty() 
//                 &&  new_ME.affected_MiniRegions.size() > 2 )
//         {
//             print_line_of_markers("WARNING! ");
//             print_line_of_markers("(");
//             std::fprintf(stderr, "WARNING: interlocking event UID %u intersects %u mini-regions of event %u\n",
//                          event_being_considered->UID,
//                          new_ME.affected_MiniRegions.size(),
//                          spawning_Event->UID);
//             new_ME.print_this_MiniEvent();     
//             //exit(1);
//             print_line_of_markers(")");
//         }




        if ( !new_ME.affected_MiniRegions.empty() )
            insert_ME_it = my_MiniEvents.insert( insert_ME_it,   type_uint__MiniEvent(new_ME.full_Event->UID, new_ME) );


    }  // end for-loop     it_interlock = my_event->local/global_interlocking_Events

    
    return MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY;

} //  end of    create_MiniEvents_and_populate_variational_positions_for_checking
























void Paired_end_read::make_partners_out_of_MiniRegions_for_each_MiniEvent
                           (const type_map_uint_to_uint_to_set_uint &MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY)
{        
// Note that the same MiniRegion may be partner to multiple MiniEvents.  This is not a problem at all, it is normal and expected - it will happen for locally interlocking MiniEvents.      


    for (type_map_uint_to_uint_to_set_uint::const_iterator it_ME = MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY.begin();
            it_ME != MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY.end();
            ++it_ME)
    {

        const type_map_uint_to_MiniEvent::iterator  it_the_ME =   my_MiniEvents.find(it_ME->first);

        //for each MiniRegion, check every other MiniRegion to see if it is a partner.
        for (type_map_uint_to_set_uint::const_iterator it_MR_A = it_ME->second.begin();
                it_MR_A != it_ME->second.end();
                ++it_MR_A)
        {                                   
            const BOOST_Interval MR_A_profile_interval(   *(it_MR_A->second.begin()) ,  *(--it_MR_A->second.end())     );
            const uint qs_of_intersecting_LCR_of_MR_A__ptr = it_the_ME->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(it_MR_A->first);
            
            for (type_map_uint_to_set_uint::const_iterator it_MR_B = it_ME->second.begin();
                    it_MR_B != it_ME->second.end();
                    ++it_MR_B)
	    {
                if ( it_MR_B->first != it_MR_A->first   //skip yourself, of course.
                     and     qs_of_intersecting_LCR_of_MR_A__ptr  !=  it_the_ME->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(it_MR_B->first) ) 
                                                    //must be on opposite LCRs in order to form a hybrid during recombination...                                                                    
		{
                    if ( BOOST_overlap(MR_A_profile_interval,   BOOST_Interval( *(it_MR_B->second.begin()), *(--it_MR_B->second.end())  )   )    )
                    {
                        it_the_ME->second.partners_for_MiniRegions_on_LCRs[it_MR_A->first].insert( it_MR_B->first );   
                        it_the_ME->second.partners_for_MiniRegions_on_LCRs[it_MR_B->first].insert( it_MR_A->first );   
                    }
		}
	    }                                                                
//                 std::fprintf(stderr, "\t\t\t\tMiniRegions %u and %u are partners, with size of intersection = %u\n",
//                                     it_MR_A->first,
//                                     it_MR_B->first,
//                                     BOOST_width_inclusive(MR_A_B_profile_brkpt_interval_intersection)   );
                                    
//                 print_set<uint>(it_MR_A->second, "it_MR_A intersected_breakpoints");                    
//                 print_set<uint>(it_MR_A->second, "it_MR_B intersected_breakpoints");                                                                                                                                                                                                                                                                                                                                                                                                      
                        
        }  // end for-loop  it_MR_A






  
        // MiniRegions have been paired up.
        
        
        
        
//         // Every MiniRegion should have AT MOST  one other partner.  Adjust if necesary                
// // Note: there are two reason why a MiniRegion may have more than one partner:                
// //              1)   There is self-similaritiy WITHIN these LCRs.  This happens in dense mosaic / tandem repeat regions, where these LCRs may actually be compositions of smaller LCRs.
// //              2)   There is a "large" gap in one LCR wrt the other.  Consider the following:
// 
// 
// //                     PER:                    mate_0                    mate_1
// //                  LCR  A:   =======================------------------========================
// //                             ...                  a                  a+1
// //                 profile:    ...                  p                 p+100 
// //                             ...                  b                 b+100
// //                  LCR  B:   =================================================================
// 
// // This result in the following MiniRegions for the PER shown:
// 
// //                     PER:                    mate_0                    mate_1
// //                  LCR  A:   =======================------------------========================
// //                            MR_0        [                                       ]
//     
//     
// //                  LCR  B:   =================================================================
// //                            MR_1       [              ]       [              ]        MR_2
// 
// //It follows from the pictures.
// //Note that nothing in the picture is discordant or unusual.  All three MiniRegions are the same length!
// // The PER appears discordant only becuase of the inclusion of the gap ( "--------" )  in LCR A.  Of course, that gap is not actually there it is only their because we are showing the alignment of LCRs A and B!!!!).  The same goes for MiniREgion 0.  It looks longer because it includes the gap.  But MR_0,1,2 all have the same amount of ABSOLUTE genomic sequence.  However, because of the gap, MR_0 happens to span a lot more of the profile (thus, MR_0 contains 100 more positons of the profile than either MR_1 or MR_2).
// 
// // Therefore, a good strategy would be to weld MR_1,MR_2 together...
//         
//         
//         
//         
        
        
        //Now we check if any MiniRegion doesn't have a partner. 
    
    
        //is anybody partner-less?
        type_set_uint::iterator insert_insertion_it = it_the_ME->second.MiniRegions_in_insertion_of_an_LCR.begin();                
        
        type_set_uint::iterator it_MR = it_the_ME->second.affected_MiniRegions.begin();        
        while (  it_MR != it_the_ME->second.affected_MiniRegions.end()  )
	{
            if ( it_the_ME->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(*it_MR) > 0 ) //of course MiniRegions in the DAR don't have a partner!!!
                ++it_MR;                               
            else if (it_the_ME->second.partners_for_MiniRegions_on_LCRs.count( *it_MR ) > 0)   // already has a parnter.
                ++it_MR;
            else
            {                             
                const type_map_uint_to_MiniRegion::const_iterator  it_the_MR = my_MiniRegions.find(*it_MR);
                
                //check if intersects the DAR:
                if (  (it_the_ME->second.full_Event->recomb_type == recomb_class__Inv  or   it_the_ME->second.full_Event->recomb_type == recomb_class__DupDel)
                      and   it_the_MR->second.chromosome_of_region == it_the_ME->second.full_Event->chromos[0]                              
                      and   BOOST_overlap( it_the_MR->second.region_interval, it_the_ME->second.full_Event->region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs)    )
                {                                               
                    
                    it_the_ME->second.MiniRegions_in_Inbetween_of_this_MiniEvent.insert( *it_MR );                                                        
                    it_the_ME->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.erase(*it_MR);
                    
    //                     std::fprintf(stderr, "\t\t\t\t\tME %u:    MR %u switched to   MiniRegions_in_Inbetween_of_this_MiniEvent.\n",
    //                                  it_ME->first, *it_MR);
                    ++it_MR;                 
                }
                else        //else, this MiniRegion must be in an insertion of an LCR.
                {
                    bool intersects_an_LCR_insertion = false;     
                    const Event *const Event_ptr = it_the_ME->second.full_Event;
                    assert( Event_ptr != NULL );
                    
                    for (uint qs=0; qs<2; ++qs)
		    {
                        if ( Event_ptr->chromos[qs] == it_the_MR->second.chromosome_of_region
                             and   BOOST_overlap( it_the_MR->second.region_interval, Event_ptr->LCRs[qs])    )
                        {
                            
                            const type_set_uint profile_pos_intersected_by_this_MR(
                                        extract_keys_of_compressed_map_and_return_as_set(
                                            get_compressed_map_inverse_of_compressed_map(  //profile to LCR
                                                    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                            Event_ptr->compressed_map_LCR_to_profile__for_each_LCR.at(qs),
                                                            it_the_MR->second.region_interval     ) )));                                                                                                
                            
                            if ( ! profile_pos_intersected_by_this_MR.empty() )                            
                            {                            
                                const BOOST_Interval profile_pos_endpoints_intersected_by_this_MR(
                                                                *profile_pos_intersected_by_this_MR.begin(),
                                                                *(--profile_pos_intersected_by_this_MR.end())    ); 
                                                            
                                for (type_set_BI::const_iterator it_large_gap = Event_ptr->large_gaps_of_profile.begin();
                                        it_large_gap != Event_ptr->large_gaps_of_profile.end();
                                        ++it_large_gap)                            
                                    if ( BOOST_overlap( *it_large_gap, profile_pos_endpoints_intersected_by_this_MR ) )
                                    {
                                        insert_insertion_it = it_the_ME->second.MiniRegions_in_insertion_of_an_LCR.insert( insert_insertion_it, *it_MR );
                                        intersects_an_LCR_insertion = true;
                                        break;                            
                                    }
                                    
                                if (intersects_an_LCR_insertion)    
                                    break;
                            }
                        }
		    }//qs
                        
                        
                    if (intersects_an_LCR_insertion)    
                        ++it_MR;
                    else  
                    {
                        //then this MiniRegion must be "on the edge" of this MiniEvent.  so delete it.  we are ignoring such edge effects.
                        it_the_ME->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.erase(*it_MR);
                        it_the_ME->second.affected_MiniRegions.erase(it_MR++);
    //                     print_line_of_markers("ERROR! ");
    //                     print_line_of_markers("(");
    //                     
    //                     std::fprintf(stderr, "\n\nThis MiniRegions allegedly intersects this MiniEvent, but does not intersect the Inbetween, niether LCR, nor an insertion in an LCR.\n\n");
    //                     
    //                     it_the_MR->second.print_this_MiniRegion();                    
    //                     it_the_ME->second.print_this_MiniEvent();                    
    //                     Event_ptr->print_this_entry();
    //                     
    //                     print_line_of_markers(")");
                    }                                                                                 
                }   // else, check if is in insertion in an LCR.
            }  // else, check DAr & insertion in LCR
	}//affected_MiniRegions

    }  //  for-loop  it_ME
    
    
    
    
    
    
    
    
    // erase any empty MiniEvents.  (may no longer be relevant after edge effects were filtered out.)
    type_map_uint_to_MiniEvent::iterator it_ME_check = my_MiniEvents.begin();
    while ( it_ME_check != my_MiniEvents.end() ) 
    {
        if (it_ME_check->second.affected_MiniRegions.empty())
            my_MiniEvents.erase(it_ME_check++);   //post-increment necessary!
        else
            ++it_ME_check;        
    }

            
                     

} // end of    make_partners_out_of_MiniRegions_for_each_MiniEvent


































//The following serves both a computational advantage and an implementation practicality.
//There may be some regions which are unaffected by any Events (although they may still have positions that are homologous to some Event's variational positions).  Thus, anytime we ever sum over all possible mappings of this PER, the probability of its mapping to these unaffected MiniRegions will always be the same.
//Thus we precompute this part of the sum in "haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents.
//Also, we will be given state vectors wherein some of the MiniRegions have hybridized.  In these cases, we will also need to know which MiniRegions did NOT hybridize.  thus we need a map of ONLY AFFECTED miniRegions.

//now separate out those which may not be affected by any Events - these MiniRegions will always contribute the same value to the sum over all mappings for this PER.

void Paired_end_read::sum_over_all_unaffected_LCRs()
{
      
    type_set_uint all_affected_MiniRegions;

    for (type_map_uint_to_MiniEvent::const_iterator it_ME = my_MiniEvents.begin();
	    it_ME != my_MiniEvents.end();
	    ++it_ME)
    {
	all_affected_MiniRegions.insert( it_ME->second.affected_MiniRegions.begin(), it_ME->second.affected_MiniRegions.end() );
    }

    //haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents;  //iniitalized to 0.00  in constructor.        
    type_map_uint_to_real::iterator it_MR_check_nohybrids = map_affected_MiniID__to__P_nohybrid.begin();  
					    //remember, this map was initialied to contain no-hybrids for ALL MiniRegions (not just for affected MiniRegions) .
    while (it_MR_check_nohybrids != map_affected_MiniID__to__P_nohybrid.end())  
    {
	if ( all_affected_MiniRegions.count(it_MR_check_nohybrids->first) == 0 )
	{     
	    switch ( my_MiniRegions.at(it_MR_check_nohybrids->first).chromosome_of_region )
	    {
		case 23:  // X
		    haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome.second += it_MR_check_nohybrids->second;      
		    ++haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome.first;
		    break;
		case 24:  // Y
		    haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome.second += it_MR_check_nohybrids->second;      
		    ++haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome.first;
		    break;       
		default:   //autosomal
		    haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal.second += it_MR_check_nohybrids->second;      
		    ++haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal.first;
		    break;                                                
	    };//switch
	    
	    map_affected_MiniID__to__P_nohybrid.erase(it_MR_check_nohybrids++);   //posti-increment is NECESSARY!!!
	}    
	else
	{  ++it_MR_check_nohybrids;  }
    }
	    
    all_affected_MiniRegions.clear();

}//sum_over_all_unaffected_LCRs














void  Paired_end_read::align_to_every_possible_hybrid_outcome_of_MiniRegions_and_Events_and_record_probabilities_and_relevant_breakpoints
			    (CUPA::GPU_alignment_workload &GPU_work,
			     PER_bundle_wrapper__space::type_map_string_to_list_PER_bundle_wrapper  &record_of_PER_algmt_jobs,
			     const bool &consider_GeneConversion_breakpoints)
{   
//This is the big one...

// For each MiniEvent:
//      1)  we figure out which (pair of) MiniRegions would recombine if this Event occurred, and which profile breakpoints are contained in these MiniRegions
//      2)  For each profile breakpoint that is contained in this pair of MiniRegions,  we construct each type of hybrid (AB and BA), align the PER to these hybrids, and calcucate the probability of the mapping (i.e. alignment) of this PER to these hybrids.
//      3)  Given the alignment of the PER to the hybrid, we check where the breakpoint of the hybrid falls wrt the alignment.  If the breakpoint is contained contained within the interval formed by the ends of the alignment, then this alignment should really be equivalernt to the alignment in the case of "no hybrid".  In this case, there would be no reason to save this breakpoint or the calculated probabilities - they are simply not relevant.


    //First, compute the probability of the PER mapping to each MiniRegion.  Note that we call this map "affected_MiniID__to__P_nohybrid", but right now it is actually all MiniIDs (those potentially affected by MiniEvents and those not).  This will be corrected in the next step.



    //first we record no-hybrid for ALL MiniRegions.  Next we will restrict this list to only affected MiniRegions.

    for (type_map_uint_to_MiniRegion::const_iterator it_mini_A = my_MiniRegions.begin();
            it_mini_A != my_MiniRegions.end();
            ++it_mini_A)
    {
        //no-hybrid case
	
	//submit work:
	const CUPA::type_list_Paired_read_alignment_bundle::const_iterator it_submitted_work__nohybrid
		    = PER_bundle_wrapper__space::format_and_submit_alignment_work_to_GPU_workload
                            (GPU_work,			     
			     mates.first, mates.second,
			     it_mini_A->second.region_sequence,
			    it_mini_A->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree,
			    it_mini_A->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile);
				
	//record wrapper:
	record_of_PER_algmt_jobs[name].push_back(PER_bundle_wrapper__space::PER_bundle_wrapper(
						    it_submitted_work__nohybrid,  this,
						    it_mini_A, it_mini_A->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree));

    } // end of no-hybrid









    
    
    

    
    
    
    //align to all hybrids        
    for (type_map_uint_to_MiniEvent::iterator MiniEvent_intersecting_minis_AB = my_MiniEvents.begin();
            MiniEvent_intersecting_minis_AB != my_MiniEvents.end();
            ++MiniEvent_intersecting_minis_AB)
    {
	
        const Event *const Event_intersecting_minis_AB = MiniEvent_intersecting_minis_AB->second.full_Event;   //readability


        //need to cycle through because some varpos can only be caught from one side (e.g. inserts relative to other mini)
        for ( type_set_uint::const_iterator it_Mini_A_MiniID = MiniEvent_intersecting_minis_AB->second.affected_MiniRegions.begin();
                it_Mini_A_MiniID != MiniEvent_intersecting_minis_AB->second.affected_MiniRegions.end();
                ++it_Mini_A_MiniID)
        {                                    
            
            if (     MiniEvent_intersecting_minis_AB->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(*it_Mini_A_MiniID) > 0 
                 or  MiniEvent_intersecting_minis_AB->second.MiniRegions_in_insertion_of_an_LCR.count( *it_Mini_A_MiniID )  > 0       )
	    {  continue;  }
                                    
            //else, has a partner and is potentially a hybrid read.                                                  
            
            const type_map_uint_to_MiniRegion::const_iterator  mini_A_it = my_MiniRegions.find(*it_Mini_A_MiniID);  //readability.

                                                                                    
            const type_map_uint_to_set_uint::const_iterator  it_to_set_of_MiniRegion_partners_for_MiniRegion_A
                                =  MiniEvent_intersecting_minis_AB->second.partners_for_MiniRegions_on_LCRs.find( *it_Mini_A_MiniID );
            
            for (type_set_uint::const_iterator it_MiniB_MiniID = it_to_set_of_MiniRegion_partners_for_MiniRegion_A->second.begin();
                    it_MiniB_MiniID != it_to_set_of_MiniRegion_partners_for_MiniRegion_A->second.end();
                    ++it_MiniB_MiniID)
            {
                
                const type_map_uint_to_MiniRegion::const_iterator mini_B_it =  MiniEvent_intersecting_minis_AB->second.parent_PER->my_MiniRegions.find(*it_MiniB_MiniID);                           
                
                    
                assert( MiniEvent_intersecting_minis_AB->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(mini_B_it->second.MiniID) == 0  );
                assert( MiniEvent_intersecting_minis_AB->second.MiniRegions_in_insertion_of_an_LCR.count(mini_B_it->second.MiniID) == 0  );



                const  uint mini_A_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID);
                const  uint mini_B_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_B_it->second.MiniID);   
                


                // get some needed maps and  some crucial, helpful characteristics:

                // For A:            
                const type_map_uint_to_uint map_MiniEvent_profile_varpos_only_to_mini_A(
                            convert_compressed_map_to_full_map(
                                    get_compressed_map_inverse_of_compressed_map(
                                        get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                Event_intersecting_minis_AB->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(mini_A_qs),
                                                mini_A_it->second.region_interval   )  ) ));

                const type_map_uint_to_uint map_mini_A_to_MiniEvent_profile(
                        convert_compressed_map_to_full_map(
                            get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                    Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at( mini_A_qs ),
                                    mini_A_it->second.region_interval  )));      
                        
                type_map_uint_to_uint map_MiniEvent_profile_to_mini_A( get_inverse_of_map(map_mini_A_to_MiniEvent_profile) );  // (see explanation for mini_B version)
                

                

                const bool mini_A_orientation_agrees_with_MiniEvent_profile
                                = (map_MiniEvent_profile_to_mini_A.begin()->second <  (++map_MiniEvent_profile_to_mini_A.begin())->second);            

                                
                const bool must_complement_mini_A_to_get_to_MiniEvent_profile =  (mini_A_qs == 1 &&  Event_intersecting_minis_AB->use_complement_of_subject );                                
            
                            

                // For B:
                const type_map_uint_to_uint map_MiniEvent_profile_to_mini_B(
                            convert_compressed_map_to_full_map(
                                    get_compressed_map_inverse_of_compressed_map(
                                            get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                  Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at( mini_B_qs ),
                                                  mini_B_it->second.region_interval  )) ) ); 
                //if the breakpoint is at a variational position of mini_A that does NOT appear in mini_B (e.g. the breakpoint is an insertion in A relative to B), then we need to know the full sequence of B in order to know what the "first" nucleotide of B would be in the hybrid, e.g. exactly how the mini sequences would hybridize.


                const bool mini_B_orientation_agrees_with_MiniEvent_profile
                                = (map_MiniEvent_profile_to_mini_B.begin()->second  <  (++map_MiniEvent_profile_to_mini_B.begin())->second);            

                                
                const bool must_complement_mini_B_to_get_to_MiniEvent_profile =  ( mini_B_qs == 1   &&   Event_intersecting_minis_AB->use_complement_of_subject );

                
                        
		// to understand the following to lines, think: "double negatives give a positive"  (if you have to reverse orientation to get to the spwaning EV profile AND you have to rteverse orientation to get to the ME profile, then the spawning profile and the ME profile must have the same orientation).
                const bool orientation_of_PER_and_ME_profile_agree
                    =  (mini_A_it->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree ==  mini_A_orientation_agrees_with_MiniEvent_profile);
                const bool must_complement_PER_to_get_to_ME_profile
                    =  (mini_A_it->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile != must_complement_mini_A_to_get_to_MiniEvent_profile);


                    
                //if miniRegions are too far apart wrt profile, then they will never bee seen together in the same PER even if there is a hybrid (they must come from different homologosu regions of the LCR).
                if ( !BOOST_overlap(
                                    BOOST_Interval( map_MiniEvent_profile_to_mini_A.begin()->first, (--map_MiniEvent_profile_to_mini_A.end())->first ),
                                    BOOST_Interval( map_MiniEvent_profile_to_mini_B.begin()->first, (--map_MiniEvent_profile_to_mini_B.end())->first )  ) )
                {
		    std::stringstream  error_strm;
                    print_line_of_markers("ERROR! ", &error_strm);
                    print_line_of_markers("(", &error_strm);
                    error_strm << "mini_A_it->second.MiniID = " << mini_A_it->second.MiniID << "   and   mini_B_it->second.MiniID = " << mini_B_it->second.MiniID << "   DO NOT INTERSECT!!!\n";
                    mini_A_it->second.print_this_MiniRegion(&error_strm);
                    mini_B_it->second.print_this_MiniRegion(&error_strm);
                    print_line_of_markers(")", &error_strm);
		    std::cerr << "\n\n" << error_strm.str() << "\n\n";
                }

                map_MiniEvent_profile_to_mini_A.clear();


                
                
                for (type_map_uint_to_uint::const_iterator it_ME_varpos_to_mini_A__low = map_MiniEvent_profile_varpos_only_to_mini_A.begin();
                        it_ME_varpos_to_mini_A__low != map_MiniEvent_profile_varpos_only_to_mini_A.end(); 
                        ++it_ME_varpos_to_mini_A__low)
                {                     
		    type_map_uint_to_uint::const_iterator  it_END_for_consideration;
		    if (consider_GeneConversion_breakpoints)
		    {  it_END_for_consideration = map_MiniEvent_profile_varpos_only_to_mini_A.end();  }
		    else
		    {
			it_END_for_consideration = it_ME_varpos_to_mini_A__low;  
			++it_END_for_consideration;
		    }
		    
		    
		    for (type_map_uint_to_uint::const_iterator it_ME_varpos_to_mini_A__high = it_ME_varpos_to_mini_A__low;
			    it_ME_varpos_to_mini_A__high != it_END_for_consideration; 
			    ++it_ME_varpos_to_mini_A__high)
		    {			
			const BOOST_Interval profile_brkpt_interval(it_ME_varpos_to_mini_A__low->first, it_ME_varpos_to_mini_A__high->first);

			const bool is_NAHR = BOOST_is_a_point(profile_brkpt_interval);
			
			//if already recorded, skip.
			if ( MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.count(profile_brkpt_interval) > 0 )
			{  continue;  }                   
			

			const type_map_uint_to_uint::const_iterator it_found_B__low = find_key_in_map_or_find_next_key_in_specified_direction_WRT_keys<uint,uint>
										    (map_MiniEvent_profile_to_mini_B,
										    it_ME_varpos_to_mini_A__low->first,
										    true);  //since is always wrt profile (regardless if B reversed or not)
										    
			const type_map_uint_to_uint::const_iterator it_found_B__high = find_key_in_map_or_find_next_key_in_specified_direction_WRT_keys<uint,uint>
										    (map_MiniEvent_profile_to_mini_B,
										    it_ME_varpos_to_mini_A__high->first,
										    true);  //since is always wrt profile (regardless if B reversed or not)										    
			// (we use "map_MiniEvent_profile_to_mini_B" and not "map_MiniEvent_profile_varpos_only_to_mini_B"   here because miniB may have other variational positions (i.e. varpos related to other events, that come before the varpos for THIS PARTICULAR MiniEvent).

			if ( it_found_B__low  == map_MiniEvent_profile_to_mini_B.end()   or   it_found_B__high == map_MiniEvent_profile_to_mini_B.end()  )
			{  continue;  }

			    
			    
			const BOOST_Interval absolute_brkpt_interval__on_A(  std::min<uint>(it_ME_varpos_to_mini_A__low->second, it_ME_varpos_to_mini_A__high->second),
									     std::max<uint>(it_ME_varpos_to_mini_A__low->second, it_ME_varpos_to_mini_A__high->second));	

			const BOOST_Interval absolute_brkpt_interval__on_B(  std::min<uint>(it_found_B__low->second, it_found_B__high->second),
									     std::max<uint>(it_found_B__low->second, it_found_B__high->second));
									     			
									    
			//so     it_map_A:  mini_coord  --->  varpos         it_found_B: same_varpos  --->  mini_B_coord
			//so our hybrid will be:  sequence A  from [beginning_of_mini_A,  it_map_A->first) NOT INCLUSIVE (since the breakpoint reflects the NEW sequence)
			// and then               sequence B  from [it_found_B->second, end_of_mini_B]
			//Note:  "beginning_of_mini_A" and "end_of_mini_B"  dpened on orientations (beginning may in fact mean end, and end may in fact mean beginning).


			//We orient ourselves WITH the MINIEVENT PROFILE.  Here, we want the lower part of the profile to be sequence A, and the upper part of the profile to be sequence B.  Depending on the orientations of A,B to the profile, we need to play some games here....
			const type_string__BI hybrid_AB__breakpoint(
						    get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
									    (mini_A_it->second,
									    absolute_brkpt_interval__on_A,    //not to be included!!
									    must_complement_mini_A_to_get_to_MiniEvent_profile,
									    !mini_A_orientation_agrees_with_MiniEvent_profile,   //always want beginning of A
									    mini_B_it->second,
									    absolute_brkpt_interval__on_B,
									    must_complement_mini_B_to_get_to_MiniEvent_profile,
									    !mini_B_orientation_agrees_with_MiniEvent_profile,    // take_this_side_of_B_sequence)   );    //CAREFUL!!!
									    is_NAHR));  // note that the breakpoint is relative to the length of the hybrid      
			

			//and the flip side....
			//We orient ourselves WITH the MINIEVENT PROFILE.  Here, we want the lower part of the profile to be sequence B, and the upper part of the profile to be sequence A.  Depending on the orientations of A,B to the profile, we need to play some games here....
			const type_string__BI hybrid_BA__breakpoint(
							get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
								    (mini_B_it->second,
								    absolute_brkpt_interval__on_B,  //profile to mini_b.   //not to be included!!
								    must_complement_mini_B_to_get_to_MiniEvent_profile,
								    !mini_B_orientation_agrees_with_MiniEvent_profile,    // tricky!
								    mini_A_it->second,
								    absolute_brkpt_interval__on_A,  //inclusive
								    must_complement_mini_A_to_get_to_MiniEvent_profile,
								    !mini_A_orientation_agrees_with_MiniEvent_profile,//take_this_side_of_A_sequence)   // //is this correct?
								    is_NAHR));  // note that the breakpoint is relative to the length of the hybrid			
			
			
			
			
						    
			if (   determine_observed_frag_length() > hybrid_AB__breakpoint.first.size()
			    or determine_observed_frag_length() > hybrid_BA__breakpoint.first.size())
			{                                                    
			    continue;
			}//suff length
			    
			    
			
			    
			//Now, we want to align the PER to these...		
			
			//submit work AB:
			const CUPA::type_list_Paired_read_alignment_bundle::const_iterator it_submitted_work__hybrid_AB
				    = PER_bundle_wrapper__space::format_and_submit_alignment_work_to_GPU_workload
					    (GPU_work,			     
					    mates.first, mates.second,
					    hybrid_AB__breakpoint.first,
					    orientation_of_PER_and_ME_profile_agree,
					    must_complement_PER_to_get_to_ME_profile);
						    
				
			//BA:

			//This hybrid string is now properly oriented and complemented so that it maps directly to "Event_intersecting_minis_AB" 's profile.
			//Thus, the orientation/complementation of the PERs should be based on "Event_intersecting_minis_AB" 's profile...
			//
			// or... maybe we can cleverly orient/complement hybrid_AB so that it will agree with the PERs........!!!!!!!!!!!!!!!!!!!!

		
			//BUT THE PERS ARE ORIENTATED/COMPLEMENTED TO FIT THE VARI_POS PROFILE, NOT THIS INTERSECTING EVENT'S PROFILE!!!!!!!!!!!!!!!!!!!!!!!

			//Finally, go through and calculate the proability according to the varpos and taking note of the breakpoint...
			// AND TAKE NOTE OF THE ORIENTATION OF A!!!!!!
			
			//submit work BA:
			const CUPA::type_list_Paired_read_alignment_bundle::const_iterator it_submitted_work__hybrid_BA
				    = PER_bundle_wrapper__space::format_and_submit_alignment_work_to_GPU_workload
					    (GPU_work,			     
					    mates.first, mates.second,
					    hybrid_BA__breakpoint.first,
					    orientation_of_PER_and_ME_profile_agree,
					    must_complement_PER_to_get_to_ME_profile);    
					
					    
					    
			//record wrapper:
			record_of_PER_algmt_jobs[name].push_back(
				PER_bundle_wrapper__space::PER_bundle_wrapper( 
							it_submitted_work__hybrid_AB, it_submitted_work__hybrid_BA,
							this,    
							mini_A_it, mini_B_it,
							MiniEvent_intersecting_minis_AB,
							profile_brkpt_interval,
							hybrid_AB__breakpoint.second, hybrid_BA__breakpoint.second,
							absolute_brkpt_interval__on_A, absolute_brkpt_interval__on_B,
							is_NAHR,
							mini_A_orientation_agrees_with_MiniEvent_profile,
							mini_B_orientation_agrees_with_MiniEvent_profile));
							
							
		    }//high
		}//low				
	    }//mini_B_id (partner of mini A)	    
	}//mini_A_id
    }  // end for-loop    it_ev



}  //  end of   align_each_PER_to_every_possible_hybrid_outcome_of_MiniRegions_and_record_probabilities_in_each_PER























type_uint__real Paired_end_read::calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints
                                                      (const bool &haploid_0_or_1,
                                                       const Sparse_map &sparse_haploid_state_vector,
                                                       const type_map_uint_to_BI &sparse_haploid_breakpoints)  const
{
    
    //first filter everything:
    Sparse_map filtered_sparse_haploid_state_vector;
    //type_map_uint_to_map_uint_to_int sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs__to__relevant_breakpoint_for_that_MR;
    type_map_uint_to_set_uint sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs;
    type_map_uint_to_set_uint sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region;
                                                //includes MiniRegions in deleted region (even though they won't be considered).
            
    
    
    filter_haploid_state_vector_through_MiniEvents
			    (filtered_sparse_haploid_state_vector,
			    sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs,
			    sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region,
			    sparse_haploid_state_vector);
    // some of the MR relevant brkpts may/will overlap.  check this! 



                          

//// ************************************************************************************************************************************** ////
//// ************************************************************************************************************************************** ////
//// ************************************************************************************************************************************** ////
//// ************************************************************************************************************************************** ////


    real sum_mPER__cond__Events_and_brkpts(0.00L);
    uint number_of_regions_contributing_to_sum = 0;
        //There are 4 types of contributions to this sum:
        //              1)  The occurring MiniEvents will result in hybrid LCRs.  We must sum over these hybridized MiniRegions.
        //              2)  The occurring MiniEvents may affect the MiniRegions lying  in the region Inbetween its two LCRs.  We must sum over these MiniRegions appropriately (according to the outcome.)
        //              3)  There are MiniRegions that are affected by MiniEvents that are NOT occurring in this state vector.  We sum over these MiniRegions, being careful not to consider a MiniRegion that was already considered in one of the previous steps (a MiniRegion may pertain to multiple different MiniEvents, e.g. exclusivity constraints).
        //              4)  There are MiniRegions that are not potentially affected by an MiniEvent (they live in "remote corners of the genome").  We always sum over these events, regardless of the state vector.


    
//// ************************************************************************************************************************************** ////    
//     1)  The occurring MiniEvents will result in hybrid LCRs.  We must sum over these hybridized MiniRegions.
           
    
    for (type_map_uint_to_set_uint::const_iterator it_ME_to_MR_on_LCRs = sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs.begin();
            it_ME_to_MR_on_LCRs != sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs.end();
            ++it_ME_to_MR_on_LCRs)
    {
        const type_uint__real returned_pair(              
              my_MiniEvents.at(it_ME_to_MR_on_LCRs->first).
                    calculate__haploid_sum___P_mPER__cond_this_Event_occurred_with_these_relevant_profile_breakpoints____ONLY_SUM_OVER_MINIREGIONS_ON_THE_LCRS
                                                    (static_cast<type_haploid_outcome>(filtered_sparse_haploid_state_vector.get_state_for_UID(it_ME_to_MR_on_LCRs->first)),
                                                     it_ME_to_MR_on_LCRs->second,
						     sparse_haploid_breakpoints.at(it_ME_to_MR_on_LCRs->first)));   
                                                     

        number_of_regions_contributing_to_sum += returned_pair.first;
        sum_mPER__cond__Events_and_brkpts += returned_pair.second;
    }//ME



    
//// ************************************************************************************************************************************** ////
//              2)  The occurring MiniEvents may affect the MiniRegions lying  in the region Inbetween its two LCRs.  We must sum over these MiniRegions appropriately (according to the outcome.)    
      
    
    for (type_map_uint_to_set_uint::const_iterator it_Inbetween = sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region.begin();
            it_Inbetween != sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region.end();
            ++it_Inbetween)
    {//ME
	switch (sparse_haploid_state_vector.get_state_for_UID(it_Inbetween->first))
        //switch (my_MiniEvents.at(it_Inbetween->first).full_Event->recomb_type)
        {
	    case hap_outcome__Del:		
		break;   //do nothing since these regions got deleted. 
		
            case hap_outcome__Dup:
            {
		for (type_set_uint::const_iterator it_inbetween_MR = it_Inbetween->second.begin();
			it_inbetween_MR != it_Inbetween->second.end();
			++it_inbetween_MR)
		{
		    sum_mPER__cond__Events_and_brkpts += (map_affected_MiniID__to__P_nohybrid.at(*it_inbetween_MR) * 2);  //each copy exists twice
		    number_of_regions_contributing_to_sum += 2; //each copy exists twice
		}                   
		break;
            }//dup                
            
            case hap_outcome__Inv:
            case hap_outcome__GeneConv_BAB:
	    case hap_outcome__GeneConv_ABA:
            {
                for (type_set_uint::const_iterator it_inbetween_MR = it_Inbetween->second.begin();
                        it_inbetween_MR != it_Inbetween->second.end();
                        ++it_inbetween_MR)
                {
                    sum_mPER__cond__Events_and_brkpts +=  map_affected_MiniID__to__P_nohybrid.at(*it_inbetween_MR);  
                    ++number_of_regions_contributing_to_sum;
                }                            
		break;
            }//inv, GeneConv 

                        
            default:
            {                
                print_line_of_markers("ERROR! ");
                print_line_of_markers("(");
                std::fprintf(stderr, "\n\nERROR!  unrecognized recomb type = %u   for MiniEvent  %u  of PER   %s\n\n",
                             my_MiniEvents.at(it_Inbetween->first).full_Event->recomb_type,
                             it_Inbetween->first,
                             name.c_str()  );
                print_map_to_set<uint, uint>(sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region,
                                             "sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region");             
                             
                std::fprintf(stderr, "\n\n\n\nMiniEvent generating this error:\n\n"); 
                my_MiniEvents.at(it_Inbetween->first).print_this_MiniEvent();
                std::fprintf(stderr, "\n\n\n\nPER generating the error:\n\n");
                print_this_PER();
                print_line_of_markers(")");
                exit(1);
		break;    
            }            
        };//switch
    }//ME


//// ************************************************************************************************************************************** ////
//              3)  There are MiniRegions that are potentially affected by MiniEvents, but are NOT occurring in this state vector.  We sum over these MiniRegions, being careful not to consider a MiniRegion that was already considered in one of the previous steps (a MiniRegion may pertain to multiple different MiniEvents due to exclusivity constraints).  
    
    
    type_map_uint_to_real MiniRegions_pertaining_to_NON_occurring_MiniEvents_that_have_not_yet_been_considered( map_affected_MiniID__to__P_nohybrid );
    
    //erase according to occurring MiniEvents (LCRs)
    for (type_map_uint_to_set_uint::const_iterator it_ME_to_MR_on_LCRs = sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs.begin();
            it_ME_to_MR_on_LCRs != sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs.end();
            ++it_ME_to_MR_on_LCRs)
    {
        for (type_set_uint::const_iterator it_MR_on_LCRs = it_ME_to_MR_on_LCRs->second.begin();
                it_MR_on_LCRs != it_ME_to_MR_on_LCRs->second.end();
                ++it_MR_on_LCRs)
	{
            MiniRegions_pertaining_to_NON_occurring_MiniEvents_that_have_not_yet_been_considered.erase(*it_MR_on_LCRs); 
	}//MR
    }//ME
    
    
        
        
    //erase according to gender
    if (gender_of_individual_is_female  or   !haploid_0_or_1)  //  female or male and X  --->  erase Y
    {
        for (type_set_uint::const_iterator it_MR_Y = set_of_MiniRegions_on__YYYY_chromosome.begin();
                it_MR_Y != set_of_MiniRegions_on__YYYY_chromosome.end();
                ++it_MR_Y)
	{  MiniRegions_pertaining_to_NON_occurring_MiniEvents_that_have_not_yet_been_considered.erase(*it_MR_Y);  }
    }
    else  //male and Y   --  erase X
    {
        for (type_set_uint::const_iterator it_MR_X = set_of_MiniRegions_on__XXXX_chromosome.begin();
                it_MR_X != set_of_MiniRegions_on__XXXX_chromosome.end();
                ++it_MR_X)
	{  MiniRegions_pertaining_to_NON_occurring_MiniEvents_that_have_not_yet_been_considered.erase(*it_MR_X);  }
    }
    
        
    //erase according to occurring MiniEvents (Inbetween)
    for (type_map_uint_to_set_uint::const_iterator it_Inbetween = sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region.begin();
            it_Inbetween != sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region.end();
            ++it_Inbetween)
    {
        for (type_set_uint::const_iterator it_inbetween_MRs = it_Inbetween->second.begin();
                it_inbetween_MRs != it_Inbetween->second.end();
                ++it_inbetween_MRs)        
	{  MiniRegions_pertaining_to_NON_occurring_MiniEvents_that_have_not_yet_been_considered.erase(*it_inbetween_MRs);  }
    }

        
    //finally, add the remaining MRs.
    for (type_map_uint_to_real::const_iterator it_remaining = MiniRegions_pertaining_to_NON_occurring_MiniEvents_that_have_not_yet_been_considered.begin();
            it_remaining != MiniRegions_pertaining_to_NON_occurring_MiniEvents_that_have_not_yet_been_considered.end();
            ++it_remaining) 
    {
        sum_mPER__cond__Events_and_brkpts += it_remaining->second; 
        ++number_of_regions_contributing_to_sum;
    }
//         if (it_remaining->second < 1.0000L)       //some MR may be unique places, i.e. have no variational positions
            

         
    
//// ************************************************************************************************************************************** ////
//              4)  There are MiniRegions that are not potentially affected by any MiniEvent (they live in "remote corners of the genome").  We always sum over these events, regardless of the state vector.    
    
    sum_mPER__cond__Events_and_brkpts += haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal.second;    
    number_of_regions_contributing_to_sum += haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal.first;
    
    if (gender_of_individual_is_female  or  !haploid_0_or_1)  //  female or male and X  --->  get X
    {
        sum_mPER__cond__Events_and_brkpts += haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome.second;    
        number_of_regions_contributing_to_sum += haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome.first;    
    }
    else // get Y
    {
        sum_mPER__cond__Events_and_brkpts += haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome.second;    
        number_of_regions_contributing_to_sum += haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome.first;              
    }
    
//     std::cerr<<  "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tafter never potentially affected MiniRegions:\t\t\t"  << sum_mPER__cond__Events_and_brkpts <<  "\n\n\n";     
    
    
//// ************************************************************************************************************************************** ////
//// ************************************************************************************************************************************** ////
//// ************************************************************************************************************************************** ////
//// ************************************************************************************************************************************** ////



    const type_uint__real number_of_regions__and__sum_mPER__cond__Events_and_brkpts(number_of_regions_contributing_to_sum, sum_mPER__cond__Events_and_brkpts);

    
    

    return number_of_regions__and__sum_mPER__cond__Events_and_brkpts;

} // end of  calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints





















  



void Paired_end_read::filter_haploid_state_vector_through_MiniEvents
                        (Sparse_map& filtered_sparse_haploid_state_vector,
                        type_map_uint_to_set_uint &sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs,
                        type_map_uint_to_set_uint sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region,
                        const Sparse_map &proposed_sparse_haploid_state_vector)    const
{
    
    type_map_uint_to_uint::iterator insert_filtered_states_it = filtered_sparse_haploid_state_vector.sparse_map_UID_to_value.begin();
    type_map_uint_to_set_uint::iterator insert_ME_to_MR_it = sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs.begin();
    type_map_uint_to_set_uint::iterator insert_Inbetween_it = sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region.begin();     
    
                              
    for (type_map_uint_to_uint::const_iterator proposed_states_it = proposed_sparse_haploid_state_vector.sparse_map_UID_to_value.begin();
	    proposed_states_it != proposed_sparse_haploid_state_vector.sparse_map_UID_to_value.end();
	    ++proposed_states_it)
    {        
	const type_map_uint_to_MiniEvent::const_iterator it_check_ME = my_MiniEvents.find(proposed_states_it->first);  
	
	if (it_check_ME != my_MiniEvents.end())
	{	
	    //record haploid event and outcome:            	
	    insert_filtered_states_it = filtered_sparse_haploid_state_vector.sparse_map_UID_to_value.insert(insert_filtered_states_it,  *proposed_states_it);                                                                                                                                                                                   
	    
	    insert_Inbetween_it = sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region
				    .insert(insert_Inbetween_it, type_uint__set_uint(it_check_ME->first, it_check_ME->second.MiniRegions_in_Inbetween_of_this_MiniEvent)); 
						
				
	    const type_set_uint MRs_on_LCRs_of_this_ME(
		    extract_keys_of_map_and_return_as_set<uint, type_set_uint>(
			    it_check_ME->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion));	    
						
	    insert_ME_to_MR_it = sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs
				    .insert(insert_ME_to_MR_it, type_uint__set_uint(it_check_ME->first, MRs_on_LCRs_of_this_ME));
	}//found ME
    }//proposed

}  // end of    filter_haploid_state_vector_and_breakpoints_through_PERs_MiniEvents



































type_map_uint_to_MiniRegion  Paired_end_read::get_MRs()  const
{
    return my_MiniRegions;
}








void Paired_end_read::add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately
                                        ( const Paired_end_read &another_PER )
{
    
    add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately(another_PER.get_MRs());
}                                        



void Paired_end_read::add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately
                                                        (const type_map_uint_to_MiniRegion &some_MRs)
                                                        //( const type_list_MiniRegion &homologous_MiniRegions )
{

    type_map_uint_to_MiniRegion::iterator insert_MR  =  my_MiniRegions.begin();// MiniID_ctr  ==  0    ?
//                                                               my_MiniRegions.begin()
//                                                               :
//                                                               --my_MiniRegions.end();
                                                              

    for (type_map_uint_to_MiniRegion::const_iterator it_new_MR = some_MRs.begin();
            it_new_MR != some_MRs.end();
            ++it_new_MR)
    {

        bool already_has_similar_interval = false;
        for (type_map_uint_to_MiniRegion::const_iterator it_old_MR = my_MiniRegions.begin();
                it_old_MR != my_MiniRegions.end();
                ++it_old_MR)        
            if (   it_old_MR->second.chromosome_of_region == it_new_MR->second.chromosome_of_region
                 and  BOOST_overlap( it_old_MR->second.region_interval, it_new_MR->second.region_interval )   )
            {
                already_has_similar_interval = true;
                
//                 //error-check that some of the properites are the same:
//                 if (  (    it_old_MR->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree
//                         != it_new_MR->orientation_of_MiniRegion_and_spawning_Event_profile_agree )
//                      ||                    
//                       ( it_old_MR->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile 
//                         != it_new_MR->MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile )
//                    )
//                 {
//                     print_line_of_markers("ERROR! ");
//                     print_line_of_markers("(");                    
//                     std::fprintf(stderr, "Two intersecting MiniRegions were found for this same PER, but they have different properties!\n");
//                     print_this_PER();
//                     std::fprintf(stderr, "\"Old\" MiniRegion:\n");
//                     it_old_MR->second.print_this_MiniRegion();
//                     std::fprintf(stderr, "\"New\" MiniRegion:\n");
//                     it_new_MR->print_this_MiniRegion();         
//                     print_line_of_markers(")");
//                 }
                break;
            }        


        if ( !already_has_similar_interval )  
        {                        
            insert_MR = my_MiniRegions.insert( insert_MR,   type_uint__MiniRegion(MiniID_ctr, it_new_MR->second) );
            insert_MR->second.MiniID = MiniID_ctr;
            ++MiniID_ctr;
            
//             if (BOOST_width_inclusive(insert_MR->second.region_interval) <=  my_Readgroup_ptr->median_insert_size)
//             {
//                 std::stringstream error_strm;
//                 print_line_of_markers("ERROR! ", &error_strm);
//                 print_line_of_markers("(", &error_strm);
//                 error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
//                 
//                 error_strm << "\n\nERROR!   in \"add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately\"\n\n"  
//                             << "BOOST_width_inclusive(insert_MR->second.region_interval) = " << BOOST_width_inclusive(insert_MR->second.region_interval)
//                             << "     <    " << my_Readgroup_ptr->median_insert_size << "  =  my_Readgroup_ptr->median_insert_size \n\n";
//                 
//                 print_line_of_markers(")", &error_strm);
//                 std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );            
//             }
        }
        
//         std::fprintf(stderr, "%s.add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions - MiniID_ctr = %u\n",
//                              name.c_str(), MiniID_ctr-1);

    } // end for-loop    homologous_MiniRegions


} // end of    add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately





void Paired_end_read::add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately
                                                        (const MiniRegion &some_MR)
                                                        //( const type_list_MiniRegion &homologous_MiniRegions )
{

    type_map_uint_to_MiniRegion::iterator insert_MR  =    my_MiniRegions.begin();// MiniID_ctr  ==  0    ?
//                                                               my_MiniRegions.begin()
//                                                               :
//                                                               --my_MiniRegions.end();

//     for (type_list_MiniRegion::const_iterator it_new_MR = homologous_MiniRegions.begin();
//             it_new_MR != homologous_MiniRegions.end();
//             ++it_new_MR)
//     {

        bool already_has_similar_interval = false;
        for (type_map_uint_to_MiniRegion::const_iterator it_old_MR = my_MiniRegions.begin();
                it_old_MR != my_MiniRegions.end();
                ++it_old_MR)  
	{
            if (   it_old_MR->second.chromosome_of_region == some_MR.chromosome_of_region
                 and  BOOST_overlap( it_old_MR->second.region_interval, some_MR.region_interval )   )
            {
                already_has_similar_interval = true;
                
//                 //error-check that some of the properites are the same:
//                 if (  (    it_old_MR->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree
//                         != it_new_MR->orientation_of_MiniRegion_and_spawning_Event_profile_agree )
//                      ||                    
//                       ( it_old_MR->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile 
//                         != it_new_MR->MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile )
//                    )
//                 {
//                     print_line_of_markers("ERROR! ");
//                     print_line_of_markers("(");                    
//                     std::fprintf(stderr, "Two intersecting MiniRegions were found for this same PER, but they have different properties!\n");
//                     print_this_PER();
//                     std::fprintf(stderr, "\"Old\" MiniRegion:\n");
//                     it_old_MR->second.print_this_MiniRegion();
//                     std::fprintf(stderr, "\"New\" MiniRegion:\n");
//                     it_new_MR->print_this_MiniRegion();         
//                     print_line_of_markers(")");
//                 }
                break;
            }     
	}


        if ( !already_has_similar_interval )  
        {            
            
            insert_MR = my_MiniRegions.insert( insert_MR,   type_uint__MiniRegion(MiniID_ctr, some_MR) );
            insert_MR->second.MiniID = MiniID_ctr;
            ++MiniID_ctr;
            
//             if (BOOST_width_inclusive(insert_MR->second.region_interval) <=  my_Readgroup_ptr->median_insert_size)
//             {
//                 std::stringstream error_strm;
//                 print_line_of_markers("ERROR! ", &error_strm);
//                 print_line_of_markers("(", &error_strm);
//                 error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
//                 
//                 error_strm << "\n\nERROR!   in \"add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately\"\n\n"  
//                             << "BOOST_width_inclusive(insert_MR->second.region_interval) = " << BOOST_width_inclusive(insert_MR->second.region_interval)
//                             << "     <    " << my_Readgroup_ptr->median_insert_size << "  =  my_Readgroup_ptr->median_insert_size \n\n";
//                 
//                 print_line_of_markers(")", &error_strm);
//                 std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );            
//             }
            
        }
        
//         std::fprintf(stderr, "%s.add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions - MiniID_ctr = %u\n",
//                              name.c_str(), MiniID_ctr-1);

//     } // end for-loop    homologous_MiniRegions


} // end of    add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately































void Paired_end_read::clear_all_MiniRegion_sequences__for_memory_considerations()
{
    
    for (type_map_uint_to_MiniRegion::iterator it_MR = my_MiniRegions.begin();
            it_MR != my_MiniRegions.end();
            ++it_MR)
    {  it_MR->second.region_sequence.clear();  }

}  // clear_all_MiniRegion_sequences__for_memory_considerations
































void  Paired_end_read::HEURISTIC__determine_which_breakpoints_are_strongly_supported_by_this_PER_for_a_given_Event____using_pseudo_posterior
                                  (const uint &UID_of_MiniEvent,
                                   const longdouble &posterior_threshold,
				   type_map_uint_to_uint &possibly_supported_breakpoints__01,
				   type_map_uint_to_uint &possibly_supported_breakpoints__10)  const
{            
        
    const type_map_uint_to_MiniEvent::const_iterator find_ME_it = my_MiniEvents.find(UID_of_MiniEvent);    
    if ( find_ME_it == my_MiniEvents.end()  or   find_ME_it->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.empty())
    {  return;  }
    
    if (gender_of_individual_is_female  and  find_ME_it->second.full_Event->chromos[0] == 24)//Y-chromo for women - impossible.
    {  return;  }
    
                  
    
    real all_possible_nohybrids__sum(0.00L);
    
    for (type_map_uint_to_real::const_iterator it_nohybrid = map_affected_MiniID__to__P_nohybrid.begin();
            it_nohybrid != map_affected_MiniID__to__P_nohybrid.end();
            ++it_nohybrid)    
    {
        //gender sensitive
        if (!gender_of_individual_is_female  or   my_MiniRegions.at(it_nohybrid->first).chromosome_of_region != 24)
	{  all_possible_nohybrids__sum += it_nohybrid->second;  }        
    }
    

    
    for (type_set_BI::const_iterator it_relevant_varpos =  find_ME_it->second.relevant_breakpoints_for_which_this_PER_actually_shows_a_hybrid_pattern.begin();
	    it_relevant_varpos != find_ME_it->second.relevant_breakpoints_for_which_this_PER_actually_shows_a_hybrid_pattern.end();
	    ++it_relevant_varpos)
    {
	if (!BOOST_is_a_point(*it_relevant_varpos))
	{  continue;  }
	
	const type_map_BI_to__real__real::const_iterator it_brkpt = find_ME_it->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.find(*it_relevant_varpos);
	
	if ((it_brkpt->second.first/(all_possible_nohybrids__sum + it_brkpt->second.first))  >=  posterior_threshold)
	{
	    ++possibly_supported_breakpoints__01[it_brkpt->first.lower()];
	}
	
	if ((it_brkpt->second.second/(all_possible_nohybrids__sum + it_brkpt->second.second))  >=  posterior_threshold)
	{
	    ++possibly_supported_breakpoints__10[it_brkpt->first.lower()];
	}
    }
    
    
}  //  HEURISTIC__determine_which_breakpoints_are_strongly_supported_by_this_PER_for_a_given_Event____using_posterior



































































bool Paired_end_read::align_to_every_possible_hybrid_outcome_AND_DISPLAY
                            (const uint &EV_UID,
                             const BOOST_Interval &desired_breakpoint,
			     const type_haploid_outcome event_outcome,
                             Star_alignment  &star_alignment_for_LCR_0,
                             Star_alignment  &star_alignment_for_LCR_1,
                             Star_alignment  &star_alignment_for_hybrid_AB,
                             Star_alignment  &star_alignment_for_hybrid_BA)
{    
    
    type_map_uint_to_MiniEvent::iterator MiniEvent_intersecting_minis_AB = my_MiniEvents.find(EV_UID);
    
    if (MiniEvent_intersecting_minis_AB == my_MiniEvents.end())
    {  return true;  }
    
    
    const bool desired_outcome_is_NAHR = test_if_haploid_outcome_is_NAHR(event_outcome);
    
    const BOOST_Interval padded_brkpts_of_interest(safely_subtract(desired_breakpoint.lower(), 15),
                                                     desired_breakpoint.upper() + 15);        
    
    
    const Event *const Event_intersecting_minis_AB = MiniEvent_intersecting_minis_AB->second.full_Event;   //readability        
            
    type_list_real  P_every_non_hybrid_location;

    
                                                                      
                                                    
    //best:
    
    
    //overall:
    real best_overall_nohybrid_probability(0.00L); 
    
    const MiniRegion *best_overall_nohybrid_MR = NULL;
    type_string__string best_PER_and_overall_no_hybrid_algmt;
    type_map_uint_to_uint  best_overall_nohybrid___MR_coords_absolute_to_algmt_profile;
    
    
    
    //LCR 0:    
    real best_LCR_0_nohybrid_probability(0.00L);
    
    const MiniRegion *best_LCR_0_nohybrid_MR = NULL;
    type_string__string best_PER_and_LCR_0_no_hybrid_algmt;
    type_map_uint_to_uint  best_LCR_0_nohybrid___MR_coords_absolute_to_algmt_profile;    
    
    
    
    //LCR 1:    
    real best_LCR_1_nohybrid_probability(0.00L);  
    
    const MiniRegion *best_LCR_1_nohybrid_MR = NULL;
    type_string__string best_PER_and_LCR_1_no_hybrid_algmt;
    type_map_uint_to_uint  best_LCR_1_nohybrid___MR_coords_absolute_to_algmt_profile;     
    
    
    
    
    
    
    type_map_uint_to_real::iterator it_nohybrid = map_affected_MiniID__to__P_nohybrid.begin();
    //first we record no-hybrid for ALL MiniRegions.  Next we will restrict this list to only affected MiniRegions.
    
    for (type_map_uint_to_MiniRegion::const_iterator it_mini_A = my_MiniRegions.begin();
            it_mini_A != my_MiniRegions.end();
            ++it_mini_A)
    {

        BOOST_Interval mini_A_algmt_endpts;
        type_string__string PER_and_mini_A_algmt;
                                    
        const real P_nohybrid_A(
                        my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                                                (mates.first,
                                                mates.second,
                                                it_mini_A->second.region_sequence,
                                                it_mini_A->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree,
                                                it_mini_A->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile,
                                                PER_and_mini_A_algmt,
                                                mini_A_algmt_endpts));      	                      
	
        // record:
        it_nohybrid = map_affected_MiniID__to__P_nohybrid.insert(it_nohybrid, type_uint__real(it_mini_A->first, P_nohybrid_A));
								  
	map_MiniID__to__nohybrid_almgt_interval_absolute.insert(
				    type_uint__BI(it_mini_A->first, mini_A_algmt_endpts + it_mini_A->second.region_interval.lower()));                                                                       
	
	
	
	
	P_every_non_hybrid_location.push_back(P_nohybrid_A);
   
        const type_map_uint_to_uint  map_mini_A_indeces_to_algmt_profile(get_relative_hybrid_coords_to_algmt_profile(PER_and_mini_A_algmt,mini_A_algmt_endpts));
        
        if ((--map_mini_A_indeces_to_algmt_profile.end())->first  >=   it_mini_A->second.region_sequence.size())
        {
            std::stringstream error_strm;
            error_strm << "\n\nERROR!   in no-hybrid displays:\n\n\n"
                        << "(--map_mini_A_indeces_to_algmt_profile.end())->first  =  "  <<  (--map_mini_A_indeces_to_algmt_profile.end())->first
                        << "   >=     " << it_mini_A->second.region_sequence.size() <<" =  it_mini_A->second.region_sequence.size()\n\n"
			<< "\nmini_A_algmt_endpts = " << mini_A_algmt_endpts
                        << "\n\nPER_and_mini_A_algmt:\nPER:\n" << PER_and_mini_A_algmt.first
                        << "\nmini:\n" << PER_and_mini_A_algmt.second
                        << "\n\n";                        
            it_mini_A->second.print_this_MiniRegion(&error_strm);            
            print_this_PER(&error_strm);            
            print_map_keys_and_values<uint,uint>(map_mini_A_indeces_to_algmt_profile, "map_mini_A_indeces_to_algmt_profile",  &error_strm,  true );                       
            
            error_message(error_strm, false);             
        }
        
        
        
        type_map_uint_to_uint  dummy_empty_map;
        type_map_uint_to_uint  map_mini_A_absolute_to_algmt_profile;
                           adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
                                        (map_mini_A_indeces_to_algmt_profile,
                                            BOOST_make_point(BOOST_width_inclusive(it_mini_A->second.region_interval)),  //i.e.  beyond the end
                                            it_mini_A->second.region_interval.lower(),
                                            false,  //A_mini_was_reversed
                                            0,  //absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid.  irrelevant
                                            false,   //B_mini_was_reversed.  irrelevant
					   true,  // is_NAHR
                                            map_mini_A_absolute_to_algmt_profile,   
                                            dummy_empty_map);           //irrelevant

					
                                        
        assert( !map_mini_A_absolute_to_algmt_profile.empty() );                                
        
                            
        
        if (P_nohybrid_A  > best_overall_nohybrid_probability)
        {
            best_overall_nohybrid_probability = P_nohybrid_A;
                
            best_overall_nohybrid_MR = &(it_mini_A->second);
            best_PER_and_overall_no_hybrid_algmt = PER_and_mini_A_algmt;
            
            best_overall_nohybrid___MR_coords_absolute_to_algmt_profile = map_mini_A_absolute_to_algmt_profile;             
        }
        
        



        int this_MR_intersects_one_of_the_Events_LCRs = -1;
                {//relevance                                                                   
                    for (uint qs =0 ; qs < 2 ; ++qs)     
		    {
                        if (it_mini_A->second.chromosome_of_region == Event_intersecting_minis_AB->chromos[qs])
                        {
                            const type_map_uint_to_uint map_brkpts_of_interest_to_Mini(
                                    convert_compressed_map_to_full_map(
                                            get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                    get_compressed_map_inverse_of_compressed_map(
                                                            get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                                    Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at(qs),
                                                                    it_mini_A->second.region_interval)),
                                                    padded_brkpts_of_interest)));                     
                                            
                            if ( !map_brkpts_of_interest_to_Mini.empty() )
                            {
                                this_MR_intersects_one_of_the_Events_LCRs = qs;
                                break;
                            }
                        }
		    }//qs
                }//relevance
                
        if (this_MR_intersects_one_of_the_Events_LCRs == 0)        
        {
            best_LCR_0_nohybrid_probability = P_nohybrid_A; 
                
            best_LCR_0_nohybrid_MR = &(it_mini_A->second);
            best_PER_and_LCR_0_no_hybrid_algmt = PER_and_mini_A_algmt;
            
            best_LCR_0_nohybrid___MR_coords_absolute_to_algmt_profile = map_mini_A_absolute_to_algmt_profile;                           
        }
        else if (this_MR_intersects_one_of_the_Events_LCRs == 1)
        {
            best_LCR_1_nohybrid_probability = P_nohybrid_A;
                
            best_LCR_1_nohybrid_MR = &(it_mini_A->second);
            best_PER_and_LCR_1_no_hybrid_algmt = PER_and_mini_A_algmt;
            
            best_LCR_1_nohybrid___MR_coords_absolute_to_algmt_profile = map_mini_A_absolute_to_algmt_profile;                           
        }
                                                                                 

    } // end of no-hybrid





        
        
    

//     if (   best_overall_nohybrid_MR != best_LCR_0_nohybrid_MR    and   best_overall_nohybrid_MR !=  best_LCR_1_nohybrid_MR )
//     {  return true;  }
    
    
    
    
    
    
    
    
    
    

    real best_hybrid_probability(0.00L);
    uint best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt;
    int  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = -1;    
    type_string__string best_PER_and_hybrid_algmt;
    type_map_uint_to_uint  best_hybrid___absolute_A_coords_to_algmt_profile;
    type_map_uint_to_uint  best_hybrid___absolute_B_coords_to_algmt_profile;
	
    
    
    
    
    
    
    
    
    
    
    



    //need to cycle through because some varpos can only be caught from one side (e.g. inserts relative to other mini)
    for ( type_set_uint::const_iterator it_Mini_A_MiniID = MiniEvent_intersecting_minis_AB->second.affected_MiniRegions.begin();
	    it_Mini_A_MiniID != MiniEvent_intersecting_minis_AB->second.affected_MiniRegions.end();
	    ++it_Mini_A_MiniID)
    {                        
	
	
	if (    MiniEvent_intersecting_minis_AB->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(*it_Mini_A_MiniID) > 0 
	    or  MiniEvent_intersecting_minis_AB->second.MiniRegions_in_insertion_of_an_LCR.count(*it_Mini_A_MiniID)  > 0)
	{  continue;   }                                                                
				
	//else, has a partner and is potentially a hybrid read.                                      
	
	
	
	const type_map_uint_to_MiniRegion::const_iterator  mini_A_it = my_MiniRegions.find(*it_Mini_A_MiniID);  //readability.

										
	
	const type_map_uint_to_set_uint::const_iterator  it_to_set_of_MiniRegion_partners_for_MiniRegion_A
			    =  MiniEvent_intersecting_minis_AB->second.partners_for_MiniRegions_on_LCRs.find( *it_Mini_A_MiniID );
	
	for (type_set_uint::const_iterator it_MiniB_MiniID = it_to_set_of_MiniRegion_partners_for_MiniRegion_A->second.begin();
		it_MiniB_MiniID != it_to_set_of_MiniRegion_partners_for_MiniRegion_A->second.end();
		++it_MiniB_MiniID)
	{
			    
	    const type_map_uint_to_MiniRegion::const_iterator mini_B_it =  MiniEvent_intersecting_minis_AB->second.parent_PER->my_MiniRegions.find(*it_MiniB_MiniID);                         

	    
	    
	    assert( MiniEvent_intersecting_minis_AB->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(mini_B_it->second.MiniID) == 0  );
	    assert( MiniEvent_intersecting_minis_AB->second.MiniRegions_in_insertion_of_an_LCR.count(mini_B_it->second.MiniID) == 0  );



	    const  uint mini_A_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID);
	    const  uint mini_B_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_B_it->second.MiniID);   
	    



	    // get some needed maps and  some crucial, helpful characteristics:

	    // For A:            
	    const type_map_uint_to_uint map_MiniEvent_profile_varpos_only_to_mini_A(
			convert_compressed_map_to_full_map(
				get_compressed_map_inverse_of_compressed_map(
				    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
					    Event_intersecting_minis_AB->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(mini_A_qs),
					    mini_A_it->second.region_interval   )  ) ));



	    const type_map_uint_to_uint map_mini_A_to_MiniEvent_profile(
		    convert_compressed_map_to_full_map(
			get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
				Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at( mini_A_qs ),
				mini_A_it->second.region_interval  )));    
		    
	    type_map_uint_to_uint map_MiniEvent_profile_to_mini_A( get_inverse_of_map(map_mini_A_to_MiniEvent_profile) );  // (see explanation for mini_B version)
	    

	    

	    const bool mini_A_orientation_agrees_with_MiniEvent_profile
			    = (map_MiniEvent_profile_to_mini_A.begin()->second <  (++map_MiniEvent_profile_to_mini_A.begin())->second);            

			    
	    const bool must_complement_mini_A_to_get_to_MiniEvent_profile =  ( mini_A_qs == 1 &&  Event_intersecting_minis_AB->use_complement_of_subject );                                
	
			

	    // For B:
	    const type_map_uint_to_uint map_MiniEvent_profile_to_mini_B(
			convert_compressed_map_to_full_map(
				get_compressed_map_inverse_of_compressed_map(
					get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
						Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at( mini_B_qs ),
						mini_B_it->second.region_interval  )) ) );  
	    //if the breakpoint is at a variational position of mini_A that does NOT appear in mini_B (e.g. the breakpoint is an insertion in A relative to B), then we need to know the full sequence of B in order to know what the "first" nucleotide of B would be in the hybrid, e.g. exactly how the mini sequences would hybridize.


	    const bool mini_B_orientation_agrees_with_MiniEvent_profile
			    = (  map_MiniEvent_profile_to_mini_B.begin()->second   <   (++map_MiniEvent_profile_to_mini_B.begin())->second  );            

			    
	    const bool must_complement_mini_B_to_get_to_MiniEvent_profile =  ( mini_B_qs == 1   &&   Event_intersecting_minis_AB->use_complement_of_subject );

	    

		    

	    // ??????????????????????????????????????????????????????????????????????????????????????????????????????????????
	    const bool orientation_of_PER_and_ev_profile_agree
		=  ( mini_A_it->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree ==  mini_A_orientation_agrees_with_MiniEvent_profile );
	    const bool must_complement_PER_to_get_to_ev_profile
		=  ( mini_A_it->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile != must_complement_mini_A_to_get_to_MiniEvent_profile  );


		
	    //if miniRegions are too far apart wrt profile, then they will never bee seen together in the same PER even if there is a hybrid (they must come from different homologosu regions of the LCR).
	    if ( !BOOST_overlap(
				BOOST_Interval( map_MiniEvent_profile_to_mini_A.begin()->first, (--map_MiniEvent_profile_to_mini_A.end())->first ),
				BOOST_Interval( map_MiniEvent_profile_to_mini_B.begin()->first, (--map_MiniEvent_profile_to_mini_B.end())->first )   )   )
	    {
		print_line_of_markers("ERROR! ");
		print_line_of_markers("(");
		std::fprintf(stderr, "mini_A_it->second.MiniID = %u   and   mini_B_it->second.MiniID = %u   DO NOT INTERSECT!!!\n",
			    mini_A_it->second.MiniID,
			    mini_B_it->second.MiniID );
		mini_A_it->second.print_this_MiniRegion();
		mini_B_it->second.print_this_MiniRegion();
		print_line_of_markers(")");                
	    }
	    
	    
	    
	    
	    
	    	    

	    
	    type_map_uint_to_uint::const_iterator it_ME_varpos_to_mini_A__low;
	    type_map_uint_to_uint::const_iterator it_ME_varpos_to_mini_A__high;
	    uint compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = 0;
		    
	    {//find brkpts
		it_ME_varpos_to_mini_A__low = map_MiniEvent_profile_varpos_only_to_mini_A.find(desired_breakpoint.lower());
		it_ME_varpos_to_mini_A__high = map_MiniEvent_profile_varpos_only_to_mini_A.find(desired_breakpoint.upper());
		
		const BOOST_Interval  profile_region_of_mini_A(make_BI(get_endpoints_of_keys_of_map<uint,uint>(map_MiniEvent_profile_to_mini_A)));		
		
		if (	(it_ME_varpos_to_mini_A__low == map_MiniEvent_profile_varpos_only_to_mini_A.end()  and  it_ME_varpos_to_mini_A__high == map_MiniEvent_profile_varpos_only_to_mini_A.end())
		    or  (it_ME_varpos_to_mini_A__low == map_MiniEvent_profile_varpos_only_to_mini_A.end()  and  BOOST_in(desired_breakpoint.lower(), profile_region_of_mini_A))
		    or  (it_ME_varpos_to_mini_A__high == map_MiniEvent_profile_varpos_only_to_mini_A.end()  and  BOOST_in(desired_breakpoint.upper(), profile_region_of_mini_A)))
		{  continue;  }
		else if (it_ME_varpos_to_mini_A__low == map_MiniEvent_profile_varpos_only_to_mini_A.end())
		{
		    it_ME_varpos_to_mini_A__low = it_ME_varpos_to_mini_A__high;
		    compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = 2;
		}
		else if (it_ME_varpos_to_mini_A__high == map_MiniEvent_profile_varpos_only_to_mini_A.end())
		{
		    it_ME_varpos_to_mini_A__high = it_ME_varpos_to_mini_A__low;
		    compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = 1;
		}	    	    
	    }//find brkpts		    
			


	    const BOOST_Interval profile_brkpt_interval(it_ME_varpos_to_mini_A__low->first, it_ME_varpos_to_mini_A__high->first);
	    
	    const bool compute_is_NAHR = BOOST_is_a_point(profile_brkpt_interval);		    		    
	    
	    
	    
	    //if already recorded, skip.
	    if (MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.count(profile_brkpt_interval) > 0)
	    {  continue;  }                    
	    

	    const type_map_uint_to_uint::const_iterator it_found_B__low = find_key_in_map_or_find_next_key_in_specified_direction_WRT_keys<uint,uint>
									(map_MiniEvent_profile_to_mini_B,
									it_ME_varpos_to_mini_A__low->first,
									true);  //since is always wrt profile (regardless if B reversed or not)
									
	    const type_map_uint_to_uint::const_iterator it_found_B__high = find_key_in_map_or_find_next_key_in_specified_direction_WRT_keys<uint,uint>
									(map_MiniEvent_profile_to_mini_B,
									it_ME_varpos_to_mini_A__high->first,
									true);  //since is always wrt profile (regardless if B reversed or not)										    
	    // (we use "map_MiniEvent_profile_to_mini_B" and not "map_MiniEvent_profile_varpos_only_to_mini_B"   here because miniB may have other variational positions (i.e. varpos related to other events, that come before the varpos for THIS PARTICULAR MiniEvent).

	    if (it_found_B__low  == map_MiniEvent_profile_to_mini_B.end()   or   it_found_B__high == map_MiniEvent_profile_to_mini_B.end())
	    {  continue;  }                  

		
		
	    const BOOST_Interval absolute_brkpt_interval__on_A(std::min<uint>(it_ME_varpos_to_mini_A__low->second, it_ME_varpos_to_mini_A__high->second),
								    std::max<uint>(it_ME_varpos_to_mini_A__low->second, it_ME_varpos_to_mini_A__high->second));	

	    const BOOST_Interval absolute_brkpt_interval__on_B(std::min<uint>(it_found_B__low->second, it_found_B__high->second),
								    std::max<uint>(it_found_B__low->second, it_found_B__high->second));

	    

								    
								    

						
	
	
	    //so     it_map_A:  mini_coord  --->  varpos         it_found_B: same_varpos  --->  mini_B_coord
	    //so our hybrid will be:  sequence A  from [beginning_of_mini_A,  it_map_A->first) NOT INCLUSIVE (since the breakpoint reflects the NEW sequence)
	    // and then               sequence B  from [it_found_B->second, end_of_mini_B]
	    //Note:  "beginning_of_mini_A" and "end_of_mini_B"  dpened on orientations (beginning may in fact mean end, and end may in fact mean beginning).


	    //We orient ourselves WITH the PROFILE.  Here, we want the lower part of the profile to be sequence A, and the upper part of the profile to be sequence B.  Depending on the orientations of A,B to the profile, we need to play some games here....
	    const type_string__BI hybrid_AB__breakpoint(
					get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
								(mini_A_it->second,
								absolute_brkpt_interval__on_A,    //not to be included!!
								must_complement_mini_A_to_get_to_MiniEvent_profile,
								!mini_A_orientation_agrees_with_MiniEvent_profile,   //always want beginning of A
								mini_B_it->second,
								absolute_brkpt_interval__on_B,
								must_complement_mini_B_to_get_to_MiniEvent_profile,
								!mini_B_orientation_agrees_with_MiniEvent_profile,  // take_this_side_of_B_sequence)   );    //CAREFUL!!!
								compute_is_NAHR));  // note that the breakpoint is relative to the length of the hybrid      
	    //and the flip side....
	    //We orient ourselves WITH the PROFILE.  Here, we want the lower part of the profile to be sequence B, and the upper part of the profile to be sequence A.  Depending on the orientations of A,B to the profile, we need to play some games here....
	    const type_string__BI hybrid_BA__breakpoint(
					    get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
							(mini_B_it->second,
							absolute_brkpt_interval__on_B,  //profile to mini_b.   //not to be included!!
							must_complement_mini_B_to_get_to_MiniEvent_profile,
							!mini_B_orientation_agrees_with_MiniEvent_profile,    // tricky!
							mini_A_it->second,
							absolute_brkpt_interval__on_A,  //inclusive
							must_complement_mini_A_to_get_to_MiniEvent_profile,
							!mini_A_orientation_agrees_with_MiniEvent_profile, //take_this_side_of_A_sequence)   // //is this correct?
							    compute_is_NAHR));  // note that the breakpoint is relative to the length of the hybrid			



	    
	    
	    if (   determine_observed_frag_length() > hybrid_AB__breakpoint.first.size()
		or determine_observed_frag_length() > hybrid_BA__breakpoint.first.size())
	    {                                                    
		continue;
	    }//suff length

	    //Now, we want to align the PER to these...


	    BOOST_Interval subset_of_hybrid_AB_that_is_aligned;
	    type_string__string PER_and_hybrid_AB_algmt;

			    
	    const real P_hybrid_AB(
				my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
							(mates.first,
							mates.second,
							hybrid_AB__breakpoint.first,
							orientation_of_PER_and_ev_profile_agree,
							must_complement_PER_to_get_to_ev_profile,                                             
							PER_and_hybrid_AB_algmt,
							subset_of_hybrid_AB_that_is_aligned)); 


	    const uint A_mini_offset = mini_A_orientation_agrees_with_MiniEvent_profile   ? 
					    mini_A_it->second.region_interval.lower()  :  mini_A_it->second.region_interval.upper();

    

	
	
	    //BA:

	    //This hybrid string is now properly oriented and complemented so that it maps directly to "Event_intersecting_minis_AB" 's profile.
	    //Thus, the orientation/complementation of the PERs should be based on "Event_intersecting_minis_AB" 's profile...
	    //
	    // or... maybe we can cleverly orient/complement hybrid_AB so that it will agree with the PERs........!!!!!!!!!!!!!!!!!!!!



	    //Now, we want to align the PER to these...

	    BOOST_Interval subset_of_hybrid_BA_that_is_aligned;
	    type_string__string PER_and_hybrid_BA_algmt;

	    const real P_hybrid_BA(
			    my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
						    (mates.first,
						    mates.second,
						    hybrid_BA__breakpoint.first,                                                  
						    orientation_of_PER_and_ev_profile_agree,
						    must_complement_PER_to_get_to_ev_profile,
						    PER_and_hybrid_BA_algmt,
						    subset_of_hybrid_BA_that_is_aligned));
                            
	    
					    
	    //BUT THE PERS ARE ORIENTATED/COMPLEMENTED TO FIT THE VARI_POS PROFILE, NOT THIS INTERSECTING EVENT'S PROFILE!!!!!!!!!!!!!!!!!!!!!!!

	    //Finally, go through and calculate the proability according to the varpos and taking note of the breakpoint...
	    // AND TAKE NOTE OF THE ORIENTATION OF A!!!!!!


	    const uint B_mini_offset = mini_B_orientation_agrees_with_MiniEvent_profile  ?
						mini_B_it->second.region_interval.lower()  :   mini_B_it->second.region_interval.upper();  			
		
		
		
		
		    
						    
						
		
	    
	    
	    //Here, we determine if a breakpoint is relevant or not.  If the breakpoint did not oiccur within the alignment of the PER to the hybrid sequence, then the alignment is (in theory) equivalent to aligning to the "no-hybrid" version of the appropriate seuqnece, which was done at the beginning of this function.  Thus, there is no need to save this probability again, and we shouldn't mark such a breakpoint as relevant.  Otherwise, we save it
	    //bool breakpoint_is_relevant__and_so_we_should_record_the_probability = false;
	
						    
	    const bool breakpoint_is_relevant_for_AB = BOOST_overlap(hybrid_AB__breakpoint.second, subset_of_hybrid_AB_that_is_aligned);
	    const bool breakpoint_is_relevant_for_BA = BOOST_overlap(hybrid_BA__breakpoint.second, subset_of_hybrid_BA_that_is_aligned);    
	    
	    


		    
	    if (breakpoint_is_relevant_for_AB or breakpoint_is_relevant_for_BA)
	    {                   		    
		//inside the "align_to_every..." loop, MiniRegions "A" and "B" have no relation to "LCR_0" and "LCR_1" whatsoever.
		// we can check though, if A really matches to LCR_0, etc..
		//but regardless, hybrids are always formed by orienting the hybrid to the ME profile.
		
		const bool mini_A_corresponds_to_LCR_0 = (MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID) == 0);
		
		type_map_uint_to_uint *const proper_map_abs_A_to_algmt  =    mini_A_corresponds_to_LCR_0   ?  
			    &best_hybrid___absolute_A_coords_to_algmt_profile  :  &best_hybrid___absolute_B_coords_to_algmt_profile;
			    
		type_map_uint_to_uint *const proper_map_abs_B_to_algmt  =    mini_A_corresponds_to_LCR_0   ?  
			    &best_hybrid___absolute_B_coords_to_algmt_profile  :  &best_hybrid___absolute_A_coords_to_algmt_profile; 
		
			    			    			    
		//"X" and "Y" are meant to break from "A" and "B".  So that whether we do AB/ABA or BA/BAB or whatever, we can just express it as a general "XY/XYX"
		type_map_uint_to_uint  map_relative_hybrid_coords_to_algmt_profile;
		BOOST_Interval  hybrid_breakpoint_indeces;
		type_uint__uint absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX;
		bool mini_X_was_reversed;
		uint absolute_Y_coord_at_breakpoint;
		bool mini_Y_was_reversed;
		bool algmt_should_be_considered_as_NAHR;
		
		
		
		if (P_hybrid_AB > P_hybrid_BA)       
		{
		    if (P_hybrid_AB > best_overall_nohybrid_probability)
		    {                     
			//record "best" values:
			best_hybrid_probability = P_hybrid_AB;
			best_PER_and_hybrid_algmt = PER_and_hybrid_AB_algmt;	
			best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt;
			
			//mapping stuff for star alignment:
			map_relative_hybrid_coords_to_algmt_profile = get_relative_hybrid_coords_to_algmt_profile(PER_and_hybrid_AB_algmt, subset_of_hybrid_AB_that_is_aligned);
			hybrid_breakpoint_indeces = hybrid_AB__breakpoint.second;
				    
			if (desired_outcome_is_NAHR)
			{//NAHR
			    best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 =   mini_A_corresponds_to_LCR_0  ?   2 : 3;
			
			    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first  =  mini_A_orientation_agrees_with_MiniEvent_profile   ?
					mini_A_it->second.region_interval.lower()  :   mini_A_it->second.region_interval.upper();
					
			    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first;
			    
			    absolute_Y_coord_at_breakpoint =   mini_B_orientation_agrees_with_MiniEvent_profile   ?
						    absolute_brkpt_interval__on_B.lower()   :   absolute_brkpt_interval__on_B.upper();
						    
			    mini_X_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
			    mini_Y_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
			    algmt_should_be_considered_as_NAHR = true;
			}//NAHR			    
			else
			{//Geneconv
			    if (compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 0)
			    {//full
				best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = mini_A_corresponds_to_LCR_0  ?   4  :  5;
			    
				if (mini_A_orientation_agrees_with_MiniEvent_profile)
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_A_it->second.region_interval.lower();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_brkpt_interval__on_A.upper();
				}
				else
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_A_it->second.region_interval.upper();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_brkpt_interval__on_A.lower();
				}
				
				absolute_Y_coord_at_breakpoint =   mini_B_orientation_agrees_with_MiniEvent_profile   ?
							absolute_brkpt_interval__on_B.lower()   :   absolute_brkpt_interval__on_B.upper();
							
				mini_X_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
				mini_Y_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
				algmt_should_be_considered_as_NAHR = false;
			    }//full
			    else
			    {//only one of the breakpoints
			    
				best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = mini_A_corresponds_to_LCR_0  ?  2  :  3;
				
				if (mini_A_orientation_agrees_with_MiniEvent_profile)
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_A_it->second.region_interval.lower();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = mini_A_it->second.region_interval.upper();
				}
				else
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_A_it->second.region_interval.upper();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = mini_A_it->second.region_interval.lower();
				}
				
				absolute_Y_coord_at_breakpoint =   mini_B_orientation_agrees_with_MiniEvent_profile   ?
							absolute_brkpt_interval__on_B.lower()   :   absolute_brkpt_interval__on_B.upper();
							
				mini_X_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
				mini_Y_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
				
				algmt_should_be_considered_as_NAHR = true;
			    }//only one of the breakpoints
			    
// 			    if (name.compare("SRR359106.64680923") == 0)
// 			    {
// 				std::stringstream debug_strm;
// 				debug_strm << "\n\nname = " << name << "\n"
// 					<< "inside decision Geneconv AB\n"
// 					<< "compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = " << compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt << "\n"
// 					<< "compute_is_NAHR = " << compute_is_NAHR << "\n"
// 					<< "mini_A_corresponds_to_LCR_0 = " << mini_A_corresponds_to_LCR_0 << "\n"
// 					<< "best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = " << best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt << "\n"
// 					<< "(P_hybrid_AB > P_hybrid_BA) = " << (P_hybrid_AB > P_hybrid_BA) << "\n\n";
// 				std::cerr << "\n\n\n" << debug_strm.str() << "\n\n\n";
// 			    }
			    
			    
			}//GeneConv
			
		    }//P_hybrid_AB > best_overall_nohybrid_probability
		}//P_hybrid_AB > P_hybrid_BA
		else
		{
		    if (P_hybrid_BA > best_overall_nohybrid_probability)
		    {
			//record "best" values:
			best_hybrid_probability = P_hybrid_BA;
			best_PER_and_hybrid_algmt = PER_and_hybrid_BA_algmt;
			best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt;
			
			//mapping stuff for star alignment:
			map_relative_hybrid_coords_to_algmt_profile = get_relative_hybrid_coords_to_algmt_profile(PER_and_hybrid_BA_algmt, subset_of_hybrid_BA_that_is_aligned);			
			hybrid_breakpoint_indeces = hybrid_BA__breakpoint.second;
			
			if (desired_outcome_is_NAHR)
			{//NAHR
			    best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = mini_A_corresponds_to_LCR_0  ?  3 : 2;
			
			    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first =  mini_B_orientation_agrees_with_MiniEvent_profile  ?
					mini_B_it->second.region_interval.lower()   :   mini_B_it->second.region_interval.upper();
			
			    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first;
			    
			    absolute_Y_coord_at_breakpoint =   mini_A_orientation_agrees_with_MiniEvent_profile   ?
						    absolute_brkpt_interval__on_A.lower()   :   absolute_brkpt_interval__on_A.upper();
			    
			    mini_X_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
			    mini_Y_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
			    algmt_should_be_considered_as_NAHR = true;
			}//NAHR
			else
			{//Geneconv
			    if (compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 0)
			    {//full
				best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = mini_A_corresponds_to_LCR_0  ?  5 : 4;
			    
				if (mini_B_orientation_agrees_with_MiniEvent_profile)
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_B_it->second.region_interval.lower();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_brkpt_interval__on_B.upper();
				}
				else
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_B_it->second.region_interval.upper();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_brkpt_interval__on_B.lower();
				}
				
				absolute_Y_coord_at_breakpoint =   mini_A_orientation_agrees_with_MiniEvent_profile   ?
							absolute_brkpt_interval__on_A.lower()   :   absolute_brkpt_interval__on_A.upper();
							
				mini_X_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
				mini_Y_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
				algmt_should_be_considered_as_NAHR = false;
			    }//full
			    else
			    {//only one of the breakpoints
			    
				//mini_A_corresponds_to_LCR_0 = true -->  "BAB" really is "BAB"
				best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5  = mini_A_corresponds_to_LCR_0  ?  3 : 2;
				
				//we se "BA" of "BAB"
				if (mini_B_orientation_agrees_with_MiniEvent_profile)
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_B_it->second.region_interval.lower();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = mini_B_it->second.region_interval.upper();
				}
				else
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_B_it->second.region_interval.upper();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = mini_B_it->second.region_interval.lower();
				}
				
				absolute_Y_coord_at_breakpoint =   mini_A_orientation_agrees_with_MiniEvent_profile   ?
							absolute_brkpt_interval__on_A.lower()   :   absolute_brkpt_interval__on_A.upper();
							
				mini_X_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
				mini_Y_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
				
				algmt_should_be_considered_as_NAHR = true;
			    }//only one of the breakpoints
			}//GeneConv
			
		    }//P_hybrid_BA > best_overall_nohybrid_probability
		} //P_hybrid_BA > P_hybrid_AB
		
		
		//record
		if (P_hybrid_AB  >  best_overall_nohybrid_probability  or   P_hybrid_BA  >  best_overall_nohybrid_probability)
		{
		    adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
				(map_relative_hybrid_coords_to_algmt_profile,
				hybrid_breakpoint_indeces,
				absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX, mini_X_was_reversed,
				absolute_Y_coord_at_breakpoint, mini_Y_was_reversed,
				algmt_should_be_considered_as_NAHR,
				*proper_map_abs_A_to_algmt, *proper_map_abs_B_to_algmt);
		}//record
		
		
// 			    if (name.compare("SRR359106.64680923") == 0)
// 			    {
// 				std::stringstream debug_strm;
// 				debug_strm << "\n\nname = " << name << "\n"
// 					<< "after \"adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates\"\n"
// 					<< "hybrid_breakpoint_indeces = " << hybrid_breakpoint_indeces << "\n"
// 					<< "absolute_brkpt_interval__on_B = " << absolute_brkpt_interval__on_B << "\n"
// 					<< "absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX = " << absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first 
// 					<< ",  " << absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second << "\n"
// 					<< "absolute_Y_coord_at_breakpoint = " << absolute_Y_coord_at_breakpoint << "\n"
// 					<< "mini_X_was_reversed = " << mini_X_was_reversed << "\n"
// 					<< "mini_Y_was_reversed = " << mini_Y_was_reversed << "\n"
// 					<< "algmt_should_be_considered_as_NAHR  = " << algmt_should_be_considered_as_NAHR  << "\n\n"
// 					<< "(P_hybrid_AB > P_hybrid_BA) = " << (P_hybrid_AB > P_hybrid_BA) << "\n\n";
// 					
// 					print_map_keys_and_values<uint,uint>(*proper_map_abs_A_to_algmt, "*proper_map_abs_A_to_algmt", &debug_strm, true);
// 					print_map_keys_and_values<uint,uint>(*proper_map_abs_B_to_algmt, "*proper_map_abs_B_to_algmt", &debug_strm, true);
// 					
// 					print_map_keys_and_values<uint,uint>(map_relative_hybrid_coords_to_algmt_profile, "map_relative_hybrid_coords_to_algmt_profile", &debug_strm, true);
// 				std::cerr << "\n\n\n" << debug_strm.str() << "\n\n\n";
// 			    }




		for (uint lu=0; lu<2; ++lu)
		{
		    MiniEvent_intersecting_minis_AB->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[mini_A_it->second.MiniID]
							    .insert(pair_at(profile_brkpt_interval,lu));    
							    
		    MiniEvent_intersecting_minis_AB->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[mini_B_it->second.MiniID]
							    .insert(pair_at(profile_brkpt_interval,lu));     
		}
		
		if ( MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID) == 0)
		{
		    MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
					.insert(type_BI__real__real( profile_brkpt_interval, type_real__real(P_hybrid_AB, P_hybrid_BA)));
		}
		else
		{
		    MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
					.insert(type_BI__real__real(profile_brkpt_interval, type_real__real(P_hybrid_BA, P_hybrid_AB)));    
		}
	    }  // relevant for both.

	}  //  partners for A

    }  // end for-loop     mini_A_ptr
    
    
    
    
    
    
    
    
    

    bool total_return_val = true;
    
    
    
    bool should_display__via_AB_ABA;
    bool should_display__via_BA_BAB;
    
    if (desired_outcome_is_NAHR)
    {//NAHR
	should_display__via_AB_ABA =  (best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 2  and  test_if_outcome_should_display__AB_ABA(event_outcome));
	should_display__via_BA_BAB =  (best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 3  and  test_if_outcome_should_display__BA_BAB(event_outcome));
    }//NAHR
    else
    {//Geneconv
	should_display__via_AB_ABA  =   (event_outcome == hap_outcome__GeneConv_ABA)
					and  (     (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 0  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 4)
						or (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 1  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 2)
						or (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 2  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 3));
					
	should_display__via_BA_BAB =    (event_outcome == hap_outcome__GeneConv_BAB)
					and  (     (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 0  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 5)
						or (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 1  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 3)  // "BA"
						or (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 2  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 2));
    }//GeneConv
    
//     assert( !(should_display__via_AB_ABA and should_display__via_BA_BAB));
    
// 	    if (name.compare("SRR017053.5539940") == 0)
// 	    {
// 		std::stringstream debug_strm;
// 		debug_strm << "\n\nname = " << name << "\n"
// 			<< "best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = " << best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt << "\n"
// 			<< "should_display__via_AB_ABA = " << should_display__via_AB_ABA << "\n"
// 			<< "should_display__via_BA_BAB = " << should_display__via_BA_BAB << "\n"
// 			<< "\nbest_PER_and_hybrid_algmt:\n" << best_PER_and_hybrid_algmt.first << "\n" << best_PER_and_hybrid_algmt.second << "\n\n"
// 			<< "\nbest_PER_and_hybrid_algmt:\n" << best_PER_and_LCR_0_no_hybrid_algmt.first << "\n" << best_PER_and_LCR_0_no_hybrid_algmt.second << "\n\n"
// 			<< "\nbest_PER_and_hybrid_algmt:\n" << best_PER_and_LCR_1_no_hybrid_algmt.first << "\n" << best_PER_and_LCR_1_no_hybrid_algmt.second << "\n\n";
// 		std::cerr << "\n\n\n" << debug_strm.str() << "\n\n\n";
// 	    }
    

						
    if (((best_hybrid_probability/best_overall_nohybrid_probability) > 1.005L)   // non-triviality condition
	  and   (should_display__via_AB_ABA  or  should_display__via_BA_BAB))
    {        
	
	const real sum_over_all_locations(sum_over_list<real>(P_every_non_hybrid_location));
	
	const type_map_uint_to_uint dummy_empty_map;
	bool hybrid_return_val = true;
		
	
	if (should_display__via_AB_ABA)   
	{
	    #pragma omp critical (add_display_alignment_to_Star_Alignment__only_one_thread_at_a_time_inside_and_DISPLAY_PERs____AB)
	    hybrid_return_val = star_alignment_for_hybrid_AB.add_pairwise_alignment_to_progressive_star_alignment
					    (best_hybrid_probability,
					    best_PER_and_hybrid_algmt,
					    best_hybrid___absolute_A_coords_to_algmt_profile,
					    best_hybrid___absolute_B_coords_to_algmt_profile,
					    true,
					    false,
					    mates,
					    sum_over_all_locations+best_hybrid_probability);
	}
	
	if (should_display__via_BA_BAB)
	{
	    #pragma omp critical (add_display_alignment_to_Star_Alignment__only_one_thread_at_a_time_inside_and_DISPLAY_PERs____BA)
	    hybrid_return_val = star_alignment_for_hybrid_BA.add_pairwise_alignment_to_progressive_star_alignment
					    (best_hybrid_probability,
					    best_PER_and_hybrid_algmt,
					    best_hybrid___absolute_B_coords_to_algmt_profile,
					    best_hybrid___absolute_A_coords_to_algmt_profile,                                            
					    true,
					    false,
					    mates,
					    sum_over_all_locations+best_hybrid_probability);
	}
					    
	bool LCR_0_return_val = true;
	bool LCR_1_return_val = true;
		
	if (best_LCR_0_nohybrid_MR != NULL)   
	{
	    #pragma omp critical (add_display_alignment_to_Star_Alignment__only_one_thread_at_a_time_inside_and_DISPLAY_PERs____LCR_0)
	    LCR_0_return_val = star_alignment_for_LCR_0.add_pairwise_alignment_to_progressive_star_alignment
						(best_LCR_0_nohybrid_probability,
						best_PER_and_LCR_0_no_hybrid_algmt,
						best_LCR_0_nohybrid___MR_coords_absolute_to_algmt_profile,
						dummy_empty_map,                                            
						best_LCR_0_nohybrid_MR->orientation_of_MiniRegion_and_spawning_Event_profile_agree,
						best_LCR_0_nohybrid_MR->MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile,
						mates,
						sum_over_all_locations+best_hybrid_probability); 
	}
						
	if (best_LCR_1_nohybrid_MR != NULL)              
	{
	    #pragma omp critical (add_display_alignment_to_Star_Alignment__only_one_thread_at_a_time_inside_and_DISPLAY_PERs____LCR_1)
	    LCR_1_return_val = star_alignment_for_LCR_1.add_pairwise_alignment_to_progressive_star_alignment
						(best_LCR_1_nohybrid_probability,
						best_PER_and_LCR_1_no_hybrid_algmt,
						best_LCR_1_nohybrid___MR_coords_absolute_to_algmt_profile,
						dummy_empty_map,                                            
						best_LCR_1_nohybrid_MR->orientation_of_MiniRegion_and_spawning_Event_profile_agree,
						best_LCR_1_nohybrid_MR->MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile,
						mates,
						sum_over_all_locations+best_hybrid_probability);
	}
						
						
	total_return_val = (hybrid_return_val  and   LCR_0_return_val   and    LCR_1_return_val);                                    
    }
    
        
    return total_return_val;  

}  //  end of   align_to_every_possible_hybrid_outcome_AND_DISPLAY


















































bool Paired_end_read::align_to_every_possible_hybrid_outcome_AND_DISPLAY_REGARDLESS_OF_NOHYBRIDS
                            (const uint &EV_UID,
                             const BOOST_Interval &desired_breakpoint,
			     const type_haploid_outcome event_outcome,
                             Star_alignment  &star_alignment_for_LCR_0,
                             Star_alignment  &star_alignment_for_LCR_1,
                             Star_alignment  &star_alignment_for_hybrid_AB,
                             Star_alignment  &star_alignment_for_hybrid_BA)
{    
    
    type_map_uint_to_MiniEvent::iterator MiniEvent_intersecting_minis_AB = my_MiniEvents.find(EV_UID);
    
    if (MiniEvent_intersecting_minis_AB == my_MiniEvents.end())
    {  return true;  }
    
    
    const bool desired_outcome_is_NAHR = test_if_haploid_outcome_is_NAHR(event_outcome);
    
    const BOOST_Interval padded_brkpts_of_interest(safely_subtract(desired_breakpoint.lower(), 15),
                                                     desired_breakpoint.upper() + 15);        
    
    
    const Event *const Event_intersecting_minis_AB = MiniEvent_intersecting_minis_AB->second.full_Event;   //readability        
            
    type_list_real  P_every_non_hybrid_location;

    
                                                                      
                                                    
    //best:
    
    
    //overall:
    real best_overall_nohybrid_probability(0.00L); 
    
    const MiniRegion *best_overall_nohybrid_MR = NULL;
    type_string__string best_PER_and_overall_no_hybrid_algmt;
    type_map_uint_to_uint  best_overall_nohybrid___MR_coords_absolute_to_algmt_profile;
    
    
    
    //LCR 0:    
    real best_LCR_0_nohybrid_probability(0.00L);
    
    const MiniRegion *best_LCR_0_nohybrid_MR = NULL;
    type_string__string best_PER_and_LCR_0_no_hybrid_algmt;
    type_map_uint_to_uint  best_LCR_0_nohybrid___MR_coords_absolute_to_algmt_profile;    
    
    
    
    //LCR 1:    
    real best_LCR_1_nohybrid_probability(0.00L);  
    
    const MiniRegion *best_LCR_1_nohybrid_MR = NULL;
    type_string__string best_PER_and_LCR_1_no_hybrid_algmt;
    type_map_uint_to_uint  best_LCR_1_nohybrid___MR_coords_absolute_to_algmt_profile;     
    
    
    
    
    
    
    type_map_uint_to_real::iterator it_nohybrid = map_affected_MiniID__to__P_nohybrid.begin();
    //first we record no-hybrid for ALL MiniRegions.  Next we will restrict this list to only affected MiniRegions.
    
    for (type_map_uint_to_MiniRegion::const_iterator it_mini_A = my_MiniRegions.begin();
            it_mini_A != my_MiniRegions.end();
            ++it_mini_A)
    {

        BOOST_Interval mini_A_algmt_endpts;
        type_string__string PER_and_mini_A_algmt;
                                    
        const real P_nohybrid_A(
                        my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                                                (mates.first,
                                                mates.second,
                                                it_mini_A->second.region_sequence,
                                                it_mini_A->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree,
                                                it_mini_A->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile,
                                                PER_and_mini_A_algmt,
                                                mini_A_algmt_endpts));      	                      
	
        // record:
        it_nohybrid = map_affected_MiniID__to__P_nohybrid.insert(it_nohybrid, type_uint__real(it_mini_A->first, P_nohybrid_A));
								  
	map_MiniID__to__nohybrid_almgt_interval_absolute.insert(
				    type_uint__BI(it_mini_A->first, mini_A_algmt_endpts + it_mini_A->second.region_interval.lower()));                                                                       
	
	
	
	
	P_every_non_hybrid_location.push_back(P_nohybrid_A);
   
        const type_map_uint_to_uint  map_mini_A_indeces_to_algmt_profile(get_relative_hybrid_coords_to_algmt_profile(PER_and_mini_A_algmt,mini_A_algmt_endpts));
        
        if ((--map_mini_A_indeces_to_algmt_profile.end())->first  >=   it_mini_A->second.region_sequence.size())
        {
            std::stringstream error_strm;
            error_strm << "\n\nERROR!   in no-hybrid displays:\n\n\n"
                        << "(--map_mini_A_indeces_to_algmt_profile.end())->first  =  "  <<  (--map_mini_A_indeces_to_algmt_profile.end())->first
                        << "   >=     " << it_mini_A->second.region_sequence.size() <<" =  it_mini_A->second.region_sequence.size()\n\n"
			<< "\nmini_A_algmt_endpts = " << mini_A_algmt_endpts
                        << "\n\nPER_and_mini_A_algmt:\nPER:\n" << PER_and_mini_A_algmt.first
                        << "\nmini:\n" << PER_and_mini_A_algmt.second
                        << "\n\n";                        
            it_mini_A->second.print_this_MiniRegion(&error_strm);            
            print_this_PER(&error_strm);            
            print_map_keys_and_values<uint,uint>(map_mini_A_indeces_to_algmt_profile, "map_mini_A_indeces_to_algmt_profile",  &error_strm,  true );                       
            
            error_message(error_strm, false);             
        }
        
        
        
        type_map_uint_to_uint  dummy_empty_map;
        type_map_uint_to_uint  map_mini_A_absolute_to_algmt_profile;
                           adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
                                        (map_mini_A_indeces_to_algmt_profile,
                                            BOOST_make_point(BOOST_width_inclusive(it_mini_A->second.region_interval)),  //i.e.  beyond the end
                                            it_mini_A->second.region_interval.lower(),
                                            false,  //A_mini_was_reversed
                                            0,  //absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid.  irrelevant
                                            false,   //B_mini_was_reversed.  irrelevant
					   true,  // is_NAHR
                                            map_mini_A_absolute_to_algmt_profile,   
                                            dummy_empty_map);           //irrelevant

					
                                        
        assert( !map_mini_A_absolute_to_algmt_profile.empty() );                                
        
                            
        
        if (P_nohybrid_A  > best_overall_nohybrid_probability)
        {
            best_overall_nohybrid_probability = P_nohybrid_A;
                
            best_overall_nohybrid_MR = &(it_mini_A->second);
            best_PER_and_overall_no_hybrid_algmt = PER_and_mini_A_algmt;
            
            best_overall_nohybrid___MR_coords_absolute_to_algmt_profile = map_mini_A_absolute_to_algmt_profile;             
        }
        
        



        int this_MR_intersects_one_of_the_Events_LCRs = -1;
                {//relevance                                                                   
                    for (uint qs =0 ; qs < 2 ; ++qs)     
		    {
                        if (it_mini_A->second.chromosome_of_region == Event_intersecting_minis_AB->chromos[qs])
                        {
                            const type_map_uint_to_uint map_brkpts_of_interest_to_Mini(
                                    convert_compressed_map_to_full_map(
                                            get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                    get_compressed_map_inverse_of_compressed_map(
                                                            get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                                    Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at(qs),
                                                                    it_mini_A->second.region_interval)),
                                                    padded_brkpts_of_interest)));                     
                                            
                            if ( !map_brkpts_of_interest_to_Mini.empty() )
                            {
                                this_MR_intersects_one_of_the_Events_LCRs = qs;
                                break;
                            }
                        }
		    }//qs
                }//relevance
                
        if (this_MR_intersects_one_of_the_Events_LCRs == 0)        
        {
            best_LCR_0_nohybrid_probability = P_nohybrid_A; 
                
            best_LCR_0_nohybrid_MR = &(it_mini_A->second);
            best_PER_and_LCR_0_no_hybrid_algmt = PER_and_mini_A_algmt;
            
            best_LCR_0_nohybrid___MR_coords_absolute_to_algmt_profile = map_mini_A_absolute_to_algmt_profile;                           
        }
        else if (this_MR_intersects_one_of_the_Events_LCRs == 1)
        {
            best_LCR_1_nohybrid_probability = P_nohybrid_A;
                
            best_LCR_1_nohybrid_MR = &(it_mini_A->second);
            best_PER_and_LCR_1_no_hybrid_algmt = PER_and_mini_A_algmt;
            
            best_LCR_1_nohybrid___MR_coords_absolute_to_algmt_profile = map_mini_A_absolute_to_algmt_profile;                           
        }
                                                                                 

    } // end of no-hybrid





        
        
    

//     if (   best_overall_nohybrid_MR != best_LCR_0_nohybrid_MR    and   best_overall_nohybrid_MR !=  best_LCR_1_nohybrid_MR )
//     {  return true;  }
    
    
    
    
    
    
    
    
    
    

    real best_hybrid_probability(0.00L);
    uint best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt;
    int  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = -1;    
    type_string__string best_PER_and_hybrid_algmt;
    type_map_uint_to_uint  best_hybrid___absolute_A_coords_to_algmt_profile;
    type_map_uint_to_uint  best_hybrid___absolute_B_coords_to_algmt_profile;
	
    
    
    
    
    
    
    
    
    
    
    



    //need to cycle through because some varpos can only be caught from one side (e.g. inserts relative to other mini)
    for ( type_set_uint::const_iterator it_Mini_A_MiniID = MiniEvent_intersecting_minis_AB->second.affected_MiniRegions.begin();
	    it_Mini_A_MiniID != MiniEvent_intersecting_minis_AB->second.affected_MiniRegions.end();
	    ++it_Mini_A_MiniID)
    {                        
	
	
	if (    MiniEvent_intersecting_minis_AB->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(*it_Mini_A_MiniID) > 0 
	    or  MiniEvent_intersecting_minis_AB->second.MiniRegions_in_insertion_of_an_LCR.count(*it_Mini_A_MiniID)  > 0)
	{  continue;   }                                                                
				
	//else, has a partner and is potentially a hybrid read.                                      
	
	
	
	const type_map_uint_to_MiniRegion::const_iterator  mini_A_it = my_MiniRegions.find(*it_Mini_A_MiniID);  //readability.

										
	
	const type_map_uint_to_set_uint::const_iterator  it_to_set_of_MiniRegion_partners_for_MiniRegion_A
			    =  MiniEvent_intersecting_minis_AB->second.partners_for_MiniRegions_on_LCRs.find( *it_Mini_A_MiniID );
	
	for (type_set_uint::const_iterator it_MiniB_MiniID = it_to_set_of_MiniRegion_partners_for_MiniRegion_A->second.begin();
		it_MiniB_MiniID != it_to_set_of_MiniRegion_partners_for_MiniRegion_A->second.end();
		++it_MiniB_MiniID)
	{
			    
	    const type_map_uint_to_MiniRegion::const_iterator mini_B_it =  MiniEvent_intersecting_minis_AB->second.parent_PER->my_MiniRegions.find(*it_MiniB_MiniID);                         

	    
	    
	    assert( MiniEvent_intersecting_minis_AB->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(mini_B_it->second.MiniID) == 0  );
	    assert( MiniEvent_intersecting_minis_AB->second.MiniRegions_in_insertion_of_an_LCR.count(mini_B_it->second.MiniID) == 0  );



	    const  uint mini_A_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID);
	    const  uint mini_B_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_B_it->second.MiniID);   
	    



	    // get some needed maps and  some crucial, helpful characteristics:

	    // For A:            
	    const type_map_uint_to_uint map_MiniEvent_profile_varpos_only_to_mini_A(
			convert_compressed_map_to_full_map(
				get_compressed_map_inverse_of_compressed_map(
				    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
					    Event_intersecting_minis_AB->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(mini_A_qs),
					    mini_A_it->second.region_interval   )  ) ));



	    const type_map_uint_to_uint map_mini_A_to_MiniEvent_profile(
		    convert_compressed_map_to_full_map(
			get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
				Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at( mini_A_qs ),
				mini_A_it->second.region_interval  )));    
		    
	    type_map_uint_to_uint map_MiniEvent_profile_to_mini_A( get_inverse_of_map(map_mini_A_to_MiniEvent_profile) );  // (see explanation for mini_B version)
	    

	    

	    const bool mini_A_orientation_agrees_with_MiniEvent_profile
			    = (map_MiniEvent_profile_to_mini_A.begin()->second <  (++map_MiniEvent_profile_to_mini_A.begin())->second);            

			    
	    const bool must_complement_mini_A_to_get_to_MiniEvent_profile =  ( mini_A_qs == 1 &&  Event_intersecting_minis_AB->use_complement_of_subject );                                
	
			

	    // For B:
	    const type_map_uint_to_uint map_MiniEvent_profile_to_mini_B(
			convert_compressed_map_to_full_map(
				get_compressed_map_inverse_of_compressed_map(
					get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
						Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at( mini_B_qs ),
						mini_B_it->second.region_interval  )) ) );  
	    //if the breakpoint is at a variational position of mini_A that does NOT appear in mini_B (e.g. the breakpoint is an insertion in A relative to B), then we need to know the full sequence of B in order to know what the "first" nucleotide of B would be in the hybrid, e.g. exactly how the mini sequences would hybridize.


	    const bool mini_B_orientation_agrees_with_MiniEvent_profile
			    = (  map_MiniEvent_profile_to_mini_B.begin()->second   <   (++map_MiniEvent_profile_to_mini_B.begin())->second  );            

			    
	    const bool must_complement_mini_B_to_get_to_MiniEvent_profile =  ( mini_B_qs == 1   &&   Event_intersecting_minis_AB->use_complement_of_subject );

	    

		    

	    // ??????????????????????????????????????????????????????????????????????????????????????????????????????????????
	    const bool orientation_of_PER_and_ev_profile_agree
		=  ( mini_A_it->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree ==  mini_A_orientation_agrees_with_MiniEvent_profile );
	    const bool must_complement_PER_to_get_to_ev_profile
		=  ( mini_A_it->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile != must_complement_mini_A_to_get_to_MiniEvent_profile  );


		
	    //if miniRegions are too far apart wrt profile, then they will never bee seen together in the same PER even if there is a hybrid (they must come from different homologosu regions of the LCR).
	    if ( !BOOST_overlap(
				BOOST_Interval( map_MiniEvent_profile_to_mini_A.begin()->first, (--map_MiniEvent_profile_to_mini_A.end())->first ),
				BOOST_Interval( map_MiniEvent_profile_to_mini_B.begin()->first, (--map_MiniEvent_profile_to_mini_B.end())->first )   )   )
	    {
		print_line_of_markers("ERROR! ");
		print_line_of_markers("(");
		std::fprintf(stderr, "mini_A_it->second.MiniID = %u   and   mini_B_it->second.MiniID = %u   DO NOT INTERSECT!!!\n",
			    mini_A_it->second.MiniID,
			    mini_B_it->second.MiniID );
		mini_A_it->second.print_this_MiniRegion();
		mini_B_it->second.print_this_MiniRegion();
		print_line_of_markers(")");                
	    }
	    
	    
	    
	    
	    
	    	    

	    
	    type_map_uint_to_uint::const_iterator it_ME_varpos_to_mini_A__low;
	    type_map_uint_to_uint::const_iterator it_ME_varpos_to_mini_A__high;
	    uint compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = 0;
		    
	    {//find brkpts
		it_ME_varpos_to_mini_A__low = map_MiniEvent_profile_varpos_only_to_mini_A.find(desired_breakpoint.lower());
		it_ME_varpos_to_mini_A__high = map_MiniEvent_profile_varpos_only_to_mini_A.find(desired_breakpoint.upper());
		
		const BOOST_Interval  profile_region_of_mini_A(make_BI(get_endpoints_of_keys_of_map<uint,uint>(map_MiniEvent_profile_to_mini_A)));		
		
		if (	(it_ME_varpos_to_mini_A__low == map_MiniEvent_profile_varpos_only_to_mini_A.end()  and  it_ME_varpos_to_mini_A__high == map_MiniEvent_profile_varpos_only_to_mini_A.end())
		    or  (it_ME_varpos_to_mini_A__low == map_MiniEvent_profile_varpos_only_to_mini_A.end()  and  BOOST_in(desired_breakpoint.lower(), profile_region_of_mini_A))
		    or  (it_ME_varpos_to_mini_A__high == map_MiniEvent_profile_varpos_only_to_mini_A.end()  and  BOOST_in(desired_breakpoint.upper(), profile_region_of_mini_A)))
		{  continue;  }
		else if (it_ME_varpos_to_mini_A__low == map_MiniEvent_profile_varpos_only_to_mini_A.end())
		{
		    it_ME_varpos_to_mini_A__low = it_ME_varpos_to_mini_A__high;
		    compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = 2;
		}
		else if (it_ME_varpos_to_mini_A__high == map_MiniEvent_profile_varpos_only_to_mini_A.end())
		{
		    it_ME_varpos_to_mini_A__high = it_ME_varpos_to_mini_A__low;
		    compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = 1;
		}	    	    
	    }//find brkpts		    
			


	    const BOOST_Interval profile_brkpt_interval(it_ME_varpos_to_mini_A__low->first, it_ME_varpos_to_mini_A__high->first);
	    
	    const bool compute_is_NAHR = BOOST_is_a_point(profile_brkpt_interval);		    		    
	    
	    
	    
	    //if already recorded, skip.
	    if (MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.count(profile_brkpt_interval) > 0)
	    {  continue;  }                    
	    

	    const type_map_uint_to_uint::const_iterator it_found_B__low = find_key_in_map_or_find_next_key_in_specified_direction_WRT_keys<uint,uint>
									(map_MiniEvent_profile_to_mini_B,
									it_ME_varpos_to_mini_A__low->first,
									true);  //since is always wrt profile (regardless if B reversed or not)
									
	    const type_map_uint_to_uint::const_iterator it_found_B__high = find_key_in_map_or_find_next_key_in_specified_direction_WRT_keys<uint,uint>
									(map_MiniEvent_profile_to_mini_B,
									it_ME_varpos_to_mini_A__high->first,
									true);  //since is always wrt profile (regardless if B reversed or not)										    
	    // (we use "map_MiniEvent_profile_to_mini_B" and not "map_MiniEvent_profile_varpos_only_to_mini_B"   here because miniB may have other variational positions (i.e. varpos related to other events, that come before the varpos for THIS PARTICULAR MiniEvent).

	    if (it_found_B__low  == map_MiniEvent_profile_to_mini_B.end()   or   it_found_B__high == map_MiniEvent_profile_to_mini_B.end())
	    {  continue;  }                  

		
		
	    const BOOST_Interval absolute_brkpt_interval__on_A(std::min<uint>(it_ME_varpos_to_mini_A__low->second, it_ME_varpos_to_mini_A__high->second),
								    std::max<uint>(it_ME_varpos_to_mini_A__low->second, it_ME_varpos_to_mini_A__high->second));	

	    const BOOST_Interval absolute_brkpt_interval__on_B(std::min<uint>(it_found_B__low->second, it_found_B__high->second),
								    std::max<uint>(it_found_B__low->second, it_found_B__high->second));

	    

								    
								    

						
	
	
	    //so     it_map_A:  mini_coord  --->  varpos         it_found_B: same_varpos  --->  mini_B_coord
	    //so our hybrid will be:  sequence A  from [beginning_of_mini_A,  it_map_A->first) NOT INCLUSIVE (since the breakpoint reflects the NEW sequence)
	    // and then               sequence B  from [it_found_B->second, end_of_mini_B]
	    //Note:  "beginning_of_mini_A" and "end_of_mini_B"  dpened on orientations (beginning may in fact mean end, and end may in fact mean beginning).


	    //We orient ourselves WITH the PROFILE.  Here, we want the lower part of the profile to be sequence A, and the upper part of the profile to be sequence B.  Depending on the orientations of A,B to the profile, we need to play some games here....
	    const type_string__BI hybrid_AB__breakpoint(
					get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
								(mini_A_it->second,
								absolute_brkpt_interval__on_A,    //not to be included!!
								must_complement_mini_A_to_get_to_MiniEvent_profile,
								!mini_A_orientation_agrees_with_MiniEvent_profile,   //always want beginning of A
								mini_B_it->second,
								absolute_brkpt_interval__on_B,
								must_complement_mini_B_to_get_to_MiniEvent_profile,
								!mini_B_orientation_agrees_with_MiniEvent_profile,  // take_this_side_of_B_sequence)   );    //CAREFUL!!!
								compute_is_NAHR));  // note that the breakpoint is relative to the length of the hybrid      
	    //and the flip side....
	    //We orient ourselves WITH the PROFILE.  Here, we want the lower part of the profile to be sequence B, and the upper part of the profile to be sequence A.  Depending on the orientations of A,B to the profile, we need to play some games here....
	    const type_string__BI hybrid_BA__breakpoint(
					    get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
							(mini_B_it->second,
							absolute_brkpt_interval__on_B,  //profile to mini_b.   //not to be included!!
							must_complement_mini_B_to_get_to_MiniEvent_profile,
							!mini_B_orientation_agrees_with_MiniEvent_profile,    // tricky!
							mini_A_it->second,
							absolute_brkpt_interval__on_A,  //inclusive
							must_complement_mini_A_to_get_to_MiniEvent_profile,
							!mini_A_orientation_agrees_with_MiniEvent_profile, //take_this_side_of_A_sequence)   // //is this correct?
							    compute_is_NAHR));  // note that the breakpoint is relative to the length of the hybrid			



	    
	    
	    if (   determine_observed_frag_length() > hybrid_AB__breakpoint.first.size()
		or determine_observed_frag_length() > hybrid_BA__breakpoint.first.size())
	    {                                                    
		continue;
	    }//suff length

	    //Now, we want to align the PER to these...


	    BOOST_Interval subset_of_hybrid_AB_that_is_aligned;
	    type_string__string PER_and_hybrid_AB_algmt;

			    
	    const real P_hybrid_AB(
				my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
							(mates.first,
							mates.second,
							hybrid_AB__breakpoint.first,
							orientation_of_PER_and_ev_profile_agree,
							must_complement_PER_to_get_to_ev_profile,                                             
							PER_and_hybrid_AB_algmt,
							subset_of_hybrid_AB_that_is_aligned)); 


	    const uint A_mini_offset = mini_A_orientation_agrees_with_MiniEvent_profile   ? 
					    mini_A_it->second.region_interval.lower()  :  mini_A_it->second.region_interval.upper();

    

	
	
	    //BA:

	    //This hybrid string is now properly oriented and complemented so that it maps directly to "Event_intersecting_minis_AB" 's profile.
	    //Thus, the orientation/complementation of the PERs should be based on "Event_intersecting_minis_AB" 's profile...
	    //
	    // or... maybe we can cleverly orient/complement hybrid_AB so that it will agree with the PERs........!!!!!!!!!!!!!!!!!!!!



	    //Now, we want to align the PER to these...

	    BOOST_Interval subset_of_hybrid_BA_that_is_aligned;
	    type_string__string PER_and_hybrid_BA_algmt;

	    const real P_hybrid_BA(
			    my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
						    (mates.first,
						    mates.second,
						    hybrid_BA__breakpoint.first,                                                  
						    orientation_of_PER_and_ev_profile_agree,
						    must_complement_PER_to_get_to_ev_profile,
						    PER_and_hybrid_BA_algmt,
						    subset_of_hybrid_BA_that_is_aligned));
                            
	    
					    
	    //BUT THE PERS ARE ORIENTATED/COMPLEMENTED TO FIT THE VARI_POS PROFILE, NOT THIS INTERSECTING EVENT'S PROFILE!!!!!!!!!!!!!!!!!!!!!!!

	    //Finally, go through and calculate the proability according to the varpos and taking note of the breakpoint...
	    // AND TAKE NOTE OF THE ORIENTATION OF A!!!!!!


	    const uint B_mini_offset = mini_B_orientation_agrees_with_MiniEvent_profile  ?
						mini_B_it->second.region_interval.lower()  :   mini_B_it->second.region_interval.upper();  			
		
		
		
		
		    
						    
						
		
	    
	    
	    //Here, we determine if a breakpoint is relevant or not.  If the breakpoint did not oiccur within the alignment of the PER to the hybrid sequence, then the alignment is (in theory) equivalent to aligning to the "no-hybrid" version of the appropriate seuqnece, which was done at the beginning of this function.  Thus, there is no need to save this probability again, and we shouldn't mark such a breakpoint as relevant.  Otherwise, we save it
	    //bool breakpoint_is_relevant__and_so_we_should_record_the_probability = false;
	
						    
	    const bool breakpoint_is_relevant_for_AB = BOOST_overlap(hybrid_AB__breakpoint.second, subset_of_hybrid_AB_that_is_aligned);
	    const bool breakpoint_is_relevant_for_BA = BOOST_overlap(hybrid_BA__breakpoint.second, subset_of_hybrid_BA_that_is_aligned);    
	    
	    


		    
	    if (breakpoint_is_relevant_for_AB or breakpoint_is_relevant_for_BA)
	    {                   		    
		//inside the "align_to_every..." loop, MiniRegions "A" and "B" have no relation to "LCR_0" and "LCR_1" whatsoever.
		// we can check though, if A really matches to LCR_0, etc..
		//but regardless, hybrids are always formed by orienting the hybrid to the ME profile.
		
		const bool mini_A_corresponds_to_LCR_0 = (MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID) == 0);
		
		type_map_uint_to_uint *const proper_map_abs_A_to_algmt  =    mini_A_corresponds_to_LCR_0   ?  
			    &best_hybrid___absolute_A_coords_to_algmt_profile  :  &best_hybrid___absolute_B_coords_to_algmt_profile;
			    
		type_map_uint_to_uint *const proper_map_abs_B_to_algmt  =    mini_A_corresponds_to_LCR_0   ?  
			    &best_hybrid___absolute_B_coords_to_algmt_profile  :  &best_hybrid___absolute_A_coords_to_algmt_profile; 
		
			    			    			    
		//"X" and "Y" are meant to break from "A" and "B".  So that whether we do AB/ABA or BA/BAB or whatever, we can just express it as a general "XY/XYX"
		type_map_uint_to_uint  map_relative_hybrid_coords_to_algmt_profile;
		BOOST_Interval  hybrid_breakpoint_indeces;
		type_uint__uint absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX;
		bool mini_X_was_reversed;
		uint absolute_Y_coord_at_breakpoint;
		bool mini_Y_was_reversed;
		bool algmt_should_be_considered_as_NAHR;
		
		

		if ((test_if_outcome_should_display__AB_ABA(event_outcome)  and  mini_A_corresponds_to_LCR_0)  or  (test_if_outcome_should_display__BA_BAB(event_outcome)  and !mini_A_corresponds_to_LCR_0))		
// 		if (P_hybrid_AB > P_hybrid_BA)       
		{
// 		    if (P_hybrid_AB > best_overall_nohybrid_probability)
		    {                     
			//record "best" values:
			best_hybrid_probability = P_hybrid_AB;
			best_PER_and_hybrid_algmt = PER_and_hybrid_AB_algmt;	
			best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt;
			
			//mapping stuff for star alignment:
			map_relative_hybrid_coords_to_algmt_profile = get_relative_hybrid_coords_to_algmt_profile(PER_and_hybrid_AB_algmt, subset_of_hybrid_AB_that_is_aligned);
			hybrid_breakpoint_indeces = hybrid_AB__breakpoint.second;
				    
			if (desired_outcome_is_NAHR)
			{//NAHR
			    best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 =   mini_A_corresponds_to_LCR_0  ?   2 : 3;
			
			    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first  =  mini_A_orientation_agrees_with_MiniEvent_profile   ?
					mini_A_it->second.region_interval.lower()  :   mini_A_it->second.region_interval.upper();
					
			    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first;
			    
			    absolute_Y_coord_at_breakpoint =   mini_B_orientation_agrees_with_MiniEvent_profile   ?
						    absolute_brkpt_interval__on_B.lower()   :   absolute_brkpt_interval__on_B.upper();
						    
			    mini_X_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
			    mini_Y_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
			    algmt_should_be_considered_as_NAHR = true;
			}//NAHR			    
			else
			{//Geneconv
			    if (compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 0)
			    {//full
				best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = mini_A_corresponds_to_LCR_0  ?   4  :  5;
			    
				if (mini_A_orientation_agrees_with_MiniEvent_profile)
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_A_it->second.region_interval.lower();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_brkpt_interval__on_A.upper();
				}
				else
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_A_it->second.region_interval.upper();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_brkpt_interval__on_A.lower();
				}
				
				absolute_Y_coord_at_breakpoint =   mini_B_orientation_agrees_with_MiniEvent_profile   ?
							absolute_brkpt_interval__on_B.lower()   :   absolute_brkpt_interval__on_B.upper();
							
				mini_X_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
				mini_Y_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
				algmt_should_be_considered_as_NAHR = false;
			    }//full
			    else
			    {//only one of the breakpoints
			    
				best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = mini_A_corresponds_to_LCR_0  ?  2  :  3;
				
				if (mini_A_orientation_agrees_with_MiniEvent_profile)
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_A_it->second.region_interval.lower();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = mini_A_it->second.region_interval.upper();
				}
				else
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_A_it->second.region_interval.upper();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = mini_A_it->second.region_interval.lower();
				}
				
				absolute_Y_coord_at_breakpoint =   mini_B_orientation_agrees_with_MiniEvent_profile   ?
							absolute_brkpt_interval__on_B.lower()   :   absolute_brkpt_interval__on_B.upper();
							
				mini_X_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
				mini_Y_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
				
				algmt_should_be_considered_as_NAHR = true;
			    }//only one of the breakpoints
			    
// 			    if (name.compare("SRR359106.64680923") == 0)
// 			    {
// 				std::stringstream debug_strm;
// 				debug_strm << "\n\nname = " << name << "\n"
// 					<< "inside decision Geneconv AB\n"
// 					<< "compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = " << compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt << "\n"
// 					<< "compute_is_NAHR = " << compute_is_NAHR << "\n"
// 					<< "mini_A_corresponds_to_LCR_0 = " << mini_A_corresponds_to_LCR_0 << "\n"
// 					<< "best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = " << best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt << "\n"
// 					<< "(P_hybrid_AB > P_hybrid_BA) = " << (P_hybrid_AB > P_hybrid_BA) << "\n\n";
// 				std::cerr << "\n\n\n" << debug_strm.str() << "\n\n\n";
// 			    }
			    
			    
			}//GeneConv
			
		    }//P_hybrid_AB > best_overall_nohybrid_probability
		}//P_hybrid_AB > P_hybrid_BA
// 		else
		else if ((test_if_outcome_should_display__AB_ABA(event_outcome)  and  !mini_A_corresponds_to_LCR_0)  or  (test_if_outcome_should_display__BA_BAB(event_outcome)  and mini_A_corresponds_to_LCR_0))				
		{
// 		    if (P_hybrid_BA > best_overall_nohybrid_probability)
		    {
			//record "best" values:
			best_hybrid_probability = P_hybrid_BA;
			best_PER_and_hybrid_algmt = PER_and_hybrid_BA_algmt;
			best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt;
			
			//mapping stuff for star alignment:
			map_relative_hybrid_coords_to_algmt_profile = get_relative_hybrid_coords_to_algmt_profile(PER_and_hybrid_BA_algmt, subset_of_hybrid_BA_that_is_aligned);			
			hybrid_breakpoint_indeces = hybrid_BA__breakpoint.second;
			
			if (desired_outcome_is_NAHR)
			{//NAHR
			    best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = mini_A_corresponds_to_LCR_0  ?  3 : 2;
			
			    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first =  mini_B_orientation_agrees_with_MiniEvent_profile  ?
					mini_B_it->second.region_interval.lower()   :   mini_B_it->second.region_interval.upper();
			
			    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first;
			    
			    absolute_Y_coord_at_breakpoint =   mini_A_orientation_agrees_with_MiniEvent_profile   ?
						    absolute_brkpt_interval__on_A.lower()   :   absolute_brkpt_interval__on_A.upper();
			    
			    mini_X_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
			    mini_Y_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
			    algmt_should_be_considered_as_NAHR = true;
			}//NAHR
			else
			{//Geneconv
			    if (compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 0)
			    {//full
				best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 = mini_A_corresponds_to_LCR_0  ?  5 : 4;
			    
				if (mini_B_orientation_agrees_with_MiniEvent_profile)
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_B_it->second.region_interval.lower();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_brkpt_interval__on_B.upper();
				}
				else
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_B_it->second.region_interval.upper();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = absolute_brkpt_interval__on_B.lower();
				}
				
				absolute_Y_coord_at_breakpoint =   mini_A_orientation_agrees_with_MiniEvent_profile   ?
							absolute_brkpt_interval__on_A.lower()   :   absolute_brkpt_interval__on_A.upper();
							
				mini_X_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
				mini_Y_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
				algmt_should_be_considered_as_NAHR = false;
			    }//full
			    else
			    {//only one of the breakpoints
			    
				//mini_A_corresponds_to_LCR_0 = true -->  "BAB" really is "BAB"
				best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5  = mini_A_corresponds_to_LCR_0  ?  3 : 2;
				
				//we se "BA" of "BAB"
				if (mini_B_orientation_agrees_with_MiniEvent_profile)
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_B_it->second.region_interval.lower();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = mini_B_it->second.region_interval.upper();
				}
				else
				{
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first = mini_B_it->second.region_interval.upper();
				    absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second = mini_B_it->second.region_interval.lower();
				}
				
				absolute_Y_coord_at_breakpoint =   mini_A_orientation_agrees_with_MiniEvent_profile   ?
							absolute_brkpt_interval__on_A.lower()   :   absolute_brkpt_interval__on_A.upper();
							
				mini_X_was_reversed = !mini_B_orientation_agrees_with_MiniEvent_profile;
				mini_Y_was_reversed = !mini_A_orientation_agrees_with_MiniEvent_profile;
				
				algmt_should_be_considered_as_NAHR = true;
			    }//only one of the breakpoints
			}//GeneConv
			
		    }//P_hybrid_BA > best_overall_nohybrid_probability
		} //P_hybrid_BA > P_hybrid_AB
		
		
		//record
// 		if (P_hybrid_AB  >  best_overall_nohybrid_probability  or   P_hybrid_BA  >  best_overall_nohybrid_probability)
		{
		    adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
				(map_relative_hybrid_coords_to_algmt_profile,
				hybrid_breakpoint_indeces,
				absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX, mini_X_was_reversed,
				absolute_Y_coord_at_breakpoint, mini_Y_was_reversed,
				algmt_should_be_considered_as_NAHR,
				*proper_map_abs_A_to_algmt, *proper_map_abs_B_to_algmt);
		}//record
		
		
// 			    if (name.compare("SRR359106.64680923") == 0)
// 			    {
// 				std::stringstream debug_strm;
// 				debug_strm << "\n\nname = " << name << "\n"
// 					<< "after \"adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates\"\n"
// 					<< "hybrid_breakpoint_indeces = " << hybrid_breakpoint_indeces << "\n"
// 					<< "absolute_brkpt_interval__on_B = " << absolute_brkpt_interval__on_B << "\n"
// 					<< "absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX = " << absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.first 
// 					<< ",  " << absolute_X_coords__hybrid_index_0___second_brkpt___XY_or_XYX.second << "\n"
// 					<< "absolute_Y_coord_at_breakpoint = " << absolute_Y_coord_at_breakpoint << "\n"
// 					<< "mini_X_was_reversed = " << mini_X_was_reversed << "\n"
// 					<< "mini_Y_was_reversed = " << mini_Y_was_reversed << "\n"
// 					<< "algmt_should_be_considered_as_NAHR  = " << algmt_should_be_considered_as_NAHR  << "\n\n"
// 					<< "(P_hybrid_AB > P_hybrid_BA) = " << (P_hybrid_AB > P_hybrid_BA) << "\n\n";
// 					
// 					print_map_keys_and_values<uint,uint>(*proper_map_abs_A_to_algmt, "*proper_map_abs_A_to_algmt", &debug_strm, true);
// 					print_map_keys_and_values<uint,uint>(*proper_map_abs_B_to_algmt, "*proper_map_abs_B_to_algmt", &debug_strm, true);
// 					
// 					print_map_keys_and_values<uint,uint>(map_relative_hybrid_coords_to_algmt_profile, "map_relative_hybrid_coords_to_algmt_profile", &debug_strm, true);
// 				std::cerr << "\n\n\n" << debug_strm.str() << "\n\n\n";
// 			    }




		for (uint lu=0; lu<2; ++lu)
		{
		    MiniEvent_intersecting_minis_AB->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[mini_A_it->second.MiniID]
							    .insert(pair_at(profile_brkpt_interval,lu));    
							    
		    MiniEvent_intersecting_minis_AB->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[mini_B_it->second.MiniID]
							    .insert(pair_at(profile_brkpt_interval,lu));     
		}
		
		if ( MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID) == 0)
		{
		    MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
					.insert(type_BI__real__real( profile_brkpt_interval, type_real__real(P_hybrid_AB, P_hybrid_BA)));
		}
		else
		{
		    MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
					.insert(type_BI__real__real(profile_brkpt_interval, type_real__real(P_hybrid_BA, P_hybrid_AB)));    
		}
	    }  // relevant for both.

	}  //  partners for A

    }  // end for-loop     mini_A_ptr
    
    
    
    
    
    
    
    
    

    bool total_return_val = true;
    
    
    
    bool should_display__via_AB_ABA;
    bool should_display__via_BA_BAB;
    
    if (desired_outcome_is_NAHR)
    {//NAHR
	should_display__via_AB_ABA =  (best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 2  and  test_if_outcome_should_display__AB_ABA(event_outcome));
	should_display__via_BA_BAB =  (best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 3  and  test_if_outcome_should_display__BA_BAB(event_outcome));
    }//NAHR
    else
    {//Geneconv
	should_display__via_AB_ABA  =   (event_outcome == hap_outcome__GeneConv_ABA)
					and  (     (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 0  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 4)
						or (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 1  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 2)
						or (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 2  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 3));
					
	should_display__via_BA_BAB =    (event_outcome == hap_outcome__GeneConv_BAB)
					and  (     (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 0  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 5)
						or (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 1  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 3)  // "BA"
						or (best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt == 2  and  best_hybrid_type___AB_2____BA_3____ABA_4____BAB_5 == 2));
    }//GeneConv
    
//     assert( !(should_display__via_AB_ABA and should_display__via_BA_BAB));
    
// 	    if (name.compare("SRR017053.5539940") == 0)
// 	    {
// 		std::stringstream debug_strm;
// 		debug_strm << "\n\nname = " << name << "\n"
// 			<< "best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt = " << best_compute_is___0_full_GeneConv___1_first_brkpt___2_second_brkpt << "\n"
// 			<< "should_display__via_AB_ABA = " << should_display__via_AB_ABA << "\n"
// 			<< "should_display__via_BA_BAB = " << should_display__via_BA_BAB << "\n"
// 			<< "\nbest_PER_and_hybrid_algmt:\n" << best_PER_and_hybrid_algmt.first << "\n" << best_PER_and_hybrid_algmt.second << "\n\n"
// 			<< "\nbest_PER_and_hybrid_algmt:\n" << best_PER_and_LCR_0_no_hybrid_algmt.first << "\n" << best_PER_and_LCR_0_no_hybrid_algmt.second << "\n\n"
// 			<< "\nbest_PER_and_hybrid_algmt:\n" << best_PER_and_LCR_1_no_hybrid_algmt.first << "\n" << best_PER_and_LCR_1_no_hybrid_algmt.second << "\n\n";
// 		std::cerr << "\n\n\n" << debug_strm.str() << "\n\n\n";
// 	    }
    

						
//     if (((best_hybrid_probability/best_overall_nohybrid_probability) > 1.005L)   // non-triviality condition
	  if   (should_display__via_AB_ABA  or  should_display__via_BA_BAB)
    {        
	
	const real sum_over_all_locations(sum_over_list<real>(P_every_non_hybrid_location));
	
	const type_map_uint_to_uint dummy_empty_map;
	bool hybrid_return_val = true;
		
	
	if (should_display__via_AB_ABA)   
	{
	    #pragma omp critical (add_display_alignment_to_Star_Alignment__only_one_thread_at_a_time_inside_and_DISPLAY_PERs____AB)
	    hybrid_return_val = star_alignment_for_hybrid_AB.add_pairwise_alignment_to_progressive_star_alignment
					    (best_hybrid_probability,
					    best_PER_and_hybrid_algmt,
					    best_hybrid___absolute_A_coords_to_algmt_profile,
					    best_hybrid___absolute_B_coords_to_algmt_profile,
					    true,
					    false,
					    mates,
					    sum_over_all_locations+best_hybrid_probability);
	}
	
	if (should_display__via_BA_BAB)
	{
	    #pragma omp critical (add_display_alignment_to_Star_Alignment__only_one_thread_at_a_time_inside_and_DISPLAY_PERs____BA)
	    hybrid_return_val = star_alignment_for_hybrid_BA.add_pairwise_alignment_to_progressive_star_alignment
					    (best_hybrid_probability,
					    best_PER_and_hybrid_algmt,
					    best_hybrid___absolute_B_coords_to_algmt_profile,
					    best_hybrid___absolute_A_coords_to_algmt_profile,                                            
					    true,
					    false,
					    mates,
					    sum_over_all_locations+best_hybrid_probability);
	}
					    
	bool LCR_0_return_val = true;
	bool LCR_1_return_val = true;
		
	if (best_LCR_0_nohybrid_MR != NULL)   
	{
	    #pragma omp critical (add_display_alignment_to_Star_Alignment__only_one_thread_at_a_time_inside_and_DISPLAY_PERs____LCR_0)
	    LCR_0_return_val = star_alignment_for_LCR_0.add_pairwise_alignment_to_progressive_star_alignment
						(best_LCR_0_nohybrid_probability,
						best_PER_and_LCR_0_no_hybrid_algmt,
						best_LCR_0_nohybrid___MR_coords_absolute_to_algmt_profile,
						dummy_empty_map,                                            
						best_LCR_0_nohybrid_MR->orientation_of_MiniRegion_and_spawning_Event_profile_agree,
						best_LCR_0_nohybrid_MR->MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile,
						mates,
						sum_over_all_locations+best_hybrid_probability); 
	}
						
	if (best_LCR_1_nohybrid_MR != NULL)              
	{
	    #pragma omp critical (add_display_alignment_to_Star_Alignment__only_one_thread_at_a_time_inside_and_DISPLAY_PERs____LCR_1)
	    LCR_1_return_val = star_alignment_for_LCR_1.add_pairwise_alignment_to_progressive_star_alignment
						(best_LCR_1_nohybrid_probability,
						best_PER_and_LCR_1_no_hybrid_algmt,
						best_LCR_1_nohybrid___MR_coords_absolute_to_algmt_profile,
						dummy_empty_map,                                            
						best_LCR_1_nohybrid_MR->orientation_of_MiniRegion_and_spawning_Event_profile_agree,
						best_LCR_1_nohybrid_MR->MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile,
						mates,
						sum_over_all_locations+best_hybrid_probability);
	}
						
						
	total_return_val = (hybrid_return_val  and   LCR_0_return_val   and    LCR_1_return_val);                                    
    }
    
        
    return total_return_val;  

}  //  end of   align_to_every_possible_hybrid_outcome_AND_DISPLAY_REGARDLESS_OF_NOHYBRIDS














































void  Paired_end_read::deep_copy_from_serialize
                                (const Paired_end_read &another_PER,
                                 const Conn_Comp *const &the_conn_comp_ptr)
{
    
    {//PER
        name = another_PER.name;
                    //  const std::string name;          --     is const, can't be changed.  
                
        //  const Readgroup_statistics *my_Readgroup_ptr;       // RECONNECT
        
        mates = another_PER.mates;
//         mate_sequences = another_PER.mate_sequences;
        
        spawning_profile_chromo_and_position = another_PER.spawning_profile_chromo_and_position;    
        
        
        my_MiniEvents = another_PER.my_MiniEvents;   // RECONNECT!!!    
        my_MiniRegions = another_PER.my_MiniRegions;
        
        set_of_MiniRegions_on__XXXX_chromosome = another_PER.set_of_MiniRegions_on__XXXX_chromosome;
        set_of_MiniRegions_on__YYYY_chromosome = another_PER.set_of_MiniRegions_on__YYYY_chromosome;
        
        map_affected_MiniID__to__P_nohybrid = another_PER.map_affected_MiniID__to__P_nohybrid;
        
        haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal
                        = another_PER.haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal;
                        
        haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome
                        = another_PER.haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome;
                        
        haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome
                        = another_PER.haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome;                        
        

        MiniID_ctr = another_PER.MiniID_ctr;
        
	failed_PER_flag = another_PER.failed_PER_flag;
	
	map_MiniID__to__nohybrid_almgt_interval_absolute = another_PER.map_MiniID__to__nohybrid_almgt_interval_absolute;
    }//PER
    
    



    
    
    reconnect_pointers_after_uploading_from_a_serialization( the_conn_comp_ptr );
    
    
    
//     type_map_uint_to_MiniEvent::iterator insert_it = my_MiniEvents.begin();
//     
//     for (type_map_uint_to_MiniEvent::const_iterator it_another_ME = another_PER.my_MiniEvents.begin();
//             it_another_ME != another_PER.my_MiniEvents.end();
//             ++it_another_ME)
//     {
//         insert_it = my_MiniEvents.insert( insert_it,    std::pair<uint, MiniEvent>
//                                                        (  it_another_ME->first,
//                                                           MiniEvent(    &(the_conn_comp.events.find(it_another_ME->first)->second),     this  ))    );
//         
//         insert_it->second.affected_MiniRegions = it_another_ME->second.affected_MiniRegions;
//         insert_it->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR = it_another_ME->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR;
//         insert_it->second.MiniRegions_in_Inbetween_of_this_MiniEvent = it_another_ME->second.MiniRegions_in_Inbetween_of_this_MiniEvent;
//         insert_it->second.partners_for_MiniRegions_on_LCRs = it_another_ME->second.partners_for_MiniRegions_on_LCRs;
//         insert_it->second.MiniRegions_in_insertion_of_an_LCR = it_another_ME->second.MiniRegions_in_insertion_of_an_LCR;
//         insert_it->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
//                              = it_another_ME->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA;
//         insert_it->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion
//                      = it_another_ME->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion;
//         insert_it->second.relevant_breakpoints_for_which_this_PER_actually_shows_a_hybrid_pattern
//                      = it_another_ME->second.relevant_breakpoints_for_which_this_PER_actually_shows_a_hybrid_pattern;                     
//     }
    
    
    
    

}  //   deep_copy_from_serialize        

























void Paired_end_read::reconnect_pointers_after_uploading_from_a_serialization
                                                (const Conn_Comp *const &the_conn_comp_ptr)
{
	
        
    //PER
	
	if (!is_simulation)
	{	
		std::string my_readgroup_name;
		{
			const int loc_of_dot = name.find('.');
			if (loc_of_dot == std::string::npos)
			{
				std::stringstream error_strm;
				print_line_of_markers("ERROR! ", &error_strm);
				print_line_of_markers("(", &error_strm);
				error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
				
				error_strm << "\n\nERROR!    could not find  \".\" in readgroup name for this read:  name = [" << name << "] !!!!\n\n";
				error_strm << "\n\tassigning arbitrary RG to it so the program doesn't crash...\n\n";
				
				print_this_PER( &error_strm );
				
				print_line_of_markers(")", &error_strm);
				std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                
				
				my_readgroup_name.assign(  Readgroup_stats.begin()->first  );
			}
			else    
			{
				my_readgroup_name.assign(   name.substr(0, loc_of_dot)   );
			}
		}
    
		if (!failed_PER_flag)    
		my_Readgroup_ptr = &Readgroup_stats.at(my_readgroup_name);     
		else
		my_Readgroup_ptr = NULL;
	}
	else // is_simulation
	{
		my_Readgroup_ptr = &(Readgroup_stats.begin()->second);
	}
    
    
    
    
    
    //MiniEvents
    for (type_map_uint_to_MiniEvent::iterator it_ME = my_MiniEvents.begin();
            it_ME != my_MiniEvents.end();
            ++it_ME)
    {
        it_ME->second.full_Event = &(the_conn_comp_ptr->events.at( it_ME->first ));
        it_ME->second.parent_PER = this;
    }




    //PER_aligner
    my_PER_aligner.my_PER_ptr = this;
        

}  //  reconnect_pointers_after_uploading_from_a_serialization                                  


















const MiniRegion *const Paired_end_read::get_ptr_to_first_MiniRegion()
{
    assert(  !my_MiniRegions.empty()   );
    return &(my_MiniRegions.begin()->second);
} // get_ptr_to_first_MiniRegion















void Paired_end_read::set_sequence_for_all_MiniRegions
                        (type_map_uint_to_list_BI__string &EV_uploaded_pieces_of_chromos)
{
    bool success_for_all_MRs = true;
    for (type_map_uint_to_MiniRegion::iterator it_MR = my_MiniRegions.begin();
            it_MR != my_MiniRegions.end();
            ++it_MR)
    {
        const bool success_found =  it_MR->second.set_region_sequence_from_preloaded_sequence(EV_uploaded_pieces_of_chromos);
        
        if (!success_found)
	{  success_for_all_MRs = false;  }
    }
    
    
    if ( !success_for_all_MRs )
        print_this_PER();
    

}//set_sequence_for_all_MiniRegions



















void Paired_end_read::extend_MRs_to_minimum_length()
{
//     const int sum_read_lengths = mates.first.Length + mates.second.Length;
//     
//     const type_int__int adv_bounds = my_PER_aligner.get_advance_bounds();
//     const int min_frag_size_this_PER = std::min<int>(  sum_read_lengths + adv_bounds.first,
//                                                         sum_read_lengths + adv_bounds.second  );    

    const int min_frag_size_this_PER = my_Readgroup_ptr->median_insert_size + my_Readgroup_ptr->median_insert_size____absolute_standard_deviation*2;
    
    for (type_map_uint_to_MiniRegion::iterator it_MR = my_MiniRegions.begin();
            it_MR != my_MiniRegions.end();
            ++it_MR)
    {
	if (BOOST_width_inclusive(it_MR->second.region_interval) <  min_frag_size_this_PER)
	{
	    const int diff_in_size = min_frag_size_this_PER - (int)BOOST_width_inclusive(it_MR->second.region_interval);
	    	    
	    it_MR->second.region_interval.set(  safe_subtract_base_1(it_MR->second.region_interval.lower(), diff_in_size + 10),
						safe_chromo_add_base_1(it_MR->second.region_interval.upper(), diff_in_size + 10, it_MR->second.chromosome_of_region));
		
	} // it_MR
    }


} // extend_MRs_to_minimum_length







































// bool Paired_end_read::align_PER_to_MAP_location
//                             (const uint &EV_UID,
//                              const BOOST_Interval &desired_breakpoint,
//                              Star_alignment  &star_alignment_for_hybrid_AB,
//                              Star_alignment  &star_alignment_for_hybrid_BA,
// 			     type_map_uint_to_list_BI__Star_alignment   &star_alignment__non_hybrids)
// {    
//     
//     uint number_of_alignments_performed = 0;   
//     
//     const BOOST_Interval padded_brkpts_of_interest( safely_subtract(desired_breakpoint.lower(),15), desired_breakpoint.upper() + 15 );    
//     
//                                                     
//     type_list_real  P_every_possible_location;          
//                                                     
//     //best:
//     
//     
//     //overall:    
//     real best_overall_probability(0.00L);
//     
//     const MiniRegion *best_overall_nohybrid_MR = NULL;
//     type_string__BI best_hybrid_sequence_and_breakpoint;      
//     type_string__string best_PER_and_overall_algmt;    
//     bool best_overall_alignment_is_hybrid_AB;
//     
//     type_map_uint_to_uint  best_hybrid___absolute_A_coords_to_algmt_profile;
//     type_map_uint_to_uint  best_hybrid___absolute_B_coords_to_algmt_profile;    
//     
//     bool best_orientation;
//     bool best_complementarity;
//     
//     
//     
//     
//     
//     
//     
//     for (type_map_uint_to_MiniRegion::const_iterator it_mini_A = my_MiniRegions.begin();
//             it_mini_A != my_MiniRegions.end();
//             ++it_mini_A)
//     {
// 
//         BOOST_Interval mini_A_algmt_endpts;
//         type_string__string PER_and_mini_A_algmt;
//                                     
//         const real P_nohybrid_A(
//                         my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
//                                                 (mates.first,
//                                                 mates.second,
//                                                 it_mini_A->second.region_sequence,
//                                                 it_mini_A->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree,
//                                                 it_mini_A->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile,
//                                                 PER_and_mini_A_algmt,
//                                                 mini_A_algmt_endpts)             );     
// 						
// 	++number_of_alignments_performed;
// 	
// 	P_every_possible_location.push_back(P_nohybrid_A);
// 						
//    
//                                                                                     
// 
//         const type_map_uint_to_uint  map_mini_A_indeces_to_algmt_profile(  get_relative_hybrid_coords_to_algmt_profile(PER_and_mini_A_algmt,mini_A_algmt_endpts)   );
//         
//         if (   (--map_mini_A_indeces_to_algmt_profile.end())->first  >=   it_mini_A->second.region_sequence.size() )
//         {
//             std::stringstream error_strm;
//             print_line_of_markers("ERROR! ", &error_strm);
//             print_line_of_markers("(", &error_strm);
//             error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
//             
//             error_strm << "\n\nERROR!   in no-hybrid displays:\n\n\n"
//                         << "(--map_mini_A_indeces_to_algmt_profile.end())->first  =  "  <<  (--map_mini_A_indeces_to_algmt_profile.end())->first
//                         << "   >=     " << it_mini_A->second.region_sequence.size() <<" =  it_mini_A->second.region_sequence.size()\n\n";
//                         
//             error_strm  << "\nmini_A_algmt_endpts = " << mini_A_algmt_endpts
//                         << "\n\nPER_and_mini_A_algmt:\nPER:\n" << PER_and_mini_A_algmt.first
//                         << "\nmini:\n" << PER_and_mini_A_algmt.second
//                         << "\n\n";
//                         
//             it_mini_A->second.print_this_MiniRegion( &error_strm );
//             
//             print_this_PER(  &error_strm  );
//             
//             print_map_keys_and_values<uint,uint>(map_mini_A_indeces_to_algmt_profile, "map_mini_A_indeces_to_algmt_profile",  &error_strm,  true );
//                         
//             
//             print_line_of_markers(")", &error_strm);
//             std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                
//         }
//         
//         
//         
//         
//         
//         type_map_uint_to_uint  dummy_empty_map;
//         type_map_uint_to_uint  map_mini_A_absolute_to_algmt_profile;
//                            adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
//                                         (  map_mini_A_indeces_to_algmt_profile,
//                                             BOOST_make_point(BOOST_width_inclusive(it_mini_A->second.region_interval)),  //i.e.  beyond the end
//                                             it_mini_A->second.region_interval.lower(),
//                                             false,  //never reversed the MR,  only possibly the PER
//                                             0,  //irrelevant
//                                             false,   //irrelevant
// 					   true,
//                                             map_mini_A_absolute_to_algmt_profile,   
//                                             dummy_empty_map     );           //irrelevant     
//                                         
//         assert(  !map_mini_A_absolute_to_algmt_profile.empty()  );                                        
//                      
//         
//         if (P_nohybrid_A  > best_overall_probability)
//         {
//             best_overall_probability = P_nohybrid_A;
//                 
//             best_overall_nohybrid_MR = &(it_mini_A->second);
//             best_PER_and_overall_algmt = PER_and_mini_A_algmt;
//             
//             best_hybrid___absolute_A_coords_to_algmt_profile = map_mini_A_absolute_to_algmt_profile;   
// 	    best_hybrid___absolute_B_coords_to_algmt_profile = dummy_empty_map;
// 
//             best_orientation = it_mini_A->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree;
//             best_complementarity = it_mini_A->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile;            
//         }
//         
// 
//     } // end of no-hybrid
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
//         
//         
//         
//         
//         
//     const type_map_uint_to_MiniEvent::iterator MiniEvent_intersecting_minis_AB = my_MiniEvents.find(EV_UID);    
//     if ( MiniEvent_intersecting_minis_AB != my_MiniEvents.end() )
//     {
//     
// 	const Event *const Event_intersecting_minis_AB = MiniEvent_intersecting_minis_AB->second.full_Event;   //readability
//                 
// 
//         //need to cycle through because some varpos can only be caught from one side (e.g. inserts relative to other mini)
//         for ( type_set_uint::const_iterator it_Mini_A_MiniID = MiniEvent_intersecting_minis_AB->second.affected_MiniRegions.begin();
//                 it_Mini_A_MiniID != MiniEvent_intersecting_minis_AB->second.affected_MiniRegions.end();
//                 ++it_Mini_A_MiniID)
//         {                        
//             
//             
//             if (     MiniEvent_intersecting_minis_AB->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(*it_Mini_A_MiniID) > 0 
//                  or  MiniEvent_intersecting_minis_AB->second.MiniRegions_in_insertion_of_an_LCR.count( *it_Mini_A_MiniID )  > 0       )
//                 continue;                                                                   
//                                     
//             //else, has a partner and is potentially a hybrid read.                                      
//             
//             
//             const type_map_uint_to_MiniRegion::const_iterator  mini_A_it = my_MiniRegions.find(*it_Mini_A_MiniID);  //readability.
//                                                                                     
//             
//             const type_set_uint *const ptr_to_set_of_MiniRegion_partners_for_MiniRegion_A 
//                         = &MiniEvent_intersecting_minis_AB->second.partners_for_MiniRegions_on_LCRs.at( *it_Mini_A_MiniID );
// 	    const type_map_uint_to_set_uint::const_iterator  it_to_set_of_MiniRegion_partners_for_MiniRegion_A
//                                 =  MiniEvent_intersecting_minis_AB->second.partners_for_MiniRegions_on_LCRs.find( *it_Mini_A_MiniID );
//             
//             for (type_set_uint::const_iterator it_MiniB_MiniID = it_to_set_of_MiniRegion_partners_for_MiniRegion_A->second.begin();
//                     it_MiniB_MiniID != it_to_set_of_MiniRegion_partners_for_MiniRegion_A->second.end();
//                     ++it_MiniB_MiniID)
//             {
//                                                          
// 		const type_map_uint_to_MiniRegion::const_iterator mini_B_it =  MiniEvent_intersecting_minis_AB->second.parent_PER->my_MiniRegions.find(*it_MiniB_MiniID);
//                 
//                 
//                 assert( MiniEvent_intersecting_minis_AB->second.MiniRegions_in_Inbetween_of_this_MiniEvent.count(mini_B_it->second.MiniID) == 0  );
//                 assert( MiniEvent_intersecting_minis_AB->second.MiniRegions_in_insertion_of_an_LCR.count(mini_B_it->second.MiniID) == 0  );
// 
// 
// 
//                 const  uint mini_A_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID);
//                 const  uint mini_B_qs = MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_B_it->second.MiniID);   
//                 
// 
// 
// 
//                 // get some needed maps and  some crucial, helpful characteristics:
// 
//                 // For A:            
//                 const type_map_uint_to_uint map_MiniEvent_profile_varpos_only_to_mini_A(
//                             convert_compressed_map_to_full_map(
//                                     get_compressed_map_inverse_of_compressed_map(
//                                         get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
//                                                 Event_intersecting_minis_AB->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(mini_A_qs),
//                                                 mini_A_it->second.region_interval   )  ) ));
// 
// 
// 
//                 const type_map_uint_to_uint map_mini_A_to_MiniEvent_profile(
//                         convert_compressed_map_to_full_map(
//                             get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
//                                     Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at( mini_A_qs ),
//                                     mini_A_it->second.region_interval  )));    
//                         
//                 type_map_uint_to_uint map_MiniEvent_profile_to_mini_A( get_inverse_of_map(map_mini_A_to_MiniEvent_profile) );  // (see explanation for mini_B version)
//                 
// 
//                 
// 
//                 const bool mini_A_orientation_agrees_with_MiniEvent_profile
//                                 = (map_MiniEvent_profile_to_mini_A.begin()->second <  (++map_MiniEvent_profile_to_mini_A.begin())->second);            
// 
//                                 
//                 const bool must_complement_mini_A_to_get_to_MiniEvent_profile =  ( mini_A_qs == 1 &&  Event_intersecting_minis_AB->use_complement_of_subject );                                
//             
//                             
// 
//                 // For B:
//                 const type_map_uint_to_uint map_MiniEvent_profile_to_mini_B(
//                             convert_compressed_map_to_full_map(
//                                     get_compressed_map_inverse_of_compressed_map(
//                                             get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
//                                                   Event_intersecting_minis_AB->compressed_map_LCR_to_profile__for_each_LCR.at( mini_B_qs ),
//                                                   mini_B_it->second.region_interval  )) ) );  
//                 //if the breakpoint is at a variational position of mini_A that does NOT appear in mini_B (e.g. the breakpoint is an insertion in A relative to B), then we need to know the full sequence of B in order to know what the "first" nucleotide of B would be in the hybrid, e.g. exactly how the mini sequences would hybridize.
// 
// 
//                 const bool mini_B_orientation_agrees_with_MiniEvent_profile
//                                 = (  map_MiniEvent_profile_to_mini_B.begin()->second   <   (++map_MiniEvent_profile_to_mini_B.begin())->second  );            
// 
//                                 
//                 const bool must_complement_mini_B_to_get_to_MiniEvent_profile =  ( mini_B_qs == 1   &&   Event_intersecting_minis_AB->use_complement_of_subject );
// 
//                 
// 
//                         
// 
//                 // ??????????????????????????????????????????????????????????????????????????????????????????????????????????????
//                 const bool orientation_of_PER_and_ev_profile_agree
//                     =  ( mini_A_it->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree ==  mini_A_orientation_agrees_with_MiniEvent_profile );
//                 const bool must_complement_PER_to_get_to_ev_profile
//                     =  ( mini_A_it->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile != must_complement_mini_A_to_get_to_MiniEvent_profile  );
// 
// 
//                     
//                 //if miniRegions are too far apart wrt profile, then they will never bee seen together in the same PER even if there is a hybrid (they must come from different homologosu regions of the LCR).
//                 if ( !BOOST_overlap(
//                                     BOOST_Interval( map_MiniEvent_profile_to_mini_A.begin()->first, (--map_MiniEvent_profile_to_mini_A.end())->first ),
//                                     BOOST_Interval( map_MiniEvent_profile_to_mini_B.begin()->first, (--map_MiniEvent_profile_to_mini_B.end())->first )   )   )
//                 {
//                     print_line_of_markers("ERROR! ");
//                     print_line_of_markers("(");
//                     std::fprintf(stderr, "mini_A_it->second.MiniID = %u   and   mini_B_it->second.MiniID = %u   DO NOT INTERSECT!!!\n",
//                                 mini_A_it->second.MiniID,
//                                 mini_B_it->second.MiniID );
//                     mini_A_it->second.print_this_MiniRegion();
//                     mini_B_it->second.print_this_MiniRegion();
//                     print_line_of_markers(")");                
//                 }
// 
//                 map_MiniEvent_profile_to_mini_A.clear();
// 
// 
//     
// 
// 
//                 
//                 for (type_map_uint_to_uint::const_iterator it_ME_varpos_to_mini_A__low = map_MiniEvent_profile_varpos_only_to_mini_A.begin();
//                         it_ME_varpos_to_mini_A__low != map_MiniEvent_profile_varpos_only_to_mini_A.end(); 
//                         ++it_ME_varpos_to_mini_A__low)
//                 {                     		    
// 		    for (type_map_uint_to_uint::const_iterator it_ME_varpos_to_mini_A__high = it_ME_varpos_to_mini_A__low;
// 			    it_ME_varpos_to_mini_A__high != map_MiniEvent_profile_varpos_only_to_mini_A.end(); 
// 			    ++it_ME_varpos_to_mini_A__high)
// 		    { 
// 			
// 			const BOOST_Interval profile_brkpt_interval(it_ME_varpos_to_mini_A__low->first, it_ME_varpos_to_mini_A__high->first);
// 			
// 			const bool is_NAHR = BOOST_is_a_point(profile_brkpt_interval);
// 			
// 			//if already recorded, skip.
// 			if ( MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA.count(profile_brkpt_interval) > 0 )
// 			    continue;                    
// 			
// 
// 			const type_map_uint_to_uint::const_iterator it_found_B__low = find_key_in_map_or_find_next_key_in_specified_direction_WRT_keys<uint,uint>
// 										    (map_MiniEvent_profile_to_mini_B,
// 										    it_ME_varpos_to_mini_A__low->first,
// 										    true);  //since is always wrt profile (regardless if B reversed or not)
// 										    
// 			const type_map_uint_to_uint::const_iterator it_found_B__high = find_key_in_map_or_find_next_key_in_specified_direction_WRT_keys<uint,uint>
// 										    (map_MiniEvent_profile_to_mini_B,
// 										    it_ME_varpos_to_mini_A__high->first,
// 										    true);  //since is always wrt profile (regardless if B reversed or not)										    
// 			// (we use "map_MiniEvent_profile_to_mini_B" and not "map_MiniEvent_profile_varpos_only_to_mini_B"   here because miniB may have other variational positions (i.e. varpos related to other events, that come before the varpos for THIS PARTICULAR MiniEvent).
// 
// 			if ( it_found_B__low  == map_MiniEvent_profile_to_mini_B.end()   or   it_found_B__high == map_MiniEvent_profile_to_mini_B.end()  )
// 			    continue;                    
// 
// // 			if ( ! BOOST_in(it_found_B->second, mini_B_it->second.region_interval) )
// // 			    continue;
// 			    
// 			    
// 			const BOOST_Interval absolute_brkpt_interval__on_A(  std::min<uint>(it_ME_varpos_to_mini_A__low->second, it_ME_varpos_to_mini_A__high->second),
// 									     std::max<uint>(it_ME_varpos_to_mini_A__low->second, it_ME_varpos_to_mini_A__high->second));	
// 
// 			const BOOST_Interval absolute_brkpt_interval__on_B(  std::min<uint>(it_found_B__low->second, it_found_B__high->second),
// 									     std::max<uint>(it_found_B__low->second, it_found_B__high->second));									     
// 									     
// 			    
// 
// 
// 			if ( !BOOST_overlap_and_share_some_bound(profile_brkpt_interval, desired_breakpoint) )
// 			{  continue;  }  
// 			
// 			
// 			//so     it_map_A:  mini_coord  --->  varpos         it_found_B: same_varpos  --->  mini_B_coord
// 			//so our hybrid will be:  sequence A  from [beginning_of_mini_A,  it_map_A->first) NOT INCLUSIVE (since the breakpoint reflects the NEW sequence)
// 			// and then               sequence B  from [it_found_B->second, end_of_mini_B]
// 			//Note:  "beginning_of_mini_A" and "end_of_mini_B"  dpened on orientations (beginning may in fact mean end, and end may in fact mean beginning).
// 
// 
// 			//We orient ourselves WITH the PROFILE.  Here, we want the lower part of the profile to be sequence A, and the upper part of the profile to be sequence B.  Depending on the orientations of A,B to the profile, we need to play some games here....
// 			const type_string__BI hybrid_AB__breakpoint(
// 						    get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
// 									    (mini_A_it->second,
// 									    absolute_brkpt_interval__on_A,    //not to be included!!
// 									    must_complement_mini_A_to_get_to_MiniEvent_profile,
// 									    !mini_A_orientation_agrees_with_MiniEvent_profile,   //always want beginning of A
// 									    mini_B_it->second,
// 									    absolute_brkpt_interval__on_B,
// 									    must_complement_mini_B_to_get_to_MiniEvent_profile,
// 									    mini_B_orientation_agrees_with_MiniEvent_profile,    // take_this_side_of_B_sequence)   );    //CAREFUL!!!
// 									    is_NAHR));  // note that the breakpoint is relative to the length of the hybrid                                
// 
// 			
// 			//and the flip side....
// 			//We orient ourselves WITH the PROFILE.  Here, we want the lower part of the profile to be sequence B, and the upper part of the profile to be sequence A.  Depending on the orientations of A,B to the profile, we need to play some games here....
// 			const type_string__BI hybrid_BA__breakpoint(
// 							get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
// 								    (mini_B_it->second,
// 								    absolute_brkpt_interval__on_B,  //profile to mini_b.   //not to be included!!
// 								    must_complement_mini_B_to_get_to_MiniEvent_profile,
// 								    !mini_B_orientation_agrees_with_MiniEvent_profile,    // tricky!
// 								    mini_A_it->second,
// 								    absolute_brkpt_interval__on_A,  //inclusive
// 								    must_complement_mini_A_to_get_to_MiniEvent_profile,
// 								    mini_A_orientation_agrees_with_MiniEvent_profile,  //take_this_side_of_A_sequence
// 								    is_NAHR));  // note that the breakpoint is relative to the length of the hybrid  
// 
// 			//This hybrid string is now properly oriented and complemented so that it maps directly to "Event_intersecting_minis_AB" 's profile.
// 			//Thus, the orientation/complementation of the PERs should be based on "Event_intersecting_minis_AB" 's profile.
// 	    
// 
// 			if (   determine_observed_frag_length() > hybrid_AB__breakpoint.first.size()
// 			    or determine_observed_frag_length() > hybrid_BA__breakpoint.first.size())
// 			{                                                    
// 			    continue;
// 			}//suff length
// 
// 
// 
// 			//Now, we want to align the PER to these...
// 
// 			BOOST_Interval subset_of_hybrid_AB_that_is_aligned;
// 			type_string__string PER_and_hybrid_AB_algmt;
// 
// 					
// 			const real P_hybrid_AB(
// 					    my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
// 								    (mates.first,
// 								    mates.second,
// 								    hybrid_AB__breakpoint.first,
// 								    orientation_of_PER_and_ev_profile_agree,
// 								    must_complement_PER_to_get_to_ev_profile,                                             
// 								    PER_and_hybrid_AB_algmt,
// 								    subset_of_hybrid_AB_that_is_aligned)             ); 
// 			++number_of_alignments_performed;
// 		    
// 			P_every_possible_location.push_back(P_hybrid_AB);
// 			
// 
// 			//Finally, go through and calculate the proability according to the varpos and taking note of the breakpoint...
// 			// AND TAKE NOTE OF THE ORIENTATION OF B!!!!!!
// 
// 
// 			const uint A_mini_offset = mini_A_orientation_agrees_with_MiniEvent_profile   ? 
// 							    mini_A_it->second.region_interval.lower()  :  mini_A_it->second.region_interval.upper();
// 															
// 
// 			//BA:
// 
// 			//This hybrid string is now properly oriented and complemented so that it maps directly to "Event_intersecting_minis_AB" 's profile.
// 			//Thus, the orientation/complementation of the PERs should be based on "Event_intersecting_minis_AB" 's profile...
// 			//
// 			// or... maybe we can cleverly orient/complement hybrid_AB so that it will agree with the PERs........!!!!!!!!!!!!!!!!!!!!
// 
// 
// 
// 			//Now, we want to align the PER to these...
// 
//     //                     type_vector_map_uint_to_uint PER_and_hybrid_BA_algmt_coords;  //these will be relative for each
// 			BOOST_Interval subset_of_hybrid_BA_that_is_aligned;
// 			type_string__string PER_and_hybrid_BA_algmt;
// 
// 
// 	//                 std::fprintf(stderr, "\nalign BA...");                
// 			const real P_hybrid_BA(
// 					my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
// 								(mates.first,
// 								mates.second,
// 								hybrid_BA__breakpoint.first,                                                  
// 								orientation_of_PER_and_ev_profile_agree,
// 								must_complement_PER_to_get_to_ev_profile,
// 								PER_and_hybrid_BA_algmt,
// 								subset_of_hybrid_BA_that_is_aligned )                 );
// 			++number_of_alignments_performed;        
// 			
// 			P_every_possible_location.push_back(P_hybrid_BA);                                                          
// 							
// 			//BUT THE PERS ARE ORIENTATED/COMPLEMENTED TO FIT THE VARI_POS PROFILE, NOT THIS INTERSECTING EVENT'S PROFILE!!!!!!!!!!!!!!!!!!!!!!!
// 
// 			//Finally, go through and calculate the proability according to the varpos and taking note of the breakpoint...
// 			// AND TAKE NOTE OF THE ORIENTATION OF A!!!!!!
// 
// 
// 
// 			const uint B_mini_offset = mini_B_orientation_agrees_with_MiniEvent_profile  ?
// 							    mini_B_it->second.region_interval.lower()  :   mini_B_it->second.region_interval.upper();
// 		
// 							    
// 							    
// 	    
//                     
//                     
//                     
// 			//Here, we determine if a breakpoint is relevant or not.  If the breakpoint did not oiccur within the alignment of the PER to the hybrid sequence, then the alignment is (in theory) equivalent to aligning to the "no-hybrid" version of the appropriate seuqnece, which was done at the beginning of this function.  Thus, there is no need to save this probability again, and we shouldn't mark such a breakpoint as relevant.  Otherwise, we save it
// 			//bool breakpoint_is_relevant__and_so_we_should_record_the_probability = false;
// 		    
// 								
// 			const bool breakpoint_is_relevant_for_AB = BOOST_overlap(hybrid_AB__breakpoint.second, subset_of_hybrid_AB_that_is_aligned);
// 			const bool breakpoint_is_relevant_for_BA = BOOST_overlap(hybrid_BA__breakpoint.second, subset_of_hybrid_BA_that_is_aligned);          
// 
// 
// 				
// 			if (breakpoint_is_relevant_for_AB   or   breakpoint_is_relevant_for_BA)
// 			{                   
// 			    
// 			    //inside the "align_to_every..." loop, MiniRegions "A" and "B" have no relation to "LCR_0" and "LCR_1" whatsoever.
// 			    // we can check though, if A really matches to LCR_0, etc..
// 			    //but regardless, hybrid AB is always formed by orienting the hybrid to the ME profile.
// 			    
// 			    
// 			    const bool AB_really_is_AB  =  (MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID) == 0);
// 			    
// 			    if (P_hybrid_AB > P_hybrid_BA)       
// 			    {
// 				if (P_hybrid_AB > best_overall_probability)
// 				{				                                   
// 				    best_overall_probability = P_hybrid_AB;
// 				    best_overall_alignment_is_hybrid_AB =  AB_really_is_AB;
// 				    
// 				    best_overall_nohybrid_MR = NULL;
// 				    
// 				    best_hybrid_sequence_and_breakpoint = hybrid_AB__breakpoint;
// 				    best_PER_and_overall_algmt = PER_and_hybrid_AB_algmt;
// 				    
// 				    best_orientation = true;
// 				    best_complementarity = false;
// 				    
// 				    type_map_uint_to_uint *const proper_map_abs_A_to_algmt  =    AB_really_is_AB   ?  
// 					    &best_hybrid___absolute_A_coords_to_algmt_profile  :  &best_hybrid___absolute_B_coords_to_algmt_profile;
// 					    
// 				    type_map_uint_to_uint *const proper_map_abs_B_to_algmt  =    AB_really_is_AB   ?  
// 					    &best_hybrid___absolute_B_coords_to_algmt_profile  :  &best_hybrid___absolute_A_coords_to_algmt_profile;                                          
// 				    
// 				    
// 				    type_uint__uint absolute_A_coords__hybrid_index_0___second_brkpt(A_mini_offset,A_mini_offset);
// 				    if (!is_NAHR)//Geneconv
// 				    {  
// 					if (!mini_A_orientation_agrees_with_MiniEvent_profile)			
// 					    absolute_A_coords__hybrid_index_0___second_brkpt.second = absolute_brkpt_interval__on_A.lower();//A was reversed
// 					else
// 					    absolute_A_coords__hybrid_index_0___second_brkpt.second = absolute_brkpt_interval__on_A.upper();//A NOT reversed			
// 				    }//GeneConv
// 				    
// 				    
// 				    const uint absolute_B_coord_at_breakpoint  		
// 						=  (!mini_B_orientation_agrees_with_MiniEvent_profile)  ?
// 							absolute_brkpt_interval__on_B.upper()  :   absolute_brkpt_interval__on_B.lower();
// 
// 					    
// 				    adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
// 						( get_relative_hybrid_coords_to_algmt_profile(PER_and_hybrid_AB_algmt, subset_of_hybrid_AB_that_is_aligned),
// 						hybrid_AB__breakpoint.second,
// 						absolute_A_coords__hybrid_index_0___second_brkpt,
// 						!mini_A_orientation_agrees_with_MiniEvent_profile,
// 						absolute_B_coord_at_breakpoint,
// 						!mini_B_orientation_agrees_with_MiniEvent_profile,
// 						is_NAHR,
// 						*proper_map_abs_A_to_algmt,
// 						*proper_map_abs_B_to_algmt);                                                                                     
// 				}
// 			    }// AB is better
// 			    else
// 			    {
// 				if (P_hybrid_BA > best_overall_probability)
// 				{				       
// 				    
// 				    best_overall_probability = P_hybrid_BA;
// 				    best_overall_alignment_is_hybrid_AB =  !AB_really_is_AB;  //now is BA
// 				    
// 				    best_overall_nohybrid_MR = NULL;
// 				    
// 				    best_hybrid_sequence_and_breakpoint = hybrid_BA__breakpoint;
// 				    best_PER_and_overall_algmt = PER_and_hybrid_BA_algmt;
// 				    
// 				    best_orientation = true;
// 				    best_complementarity = false;				
// 				    
// 				    type_map_uint_to_uint *const proper_map_abs_A_to_algmt  =    AB_really_is_AB   ?  
// 					    &best_hybrid___absolute_A_coords_to_algmt_profile  :  &best_hybrid___absolute_B_coords_to_algmt_profile;
// 					    
// 				    type_map_uint_to_uint *const proper_map_abs_B_to_algmt  =    AB_really_is_AB   ?  
// 					    &best_hybrid___absolute_B_coords_to_algmt_profile  :  &best_hybrid___absolute_A_coords_to_algmt_profile;    
// 					    
// 					    					    					    
// 					    
// 				    type_uint__uint absolute_B_coords__hybrid_index_0___second_brkpt(B_mini_offset,B_mini_offset);
// 				    if (!is_NAHR)//Geneconv
// 				    {  
// 					if (!mini_B_orientation_agrees_with_MiniEvent_profile)			
// 					    absolute_B_coords__hybrid_index_0___second_brkpt.second = absolute_brkpt_interval__on_B.lower();//B was reversed
// 					else
// 					    absolute_B_coords__hybrid_index_0___second_brkpt.second = absolute_brkpt_interval__on_B.upper();//B NOT reversed			
// 				    }//GeneConv				    
// 				    
// 				    const uint absolute_A_coord_at_breakpoint  		
// 						=  (!mini_A_orientation_agrees_with_MiniEvent_profile)  ?
// 							absolute_brkpt_interval__on_A.upper()  :   absolute_brkpt_interval__on_A.lower();					    
// 				    
// 				    adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
// 						( get_relative_hybrid_coords_to_algmt_profile(PER_and_hybrid_BA_algmt, subset_of_hybrid_BA_that_is_aligned),
// 						hybrid_BA__breakpoint.second,
// 						absolute_B_coords__hybrid_index_0___second_brkpt,
// 						!mini_B_orientation_agrees_with_MiniEvent_profile,
// 						absolute_A_coord_at_breakpoint,
// 						!mini_A_orientation_agrees_with_MiniEvent_profile,
// 						is_NAHR,
// 						*proper_map_abs_B_to_algmt,
// 						*proper_map_abs_A_to_algmt  );                                                                                                                                         
// 				}
// 			    } // BA is better
// 			}
// 			    
// 			    
// 				
// 			
// 			    
// 								
// 			if ( breakpoint_is_relevant_for_AB   or   breakpoint_is_relevant_for_BA )
// 			{
// 
// 			    for (uint lu=0; lu<2; ++lu)
// 			    {
// 				MiniEvent_intersecting_minis_AB->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[mini_A_it->second.MiniID]
// 									.insert(pair_at(profile_brkpt_interval,lu));    
// 									
// 				MiniEvent_intersecting_minis_AB->second.map_MiniRegions_on_MiniEvent_LCRs__to__relevant_profile_breakpoints_for_that_MiniRegion[mini_B_it->second.MiniID]
// 									.insert(pair_at(profile_brkpt_interval,lu));     
// 			    }                                                         
// 																
// 			    
// 			    if ( MiniEvent_intersecting_minis_AB->second.map_MiniRegions_ON_LCRs__to__intersecting_LCR.at(mini_A_it->second.MiniID) == 0 )
// 				MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
// 						    .insert(type_BI__real__real( profile_brkpt_interval,
// 												type_real__real(P_hybrid_AB, P_hybrid_BA)
// 											    ));
// 			    else
// 				MiniEvent_intersecting_minis_AB->second.map_relevant_profile_breakpoint_to_____P_hybrid_AB___P_hybrid_BA
// 						    .insert(type_BI__real__real( profile_brkpt_interval,
// 												type_real__real(P_hybrid_BA, P_hybrid_AB)
// 											    ));    
// 											    
// 
// 											    
// 			}  // relevant for both.
// 
// 		    }//varpos upp
// 		}//varpos low    
// 		
//             }  //  partners for A
//         }  // end for-loop     mini_A_ptr
// 
//     }//found MiniEvent
//     
//     
//     
//     
//     
//     Star_alignment *appropriate_Star;
//     
//     if (best_overall_nohybrid_MR == NULL)
//     {
// 	appropriate_Star =   best_overall_alignment_is_hybrid_AB  ?  &star_alignment_for_hybrid_AB  :  &star_alignment_for_hybrid_BA;		        
//     }
//     else
//     {
// 	for (type_map_uint_to_list_BI__Star_alignment::iterator it_chr = star_alignment__non_hybrids.begin();
// 		it_chr != star_alignment__non_hybrids.end();
// 		++it_chr)
// 	{
// 	    if (it_chr->first != best_overall_nohybrid_MR->chromosome_of_region)
// 		continue;
// 	    
// 	    for (type_list_BI__Star_alignment::iterator  it_st = it_chr->second.begin();
// 		    it_st != it_chr->second.end();	    
// 		    ++it_st)
// 	    {
// 		if ( BOOST_overlap(it_st->first, best_overall_nohybrid_MR->region_interval)  )
// 		{
// 		    appropriate_Star = &(it_st->second);
// 		    break;		
// 		}	    
// 	    }//bi	
// 	}//chr
//     }//MR
// 	
//     
// 
//     
//     appropriate_Star->add_pairwise_alignment_to_progressive_star_alignment
// 			    (   best_overall_probability,
// 				best_PER_and_overall_algmt,
// 				best_hybrid___absolute_A_coords_to_algmt_profile,
// 				best_hybrid___absolute_B_coords_to_algmt_profile,  				
// 				best_orientation,
// 				best_complementarity,
// 				mates,
// 				sum_over_list<real>(P_every_possible_location)   );    
// 
// 				
//     return true;				
//         
// }  //  end of   align_PER_to_MAP_location
// 
// 
// 
// 

































void  set_up_and_align_PERs_using_GPU____OMP
			(type_map_string_to_PER &PERs_for_alignment,
			 const Event *const spawning_Event,
			 const bool &create_and_consider_MiniEvents,
			const bool &consider_GeneConversion_breakpoints,  //irrelevant if "create_and_consider_MiniEvents" = false
			const bool &sum_over_all_unaffected_LCRs,  //"false" for heuristic
			const bool &special_save_breakpoints_actually_captured_in_alignment, //"true" for heuristic
			const bool &do_PER_orient_formatting)  // orthogonal to all other options
{
         
    std::cerr<<" inside  \"set_up_and_align_PERs_using_GPU_and_other_parallelism\"\n";    
    const double time_begin = omp_get_wtime();     
    
    
    const std::vector<type_map_string_to_PER::iterator> PERs___loop(get_loopable_iterators<type_map_string_to_PER>(PERs_for_alignment));
					
    //align on GPU
    CUPA::GPU_alignment_workload  GPU_work;
    PER_bundle_wrapper__space::type_map_string_to_list_PER_bundle_wrapper  record_of_PER_work;
    
    const double time_begin_per_PERs = omp_get_wtime();
    
    
            
    
    
    #pragma omp parallel
    {//parallel
	CUPA::GPU_alignment_workload   GPU_work___thread;
	PER_bundle_wrapper__space::type_map_string_to_list_PER_bundle_wrapper  record_of_PER_work___thread;
	
	
	
	if (do_PER_orient_formatting)
	{
	    #pragma omp for schedule (dynamic,100)	
	    for (uint per_ctr=0; per_ctr < PERs___loop.size(); ++per_ctr)
	    {
		PERs___loop.at(per_ctr)->second.submit_PER_orienting_work_to_GPU(
						    PERs___loop.at(per_ctr)->second.my_MiniRegions.begin()->second,
						     GPU_work___thread, record_of_PER_work___thread);
	    }//per_ctr
	}//do_PER_orient_formatting
	else
	{	
	    #pragma omp for schedule (dynamic,100)	
	    for (uint per_ctr=0; per_ctr < PERs___loop.size(); ++per_ctr)
	    {
		if (create_and_consider_MiniEvents)
		{
		    const type_map_uint_to_uint_to_set_uint MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY(
									PERs___loop.at(per_ctr)->second.create_MiniEvents_and_identify_gender_of_MRs(spawning_Event));

		    PERs___loop.at(per_ctr)->second.make_partners_out_of_MiniRegions_for_each_MiniEvent
									(MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY);
		}//create_and_consider_MiniEvents
		
		
		PERs___loop.at(per_ctr)->second.
		    align_to_every_possible_hybrid_outcome_of_MiniRegions_and_Events_and_record_probabilities_and_relevant_breakpoints(
			    GPU_work___thread, record_of_PER_work___thread,
			    consider_GeneConversion_breakpoints);
				    
		//delete information that will no longer be used.
		PERs___loop.at(per_ctr)->second.clear_all_MiniRegion_sequences__for_memory_considerations();
	    }//per_ctr
	}//alignments
	
	
	
	//combine across threads:
	
	#pragma omp critical (splice_together_GPU_workloads_across_threads__ie___absorb_another_workload)
	{
	    GPU_work.absorb_another_workload(GPU_work___thread);
	}//absorb GPU
	
	#pragma omp critical (splice_together_record_of_GPU_work_submission_across_threads)
	{
	    for (PER_bundle_wrapper__space::type_map_string_to_list_PER_bundle_wrapper::iterator it_rec_thread = record_of_PER_work___thread.begin();
		    it_rec_thread != record_of_PER_work___thread.end();
		    ++it_rec_thread)
	    {
		const PER_bundle_wrapper__space::type_map_string_to_list_PER_bundle_wrapper::iterator insert_rec_it 
		    = record_of_PER_work.insert(std::pair<std::string, PER_bundle_wrapper__space::type_list_PER_bundle_wrapper>
						(it_rec_thread->first, PER_bundle_wrapper__space::type_list_PER_bundle_wrapper())).first;
						
		insert_rec_it->second.splice(insert_rec_it->second.end(), it_rec_thread->second);	    
	    }//record thread
	}//absorb records

	mpfr_free_cache();
    }//parallel
    
    std::cerr << "\n\n\ttime to prepare PERs for GPU workload: " << omp_get_wtime() - time_begin_per_PERs << "  s.\n\n";
    
    
    
    
    

    
    //execute!!!
    const CUPA::type_list_GPU_timing_report  timings_for_GPU_kernels = GPU_work.perform_alignments_on_GPU();


    
    
    
				    
				    
    //post process
    std::cerr << "process_results...\n\\n";    
				    
    const double time_begin_postprocess = omp_get_wtime();
    
    const std::vector<PER_bundle_wrapper__space::type_map_string_to_list_PER_bundle_wrapper::const_iterator> PER_results___loop(
		    get_loopable_constiterators<PER_bundle_wrapper__space::type_map_string_to_list_PER_bundle_wrapper>(record_of_PER_work));
    
    #pragma omp parallel
    {//parallel
    
	if (do_PER_orient_formatting)
	{
	    #pragma omp for schedule (dynamic,100)
	    for (uint per_ctr = 0; per_ctr < PER_results___loop.size(); ++per_ctr)
	    {	    
		PER_bundle_wrapper__space::process_orienting_of_an_illdefined_PER(PER_results___loop.at(per_ctr));
	    }//per_ctr		
	}
	else
	{    
	    #pragma omp for schedule (dynamic,100)
	    for (uint per_ctr = 0; per_ctr < PER_results___loop.size(); ++per_ctr)
	    {
		for (PER_bundle_wrapper__space::type_list_PER_bundle_wrapper::const_iterator it_job = PER_results___loop.at(per_ctr)->second.begin();
			it_job != PER_results___loop.at(per_ctr)->second.end();
			++it_job)
		{		
		    it_job->process_results(special_save_breakpoints_actually_captured_in_alignment);
		}//job
	    }//per_ctr
	}//alignments
        
	mpfr_free_cache();
    }//parallel
    
    std::cerr << "\n\ntime to post process alignment results:  " << omp_get_wtime() - time_begin_postprocess << "  s.\n\n";
    
    
    
//     if (PERs_for_alignment.count("ERR006198.4865334") > 0)
//     {
// 	PERs_for_alignment.at("ERR006198.4865334").print_this_PER();		
//     }
    
    
    
    if (sum_over_all_unaffected_LCRs)
    {	
	const double time_begin_sum_unaff = omp_get_wtime();
	
	#pragma omp parallel
	{//parallel
	
	    #pragma omp for schedule (dynamic,100)
	    for (uint per_ctr = 0; per_ctr < PERs___loop.size(); ++per_ctr)
	    {
		PERs___loop.at(per_ctr)->second.sum_over_all_unaffected_LCRs();
	    }//per_ctr
	    
	    mpfr_free_cache();
	}//parallel
	
	std::cerr << "\n\ntime to sum over all unaffected LCRs:  " << omp_get_wtime() - time_begin_sum_unaff << "  s.\n\n";
	
    }//sum_over_all_unaffected_LCRs
    
    
    
    //print times:
    CUPA::GPU_timing_report::print_total_times_across_multiple_reports(timings_for_GPU_kernels, true, NULL);
			
					    
    std::cerr << "\n\nDONE  with  \"set_up_and_align_PERs_using_GPU____OMP\".\n\n";    

}//set_up_and_align_PERs_using_GPU____OMP