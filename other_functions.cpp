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

#include <limits.h>


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

#include <boost/math/distributions/normal.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>


#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/cauchy.hpp>



#include <gmpfrxx.h>
#include <boost/math/bindings/mpfr.hpp>


#include <mpreal.h>
#include <mpfr.h>


#include <boost/random.hpp>
#include <boost/graph/graph_concepts.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>


#include <api/BamReader.h>
#include <api/BamAux.h>
#include <api/BamAlignment.h>


#include <general_typedefs.h>
#include <translations_through_alignment.h>
#include <Readgroup_statistics.h>

#include "Event.h"
#include "other_functions.h"
#include "globals.h"
#include "MiniRegion.h"
#include "Paired_end_read.h"

#include "Conn_Comp.h"
#include "io_functions.h"
#include "Call_set.h"
#include "Visitation_object.h"

#include "Star_alignment.h"








void remove_degenerate_Events
		    (type_map_uint_to_CC &Conn_Comps)
{
    std::cerr << "\n\nchecking degenerate events...\n\n";
    uint number_of_degenerate_Events = 0;
    uint number_of_CCs_with_degenerate_Events = 0;
    uint number_of_events_removed_related_to_degeneracy = 0;

    type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
    while ( it_cc != Conn_Comps.end() )
    {
	bool contains_degenerate_Event = false;
	for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
		it_ev != it_cc->second.events.end();
		++it_ev)
	{
	    if (BOOST_overlap(it_ev->second.LCRs[0],it_ev->second.LCRs[1]))
	    {
// 		it_ev->second.print_this_entry();
		contains_degenerate_Event = true;	
		++number_of_degenerate_Events;
	    }	    
	}//ev
	
	if (contains_degenerate_Event)
	{
	    ++number_of_CCs_with_degenerate_Events;
	    number_of_events_removed_related_to_degeneracy += it_cc->second.events.size();
	    Conn_Comps.erase(it_cc++);	    
	}
	else
	{  ++it_cc;  }
    }//cc    
    
    std::cerr << "\n\n\nRemoval of degenerate events report:"
	    << "\n\t\tnumber_of_degenerate_Events = " << number_of_degenerate_Events
	    << "\n\t\tnumber_of_CCs_with_degenerate_Events = " << number_of_CCs_with_degenerate_Events
	    << "\n\t\tnumber_of_events_removed_related_to_degeneracy = " << number_of_events_removed_related_to_degeneracy
	    << "\n\n";	    

}//remove_degenerate_Events









void remove_undesirable_Connected_Components_from_set_of_Connected_Components
		    (type_map_uint_to_CC &Conn_Comps,
			const type_set_uint &targeted_CC_IDs,
			const type_set_uint &display_Event_and_breakpoint_range,
			const BOOST_Interval &acceptable_Connected_Component_size_range,
			const BOOST_Interval &acceptable_Connected_Component_compute_range,
			const int &last_completed_CC)

{       
    
//     remove_degenerate_Events(Conn_Comps);
       
    type_map_uint_to_CC::iterator cc_remove_undesirables_it = Conn_Comps.begin();
    while ( cc_remove_undesirables_it != Conn_Comps.end() )                                                   
    {
	if ( !targeted_CC_IDs.empty() )      
	{
	    if (targeted_CC_IDs.count(cc_remove_undesirables_it->first) == 0)
		Conn_Comps.erase(cc_remove_undesirables_it++);   //must be post-increment!
	    else
		++cc_remove_undesirables_it;            
	}
	else if ( !display_Event_and_breakpoint_range.empty() )
	{
	    bool contains_a_display_event = false;
	    for (type_set_uint::const_iterator it_disp = display_Event_and_breakpoint_range.begin();
		    it_disp != display_Event_and_breakpoint_range.end();
		    ++it_disp)
	    {
		if ( cc_remove_undesirables_it->second.events.count(*it_disp) > 0 )
		{
		    contains_a_display_event = true;
		    break;                            
		}
	    }
		
	    if (contains_a_display_event)
		++cc_remove_undesirables_it;
	    else
		Conn_Comps.erase(cc_remove_undesirables_it++);   //must be post-increment!                                                        
	}                                                                
	else
	{                        
	    if ( !BOOST_empty(  acceptable_Connected_Component_size_range  )  )
	    {                                                
		const uint cc_size = cc_remove_undesirables_it->second.events.size();
		if (   ! BOOST_in( cc_size, acceptable_Connected_Component_size_range) )        
		    Conn_Comps.erase(cc_remove_undesirables_it++);   //must be post-increment!
		else
		    ++cc_remove_undesirables_it;
	    }
	    else if (  !BOOST_empty(   acceptable_Connected_Component_compute_range)    )
	    {                 
		const int total_diploid_compute_size = read_total_diploid_compute_size_from_file_for_given_Conn_Comp(  cc_remove_undesirables_it->second  );
		if (  total_diploid_compute_size == -1  or   !BOOST_in(  (uint)total_diploid_compute_size, acceptable_Connected_Component_compute_range ) )        
		    Conn_Comps.erase(cc_remove_undesirables_it++);   //must be post-increment!
		else
		    ++cc_remove_undesirables_it;
	    }                            
	}   
    }//cc
	
	
	
	
    if (last_completed_CC  >=   0)
    {
	if (my_MPI_rank == 0)
	    std::cerr << "\n\nremoving connected components until where we left off...  last_completed_CC = " << last_completed_CC << "\n\n";
	
	type_map_uint_to_CC::iterator cc_it = Conn_Comps.begin();
	while ( cc_it != Conn_Comps.end()   and   cc_it->first != (uint)last_completed_CC )            
	    Conn_Comps.erase(cc_it++);  //post-increment necessary!
	    
	if ( cc_it != Conn_Comps.end() )
	{
	    assert(  Conn_Comps.begin()->first    ==    (uint)last_completed_CC  );
	    Conn_Comps.erase( Conn_Comps.begin() );  // this is the last CC!
	}        
    } // pick up where we left off                    
    
    if (my_MPI_rank == 0)
        std::fprintf(stderr, "\nFinished eliminating undesirable Conn Comps.  Number of remaining Conn Conmps = %u\n", Conn_Comps.size());      
	
} // remove_undesirable_Connected_Components_from_set_of_Connected_Components







    
    
    
    
    
    
    
    
    
    
void count_decent_reads_with_left_endpoints_in_set_region
	(BamTools::BamReader &a_BAM_reader,  //region already set
	const BOOST_Interval &absolute_genome_interval,
	type_map_uint_to_uint &counts_per_position)
{  
    BamTools::BamAlignment temp_core_Balgmt;
    while (a_BAM_reader.GetNextAlignmentCore(temp_core_Balgmt))    
    {
        if ( !temp_core_Balgmt.IsFailedQC())
	{
            if ( temp_core_Balgmt.IsProperPair() )
            {
                if (BOOST_in((uint)(temp_core_Balgmt.Position+1),absolute_genome_interval)  //genome is 1-based, bam-files are 0-based
                     and temp_core_Balgmt.Position < temp_core_Balgmt.MatePosition)  //i.e. is first mate
		{  ++counts_per_position[(uint)(temp_core_Balgmt.Position+1)];  }
            }
            else if (temp_core_Balgmt.IsPaired())    
	    {
                if (temp_core_Balgmt.IsMapped() and  BOOST_in((uint)(temp_core_Balgmt.Position+1), absolute_genome_interval))
                {  ++counts_per_position[(uint)(temp_core_Balgmt.Position+1)];  }
	    }
	}//qc
    }//while
    //Recall that a read is counted only if it STARTS at the considered genome position (i.e. left-endpoint of read lands on that position).  Other reads that overlap are ignored. 

}//count_decent_reads_with_left_endpoints_in_set_region
    
    
    
    
uint count_decent_reads_with_left_endpoints_in_set_region
	(BamTools::BamReader &a_BAM_reader,  //region already set
	const BOOST_Interval &absolute_genome_interval)   
{
    uint Balgmt_ctr = 0;
    BamTools::BamAlignment temp_core_Balgmt;
    while (a_BAM_reader.GetNextAlignmentCore(temp_core_Balgmt))    
    {
        if ( !temp_core_Balgmt.IsFailedQC())
	{
            if ( temp_core_Balgmt.IsProperPair() )
            {
                if (BOOST_in((uint)(temp_core_Balgmt.Position+1),absolute_genome_interval)  //genome is 1-based, bam-files are 0-based
                     and temp_core_Balgmt.Position < temp_core_Balgmt.MatePosition)  //i.e. is first mate
		{  ++Balgmt_ctr;  }
            }
            else if (temp_core_Balgmt.IsPaired())    
	    {
                if (temp_core_Balgmt.IsMapped() and  BOOST_in((uint)(temp_core_Balgmt.Position+1), absolute_genome_interval))
                {  ++Balgmt_ctr;  }
	    }
	}//qc
    }//while
    //Recall that a read is counted only if it STARTS at the considered genome position (i.e. left-endpoint of read lands on that position).  Other reads that overlap are ignored. 
    
    return  Balgmt_ctr;

}//count_decent_reads_with_left_endpoints_in_set_region
    
    
    





// void count_reads_whose_left_endpoint_lies_in_interval
// 			(const uint &chromosome_value,  //NOT index 
// 			const BOOST_Interval &absolute_genome_interval,
// 			BamTools::BamReader &a_BAM_reader,
// 			type_map_uint_to_uint &counts_per_position,
// 			 const bool &use_MPI)
// {
//     if (BOOST_empty(absolute_genome_interval))
//     {
// 	error_message("inside  \"count_reads_whose_left_endpoint_lies_in_interval\",  absolute_genome_interval is empty!!!!", false);
// 	return;    
//     }
//         
//     
// //     //assume  "counts_per_position"  has been pre-allocated  for  "absolute_genome_interval"
// //     if (counts_per_position.count(absolute_genome_interval.lower()) == 0 
// // 	or  counts_per_position.count(absolute_genome_interval.upper()) == 0
// // 	or  counts_per_position.size() < BOOST_width_inclusive(absolute_genome_interval))
// //     {  
// // 	initialize_map_across_interval__via_insert<uint,uint>(counts_per_position, absolute_genome_interval, 0);  
// //     }                
//     const bool success_SetRegion =  a_BAM_reader.SetRegion(  map_chromosome_value_to_BAM_Ref_IDs.at(chromosome_value),
// 							     safe_subtract_base_1(absolute_genome_interval.lower(),1),
// 							    map_chromosome_value_to_BAM_Ref_IDs.at(chromosome_value),
// 							     absolute_genome_interval.upper()-1);                                                                 
// 		    if (!success_SetRegion)                                                         
// 		    {
// 			std::stringstream error_strm;                            
// 			error_strm << "ERROR:  in  \"count_reads_whose_left_endpoint_lies_in_interval\"   SetRegion FAILED!!!!\n\n\t\t\tSetRegion("
// 				<< map_chromosome_value_to_BAM_Ref_IDs.at(chromosome_value) <<"," << absolute_genome_interval<<",\n\t\t\t\t" 
// 				<< map_chromosome_value_to_BAM_Ref_IDs.at(chromosome_value) <<", " << absolute_genome_interval<<")\n\n\n";  
// 			error_message(error_strm, false);
// 			return; 
// 		    }                                                   
//     //NOTE:  BAMTools indexes chromsomes from 0!!!!!.   e.g. chromosome 5 has index ("RefID") 4!!!!
//     //, also, BAM files index the GENOME starting at 0, while SAM files and the Reference start at 1.
//     //Thus we must subtract one from our chromosome value when looking up stuff in the BAM files.    
//     
//     const double time_begin__count = omp_get_wtime();
//     
//     if (!use_MPI)
//     {           
// 	count_decent_reads_with_left_endpoints_in_set_region(a_BAM_reader, absolute_genome_interval, counts_per_position);    
//     }
//     else
//     {		
// 	const type_vector_BI  subintervals(split_interval_into_uniformly_sized_subintervals(absolute_genome_interval,(uint)size_MPI_WORLD, 1000));
// 	for (int sub=my_MPI_rank; sub < subintervals.size(); sub += size_MPI_WORLD)
// 	{
// 	    count_decent_reads_with_left_endpoints_in_set_region(a_BAM_reader, subintervals.at(sub), counts_per_position);	
// 	}
// 	
// 	type_vector_map_uint_to_uint  gathered_counts_per_pos;
// 	boost::mpi::gather<type_map_uint_to_uint>(*world_MPI_boost_ptr, counts_per_position, gathered_counts_per_pos, 0);
// 	
// 	//don't add yourself!!!
// 	for (int rank=1; rank < size_MPI_WORLD; ++rank)
// 	{
// 	    add_map_two_maps<uint,uint>(counts_per_position, gathered_counts_per_pos.at(rank));	    	
// 	}		
//     }//MPI
//     
//     if (my_MPI_rank == 0)
//     {
// 	std::cerr << "\n\ntime for \"count_reads_whose_left_endpoint_lies_in_interval\" = " << (omp_get_wtime() - time_begin__count)/60.00 << "  min.\n\n";
//     }
//         
// }//count_reads_whose_left_endpoint_lies_in_interval









// uint count_reads_whose_left_endpoint_falls_at_position
//                                     ( const uint &chromosome_value,  //NOT index 
//                                       const uint &absolute_genome_pos,
//                                       BamTools::BamReader &a_BAM_reader)                                                     
// {   
//     const bool success_SetRegion =  a_BAM_reader.SetRegion(  map_chromosome_value_to_BAM_Ref_IDs.at(chromosome_value), (int)absolute_genome_pos - 1 - 3,
//                                                              map_chromosome_value_to_BAM_Ref_IDs.at(chromosome_value), (int)absolute_genome_pos - 1 + 3   );                                                                 
//                         if (!success_SetRegion)                                                         
//                         {
//                             std::stringstream error_strm;                            
//                             error_strm << "ERROR:  in  \"count_reads_whose_left_endpoint_falls_at_position\",   SetRegion FAILED!!!!\n\n\t\t\tSetRegion("
//                                     << map_chromosome_value_to_BAM_Ref_IDs.at(chromosome_value) <<"," << (absolute_genome_pos - 1 - 1) <<",\n\t\t\t\t" 
//                                     << map_chromosome_value_to_BAM_Ref_IDs.at(chromosome_value) <<", " << (absolute_genome_pos - 1 + 1) <<")\n\n\n";  
// 			    error_message(error_strm, false);
// 			    return 0; 
//                         }                                                                                          
//         //NOTE:  BAMTools indexes chromsomes from 0!!!!!.   e.g. chromosome 5 has index ("RefID") 4!!!!
//         //, also, BAM files index the GENOME starting at 0, while SAM files and the Reference start at 1.
//         //Thus we must subtract one from our chromosome value when looking up stuff in the BAM files.
//                             
//                             
//     uint Balgmts_ctr = 0;   
//     BamTools::BamAlignment temp_core_Balgmt;
//     
// 
//     while ( a_BAM_reader.GetNextAlignmentCore(temp_core_Balgmt) )    
//     {
// //         if ( !temp_core_Balgmt.IsFailedQC())
//             if ( temp_core_Balgmt.IsProperPair() )
//             {
//                 if ( temp_core_Balgmt.Position  == (int)absolute_genome_pos - 1  //genome is 1-based, bam-files are 0-based
//                      and temp_core_Balgmt.Position < temp_core_Balgmt.MatePosition )  //i.e. is first mate
//                     ++Balgmts_ctr;                 
//             }
//             else if ( temp_core_Balgmt.IsPaired() )    
// 	    {
//                 if (  temp_core_Balgmt.IsMapped() and  temp_core_Balgmt.Position  ==  (int)absolute_genome_pos - 1)
//                     ++Balgmts_ctr;                         
// 	    }
//     }
//     //Recall that a read is counted only if it STARTS at the considered genome position (i.e. left-endpoint of read lands on that position).  Other reads that overlap are ignored.   
//         
//     
//     return Balgmts_ctr;
// 
// }   //  end of count_reads_whose_left_endpoint_falls_at_position


























//This way is for Mini_regions
type_uint__uint create_list_of_unformatted_Paired_end_reads_from_list_of_Balgmts
                              ( type_list_Balgmts &Balgmts_from_a_Mini_region,
                                BamTools::BamReader &my_BAM_reader,
				type_list_PER &created_PERs)
{
//WARNING: this function erases "Balgmts_from_a_Mini_region".  It should not be used upon return!

    uint number_of_discarded_mates = 0;
        
    type_list_Balgmts::iterator it_begin = Balgmts_from_a_Mini_region.begin();    //always start at begin since deleting pairs of Balgmts will invalidate all other iterators.
    
    
    //go through and find pairs of Balgmts (have same name), add them to  "created_PERs", and then delete them from  "Balgmts_from_a_Mini_region"
    while ( !Balgmts_from_a_Mini_region.empty() )
    {                
        bool found_mate = false;                            
        
        //search through remaining alignments to find mate  (mate for first algmt, i.e. at .begin)
        type_list_Balgmts::iterator it_mate = ++Balgmts_from_a_Mini_region.begin();
            
        while (it_mate != Balgmts_from_a_Mini_region.end())   
	{
            if (    it_mate->Name.compare(it_begin->Name) == 0  )
                // &&  it_mate->AlignmentFlag != it_begin->AlignmentFlag )   //unnecessary.
            {
                found_mate = true;
                break;
            }     
	    ++it_mate;
	}
    


    
        if ( !found_mate )  
        {
            if ( it_begin->IsProperPair() )
            {
                BamTools::BamAlignment Balgmt;
                const bool set_sucess 
                     =  my_BAM_reader.SetRegion( (int)it_begin->MateRefID, (int)it_begin->MatePosition       -        3,      
                                                (int)it_begin->MateRefID, (int)it_begin->MatePosition        +        3       );      
           
                if (!set_sucess)                                
                {
                    std::stringstream warning_strm;
                    warning_strm << "\n\nWARNING!    unable to set region when looking for extra mate (proper pair):  it_begin->MateRefID  =  "  << it_begin->MateRefID  
                            << " ,     it_begin->MatePosition  ="  <<  it_begin->MatePosition  << "   +- 5"
                            << "        for   read:   "  << it_begin->Name <<  "\n\n";
                    warning_message(warning_strm, false);              
                }
                else
                {
                    while( my_BAM_reader.GetNextAlignmentCore(Balgmt) )
		    {
                        if (  Balgmt.Position == it_begin->MatePosition )
                        {
                            Balgmt.BuildCharData();                          
                            if ( Balgmt.Name.compare(it_begin->Name) == 0 )
                            {
                                it_mate = Balgmts_from_a_Mini_region.insert( Balgmts_from_a_Mini_region.end(),  Balgmt);    
                                found_mate = true;
                                break;                                                                            
                            }                                                                    
                        } 
		    }
                }
            }//is proper pair             
            else   //is not proper pair.
            {
                if ( it_begin->IsMateMapped() )
                {
                    BamTools::BamAlignment Balgmt;
                    const bool set_sucess 
                            =  my_BAM_reader.SetRegion(  (int)it_begin->MateRefID, ((int)it_begin->MatePosition)      -         3,        // - 5 !!!!!!!!!
                                                        (int)it_begin->MateRefID, ((int)it_begin->MatePosition)      +         3       );  // + 5 !!!!!!!!!            
                    
                    if (set_sucess)  
                    {
                        while( my_BAM_reader.GetNextAlignmentCore(Balgmt) ) 
			{
                            if (  Balgmt.Position == it_begin->MatePosition )
                            {                                                
                                Balgmt.BuildCharData();                                
                                if ( Balgmt.Name.compare(it_begin->Name) == 0 )
                                {
                                    it_mate = Balgmts_from_a_Mini_region.insert( Balgmts_from_a_Mini_region.end(), Balgmt);    
                                    found_mate = true;
                                    break;                                                                            
                                }                                                                    
                            }  
			}
                    }
                    else
                    {
			
                        std::stringstream warning_strm;                        
                        warning_strm << "\n\nWARNING!    unable to set region when looking for extra mate:  it_begin->MateRefID  =  "  << it_begin->MateRefID  
                                << " ,     it_begin->MatePosition  ="  <<  it_begin->MatePosition  << "   +- 5"
                                << "        for   read:   "  << it_begin->Name <<  "\n\n"
                                <<  "my_BAM_reader.IsOpen()  =   "  << my_BAM_reader.IsOpen()
                                << "\n\nmy_BAM_reader.GetFilename() =  " << my_BAM_reader.GetFilename()
                                << "\n\nmy_BAM_reader.GetErrorString() = "  <<    my_BAM_reader.GetErrorString() << "\n\n";
			print_Balgmt(*it_begin, &warning_strm);
			warning_message(warning_strm, false);                  
                    }                                
                }
                else                    
                {                    
                    std::stringstream error_strm;                    
                    error_strm << "\n\nERROR!  Should never reach the case  \"!.isProperPair()\"  but \".isPaired()\"  but  \"!.IsMateMapped() \"\n\n";
                    print_Balgmt(*it_begin, &error_strm);
		    error_message(error_strm, false);
                    
//                     exit(1);
//                     BamTools::BamAlignment Balgmt;
//                     my_BAM_reader.SetRegion(it_begin->MateRefID ,it_begin->Position      -        3,        // - 5 !!!!!!!!!
//                                             it_begin->MateRefID ,it_begin->Position     +         3       );  // + 5 !!!!!!!!!
                                
//                     while( my_BAM_reader.GetNextAlignmentCore(Balgmt) )                
//                         if (  Balgmt.Position == it_begin->Position    &&    ( it_begin->IsMapped() != Balgmt.IsMapped() )    )
//                         {                                                
//                             Balgmt.BuildCharData();
//                             
//                             if ( Balgmt.Name.compare(it_begin->Name) == 0 )
//                             {
//                                 it_mate = Balgmts_from_a_Mini_region.insert( Balgmts_from_a_Mini_region.end(),
//                                                                             Balgmt);    
//                                 found_mate = true;
//                                 break;                                                                            
//                             }                                                                    
//                         }
                }   
            }//not proper pair

        }//if !found_mate
            
        
        if (found_mate)       
        {    
            created_PERs.push_back(   Paired_end_read( *it_begin, *it_mate )    );  //thus, the first mate is, by default, guaranteed to be in the search region.
//             else                
//                 created_PERs.push_back(   Paired_end_read( *it_mate, *it_begin )    );
            Balgmts_from_a_Mini_region.erase( it_mate );
            it_begin = Balgmts_from_a_Mini_region.erase( it_begin );
        }
        else
        {
	    std::stringstream error_strm;
            error_strm << "\n\n\n\nERROR:  unable to find mate, even after searching.  DISCARDING   read   "  << it_begin->Name  << "\n\n";	    
            print_Balgmt(*it_begin, &error_strm);
	    std::cerr << error_strm.str() << "\n\n";
	    
            it_begin = Balgmts_from_a_Mini_region.erase(it_begin);
	    ++number_of_discarded_mates;
        }    
    }  // end while     !Balgmts_from_a_Mini_region.empty()
    
    
    
    type_list_PER::iterator it_check_errors = created_PERs.begin();
    uint number_failures = 0;
    while (it_check_errors != created_PERs.end())
    {
	if (it_check_errors->failed_PER_flag)
	{
	    it_check_errors = created_PERs.erase(it_check_errors);
	    ++number_failures;
	}
	else
	{  ++it_check_errors;  }
    }        
    
    return type_uint__uint(number_failures,number_of_discarded_mates);
  
}  // end of   create_list_of_unformatted_Paired_end_reads_from_listof_Balgmts

































BOOST_Interval get_padded_region_below_and_above_being_wary_of_chromosome_centromere_and_endpoints
                                  (const uint &chromosome_of_position_to_be_padded,
                                   const uint &position_to_be_padded_below_and_above,
                                   const uint &pad_below_amount,
                                   const uint &pad_above_amount)
{

    //assume p-arm:
    uint extreme_lower_bound = 1;
    uint extreme_upper_bound = Event::centromere_coordinates.at(chromosome_of_position_to_be_padded).first;


    if ( position_to_be_padded_below_and_above > Event::centromere_coordinates.at(chromosome_of_position_to_be_padded).second ) //q-arm
    {
        extreme_lower_bound = Event::centromere_coordinates.at(chromosome_of_position_to_be_padded).second;  //top of centromere  
        extreme_upper_bound = Event::chromosome_lengths.at(chromosome_of_position_to_be_padded);  //top of chromosome
    }
    else if ( position_to_be_padded_below_and_above > Event::centromere_coordinates.at(chromosome_of_position_to_be_padded).first )  //inside centromere
    {        
        std::fprintf(stderr, "You asked to pad a genome position that is in the centromere!!!\n");
        std::fprintf(stderr, "chromosome_of_position_to_be_padded = %u, position_to_be_padded_below_and_above = %u\npad_below_amount = %u, pad_above_amount = %u\nin function  \"get_padded_region_below_and_above_being_wary_of_chromosome_centromere_and_endpoints\"\n\n\n",
                        chromosome_of_position_to_be_padded,
                        position_to_be_padded_below_and_above,
                        pad_below_amount,
                        pad_above_amount);
        exit(1);
    }
  
    if (  (int)position_to_be_padded_below_and_above - (int)pad_below_amount   <   (int)extreme_lower_bound  
            or     (int)position_to_be_padded_below_and_above + (int)pad_above_amount   >  (int)extreme_upper_bound )
    {
        std::stringstream error_strm;
        print_line_of_markers("WARNING! ", &error_strm);
        print_line_of_markers("(", &error_strm);
        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";      
        
        error_strm << "\n\nERROR!    position_to_be_padded_below_and_above =  "  << position_to_be_padded_below_and_above  << " does not fit cleanly in mapped parts of chromosome in  \"get_padded_region_below_and_above_being_wary_of_chromosome_centromere_and_endpoints\".\n\n";
        error_strm << "\tpad_below_amount  =  " << pad_below_amount  <<  "\n\tpad_above_amount  =  " << pad_above_amount << "\n\n";
        error_strm << "\tchromosome_of_position_to_be_padded = " << chromosome_of_position_to_be_padded << " extreme_lower_bound = " << extreme_lower_bound
                    << "    extreme_upper_bound =  " << extreme_upper_bound 
                    << "Event::centromere_coordinates.at(chromosome_of_position_to_be_padded)  =  <" <<
                    Event::centromere_coordinates.at(chromosome_of_position_to_be_padded).first
                    << ",   " << Event::centromere_coordinates.at(chromosome_of_position_to_be_padded).second  <<  ">\n\n";
        error_strm  <<  "Event::chromosome_lengths.at(chromosome_of_position_to_be_padded)  =  "  << Event::chromosome_lengths.at(chromosome_of_position_to_be_padded);     
                    
        
        print_line_of_markers(")", &error_strm);
        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
    }
      
  
  
    return BOOST_Interval(  (uint)std::max<int>((int)position_to_be_padded_below_and_above - (int)pad_below_amount, (int)extreme_lower_bound ),
                            (uint)std::min<int>((int)position_to_be_padded_below_and_above + (int)pad_above_amount, (int)extreme_upper_bound)     );  

}  // end of   get_padded_region_below_and_above_being_wary_of_chromosome_centromere_and_endpoints






















void adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
                        (const type_map_uint_to_uint &hybrid_indeces_to_algmt_profile,
                         const BOOST_Interval &hybrid_breakpoints,//&index_of_hybrid_AB_that_is_first_B_nucleotide_IE_the_breakpoint_varpos,
                         const uint &absolute_genome_position_of_A_mapping_to_index_0_of_hybrid,  //Note, index 0 of hybrid is NOT  necessarily the first index in the algmt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         const bool &A_mini_was_reversed,
                         const uint &absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid,
                         const bool &B_mini_was_reversed,
			 const bool &is_NAHR,
                         type_map_uint_to_uint &absolute_A_coords_to_algmt_profile,
                         type_map_uint_to_uint &absolute_B_coords_to_algmt_profile)
{
    assert(BOOST_is_a_point(hybrid_breakpoints));
    
    adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
                        ( hybrid_indeces_to_algmt_profile,
                          hybrid_breakpoints,
                         type_uint__uint(absolute_genome_position_of_A_mapping_to_index_0_of_hybrid,
					 absolute_genome_position_of_A_mapping_to_index_0_of_hybrid),
                          A_mini_was_reversed,
                          absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid,
                          B_mini_was_reversed,
			  is_NAHR,
                         absolute_A_coords_to_algmt_profile,
                         absolute_B_coords_to_algmt_profile);         
}//adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates                        








//"A" and "B" do NOT necessarily match "0" and "1".
//Thus, NAHR is always of the form "AB", and GeneConverison always of the form "ABA". 
//So if you want "BA or "BAB", then you have to "trick" this function by swapping the coordinates.
void adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
                        (const type_map_uint_to_uint &hybrid_indeces_to_algmt_profile,
                         const BOOST_Interval &hybrid_breakpoints,//&index_of_hybrid_AB_that_is_first_B_nucleotide_IE_the_breakpoint_varpos,
                         const type_uint__uint &absolute_genome_position_of_A_mapping_to_index_0_of_hybrid___and___to_second_GeneConv_breakpoint_on_hybrid,  //Note, index 0 of hybrid is NOT  necessarily the first index in the algmt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         const bool &A_mini_was_reversed,
                         const uint &absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid,
                         const bool &B_mini_was_reversed,
			 const bool &is_NAHR,
                         type_map_uint_to_uint &absolute_A_coords_to_algmt_profile,
                         type_map_uint_to_uint &absolute_B_coords_to_algmt_profile)                        
{
                     
    absolute_A_coords_to_algmt_profile.clear();
    absolute_B_coords_to_algmt_profile.clear();
      //absolute genome to position    --->     algmt_coords of PER and hybrid.   
      
    if (is_NAHR)
    {//NAHR
	//"hybrid_breakpoints" should be a point.
	{//A	
	    const type_map_uint_to_uint  map_hybrid_to_algmt__on_A(
		    get_map_intersection_of_map_keys_and_interval<uint,uint>(hybrid_indeces_to_algmt_profile, 
									    BOOST_Interval(0, safely_subtract(hybrid_breakpoints.lower(),1))));

	    type_map_uint_to_uint::iterator insert_it_adjusted_coords = absolute_A_coords_to_algmt_profile.begin();
	    for (type_map_uint_to_uint::const_iterator it_coords = map_hybrid_to_algmt__on_A.begin();
		    it_coords != map_hybrid_to_algmt__on_A.end();
		    ++it_coords)
	    {
		const uint adjusted_key  =  A_mini_was_reversed    ?
		    absolute_genome_position_of_A_mapping_to_index_0_of_hybrid___and___to_second_GeneConv_breakpoint_on_hybrid.first - it_coords->first
		    : absolute_genome_position_of_A_mapping_to_index_0_of_hybrid___and___to_second_GeneConv_breakpoint_on_hybrid.first + it_coords->first;            
		
		insert_it_adjusted_coords = absolute_A_coords_to_algmt_profile.insert(insert_it_adjusted_coords, type_uint__uint(adjusted_key, it_coords->second));                                                                
	    }		
	}//A
	
	{//B	    		    
	    const type_map_uint_to_uint  map_hybrid_to_algmt__on_B(
		    get_map_intersection_of_map_keys_and_interval<uint,uint>(hybrid_indeces_to_algmt_profile,
									    BOOST_Interval(hybrid_breakpoints.upper(), std::numeric_limits<int>::max())));	
									    
	    type_map_uint_to_uint::iterator insert_it_adjusted_coords = absolute_B_coords_to_algmt_profile.begin();
	    for (type_map_uint_to_uint::const_iterator it_coords = map_hybrid_to_algmt__on_B.begin();
		    it_coords != map_hybrid_to_algmt__on_B.end();
		    ++it_coords)
	    {		
		const uint adjusted_key = B_mini_was_reversed  ? 
					absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid - (it_coords->first - hybrid_breakpoints.upper())
				    :   absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid + (it_coords->first - hybrid_breakpoints.upper());
					
		insert_it_adjusted_coords = absolute_B_coords_to_algmt_profile.insert(insert_it_adjusted_coords, type_uint__uint(adjusted_key, it_coords->second));                                                                
	    }   									    
	}//B        
    }//NAHR
    else
    {//GeneConv
	{//first A
	    const type_map_uint_to_uint  map_hybrid_to_algmt__on_A(
		    get_map_intersection_of_map_keys_and_interval<uint,uint>(hybrid_indeces_to_algmt_profile, 
									    BOOST_Interval(0, safely_subtract(hybrid_breakpoints.lower(),1) )));

	    type_map_uint_to_uint::iterator insert_it_adjusted_coords = absolute_A_coords_to_algmt_profile.begin();
	    for (type_map_uint_to_uint::const_iterator it_coords = map_hybrid_to_algmt__on_A.begin();
		    it_coords != map_hybrid_to_algmt__on_A.end();
		    ++it_coords)
	    {
		const uint adjusted_key  =  A_mini_was_reversed    ?
		    absolute_genome_position_of_A_mapping_to_index_0_of_hybrid___and___to_second_GeneConv_breakpoint_on_hybrid.first - it_coords->first
		    : absolute_genome_position_of_A_mapping_to_index_0_of_hybrid___and___to_second_GeneConv_breakpoint_on_hybrid.first + it_coords->first;            
		
		insert_it_adjusted_coords = absolute_A_coords_to_algmt_profile.insert(insert_it_adjusted_coords, type_uint__uint(adjusted_key, it_coords->second));                                                                
	    }		
	}//first A
	
	//gene conversion tract
	{//B	    		    
	    const type_map_uint_to_uint  map_hybrid_to_algmt__on_B(
		    get_map_intersection_of_map_keys_and_interval<uint,uint>(
							hybrid_indeces_to_algmt_profile,
							BOOST_Interval(hybrid_breakpoints.lower(), safely_subtract(hybrid_breakpoints.upper(),1))));
									    
	    type_map_uint_to_uint::iterator insert_it_adjusted_coords = absolute_B_coords_to_algmt_profile.begin();
	    for (type_map_uint_to_uint::const_iterator it_coords = map_hybrid_to_algmt__on_B.begin();
		    it_coords != map_hybrid_to_algmt__on_B.end();
		    ++it_coords)
	    {		
		const uint adjusted_key = B_mini_was_reversed  ? 
					absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid - (it_coords->first - hybrid_breakpoints.lower())
				    :   absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid + (it_coords->first - hybrid_breakpoints.lower());
					
		insert_it_adjusted_coords = absolute_B_coords_to_algmt_profile.insert(insert_it_adjusted_coords, type_uint__uint(adjusted_key, it_coords->second));                                                                
	    }   									    
	}//B    
		
	{//second A
	    const type_map_uint_to_uint  map_hybrid_to_algmt__on_A2(
		    get_map_intersection_of_map_keys_and_interval<uint,uint>(hybrid_indeces_to_algmt_profile, 
									     BOOST_Interval(hybrid_breakpoints.upper(), std::numeric_limits<int>::max())));	
									    
	    type_map_uint_to_uint::iterator insert_it_adjusted_coords = absolute_A_coords_to_algmt_profile.begin();
	    for (type_map_uint_to_uint::const_iterator it_coords = map_hybrid_to_algmt__on_A2.begin();
		    it_coords != map_hybrid_to_algmt__on_A2.end();
		    ++it_coords)
	    {
		const uint adjusted_key  =  A_mini_was_reversed    ?
				absolute_genome_position_of_A_mapping_to_index_0_of_hybrid___and___to_second_GeneConv_breakpoint_on_hybrid.second
										- (it_coords->first - hybrid_breakpoints.upper())
				: absolute_genome_position_of_A_mapping_to_index_0_of_hybrid___and___to_second_GeneConv_breakpoint_on_hybrid.second
										+ (it_coords->first - hybrid_breakpoints.upper());            
		
		insert_it_adjusted_coords = absolute_A_coords_to_algmt_profile.insert(insert_it_adjusted_coords,
											type_uint__uint(adjusted_key, it_coords->second));                                                                
	    }
	}//second A        
    }//GeneConv            

}//adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates

























//"A" and "B" are abstract - they do not necesssariyl map as "A" --> "LCR 0",  "B" --> "LCR 1".   (it may be reversed!!!)
//This function always creates hybrids of the form "AB" or "ABA"
//Thus, if you want a "BA" or "BAB" hybrid, then you have to "trick" the function by swapping the coordinates.
//"breakpoints_on_A/B_mini_region_absolute" mark the actual breakpoint on each sequence.  Inclusion of the breakpoint depends, of course on NAHR vs. GeneConv.
void  get_variational_indeces_mapped_onto_hybrid_sequence
			(const MiniRegion &A_mini,
			const BOOST_Interval &breakpoints_on_A_mini_region_absolute,  //not to be included!!
			const bool &A_mini_must_be_reversed,   //  0 = lower() is boundary,  1= upper() is boundary
			const MiniRegion &B_mini,
			const BOOST_Interval &breakpoints_on_B_mini_region_absolute, //this one is included!!
			const bool &B_mini_must_be_reversed,     //  0 = lower() is boundary,  1= upper() is boundary
			 const type_set_uint &absolute_variational_positions__on_A,
			 const type_set_uint &absolute_variational_positions__on_B,
			 const bool &is_NAHR,
			 type_set_uint &hybrid_varpos_indeces__A,
			 type_set_uint &hybrid_varpos_indeces__B)
{

    hybrid_varpos_indeces__A.clear();
    hybrid_varpos_indeces__B.clear();
    
    //Gene Conversion:
    if (!is_NAHR)
    {	
	//always of the form "ABA".
	
	const uint length_of_first_A_tract =  A_mini_must_be_reversed  ?
						A_mini.region_interval.upper() - breakpoints_on_A_mini_region_absolute.upper()  //not inclusive
						:
						breakpoints_on_A_mini_region_absolute.lower() - A_mini.region_interval.lower();  //not inclusive
	
	const uint length_of_B_tract =   BOOST_width_inclusive(breakpoints_on_B_mini_region_absolute) - 1; //only includes one of its two endpoints, whether reversed or not.
	
	
	{//first A
	    if (!A_mini_must_be_reversed)
	    {
		type_set_uint::iterator insert_index_A_it = hybrid_varpos_indeces__A.begin();
		
		for (type_set_uint::const_iterator it_varpos_A = absolute_variational_positions__on_A.begin();
			it_varpos_A != absolute_variational_positions__on_A.end();
			++it_varpos_A)
		{
		    if (*it_varpos_A >= breakpoints_on_A_mini_region_absolute.lower())
		    {  break;  }
		    else if (*it_varpos_A >= A_mini.region_interval.lower())
		    {
			insert_index_A_it = hybrid_varpos_indeces__A.insert(insert_index_A_it, *it_varpos_A - A_mini.region_interval.lower()); 
		    }		
		}
	    }//DO NOT reverse
	    else
	    {
		type_set_uint::iterator insert_index_A_it = hybrid_varpos_indeces__A.begin();
		
		for (type_set_uint::const_reverse_iterator it_varpos_A = absolute_variational_positions__on_A.rbegin();
			it_varpos_A != absolute_variational_positions__on_A.rend();
			++it_varpos_A)
		{
		    if (*it_varpos_A <= breakpoints_on_A_mini_region_absolute.upper())
		    {  break;  }
		    else if (*it_varpos_A <= A_mini.region_interval.upper())
		    {
			insert_index_A_it = hybrid_varpos_indeces__A.insert(insert_index_A_it, A_mini.region_interval.upper() - (*it_varpos_A)); 
		    }
		}
	    }//reverse
	}//first A
	
	{//B tract
	
	    if (!B_mini_must_be_reversed)
	    {
		//lower is included, upper excluded
		type_set_uint::iterator insert_index_B_it = hybrid_varpos_indeces__B.begin();
		
		for (type_set_uint::const_iterator it_varpos_B = absolute_variational_positions__on_B.begin();
			it_varpos_B != absolute_variational_positions__on_B.end();
			++it_varpos_B)
		{
		    if (*it_varpos_B >= breakpoints_on_B_mini_region_absolute.upper())
		    {  break;  }
		    else if (*it_varpos_B >= breakpoints_on_B_mini_region_absolute.lower())
		    {
			insert_index_B_it = hybrid_varpos_indeces__B.insert(insert_index_B_it, *it_varpos_B - breakpoints_on_B_mini_region_absolute.lower() + length_of_first_A_tract);
		    }		    
		}//varpos
	    }//DO NOT reverse
	    else
	    {//reverse
		//lower excluded, upper included
		type_set_uint::iterator insert_index_B_it = hybrid_varpos_indeces__B.begin();
		
		for (type_set_uint::const_iterator it_varpos_B = absolute_variational_positions__on_B.begin();
			it_varpos_B != absolute_variational_positions__on_B.end();
			++it_varpos_B)
		{
		    if (*it_varpos_B > breakpoints_on_B_mini_region_absolute.upper())
		    {  break;  }
		    else if (*it_varpos_B > breakpoints_on_B_mini_region_absolute.lower())
		    {
			insert_index_B_it = hybrid_varpos_indeces__B.insert(insert_index_B_it, breakpoints_on_B_mini_region_absolute.upper() - (*it_varpos_B) + length_of_first_A_tract);
		    }    
		}//varpos
	    }//reverse
	}//B tract
	
	{//second A
	    if (!A_mini_must_be_reversed)
	    {
		type_set_uint::iterator insert_index_A_it = hybrid_varpos_indeces__A.begin();
		
		for (type_set_uint::const_reverse_iterator it_varpos_A = absolute_variational_positions__on_A.rbegin();
			it_varpos_A != absolute_variational_positions__on_A.rend();
			++it_varpos_A)
		{
		    if (*it_varpos_A < breakpoints_on_A_mini_region_absolute.upper())
		    {  break;  }
		    else if (*it_varpos_A <= A_mini.region_interval.upper())
		    {
			insert_index_A_it = hybrid_varpos_indeces__A.insert(insert_index_A_it, 
									    *it_varpos_A - breakpoints_on_A_mini_region_absolute.upper() + length_of_first_A_tract + length_of_B_tract); 
		    }		
		}
	    }//DO NOT reverse
	    else
	    {
		type_set_uint::iterator insert_index_A_it = hybrid_varpos_indeces__A.begin();
		
		for (type_set_uint::const_iterator it_varpos_A = absolute_variational_positions__on_A.begin();
			it_varpos_A != absolute_variational_positions__on_A.end();
			++it_varpos_A)
		{
		    if (*it_varpos_A > breakpoints_on_A_mini_region_absolute.lower())
		    {  break;  }
		    else if (*it_varpos_A >= A_mini.region_interval.lower())
		    {
			insert_index_A_it = hybrid_varpos_indeces__A.insert(insert_index_A_it,
									    breakpoints_on_A_mini_region_absolute.lower() - (*it_varpos_A) + length_of_first_A_tract + length_of_B_tract);
		    }
		}
	    }//reverse
	}//second A
	
    }//gene conversion       
    else
    {//NAHR
	//always of the form "AB"
	
	const uint length_of_A_sequence =  A_mini_must_be_reversed  ?  
						    A_mini.region_interval.upper() - breakpoints_on_A_mini_region_absolute.upper()  //not inclusive
						    :
						    breakpoints_on_A_mini_region_absolute.lower() - A_mini.region_interval.lower();  //not inclusive
						    
	{//A
	    if (!A_mini_must_be_reversed)
	    {
		type_set_uint::iterator insert_index_A_it = hybrid_varpos_indeces__A.begin();
		
		for (type_set_uint::const_iterator it_varpos_A = absolute_variational_positions__on_A.begin();
			it_varpos_A != absolute_variational_positions__on_A.end();
			++it_varpos_A)
		{
		    if (*it_varpos_A >= breakpoints_on_A_mini_region_absolute.lower())
		    {  break;  }
		    else if (*it_varpos_A >= A_mini.region_interval.lower())
		    {
			insert_index_A_it = hybrid_varpos_indeces__A.insert(insert_index_A_it, *it_varpos_A - A_mini.region_interval.lower()); 
		    }		
		}
	    }//DO NOT reverse
	    else
	    {
		type_set_uint::iterator insert_index_A_it = hybrid_varpos_indeces__A.begin();
		
		for (type_set_uint::const_reverse_iterator it_varpos_A = absolute_variational_positions__on_A.rbegin();
			it_varpos_A != absolute_variational_positions__on_A.rend();
			++it_varpos_A)
		{
		    if (*it_varpos_A <= breakpoints_on_A_mini_region_absolute.upper())
		    {  break;  }
		    else if (*it_varpos_A <= A_mini.region_interval.upper())
		    {
			insert_index_A_it = hybrid_varpos_indeces__A.insert(insert_index_A_it, A_mini.region_interval.upper() - (*it_varpos_A)); 
		    }
		}
	    }//reverse
	}//A
	
	{//B
	    if (!B_mini_must_be_reversed)
	    {   
		for (type_set_uint::const_reverse_iterator it_varpos_B = absolute_variational_positions__on_B.rbegin();
			it_varpos_B != absolute_variational_positions__on_B.rend();
			++it_varpos_B)
		{
		    if (*it_varpos_B < breakpoints_on_B_mini_region_absolute.upper())
		    {  break;  }
		    else if (*it_varpos_B <= B_mini.region_interval.upper())
		    {
			hybrid_varpos_indeces__B.insert(*it_varpos_B - breakpoints_on_B_mini_region_absolute.upper() + length_of_A_sequence);
		    }
		}
	    }//DO NOT reverse
	    else
	    {    
		for (type_set_uint::const_iterator it_varpos_B = absolute_variational_positions__on_B.begin();
			it_varpos_B != absolute_variational_positions__on_B.end();
			++it_varpos_B)
		{
		    if (*it_varpos_B > breakpoints_on_B_mini_region_absolute.lower())
		    {  break;  }
		    else if (*it_varpos_B >= B_mini.region_interval.lower())
		    {
			hybrid_varpos_indeces__B.insert(breakpoints_on_B_mini_region_absolute.lower() - (*it_varpos_B) + length_of_A_sequence); 
		    }
		}
	    }//reverse
	}//B
	
    }//NAHR
        
}//get_variational_indeces_mapped_onto_hybrid_sequence























//"A" and "B" are abstract - they do not necesssariyl map as "A" --> "LCR 0",  "B" --> "LCR 1".   (it may be reversed!!!)
//This function always creates hybrids of the form "AB" or "ABA"
//Thus, if you want a "BA" or "BAB" hybrid, then you have to "trick" the function by swapping the coordinates.
//"breakpoints_on_A/B_mini_region_absolute" mark the actual breakpoint on each sequence.  Inclusion of the breakpoint depends, of course on NAHR vs. GeneConv.
type_string__BI  get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
                                    (const MiniRegion &A_mini,
                                     const BOOST_Interval &breakpoints_on_A_mini_region_absolute,  //not to be included!!
                                     const bool &must_complement_mini_A_to_get_to_profile,
                                     const bool &A_mini_must_be_reversed,   //  0 = lower() is boundary,  1= upper() is boundary
                                     const MiniRegion &B_mini,
                                     const BOOST_Interval &breakpoints_on_B_mini_region_absolute, //this one is included!!
                                     const bool &must_complement_mini_B_to_get_to_profile,                                     
                                     const bool &B_mini_must_be_reversed, //0 = lower() is boundary,  1= upper() is boundary
				     const bool &is_NAHR) //else, GeneConv
{
    
    type_string__BI  hybrid_seq___breakpoint;    
    hybrid_seq___breakpoint.first.reserve(std::min<uint>(A_mini.region_sequence.size(), B_mini.region_sequence.size()));
    
    if (is_NAHR)
    {//NAHR
	//always of the form "AB"
    
	{//A
	    std::string subseq_A;
	
	    if (!A_mini_must_be_reversed)
	    {
		const int length_of_subseq = (int)breakpoints_on_A_mini_region_absolute.lower() - A_mini.region_interval.lower();  //not inclusive
		subseq_A = A_mini.region_sequence.substr(0, length_of_subseq);
	    }//DO NOT reverse
	    else
	    {//reverse
		const int start_index = breakpoints_on_A_mini_region_absolute.upper() - A_mini.region_interval.lower() + 1; // not inclusive
		subseq_A = get_reverse_of_sequence(A_mini.region_sequence.substr(start_index));
	    }//reverse
	    
	    if (must_complement_mini_A_to_get_to_profile)
	    {
		subseq_A = get_complement_of_sequence(subseq_A);
	    }
	    
	    hybrid_seq___breakpoint.first = subseq_A;
	    hybrid_seq___breakpoint.second = BOOST_make_point(subseq_A.size());
	}//A
	
	{//B
	    std::string subseq_B;
	    
	    if (!B_mini_must_be_reversed)
	    {
		const int start_index = (int)breakpoints_on_B_mini_region_absolute.upper() - B_mini.region_interval.lower();  // inclusive
		subseq_B = B_mini.region_sequence.substr(start_index);
	    }//DO NOT reverse
	    else
	    {//reverse
		const int length_of_subseq = (int)breakpoints_on_B_mini_region_absolute.lower() - B_mini.region_interval.lower() + 1;//inclusive
		subseq_B = get_reverse_of_sequence(B_mini.region_sequence.substr(0, length_of_subseq));		
	    }//reverse
	    
	    if (must_complement_mini_B_to_get_to_profile)
	    {
		subseq_B = get_complement_of_sequence(subseq_B);
	    }		    
	    
	    hybrid_seq___breakpoint.first.append(subseq_B);
	}//B
    
    }//NAHR
    else
    {//GeneConv
	//always of the form  "ABA"
	
	uint first_breakpoint_index_on_hybrid;
	{//first A
	    std::string subseq_A;
	    
	    if (!A_mini_must_be_reversed)
	    {
		const int length_of_subseq = (int)breakpoints_on_A_mini_region_absolute.lower() - A_mini.region_interval.lower(); // not inclusive
		subseq_A = A_mini.region_sequence.substr(0, length_of_subseq);
	    }//DO NOT reverse
	    else
	    {//reverse
		const int start_index = (int)breakpoints_on_A_mini_region_absolute.upper() - A_mini.region_interval.lower() + 1; // not inclusive
		subseq_A = get_reverse_of_sequence(A_mini.region_sequence.substr(start_index));
	    }//reverse
	    
	    if (must_complement_mini_A_to_get_to_profile)
	    {
		subseq_A = get_complement_of_sequence(subseq_A);
	    }
	    
	    hybrid_seq___breakpoint.first = subseq_A;
	    first_breakpoint_index_on_hybrid = subseq_A.size();
	    //hybrid_seq___breakpoint.second.set(subseq_A.size(), std::numeric_limits<int>::max());
	}//first A
	
	{//B tract
	    std::string subseq_B;
	    
	    if (!B_mini_must_be_reversed)
	    {
		const int start_index = (int)breakpoints_on_B_mini_region_absolute.lower() - B_mini.region_interval.lower(); // inclusive
		const int length_of_tract = (int)breakpoints_on_B_mini_region_absolute.upper() - breakpoints_on_B_mini_region_absolute.lower(); // not inclusive
		subseq_B = B_mini.region_sequence.substr(start_index, length_of_tract);
	    }//DO NOT reverse
	    else
	    {
		const int start_index = (int)breakpoints_on_B_mini_region_absolute.lower() - B_mini.region_interval.lower() + 1; // not inclsuive
		const int length_of_tract = (int)breakpoints_on_B_mini_region_absolute.upper() - breakpoints_on_B_mini_region_absolute.lower(); // not inclusive
		subseq_B = get_reverse_of_sequence(B_mini.region_sequence.substr(start_index, length_of_tract));
	    }//reverse
	    
	    if (must_complement_mini_B_to_get_to_profile)
	    {
		subseq_B = get_complement_of_sequence(subseq_B);
	    }
	    
	    hybrid_seq___breakpoint.first.append(subseq_B);
	    hybrid_seq___breakpoint.second.set(first_breakpoint_index_on_hybrid, hybrid_seq___breakpoint.first.size());
	}//B tract
	
	{//second A
	    std::string subseq_A;
	    
	    if (!A_mini_must_be_reversed)
	    {
		const int start_index = (int)breakpoints_on_A_mini_region_absolute.upper() - A_mini.region_interval.lower(); // inclusive
		subseq_A = A_mini.region_sequence.substr(start_index);
	    }//DO NOT reverse
	    else
	    {//reverse
		const int length_of_subseq = (int)breakpoints_on_A_mini_region_absolute.lower() - A_mini.region_interval.lower() + 1;  // inclusive
		subseq_A = get_reverse_of_sequence(A_mini.region_sequence.substr(0, length_of_subseq));
	    }//reverse
	    
	    if (must_complement_mini_A_to_get_to_profile)
	    {
		subseq_A = get_complement_of_sequence(subseq_A);
	    }
	    
	    hybrid_seq___breakpoint.first.append(subseq_A);
	}//second A            
    }//GeneConv
        
    
    return  hybrid_seq___breakpoint;        
    
}//get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint




















//"0" and "1" stand for the TRUE LCR sides.
type_vector_int  label_absolute_coordinates_on_hybrid
		    (const MiniRegion &mini_LCR_0,
		     const BOOST_Interval &absolute_breakpoint_region__LCR_0,
		     const MiniRegion &mini_LCR_1,
		     const BOOST_Interval &absolute_breakpoint_region__LCR_1,
		     const bool &is_NAHR,  //else, GeneConv
		     const bool &AB_ABA__false_____BA_BAB__true,
		     const type_recomb_class  &recomb_type_of_Event)
{

    type_vector_int hybrid_indeces_labeled_as_absolute_coords;
    hybrid_indeces_labeled_as_absolute_coords.reserve(maximum_average_fragment_length);
    
    if (!AB_ABA__false_____BA_BAB__true)    
    {
	for (int j = (int)mini_LCR_0.region_interval.lower(); j < (int)absolute_breakpoint_region__LCR_0.lower(); ++j)
	{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
	
	if (is_NAHR)
	{//NAHR
	    if (recomb_type_of_Event == recomb_class__DupDel)
	    {
		for (int j = (int)absolute_breakpoint_region__LCR_1.upper(); j <= (int)mini_LCR_1.region_interval.upper();   ++j)
		{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }		    
	    }
	    else
	    {
		for (int j = (int)absolute_breakpoint_region__LCR_1.lower(); j >= (int)mini_LCR_1.region_interval.lower();   --j)
		{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }						    
	    }
	}//NAHR
	else
	{//GeneConv
	    if (recomb_type_of_Event == recomb_class__DupDel)
	    {
		for (int j = (int)absolute_breakpoint_region__LCR_1.lower(); j < (int)absolute_breakpoint_region__LCR_1.upper();   ++j)
		{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
		
		for (int j = (int)absolute_breakpoint_region__LCR_0.upper(); j <= (int)mini_LCR_0.region_interval.upper();   ++j)
		{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
	    }
	    else
	    {
		for (int j = (int)absolute_breakpoint_region__LCR_1.upper(); j > (int)absolute_breakpoint_region__LCR_1.lower();   --j)
		{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
		
		for (int j = (int)absolute_breakpoint_region__LCR_0.upper(); j <= (int)mini_LCR_0.region_interval.upper();   ++j)
		{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
	    }		    		
	}//Geneconv
    }//AB__ABA
    else
    {//BA__BAB
    	
	if (recomb_type_of_Event == recomb_class__DupDel)
	{
	    for (int j = (int)mini_LCR_1.region_interval.lower(); j < (int)absolute_breakpoint_region__LCR_1.lower(); ++j)
	    {  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
	}
	else
	{
	    for (int j = (int)mini_LCR_1.region_interval.upper(); j > (int)absolute_breakpoint_region__LCR_1.upper(); --j)
	    {  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
	}
	
		
	if (is_NAHR)
	{//NAHR		    
	    for (int j = (int)absolute_breakpoint_region__LCR_0.upper(); j <= (int)mini_LCR_0.region_interval.upper();   ++j)
	    {  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }		    
	}//NAHR
	else
	{//GeneConv
			
	    for (int j = (int)absolute_breakpoint_region__LCR_0.lower(); j < (int)absolute_breakpoint_region__LCR_0.upper();   ++j)
	    {  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
	
	
	    if (recomb_type_of_Event == recomb_class__DupDel)
	    {
		for (int j = (int)absolute_breakpoint_region__LCR_1.upper(); j <= (int)mini_LCR_1.region_interval.upper(); ++j)
		{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
	    }
	    else
	    {
		for (int j = (int)absolute_breakpoint_region__LCR_1.lower(); j >= (int)mini_LCR_1.region_interval.lower(); --j)
		{  hybrid_indeces_labeled_as_absolute_coords.push_back(j);  }
	    }		    		
	}//Geneconv
    }//BA__BAB


    return  hybrid_indeces_labeled_as_absolute_coords;

}//label_absolute_coordinates_on_hybrid











































void set_map_from_chromosome_values_to_RefID_in_BAM_file
                                    (BamTools::BamReader &a_BAM_reader)
{

    BamTools::RefVector the_refs = a_BAM_reader.GetReferenceData();  

    char *name_of_chromo = new char[option_size];

			std::cerr << "the_refs.size() = " << the_refs.size() << "\n";
		   	for (uint k=0; k<the_refs.size(); ++k)
			{		
				std:: cerr << "\t\t[" << the_refs.at(k).RefName.c_str() << "]\n\n";		
			}

	//determine naming convention
	int chrnametype = 0;
	bool found_chr_name = false;
	while (!found_chr_name and chrnametype < 4)
	{	
	        if (chrnametype == 0)
	            {  std::sprintf(name_of_chromo, "chr1");  }
	        else if (chrnametype == 1)
	            {  std::sprintf(name_of_chromo, "1");  }
		else
		{
			std::cerr << "\n\n\nERROR!  cannot find chrnametype in setting bam ref id!!!!!\n\n\n";
			std::cerr << "the_refs.size() = " << the_refs.size() << "\n";
		   	for (uint k=0; k<the_refs.size(); ++k)
			{		
				std:: cerr << "\t\t[" << the_refs.at(k).RefName.c_str() << "]\n\n";		
			}
			exit(1);
		}

	        for (uint k=0; k<the_refs.size(); ++k)
		{
	            if ( strcasecmp(name_of_chromo, the_refs.at(k).RefName.c_str()) == 0 )
	            {
			found_chr_name = true;
	                break;
	            }            
		}
		if (!found_chr_name)
		{	++chrnametype;  }
	}//while chrnametype

std::cerr << "determined   chrnametype = " << chrnametype << "\n\n";
    
    // for chromos 1,...,22
    for (uint j=1; j<=22; ++j)
    {

	        if (chrnametype == 0)
	            std::sprintf(name_of_chromo, "chr%u", j);
	        else if (chrnametype == 1)
	            std::sprintf(name_of_chromo, "%u", j);
		else
		{
			std::cerr << "\n\n\nERROR!  cannot find chrnametype in setting bam ref id!!!!!\n\n\n";
			exit(1);
		}
	        
	        for (uint k=0; k<the_refs.size(); ++k)
		{
	            if ( strcasecmp(name_of_chromo, the_refs.at(k).RefName.c_str()) == 0 )
	            {
	                map_chromosome_value_to_BAM_Ref_IDs.insert( type_uint__uint(j,k) );
	                break;
	            }            
		}
    }
    
    
    // for X:
    if (chrnametype == 0)
        std::sprintf(name_of_chromo, "chrX");
    else if (chrnametype == 1)
        std::sprintf(name_of_chromo, "X");
		else
		{
			std::cerr << "\n\n\nERROR!  cannot find chrnametype in setting bam ref id!!!!!\n\n\n";
			exit(1);
		}
    
    for (uint k=0; k<the_refs.size(); ++k)
    {
        if ( strcasecmp(name_of_chromo, the_refs.at(k).RefName.c_str()) == 0 )
        {            
            map_chromosome_value_to_BAM_Ref_IDs.insert( type_uint__uint(23,k) );
            break;
        }   
    }
        
        
    // for Y:
    if (chrnametype == 0)
        std::sprintf(name_of_chromo, "chrY");
    else if (chrnametype == 1)
        std::sprintf(name_of_chromo, "Y");   
 		else
		{
			std::cerr << "\n\n\nERROR!  cannot find chrnametype in setting bam ref id!!!!!\n\n\n";
			exit(1);
		}
    
    for (uint k=0; k<the_refs.size(); ++k)
    {
        if ( strcasecmp(name_of_chromo, the_refs.at(k).RefName.c_str()) == 0 )
        {
            map_chromosome_value_to_BAM_Ref_IDs.insert( type_uint__uint(24,k) );
            break;
        }           
    }
        
    delete[] name_of_chromo;
        
} // set_map_from_chromosome_values_to_RefID_in_BAM_file
























































































type_map_uint_to_uint get_relative_hybrid_coords_to_algmt_profile
                        (const type_string__string &PER_and_hybrid_algmt,
                         const BOOST_Interval &endpoints_of_hybrid_on_algmt)
{
    if (BOOST_empty(endpoints_of_hybrid_on_algmt))
    {
        std::stringstream error_strm;        
        error_strm << "\n\nERROR!   in  \"get_relative_hybrid_coords_to_algmt_profile()\":" <<
                      "endpoints_of_hybrid_on_algmt is empt!!!!!!\n\n\n";                      
        error_strm << "PER_and_hybrid_algmt:\n\n";        
        error_strm << PER_and_hybrid_algmt.first << "\n" << PER_and_hybrid_algmt.second << "\n\n";        
        error_message(error_strm, false);          
    }        
    
    type_map_uint_to_uint Hybrid_indeces_to_algmt_profile;
    type_map_uint_to_uint::iterator insert_hybrid_it = Hybrid_indeces_to_algmt_profile.begin();        
    
    
    
    uint hybrid_indx = endpoints_of_hybrid_on_algmt.lower();
    for (uint algmt_indx = 0; algmt_indx < PER_and_hybrid_algmt.second.size(); ++algmt_indx)
    {
        if (PER_and_hybrid_algmt.second.at(algmt_indx) != '-')
        {
            insert_hybrid_it = Hybrid_indeces_to_algmt_profile.insert(insert_hybrid_it, type_uint__uint(hybrid_indx, algmt_indx));
            ++hybrid_indx;                                                
        }
    }
    
    
    return Hybrid_indeces_to_algmt_profile;            

}  // get_relative_hybrid_coords_to_algmt_profile
























































                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
bool check_that_relevant_regions_are_well_defined__and__hide_undefined_regions_if_necessary
			(const type_map_uint_to_Event::iterator  &event_being_summed_out,
			 type_map_uint_to_Event &mutable_events_of_this_CC)
{
    
//check that this region of the genome is well-defined - it does not contain unsequenced/unmapped/unlocalized regions full of "N".
// if_not__hide_small_region_around_undefined_regions_and_remove_from_natural_poisson_intervals                      

    
    if (my_MPI_rank == 0)
	std::cerr << "\n\n\ncheck that this region of the genome is well-defined - it does not contain unsequenced/unmapped/unlocalized regions full of \"N\"...\n\n\n";
    
    
    type_map_uint_to_set_uint  map_chr_to_pos_of_N_nucleotides;
		    
    for (type_map_uint_to_list_BI__string::const_iterator it_check = event_being_summed_out->second.uploaded_pieces_of_chromos.begin();
	    it_check != event_being_summed_out->second.uploaded_pieces_of_chromos.end();
	    ++it_check)  
    {
	for (type_list__BI__string::const_iterator it_reg = it_check->second.begin();
		it_reg != it_check->second.end();
		++it_reg)	
	{
	    for (uint j = 0; j < it_reg->second.size(); ++j)
	    {
		if (it_reg->second.at(j) == 'N'   or   it_reg->second.at(j) == 'n')                            
		    map_chr_to_pos_of_N_nucleotides[it_check->first].insert(  it_reg->first.lower() + j  );    
	    }	
	}
    }
				
		
		
    uint number_of_N = 0;
    for (type_map_uint_to_set_uint::const_iterator it_chr = map_chr_to_pos_of_N_nucleotides.begin();
	    it_chr != map_chr_to_pos_of_N_nucleotides.end();
	    ++it_chr)
    {
	number_of_N += it_chr->second.size();
    }
    
    
    
    if (number_of_N > 50)
    {
	if (my_MPI_rank == 0)
	{
	    std::stringstream warning_strm;
	    print_line_of_markers("WARNING! ", &warning_strm);
	    print_line_of_markers("(", &warning_strm);
	    warning_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
	    
	    warning_strm << "\n\nWARNING!     number_of_N  = "  << number_of_N 
			    <<  "  for event "  << event_being_summed_out->second.UID 
			    << "   of CC " << event_being_summed_out->second.my_Conn_Comp->Conn_Comp_ID << "\n\n";
	    
	    warning_strm << "\n\nscanning for regions of high concentrations of \"N\"  nucleotides in reference for this Event...\n\n";                
			    
	    print_line_of_markers(")", &warning_strm);
	    std::fprintf(stderr, "\n\n%s\n\n", warning_strm.str().c_str() );                                                                                       
	}
	
	
	const type_map_uint_to_list_BI regions_near_N_nucleotides( 
			scan_region_for_unidentified__N__nucleotides_and_merge_for_ignoring(  
					    map_chr_to_pos_of_N_nucleotides, 20, 50));
						 
					    
	if (!regions_near_N_nucleotides.empty())
	{                 
	    if (my_MPI_rank == 0)
	    {
		std::cerr << "\n\nERASING regions of high concentrations of \"N\"  nucleotides in reference...\n\n";                        
		
		print_map_to_list<uint, BOOST_Interval>(regions_near_N_nucleotides, "regions_near_N_nucleotides", NULL, true);                            
	    }
	    
	    
	    
	    type_map_uint_to_uint_to_uint_to_longdouble::iterator it_nat = event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.begin();
	    while (it_nat != event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.end())
	    {
		std::vector<type_map_uint_to_uint_to_longdouble::iterator> nat_chromos_with_ignored_regions;
		
		for (type_map_uint_to_uint_to_longdouble::iterator it_nat_chr = it_nat->second.begin();
			it_nat_chr != it_nat->second.end();
			++it_nat_chr)                            
		{
		    if (regions_near_N_nucleotides.count(it_nat_chr->first) > 0)
		    {  nat_chromos_with_ignored_regions.push_back(it_nat_chr);  }
		}
					    
		const uint number_relevant_chromos = nat_chromos_with_ignored_regions.size();                           
		#pragma omp parallel
		{
		    #pragma omp for schedule(dynamic,1)
		    for (uint nat_chr_ctr = 0; nat_chr_ctr < number_relevant_chromos; ++nat_chr_ctr)
		    {                                 
			for (type_list_BI::const_iterator it_bi_ignore 
					= regions_near_N_nucleotides.at(nat_chromos_with_ignored_regions.at(nat_chr_ctr)->first).begin();
				it_bi_ignore != regions_near_N_nucleotides.at(nat_chromos_with_ignored_regions.at(nat_chr_ctr)->first).end();
				++it_bi_ignore)
			{                                    
			    type_map_uint_to_longdouble::iterator it_nat_pos = nat_chromos_with_ignored_regions.at(nat_chr_ctr)->second.begin();
			    while (it_nat_pos  !=  nat_chromos_with_ignored_regions.at(nat_chr_ctr)->second.end())
			    {
				if (BOOST_in(it_nat_pos->first, *it_bi_ignore))
				    nat_chromos_with_ignored_regions.at(nat_chr_ctr)->second.erase(it_nat_pos++);
				else
				    ++it_nat_pos;
			    }
			}                                
		    }                                
					    
		    mpfr_free_cache();
		}//parallel  
						
		type_map_uint_to_uint_to_longdouble::iterator it_chr = it_nat->second.begin();
		while (it_chr != it_nat->second.end())
		{
		    if (it_chr->second.empty())
		    {  it_nat->second.erase(it_chr++);  }
		    else
		    {  ++it_chr;  }
		}
		
		if (it_nat->second.empty())
		{  event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.erase(it_nat++);  }
		else
		{  ++it_nat;  }
	    }//it_nat      
	    
	    	    
	    
		if (my_MPI_rank == 0)
		{
		    std::cerr << "\n\nnumber of natural poisson intervals remaining after erasing \"N\" positions = " 
				<< event_being_summed_out->second.natural_poisson_intervals___gender_INsensitive.size() << "\n\n";                        
		}
	    
		
	    
	    
	    
	    
	    const uint number_neighbors_and_self = event_being_summed_out->second.local_interlocking_Events.size() 
					    + event_being_summed_out->second.global_interlocking_Events.size() + 1;
					    
	    std::vector< type_map_uint_to_Event::iterator > neighbors_and_self;
	    neighbors_and_self.reserve(number_neighbors_and_self);
			{//scope                                    
			    neighbors_and_self.push_back( event_being_summed_out );
			    
			    for (type_map_uint_to_Interlock::const_iterator it_neighb = event_being_summed_out->second.local_interlocking_Events.begin();
				    it_neighb != event_being_summed_out->second.local_interlocking_Events.end();
				    ++it_neighb)
			    {  neighbors_and_self.push_back(mutable_events_of_this_CC.find(it_neighb->second.entry->UID));  }  
			    
			    for (type_map_uint_to_Interlock::const_iterator it_neighb = event_being_summed_out->second.global_interlocking_Events.begin();
				    it_neighb != event_being_summed_out->second.global_interlocking_Events.end();
				    ++it_neighb)
			    {  neighbors_and_self.push_back(mutable_events_of_this_CC.find(it_neighb->second.entry->UID));  }
			}//scope
				    
	    #pragma omp parallel
	    {
		#pragma omp for schedule(dynamic,1)                        
		for (uint ev_ctr = 0; ev_ctr < number_neighbors_and_self; ++ev_ctr)    
		{  neighbors_and_self.at(ev_ctr)->second.ignore_regions_of_the_genome_from_consideration_by_this_Event(regions_near_N_nucleotides);  }
		
		mpfr_free_cache();
	    }
	    
	    if (my_MPI_rank == 0)
	    {
		std::cerr << "\n\nDone erasing regions of high concentrations of \"N\"  nucleotides in reference...\n\n";
		save_Events_affected_by_Undefined_regions_to_file(mutable_events_of_this_CC, regions_near_N_nucleotides);
	    }//0
	    
	    return true;
	}// !regions_near_N_nucleotides.empty()                                                                                                                                                                                                  
	else
	{                                        
	    if (my_MPI_rank == 0)
		std::cerr << "\n\nnot enough \"N\"  nucleotides in reference for this event to worry about it.  nothing changed...\n\n";                    
	}
    }// N > 50       
    
    return false;
	    
}//check_that_relevant_regions_are_well_defined__and__hide_undefined_regions_if_necessary

            
                                
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
void eliminate_GeneConversion_event_outcomes
	    (std::vector<type_vector_Sparse_map::const_iterator>  &iterators_to_haploid_state_vectors)
{
    std::vector<type_vector_Sparse_map::const_iterator>  NAHR_only__state_vectors;
    NAHR_only__state_vectors.reserve(iterators_to_haploid_state_vectors.size());   
    
    for (uint j=0; j < iterators_to_haploid_state_vectors.size(); ++j)
    {
	const type_map_uint_to_uint map_outcome_to_Event(get_inverse_of_map(iterators_to_haploid_state_vectors.at(j)->sparse_map_UID_to_value));
	if (map_outcome_to_Event.count(hap_outcome__GeneConv_ABA) == 0  and  map_outcome_to_Event.count(hap_outcome__GeneConv_BAB) == 0)
	{  NAHR_only__state_vectors.push_back(iterators_to_haploid_state_vectors.at(j));  }
    }//j            
    
    iterators_to_haploid_state_vectors = NAHR_only__state_vectors;
//     iterators_to_haploid_state_vectors.shrink_to_fit();

}//eliminate_GeneConversion_event_outcomes

                    
                    
                    
                    
                    
                    
void  make_gender_sensitive_iterators_for_state_vectors
					(std::vector< type_vector_Sparse_map::const_iterator >  &iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors,
					 std::vector< type_vector_Sparse_map::const_iterator >  &iterators_to_haploid_1_X_or_Y_sensitive_state_vectors,
					 const type_map_uint_to_Event &events_of_some_CC,
					 const type_list_Visitation_object::const_iterator &it_Visit)
{

    iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.clear();
    iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.clear();
    
    type_set_uint  events_that_occurr_on_XXXX_chromosome;
    type_set_uint  events_that_occurr_on_YYYY_chromosome;
    
    for (type_map_uint_to_Event::const_iterator it_ev = events_of_some_CC.begin();
	    it_ev != events_of_some_CC.end();
	    ++it_ev)
    {
	if (it_ev->second.chromos[0] == 23)
	{  events_that_occurr_on_XXXX_chromosome.insert(it_ev->first);  }
	else if (it_ev->second.chromos[0] == 24)
	{  events_that_occurr_on_YYYY_chromosome.insert(it_ev->first);  }
    }
    
    
    const uint haploid_state_vec_size = it_Visit->permissible_haploid_state_vectors_sparse.size();
    iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.reserve(haploid_state_vec_size);
    iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.reserve(haploid_state_vec_size);        
    
    
            
    for (type_vector_Sparse_map::const_iterator  it_state_vec =  it_Visit->permissible_haploid_state_vectors_sparse.begin();
	    it_state_vec != it_Visit->permissible_haploid_state_vectors_sparse.end();
	    ++it_state_vec)
    {
	bool contains_an_XXXX_event = false;
	bool contains_a_YYYY_event = false;
    
	const type_set_uint occurring_events_in_this_state_vec(  extract_keys_of_map_and_return_as_set<uint,uint>( it_state_vec->sparse_map_UID_to_value )  );
	
	for (type_set_uint::const_iterator it_occ_ev = occurring_events_in_this_state_vec.begin();
		it_occ_ev != occurring_events_in_this_state_vec.end();
		++it_occ_ev)
	{
	    if (events_that_occurr_on_XXXX_chromosome.count(*it_occ_ev) > 0)                    
	    {  contains_an_XXXX_event = true;  }          
	    if (events_that_occurr_on_YYYY_chromosome.count(*it_occ_ev) > 0)
	    {  contains_a_YYYY_event = true;  }
	}
	
	
	if (contains_an_XXXX_event  and  contains_a_YYYY_event)
	{  continue;  }
	else if (contains_an_XXXX_event)
	{
	    iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.push_back(it_state_vec);
	    if (gender_of_individual_is_female)
	    {  iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.push_back(it_state_vec);  }
	}
	else if (contains_a_YYYY_event)
	{
	    if (!gender_of_individual_is_female)
	    {  iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.push_back(it_state_vec);  }               
	}
	else // has nothing to do with X or Y
	{
	    iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.push_back(it_state_vec);
	    iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.push_back(it_state_vec);                
	}                                                    
    }//it_state_vecs   
    
    
//     iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors.shrink_to_fit();
//     iterators_to_haploid_1_X_or_Y_sensitive_state_vectors.shrink_to_fit();    
            
}//make gender sensitive iterators for states                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
		
                            
                
                    
        







































                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                
                
                
                
                
                
                
                    
void  write_line_of__Call_vs_Val__of_table__to_ostream
                (   const Call_set *const &ptr_to_call,
                    const Validated_rearrangement *const &ptr_to_val,
                    std::ostream &some_ostream)
{

    if (ptr_to_val != NULL)
    {
        std::string val_outcome_str("none");
        if (ptr_to_val->recomb_type == 1)
            val_outcome_str = "inversion";
        else //dupdel
	{
            if (ptr_to_val->outcome == 1)
                val_outcome_str = "deletion";
            else
                val_outcome_str = "duplication";
	}


        some_ostream    << ptr_to_val->variant_accession_name << "\t"
                        << ptr_to_val->NA_subject << "\t"
                        << "chr" << ptr_to_val->chromosome << "\t"
                        << ptr_to_val->breakpoints.lower() << "\t"
                        << ptr_to_val->breakpoints.upper() << "\t"
                        << ptr_to_val->published_study << "\t"
                        << val_outcome_str << "\t";
    }
    else
    {
        some_ostream << " \t"
                        " \t"
                        " \t"
                        " \t"
                        " \t"
                        " \t"
                        " \t" ;      
    }
   

    if (ptr_to_call != NULL)
    {
                                        
        const std::string call_outcome_str  =    ptr_to_call->potentially_gene_conversion    ?    
                                                    "gene conversion"   :   ptr_to_call->var_outcome;
         
        if (   BOOST_empty(ptr_to_call->brkpts)   )    
            some_ostream    << "-\t-\t";                                   
        else                                
            some_ostream    << ptr_to_call->brkpts.lower() << "\t"
                            << ptr_to_call->brkpts.upper() << "\t"
                            << call_outcome_str << "\t"
                            << ptr_to_call->cc_of_event << "\t"
                            << ptr_to_call->event_of_var << "\n";
    }
    else
    {
        some_ostream  <<  " \t"
                    << " \t"
                    << " \t"
                    << " \t"
                    << " \n";                
    }


} //write_line_of__Call_vs_Val__of_table__to_ostream
                    






















































type_map_uint_to_list_BI  scan_region_for_unidentified__N__nucleotides_and_merge_for_ignoring
                                        (const type_map_uint_to_set_uint &map_chromo_to_pos,
                                         const int &min_allowable_gap_size_between_consecutive_N_nucleotides,
                                         const int &min_size_for_ignoring)
{
    
    type_map_uint_to_list_BI  regions_to_ignore;
    
    
    for (type_map_uint_to_set_uint::const_iterator it_chr = map_chromo_to_pos.begin();
            it_chr != map_chromo_to_pos.end();
            ++it_chr)
    {        
        int beginning_of_this_N_region = (int)(*(it_chr->second.begin()));
        int last_pos = beginning_of_this_N_region;
        
        for (type_set_uint::const_iterator it_pos = ++it_chr->second.begin();
                it_pos != it_chr->second.end();
                ++it_pos)        
	{
            if ( (int)*it_pos - last_pos  >  (int)min_allowable_gap_size_between_consecutive_N_nucleotides )
            {
                if (last_pos - beginning_of_this_N_region + 1  >=  min_size_for_ignoring )
                {
                    const uint lb_for_new_BI = (uint)std::max<int>(1, ((int)beginning_of_this_N_region) 
                                                                                    - ((int)maximum_average_fragment_length) 
                                                                                    - maximum_frag_length_absolute_dev);
                    const uint ub_for_new_BI =  (uint)last_pos + (uint)maximum_average_fragment_length + (uint)maximum_frag_length_absolute_dev;                    
                            if (lb_for_new_BI > ub_for_new_BI)
                            {
                                std::stringstream error_strm;                                
                                error_strm << "\n\nERROR!   in   \"scan_region_for_unidentified__N__nucleotides_and_merge_for_ignoring\"   somehow    lb_for_new_BI  =  "
                                            << lb_for_new_BI << "   >   " 
                                            <<  ub_for_new_BI  << "  =  ub_for_new_BI  !!!!"
                                            << "\n\t\tlast_pos = " << last_pos 
                                            << ",  beginning_of_this_N_region = " <<  beginning_of_this_N_region 
                                            << ",   min_size_for_ignoring = " << min_size_for_ignoring 
                                            << ",   *it_pos = " <<  *it_pos 
                                            << ",    it_chr->first   =  "  << it_chr->first
                                            << ",   maximum_average_fragment_length = "  <<  maximum_average_fragment_length
                                            << ",   maximum_frag_length_absolute_dev = " << maximum_frag_length_absolute_dev
                                            << "\n\n\n";
                                error_message(error_strm,false);                                 
                            }
                    
                    regions_to_ignore[it_chr->first].push_back(BOOST_Interval(lb_for_new_BI, ub_for_new_BI));    
                }
                
                beginning_of_this_N_region = (int)*it_pos;
                last_pos =  (int)*it_pos;
            }
            else
	    {  last_pos = (int)*it_pos;  }
	}
                
        //last time
        if (last_pos - beginning_of_this_N_region + 1  >=  min_size_for_ignoring)
        {
            const uint lb_for_new_BI = (uint)std::max<int>(1, ((int)beginning_of_this_N_region) 
                                                                            - ((int)maximum_average_fragment_length) 
                                                                            - maximum_frag_length_absolute_dev);
            const uint ub_for_new_BI =  (uint)last_pos + (uint)maximum_average_fragment_length + (uint)maximum_frag_length_absolute_dev;            
                    if (lb_for_new_BI > ub_for_new_BI)
                    {
                        std::stringstream error_strm;                        
                        error_strm << "\n\nERROR!   in   \"scan_region_for_unidentified__N__nucleotides_and_merge_for_ignoring\"  (last time)   somehow    lb_for_new_BI  =  "
                                    << lb_for_new_BI << "   >   " 
                                    <<  ub_for_new_BI  << "  =  ub_for_new_BI  !!!!"
                                    << "\n\t\tlast_pos = " << last_pos  
                                    << ",  beginning_of_this_N_region = " <<  beginning_of_this_N_region 
                                    << ",   min_size_for_ignoring = " << min_size_for_ignoring 
                                    << ",    it_chr->first   =  "  << it_chr->first
                                    << ",   maximum_average_fragment_length = "  <<  maximum_average_fragment_length
                                    << ",   maximum_frag_length_absolute_dev = " << maximum_frag_length_absolute_dev
                                    << "\n\n\n";
                        
                        error_message(error_strm,false);                                      
                    }            
                        
            regions_to_ignore[it_chr->first].push_back(BOOST_Interval(lb_for_new_BI, ub_for_new_BI));  
        }
    }
    
    
    
    merge_intervals_on_each_chromosome(regions_to_ignore);
        
            

    return regions_to_ignore;
    
} // scan_region_for_unidentified__N__nucleotides_and_merge_for_ignoring




































void make_MAP_calls_from_marginal_posterior_distributions_per_Event
                                      (const type_map_uint_to_Marginal_Event_posterior_distribution  &marginal_distributions,
                                       const Conn_Comp &the_Conn_Comp,
                                       std::stringstream &output_marginal_Positive_calls_strm,
                                        std::stringstream &output_marginal_Negative_calls_strm,
                                        std::stringstream &output_marginal_all_calls_strm)
{
    
    std::cerr << "\nmake_MAP_calls_from_marginal_posterior_distributions_per_Event\n";         
    
    for (type_map_uint_to_Event::const_iterator it_ev = the_Conn_Comp.events.begin();
	    it_ev != the_Conn_Comp.events.end();
	    ++it_ev)
    {
	if (it_ev->second.compressed_map_LCR_to_profile__for_each_LCR.size() != 2)
	{
	    std::stringstream error_strm;
	    print_line_of_markers("ERROR! ", &error_strm);
	    print_line_of_markers("(", &error_strm);
	    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
	    
	    error_strm << "\n\nERROR!   ev"  <<  it_ev->second.UID   << "   has empty compressed_map_LCR_to_profile__for_each_LCR  in  \"make_MAP_calls_from_marginal_posterior_distributions_per_Event\"!!!!!!!!!! \n\n";
	    error_strm << "\n\nit_ev->second.compressed_map_LCR_to_profile__for_each_LCR.size() = " << it_ev->second.compressed_map_LCR_to_profile__for_each_LCR.size() << "\n\n";
	    
	    print_line_of_markers(")", &error_strm);
	    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );	
	}
    }
    
    
    
    
    std::cerr << "\n\t\t\tget MAP calls\n";    
    type_map_uint_to_Sampled_diploid_Event_data  diploid_marginal_MAP_calls;
            
    
    for (type_map_uint_to_Marginal_Event_posterior_distribution::const_iterator it_ev = marginal_distributions.begin();
            it_ev != marginal_distributions.end();
            ++it_ev)
    {                                
        const type_map_uint_to_Event::const_iterator occurring_Event = the_Conn_Comp.events.find( it_ev->first );
	assert( occurring_Event !=  the_Conn_Comp.events.end() );
        
        
        type_map__hap_outcome__hap_outcome__to__longdouble::const_iterator it_best_marginal_outcome__this_event = it_ev->second.map_diploid_outcomes_to_posterior_probability.begin();        
        longdouble outcome_normalization_constant = 0;
		//necessary because of the haploid neutral stuff.
	
        for (type_map__hap_outcome__hap_outcome__to__longdouble::const_iterator it_outcome = it_ev->second.map_diploid_outcomes_to_posterior_probability.begin();
                it_outcome != it_ev->second.map_diploid_outcomes_to_posterior_probability.end();
                ++it_outcome)
        {
            outcome_normalization_constant += it_outcome->second;            
            if (it_outcome->second  >  it_best_marginal_outcome__this_event->second)
	    {  it_best_marginal_outcome__this_event = it_outcome;  }
        }                                
        
        
	diploid_marginal_MAP_calls[it_ev->first].associated_genome = genome_name;          
        diploid_marginal_MAP_calls[it_ev->first].cc_of_event = it_ev->second.CC_ID;
	diploid_marginal_MAP_calls[it_ev->first].event_UID = it_ev->second.ev_UID;
	diploid_marginal_MAP_calls[it_ev->first].the_diploid_Event_outcome = it_best_marginal_outcome__this_event->first;
	diploid_marginal_MAP_calls[it_ev->first].P_diploid_outcome =  it_best_marginal_outcome__this_event->second / outcome_normalization_constant;
        
		
	type_map__BI__BI__to__longdouble::const_iterator it_best_brkpts__this_event 
		    = it_ev->second.map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability.at(it_best_marginal_outcome__this_event->first).begin();
    
	longdouble brkpts_normalization_constant = 0;
	for (type_map__BI__BI__to__longdouble::const_iterator it_brkpts 
			    = it_ev->second.map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability.at(it_best_marginal_outcome__this_event->first).begin();
		it_brkpts != it_ev->second.map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability.at(it_best_marginal_outcome__this_event->first).end();
		++it_brkpts)
	{
	    brkpts_normalization_constant += it_brkpts->second;
	    if ( it_brkpts->second  >  it_best_brkpts__this_event->second )
		it_best_brkpts__this_event = it_brkpts;
	}   
	
	diploid_marginal_MAP_calls[it_ev->first].the_diploid_profile_brkpts = it_best_brkpts__this_event->first;
	diploid_marginal_MAP_calls[it_ev->first].P_diploid_breakpoint = it_best_brkpts__this_event->second / brkpts_normalization_constant;                     
    }//ev        

    
    
    append_diploid_Sampled_Events_to_file(diploid_marginal_MAP_calls, output_dir);    
                                                
// //format:  
// //	genome_name
// //	population 
// //	recomb_outcome
// //	P(outcome)
// //	chr_lower
// //	lower_brkpt
// //	chr_upper
// //	upper_brkpt
// //	P(brkpt)
// //	amount_affected
// //	CC_ID
// //	EV_ID
// //	breakpoint lower
// //	breakpoint upper               
// 
//     std::cerr << "\n\t\t\twrite to streams\n";   
//     const type_haploid_outcome__haploid_outcome  diploid_no_event__outcome(hap_outcome__None, hap_outcome__None);
//     
//     for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_ev = diploid_marginal_MAP_calls.begin();
//             it_ev != diploid_marginal_MAP_calls.end();
//             ++it_ev)
//     {     	
// 	const type_map_uint_to_Event::const_iterator occurring_Event = the_Conn_Comp.events.find( it_ev->first );
// 	
// 	for (uint hap=0; hap<2; ++hap)
// 	{	    	
// 	    if (pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap) == hap_outcome__None)
// 	    {
// 		std::stringstream temp_osstrm;
// 		temp_osstrm.precision(10);
// 		temp_osstrm
// 			<< genome_name
// 			<< "\t" << population_abbreviation
// 			<< "\tnone"
// 			<< "\t" << it_ev->second.P_diploid_outcome
// 			<< "\tchr " <<  occurring_Event->second.chromos[0]
// 			<< ": " <<  occurring_Event->second.region_between_and_including_the_LCRs_themselves.lower()
// 			<< "\tchr " <<  occurring_Event->second.chromos[0]
// 			<< ": " <<  occurring_Event->second.region_between_and_including_the_LCRs_themselves.upper()
// 			<< "\t1"
// 			<< "\t" <<  BOOST_width_inclusive( occurring_Event->second.region_between_and_including_the_LCRs_themselves )                                                   
// 			<< "\t" << the_Conn_Comp.Conn_Comp_ID
// 			<< "\t" << occurring_Event->first
// 			<< "\t-\t-\n";
// 						    
// 		    output_marginal_Negative_calls_strm << temp_osstrm.str();
// 		    output_marginal_all_calls_strm << temp_osstrm.str();
// 	    }
// 	    else
// 	    {   		
// 		BOOST_Interval affected_absolute_brkpts;		
// 		{
// 		    const type_BI__BI absolute_brkpts(
// 					    convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
// 						pair_at<BOOST_Interval>(it_ev->second.the_diploid_profile_brkpts,hap),
// 						occurring_Event->second.compressed_map_LCR_to_profile__for_each_LCR,
// 						pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap)));					
// 
// 		    if (test_if_haploid_outcome_is_NAHR(pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap)))
// 		    {
// 			affected_absolute_brkpts = absolute_brkpts.first;
// 		    }
// 		    else
// 		    {
// 			if (pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap) == hap_outcome__GeneConv_ABA)
// 			{
// 			    affected_absolute_brkpts = absolute_brkpts.first;
// 			}
// 			else
// 			{
// 			    affected_absolute_brkpts = absolute_brkpts.second;
// 			}
// 		    }//GeneConv
// 		}
// 										   	    
// 		std::stringstream table_entry_this_outcome; 
// 		table_entry_this_outcome.precision(10);		
// 		table_entry_this_outcome 
// 		    <<  genome_name
// 		    <<  "\t" << population_abbreviation
// 		    <<  "\t" << convert_haploid_outcome_to_string(pair_at<type_haploid_outcome>(it_ev->second.the_diploid_Event_outcome,hap)) //recomb outcome                                      
// 		    << "\t" << it_ev->second.P_diploid_outcome                         // P outcome                  
// 		    << "\tchr " << occurring_Event->second.chromos[0]                           // chr lower brkpt
// 		    << ": " << affected_absolute_brkpts.lower()               // lower brkpt
// 		    << "\tchr " << occurring_Event->second.chromos[0]                           // chr upper brkpt
// 		    << ": " << affected_absolute_brkpts.upper()               // upper brkpt
// 		    << "\t" << it_ev->second.P_diploid_breakpoint   // P brkpt                                                                                                  
// 		    << "\t" << BOOST_width_inclusive(  affected_absolute_brkpts  )    // amount of DNA affected                            
// 		    << "\t" << the_Conn_Comp.Conn_Comp_ID
// 		    << "\t" << occurring_Event->first
// 		    << "\t" << pair_at<BOOST_Interval>(it_ev->second.the_diploid_profile_brkpts,hap).lower()
// 		    << "\t" << pair_at<BOOST_Interval>(it_ev->second.the_diploid_profile_brkpts,hap).upper()
// 		    << "\n";		    		
// 				    		    
// 		output_marginal_Positive_calls_strm << table_entry_this_outcome.str();
// 		output_marginal_all_calls_strm << table_entry_this_outcome.str();								    
// 	    }//outcome != 0                
// 	}//hap                                      
//     }//ev
                
    
} // make_MAP_calls_from_marginal_posterior_distributions_per_Event































































type_BI__BI  convert_each_profile_breakpoint_to_absolute_breakpoints
				(const BOOST_Interval &profile_brkpts,
				 const type_vector_vector_vector_int &compressed_maps_LCRs_to_profile__of_some_Event)
{
    return  type_BI__BI(
		    convert_profile_breakpoint_to_absolute_breakpoints
				(profile_brkpts.lower(),
				 compressed_maps_LCRs_to_profile__of_some_Event),
		    convert_profile_breakpoint_to_absolute_breakpoints
				(profile_brkpts.upper(),
				 compressed_maps_LCRs_to_profile__of_some_Event)  );
}				 



BOOST_Interval  convert_profile_breakpoint_to_absolute_breakpoints
				(const uint &profile_brkpt_index,
				 const type_vector_vector_vector_int &compressed_maps_LCRs_to_profile__of_some_Event)
{        
    //get absolute brkpts
    type_uint__uint absolute_brkpts;       
    
    for (uint qs = 0; qs < 2; ++qs)
    {    
	const BOOST_Interval desired_upper_interval(profile_brkpt_index, std::numeric_limits<int>::max()); // NOT UINT,  but we want MAX for INT !!!!!  (when comparing vlaues inside of "compressed map" routines, we cast to "int" for many comparisons.  Hence, the numeric max of "uint" would be out of range, giving us bizarre answers!!!!
	
	const type_vector_vector_int compressed_map__profile_to_LCR(
				    get_compressed_map_inverse_of_compressed_map(compressed_maps_LCRs_to_profile__of_some_Event.at(qs)));
	
	const type_vector_vector_int  compressed_map__profile_to_LCR__intersecton(
		    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(compressed_map__profile_to_LCR, desired_upper_interval));
		    
	const type_vector_vector_int *const  ptr__compressed_map_profile =   compressed_map__profile_to_LCR__intersecton[0].empty()  ?
						&compressed_map__profile_to_LCR   :    &compressed_map__profile_to_LCR__intersecton;
		
						
	const type_int__int  closest_profile_index__to_abs_coord = get_extremal_entry_in_compressed_map(*ptr__compressed_map_profile, false);
	pair_at<uint>(absolute_brkpts,qs) = (uint)closest_profile_index__to_abs_coord.second; 	    
	
    }   
    
//     assert(absolute_brkpts.first <= absolute_brkpts.second);
    
    return BOOST_Interval(absolute_brkpts.first,  absolute_brkpts.second);

} // convert_profile_breakpoint_to_absolute_breakpoints




type_BI__BI convert_absolute_breakpoints_to_outcome_meaningful_intervals
		(const type_BI__BI &absolute_brkpts_each_LCR,
		 const type_haploid_outcome &corresponding_hap_outcome)
{
    if (corresponding_hap_outcome == hap_outcome__None)
    {  return empty_BI__BI;  }
    else if (test_if_haploid_outcome_is_GeneConversion(corresponding_hap_outcome))
    {  return absolute_brkpts_each_LCR;  }
    else //if (corresponding_hap_outcome == hap_outcome__Del or corresponding_hap_outcome == hap_outcome__Dup)   or  INV  or TRANSLOC
    {  
	const BOOST_Interval whole_reg(absolute_brkpts_each_LCR.first.lower(), absolute_brkpts_each_LCR.second.upper());
	return type_BI__BI(whole_reg,whole_reg);
    }   

}//convert_absolute_breakpoints_and_outcomes_to_meaningful_intervals






//return value:  if  GeneConv: ".first" = breakpoint interval on LCR 0.  ".second" on LCR 1
//		else :  ".first" = ".second" = hull of breakpoints.  (breakpoint LCR 0 through breakpoint LCR 1)
type_BI__BI convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals
				(const BOOST_Interval &profile_brkpts,
				 const type_vector_vector_vector_int &compressed_maps_LCRs_to_profile__of_some_Event,
				const type_haploid_outcome &corresponding_hap_outcome)
{     
    return  convert_absolute_breakpoints_to_outcome_meaningful_intervals(
			    BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(
					convert_each_profile_breakpoint_to_absolute_breakpoints(profile_brkpts, compressed_maps_LCRs_to_profile__of_some_Event)),
			    corresponding_hap_outcome);
    
}//convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals


























// i.e. double these numbers to get the diploid rate
void  calculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP
		    (type_map_uint_to_list_BI__string &pre_uploaded_regions_across_all_chromosomes,
		     const uint &chromosome_value,
		     type_map_uint_to_longdouble &positions_along_chromosome,
		     const uint &CC_ID,
		     const uint &EV_ID,
		     const uint &color_of_interval)
{

    type_map_uint_to_list_BI__string::const_iterator  it_to_uploaded_chromo = pre_uploaded_regions_across_all_chromosomes.find(chromosome_value);
    if (it_to_uploaded_chromo == pre_uploaded_regions_across_all_chromosomes.end())
    {
	std::stringstream error_strm;
	error_strm << "\n\nERROR!   unable to find chromo  " << chromosome_value 
		    << " of natural poisson interval color " << color_of_interval
		    << "  of event " << EV_ID  << "  of connected component  " << CC_ID
		    << "  among uploaded chromosomes!!!\n\t\t\tperforming an emergency upload of the entire chromosome!!!\n\n";
		    
	print_set<uint>(  extract_keys_of_map_and_return_as_set<uint, type_list__BI__string>(pre_uploaded_regions_across_all_chromosomes),
			    "chromosomes in \"pre_uploaded_regions_across_all_chromosomes\"", &error_strm, true);
				    
	if (Event::chromosome_lengths.count(chromosome_value) > 0)
	{
	    std::stringstream warning_strm;
	    warning_strm << "\n\nDOUBLE ERROR!     UNABLE TO identify  chromosome  "  << chromosome_value << " among chromosome lengths!!!\n\n";	    
	    print_map_keys_and_values<uint,uint>(Event::chromosome_lengths, "Event::chromosome_lengths", &warning_strm, true);	    
	    error_message(warning_strm,true);          
	}
	
	error_message(error_strm,false);         
	
	
	
	const BOOST_Interval whole_chromo_cooords(1, Event::chromosome_lengths.at(chromosome_value));
			    
	type_map_uint_to_list_BI__string::iterator it_up_chro
			=  pre_uploaded_regions_across_all_chromosomes.insert(pre_uploaded_regions_across_all_chromosomes.end(),   std::pair<uint,type_list__BI__string>
												    (chromosome_value,  type_list__BI__string()));
											    
	it_up_chro->second.push_back(type_BI__string(whole_chromo_cooords, get_specific_DNA_sequence_from_Reference_genome(chromosome_value, whole_chromo_cooords))); 								    
	it_to_uploaded_chromo = it_up_chro;                                                                       
    
    }//unalbe to find chromo
    
    
    
    
    
    const std::vector< type_map_uint_to_longdouble::iterator>  pos____loop(get_loopable_iterators<type_map_uint_to_longdouble>(positions_along_chromosome));
    const uint number_of_pos_this_chromo_this_nat = pos____loop.size();
    
    
//     if (!pos____loop.empty()  and   !BOOST_in(pos____loop[number_of_pos_this_chromo_this_nat-1]->first + maximum_average_fragment_length, ) )
					
					
    #pragma omp parallel
    {//parallel            

	type_map_uint_to_list_BI__string  emergency_uploaded_regions___thread_specific;
		
	#pragma omp for  schedule(dynamic,100)
	for (uint pos_ctr = 0; pos_ctr < number_of_pos_this_chromo_this_nat; ++pos_ctr)
	{
	    
	    //get the right interval and string.
	    type_list__BI__string::const_iterator  it_relevant_piece_of_chr = it_to_uploaded_chromo->second.end();
			for (type_list__BI__string::const_iterator it_piece = it_to_uploaded_chromo->second.begin();
				it_piece != it_to_uploaded_chromo->second.end();
				++it_piece)
			{
			    if ( BOOST_in(pos____loop.at(pos_ctr)->first, it_piece->first))
			    {
				it_relevant_piece_of_chr = it_piece;
				break;                    
			    }
			}
	    
	    
												    
	    if (it_relevant_piece_of_chr ==  it_to_uploaded_chromo->second.end()
		or   !BOOST_in(std::min<uint>(Event::chromosome_lengths.at( chromosome_value), pos____loop.at(pos_ctr)->first + maximum_average_fragment_length),
				it_relevant_piece_of_chr->first))
	    {
		bool must_do_emergency_upload = false;
		
		type_map_uint_to_list_BI__string::iterator it_emergency_chr = emergency_uploaded_regions___thread_specific.find( chromosome_value );
		
		if (it_emergency_chr == emergency_uploaded_regions___thread_specific.end())
		{    must_do_emergency_upload = true;  }
		else
		{                                                
				for (type_list__BI__string::const_iterator it_piece = it_emergency_chr->second.begin();
					it_piece != it_emergency_chr->second.end();
					++it_piece)
				{
				    if (BOOST_in(pos____loop.at(pos_ctr)->first, it_piece->first)  and  BOOST_in(pos____loop.at(pos_ctr)->first + maximum_average_fragment_length, it_piece->first))
				    {
					it_relevant_piece_of_chr = it_piece;
					break;                    
				    }    
				}
				   
				    
		    if (it_relevant_piece_of_chr ==  it_to_uploaded_chromo->second.end())  //still!        
		    {  must_do_emergency_upload = true;  }
		}//found emergency chr
		
		
		if (must_do_emergency_upload)
		{                                                        
		    std::stringstream error_strm;
		    error_strm << "\n\nERROR!    unable to find relevant pre_uploaded chromos or emergency chromos in \"calculate_total_GC_rate_for_every_base\" !!!\n\n";
		    error_strm << "\n\tnatural poisson color = " << color_of_interval
				<< "\n\tchromo = " << chromosome_value
				<< "\n\tposition = " << pos____loop.at(pos_ctr)->first;
		    error_strm << "\n\nuploaded chromo regions for this chromo...\n\n";
		    for (type_list__BI__string::const_iterator it_piece = it_to_uploaded_chromo->second.begin();
			    it_piece != it_to_uploaded_chromo->second.end();
			    ++it_piece) 
			error_strm << "\t" << it_piece->first << "\n";       
		    
		    
		    
		    const BOOST_Interval spontaneous_interval( safe_subtract_base_1(pos____loop.at(pos_ctr)->first, 30000),
							       safe_chromo_add_base_1(pos____loop.at(pos_ctr)->first, 30000, chromosome_value));
			
		    error_strm << "\n\nemergency uploading region   chr: " << chromosome_value 
				<<  "  " << spontaneous_interval  << "    for this thread.\n\n\n";                                                                                                                                                    
		    
		    error_message(error_strm, false);
		    
		    emergency_uploaded_regions___thread_specific[chromosome_value]
			    .push_back(  type_BI__string( spontaneous_interval,
							    get_specific_DNA_sequence_from_Reference_genome( chromosome_value, spontaneous_interval )   ));  
							    
		    it_relevant_piece_of_chr = --emergency_uploaded_regions___thread_specific.at(chromosome_value).end();                                         
		}//emergency upload
	    }
					    
			
	    
	    longdouble aggregate_RG_GC_fragmentation_rate = 0.00L;
	    	    
	    
	    //go through each read group and count the GC content according to the appropriate fragment length for that readgroup.  Then get the ratio for the read group.
	    for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
		    it_rg != Readgroup_stats.end();
		    ++it_rg)
	    {		
		uint GC_count = 0;
		const uint relative_start_pos = pos____loop.at(pos_ctr)->first - it_relevant_piece_of_chr->first.lower();
		const uint relative_end_pos = std::min<uint>(pos____loop.at(pos_ctr)->first + (uint)it_rg->second.median_insert_size,
							     Event::chromosome_lengths.at(chromosome_value)+1)
							- it_relevant_piece_of_chr->first.lower();
		
		for (uint j = relative_start_pos; j < relative_end_pos; ++j)
		{
		    if (it_relevant_piece_of_chr->second.at(j) == 'G'  or  it_relevant_piece_of_chr->second.at(j) == 'C')
		    {  ++GC_count;  }
		}
		    
		const type_map_uint_to_longdouble::const_iterator it_find_GC = it_rg->second.haploid_fragmentation_rate__per_GC_count.find(GC_count);
		
		if ( it_find_GC == it_rg->second.haploid_fragmentation_rate__per_GC_count.end() )
		{
		    std::stringstream error_strm;
		    error_strm << "\n\nERROR!   could not find GC_count = "  <<  GC_count
				<< " in Read group  " << it_rg->first << " for chr " << chromosome_value <<" : " <<  pos____loop.at(pos_ctr)->first 
				<< "\n\t\tusing some approximation based on nearby neighbors instead (if is near beginning or near end)\n\n\n";
		    it_rg->second.print_this_Readgroup(&error_strm);                                                                        
		    
		    error_message(error_strm,false);		    
		    
		    if (GC_count < it_rg->second.haploid_fragmentation_rate__per_GC_count.begin()->first)
			aggregate_RG_GC_fragmentation_rate  +=  it_rg->second.haploid_fragmentation_rate__per_GC_count.begin()->second;
		    else if (GC_count >   (--it_rg->second.haploid_fragmentation_rate__per_GC_count.end())->first     )
			aggregate_RG_GC_fragmentation_rate  +=   (--it_rg->second.haploid_fragmentation_rate__per_GC_count.end())->second;
		    else 
			aggregate_RG_GC_fragmentation_rate += it_rg->second.haploid_fragmentation_rate__per_GC_count.at(GC_count);
		}
		else
		{                                                                        
		    aggregate_RG_GC_fragmentation_rate += it_find_GC->second;                      
		    
		    if (it_find_GC->second > 3.00L)
		    {
			std::stringstream error_strm;			
			error_strm << "\n\nERROR!   GC_count  " <<  GC_count <<" has rate =   " <<  it_find_GC->second
				    <<  " for            chr " << chromosome_value << "  pos   " << pos____loop.at(pos_ctr)->first <<  "\n\n";
			
			error_message(error_strm,false);                     
		    }
		}
	    }// rg
			    
	    pos____loop.at(pos_ctr)->second = aggregate_RG_GC_fragmentation_rate;  	    
		
	}  //pos 
	
	
	
	#pragma omp critical (add_novel_emergency_uploaded_regions_to_preuploaded_regions_for_this_Event)
	{//critical
	    for (type_map_uint_to_list_BI__string::const_iterator  it_emerg_chr  =  emergency_uploaded_regions___thread_specific.begin();
		    it_emerg_chr !=  emergency_uploaded_regions___thread_specific.end();
		    ++it_emerg_chr)
	    {
		for (type_list__BI__string::const_iterator it_emerg_reg = it_emerg_chr->second.begin();
			it_emerg_reg != it_emerg_chr->second.end();
			++it_emerg_reg)
		{
		    const type_map_uint_to_list_BI__string::iterator it_find_exist_chr 
								    =  pre_uploaded_regions_across_all_chromosomes.find( it_emerg_chr->first );
		    if (it_find_exist_chr == pre_uploaded_regions_across_all_chromosomes.end() )
			pre_uploaded_regions_across_all_chromosomes[it_emerg_chr->first].push_back( *it_emerg_reg );
		    else           
		    {
			bool has_a_significantly_overlapping_reg = false;
			for (type_list__BI__string::iterator it_exist_reg = it_find_exist_chr->second.begin();
				it_exist_reg != it_find_exist_chr->second.end();
				++it_exist_reg)
			{
			    if (    (double)BOOST_width_inclusive( BOOST_intersect(it_exist_reg->first, it_emerg_reg->first ) ) 
				    /  (double)BOOST_width_inclusive(it_emerg_reg->first)   >   0.75    )
			    {
				has_a_significantly_overlapping_reg = true;
				break;
			    }
			}
			    
			if (!has_a_significantly_overlapping_reg)
			    it_find_exist_chr->second.push_back(*it_emerg_reg);
		    }
		}//emerg_reg
	    }//emerg_chr
	}//critical
			    
			    
	mpfr_free_cache();
    }//parallel    



} //calculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP

























// void  display_MAP_PER_space_in_some_unique_region_of_the_genome
// 		(const uint &chromosome,
// 		 const uint &some_genome_position___base_1,
// 		 const uint &padding_amount_for_position)
// {
// 
//     std::cerr << "\n\n\n\n\"display_MAP_PER_space_in_some_unique_region_of_the_genome\"...\n\n";
//     std::cerr << "\n\t\tchromosome = " << chromosome
// 	    << "\n\t\tsome_genome_position___base_1 = " << some_genome_position___base_1 
// 	    << "\n\t\tpadding_amount_for_position = " << padding_amount_for_position << "\n\n\n";
//     
//     const BOOST_Interval search_region( some_genome_position___base_1 - 1 - padding_amount_for_position,
// 					 some_genome_position___base_1 - 1 + padding_amount_for_position);
// 					 
//     const BOOST_Interval broader_region( some_genome_position___base_1 - 1 - padding_amount_for_position*10,
// 					 some_genome_position___base_1 - 1 + padding_amount_for_position*10);
//     					 
//     
//     type_list_PER found_PERs;
//     
//     //raw upload all PERs, pair them, and give them spawning MiniRegions.
//     
//     bool success_raw_upload = true;
//     
//     #pragma omp parallel
//     {
//         BamTools::BamReader *const ptr_to_thread_BAM_Reader = &(*(my_BAM_readers.begin() + (uint)omp_get_thread_num())); //ptr arithmetic  
//         
//                        
// 	type_list_Balgmts found_Balgmts;            
// 	
// 	const bool success_set_region 
// 	    = ptr_to_thread_BAM_Reader->SetRegion(
// 				map_chromosome_value_to_BAM_Ref_IDs.at(chromosome),
// 				search_region.lower(),                 //BAMfiles are 0-based
// 				map_chromosome_value_to_BAM_Ref_IDs.at(chromosome),
// 				search_region.upper()      );  //BAMfiles are 0-based
// 	    
// 	if (!success_set_region)
// 	{
// 	    std::stringstream error_strm;
// 	    print_line_of_markers("ERROR! ", &error_strm);
// 	    print_line_of_markers("(", &error_strm);
// 	    error_strm << "success_set_region  failed in  \"display_MAP_PER_space_in_some_unique_region_of_the_genome\"\n\n";
// 	    
// 	    print_line_of_markers(")", &error_strm);
// 	    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                      
// 	}                           
// 	else
// 	{                
// 	    BamTools::BamAlignment Balgmt;
// 	    
// 	    while ( ptr_to_thread_BAM_Reader->GetNextAlignmentCore(Balgmt) )  
// 	    {
// 		if ( !Balgmt.IsFailedQC() )  
// 		{
// 		    if ( Balgmt.IsProperPair()     or     ( Balgmt.IsPaired()  and  Balgmt.IsMapped()  and   Balgmt.IsMateMapped()   )       )  
// 		    {
// 			Balgmt.BuildCharData();
// 			found_Balgmts.push_back( Balgmt );                    
// 		    }
// 		}
// 	    }//while
// 	}//success
// 	
// 
//     
// 	//Now pair up these reads to form Paired_end_reads;
// 	found_PERs =  create_list_of_unformatted_Paired_end_reads_from_list_of_Balgmts(found_Balgmts, ptr_to_thread_BAM_Reader);
// 	found_Balgmts.clear(); //should already be erased anyway!!  (just for emphasis).
// 
//     
//         mpfr_free_cache();
//     }//omp parallel
//         
//         
//         
//         
//         
//     {
//         std::stringstream  diag_ss;
//         diag_ss   << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() 
//                     << "  inside \"display_MAP_PER_space_in_some_unique_region_of_the_genome...\"  -    uploaded raw PERs.   Now adding miniregions...\n\n";
//         std::fprintf(stderr, "\n\n%s\n\n", diag_ss.str().c_str() );
//     }
//         
//         
//         
//         
//     
//         
//         
//     type_map_uint_to_list_BI__string uploaded_pieces_of_chromos;
//     uploaded_pieces_of_chromos[chromosome].push_back(
// 			type_BI__string(broader_region, get_specific_DNA_sequence_from_Reference_genome(chromosome, broader_region)));
//         
//         
//     type_vector_int  foundation_coords;
//     foundation_coords.reserve( BOOST_width_inclusive(broader_region) );
//     for (uint pos=broader_region.lower(); pos <= broader_region.upper(); ++pos)
// 	foundation_coords.push_back(pos);
// 	            
//     char region_name[option_size];
//     std::sprintf(region_name, "chr %u:  [%u,  %u]", chromosome, broader_region.lower(), broader_region.upper());
//     
//     Star_alignment unique_region_Star
// 			(foundation_coords,
// 			 get_specific_DNA_sequence_from_Reference_genome(chromosome, broader_region),
// 			 region_name,
// 			 (int)BOOST_width_inclusive(broader_region));
//                 
//     for (type_list_PER::iterator it_found_per = found_PERs.begin();
// 	    it_found_per != found_PERs.end();
// 	    ++it_found_per)
//     {	
// 	it_found_per->add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately(
// 			MiniRegion(chromosome, broader_region, true, false ));			
//             
// 	it_found_per->set_sequence_for_all_MiniRegions( uploaded_pieces_of_chromos );
// 	it_found_per->orient_and_complement_each_mate_according_to_its_spawning_Event_profile(NULL);	
// 	
// 	//align:
// 	for (type_map_uint_to_MiniRegion::const_iterator it_mini_A = it_found_per->my_MiniRegions.begin();
// 		it_mini_A != it_found_per->my_MiniRegions.end();
// 		++it_mini_A)
// 	{
// 
// 	    std::stringstream  nohybrid_ss;
// 	    BOOST_Interval mini_A_algmt_endpts;
// 	    type_string__string PER_and_mini_A_algmt;
// 					
// 	    const real P_nohybrid_A(
// 			    it_found_per->my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
// 						    (it_found_per->mates.first,
// 						    it_found_per->mates.second,
// 						    it_mini_A->second.region_sequence,
// 						    it_mini_A->second.orientation_of_MiniRegion_and_spawning_Event_profile_agree,
// 						    it_mini_A->second.MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile,
// 						    PER_and_mini_A_algmt,
// 						    mini_A_algmt_endpts)             ); 
// 						    
//     
// 	    nohybrid_ss << "\n" << it_mini_A->second.chromosome_of_region 
// 				<< ":  [" << it_mini_A->second.region_interval.lower() 
// 				<< ", " << it_mini_A->second.region_interval.upper() << "]\n";
// 											
// 
// 	    const type_map_uint_to_uint  map_mini_A_indeces_to_algmt_profile(  get_relative_hybrid_coords_to_algmt_profile(PER_and_mini_A_algmt,mini_A_algmt_endpts)   );
// 	    
// 	    if (   (--map_mini_A_indeces_to_algmt_profile.end())->first  >=   it_mini_A->second.region_sequence.size() )
// 	    {
// 		std::stringstream error_strm;
// 		print_line_of_markers("ERROR! ", &error_strm);
// 		print_line_of_markers("(", &error_strm);
// 		error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
// 		
// 		error_strm << "\n\nERROR!   in no-hybrid displays:\n\n\n"
// 			    << "(--map_mini_A_indeces_to_algmt_profile.end())->first  =  "  <<  (--map_mini_A_indeces_to_algmt_profile.end())->first
// 			    << "   >=     " << it_mini_A->second.region_sequence.size() <<" =  it_mini_A->second.region_sequence.size()\n\n";
// 			    
// 		error_strm  << "\nmini_A_algmt_endpts = " << mini_A_algmt_endpts
// 			    << "\n\nPER_and_mini_A_algmt:\nPER:\n" << PER_and_mini_A_algmt.first
// 			    << "\nmini:\n" << PER_and_mini_A_algmt.second
// 			    << "\n\n";
// 			    
// 		it_mini_A->second.print_this_MiniRegion( &error_strm );
// 		
// 		it_found_per->print_this_PER(  &error_strm  );
// 		
// 		print_map_keys_and_values<uint,uint>(map_mini_A_indeces_to_algmt_profile, "map_mini_A_indeces_to_algmt_profile",  &error_strm,  true );
// 			    
// 		
// 		print_line_of_markers(")", &error_strm);
// 		std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                
// 	    }
//         
//         
//         
//         
// 	    
// 	    type_map_uint_to_uint  dummy_empty_map;
// 	    type_map_uint_to_uint  map_mini_A_absolute_to_algmt_profile;
// 			    adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
// 					    (  map_mini_A_indeces_to_algmt_profile,
// 						BOOST_width_inclusive(  it_mini_A->second.region_interval  ),  //i.e.  beyond the end
// 						it_mini_A->second.region_interval.lower(),
// 						false,  //never reversed the MR,  only possibly the PER
// 						0,  //irrelevant
// 						false,   //irrelevant
// 						map_mini_A_absolute_to_algmt_profile,   
// 						dummy_empty_map     );           //irrelevant     
// 					    
// 	    assert(  !map_mini_A_absolute_to_algmt_profile.empty()  );                                
//         
//                             
// 	    it_found_per->calculate_probability_of_PER_given_hybrid_AND_DISPLAY
// 					    (   map_mini_A_indeces_to_algmt_profile,  //get_relative_hybrid_coords_to_algmt_profile(PER_and_mini_A_algmt,mini_A_algmt_endpts),
// 						PER_and_mini_A_algmt,
// 						it_mini_A->second.region_sequence.size(),  //will never reach this... (this is the intention)
// 						it_mini_A->second.region_interval.lower(), 
// 						false,
// 						it_mini_A->second.chromosome_of_region,
// 						0,
// 						false,
// 						0,
// 					    nohybrid_ss,
// 					    false,
// 					    NULL);         // be careful here!! (we are.)	
// 					    
// 	    unique_region_Star.add_pairwise_alignment_to_progressive_star_alignment
//                                                     ( P_nohybrid_A,
//                                                     PER_and_mini_A_algmt,
//                                                     map_mini_A_absolute_to_algmt_profile,
//                                                     dummy_empty_map,                                            
//                                                     true,
//                                                     false,
// 						      it_found_per->mates,
// 						      P_nohybrid_A);  	
// 	}
//     }//found_PERs
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
//     {//save star alignments                       
//         char star_algmts_filename[option_size];        
//         std::sprintf(star_algmts_filename, "%s/unique_star_alignments___chr_%u___%u___%s.html",
//                                         output_dir.c_str(),					
//                                         chromosome, some_genome_position___base_1,
// 					genome_name.c_str()  );       
//                             
//         std::ofstream  save_star_algmts_filestream( star_algmts_filename );        
//                                     if (  !save_star_algmts_filestream.is_open()  )
//                                     {
//                                         std::fprintf(stderr, "\n\n\n\n\nERROR: Command [%s] could not be opened for \"save group of star_alignments_to_file_as_html\".\n\n\n\n\n\n", 
//                                                             star_algmts_filename);
//                                         return;
//                                     }                                           
//                                     
//                                     
//         save_star_algmts_filestream << "\n<pre>\n<font size=1>\n\n";            
//         save_star_algmts_filestream << "\n<span style=\"font-size: 24pt\">Genome:  "  << genome_name << "    chromosome  " <<  chromosome
//                                     << ",     position  " <<  some_genome_position___base_1
//                                     << "</span>\n\n";    				     
// 		
// 	type_vector_string  html_output;
// 	const bool success_create_html =  unique_region_Star.create_HTML_output_for_star_alignment(   
// 								    (int)chromosome, -1,
// 								    0,
// 								    true,
// 								    html_output    );
// 	
// 	if (!success_create_html)
// 	{
// 	    save_star_algmts_filestream.close();
// 	    boost::filesystem::remove(    boost::filesystem::path( star_algmts_filename )    );        
// 	    return;
// 	}
// 									    
// 	for (uint row=0; row < html_output.size(); ++row)
// 	    save_star_algmts_filestream << html_output.at(row) << "\n";       
// 	
// 	
//         save_star_algmts_filestream << "\n\n</font>\n</pre>\n";
//         save_star_algmts_filestream.close();   	    	        				    				    
//     }//save star alignments 
// 
// 
//     std::cerr << "\n\n\n\nDONE  \"display_MAP_PER_space_in_some_unique_region_of_the_genome\"...\n\n";
// 
// }//display_MAP_PER_space_in_some_unique_region_of_the_genome






















// type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble  construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region
// 		    (const Conn_Comp &the_CC,
// 		     const Sampled_diploid_Event_data &some_Event_result)
// {   
//     return  construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region(
// 				    the_CC.events.at(some_Event_result.event_UID).chromos[0],
// 				    the_CC.events.at(some_Event_result.event_UID).region_between_and_including_the_LCRs_themselves,
// 				    100000);
// }//construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region
// 
// 
// 
// // type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble  construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region
// // 		    (const Pseudo_Read_Depth_display_container &pseudo_RD_display_region)
// // {
// //     return  construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region(
// // 				    pseudo_RD_display_region.chromo_of_region,
// // 				    pseudo_RD_display_region.affected_regions_of_interest,
// // 				    pseudo_RD_display_region.padding_amount_for_neighboring_unique_regions);
// // }//construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region










void subtract_or_add_map_values_from_another_map__within_affected_region
			(const type_map_uint_to_longdouble &map_pos_to_adjustment_value,
			 const BOOST_Interval &affected_region,
			const bool &subtract__false____add_true,
			type_map_uint_to_longdouble &updated_map,
			type_map_uint_to_int  &contributions_per_position)
			
{
    const type_map_uint_to_longdouble restricted_adjustments(
		get_map_intersection_of_map_keys_and_interval<uint,longdouble>(map_pos_to_adjustment_value, affected_region));
    
    if (restricted_adjustments.empty())
    {  
	return;
    }        
    else
    {    	
	if (subtract__false____add_true)
	{//add
	    for (type_map_uint_to_longdouble::const_iterator it_adj = restricted_adjustments.begin();
		    it_adj != restricted_adjustments.end();
		    ++it_adj)
	    {
		updated_map.at(it_adj->first) += it_adj->second;
		++contributions_per_position.at(it_adj->first);
	    }
	}//add
	else
	{//subtract
	    for (type_map_uint_to_longdouble::const_iterator it_adj = restricted_adjustments.begin();
		    it_adj != restricted_adjustments.end();
		    ++it_adj)
	    {
		updated_map.at(it_adj->first) -= it_adj->second;
		--contributions_per_position.at(it_adj->first);
	    }	
	}//subtract
    }
    
}//subtract_or_add_map_values_from_another_map__within_affected_region








void  raw_add_expected_fragmentation_rate_of_homologous_regions__raw_haploid_fragmentation_rate
			    (const Event &the_Event,
			    const type_map_uint_to_uint_to_longdouble  &map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
			    type_map_uint_to_longdouble &updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    type_map_uint_to_int &contributions_per_position_to_updated_expected_fragmentation_rate) 
{
    
    //add from repetitive locations
    for (type_list_Coords_and_homologous_region::const_iterator it_homol = the_Event.regions_homologous_to_directly_affected_region.begin();
	    it_homol != the_Event.regions_homologous_to_directly_affected_region.end();
	    ++it_homol)
    {
	const type_map_uint_to_uint map_affected_region_of_interest_to_homol_reg(
		convert_compressed_map_to_full_map(
		    get_compressed_map_inverse_of_compressed_map(it_homol->compressed_map_homologous_region_to_profile))); //DAR to homol reg	     

	const type_map_uint_to_uint_to_longdouble::const_iterator it_chr_homol_reg_to_raw_rate 
	= map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate.find(it_homol->chromosome_of_homologous_region);
	
	for (type_map_uint_to_uint::const_iterator it_aff = map_affected_region_of_interest_to_homol_reg.begin();
		it_aff != map_affected_region_of_interest_to_homol_reg.end();
		++it_aff)
	{
	    updated_repeat_aware_expected_fragmentation_rate_per_positon.at(it_aff->first) += it_chr_homol_reg_to_raw_rate->second.at(it_aff->second);
	    ++contributions_per_position_to_updated_expected_fragmentation_rate.at(it_aff->first);
	}//it_aff
		    
    }//homol
		      
}//raw_add_expected_fragmentation_rate_of_homologous_regions__raw_haploid_fragmentation_rate
		      
		      











void  adjust_expected_fragmentation_rate_by_homologous_regions_raw_haploid_fragmentation_rate
			    (const Event &the_Event,
			    const type_map_uint_to_uint_to_longdouble  &map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
			    const BOOST_Interval &affected_region_on_chromosome_of_Event,
			    const bool &subtract_false____add_true,
			    type_map_uint_to_longdouble &updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    type_map_uint_to_int &contributions_per_position_to_updated_expected_fragmentation_rate) 
{
    
    //add from repetitive locations
    for (type_list_Coords_and_homologous_region::const_iterator it_homol = the_Event.regions_homologous_to_directly_affected_region.begin();
	    it_homol != the_Event.regions_homologous_to_directly_affected_region.end();
	    ++it_homol)
    {
	if (it_homol->chromosome_of_homologous_region == the_Event.chromos[0])
	{	
	    const type_map_uint_to_uint map_affected_region_of_interest_to_homol_reg(
			    convert_compressed_map_to_full_map(
				get_compressed_map_inverse_of_compressed_map( // DAR to homol reg
				    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map( // expanded reg to homol reg
					    it_homol->compressed_map_homologous_region_to_profile, //homol reg to DAR
					    affected_region_on_chromosome_of_Event))));    
	
	    if (map_affected_region_of_interest_to_homol_reg.empty())
	    {  continue;  }
	    
	    
	    if (subtract_false____add_true)
	    {
		const type_map_uint_to_uint_to_longdouble::const_iterator it_chr_homol_reg_to_raw_rate 
				    = map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate.find(it_homol->chromosome_of_homologous_region);
		
		for (type_map_uint_to_uint::const_iterator it_aff = map_affected_region_of_interest_to_homol_reg.begin();
			it_aff != map_affected_region_of_interest_to_homol_reg.end();
			++it_aff)
		{
		    updated_repeat_aware_expected_fragmentation_rate_per_positon.at(it_aff->first) += it_chr_homol_reg_to_raw_rate->second.at(it_aff->second);
		    ++contributions_per_position_to_updated_expected_fragmentation_rate.at(it_aff->first);
		}//it_aff				
	    }
	    else
	    {
		const type_map_uint_to_uint_to_longdouble::const_iterator it_chr_homol_reg_to_raw_rate 
				    = map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate.find(it_homol->chromosome_of_homologous_region);
		
		for (type_map_uint_to_uint::const_iterator it_aff = map_affected_region_of_interest_to_homol_reg.begin();
			it_aff != map_affected_region_of_interest_to_homol_reg.end();
			++it_aff)
		{
		    updated_repeat_aware_expected_fragmentation_rate_per_positon.at(it_aff->first) -= it_chr_homol_reg_to_raw_rate->second.at(it_aff->second);
		    --contributions_per_position_to_updated_expected_fragmentation_rate.at(it_aff->first);
		}//it_aff				
	    }		
	}//chromosome
    }//homol
        
}//adjust_expected_fragmentation_rate_by_homologous_regions_raw_haploid_fragmentation_rate














void  adjust_expected_fragmentation_rate_according_to_Event_outcome_and_breakpoints__on_repeat_and_raw_region
				(const Event &the_Event,
				 const type_haploid_outcome &the_haploid_outcome,				    
				const BOOST_Interval &the_haploid_breakpoints,
				const type_map_uint_to_uint_to_longdouble  &map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
				const type_map_uint_to_longdouble &raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
				type_map_uint_to_longdouble &updated_repeat_aware_expected_fragmentation_rate_per_positon,
				type_map_uint_to_int &contributions_per_position_to_updated_expected_fragmentation_rate)
{
    if (the_haploid_outcome == hap_outcome__None)
    {  return;  }
    
    
    
    const type_BI__BI meaningful_haploid_absolute_breakpoints
			    = convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
								the_haploid_breakpoints,
								the_Event.compressed_map_LCR_to_profile__for_each_LCR,
								the_haploid_outcome);    
			    
    
    if (test_if_haploid_outcome_is_GeneConversion(the_haploid_outcome))
    {//GeneConv

	const BOOST_Interval  affected_region__LCR_0(meaningful_haploid_absolute_breakpoints.first.lower(), meaningful_haploid_absolute_breakpoints.first.upper()-1);
	const BOOST_Interval  affected_region__LCR_1(meaningful_haploid_absolute_breakpoints.second.lower(), meaningful_haploid_absolute_breakpoints.second.upper()-1);    

	{//erased	    
	    const BOOST_Interval *const erased_tract__ptr 
				=   (the_haploid_outcome == hap_outcome__GeneConv_ABA)    ?
					&affected_region__LCR_0   :   &affected_region__LCR_1;
		
	    //raw adjust on exapnded region of interest
	    subtract_or_add_map_values_from_another_map__within_affected_region(
				    raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
				    *erased_tract__ptr,
				    false,
				    updated_repeat_aware_expected_fragmentation_rate_per_positon,
				    contributions_per_position_to_updated_expected_fragmentation_rate);
	    
	    //repeat-adjust from hoomologous positions
	    adjust_expected_fragmentation_rate_by_homologous_regions_raw_haploid_fragmentation_rate
			    (the_Event,
			    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
			    *erased_tract__ptr,
			    false,
			    updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    contributions_per_position_to_updated_expected_fragmentation_rate);
	    	    
	}//erased
	
	{//copied	  
	    const BOOST_Interval *const copied_tract__ptr 
				    =   (the_haploid_outcome == hap_outcome__GeneConv_ABA)    ?
					    &affected_region__LCR_1   :  &affected_region__LCR_0;    
	
	    //raw adjust
	    subtract_or_add_map_values_from_another_map__within_affected_region(
				    raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
				    *copied_tract__ptr,
				    false,
				    updated_repeat_aware_expected_fragmentation_rate_per_positon,
				    contributions_per_position_to_updated_expected_fragmentation_rate);
	    
	    //repeat-adjust from hoomologous positions
	    adjust_expected_fragmentation_rate_by_homologous_regions_raw_haploid_fragmentation_rate
			    (the_Event,
			    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
			    *copied_tract__ptr,
			    true,
			    updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    contributions_per_position_to_updated_expected_fragmentation_rate);
	}//copied
    
    }//GeneConv
    else if (test_if_haploid_outcome_is_NAHR(the_haploid_outcome))
    {//NAHR - dup or del
    
	//assert(BOOST_is_a_point(proposed_haploid_breakpoints.at(occurring_interlocking_event_it->first)));  // not necessarily a point, as we may have skipped some "hidden" brkpts that are not being summed out at the moment...
    
	const BOOST_Interval affected_region(meaningful_haploid_absolute_breakpoints.first.lower(), meaningful_haploid_absolute_breakpoints.first.upper()-1);
		//yes, you do this for both dups and dels.
		//carefully draw some pictures and think hard.
		//recall that you don't even get to this point in the code for inversions (or translocs).               	
			
	//raw adjust
	subtract_or_add_map_values_from_another_map__within_affected_region(
				    raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
				    affected_region,
				    (the_haploid_outcome == hap_outcome__Dup),
				    updated_repeat_aware_expected_fragmentation_rate_per_positon,
				    contributions_per_position_to_updated_expected_fragmentation_rate);
	
	//repeat-adjust from hoomologous positions
	adjust_expected_fragmentation_rate_by_homologous_regions_raw_haploid_fragmentation_rate
			    (the_Event,
			    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
			    affected_region,
			    (the_haploid_outcome == hap_outcome__Dup),
			    updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    contributions_per_position_to_updated_expected_fragmentation_rate);	
	
	
    }//NAHR    
	
}//adjust_expected_fragmentation_rate_according_to_Event_outcome_and_breakpoints__on_repeat_and_raw_region




















void appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
				(const Event &the_Event,
				 const type_haploid_outcome &the_haploid_outcome,
				 const BOOST_Interval &the_haploid_breakpoints,
				const type_map_uint_to_longdouble &raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
				const BOOST_Interval &expanded_region_of_interest,
				const type_map_uint_to_uint_to_longdouble  &map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
				type_map_uint_to_longdouble &updated_repeat_aware_expected_fragmentation_rate_per_positon,
				type_map_uint_to_int &contributions_per_position_to_updated_expected_fragmentation_rate)
{    
        
    //add raw NULL rate for region of interest
    subtract_or_add_map_values_from_another_map__within_affected_region(
			    raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
			    expanded_region_of_interest,
			    true,
			    updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    contributions_per_position_to_updated_expected_fragmentation_rate);
    
    
    //add repeat-aware NULL rate for region of interest
    raw_add_expected_fragmentation_rate_of_homologous_regions__raw_haploid_fragmentation_rate(
			    the_Event,
			    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
			    updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    contributions_per_position_to_updated_expected_fragmentation_rate);    
      

    adjust_expected_fragmentation_rate_according_to_Event_outcome_and_breakpoints__on_repeat_and_raw_region(
			    the_Event,
			    the_haploid_outcome,				    
			    the_haploid_breakpoints,
			    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
			    raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
			    updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    contributions_per_position_to_updated_expected_fragmentation_rate);
       
    
}//appropriately_adjust_expected_fragmentation_rate_and_numbeR_of_contributions_according_to_haploid_rate_and_event_outcome














type_map_uint_to_longdouble  divide_sum_by_contributions
					(const type_map_uint_to_longdouble &the_map_of_sums,
					const type_map_uint_to_int &the_map_of_number_ofcontributions_to_the_sums)    
{
    type_map_uint_to_longdouble the_map_of_rates(the_map_of_sums);
    
    type_map_uint_to_int::const_iterator it_contrib = the_map_of_number_ofcontributions_to_the_sums.begin();
    for (type_map_uint_to_longdouble::iterator it_rate = the_map_of_rates.begin();
	    it_rate != the_map_of_rates.end();
	    ++it_rate)
    {
	if (it_contrib->second > 0)
	{
	    it_rate->second /= it_contrib->second;
	}
	else
	{
	    it_rate->second = 0;	 
	}
	
	++it_contrib;
    }//rate
    
    
    return the_map_of_rates;
    
}//divide_sum_by_contributions













type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble  construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region
		    (const uint &chromo_of_region,
		     const BOOST_Interval &expanded_region_of_interest,
		    const type_haploid_outcome__haploid_outcome  &the_diploid_outcome,
		    const type_BI__BI &profile_breakpoints,
		    Event &the_Event)
{	        

    std::cerr << "\n\n\t\tinside \"construct_pseudo_RD_expected_Read_Depth_counts_under_NULL_hypothesis_for_desired_display_region\"!!!!\n";             
				    
    std::cerr << "\t\t\tpre_uploaded_regions_across_all_chromosomes...\n";
    //pre-upload relevant parts:
    type_map_uint_to_list_BI__string pre_uploaded_regions_across_all_chromosomes;
		{    
		    //upload extra so there are no problems with border-case PERs.		    
		    const uint pad_amount = std::max<uint>(40000, maximum_average_fragment_length*10);
		    const BOOST_Interval upload_region(safe_subtract_base_1(expanded_region_of_interest.lower(), pad_amount),
						    safe_chromo_add_base_1(expanded_region_of_interest.upper(), pad_amount, chromo_of_region)  );
		    
		    std::cerr << "\n\n\tupload_region = " << upload_region << "\nexpanded_region_of_interest = " << expanded_region_of_interest << "\n\n";
		    
		                                        
		    const std::string desired_sequence(get_specific_DNA_sequence_from_Reference_genome(chromo_of_region, upload_region));                                                                                                 
												    
		    pre_uploaded_regions_across_all_chromosomes[chromo_of_region].push_back(type_BI__string(upload_region, desired_sequence));     
		}    
    
    
    
    
    
    
    type_map_uint_to_longdouble expected_fragmentation_rate___HAPLOID_NULL;
    initialize_map_across_interval__via_insert<uint,longdouble>(expected_fragmentation_rate___HAPLOID_NULL, expanded_region_of_interest, 0.00L);

    std::cerr << "\t\t\tcalculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP...\n";
    calculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP
		    (pre_uploaded_regions_across_all_chromosomes,
		     chromo_of_region,
		     expected_fragmentation_rate___HAPLOID_NULL,
		     0, 0, 0);
			    
    
    
    
    
    type_map_uint_to_uint_to_longdouble map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate;
    {//raw homol reg rates    
    
	std::cerr << "\n\tpre_upload_homologous_regions_from_the_Reference_genome...\n";
	the_Event.pre_upload_homologous_regions_from_the_Reference_genome();
	
	std::cerr << "\n\tcalculate_GC_sensitive_per_position for Event pre-uploaded regions...\n";
	
	for (type_list_Coords_and_homologous_region::const_iterator it_homol = the_Event.regions_homologous_to_directly_affected_region.begin();
		it_homol != the_Event.regions_homologous_to_directly_affected_region.end();
		++it_homol)	     
	{
	    const BOOST_Interval homol_reg(get_endpoints_of_keys_of_compressed_map(it_homol->compressed_map_homologous_region_to_profile));
	    	    
	    type_map_uint_to_longdouble hap_frag_rate;
	    initialize_map_across_interval__via_insert<uint,longdouble>(hap_frag_rate, homol_reg, 0);
	    
	    calculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP
			(the_Event.uploaded_pieces_of_chromos,
			it_homol->chromosome_of_homologous_region,
			hap_frag_rate,
			0, 0, 0);	
			
	    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate.insert(std::pair<uint,type_map_uint_to_longdouble>(it_homol->chromosome_of_homologous_region, type_map_uint_to_longdouble()));
	    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate.at(it_homol->chromosome_of_homologous_region).insert(hap_frag_rate.begin(), hap_frag_rate.end());	    
	}//homol

	the_Event.uploaded_pieces_of_chromos.clear();
	
    }//raw homol reg rates
    
    
    
    
    
    
    

				
    
    
    
    
    
    
    std::cerr << "\nthe_Event.regions_homologous_to_directly_affected_region.size() = " 
	    <<  the_Event.regions_homologous_to_directly_affected_region.size()  << "\n";
    
	
    type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble  map_diploid_Event_outcome_to_absolute_position_to_expected_read_depths;
    
 
    type_map_uint_to_int contributing_locations_per_position_in_haploid_reference_genome;    
    initialize_map_across_interval__via_insert<uint,int>(contributing_locations_per_position_in_haploid_reference_genome, expanded_region_of_interest, 0);
    
    {//NULL
	std::cerr << "\t\t\tcreating NULL...\n";
	    
	type_map_uint_to_longdouble expected_fragmentation_sum___NULL;
	initialize_map_across_interval__via_insert<uint,longdouble>(expected_fragmentation_sum___NULL, expanded_region_of_interest, 0);
	
	appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
				    (the_Event,
				    hap_outcome__None,
				    empty_BI,
				    expected_fragmentation_rate___HAPLOID_NULL,
				    expanded_region_of_interest,
				    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
				    expected_fragmentation_sum___NULL,
				    contributing_locations_per_position_in_haploid_reference_genome);    
				    
	if (chromo_of_region < 23  or  (chromo_of_region == 23  and  gender_of_individual_is_female))  //diploid
	{
	    std::cerr << "\t\t\t\tmaking NULL diploid...\n";
	    multiply_values_of_map_by_constant<uint,longdouble>(expected_fragmentation_sum___NULL, 2);
	}
	
	map_diploid_Event_outcome_to_absolute_position_to_expected_read_depths[type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None)] 
					= divide_sum_by_contributions(expected_fragmentation_sum___NULL, contributing_locations_per_position_in_haploid_reference_genome);
    }//NULL
    
    {//default alternatives
    
	{//del
	    type_map_uint_to_longdouble expected_fragmentation_sum__diploid_alternative;
	    initialize_map_across_interval__via_insert<uint,longdouble>(expected_fragmentation_sum__diploid_alternative, expanded_region_of_interest, 0);
	    
	    type_map_uint_to_int dummy_contributions_alternative;    
	    initialize_map_across_interval__via_insert<uint,int>(dummy_contributions_alternative, expanded_region_of_interest, 0);
	    
	    appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
				    (the_Event,
				    hap_outcome__None,
				    empty_BI,
				    expected_fragmentation_rate___HAPLOID_NULL,
				    expanded_region_of_interest,
				    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
				    expected_fragmentation_sum__diploid_alternative,
				    dummy_contributions_alternative);   
				 
	    BOOST_Interval dummy_breakpoint;
	    if (!the_Event.uploaded_variational_positions_of_profile.empty())
			dummy_breakpoint.set(*the_Event.uploaded_variational_positions_of_profile.begin(), *the_Event.uploaded_variational_positions_of_profile.begin());
	    else
			dummy_breakpoint.set(the_Event.last_profile_nongap_index, the_Event.last_profile_nongap_index);
	    
	    appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
					(the_Event,
					hap_outcome__Del,
					dummy_breakpoint,
					expected_fragmentation_rate___HAPLOID_NULL,
					expanded_region_of_interest,
					map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
					expected_fragmentation_sum__diploid_alternative,
					dummy_contributions_alternative); 
					
	    map_diploid_Event_outcome_to_absolute_position_to_expected_read_depths[type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__Del)] 
	    = divide_sum_by_contributions(expected_fragmentation_sum__diploid_alternative, contributing_locations_per_position_in_haploid_reference_genome);    
	}//del
	
	{//dup
	    type_map_uint_to_longdouble expected_fragmentation_sum__diploid_alternative;
	    initialize_map_across_interval__via_insert<uint,longdouble>(expected_fragmentation_sum__diploid_alternative, expanded_region_of_interest, 0);
	    
	    type_map_uint_to_int dummy_contributions_alternative;    
	    initialize_map_across_interval__via_insert<uint,int>(dummy_contributions_alternative, expanded_region_of_interest, 0);
	    
	    appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
				    (the_Event,
				    hap_outcome__None,
				    empty_BI,
				    expected_fragmentation_rate___HAPLOID_NULL,
				    expanded_region_of_interest,
				    map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
				    expected_fragmentation_sum__diploid_alternative,
				    dummy_contributions_alternative);   
				 
	    BOOST_Interval dummy_breakpoint;
	    if (!the_Event.uploaded_variational_positions_of_profile.empty())
			dummy_breakpoint.set(*the_Event.uploaded_variational_positions_of_profile.begin(), *the_Event.uploaded_variational_positions_of_profile.begin());
	    else
			dummy_breakpoint.set(the_Event.last_profile_nongap_index, the_Event.last_profile_nongap_index);
	    
	    appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
					(the_Event,
					hap_outcome__Dup,
					dummy_breakpoint,
					expected_fragmentation_rate___HAPLOID_NULL,
					expanded_region_of_interest,
					map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
					expected_fragmentation_sum__diploid_alternative,
					dummy_contributions_alternative); 
					
	    map_diploid_Event_outcome_to_absolute_position_to_expected_read_depths[type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__Dup)] 
	    = divide_sum_by_contributions(expected_fragmentation_sum__diploid_alternative, contributing_locations_per_position_in_haploid_reference_genome);    
	}//dup
	
    }//default alternatives
    
    
    if (the_diploid_outcome != type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None))
    {//alternative
		std::cerr << "\t\t\tcreating alternative...\n";
		
		type_map_uint_to_longdouble expected_fragmentation_sum__diploid_alternative;
		initialize_map_across_interval__via_insert<uint,longdouble>(expected_fragmentation_sum__diploid_alternative, expanded_region_of_interest, 0);
		
		type_map_uint_to_int dummy_contributions_alternative;    
		initialize_map_across_interval__via_insert<uint,int>(dummy_contributions_alternative, expanded_region_of_interest, 0);
		
		appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
						(the_Event,
						the_diploid_outcome.first,
						profile_breakpoints.first,
						expected_fragmentation_rate___HAPLOID_NULL,
						expanded_region_of_interest,
						map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
						expected_fragmentation_sum__diploid_alternative,
						dummy_contributions_alternative);    
			
		appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
						(the_Event,
						the_diploid_outcome.second,
						profile_breakpoints.second,
						expected_fragmentation_rate___HAPLOID_NULL,
						expanded_region_of_interest,
						map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
						expected_fragmentation_sum__diploid_alternative,
						dummy_contributions_alternative);    
						
		map_diploid_Event_outcome_to_absolute_position_to_expected_read_depths[the_diploid_outcome] 
		= divide_sum_by_contributions(expected_fragmentation_sum__diploid_alternative, contributing_locations_per_position_in_haploid_reference_genome);
    
    }//alternative
            
    
    
    
    std::cerr << "\t\tDone with \"construct_pseudo_RD_expected_Read_Depth_counts_under_NULL_hypothesis_for_desired_display_region\"!!!!\n";
    
    return  map_diploid_Event_outcome_to_absolute_position_to_expected_read_depths;
    
}//construct_pseudo_RD_expected_Read_Depth_counts_under_NULL_hypothesis_for_desired_display_region








// void  display__GC_rate_sensitive__LCRmasked__Read_Depth_sliding_window__with__unique_padding__for_a_given_region_of_genome
// 		    (std::map<uint, Conn_Comp> &all_CCs,
// 		     const Pseudo_Read_Depth_display_container &pseudo_RD_display_region) 
// {	
//     
//     std::cerr << "\n\n inside \"display__GC_rate_sensitive__LCRmasked__Read_Depth_sliding_window__with__unique_padding__for_a_given_region_of_genome\"!!!!\n";
//     
//     const BOOST_Interval expanded_region_of_interest(  safe_subtract_base_1(pseudo_RD_display_region.affected_region_of_interest.lower(),
// 									    pseudo_RD_display_region.padding_amount_for_neighboring_unique_regions),
// 						       safe_chromo_add_base_1(pseudo_RD_display_region.affected_region_of_interest.upper(),
// 									      pseudo_RD_display_region.padding_amount_for_neighboring_unique_regions,
// 									      pseudo_RD_display_region.chromo_of_region)  );
// 	
//     std::cerr << "\n\t\tpseudo_RD_display_region.affected_region_of_interest = " << pseudo_RD_display_region.affected_region_of_interest << "\n\n";
//     std::cerr << "\n\t\texpanded_region_of_interest = " << expanded_region_of_interest << "\n\n";
// 			    
// 	
//     
// 			    
// 			    
//     std::cerr << "\t\t\tpre_uploaded_regions_across_all_chromosomes...\n\n";
//     //pre-upload relevant parts:
//     type_map_uint_to_list_BI__string pre_uploaded_regions_across_all_chromosomes;
// 		{    
// 		    //upload extra so there are no problems with border-case PERs.
// 		    const BOOST_Interval upload_region(  safe_subtract_base_1(expanded_region_of_interest.lower(), 30000),
// 							 safe_chromo_add_base_1(expanded_region_of_interest.upper(), 30000,
// 										pseudo_RD_display_region.chromo_of_region)  );
// 						
// 		    std::cerr << "\n\n\t\tpreuploading region:  " << upload_region << "\n\n";
// 		    
// 		    const std::string desired_sequence(
// 						get_specific_DNA_sequence_from_Reference_genome(
// 								pseudo_RD_display_region.chromo_of_region, 
// 								upload_region)  );                                                                                                 
// 												    
// 		    pre_uploaded_regions_across_all_chromosomes[ pseudo_RD_display_region.chromo_of_region ]
// 						.push_back( type_BI__string(upload_region, desired_sequence) );     
// 		}    
//     
//     
//     
//     
//     
//     
//     
//     
//     type_map_uint_to_longdouble positions_along_chromosome___NO_EVENT;
// 	    {
// 		type_map_uint_to_longdouble::iterator it_init = positions_along_chromosome___NO_EVENT.begin();
// 		for (uint pos = expanded_region_of_interest.lower();  pos <=  expanded_region_of_interest.upper();  ++pos)		
// 		    it_init = positions_along_chromosome___NO_EVENT.insert(it_init, type_uint__longdouble(pos, 0.00L)   );		
// 	    }//init
// 	                
//     
//     
//     std::cerr << "\t\t\tcalculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP...\n\n";
//     calculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP
// 		    (pre_uploaded_regions_across_all_chromosomes,
// 		     pseudo_RD_display_region.chromo_of_region,
// 		     positions_along_chromosome___NO_EVENT,
// 		     0, 0, 0);
// 	
//     if (  pseudo_RD_display_region.chromo_of_region < 23 
// 	 or  (pseudo_RD_display_region.chromo_of_region == 23  and  gender_of_individual_is_female)    )
//     {		     
// 	for (type_map_uint_to_longdouble::iterator it_no = positions_along_chromosome___NO_EVENT.begin();
// 		it_no != positions_along_chromosome___NO_EVENT.end();
// 		++it_no)
// 	{
// 	    it_no->second *= 2.00L;    
// 	}//make diploid
//     }
//     else if (pseudo_RD_display_region.chromo_of_region == 24  and  gender_of_individual_is_female)
//     {
// 	std::stringstream error_strm;
// 	print_line_of_markers("ERROR! ", &error_strm);
// 	print_line_of_markers("(", &error_strm);
// 	error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";	
// 	error_strm << "\n\nERROR!    pseudo_RD_display_region.chromo_of_region == 24  and  gender_of_individual_is_female  !!!!!!!!!   THIS IS IMPOSSIBLE!!!!\n\n";	
// 	print_line_of_markers(")", &error_strm);
// 	std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
// 	
// 	for (type_map_uint_to_longdouble::iterator it_no = positions_along_chromosome___NO_EVENT.begin();
// 		it_no != positions_along_chromosome___NO_EVENT.end();
// 		++it_no)
// 	{
// 	    it_no->second = 0.00L;    
// 	}//make diploid    
//     }	
//     
// 	
// 	
// 	
// 		    
//     type_map_uint_to_longdouble   positions_along_chromosome___called_Event_outcome(positions_along_chromosome___NO_EVENT);
//     if (  pseudo_RD_display_region.chromo_of_region < 23 
// 	 or  (pseudo_RD_display_region.chromo_of_region == 23  and  gender_of_individual_is_female)    )
//     {//adjust    
// 	longdouble mult_factor;
// 	if (pseudo_RD_display_region.call_type == 1) // hap del
// 	    mult_factor = 0.5L;
// 	else if (pseudo_RD_display_region.call_type == 11) // diploid del
// 	    mult_factor = 0.00L;
// 	else if (pseudo_RD_display_region.call_type == 2)  //hap dup
// 	    mult_factor = 1.5L;
// 	else if (pseudo_RD_display_region.call_type == 22)  //diploid dup
// 	    mult_factor = 2.00L;
// 	else 
// 	{
// 	    std::stringstream error_strm;
// 	    print_line_of_markers("ERROR! ", &error_strm);
// 	    print_line_of_markers("(", &error_strm);
// 	    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
// 	    
// 	    error_strm << "\n\nERROR!   unrecognize call type = "
// 		    << pseudo_RD_display_region.call_type 
// 		    <<  "   in  \"display__GC_rate_sensitive__LCRmasked__Read_Depth_sliding_window__with__unique_padding__for_a_given_region_of_genome\"!!!\n\n";
// 	    
// 	    print_line_of_markers(")", &error_strm);
// 	    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );	    
// 	    exit(1);
// 	}
// 
// 
// 	for (type_map_uint_to_longdouble::iterator it_pos = positions_along_chromosome___called_Event_outcome.begin();
// 		it_pos != positions_along_chromosome___called_Event_outcome.end();
// 		++it_pos)
// 	    it_pos->second *= mult_factor;
//     }//adjust
//     else
//     {//X and male, or Y and M/F
// 	longdouble mult_factor;
// 	if (pseudo_RD_display_region.call_type == 1) // hap del
// 	    mult_factor = 0.00L;
// 	else if (pseudo_RD_display_region.call_type == 2)  //hap dup
// 	    mult_factor = 2.00L;
// 	else
// 	{
// 	    std::stringstream error_strm;
// 	    print_line_of_markers("ERROR! ", &error_strm);
// 	    print_line_of_markers("(", &error_strm);
// 	    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
// 	    
// 	    error_strm << "\n\nERROR!   invalid call type = "
// 		    << pseudo_RD_display_region.call_type 
// 		    << "\n on chromo " << pseudo_RD_display_region.chromo_of_region
// 		    << "\n\t with gender = " << gender_of_individual_is_female
// 		    << "\n"
// 		    <<  "   in  \"display__GC_rate_sensitive__LCRmasked__Read_Depth_sliding_window__with__unique_padding__for_a_given_region_of_genome\"!!!\n\n";
// 	    
// 	    print_line_of_markers(")", &error_strm);
// 	    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );	    
// 	    exit(1);
// 	}
//     }//X and male, or Y and M/F
//     
//     
//     
//     
// 		    
// 		    
// 		    
// 		    
//     type_map_uint_to_uint  counts_for_these_positions;
// 	    {
// 		type_map_uint_to_uint::iterator it_init = counts_for_these_positions.begin();
// 		for (type_map_uint_to_longdouble::const_iterator it_pos = positions_along_chromosome___called_Event_outcome.begin();
// 			it_pos != positions_along_chromosome___called_Event_outcome.end();
// 			++it_pos)		
// 		    it_init = counts_for_these_positions.insert(it_init, type_uint__uint(it_pos->first, 0)   );
// 		
// 	    }//init  
// 	    
//     std::cerr << "\t\t\tcount_reads_whose_left_endpoint_falls_at_position...\n\n";
//     {//count 	
// 	const uint number_of_positions_to_count = counts_for_these_positions.size();
// 	
// 	std::vector< type_map_uint_to_uint::iterator > positions_to_count____loop;
// 	positions_to_count____loop.reserve( number_of_positions_to_count );         
// 						    {//scope
// 							for (type_map_uint_to_uint::iterator it_pos = counts_for_these_positions.begin();
// 								it_pos != counts_for_these_positions.end();
// 								++it_pos)
// 							    positions_to_count____loop.push_back( it_pos );
// 						    }//scope
// 										
// // 	#pragma omp parallel  shared(count_for_this_chromosome)
// // 	{//parallel
// 	    BamTools::BamReader *const ptr_to_thread_BAM_Reader = &(*(my_BAM_readers.begin() + (uint)omp_get_thread_num())); //ptr arithmetic            
// 	    
// // 	    #pragma omp for schedule(dynamic,100)  reduction(+:count_for_this_chromosome)
// 	    for (int j = 0; j < number_of_positions_to_count; ++j)
// 		positions_to_count____loop[j]->second = count_reads_whose_left_endpoint_falls_at_position(
// 											pseudo_RD_display_region.chromo_of_region,
// 											positions_to_count____loop[j]->first,
// 											*ptr_to_thread_BAM_Reader);       
// 	    mpfr_free_cache();                                                                                                                                                                                         
// // 	}//parallel                 
//     }//count		    
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
//     
//     
//     
//     {//mask non-unique regions
// 	std::cerr << "\t\t\tidentifying nonuniqueness in genome...\n\n";
// 	//non-uniqueness
// 	const type_map_uint_to_list_BI non_uniq_regions_of_genome(  load_Non_unique_regions_of_genome()  );	
// 	
// 	type_map_uint_to_list_BI::const_iterator it_chr_non = non_uniq_regions_of_genome.find(  pseudo_RD_display_region.chromo_of_region  );	
// 	if (it_chr_non != non_uniq_regions_of_genome.end() )
// 	{
// 	    for (type_list_BI::const_iterator it_non = it_chr_non->second.begin();
// 		    it_non != it_chr_non->second.end();
// 		    ++it_non)
// 	    {
// 		for (uint pos = it_non->lower(); pos <= it_non->upper(); ++pos)
// 		{
// 		    type_map_uint_to_longdouble::iterator it_gc = positions_along_chromosome___NO_EVENT.find(pos);
// 		    if (it_gc != positions_along_chromosome___NO_EVENT.end())
// 		    {
// 			it_gc->second = -99999.99L; // dummy value
// 		    }
// 		}
// 	    }		
// 	}//found chr	
//     }//mask non-unique regions
//     
//     
//     
//     
//     
//     
//     {//RD hypothesis test    
// 	longdouble poisson_parameter__null_hypothesis = 0.00L;
// 	longdouble poisson_parameter__alternative_hypothesis = 0.00L;
// 	
// 	std::map<longdouble, uint> unique_pos_along_chromosome___NO_EVENT__flipped(get_inverse_of_map__temp<uint,longdouble>(positions_along_chromosome___NO_EVENT));
// 	unique_pos_along_chromosome___NO_EVENT__flipped.erase(1.00L);
// 	
// 	const type_set_uint acceptable_positions(
// 		extract_keys_of_map_and_return_as_set<uint,longdouble>(
// 		    get_map_intersection_of_map_keys_and_interval<uint,longdouble>(
// 			get_inverse_of_map__temp<longdouble, uint>(unique_pos_along_chromosome___NO_EVENT__flipped),
// 			pseudo_RD_display_region.affected_region_of_interest)));	
// 	
// 	const longdouble poisson_mean__null_hypothesis	
// 		=  sum_over_map_values<uint,longdouble>(
// 			get_map_intersection_of_map_keys_and_set<uint,longdouble>(positions_along_chromosome___NO_EVENT, acceptable_positions));
// 				
// 	const longdouble poisson_mean__alternative_hypothesis	
// 		=  sum_over_map_values<uint,longdouble>(
// 			get_map_intersection_of_map_keys_and_set<uint,longdouble>(positions_along_chromosome___called_Event_outcome, acceptable_positions));
// 				
// 	const uint observed_data__read_count_inside_called_region	
// 		=  sum_over_map_values<uint,uint>(
// 			get_map_intersection_of_map_keys_and_set<uint, uint>(counts_for_these_positions, acceptable_positions));
// 	
// 	//poisson:
// 	boost::math::poisson_distribution<mpfr_class>  my_poisson__null( (double)poisson_mean__null_hypothesis + 0.0000001);
// 	boost::math::poisson_distribution<mpfr_class>  my_poisson__alternative( (double)poisson_mean__alternative_hypothesis + 0.0000001);
// 	
// 	const real P_observed_cond_null__poisson(boost::math::pdf<mpfr_class>(my_poisson__null, 
// 									      mpfr_class(observed_data__read_count_inside_called_region)).get_mpfr_t());
// 	const real P_observed_cond_alternative__poisson(boost::math::pdf<mpfr_class>(my_poisson__alternative,
// 										     mpfr_class(observed_data__read_count_inside_called_region)).get_mpfr_t());
// 	
// 	
// 	
// 	
// 	//negative binomial:   copied from   Event::calculate_probability_of_observed_read_depth_given_expected_read_depth_over_natural_poisson_intervals   
// 	
// 	const mpfr_class the_number_one(1.00);		
// 	const mpfr_class expected_lambda__null((double)poisson_mean__null_hypothesis + 0.0000001);
// 	const mpfr_class expected_lambda__alternative((double)poisson_mean__alternative_hypothesis + 0.0000001);
// 		
//         // CAREFULLY DERIVED VERSION:                                       
//         const mpfr_class parameter_BOOST_r(the_number_one / (double)RD_variance_scale___alpha_squared);
//                 //BOOST_r  <==>  my_k  <==> Gamma shape parameter k.                                                                
//         const mpfr_class parameter_BOOST_p___null(the_number_one / (the_number_one + (expected_lambda__null*(double)RD_variance_scale___alpha_squared)));
// 	const mpfr_class parameter_BOOST_p___alternative(the_number_one / (the_number_one + (expected_lambda__alternative*(double)RD_variance_scale___alpha_squared)));
//                 
//         boost::math::negative_binomial_distribution<mpfr_class> my_negative_binomial__null(parameter_BOOST_r, parameter_BOOST_p___null);
// 	boost::math::negative_binomial_distribution<mpfr_class> my_negative_binomial__alternative(parameter_BOOST_r, parameter_BOOST_p___alternative);
// 	
// 	const real P_observed_cond_null__negativebinomial(
// 			boost::math::pdf<mpfr_class>(my_negative_binomial__null, mpfr_class(observed_data__read_count_inside_called_region)).get_mpfr_t()); 
// 	const real P_observed_cond_alternative__negativebinomial(
// 			boost::math::pdf<mpfr_class>(my_negative_binomial__alternative, mpfr_class(observed_data__read_count_inside_called_region)).get_mpfr_t());
// 				
// 	
// 	std::cerr << "\n\nobserved_data__read_count_inside_called_region = " << observed_data__read_count_inside_called_region
// 		    << "\npoisson_mean__null_hypothesis = " << poisson_mean__null_hypothesis
// 		    << "\npoisson_mean__alternative_hypothesis = " << poisson_mean__alternative_hypothesis
// 		    << "\nPoisson( data  ;  null ) = " << P_observed_cond_null__poisson.toString(15)
// 		    << "\nPoisson( data  ;  alternative ) = " << P_observed_cond_alternative__poisson.toString(15)
// 		    << "\n\nNegativeBinomial( data  ;  null ) = " << P_observed_cond_null__negativebinomial.toString(15)
// 		    << "\nNegativeBinomial( data  ;  alternative ) = " << P_observed_cond_alternative__negativebinomial.toString(15)		    
// 		    << "\n\n";
// 		    
// 	{//save to file			    	
// 	    char out_RD_hypothesis_tests[option_size];	
// 	    std::sprintf(out_RD_hypothesis_tests, "%s/Read_depth_for_unique_regions___chr__%u__%u__%u__call__%u__%s__%s__tests", 
// 					pseudo_RD_display_region.some_dir_for_saving_data.c_str(),
// 					pseudo_RD_display_region.chromo_of_region,
// 					pseudo_RD_display_region.affected_region_of_interest.lower(),
// 					pseudo_RD_display_region.affected_region_of_interest.upper(),
// 					pseudo_RD_display_region.call_type,
// 					genome_name.c_str(),
// 					population_abbreviation.c_str()  );
// 					
// 	    std::ofstream  outfs(out_RD_hypothesis_tests);
// 				if (!outfs.good())
// 				{
// 				    std::stringstream error_strm;
// 				    print_line_of_markers("ERROR! ", &error_strm);
// 				    print_line_of_markers("(", &error_strm);
// 				    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";	    
// 				    error_strm << "\n\nERROR!  unable to open file out_RD_hypothesis_tests = ["  <<   out_RD_hypothesis_tests
// 						<< "] in \"display__GC_rate_sensitive__LCRmasked__Read_Depth_sliding_window__with__unique_padding__for_a_given_region_of_genome\".\n\n";	    
// 				    print_line_of_markers(")", &error_strm);
// 				    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() ); 	
// 				    exit(1);
// 				}
// 				
// 	    outfs << observed_data__read_count_inside_called_region
// 		  << "\n" << poisson_mean__null_hypothesis
// 		  << "\n" << poisson_mean__alternative_hypothesis
// 		  << "\n" << P_observed_cond_null__poisson.toString(10)
// 		  << "\n" << P_observed_cond_alternative__poisson.toString(10)
// 		  << "\n" << P_observed_cond_null__negativebinomial.toString(10)
// 		  << "\n" << P_observed_cond_alternative__negativebinomial.toString(10);
// 		  
// 	    outfs.close();
// 	}//save to file    
//     
//     }//RD hypothesis test
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
//     std::cerr << "\t\t\twrite to stream...\n\n";
//     //write to stream
//     std::stringstream out_ss;
//     out_ss.precision(15);
//     
//     {//write to stream
// 	type_map_uint_to_longdouble::const_iterator it_gc_NO = positions_along_chromosome___NO_EVENT.begin();
// 	type_map_uint_to_longdouble::const_iterator it_gc_yes = positions_along_chromosome___called_Event_outcome.begin();
// 	type_map_uint_to_uint::const_iterator it_count = counts_for_these_positions.begin();    
// 	
// 	
// 	while (  it_gc_NO != positions_along_chromosome___NO_EVENT.end()  )
// 	{
// 	    if (it_gc_NO->second > -1.00L)  // safeguard against dummy values
// 	    {
// 		out_ss 
// 		    << it_gc_NO->first << "\t" 
// 		    << it_gc_NO->second << "\t"
// 		    << it_gc_yes->second << "\t"
// 		    << it_count->second << "\n";
// 	    }
// 	    else
// 	    {
// 		out_ss 
// 		    << it_gc_NO->first << "\tNaN\tNaN\tNaN\n";	    	
// 	    }	
// 	    
// 	    ++it_gc_NO;
// 	    ++it_gc_yes;
// 	    ++it_count;
// 	}//pos    	
//     }//write to stream
//     
//     
//     
//     
//     std::cerr << "\t\t\tsave to file...\n\n";
//     {//save to file		  
// 	char out_RD_fname[option_size];	
// 	std::sprintf(out_RD_fname, "%s/Read_depth_for_unique_regions___chr__%u__%u__%u__call__%u__%s__%s", 
// 				    pseudo_RD_display_region.some_dir_for_saving_data.c_str(),
// 				    pseudo_RD_display_region.chromo_of_region,
// 				    pseudo_RD_display_region.affected_region_of_interest.lower(),
// 				    pseudo_RD_display_region.affected_region_of_interest.upper(),
// 				    pseudo_RD_display_region.call_type,
// 				    genome_name.c_str(),
// 				    population_abbreviation.c_str()  );
// 				    
// 	FILE *outfile = std::fopen(out_RD_fname, "w");
// 			    if (outfile == NULL)
// 			    {
// 				std::stringstream error_strm;
// 				print_line_of_markers("ERROR! ", &error_strm);
// 				print_line_of_markers("(", &error_strm);
// 				error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";	    
// 				error_strm << "\n\nERROR!  unable to open file out_RD_fname = ["  <<   out_RD_fname
// 					    << "] in \"display__GC_rate_sensitive__LCRmasked__Read_Depth_sliding_window__with__unique_padding__for_a_given_region_of_genome\".\n\n";	    
// 				print_line_of_markers(")", &error_strm);
// 				std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() ); 	
// 				exit(1);
// 			    }		
// 		    
// 	std::fprintf(outfile, "%s", out_ss.str().c_str());
// 	
// 	std::fclose(outfile);		    
//     }//save to file          
// 		    
// 		    
// 		    
// } // display__GC_rate_sensitive__LCRmasked__Read_Depth_sliding_window__with__unique_padding__for_a_given_region_of_genome
// 
// 



































void prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population
		    (const std::string &BAM_filename,
		     const std::string &gender_and_population_filename,
		    BamTools::BamReader &my_BAM_reader)
{
    
    std::cerr << "\n\ninside \"prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population\"\n\nmax_number_of_OMP_threads = " << max_number_of_OMP_threads << "\n\nBAM_filename = [" << BAM_filename << "]\n\n";
    

    const bool success_open_BAM = my_BAM_reader.Open(BAM_filename);  //sorted by position    

    if (!success_open_BAM)
    {        
	std::cerr << "\nERROR: unable to open BAM file [" << BAM_filename << "]  in \"prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population\".\n\n";   
	std::cerr.flush();
	exit(1);
    }
    else
    {  std::cerr << "successfully opened    BAM_index_filename  =   ["  << BAM_filename << "]\n\n";  }
    
    	    
    const bool success_open_index = my_BAM_reader.LocateIndex();
    
    if (!success_open_index)
    {
	bool create_indx_success;
	
	world_MPI_boost_ptr->barrier();
	if (my_MPI_rank == 0)
	{      
	    std::fprintf(stderr, "ERROR: unable to Locate BAM file index [%s].\n\n", BAM_filename.c_str());
	    
	    std::fprintf(stderr, "creating index, this may take a while...\n");
	    create_indx_success = my_BAM_reader.CreateIndex();
	}	
	world_MPI_boost_ptr->barrier();
		    
	if (my_MPI_rank > 0)
	{  create_indx_success = my_BAM_reader.LocateIndex();  }//OpenIndex(BAM_filename);//LocateIndex();
	
	if (!create_indx_success)
	{
	    print_line_of_markers("ERROR! ");
	    print_line_of_markers("(");
	    std::fprintf(stderr, "ERROR: unable to build index for BAM file index [%s].\n\n", BAM_filename.c_str());
	    print_line_of_markers(")");    
	    exit(1);
	}                        
    }    
    
    
    world_MPI_boost_ptr->barrier();
    



    if (my_MPI_rank == 0)
    {  std::fprintf(stderr, "BAMfile and index successfully loaded.\n\n");  }
    






    
    
    set_map_from_chromosome_values_to_RefID_in_BAM_file(my_BAM_reader);    
    if (my_MPI_rank == 0)
    {  print_map_keys_and_values<uint, uint>(map_chromosome_value_to_BAM_Ref_IDs, "map_chromosome_value_to_BAM_Ref_IDs");  }
    
    
            


    
    
// bas
    read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda(BAM_filename, true, &my_BAM_reader);
//    DEBUG simulation change by MATT  .  uncomment above line to remove the debug change!!!!!

	    

    
    world_MPI_boost_ptr->barrier();       
    
    
    
    if (!gender_and_population_filename.empty())
    {  read_gender_and_population_from_file(gender_and_population_filename);  }
    

    
} // prepare_BAM_files_for_IO_operations__and__load_BAS_data



















void read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda
		    (const std::string &BAM_filename,
		    const bool &display_info,
		    BamTools::BamReader *const &my_BAM_reader_ptr)		    
		    
{
    
    if (my_MPI_rank == 0)
    {  std::fprintf(stderr, "\nreading BAS file...\n");  }       
    
    Readgroup_stats = Readgroup_statistics::read_bas_file(BAM_filename.c_str());
    
    
/*  DEBUG BY MATT BEGIN
    {//eliminate those of size 0  -  we ONLY  consider Paired end Reads.
	type_map_string_to_Readgroup_stats::iterator it_rg =  Readgroup_stats.begin();
	while (it_rg != Readgroup_stats.end())
	{
	    if (it_rg->second.mean_insert_size == 0   or    it_rg->second.number_mapped_paired_reads == 0)
		Readgroup_stats.erase(it_rg++);   //post-increment necessary!
	    else
		++it_rg;       
	}
	    
	    
	if (  Readgroup_stats.empty()  )
	{
	    std::stringstream error_strm;
	    print_line_of_markers("ERROR! ", &error_strm);
	    print_line_of_markers("(", &error_strm);
	    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
	    
	    error_strm << "\n\nERROR!    Readgroup_stats   is empty after removing   readgroups with  mean_insert_size == 0   or   number_mapped_paired_reads == 0\n\n";
	    
	    print_line_of_markers(")", &error_strm);
	    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
	    exit(1);                
	}
    }//eliminate those of size 0  -  we ONLY  consider Paired end Reads.        
  DEBUG BY MATT END      */
	
	
    
    if (my_MPI_rank == 0)
    {  std::fprintf(stderr, "\nreading GC count rates...\n");  }
    
    for (type_map_string_to_Readgroup_stats::iterator it_rg =  Readgroup_stats.begin();
	    it_rg != Readgroup_stats.end();
	    ++it_rg)
    {
	it_rg->second.median_insert_size = 299; // DEBUG BY MATT
	it_rg->second.read_diploid_GC_rates_from_file__and_make_them_haploid_rates(BAM_filename.c_str());	
	
	bool has_nonpositive_values = false;	    
	for (type_map_uint_to_longdouble::const_iterator it_gc = it_rg->second.haploid_fragmentation_rate__per_GC_count.begin();
	    it_gc != it_rg->second.haploid_fragmentation_rate__per_GC_count.end();
	    ++it_gc)	    
	{
	    if (it_gc->second <= 0)
	    {
		has_nonpositive_values = true;
		break;
	    }
	}		
	
	if (has_nonpositive_values)	
	{//adjust negative values
	    std::cerr << "adjust negative values"; 
	    
	    longdouble smallest_pos_value = std::numeric_limits<longdouble>::max();
	
	    for (type_map_uint_to_longdouble::const_iterator it_gc = it_rg->second.haploid_fragmentation_rate__per_GC_count.begin();
		it_gc != it_rg->second.haploid_fragmentation_rate__per_GC_count.end();
		++it_gc)
	    {
		if (it_gc->second > 0 and  it_gc->second < smallest_pos_value)
		{  smallest_pos_value = it_gc->second;  }	    	    	    
	    }	    
	    
	    for (type_map_uint_to_longdouble::iterator it_gc = it_rg->second.haploid_fragmentation_rate__per_GC_count.begin();
		it_gc != it_rg->second.haploid_fragmentation_rate__per_GC_count.end();
		++it_gc)
	    {
		if (it_gc->second <= 0)
		{  it_gc->second = smallest_pos_value/10;  }
	    }	    	    
	}//adjust negative values
	
	
	    
	if (my_MPI_rank == 0 and display_info)
	{
	    std::stringstream out_strm;
	    out_strm << "\n\n" << it_rg->first << "\n";
	    print_map_keys_and_values<uint, longdouble>(it_rg->second.haploid_fragmentation_rate__per_GC_count, "haploid_fragmentation_rate__per_GC_count", &out_strm);	    
	    std::cerr << out_strm.str() << "\n\n";
	}			
	
    }//it_rg
    
    
    
    
    
    
    
    
    int total_number_of_mapped_Proper_paired_MATES = 0;
    maximum_frag_length_absolute_dev = 1;
    
    haploid_coverage = 0.00L;     // this is the coverage for EACH haploid separately.  So if at location x, there is a SNP which is "A" in the first haploid and "G" in the second haploid, then we expect to see "haploid_coverage" paired-end reads which show "A" and "haploid_coverage" paired-end reas which show "G".

    for (type_map_string_to_Readgroup_stats::const_iterator it_rg =  Readgroup_stats.begin();
	    it_rg != Readgroup_stats.end();
	    ++it_rg)
    {
	total_number_of_mapped_Proper_paired_MATES += it_rg->second.number_mapped_Proper_Pairs_reads;	    

	if ( (uint)it_rg->second.median_insert_size   >  maximum_average_fragment_length)
	{  maximum_average_fragment_length = (uint)it_rg->second.median_insert_size;  }
		    
	if (it_rg->second.median_insert_size____absolute_standard_deviation > maximum_frag_length_absolute_dev)
	{  maximum_frag_length_absolute_dev = it_rg->second.median_insert_size____absolute_standard_deviation;  }
	
	
	longdouble lambda_0___group = it_rg->second.number_mapped_Proper_Pairs_reads;  // diploid
	lambda_0___group /= 2;  // we want paired-reads, not mates.
	lambda_0___group /= (length_of_haploid_genome*2);   //this is the contribution to lambda_0 just from this read group.	
	
	haploid_coverage +=  (lambda_0___group * it_rg->second.median_insert_size);
	
// 	coverage_by_RG[it_rg->first] =  lambda_0___group * it_rg->second.median_insert_size;
    }
    
    
    //This works too.        
    //lambda_0 = total_number_of_mapped_Paired_reads;
    lambda_0 = total_number_of_mapped_Proper_paired_MATES;
    lambda_0 /= 2;  //because we want paired end reads, not mates.
    lambda_0 /= length_of_haploid_genome;
    
    
    std::cerr << "\n\ncalculating mate lengths...\n\n";
    for (type_map_string_to_Readgroup_stats::iterator it_rg = Readgroup_stats.begin();
	    it_rg != Readgroup_stats.end();
	    ++it_rg)
    {
	it_rg->second.average_mate_length = 37;
    }
// 	if (it_rg->second.median_insert_size == 0  or   it_rg->second.number_mapped_Proper_Pairs_reads  < 1000)
// 	{  continue;  }
// 	
// 	std::cerr << "\n\ncalculating mate lengths  for  it_rg = " << it_rg->first << "...\n\n";
// 	
// 	my_BAM_reader_ptr->Jump(map_chromosome_value_to_BAM_Ref_IDs.at(1), 152769227); // any old location.
// 	uint number_mates_from_this_RG = 0;
// 	uint sum_mate_lengths = 0;
// 	
// 	uint number_attempts = 0;
// 	BamTools::BamAlignment Balgmt;
// 	while (my_BAM_reader_ptr->GetNextAlignmentCore(Balgmt) and number_mates_from_this_RG < 2  and  number_attempts < 500)
// 	{
// 	    ++number_attempts;
// 	    const type_map_string_to_Readgroup_stats::const_iterator it_rg_found = identify_Readgroup_from_name_of_PER(Balgmt.Name);
// 	    if (it_rg_found != Readgroup_stats.end()
// 		and  it_rg_found->first.compare(it_rg->first) == 0)
// 	    {
// 		++number_mates_from_this_RG;
// 		sum_mate_lengths += Balgmt.Length;	    
// 	    }	
// 	}
// 	
// 	if (number_mates_from_this_RG < 2 and number_attempts >= 500)
// 	{
// 	    std::cerr << "\nunable to estimate mate length for readgroup [" << it_rg->first << "].  giving default = 76\n";
// 	     it_rg->second.average_mate_length = 76;	    
// 	}
// 	else
// 	{
// 	    it_rg->second.average_mate_length = (int)std::floor((double)sum_mate_lengths/number_mates_from_this_RG);
// 	    std::cerr << "\tReadgroup [" << it_rg->second.readgroup_name << "], average mate length = " << it_rg->second.average_mate_length << "\n";
// 	}
//     }//rg
//     my_BAM_reader_ptr->Rewind();
    
	    
    
    if (my_MPI_rank == 0 and display_info)
    {
	std::cerr << "\n\n\n\n\nReadGroup summary:\n\t\ttotal_number_of_mapped_Proper_paired_MATES = " << total_number_of_mapped_Proper_paired_MATES
		<< "\n\t\tmaximum_average_fragment_length = " << maximum_average_fragment_length
		<< "\n\t\tmaximum_frag_length_standard_dev = " << maximum_frag_length_absolute_dev
		<< "\n\t\tlambda_0 = " << lambda_0
		<< "\n\t\thaploid_coverage = " << haploid_coverage
		<< "\n\n\n\n\n";
		    
		    
		    
	std::fprintf(stderr, "\n\n\nReadgroup stats:\n\n");
	for (type_map_string_to_Readgroup_stats::const_iterator it_rg =  Readgroup_stats.begin();
		it_rg != Readgroup_stats.end();
		++it_rg)              
	{
	    std::fprintf(stderr, "group_name = %s\n\t\tmean_insert_size = %d\n\t\tmean_insert_size____standard_deviation = %d\n\t\tmedian_insert_size = %d\n\t\tmedian_insert_size____absolute_standard_deviation = %d\n\t\tnumber_mapped_paired_reads = %d\n\t\tnumber_mapped_Proper_Pairs_reads = %d",
			it_rg->second.readgroup_name.c_str(),
			it_rg->second.mean_insert_size,
			it_rg->second.mean_insert_size____standard_deviation,
			it_rg->second.median_insert_size,
			it_rg->second.median_insert_size____absolute_standard_deviation,
			it_rg->second.number_mapped_paired_reads,
			it_rg->second.number_mapped_Proper_Pairs_reads   );
	    std::cerr << "\n\t\taverage_mate_length = " <<  it_rg->second.average_mate_length << "\n\n\n";
			
	}
	std::fprintf(stderr, "\n\n\n\n\n\n");			
    }//0
    
    if (maximum_average_fragment_length == 0)
    {
	std::stringstream error_strm;
	print_line_of_markers("ERROR! ", &error_strm);
	print_line_of_markers("(", &error_strm);
	error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
	
	error_strm << "\n\nERROR!   maximum_average_fragment_length == 0 !!!!!\n\n";
	
	print_line_of_markers(")", &error_strm);
	std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
	exit(1);
    }    
    
}//read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda























uint  calculate_NON_uniqueness_of_region
	    (const uint &chr_of_reg,
	     const BOOST_Interval &the_reg,
	     const type_map_uint_to_list_BI &NON_unique_regions_of_genome)
{
    uint amount_of_affected_region_that_is_NON_unique = 0;
    
    const type_map_uint_to_list_BI::const_iterator  it_nonuniq_chr = NON_unique_regions_of_genome.find(  chr_of_reg  );	    
    if (it_nonuniq_chr != NON_unique_regions_of_genome.end())
    {
	for (type_list_BI::const_iterator it_bi = it_nonuniq_chr->second.begin();
		it_bi != it_nonuniq_chr->second.end();
		++it_bi)
	{
	    amount_of_affected_region_that_is_NON_unique += BOOST_width_inclusive(BOOST_intersect(*it_bi, the_reg));	    
	}	
    }//found chr
        
    
    return amount_of_affected_region_that_is_NON_unique;

} // calculate_NON_uniqueness_of_region










































uint  copy_contents_from_previous_output_directory_with_given_jobID_to_current_output_directory__and__return_last_CC_analyzed
			(const std::string  &previous_jobID)
{
            
    std::string output_root_dir(output_dir);
    std::string scratch_root_dir(scratch_job_output_dir);
    {
	const std::string current_job_id = output_dir.substr( output_dir.find_last_of('/')+1 );
	output_root_dir.erase(output_dir.find_last_of('/'),
				current_job_id.size()+1);    
	scratch_root_dir.erase(scratch_job_output_dir.find_last_of('/'),
				current_job_id.size()+1);				
    }
    
    
    std::string previous_job_output_dir(output_root_dir);
    previous_job_output_dir.push_back('/');
    previous_job_output_dir.append(previous_jobID);
    
    std::cerr << "\n\ninside \"copy_contents_from_previous_output_directory_with_given_jobID_to_current_output_directory__and__return_last_CC_analyzed\"\n\n";
    std::cerr << "copying from previous jobID = [" << previous_jobID << "],  old output dir = [" << previous_job_output_dir << "]\n\n\n";
    
    const boost::filesystem::recursive_directory_iterator  it_dir_END;
    
    for (boost::filesystem::recursive_directory_iterator it_old_dir(previous_job_output_dir);
	    it_old_dir != it_dir_END;
	    ++it_old_dir)
    {
	if (  ! boost::filesystem::is_directory((*it_old_dir).path())  )
	{                 
	    
	    std::string old_extension(  (*it_old_dir).path().string()  );
	    old_extension.assign(  old_extension.substr(old_extension.find(previous_jobID) + previous_jobID.size())  );
	    
	    std::string dest_path(   output_dir  );
	    dest_path.append(  old_extension );
	    
	    boost::filesystem::copy(   (*it_old_dir).path(),
					boost::filesystem::path(dest_path)   );                         
	}  		    
    }//it_old_dir
    
    
    
    std::string previous_scratch_lastanalyzed(scratch_root_dir);    
    previous_scratch_lastanalyzed.push_back('/');
    previous_scratch_lastanalyzed.append(previous_jobID);
    previous_scratch_lastanalyzed.push_back('/');
    previous_scratch_lastanalyzed.append("last_cc_analyzed");
    
    
    std::cerr << "\n\n\t\tDetecting last analyzed CC_ID from scratch file [" << previous_scratch_lastanalyzed << "]\n";   
    
    
    uint last_cc_analyzed;
    
    std::ifstream  ifs_lastcc(previous_scratch_lastanalyzed);
    if ( !ifs_lastcc.good() )
    {
	std::stringstream error_strm;
	print_line_of_markers("ERROR! ", &error_strm);
	print_line_of_markers("(", &error_strm);
	error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
	
	error_strm << "\n\nERROR!   failed to open ifstream   filename = ["  << previous_scratch_lastanalyzed   << "]\n\n";
	
	print_line_of_markers(")", &error_strm);
	std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
    
    }
    ifs_lastcc >> last_cc_analyzed;
    ifs_lastcc.close();
    
    std::cerr << "\n\t\tlast_cc_analyzed = " << last_cc_analyzed << "\n\n";
    
    
    std::cerr << "\n\nDONE with \"copy_contents_from_previous_directory_with_given_jobID__and__return_last_CC_analyzed\"!!!!\n\n";
    
    
    return last_cc_analyzed;

}//copy_contents_from_previous_directory_with_given_jobID__and__return_last_CC_analyzed











uint safe_subtract_base_1
	(const uint &some_value__base_1,
	 const uint &subtract_this_amount)
{
    return safe_subtract_base_1(some_value__base_1, (int)subtract_this_amount);
}


uint safe_chromo_add_base_1
	(const uint &some_value__base_1,
	 const uint &add_this_amount,
	 const uint &chromo_val)
{
    return safe_chromo_add_base_1(some_value__base_1, (int)add_this_amount, chromo_val);
}




uint safe_subtract_base_1
	(const uint &some_value__base_1,
	 const int  &subtract_this_amount)
{
    return (uint)std::max<int>(((int)some_value__base_1) - subtract_this_amount, 1);
}



uint safe_chromo_add_base_1
	(const uint &some_value__base_1,
	 const int &add_this_amount,
	 const uint &chromo_val)
{
    return  (uint)std::min<int>(((int)some_value__base_1) + add_this_amount, (int)Event::chromosome_lengths.at(chromo_val));
}

















    
    
void check_informativeness_of_every_potential_NAHR_breakpoint
	    (type_list__uint__uint  &informativeness_checking_parameters,
	     const std::map<uint, Conn_Comp> &Conn_Comps,
	     const std::string &some_outputdir,
	     const std::string &BAM_dir)
{
    std::cerr << "\n\n\n\n\tchecking overlapping Event LCRs...\n";
    for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
	    it_cc != Conn_Comps.end();
	    ++it_cc)
    {
	for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
		it_ev != it_cc->second.events.end();
		++it_ev)    
	{	
	    if ( BOOST_overlap(it_ev->second.LCRs[0],it_ev->second.LCRs[1]) )
		it_ev->second.print_this_entry();		
	}	
    }
    std::cerr << "\n\n\t\tDONE checking overlapping Event LCRs!!!!!\n";
    
    std::cerr << "\n\ninside \"check_informativeness_of_every_potential_breakpoint\"\n\n";
    typedef  std::map<uint, std::string>  type_map_uint_to_string;
    type_map_uint_to_string   map_fraglen_to_GC_rates_name;
		{
		    std::string GC_rates_by_frag_size_fname(BAM_dir);
		    GC_rates_by_frag_size_fname.append("/GC_rates__Readgroups/GC_rates_by_frag_size_library");    
		    std::cerr << "\n\treading from [" << GC_rates_by_frag_size_fname << "]...\n";
		    
		    map_fraglen_to_GC_rates_name = read_map_from_file<uint, std::string>(GC_rates_by_frag_size_fname, false);
		    
		    for (type_map_uint_to_string::iterator it_frag = map_fraglen_to_GC_rates_name.begin();
			    it_frag != map_fraglen_to_GC_rates_name.end();
			    ++it_frag)
		    {
			char ratesname[option_size];
			std::sprintf(ratesname, "%s/GC_rates__Readgroups/%s",
				     BAM_dir.c_str(),
				     it_frag->second.c_str()  );
				     
			it_frag->second.assign(ratesname);
		    }
		}
		
    
    std::cerr << "\n\tvectorizing loopable_Events...\n"; 
    std::vector<type_map_uint_to_Event::const_iterator>  loopable_Events;
			    {//init
				uint ctr=0;
				loopable_Events.reserve(Conn_Comps.size()*5);
				for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
					it_cc != Conn_Comps.end();
					++it_cc)
				{
				    std::cerr << "\n\t\tctr = " << ctr << ",  cc = " << it_cc->first << ", event.size() = " << it_cc->second.events.size();
				    for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
					    it_ev != it_cc->second.events.end();
					    ++it_ev)
				    {  loopable_Events.push_back(it_ev);  }
				    ++ctr;
				}	
			    }//init
			    
			    
			    
    std::cerr << "\n\n\nuploaded_pieces_of_chromos__all_Events.....!\n";
    type_map_uint_to_list_BI__string  uploaded_pieces_of_chromos__all_Events;
    {//upload    
	std::cerr << "\n\tgathering relevant genome...";
	type_map_uint_to_list_BI  map_chromo_to_homol_regs;	
	for (uint ev_ctr = 0; ev_ctr < loopable_Events.size(); ++ev_ctr)
	{
	    for (uint qs=0; qs<2; ++qs)
	    {
		map_chromo_to_homol_regs[loopable_Events[ev_ctr]->second.chromos[qs]]
				.push_back( BOOST_Interval(
						    safe_subtract_base_1(loopable_Events[ev_ctr]->second.LCRs[qs].lower(), 15000),
						    safe_chromo_add_base_1(loopable_Events[ev_ctr]->second.LCRs[qs].upper(), 15000, loopable_Events[ev_ctr]->second.chromos[qs]))  );		
	    }//qs   
	}//ev	
	std::cerr << "merging...";
	merge_intervals_on_each_chromosome(map_chromo_to_homol_regs);                             		    
	
	std::cerr << "uploading...";
	for (type_map_uint_to_list_BI::const_iterator it_chr = map_chromo_to_homol_regs.begin();
		it_chr != map_chromo_to_homol_regs.end();
		++it_chr)
	{
	    for (type_list_BI::const_iterator it_BI = it_chr->second.begin();
		    it_BI != it_chr->second.end();
		    ++it_BI)
	    {
		const std::string desired_sequence(
					    get_specific_DNA_sequence_from_Reference_genome(it_chr->first, *it_BI)  );                                                                                                 
												
		uploaded_pieces_of_chromos__all_Events[ it_chr->first ]
					    .push_back( type_BI__string(*it_BI, desired_sequence) );     
	    }//bi	
	}//chr
	std::cerr << "DONE.\n\n";
    }//upload
			    
				    
				    
				    
				
				
				
				
				
				    
    //construct miniRegions:
    std::cerr << "\n\t\tconstruct miniRegions...";
    
    std::vector< std::pair<MiniRegion,MiniRegion> > loopable_MR_pairs(loopable_Events.size(), std::pair<MiniRegion,MiniRegion>() );

    #pragma omp parallel for schedule(dynamic,5)
    for (uint ev_ctr = 0; ev_ctr < loopable_Events.size(); ++ev_ctr)
    {
	MiniRegion mini_LCR[2];
	for (uint qs=0; qs<2; ++qs)
	{
	    mini_LCR[qs].chromosome_of_region = loopable_Events[ev_ctr]->second.chromos[qs];
	    mini_LCR[qs].region_interval = BOOST_Interval(safe_subtract_base_1(loopable_Events[ev_ctr]->second.LCRs[qs].lower(), 10000),
							    safe_chromo_add_base_1(loopable_Events[ev_ctr]->second.LCRs[qs].upper(), 10000, loopable_Events[ev_ctr]->second.chromos[qs]));
	    mini_LCR[qs].orientation_of_MiniRegion_and_spawning_Event_profile_agree
				    = (qs == 0 or  loopable_Events[ev_ctr]->second.recomb_type == recomb_class__DupDel );
	    mini_LCR[qs].MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile 
				    = (qs == 0 or  loopable_Events[ev_ctr]->second.recomb_type == recomb_class__DupDel );
	    
	    mini_LCR[qs].set_region_sequence_from_preloaded_sequence(uploaded_pieces_of_chromos__all_Events);
	}
	
	loopable_MR_pairs.at(ev_ctr) = std::pair<MiniRegion,MiniRegion>(mini_LCR[0],mini_LCR[1]);
    }//ev_ctr    
    
    
    
    
    

				    
				    
				    
    
    for (type_list__uint__uint::const_iterator it_info_param = informativeness_checking_parameters.begin();
	    it_info_param != informativeness_checking_parameters.end();
	    ++it_info_param)
    {
	std::cerr << "\n\tit_info_param = [" << it_info_param->first << ",  " << it_info_param->second << "]\n\n";
	Readgroup_stats.clear();
	
	//set up dummy Readgroup
	std::string rgname;
	{
	    char rgchar[option_size];
	    std::sprintf(rgchar, "dummy_RG__frag_%u__mate_%u", it_info_param->first, it_info_param->second);
	    rgname.assign(rgchar);
	}
	
	Readgroup_stats[rgname] =  Readgroup_statistics(rgname, 1, 1, it_info_param->first, 50, it_info_param->first, 50);
	Readgroup_stats[rgname].average_mate_length =  it_info_param->second;
	
	std::cerr << "\n\t\tfinding dummy gc rate\n\n";
	{//dummy GC rate
	    type_map_uint_to_string  map_diff_fraglen_to_GC_rates_name;
	    for (type_map_uint_to_string::const_iterator  it_frag = map_fraglen_to_GC_rates_name.begin();
		    it_frag != map_fraglen_to_GC_rates_name.end();
		    ++it_frag)
	    {  map_diff_fraglen_to_GC_rates_name[(uint)std::abs((double)it_frag->first - (double)it_info_param->first)] = it_frag->second;  }
	    
	
// 		type_map_uint_to_string::const_iterator  it_closest = map_fraglen_to_GC_rates_name.end();
// 		uint leniancy = 5;	  
// 		uint closest_diff = 1000;
// 		while (it_closest == map_fraglen_to_GC_rates_name.end()  and   leniancy < 10000)
// 		{	    
// 		    for (type_map_uint_to_string::const_iterator  it_frag = map_fraglen_to_GC_rates_name.begin();
// 			    it_frag != map_fraglen_to_GC_rates_name.end();
// 			    ++it_frag)
// 		    {  
// 			
// 			if (std::abs((double)it_frag->first - (double)it_info_param->first) <= leniancy  )
// 			{
// 			    it_closest = it_frag;
// 			    closest_diff = std::abs((double)it_frag->first - (double)it_info_param->first);
// 			    break;			
// 			}
// 		    }
// 		    leniancy += 5;
// 		}
	    type_map_uint_to_string::const_iterator  it_closest = map_diff_fraglen_to_GC_rates_name.end();	    
	    for (type_map_uint_to_string::const_iterator it_potential = map_diff_fraglen_to_GC_rates_name.begin();
		    it_potential != map_diff_fraglen_to_GC_rates_name.end();
		    ++it_potential)
	    {
		std::cerr << "\n\t\treading map from [" << it_potential->second << "]...";
		Readgroup_stats[rgname].haploid_fragmentation_rate__per_GC_count 
				= read_map_from_file<uint, longdouble>(it_potential->second, true);
				
		if (std::abs((double)Readgroup_stats[rgname].haploid_fragmentation_rate__per_GC_count.size() - ((double)it_info_param->first+1))  ==  it_potential->first  )
		{//check that is not degenerate
		    it_closest = it_potential;
		    break;		
		}	    
	    }
	    
	    type_map_uint_to_string::const_iterator it_closest_actual_val = map_fraglen_to_GC_rates_name.end();
	    for (type_map_uint_to_string::const_iterator  it_frag = map_fraglen_to_GC_rates_name.begin();
		    it_frag != map_fraglen_to_GC_rates_name.end();
		    ++it_frag)
	    {  
		if (it_frag->second.compare(it_closest->second) == 0)
		{
		    it_closest_actual_val = it_frag;
		    break;
		}
	    }


	    Readgroup_stats[rgname].median_insert_size = it_closest_actual_val->first;
	    Readgroup_stats[rgname].mean_insert_size = it_closest_actual_val->first;
	    
	    
	    std::cerr << "\n\n\n\t\tdesired frag length = " << it_info_param->first
		    << " , closest frag length difference = " << it_closest->first
		    << ",   actual closest frag length = " << (uint)std::abs((double)it_closest->first - (double)it_info_param->first) 
		    << "\n\t\tdummy RG = [" << it_closest->second << "]\n\t\t"
		    << "Readgroup_stats[rgname].haploid_fragmentation_rate__per_GC_count.size() = " << Readgroup_stats[rgname].haploid_fragmentation_rate__per_GC_count.size() << "\n\n";		    	   
	    std::cerr.flush();
	    
	    std::cerr << "\n\n\n\t\thalve...\n\n";
	    for (type_map_uint_to_longdouble::iterator it_rate = Readgroup_stats.at(rgname).haploid_fragmentation_rate__per_GC_count.begin();
		    it_rate != Readgroup_stats.at(rgname).haploid_fragmentation_rate__per_GC_count.end();
		    ++it_rate)
	    {  it_rate->second /= 2;  }
	}//dummy GC rate
	
	
	print_map_keys_and_values<uint, longdouble>(Readgroup_stats.at(rgname).haploid_fragmentation_rate__per_GC_count, "haploid_fragmentation_rate__per_GC_count");
	
	
	//normalize by number of reads.    This requires the reference
	
	//GC-rate from file looks like:
	//
	//	lambda_GC =   Number-of-reads-mapping-to-region-with-GC-content-equals-g  /   Number-of-positons-in-Reference-with-GC-content-of-average-frag-length-equal-to-g
	//		  :=   R(g) / F(g)
	//
	//  If we divide this by the entire number of reads in the sample, then we obtain a normalized rate.  i.e. " the percent of reads in the sample which map to this GC count".
	//  i.e., if T_R = total number of reads in sample
	//
	//	p(g)  =  R(g) / ( F(g) * T_R ) 
	//
	//  Thereafter, if somebody wants to sequence with coverage "T" (assume same sequencing bias), then the expected number of reads at a Refernece location with GC-content = g
	//  	would be given by:    p(g) * T
	//
	//  Now, for simplicity we want to express values in terms of coverage.   lambda = coverage :=  (T / L_F) * L_R
	//	where  T = number of reads,  L_F = length of Reference genome,  and L_R = length of read, i.e. fragment length.
	//  We want to do "something-times-coverage = expected number of reads displaying breakpoint".  Note that T = lambda * L_F / L_R
	//  Thus, we should do:
	//	
	//		rho(g)   :=   p(g) * L_F / L_R
	//
	//   Thus, when we multiply   " rho(g)  * lambda "   ,  we will get the expected number of reads at the position, given the fragment length and bias.
	//
	//   In summary, we want to record:
	//			(  R(g) / F(g)  )  *   (1 / T_R)  *   (L_F / L_R)
	{//adjust
	    char exact_ref_count_name[option_size];
	    std::sprintf(exact_ref_count_name, "%s/GC_rates__Reference/exact_ref_count_per_GC_composition_for_fragment_length_%u",
						BAM_dir.c_str(),
						Readgroup_stats[rgname].median_insert_size);
						
	    std::cerr << "\n\n\t\treading exact from   exact_ref_count_name = [" << exact_ref_count_name << "]\n";
						
	    const type_map_uint_to_uint  exact_Ref_counts(  read_map_from_file<uint,uint>(exact_ref_count_name, false)  );
	    
	    //get total number of reads:
	    longdouble total_number_of_reads_in_sample = 0;	    
	    for (type_map_uint_to_longdouble::const_iterator it_rate = Readgroup_stats.at(rgname).haploid_fragmentation_rate__per_GC_count.begin();
		    it_rate != Readgroup_stats.at(rgname).haploid_fragmentation_rate__per_GC_count.end();
		    ++it_rate)
	    {  total_number_of_reads_in_sample +=   it_rate->second * exact_Ref_counts.at(it_rate->first);  }
	    
	    std::cerr << "\n\n\t\tadjusting...\n";
	    //Now we can adjust:
	    for (type_map_uint_to_longdouble::iterator it_rate = Readgroup_stats.at(rgname).haploid_fragmentation_rate__per_GC_count.begin();
		    it_rate != Readgroup_stats.at(rgname).haploid_fragmentation_rate__per_GC_count.end();
		    ++it_rate)
	    {
		it_rate->second = (it_rate->second  /  total_number_of_reads_in_sample)  
				    *  length_of_haploid_genome  /  Readgroup_stats.at(rgname).median_insert_size;
	    }	
	    
	    print_map_keys_and_values<uint, longdouble>(Readgroup_stats.at(rgname).haploid_fragmentation_rate__per_GC_count,
							"adjusted haploid_fragmentation_rate__per_GC_count",
							NULL, true);
	}//adjust
	
	
	
	
	
	
	
	//save to here:
	type_map_uint_to_uint_to_longdouble   expected_haploid_frag_rate_per_Event_per_breakpoint;
	
	for (uint ev_ctr = 0; ev_ctr < loopable_Events.size(); ++ev_ctr)
	{  expected_haploid_frag_rate_per_Event_per_breakpoint[loopable_Events.at(ev_ctr)->first] = type_map_uint_to_longdouble();  }
	
	
	
	std::cerr << "\n\t\tcalculating expected_haploid_frag_rate_per_Event...";	
	std::cerr.flush();
	
	//MPI
	for (int ev_ctr = my_MPI_rank; ev_ctr < loopable_Events.size(); ev_ctr += size_MPI_WORLD)
	{	    
// 	    if (ev_ctr % 500 == 0)
// 	    {  std::cerr << "\n\t\tev_ctr = " << ev_ctr;  }
	    
	    const std::vector<type_set_uint::const_iterator> loopable_varpos 
				= get_loopable_constiterators<type_set_uint>(loopable_Events[ev_ctr]->second.variational_positions_of_profile);
	    
	    //allocate
	    for (uint vp_ctr = 0; vp_ctr < loopable_varpos.size(); ++vp_ctr)
	    {  expected_haploid_frag_rate_per_Event_per_breakpoint.at(loopable_Events[ev_ctr]->second.UID)[*loopable_varpos.at(vp_ctr)] = 0;  }
		
	    //compute
	    #pragma omp parallel for schedule(dynamic,10)
	    for (uint vp_ctr = 0; vp_ctr < loopable_varpos.size(); ++vp_ctr)
	    {			
		expected_haploid_frag_rate_per_Event_per_breakpoint.at(loopable_Events[ev_ctr]->second.UID).at(*loopable_varpos.at(vp_ctr))
				= loopable_Events[ev_ctr]->second.calculate_expected_number_of_PERs_at_breakpoint(
									    BOOST_make_point(*loopable_varpos.at(vp_ctr)),
									    loopable_MR_pairs.at(ev_ctr).first,
									    loopable_MR_pairs.at(ev_ctr).second,
									    false).first;
	    }//brkpt	    
	}//ev_ctr
	
// 	std::cerr << "\n\nrank " << my_MPI_rank << "  waiting at gather...";
	
	std::vector<type_map_uint_to_uint_to_longdouble> expected_haploid_frag_rate_per_event___by_rank;
	boost::mpi::gather<type_map_uint_to_uint_to_longdouble>( 
					    *world_MPI_boost_ptr,
					    expected_haploid_frag_rate_per_Event_per_breakpoint,
					    expected_haploid_frag_rate_per_event___by_rank, 0);
					    
// 	std::cerr << "\n\nrank " << my_MPI_rank << " past gather.";		
					    
	if (my_MPI_rank == 0)
	{
	    //combine info across ranks
	    for (int rank=1; rank < size_MPI_WORLD; ++rank)
	    {
		for (int ev_rank_ctr = rank; ev_rank_ctr < loopable_Events.size();  ev_rank_ctr += size_MPI_WORLD)  //just like above
		{
		    expected_haploid_frag_rate_per_Event_per_breakpoint.at( loopable_Events[ev_rank_ctr]->second.UID ) 
			    = expected_haploid_frag_rate_per_event___by_rank.at(rank).at( loopable_Events[ev_rank_ctr]->second.UID );
		}//ev_rank_ctr
	    }//rank
		
	    //save
	    std::cerr << "\n\t\twriting to file\n\n";
	    char outfname[option_size];
	    std::sprintf(outfname, "%s/expected_haploid_frag_rate_per_Event__frag_%u__mate_%u",
				    some_outputdir.c_str(),
				    it_info_param->first,
				    it_info_param->second  );				
				
	    std::ofstream  outfs(outfname);			
			if ( !outfs.good()  )
			{
			    std::stringstream error_strm;
			    print_line_of_markers("ERROR! ", &error_strm);
			    print_line_of_markers("(", &error_strm);
			    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
			    
			    error_strm << "\n\nERROR!   not good ofstream   in  \"check_informativeness_of_every_potential_breakpoint\",\n\t\t outfname =    "  
			    <<  outfname  << "\n\n";
			    
			    print_line_of_markers(")", &error_strm);
			    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );  
			    return;
			}                             	
	
	    for (type_map_uint_to_uint_to_longdouble::const_iterator it_ev = expected_haploid_frag_rate_per_Event_per_breakpoint.begin();
		    it_ev != expected_haploid_frag_rate_per_Event_per_breakpoint.end();
		    ++it_ev)
	    {
		for (type_map_uint_to_longdouble::const_iterator it_brkpt = it_ev->second.begin();
			it_brkpt != it_ev->second.end();
			++it_brkpt)
		{
		    outfs << it_ev->first << "\t" << it_brkpt->first << "\t" << it_brkpt->second << "\n";
		}	
	    }//ev		
	    
	    outfs.close();	
	}//0
	
    }//info	    
    
}//check_informativeness_of_every_potential_NAHR_breakpoint
    
    
    
    
    
    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

// void validate_using_probe_intensity_comparison
// 			(const std::string &probe_intensity_filename,
// 			 const uint &chr_of_call,
// 			 const BOOST_Interval &brkpt_of_call,
// 			 const uint &some_aCGH_data.padding_amount,
// 			 const type_set_string &individuals_with_call,
// 			 const std::string &some_aCGH_data.acceptable_genomes_filename,
// 			 const std::string &some_aCGH_data.save_prefix)
void validate_using_probe_intensity_comparison
	    (const aCGH_validation_container &some_aCGH_data)
{
    
    std::cerr << "\n\ninside \"validate_using_probe_intensity_comparison\"...\n\n";
    // Three types of informaton is extracted and saved to file for further analysis:
    //		1: per-position sum of probe intensities (and number of contributing probes) for each individual.  This is for making "read-depth" plots.
    //		2: a list of probe intensities falling withing the rearrangement, for each individual.  This is for rank-sum tests.
    //		3: a list of probe intensities falling outside the rearrangement, for each individual.  For comparing inside earrangemnet vs. outside rearrangments rank-sum test.
    // For each type of information above, individuals are split into two groups (and files are saved as such): individuals with call vs. those without.
    
    const BOOST_Interval  padded_brkpt_region( safe_subtract_base_1(some_aCGH_data.brkpt_of_call.lower(), some_aCGH_data.padding_amount),
					       safe_chromo_add_base_1(some_aCGH_data.brkpt_of_call.upper(), some_aCGH_data.padding_amount, some_aCGH_data.chr_of_call)   );
					       
    std::cerr << "\n\t\tpadded_brkpt_region = " << padded_brkpt_region << "\n";
					       
    type_map_string_to_uint_to__uint__uint  map_individual_to_position_to__probe_sum__probe_count; 
    type_map_string_to_list_uint  map_individual_to_probes_within_region_of_call;
    type_map_string_to_list_uint  map_individual_to_probes_in_padded_region_outside_of_call;
    type_map_uint_to_string  map_column_to_individualname;
					       				         

    FILE *in_probe_intensities_file = std::fopen(some_aCGH_data.probe_intensity_filename.c_str(), "r");
                                    if (in_probe_intensities_file == NULL)
                                    {
                                        std::fprintf(stderr, "ERROR: unable to open file [%s] in \"validate_using_probe_intensity_comparison\".\n\n",
                                                    some_aCGH_data.probe_intensity_filename.c_str() );
                                        exit(1);
                                    }          
                                    
    //get header and use it to intiailize names and column<->individual correspondence
    {//header
	char header[option_size*10];    
	std::fscanf(in_probe_intensities_file, "%[^\n]\n", header);    
	std::stringstream editstrm;
	editstrm.str(header);
            
// 	type_map_uint_to__uint__uint dummy_empty;
// 					{//init
// 					    type_map_uint_to__uint__uint::iterator insert_it = dummy_empty.begin();
// 					    for (uint j=padded_brkpt_region.lower(); j <= padded_brkpt_region.upper(); ++j)
// 					    {  insert_it = dummy_empty.insert( insert_it, std::pair<uint, type_uint__uint>(j, type_uint__uint(0,0)) );  }
// 					}//init
	        
	{
	    std::string individual_name;
	    uint ctr = 0;
	    while (editstrm >> individual_name)
	    {
		map_individual_to_position_to__probe_sum__probe_count[individual_name] = type_map_uint_to__uint__uint();//dummy_empty;			   
		map_individual_to_probes_within_region_of_call[individual_name] = type_list_uint();
		map_individual_to_probes_in_padded_region_outside_of_call[individual_name] = type_list_uint();
		map_column_to_individualname[ctr] = individual_name; 
		++ctr;
	    }
	}
    }//header
    
    
//     print_map_keys_and_values<uint,std::string>(map_column_to_individualname, "map_column_to_individualname", NULL, false);
    
    
    char chr_of_call_char;
    if (some_aCGH_data.chr_of_call  == 23)
    {  chr_of_call_char = 'X';  }
    else if  (some_aCGH_data.chr_of_call  == 24)
    {  chr_of_call_char = 'Y';  }	
    else
    {  chr_of_call_char = (char)some_aCGH_data.chr_of_call;  }
    
    
    std::cerr << "\n\t\tscanning probe intensity file...";
    
    uint number_good_lines = 0;
    char in_probe_name[option_size];
    char in_chr_val[option_size];
    char in_begin[option_size];
    char in_end[option_size];   
    while (std::fscanf(in_probe_intensities_file, "%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t",
						    in_probe_name,
						    in_chr_val,
						    in_begin,
						    in_end)   !=   EOF)
    {	
	const BOOST_Interval probe_region((uint)atoi(in_begin), (uint)atoi(in_end));
	if (  (char)atoi(in_chr_val) == chr_of_call_char
	      and  BOOST_subset(probe_region, padded_brkpt_region)   )
	{	    
	    ++number_good_lines;	    
	    char in_rest_of_line[option_size*10];
	    std::fscanf(in_probe_intensities_file, "%[^\n]\n", in_rest_of_line);
	    
	    	    
	    type_vector_string  in_intensities;
	    const std::string in_line_str(in_rest_of_line);
	    boost::split(in_intensities, in_line_str, boost::is_any_of("\t"), boost::algorithm::token_compress_on);
	    
	    for (uint i_ctr=0; i_ctr < in_intensities.size(); ++i_ctr)
	    {
		if (BOOST_overlap(probe_region, some_aCGH_data.brkpt_of_call))
		{  
		    map_individual_to_probes_within_region_of_call.at(map_column_to_individualname.at(i_ctr))
								    .push_back((uint)atoi(in_intensities.at(i_ctr).c_str())); 
		}	
		else if (probe_region.upper() < some_aCGH_data.brkpt_of_call.lower()  or some_aCGH_data.brkpt_of_call.upper() < probe_region.lower())
		{ 
		    map_individual_to_probes_in_padded_region_outside_of_call.at(map_column_to_individualname.at(i_ctr))
								    .push_back((uint)atoi(in_intensities.at(i_ctr).c_str()));
		}
		
		for (uint pos=(uint)atoi(in_begin); pos <= (uint)atoi(in_end); ++pos)		
		{
		    map_individual_to_position_to__probe_sum__probe_count.at(map_column_to_individualname.at(i_ctr))[pos].first
								    += (uint)atoi(in_intensities.at(i_ctr).c_str());
		    ++map_individual_to_position_to__probe_sum__probe_count.at(map_column_to_individualname.at(i_ctr))[pos].second;		    
		}					    
	    }//intensities			    
	    
// 	    std::stringstream instrm;
// 	    instrm.str(in_rest_of_line);
// 	    
// 	    uint in_probe_intensity;
// 	    uint ctr=0;
// 	    while (instrm >> in_probe_intensity)
// 	    {
// 		if (BOOST_subset(probe_region, some_aCGH_data.brkpt_of_call))
// 		{  
// 		    map_individual_to_probes_within_region_of_call.at(map_column_to_individualname.at(ctr)).push_back(in_probe_intensity); 
// 		}	
// 		else if (probe_region.upper() < some_aCGH_data.brkpt_of_call.lower()  or some_aCGH_data.brkpt_of_call.upper() < probe_region.lower())
// 		{ 
// 		    map_individual_to_probes_in_padded_region_outside_of_call.at(map_column_to_individualname.at(ctr)).push_back(in_probe_intensity);
// 		}
// 		
// 		for (uint j=(uint)atoi(in_begin); j <= (uint)atoi(in_end); ++j)		
// 		{
// 		    map_individual_to_position_to__probe_sum__probe_count.at(map_column_to_individualname.at(ctr))[j].first += in_probe_intensity;
// 		    ++map_individual_to_position_to__probe_sum__probe_count.at(map_column_to_individualname.at(ctr))[j].second;		    
// 		}
// 		
// 		++ctr;
// 	    }
	}//same chr
	else if (       ((char)atoi(in_chr_val) == chr_of_call_char  and  padded_brkpt_region.upper() < (uint)atoi(in_begin))  
		    or  (chr_of_call_char != 'X' and chr_of_call_char != 'Y'
			    and  strcasecmp(in_chr_val,"X") != 0   and     strcasecmp(in_chr_val,"Y") != 0
			    and  (int)chr_of_call_char < atoi(in_chr_val))     )
	{  break;  }
	else
	{
	    std::fscanf(in_probe_intensities_file, "%*[^\n]");
	    std::fscanf(in_probe_intensities_file, "\n");	
	}
    }//while
    
    
    std::fclose(in_probe_intensities_file);
    
    
    
    std::cerr << "DONE!!!\n\nnumber of good lines read in from probe intensity file = " << number_good_lines << "\n\n";
    
    
    
    
    
    
    //save to file:                    
    
    const type_set_string acceptable_individuals( read_genome_list_from_file(some_aCGH_data.acceptable_genomes_filename) );
    
    
    std::cerr << "\n\t\tsaving map_positions_to_uniqueness...\n";
    {//map_positions_to_uniqueness	
	type_map_uint_to_uint map_positions_to_uniqueness;
    
	type_map_uint_to_uint::iterator insert_it = map_positions_to_uniqueness.begin();
	for (uint j=padded_brkpt_region.lower(); j <= padded_brkpt_region.upper(); ++j)
	{  insert_it = map_positions_to_uniqueness.insert(insert_it, type_uint__uint(j,1));  }
	
	const type_map_uint_to_list_BI non_uniq_regions_of_genome(  load_Non_unique_regions_of_genome()  );
	if (non_uniq_regions_of_genome.count(some_aCGH_data.chr_of_call) > 0)
	{
	    for (type_list_BI::const_iterator it_bi = non_uniq_regions_of_genome.at(some_aCGH_data.chr_of_call).begin();
		    it_bi != non_uniq_regions_of_genome.at(some_aCGH_data.chr_of_call).end();
		    ++it_bi)
	    {
		const BOOST_Interval nonuniq_and_pos(BOOST_intersect(*it_bi, padded_brkpt_region));
		if (!BOOST_empty(nonuniq_and_pos))
		{
		    for (uint j=nonuniq_and_pos.lower(); j <= nonuniq_and_pos.upper(); ++j)
		    {  map_positions_to_uniqueness.at(j) = 0; }
		}	
	    }			    
	}
	
	char save_name[option_size];
	std::sprintf(save_name, "%s__uniqueness_by_position", some_aCGH_data.save_prefix.c_str());		
	write_map_to_file<uint,uint>(save_name, map_positions_to_uniqueness);
    }//map_positions_to_uniqueness




    std::cerr << "\n\t\tsaving map_individual_to_position_to__probe_sum__probe_count...\n";
    {//map_individual_to_position_to__probe_sum__probe_count
	char save_name[option_size];
	
	std::sprintf(save_name, "%s__person_to_probe_sum_per_position__positive", some_aCGH_data.save_prefix.c_str());
	std::ofstream outfs_sum_pos(save_name);
	check_fstream_for_goodness(outfs_sum_pos, save_name, "validate_using_probe_intensity_comparison", true);		
	
	std::sprintf(save_name, "%s__person_to_probe_sum_per_position__negative", some_aCGH_data.save_prefix.c_str());
	std::ofstream outfs_sum_neg(save_name);
	check_fstream_for_goodness(outfs_sum_neg, save_name, "validate_using_probe_intensity_comparison", true);
		
	std::sprintf(save_name, "%s__person_to_probe_sum_per_position__all", some_aCGH_data.save_prefix.c_str());
	std::ofstream outfs_sum_all(save_name);
	check_fstream_for_goodness(outfs_sum_all, save_name, "validate_using_probe_intensity_comparison", true);
	
	
	for (type_map_string_to_uint_to__uint__uint::const_iterator it_ind = map_individual_to_position_to__probe_sum__probe_count.begin();
		it_ind != map_individual_to_position_to__probe_sum__probe_count.end();
		++it_ind)
	{
// 	    //always write to all
// 	    {//all
// 		type_map_uint_to__uint__uint::const_iterator it_last = --it_ind->second.end();
// 		for (type_map_uint_to__uint__uint::const_iterator it_pos = it_ind->second.begin();
// 			it_pos != it_last;
// 			++it_pos)
// 		{
// 		    outfs_sum_all << it_pos->second.first << "\t";//tab (not linebreak).
// 		}	    
// 		outfs_sum_all << it_last->second.first << "\n";//linebreak.
// 		
// 		for (type_map_uint_to__uint__uint::const_iterator it_pos = it_ind->second.begin();
// 			it_pos != it_last;
// 			++it_pos)	  
// 		{
// 		    outfs_sum_all<< it_pos->second.second << "\t";//linebreak.
// 		}		    
// 		outfs_sum_all  << it_last->second.second << "\n";//linebreak.	    	    
// 	    }//all
	    
	    
	    
	    //write to subsets, if appropriate:	    	    
	    std::ofstream *ptr_correct_strm = NULL;
	    if (some_aCGH_data.individuals_with_call.count(it_ind->first) > 0)
		ptr_correct_strm = &outfs_sum_pos;
	    else if (acceptable_individuals.count(it_ind->first) > 0)
		ptr_correct_strm = &outfs_sum_neg;
			    
	    if (ptr_correct_strm != NULL)
	    {	    
		type_map_uint_to__uint__uint::const_iterator it_last = --it_ind->second.end();
		for (type_map_uint_to__uint__uint::const_iterator it_pos = it_ind->second.begin();
			it_pos != it_last;
			++it_pos)
		{
		    (*ptr_correct_strm) << it_pos->second.first << "\t";//tab (not linebreak).
		}	    
		(*ptr_correct_strm) << it_last->second.first << "\n";//linebreak.
		
		for (type_map_uint_to__uint__uint::const_iterator it_pos = it_ind->second.begin();
			it_pos != it_last;
			++it_pos)	  
		{
		    (*ptr_correct_strm) << it_pos->second.second << "\t";//linebreak.
		}		    
		(*ptr_correct_strm) << it_last->second.second << "\n";//linebreak.
	    }	    	    
	}//ind	
	
	outfs_sum_pos.close();
	outfs_sum_neg.close();
	outfs_sum_all.close();
	
    }//map_individual_to_position_to__probe_sum__probe_count
    
    
    
    
    
    std::cerr << "\n\t\tsaving map_individual_to_probes_within_region_of_call...\n";
    {//map_individual_to_probes_within_region_of_call
	char save_name[option_size];
	
	std::sprintf(save_name, "%s__person_to_probes__within_call__positive", some_aCGH_data.save_prefix.c_str());	
	std::ofstream outfs_with_call(save_name);
	check_fstream_for_goodness(outfs_with_call, save_name, "validate_using_probe_intensity_comparison", true);
	
	std::sprintf(save_name, "%s__person_to_probes__within_call__negative", some_aCGH_data.save_prefix.c_str());
	std::ofstream outfs_without_call(save_name);
	check_fstream_for_goodness(outfs_without_call, save_name, "validate_using_probe_intensity_comparison", true);
	
	std::sprintf(save_name, "%s__person_to_probes__within_call__all", some_aCGH_data.save_prefix.c_str());
	std::ofstream outfs_all(save_name);
	check_fstream_for_goodness(outfs_all, save_name, "validate_using_probe_intensity_comparison", true);	
		
				
	for (type_map_string_to_list_uint::const_iterator it_ind = map_individual_to_probes_within_region_of_call.begin();
		it_ind != map_individual_to_probes_within_region_of_call.end();
		++it_ind)
	{
// 	    //always write to all:
// 	    print_singlecontainer_to_stream<std::ofstream, type_list_uint>(outfs_all, it_ind->second, true);
	    
	    if (some_aCGH_data.individuals_with_call.count(it_ind->first) > 0)
	    {  
		print_singlecontainer_to_stream<std::ofstream, type_list_uint>(outfs_with_call, it_ind->second, true);
	    }
	    else if (acceptable_individuals.count(it_ind->first) > 0)
	    {  
		print_singlecontainer_to_stream<std::ofstream, type_list_uint>(outfs_without_call, it_ind->second, true);
	    }						
	}//it_ind
	
	outfs_with_call.close();
	outfs_without_call.close();
	outfs_all.close();
				
    }//map_individual_to_probes_within_region_of_call
    
    
    std::cerr << "\n\t\tsaving map_individual_to_probes_in_padded_region_outside_of_call...\n";
    {//map_individual_to_probes_in_padded_region_outside_of_call
	char save_name[option_size];
	
	std::sprintf(save_name, "%s__person_to_probes__padded_region_outside_call__positive", some_aCGH_data.save_prefix.c_str());	
	std::ofstream outfs_with_call(save_name);
	check_fstream_for_goodness(outfs_with_call, save_name, "validate_using_probe_intensity_comparison", true);
	
	std::sprintf(save_name, "%s__person_to_probes__padded_region_outside_call__negative", some_aCGH_data.save_prefix.c_str());
	std::ofstream outfs_without_call(save_name);
	check_fstream_for_goodness(outfs_without_call, save_name, "validate_using_probe_intensity_comparison", true);
	
	std::sprintf(save_name, "%s__person_to_probes__padded_region_outside_call__all", some_aCGH_data.save_prefix.c_str());
	std::ofstream outfs_all(save_name);
	check_fstream_for_goodness(outfs_all, save_name, "validate_using_probe_intensity_comparison", true);	
		
				
	for (type_map_string_to_list_uint::const_iterator it_ind = map_individual_to_probes_in_padded_region_outside_of_call.begin();
		it_ind != map_individual_to_probes_in_padded_region_outside_of_call.end();
		++it_ind)
	{
// 	    //always write to all:
// 	    print_singlecontainer_to_stream<std::ofstream, type_list_uint>(outfs_all, it_ind->second, true);
	    
	    if (some_aCGH_data.individuals_with_call.count(it_ind->first) > 0)
	    {  
		print_singlecontainer_to_stream<std::ofstream, type_list_uint>(outfs_with_call, it_ind->second, true);
	    }
	    else if (acceptable_individuals.count(it_ind->first) > 0)
	    {  
		print_singlecontainer_to_stream<std::ofstream, type_list_uint>(outfs_without_call, it_ind->second, true);
	    }	    
	}//it_ind
	
	
	outfs_with_call.close();
	outfs_without_call.close();
	outfs_all.close();
	
    }//map_individual_to_probes_in_padded_region_outside_of_call    
        
        
        
    std::cerr << "\n\t\tDONE  with \"validate_using_probe_intensity_comparison\"!!!!\n";
                     
}//validate_using_probe_intensity_comparison


























// void remove_points_from_set_of_BI
// 		    (type_set_BI &some_set)
// {
//     type_set_BI::iterator  it_bi = some_set.begin();
//     while (it_bi != some_set.end())
//     {
// 	if (BOOST_is_a_point(*it_bi))
// 	    some_set.erase(it_bi++);
// 	else
// 	    ++it_bi;    
//     }    
// 
// }//remove_points_from_set_of_BI


// void remove_nondegenerate_intervals_from_set_of_BI
// 		    (type_set_BI &some_set)
// {
//     type_set_BI::iterator  it_bi = some_set.begin();
//     while (it_bi != some_set.end())
//     {
// 	if (BOOST_is_a_point(*it_bi))
// 	    ++it_bi;
// 	else
// 	    some_set.erase(it_bi++);    
//     }    
// 
// }//remove_nondegenerate_intervals_from_set_of_BI






uint  create_haploid_breakpoint_space
		(type_map_uint_to_vector_BI &brkpt_space_per_Event,
		 const Sparse_map &haploid_outcomes_per_Event,
		type_map_uint_to_Breakpoint_complex &effective_Breakpoint_complex)
{
    //brkpt_space:  [event][brkpt]
    
    type_list_uint  space_product;
    space_product.push_back(1);  // last entry will be irrelevant.
    
    for (type_map_uint_to_uint::const_iterator it_ev = haploid_outcomes_per_Event.sparse_map_UID_to_value.begin();
	    it_ev != haploid_outcomes_per_Event.sparse_map_UID_to_value.end();
	    ++it_ev)
    {
	switch(static_cast<type_haploid_outcome>(it_ev->second))
	{
	    case hap_outcome__Del:
		space_product.push_back(effective_Breakpoint_complex.at(it_ev->first).breakpoints_NAHR__AB.size());
		break;
	    case hap_outcome__Dup:
		space_product.push_back(effective_Breakpoint_complex.at(it_ev->first).breakpoints_NAHR__BA.size());
		break;
	    case hap_outcome__Inv:
		space_product.push_back(effective_Breakpoint_complex.at(it_ev->first).breakpoints_NAHR__AB.size()); // either one
		break;
	    case hap_outcome__GeneConv_ABA:
		space_product.push_back(effective_Breakpoint_complex.at(it_ev->first).breakpoints_GeneConversion___ABA.size());
		break;
	    case hap_outcome__GeneConv_BAB:
		space_product.push_back(effective_Breakpoint_complex.at(it_ev->first).breakpoints_GeneConversion___BAB.size());
		break;		
	    default:
		std::stringstream error_strm;
		error_strm << "\n\n\nERROR in \"create_haploid_breakpoint_space\"!!!  unrecognized outcome = " 
			    << convert_haploid_outcome_to_string(static_cast<type_haploid_outcome>(it_ev->second)) << "\n";
		error_message(error_strm,true);
		break;
	}//switch
    }//ev
    
    const uint total_number_of_breakpoints = product_over_container<type_list_uint, uint>(space_product);
    
    

    
    
    //get cumulative product
    for (type_list_uint::iterator it_prev_prod = space_product.begin(),
				  it_cum_prod = ++space_product.begin();
	    it_cum_prod != space_product.end();
	    ++it_prev_prod, ++it_cum_prod)
    {
	*it_cum_prod *= *it_prev_prod;
    }    
    
    
    

    
    brkpt_space_per_Event.clear();
    type_map_uint_to_vector_BI::iterator insert_it_ev = brkpt_space_per_Event.begin();    
			
    type_list_uint::const_iterator it_cumprod_this_ev = space_product.begin();    
    for (type_map_uint_to_uint::const_iterator it_ev = haploid_outcomes_per_Event.sparse_map_UID_to_value.begin();
	    it_ev != haploid_outcomes_per_Event.sparse_map_UID_to_value.end();
	    ++it_ev)
    {	
	type_vector_BI brkpt_space(total_number_of_breakpoints, empty_BI);
		
	const type_set_BI* it_appropriate_breakpoints = NULL;
	
	switch(static_cast<type_haploid_outcome>(it_ev->second))
	{
	    case hap_outcome__Del:
		it_appropriate_breakpoints = &effective_Breakpoint_complex.at(it_ev->first).breakpoints_NAHR__AB;
		break;
	    case hap_outcome__Dup:
		it_appropriate_breakpoints = &effective_Breakpoint_complex.at(it_ev->first).breakpoints_NAHR__BA;
		break;
	    case hap_outcome__Inv:
		it_appropriate_breakpoints = &effective_Breakpoint_complex.at(it_ev->first).breakpoints_NAHR__AB; // either one
		break;
	    case hap_outcome__GeneConv_ABA:
		it_appropriate_breakpoints = &effective_Breakpoint_complex.at(it_ev->first).breakpoints_GeneConversion___ABA;
		break;
	    case hap_outcome__GeneConv_BAB:
		it_appropriate_breakpoints = &effective_Breakpoint_complex.at(it_ev->first).breakpoints_GeneConversion___BAB;
		break;		
	    default:
		std::stringstream error_strm;
		error_strm << "\n\n\nERROR in \"create_haploid_breakpoint_space\"!!!  for it_appropriate_breakpoints.  unrecognized outcome = " 
			    << convert_haploid_outcome_to_string(static_cast<type_haploid_outcome>(it_ev->second)) << "\n";
		error_message(error_strm,true);
		break;
	}//switch
	
	

	uint total_brkpt_ctr = 0;
	for (uint q = 0; q < total_number_of_breakpoints / (it_appropriate_breakpoints->size()*(*it_cumprod_this_ev)); ++q)
	{
	    for (type_set_BI::const_iterator it_bi = it_appropriate_breakpoints->begin();
		    it_bi != it_appropriate_breakpoints->end();
		    ++it_bi)
	    {
		for (uint p=0; p < *it_cumprod_this_ev; ++p)
		{
		    brkpt_space.at(total_brkpt_ctr) = *it_bi;
		    ++total_brkpt_ctr;
		}//p
	    }//bi
	}//q
	
	insert_it_ev = brkpt_space_per_Event.insert(insert_it_ev, std::pair<uint,type_vector_BI>(it_ev->first, brkpt_space));
	++it_cumprod_this_ev;
    }//ev        
    
    
    return total_number_of_breakpoints;

}//create_haploid_breakpoint_space










void  set_breakpoints_using_breakpoint_space
	(type_map_uint_to_BI &haploid_breakpoint_instance,
	 const type_map_uint_to_vector_BI &haploid_breakpoint_space,
	 const uint &brkpt_ctr)
{
    for (type_map_uint_to_BI::iterator it_ev = haploid_breakpoint_instance.begin();
	    it_ev != haploid_breakpoint_instance.end();
	    ++it_ev)
    {
	it_ev->second = haploid_breakpoint_space.at(it_ev->first).at(brkpt_ctr);        
    }//ev

}//set_breakpoints_using_breakpoint_space













uint  calculate_GC_content_of_sequence
	    (const std::string &some_sequence)	     
{
    uint gc_count = 0;
    for (uint pos=0; pos < some_sequence.size(); ++pos)
    {
		if (some_sequence[pos] == 'G'  or  some_sequence[pos] == 'C')
		{  ++gc_count;  }
    }    
    
    return gc_count;
    
}//calculate_GC_content_of_sequence













type_map_string_to_Readgroup_stats::const_iterator  identify_Readgroup_from_name_of_PER
			(const std::string &name_of_PER)
{
	if (is_simulation)
	{
		return  Readgroup_stats.begin();
	}



    const std::string rg_name(name_of_PER.substr(0, name_of_PER.find_last_of('.')));            
    type_map_string_to_Readgroup_stats::const_iterator it_find_rg = Readgroup_stats.find(rg_name);        
    if (it_find_rg == Readgroup_stats.end())
    {
	it_find_rg = Readgroup_stats.find(name_of_PER.substr(0, name_of_PER.find_first_of('.')));
	if (it_find_rg == Readgroup_stats.end())
	{
	    for (type_map_string_to_Readgroup_stats::const_iterator it_rg_all = Readgroup_stats.begin();
		    it_rg_all != Readgroup_stats.end();
		    ++it_rg_all)
	    {
		if (name_of_PER.find(it_rg_all->first) != std::string::npos)
		{
		    it_find_rg = it_rg_all;
		    break;		    
		}		
	    }	    
	}	    
    }
    
    return it_find_rg;    

}//identify_Readgroup_from_name_of_PER































// void determine_graphability_of_every_Event
// 		    (const type_map_uint_to_CC &the_Conn_Comps)
// {
//     type_map_uint_to_uint  map_EV_ID_to_graphability;
//     
// 
//     std::cerr << "\t\t\tidentifying nonuniqueness in genome...\n\n";
//     //non-uniqueness
//     const type_map_uint_to_list_BI non_uniq_regions_of_genome(load_Non_unique_regions_of_genome());	
// 	
//     
//     
//     
//     for (type_map_uint_to_CC::const_iterator it_cc = the_Conn_Comps.begin();
// 	    it_cc != the_Conn_Comps.end();
// 	    ++it_cc)
//     {
// 	for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.begin();
// 		it_ev != it_cc->second.end();
// 		++it_ev)
// 	{
// 	    type_set_uint affected_reg_as_set;
// 	    for (uint j = it_ev->second.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs.lower();
// 		     j <= it_ev->second.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs.upper(); ++j )
// 	    
// 	    const type_map_uint_to_list_BI::const_iterator it_chr_non = non_uniq_regions_of_genome.find(it_ev->second.chromos[0]);	
// 	    if (it_chr_non != non_uniq_regions_of_genome.end() )
// 	    {
// 		for (type_list_BI::const_iterator it_non = it_chr_non->second.begin();
// 			it_non != it_chr_non->second.end();
// 			++it_non)
// 		{
// 		    if ()
// 		    
// 		    for (uint pos = it_non->lower(); pos <= it_non->upper(); ++pos)
// 		    {
// 			type_map_uint_to_longdouble::iterator it_gc = positions_along_chromosome___NO_EVENT.find(pos);
// 			if (it_gc != positions_along_chromosome___NO_EVENT.end())
// 			{
// 			    it_gc->second = -99999.99L; // dummy value
// 			}
// 		    }
// 		}		
// 	    }//found chr		    
// 		
// 	
// 	}//it_ev      
//     }//it_cc       
//     
// }//determine_graphability_of_every_Event
// 








real  calculate_standard_normal_tail
		(const longdouble &standardized_Z_score)
{
    boost::math::normal_distribution<mpfr_class>  my_standard_normal_dist(0,1);
    
    return real(boost::math::cdf<mpfr_class>(my_standard_normal_dist, mpfr_class((double)standardized_Z_score)).get_mpfr_t());        
}



real  calculate_normal_pdf_or_cdf
		(const longdouble &observation,
		const longdouble &mean,
		 const longdouble &sigma,
		 const bool &pdf__false__cdf__true)
{
    boost::math::normal_distribution<mpfr_class>  my_normal_dist((double)mean, (double)sigma);
        
    if (!pdf__false__cdf__true)
	return real(boost::math::pdf<mpfr_class>(my_normal_dist, mpfr_class((double)observation)).get_mpfr_t());
    else
	return real(boost::math::cdf<mpfr_class>(my_normal_dist, mpfr_class((double)observation)).get_mpfr_t());
}//calculate_pdf_or_cdf_for_poisson_distribution




real  calculate_poisson_pdf_or_cdf
		(const longdouble &observation,
		const longdouble &poisson_parameter,
		 const longdouble &pseudo_amount,
		 const bool &pdf__false__cdf__true)
{
    boost::math::poisson_distribution<mpfr_class>  my_poisson_dist((double)(poisson_parameter + pseudo_amount));
        
    if (!pdf__false__cdf__true)
	return real(boost::math::pdf<mpfr_class>(my_poisson_dist, mpfr_class((double)observation)).get_mpfr_t());
    else
	return real(boost::math::cdf<mpfr_class>(my_poisson_dist, mpfr_class((double)observation)).get_mpfr_t());
}//calculate_pdf_or_cdf_for_poisson_distribution



real  calculate_negative_binomial_pdf_or_cdf__transforming_from_poisson_parameter
		(const longdouble &observation,
		const longdouble &poisson_parameter,
		 const longdouble &pseudo_amount,
		 const bool &pdf__false__cdf__true) 
{
    
    
    {//normal approx
	const mpfr_class mu((double)poisson_parameter);
	const mpfr_class sig((double)std::max<longdouble>(3,(sqrt(poisson_parameter*(1.00 + sqrt(RD_variance_scale___alpha_squared))))));
	const boost::math::normal_distribution<mpfr_class>  my_normal(mu,sig);    
		
	if (!pdf__false__cdf__true)    
	    return   real(boost::math::pdf<mpfr_class>(my_normal, (double)observation).get_mpfr_t());
	else
	    return   real(boost::math::cdf<mpfr_class>(my_normal, (double)observation).get_mpfr_t());	
    }//normal approx    
    
    
    //negative binomial:   copied from   Event::calculate_probability_of_observed_read_depth_given_expected_read_depth_over_natural_poisson_intervals               
    //NEGATIVE BINOMIAL RETURNS NEGATIVE VALUES!!!!
    // e.g.  
// // //         Observations:
// // //                 observed_counts__inferred = 19558.54028
// // //                 observed_counts__UNIQUE_ONLY = 15627.99999
// // // 
// // //         Parameters:
// // //                 poisson_parameter__null_hypothesis__inferred_counts = 20183.90696
// // //                 poisson_parameter__alternative_hypothesis__inferred_counts = 20171.4506
// // //                 poisson_parameter__null_hypothesis__UNIQUE_ONLY = 16210.80215
// // //                 poisson_parameter__alternative_hypothesis__UNIQUE_ONLY = 16210.80215
// // // 
// // //         Likelihoods:
// // //                 poisson_likelihood__null_hypothesis__inferred_counts = 1.578555664161381251138685e-21
// // //                 poisson_likelihood__alternative_hypothesis__inferred_counts = 2.313416091357746366358837e-21
// // //                 negative_binomial_likelihood__null_hypothesis__inferred_counts = -2.554321848802513613565818e-18
// // //                 negative_binomial_likelihood__alternative_hypothesis__inferred_counts = -2.573370350779187564020547e-18
// // // 
// // //                 poisson_likelihood__null_hypothesis__UNIQUE_ONLY = 2.318748286475510302095533e-22
// // //                 poisson_likelihood__alternative_hypothesis__UNIQUE_ONLY = 2.318748286475510302095533e-22
// // //                 negative_binomial_likelihood__null_hypothesis__UNIQUE_ONLY = -1.335090308542963771063155e-19
// // //                 negative_binomial_likelihood__alternative_hypothesis__UNIQUE_ONLY = -1.335090308542963771063155e-19
// // // 
// // //         Tails:
// // //                 poisson_tail__null_hypothesis__inferred_counts = 9.743201588094714647703837e-20
// // //                 poisson_tail__alternative_hypothesis__inferred_counts = 1.453598128552905779459525e-19
// // //                 negative_binomial_tail__null_hypothesis__inferred_counts = -4.191186949747015302396115e-15
// // //                 negative_binomial_tail__alternative_hypothesis__inferred_counts = -4.253138012666306516661584e-15
// // // 
// // //                 poisson_tail__null_hypothesis__UNIQUE_ONLY = 1.237464191208852158145896e-20
// // //                 poisson_tail__alternative_hypothesis__UNIQUE_ONLY = 1.237464191208852158145896e-20
// // //                 negative_binomial_tail__null_hypothesis__UNIQUE_ONLY = -1.656467622702920409203024e-16
// // //                 negative_binomial_tail__alternative_hypothesis__UNIQUE_ONLY = -1.656467622702920409203024e-16    
// RD_variance_scale___alpha_squared was about 0.05, I'm pretty sure.
//  pseudo_amount = longdouble base_noise_read_depth_percentage_of_normal_Reference_read_depth = 0.0005; i'm pretty sure 

//     const mpfr_class the_number_one(1.00);		
//     const mpfr_class expected_lambda( (double)(poisson_parameter + pseudo_amount) );
// 	    
//     // CAREFULLY DERIVED VERSION:                                       
//     const mpfr_class parameter_BOOST_r(the_number_one / (double)RD_variance_scale___alpha_squared);
// 	    //BOOST_r  <==>  my_k  <==> Gamma shape parameter k.                                                                
//     const mpfr_class parameter_BOOST_p(the_number_one / (the_number_one + (expected_lambda*(double)RD_variance_scale___alpha_squared)));
// 	    
//     boost::math::negative_binomial_distribution<mpfr_class> my_negative_binomial(parameter_BOOST_r, parameter_BOOST_p);
//     
//     if (!pdf__false__cdf__true)    
// 	return real(boost::math::pdf<mpfr_class>(my_negative_binomial, mpfr_class((double)observation)).get_mpfr_t());
//     else
// 	return real(boost::math::cdf<mpfr_class>(my_negative_binomial, mpfr_class((double)observation)).get_mpfr_t());
}//calculate_negative_binomial_pdf_or_cdf__transforming_from_poisson_parameter




















void  subdivide_intervals
	    (type_map_uint_to_list_BI  &map_chromo_to_intervals,
	     const uint &subdivison_size)
{
    if (subdivison_size > 0)
    {
        for (type_map_uint_to_list_BI::iterator it_chr = map_chromo_to_intervals.begin();
                it_chr != map_chromo_to_intervals.end();
                ++it_chr)        
        {
            type_list_BI new_subdivison;
            
            for (type_list_BI::iterator it_BI = it_chr->second.begin();
                    it_BI != it_chr->second.end();
                    ++it_BI)
            {
                const uint search_width = BOOST_width_inclusive(*it_BI);
                if (  search_width <= subdivison_size )
                    new_subdivison.push_back( *it_BI );
                else
                {
                    const uint number_pieces = (uint)ceil(   ((float)search_width) / ((float)subdivison_size) );
                    const uint piece_size =   (uint)ceil(   ((float)search_width) / ((float)number_pieces) );
                    
                    uint next_lower_bi_ctr = it_BI->lower();
                    
                    while ( next_lower_bi_ctr <= it_BI->upper() )
                    {
                        const uint next_upper_bd =  std::min<uint>(  next_lower_bi_ctr + piece_size - 1, it_BI->upper()  );  // -1 because inclusive
                        new_subdivison.push_back(  BOOST_Interval(next_lower_bi_ctr, next_upper_bd)  );
                        next_lower_bi_ctr = next_upper_bd + 1;                    
                    }                                                        
                }                                                                
            }//bi
            
            it_chr->second = new_subdivison;
        }//chr
                               
    } // subdivide    

}//subdivide_intervals













type_map_uint_to_list_BI__string  upload_list_of_chromosome_regions_from_file
				    (const type_map_uint_to_list_BI  &map_chr_to_list_of_regions)
{
    type_map_uint_to_list_BI__string  desired_pieces_of_chromos;
    
    for (type_map_uint_to_list_BI::const_iterator it_chr = map_chr_to_list_of_regions.begin();
	    it_chr != map_chr_to_list_of_regions.end();
	    ++it_chr)
    {
	    for (type_list_BI::const_iterator it_BI = it_chr->second.begin();
		    it_BI != it_chr->second.end();
		    ++it_BI)
	    {
		const std::string desired_sequence(
					get_specific_DNA_sequence_from_Reference_genome(it_chr->first, *it_BI)  );                                                                                                 
											    
		desired_pieces_of_chromos[ it_chr->first ]
					.push_back( type_BI__string(*it_BI, desired_sequence) );     
	    }
    }
    
    
    return desired_pieces_of_chromos;    
}//upload_list_of_chromosome_regions_from_file














type_map_uint_to_list_BI  pad_all_regions
		(const type_map_uint_to_list_BI &some_map_chr_to_regions,
		 const uint &padding_amount)
{		
    type_map_uint_to_list_BI expanded_map(some_map_chr_to_regions);
    for (type_map_uint_to_list_BI::iterator it_chr = expanded_map.begin();
	    it_chr != expanded_map.end();
	    ++it_chr)
    {
	for (type_list_BI::iterator it_bi = it_chr->second.begin();
		it_bi != it_chr->second.end();
		++it_bi)
	{
	    it_bi->set(safe_subtract_base_1(it_bi->lower(), padding_amount),
		       safe_chromo_add_base_1(it_bi->upper(), padding_amount, it_chr->first));
	}//bi
    }//chr
    
    return  expanded_map;
    
}//pad_all_regions

















type_vector_BI split_interval_into_uniformly_sized_subintervals
		(const BOOST_Interval &big_interval,
		 const uint &number_of_subintervals,
		 const uint &minimum_sublength)
{        
    if (number_of_subintervals == 0)
    {  return type_vector_BI(1,big_interval);  }
    
    type_vector_BI  ordered_subintervals(number_of_subintervals, empty_BI);    
    const uint sublength = BOOST_width_inclusive(big_interval)/number_of_subintervals;
    
    if (sublength == 0  or    (minimum_sublength > 0  and  sublength < minimum_sublength))
    {  return type_vector_BI(1,big_interval);  }
    
    
    uint last_lb = big_interval.lower();
    for (uint sub=0; sub < number_of_subintervals-1; ++sub)
    {
	ordered_subintervals.at(sub).set(big_interval.lower() + sublength*sub,
					 big_interval.lower() + sublength*(sub+1)-1);    
    }
    
    //last
    ordered_subintervals.at(number_of_subintervals-1).set(big_interval.lower() + sublength*(number_of_subintervals-1),
							  big_interval.upper());
    
    return  ordered_subintervals;
}//split_interval_into_uniformly_sized_subintervals


























BOOST_Interval  get_positions_of_hybrid_lying_BELOW_LCR
		    (const type_map_uint_to_uint  &absolute_positions_to_hybrid_indeces,
		     const BOOST_Interval  &some_LCR_coordinates)
{
    BOOST_Interval hybrid_indeces_below_LCR_0;
       
    const type_set_uint  set_hybrid_indeces_below_LCR_0(
		    extract_values_of_map_and_return_as_set<uint,uint>(
			get_map_intersection_of_map_keys_and_interval<uint,uint>(absolute_positions_to_hybrid_indeces,
										BOOST_Interval(0, safe_subtract_base_1(some_LCR_coordinates.lower(),1))  )));
										    
    if (set_hybrid_indeces_below_LCR_0.empty())
	hybrid_indeces_below_LCR_0 = empty_BI;
    else
	hybrid_indeces_below_LCR_0.set(*set_hybrid_indeces_below_LCR_0.begin(), *--set_hybrid_indeces_below_LCR_0.end());
      
    
    return  hybrid_indeces_below_LCR_0;
    
}//get_positions_of_hybrid_lying_BELOW_LCR







BOOST_Interval  get_positions_of_hybrid_lying_ABOVE_LCR
		    (const type_map_uint_to_uint  &absolute_positions_to_hybrid_indeces,
		     const BOOST_Interval  &some_LCR_coordinates,
		     const uint &chromosome)
{    
    
    BOOST_Interval hybrid_indeces_above_LCR_1;
        
    const type_set_uint  set_hybrid_indeces_above_LCR_1(
		    extract_values_of_map_and_return_as_set<uint,uint>(
			get_map_intersection_of_map_keys_and_interval<uint,uint>(absolute_positions_to_hybrid_indeces,
										BOOST_Interval(some_LCR_coordinates.upper()+1,
											       Event::chromosome_lengths.at(chromosome)))));
										    
    if (set_hybrid_indeces_above_LCR_1.empty())
	hybrid_indeces_above_LCR_1 = empty_BI;
    else
	hybrid_indeces_above_LCR_1.set(*set_hybrid_indeces_above_LCR_1.begin(), *--set_hybrid_indeces_above_LCR_1.end());
        
        
    return  hybrid_indeces_above_LCR_1;
   
}//get_positions_of_hybrid_lying_ABOVE_LCR	     






    

BOOST_Interval  get_positions_of_hybrid_lying_INBETWEEN_LCRs
		    (const type_map_uint_to_uint  &absolute_positions_to_hybrid_indeces,
		     const BOOST_Interval  &region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs)
{
    
    BOOST_Interval hybrid_indeces_inbetween_LCRs;
       
    const type_set_uint  set_hybrid_indeces_inbetween_LCRs(
		    extract_values_of_map_and_return_as_set<uint,uint>(
			get_map_intersection_of_map_keys_and_interval<uint,uint>(absolute_positions_to_hybrid_indeces,
										region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs)));
										    
    if (set_hybrid_indeces_inbetween_LCRs.empty())
	hybrid_indeces_inbetween_LCRs = empty_BI;
    else
	hybrid_indeces_inbetween_LCRs.set(*set_hybrid_indeces_inbetween_LCRs.begin(), *--set_hybrid_indeces_inbetween_LCRs.end());
       
    
    return  hybrid_indeces_inbetween_LCRs;            

}//get_positions_of_hybrid_lying_INBETWEEN_LCRs
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
type_list_BI convert_set_to_list_of_contiguous_intervals
			(const type_set_uint &some_set_of_coords)
{
    
    if (some_set_of_coords.empty())
    {  return type_list_BI();  }
    
        
    type_list_BI list_of_BI;
    
    type_set_uint::const_iterator it_check = some_set_of_coords.begin();
    
    uint last_lower = *it_check;
    uint last_upper = *it_check;
    ++it_check;
    
    
    while (it_check != some_set_of_coords.end())
    {
	if (*it_check - last_upper > 1)
	{
	    list_of_BI.push_back(BOOST_Interval(last_lower, last_upper));
	    last_lower = *it_check;
	    last_upper = *it_check;
	    ++it_check;
	}
	else
	{
	    last_upper = *it_check;
	    ++it_check;
	}
    }//it_check
    
    
    //last
    list_of_BI.push_back(BOOST_Interval(last_lower, last_upper));
    
    return  list_of_BI;
    
}//convert_set_to_list_of_contiguous_intervals









BamTools::BamRegion create_BAM_region
				(const uint &some_chromo_value,
				const BOOST_Interval &some_interval,
				const uint interval_padding_amount)
{
    return BamTools::BamRegion( map_chromosome_value_to_BAM_Ref_IDs.at(some_chromo_value),
				safe_subtract_base_1(some_interval.lower(), 1 + interval_padding_amount),   // BAMfiles are 0-based
				map_chromosome_value_to_BAM_Ref_IDs.at(some_chromo_value),
				safe_chromo_add_base_1(safe_subtract_base_1(some_interval.upper(),1), interval_padding_amount, some_chromo_value));  // BAMfiles are 0-based
    
}//create_BAM_region
				
				
				
				
				
				
				
				
				
				
				

				
				
				
				

				

















longdouble get_adjustment_value_from_cumulative_map
			    (const BOOST_Interval  &affected_region,
			     const type_map_uint_to_longdouble &map_positions_to_cumulative_sum)
{    
    if (map_positions_to_cumulative_sum.empty())
    {  return 0.00L;  }

    type_map_uint_to_longdouble::const_iterator LB_cum = map_positions_to_cumulative_sum.lower_bound(affected_region.lower());
    type_map_uint_to_longdouble::const_iterator UB_cum = map_positions_to_cumulative_sum.lower_bound(affected_region.upper());
    //lower_bound()  "Returns an iterator pointing to the first element in the container whose key is not considered to go before k "
    //recall: "lower_bound(K)" returns "end()" if all of map is < K.
    
    
    
    if (LB_cum == map_positions_to_cumulative_sum.end()  // all keys in map  <  affected_region
	or   LB_cum->first > affected_region.upper())  // map-keys are not affected.
    {
	return  0.00L;  
    }	
    else
    {
	//Denote "affected_region" = [A,B]	
	//Then we have:   A  <=  L  <=  B
	
	if (UB_cum == map_positions_to_cumulative_sum.end())
	{  	    	    
	    UB_cum = --map_positions_to_cumulative_sum.end();

	    //Then we now have:
	    //	    A  <=   L   <=  U   <=  B
	    
	    if (LB_cum != map_positions_to_cumulative_sum.begin())
	    {
		--LB_cum;
		//Now:      L  <  A  <=  U  <=  B
		return  UB_cum->second - LB_cum->second;
	    }
	    else
	    {
		//everything below is affected, so subtract all of it.
		return  UB_cum->second;
	    }
	}
	else
	{
	    //Then we now have:
	    //	    A  <=   L   <=  B   <=  U
	    
	    if (UB_cum->first > affected_region.upper())
	    {  
		//Then:    A  <=  L  <=  B  <  U
		--UB_cum; 
		//Now:     A  <=  L  <=  U  <=  B
		
	    }
	    
	    //Now we have:     A  <=  L  <=  U  <=  B	    
	    if (LB_cum != map_positions_to_cumulative_sum.begin())
	    {
		--LB_cum;
		//Now we have:    L  <  A  <=  U  <= B
		return  UB_cum->second - LB_cum->second;		
	    }
	    else
	    {
		//everything below is affected, so subtract all of it.
		return  UB_cum->second;
	    }
	}
    }
    
}//get_adjustment_value_from_cumulative_map



























bool  test_if_outcome_should_display__AB_ABA
			    (const type_haploid_outcome &some_event_outcome)
{
    if (some_event_outcome == hap_outcome__Inv
	    or  some_event_outcome == hap_outcome__Del
	    or  some_event_outcome == hap_outcome__GeneConv_ABA
	    or  some_event_outcome == hap_outcome__Transloc_A)
	return true;
    else
	return false;
    
}//test_if_outcome_should_display__AB_ABA

bool  test_if_outcome_should_display__BA_BAB
			    (const type_haploid_outcome &some_event_outcome)
{
    if (some_event_outcome == hap_outcome__Inv
	    or  some_event_outcome == hap_outcome__Dup
	    or  some_event_outcome == hap_outcome__GeneConv_BAB
	    or  some_event_outcome == hap_outcome__Transloc_B)
	return true;
    else
	return false;
    
}//test_if_outcome_should_display__AB_ABA




























bool test_if_diploid_state_vectors_contain_non_identifiable_diploid_outcomes
		(const Sparse_map &state_vector_hap_0,
		const Sparse_map &state_vector_hap_1)
{
    for (type_map_uint_to_uint::const_iterator it_state_0 = state_vector_hap_0.sparse_map_UID_to_value.begin();
	    it_state_0 != state_vector_hap_0.sparse_map_UID_to_value.end();
	    ++it_state_0)
    {
	const type_haploid_outcome outcome_state_0 = static_cast<type_haploid_outcome>(it_state_0->second);
	const type_haploid_outcome outcome_state_1 = static_cast<type_haploid_outcome>(state_vector_hap_1.get_state_for_UID(it_state_0->first));
	if (    (outcome_state_0 == hap_outcome__Del  and   outcome_state_1 == hap_outcome__Dup)
	    or  (outcome_state_0 == hap_outcome__Dup  and   outcome_state_1 == hap_outcome__Del))
	{  return true;  }		
    }//state 0    
    
    return  false;
    
}//test_if_diploid_state_vectors_contain_non_identifiable_diploid_outcomes





















void set_prior_Event_outcome_probabilities
		    (const real &desired_probability_of_NO_outcome,
		    const real &fraction_of_all_occurring_amount_that_goes_to_GeneConversion_instead_of_NAHR)  
{
    const real the_number_one(1);
    
    P_E_absent_cond_theta = desired_probability_of_NO_outcome;
    
    P_E_GeneConv_onetype_cond_theta = ((the_number_one - P_E_absent_cond_theta)*fraction_of_all_occurring_amount_that_goes_to_GeneConversion_instead_of_NAHR)/2;
    
    P_E_inv_cond_theta = the_number_one - (P_E_absent_cond_theta + P_E_GeneConv_onetype_cond_theta*2);
    
    P_E_dup_or_del_cond_theta = P_E_inv_cond_theta / 2;   
    
    
    std::cerr << "\n\nset prior probabilities to:\n"
	    << "\tP(no event) = " << P_E_absent_cond_theta.toString(15) << "\n"
	    << "\tP(gene conversion \"ABA\") = " << P_E_GeneConv_onetype_cond_theta.toString(15) << "\n"
	    << "\tP(inversion) = " << P_E_inv_cond_theta.toString(15) << "\n"
	    << "\tP(del) = " << P_E_dup_or_del_cond_theta.toString(15) << "\n\n";
    
}//set_prior_Event_outcome_probabilities


































type_multimap_uint_to_uint  get_sparseness_of_variational_positions
				    (const type_set_uint  &some_set_of_vp)
{    
    
    type_multimap_uint_to_uint multimap_sparseness_to_vp;
    
    if (some_set_of_vp.empty())
    {  return  multimap_sparseness_to_vp;  }    
    
    //first
    multimap_sparseness_to_vp.insert(type_uint__uint(*some_set_of_vp.begin(),  *some_set_of_vp.begin()));  //from first vp to beginning of LCR - the distance is = to the index itself!!!
    
    for (type_set_uint::const_iterator it_vp_below = some_set_of_vp.begin(),
					    it_vp = ++some_set_of_vp.begin();              
	    it_vp != some_set_of_vp.end();
	    ++it_vp_below, ++it_vp)
    {
	multimap_sparseness_to_vp.insert(type_uint__uint(*it_vp - *it_vp_below, *it_vp));
    }
    
	
    return  multimap_sparseness_to_vp;
    
}//get_sparseness_of_variational_positions








type_multimap_uint_to_uint::const_iterator  find_sparsest_vp_less_or_greater_than_value
						    (const type_multimap_uint_to_uint &multimap_sparseness_to_vp,
						    const uint &some_vp,
						    const bool &less_than__false____greater__true)
{        
    
    if (less_than__false____greater__true)
    {
	for (type_multimap_uint_to_uint::const_reverse_iterator it_sp_to_vp = multimap_sparseness_to_vp.rbegin();
		it_sp_to_vp != multimap_sparseness_to_vp.rend();
		++it_sp_to_vp)
	{
	    if (it_sp_to_vp->second > some_vp)
	    {
		return   --(it_sp_to_vp.base());	
	    }	    	    
	}	
    }//greater
    else 
    {
	for (type_multimap_uint_to_uint::const_reverse_iterator it_sp_to_vp = multimap_sparseness_to_vp.rbegin();
		it_sp_to_vp != multimap_sparseness_to_vp.rend();
		++it_sp_to_vp)
	{
	    if (it_sp_to_vp->second < some_vp)
	    {
		return   --(it_sp_to_vp.base());
	    }
	}			
    }//less
    
    
    return  multimap_sparseness_to_vp.end();
    
    
}//find_sparsest_vp_less_or_greater_than_value






void   add_sparse_vp_to_form_Geneconversion_intervals
		(const Event &the_Event,
		const uint &maximum_number_of_allowable_breakpoints,
		const type_map_uint_to_uint &first_brkpt_to_count,
		const type_map_uint_to_uint &second_brkpt_to_count,
		type_set_BI  &set_of_GeneConv_breakpoints)
{
    
    if (set_of_GeneConv_breakpoints.size() >= maximum_number_of_allowable_breakpoints)
    {  return;  }
    
    
    const type_multimap_uint_to_uint  multimap_sparseness_to_vp(
			get_sparseness_of_variational_positions(the_Event.uploaded_variational_positions_of_profile));
    
    
    type_multimap_uint_to_uint multimap_count_to_vp(invert_map_into_multimap<uint,uint>(first_brkpt_to_count));
    {
	type_multimap_uint_to_uint mm_second_brkpt(invert_map_into_multimap<uint,uint>(second_brkpt_to_count));	
	multimap_count_to_vp.insert(mm_second_brkpt.begin(), mm_second_brkpt.end());
    }
    
    for (type_multimap_uint_to_uint::const_reverse_iterator it_count_to_vp = multimap_count_to_vp.rbegin();
	    it_count_to_vp != multimap_count_to_vp.rend();
	    ++it_count_to_vp)
    {
	const type_multimap_uint_to_uint::const_iterator it_sparsest_to_vp
			    =  find_sparsest_vp_less_or_greater_than_value(
						multimap_sparseness_to_vp,
						it_count_to_vp->second,
						(first_brkpt_to_count.count(it_count_to_vp->second) > 0) );
			    			    
	if (it_sparsest_to_vp != multimap_sparseness_to_vp.end())
	{
	    set_of_GeneConv_breakpoints.insert(BOOST_Interval(std::min<uint>(it_sparsest_to_vp->second, it_count_to_vp->second),
							      std::max<uint>(it_sparsest_to_vp->second, it_count_to_vp->second)  ));
	    
	    if (set_of_GeneConv_breakpoints.size() >= maximum_number_of_allowable_breakpoints)
	    {  break;  }	    
	}
    }
    
    
    if (set_of_GeneConv_breakpoints.empty())
    {
	set_of_GeneConv_breakpoints.insert(BOOST_Interval(0, the_Event.last_profile_nongap_index));
    }//empty        
    
    
}//add_sparse_vp_to_form_Geneconversion_intervals







void  preliminary_refine_breakpoints_for_inversion_heuristic
		(const Event &the_Event,
		const uint &maximum_number_of_allowable_breakpoints,
		const type_map_uint_to_uint &map_varpos_to_number_of_decent_reads__SAVE,
		type_set_uint  &refined_varpos,
		type_map_uint_to_uint  &refined_map_varpos_to_count)
		
{
    refined_map_varpos_to_count = remove_values_less_than_or_equal_to_threshold_from_map<uint,uint>(map_varpos_to_number_of_decent_reads__SAVE,heuristic_posterior_COUNT_threshold);        
    refined_varpos = extract_keys_of_map_and_return_as_set<uint,uint>(refined_map_varpos_to_count);    
    
    prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained<uint,uint>(
				    refined_varpos,
				    invert_map_into_multimap<uint,uint>(refined_map_varpos_to_count),
				    maximum_number_of_allowable_breakpoints,
				    refined_map_varpos_to_count);
    
    
}//preliminary_refine_breakpoints_for_inversion_heuristic










type_set_uint refine_breakpoints_using_heuristic
		(const Event &the_Event,
		const uint &maximum_number_of_allowable_breakpoints,
		const type_map_uint_to_uint &map_varpos_to_number_of_decent_reads__SAVE)
{
    type_map_uint_to_uint  map_varpos_to_number_of_decent_reads = remove_values_less_than_or_equal_to_threshold_from_map<uint,uint>(map_varpos_to_number_of_decent_reads__SAVE,heuristic_posterior_COUNT_threshold);
    
    if (map_varpos_to_number_of_decent_reads.empty())
    {
	map_varpos_to_number_of_decent_reads = remove_zero_values_from_map<uint,uint>(map_varpos_to_number_of_decent_reads__SAVE);
    }
    
    type_set_uint supported_vp = extract_keys_of_map_and_return_as_set<uint,uint>(map_varpos_to_number_of_decent_reads);
		    
    if (supported_vp.size() > maximum_number_of_allowable_breakpoints)
    {      		
	prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained<uint,uint>(
					supported_vp,
					invert_map_into_multimap<uint,uint>(map_varpos_to_number_of_decent_reads),
					maximum_number_of_allowable_breakpoints,
					map_varpos_to_number_of_decent_reads);
    }
    else if (supported_vp.size() < 2)
    {//add some sparse varpos
    
	if (the_Event.my_original_Breakpoint_complex.variational_positions__ALL.size() == 0)
	{
	    supported_vp.insert(0);
	    supported_vp.insert(the_Event.last_profile_nongap_index);
	}
	else if (the_Event.my_original_Breakpoint_complex.variational_positions__ALL.size() == 1)
	{
	    supported_vp = the_Event.my_original_Breakpoint_complex.variational_positions__ALL;
	    
	    if (the_Event.large_gaps_of_profile.size() == 1)
	    {
		supported_vp.insert(the_Event.large_gaps_of_profile.begin()->lower());
		
		if (supported_vp.size() == 1)
		{  supported_vp.insert(the_Event.last_profile_nongap_index);  }			
		if (supported_vp.size() == 1)
		{  supported_vp.insert(0);  }
	    }
	    else
	    {             
		supported_vp.insert(the_Event.last_profile_nongap_index);
		if (supported_vp.size() == 1)
		{  supported_vp.insert(0);  }
	    }                
	}//the_Event.my_original_Breakpoint_complex.variational_positions__ALL.size() == 1
	else if (the_Event.my_original_Breakpoint_complex.variational_positions__ALL.size() == 2) // >= 2 original vp
	{
	    supported_vp = the_Event.my_original_Breakpoint_complex.variational_positions__ALL;
	}
	else  // the_Event.my_original_Breakpoint_complex.variational_positions__ALL.size() >= 3
	{                     
	    type_multimap_uint_to_uint multimap_sparseness_to_vp(	    
					    get_sparseness_of_variational_positions(
						    the_Event.my_original_Breakpoint_complex.variational_positions__ALL));	    	    
	    
	    type_multimap_uint_to_uint::const_reverse_iterator it_most_sparse_vp = multimap_sparseness_to_vp.rbegin();
	    assert(multimap_sparseness_to_vp.size() > 2);
	    
	    while (supported_vp.size() < 2  and  it_most_sparse_vp != multimap_sparseness_to_vp.rend())
	    {
		supported_vp.insert(it_most_sparse_vp->second);
		++it_most_sparse_vp;
	    }		        	    
						    
	}//add some sparse varpos                 
    }//varpos
    
    
    //just in case I'm incompetent...
    if (supported_vp.size() < 2)
    {
	std::cerr << "\n\nERROR in \"refine_breakpoints_using_heuristic\" - made it past with < 2 supported_vp!!!\n";
	supported_vp.insert(0);
	supported_vp.insert(the_Event.last_profile_nongap_index);
    }
        
    
    
    return  supported_vp;
    
    
}//refine_breakpoints_using_heuristic	    

























type_set_BI  group_breakpoints_according_to_relevant_variational_position_constraints
					    (const type_set_uint &relevant_variational_positions,
					    const type_set_BI &the_space_of_potential_breakpoints)
{
    //When summing out Event X, not all of it's neighbors' breakpoints will be relevant for computing the likelihood of the data relevant to Event X (i.e., the portions of the genome that must be considered when summing out Event X).  Given a set of relevant variational positions, we take the set of ALL possible variational positions, and form sets ("groups", intervals) of variational positons according to the relevancy in considering the data associated with Event X.
    //For example, suppose Event Y has breakpoints  b1, ..., b5.  Suppose only b1 and b4 must be considered when summing out Event X.  Then we formt eh intervals [b1,b1], [b2,b4], [b5,b5].  As far as tghe data for Event X is concerned, a breakpoint in Event Y at b2,b3 or b4 is indistinguishable.  When we later sum over Event Y, then of course b2, b3, b4 will actually have distinguishable effects on the data for Event Y, and so we must "untangle"/"syntheisze"/"resolve" the breakpoint intervals from Event X with Event Y.  (this happens in Partial_sum_function).             
    
    
   //Strategy:
    //For every original interval, just "grow" the set to its maximum potential.  That is:
    //		For the right endpoint:
    //			1) if the right endpoint is already a relevant varpos, then stop.
    //			2) if the right endpoint is not already a relevant varpos, then move it as far to the right as possible without hitting a relevant varpos.
    //		Analogous for the left endpoint.
    //		If neither endpoints is a relevant varpo,s then grow on both sides.
    
    
    
    type_set_BI  efficient_breakpoints_for_computation;
    
    //doens't matter if the variational position is meaningful to the type of outcome being considered, or if is Gene Converison left or right endpoint - rememebr, it is how it affects the genome!!!  
    
    for (type_set_BI::const_iterator it_potential = the_space_of_potential_breakpoints.begin();
	    it_potential != the_space_of_potential_breakpoints.end();
	    ++it_potential)
    {
	const type_set_uint rel_vp_intersected_by_potential(get_intersection_of_set_and_interval(relevant_variational_positions,*it_potential));
	BOOST_Interval grown_brkpt_BI(*it_potential);
	
	if (relevant_variational_positions.count(it_potential->lower()) == 0  and  relevant_variational_positions.count(it_potential->upper()) == 0)
	{//neither
	
	    type_set_BI merges_space(the_space_of_potential_breakpoints);	    	    	    
	    type_set_BI::const_iterator it_cand = merges_space.begin();
	    
	    while (it_cand != merges_space.end())
	    {
		if (BOOST_overlap(grown_brkpt_BI, *it_cand)
		    and  (relevant_variational_positions.count(it_cand->lower()) == 0  and  relevant_variational_positions.count(it_cand->upper()) == 0)
		    and  (get_intersection_of_set_and_interval(relevant_variational_positions, *it_cand) == rel_vp_intersected_by_potential))
		{
		    grown_brkpt_BI = BOOST_hull(grown_brkpt_BI, *it_cand);    
		    merges_space.erase(it_cand);
		    it_cand = merges_space.begin();
		}				
		else
		{  ++it_cand;  }
	    }//2
	}//neither
	else if (relevant_variational_positions.count(it_potential->lower()) > 0   and    relevant_variational_positions.count(it_potential->upper()) == 0)
	{//left is fixed - right can be moved
	
	    type_set_BI merges_space(the_space_of_potential_breakpoints);
	    type_set_BI::const_iterator it_cand = merges_space.begin();
	    
	    while (it_cand != merges_space.end())
	    {
		if (BOOST_overlap(grown_brkpt_BI, *it_cand)
		    and  (it_cand->lower() == it_potential->lower()   and   relevant_variational_positions.count(it_cand->upper()) == 0)
		    and  (get_intersection_of_set_and_interval(relevant_variational_positions, *it_cand) == rel_vp_intersected_by_potential))
		{
		    grown_brkpt_BI = BOOST_hull(grown_brkpt_BI, *it_cand);    
		    merges_space.erase(it_cand);
		    it_cand = merges_space.begin();
		}				
		else
		{  ++it_cand;  }
	    }//2
	}//left is fixed - right can be moved
	else if (relevant_variational_positions.count(it_potential->lower()) > 0   and    relevant_variational_positions.count(it_potential->upper()) == 0)
	{//right is fixed - left can be moved
	
	    type_set_BI merges_space(the_space_of_potential_breakpoints);
	    type_set_BI::const_iterator it_cand = merges_space.begin();
	    
	    while (it_cand != merges_space.end())
	    {
		if (BOOST_overlap(grown_brkpt_BI, *it_cand)
		    and  (relevant_variational_positions.count(it_cand->lower()) == 0   and   it_cand->upper() == it_potential->upper())
		    and  (get_intersection_of_set_and_interval(relevant_variational_positions, *it_cand) == rel_vp_intersected_by_potential))
		{
		    grown_brkpt_BI = BOOST_hull(grown_brkpt_BI, *it_cand);    
		    merges_space.erase(it_cand);
		    it_cand = merges_space.begin();
		}				
		else
		{  ++it_cand;  }
	    }//2
	}//right is fixed - left can be moved
	
	
	//add
	efficient_breakpoints_for_computation.insert(grown_brkpt_BI);		
	
    }//the_space_of_potential_breakpoints
    
    
    
    return  efficient_breakpoints_for_computation;
    
    
}//group_breakpoints_according_to_relevant_variational_position_constraints



























void determine_relevant_NAHR_and_GeneConversion_breakpoints_for_each_interlocking_event_and_self_according_to_relevant_variational_positions
                                                        (const Conn_Comp &the_Conn_Comp,
							 const type_map_uint_to_set_uint &relevant_varpos_per_Event,
							type_map_uint_to_Breakpoint_complex &relevant_breakpoints_per_Event)
{
    //When summing out Event X, not all of it's neighbors' breakpoints will be relevant for computing the likelihood of the data relevant to Event X (i.e., the portions of the genome that must be considered when summing out Event X).  Given a set of relevant variational positions, we take the set of ALL possible variational positions, and form sets ("groups", intervals) of variational positons according to the relevancy in considering the data associated with Event X.
    //For example, suppose Event Y has breakpoints  b1, ..., b5.  Suppose only b1 and b4 must be considered when summing out Event X.  Then we formt eh intervals [b1,b1], [b2,b4], [b5,b5].  As far as tghe data for Event X is concerned, a breakpoint in Event Y at b2,b3 or b4 is indistinguishable.  When we later sum over Event Y, then of course b2, b3, b4 will actually have distinguishable effects on the data for Event Y, and so we must "untangle"/"syntheisze"/"resolve" the breakpoint intervals from Event X with Event Y.  (this happens in Partial_sum_function).        

    for (type_map_uint_to_set_uint::const_iterator it_ev = relevant_varpos_per_Event.begin();
            it_ev != relevant_varpos_per_Event.end();
            ++it_ev)
    {
        const type_map_uint_to_Event::const_iterator Event_being_considered = the_Conn_Comp.events.find(it_ev->first);  //readability.
	
	//NAHR
	relevant_breakpoints_per_Event[Event_being_considered->first].breakpoints_NAHR__AB
			    = group_breakpoints_according_to_relevant_variational_position_constraints(
				    it_ev->second, Event_being_considered->second.my_original_Breakpoint_complex.breakpoints_NAHR__AB);
			    
	relevant_breakpoints_per_Event[Event_being_considered->first].breakpoints_NAHR__BA
			    = group_breakpoints_according_to_relevant_variational_position_constraints(
				    it_ev->second, Event_being_considered->second.my_original_Breakpoint_complex.breakpoints_NAHR__BA);
                                                                                        
	//GeneConversion
	relevant_breakpoints_per_Event[Event_being_considered->first].breakpoints_GeneConversion___ABA
			    = group_breakpoints_according_to_relevant_variational_position_constraints(
				    it_ev->second, Event_being_considered->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA);
			    
	relevant_breakpoints_per_Event[Event_being_considered->first].breakpoints_GeneConversion___BAB
			    = group_breakpoints_according_to_relevant_variational_position_constraints(
				    it_ev->second, Event_being_considered->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB);
                                                                                                                                                     
    }  // end for-loop   relevant_varpos_per_Event    
                
}//determine_relevant_NAHR_and_GeneConversion_breakpoints_for_each_interlocking_event_and_self_according_to_relevant_variational_positions































void redo_all_marginal_MAP_calls_using_saved_samples_from_posterior	
		(const std::string &basedir,
		const type_map_uint_to_CC &Conn_Comps)
{
    
    const boost::filesystem::directory_iterator it_name_end;    
    for (boost::filesystem::directory_iterator it_name(basedir);
	    it_name != it_name_end;
	    ++it_name)
    {
	if (!boost::filesystem::is_directory(*it_name))
	{  continue;  }
	
	output_dir = boost::filesystem::complete(it_name->path()).string();
	
	{
	    type_vector_string splitted_name;	    
	    boost::split(splitted_name, it_name->path().filename().string(), boost::is_any_of("_"), boost::algorithm::token_compress_on);
	    genome_name = splitted_name.at(1);
	}
	
	std::cerr << "\n\t\t\toutput_dir = [" << output_dir << "]\n";
	std::cerr << "\n\t\t\tgenome_name = [" << genome_name << "]\n";
	
	{//rename
	    std::string original_raw_results(output_dir);
	    original_raw_results.append("/raw_posterior_marginal_MAP_calls");
	    std::string old_raw_results(output_dir);
	    old_raw_results.append("/OLD_raw_posterior_marginal_MAP_calls");
	    
	    if (boost::filesystem::exists(original_raw_results))
	    {  boost::filesystem::rename(original_raw_results, old_raw_results);  }
	}//rename
	
	std::string samplesdir(output_dir);
	samplesdir.append("/Samples_from_posterior");	
	boost::filesystem::path samples_path(samplesdir);
	
	const boost::filesystem::directory_iterator it_samp_end;
	for (boost::filesystem::directory_iterator it_sample(samples_path);
		it_sample != it_samp_end;
		++it_sample)
	{
	    
	    type_map_uint_to_CC::const_iterator the_cc_it;	    
	    {
		//ex: "samples_from_posterior_CC_29041__number_of_samples_1000"
		type_vector_string splitted_name;	    
		boost::split(splitted_name, it_sample->path().filename().string(), boost::is_any_of("_"), boost::algorithm::token_compress_on);
		
		const std::string cc_as_str(splitted_name.at(4));	    
		const int in_cc_id = atoi(cc_as_str.c_str());
		std::cerr << "\t\t\t\t\tin_cc_id = " << in_cc_id << "\n";
		
		
		the_cc_it = Conn_Comps.find(in_cc_id);
		assert(the_cc_it != Conn_Comps.end());
	    }
	    
	    
	    std::cerr << "\t\t" << boost::filesystem::complete(*it_sample).string() << "\n";
	    	    
	    const type_vector_map_uint_to_Sampled_diploid_Event_data   loaded_samples_from_posterior(load_posterior_samples_from_file(boost::filesystem::complete(*it_sample).string()));
	    
	    std::cerr << "\tmarginals...\n";
	    const type_map_uint_to_Marginal_Event_posterior_distribution marginal_distributions
							(the_cc_it->second.compute_marginal_probabilities_from_posterior_sample(loaded_samples_from_posterior));      
		
	    std::stringstream output_marginal_Positive_calls_strm;
	    std::stringstream output_marginal_Negative_calls_strm;
	    std::stringstream output_marginal_all_calls_strm;
	    
	    std::cerr << "\tmake MAP calls...\n";
	    make_MAP_calls_from_marginal_posterior_distributions_per_Event
					(marginal_distributions,
					the_cc_it->second,
					output_marginal_Positive_calls_strm,
					output_marginal_Negative_calls_strm,
					output_marginal_all_calls_strm);
	}//sample
    }//name
    
    
}//redo_marginal_MAP_calls_using_saved_samples_from_posterior





















longdouble get_standard_deviation_from_mean_for_RD
				(const longdouble &mean)
{
    return  std::max<longdouble>(3,(sqrt(mean*(1.00 + sqrt(RD_variance_scale___alpha_squared)))));    
}//get_standard_deviation_from_mean_for_RD































void redo_all_inferred_RD_tests_and_make_ratios    
				(const std::string &postdatatdir,
				 type_map_uint_to_CC &Conn_Comps)
{
    
    print_line_of_markers("==");
    std::cerr << " \n\n\n\n\n\nredo_all_inferred_RD_tests_and_make_ratios:\n\n";
	
	
    std::string calls_fname(postdatatdir);
    calls_fname.append("/all_calls");
    
    type_map_string_to_uint_to_Sampled_diploid_Event_data map_genome_to_event_to_sampled_Event(read_Sampled_Events_from_file(calls_fname));
    
    for (type_map_string_to_uint_to_Sampled_diploid_Event_data::reverse_iterator it_gen = map_genome_to_event_to_sampled_Event.rbegin();
	    it_gen != map_genome_to_event_to_sampled_Event.rend();
	    ++it_gen)
    {
	
	genome_name = it_gen->first;
	
	std::string full_job_dirname;
	{//full_job_dirname
	    const boost::filesystem::directory_iterator it_end;  
	    for (boost::filesystem::directory_iterator it_jobdir(postdatatdir);
		    it_jobdir != it_end;
		    ++it_jobdir)
	    {		
		if (boost::filesystem::is_directory(it_jobdir->path()))
		{
		    if (it_jobdir->path().filename().string().find(it_gen->first) != std::string::npos)
		    {
			full_job_dirname = boost::filesystem::complete(it_jobdir->path()).string();
			break;			
		    }		    		    
		}//dir				
	    }//it_jobdir	    	    	    
	}//full_job_dirname
	
	
	std::cerr << "\n\t\tfull_job_dirname = [" << full_job_dirname << "]\n";
	
	output_dir = full_job_dirname;
	
	
	
	std::string bam_filename;	
	{//full_job_dirname
	    const boost::filesystem::directory_iterator it_end;  
	    for (boost::filesystem::directory_iterator it_bamdir("/gpfs/scratch/mmparks/SAMs_and_BAMs");
		    it_bamdir != it_end;
		    ++it_bamdir)
	    {		
		if (!boost::filesystem::is_directory(it_bamdir->path()))
		{
		    if (it_bamdir->path().filename().string().find(it_gen->first) != std::string::npos
			and  it_bamdir->path().filename().string().find(".bas") == std::string::npos
			and  it_bamdir->path().filename().string().find(".bai") == std::string::npos
			and  it_bamdir->path().filename().string().find(".bam") == it_bamdir->path().filename().string().size() - 4)			
		    {
			bam_filename = boost::filesystem::complete(it_bamdir->path()).string();
			break;			
		    }		    		    
		}//file
	    }//it_jobdir	    	    	    
	}//full_job_dirname	
	
	std::cerr << "\t\tbam_filename = [" << bam_filename << "]\n\n";
	
	
	
	assert(!full_job_dirname.empty());
	assert(!bam_filename.empty());			
	
	
	
	read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda(bam_filename, false, NULL);
	
	
	const std::vector<type_map_uint_to_Sampled_diploid_Event_data::iterator>  myev__loop(
					get_loopable_iterators<type_map_uint_to_Sampled_diploid_Event_data>(it_gen->second));
	
	#pragma omp parallel for schedule(dynamic,10)
	for (int myev_ctr = my_MPI_rank; myev_ctr < myev__loop.size(); myev_ctr += size_MPI_WORLD)	
// 	for (type_map_uint_to_Sampled_diploid_Event_data::iterator it_myev = it_gen->second.begin();
// 		it_myev != it_gen->second.end();
// 		++it_myev)
	{    
	    const type_map_uint_to_Sampled_diploid_Event_data::iterator it_myev = myev__loop.at(myev_ctr);
	    	    
	    it_myev->second.print_this_Sampled_diploid_Event_data();	    
	    
	    const type_map_uint_to_Event::iterator result_Ev_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID); 
	    
	    
	    if (result_Ev_it->second.chromos[0] == 24  and   gender_of_individual_is_female)
	    {
		std::cerr << "Event is on Y chromosome, but individual is female.  Skipping read-depth graph in in \"create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights\".\n";
		continue;
	    }
	    
	    const BOOST_Interval expanded_region_of_interest( 
				    safe_subtract_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.lower(), 100000),
				    safe_chromo_add_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.upper(), 100000, result_Ev_it->second.chromos[0]));        
				    
	    //identify all events intersecting this region:
	    type_map_uint_to_Event_ptr  all_events_near_expanded_region;
	    for (type_map_uint_to_Event::iterator it_ev = Conn_Comps.at(it_myev->second.cc_of_event).events.begin();
		    it_ev != Conn_Comps.at(it_myev->second.cc_of_event).events.end();
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
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    char out_RD_core_fname[option_size];
	    
	    std::sprintf(out_RD_core_fname, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s_%s___%s", 
					output_dir.c_str(),
					it_myev->second.cc_of_event, result_Ev_it->second.UID,
					convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first).c_str(),
					convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second).c_str(),
					genome_name.c_str());		    
			    
	    
	    //if already computed, skip)
	    {//check already computed
		std::string testname(out_RD_core_fname);
		testname.append(".tests");
		
		if (boost::filesystem::exists(testname) != 0)
		{
		    std::ifstream infs(testname);
		    
		    std::string inteststr(option_size*10, '\0');
		    char intestfile[option_size*10];
		    std::sprintf(intestfile, "%s", inteststr.c_str());
		    		    
		    infs.read(intestfile, option_size*10);
		    infs.close();
		    
		    inteststr.assign(intestfile);
		    const uint num_newlines = std::count(inteststr.begin(), inteststr.end(), '\n');		    
		    if (num_newlines >= 31)
		    {
			std::cerr << "\n\tnum_newlines = " << num_newlines << " !!!  continue!\n";
			continue;			
		    }
		}//does not exit
	    }//check already computed
	    
	    
	    
	    
	    std::string observed_fname(out_RD_core_fname);
	    observed_fname.append(".data.observed");	    
	    
	    if (boost::filesystem::exists(observed_fname) == 0)
	    {//search
		std::string searchdir(output_dir);
		searchdir.append("/Read_Depth");
						
		std::stringstream cc_strm;
		cc_strm << "cc_" << it_myev->second.cc_of_event << "_";
		std::stringstream ev_strm;
		ev_strm << "ev_" << it_myev->second.event_UID << "_";
		
		char observed_fname_old[option_size];
		
		bool found_something = false;
		const boost::filesystem::directory_iterator it_end;  
		for (boost::filesystem::directory_iterator it_rdfile(searchdir);
			it_rdfile != it_end;
			++it_rdfile)
		{		
		    if (!boost::filesystem::is_directory(it_rdfile->path()))
		    {
			if (it_rdfile->path().filename().string().find(cc_strm.str()) != std::string::npos
			    and  it_rdfile->path().filename().string().find(ev_strm.str()) != std::string::npos
			    and  it_rdfile->path().filename().string().find(genome_name) != std::string::npos
			    and it_rdfile->path().filename().string().find("RD_hypothesis_test__cc_") != std::string::npos)
			{
			    
			   
			    char scanformat[option_size];
			    std::sprintf(scanformat, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%%s", 
						    output_dir.c_str(),
						    it_myev->second.cc_of_event,
						    it_myev->second.event_UID);	
			    
			    char in_calls_and_genome[option_size];		    
			    			    
// 			    std::cerr << "\n\t\t\t\ttry name: [" << boost::filesystem::complete(it_rdfile->path()).string().c_str() << "\n\t\t\t\tscanformat = [ " << scanformat << "]\n";
			    std::sscanf(boost::filesystem::complete(it_rdfile->path()).string().c_str(), scanformat,
					    in_calls_and_genome);
			    
			    std::string calls_CHOPPED(in_calls_and_genome);
			    calls_CHOPPED = calls_CHOPPED.substr(0,calls_CHOPPED.find(genome_name));
			    
// 			    std::cerr << "\n\t\t\t\tcalls_CHOPPED = [" << calls_CHOPPED << "]\n";
			    
			    std::sprintf(observed_fname_old, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s%s.data.observed", 
					output_dir.c_str(),
					it_myev->second.cc_of_event, result_Ev_it->second.UID,
					calls_CHOPPED.c_str(),
					genome_name.c_str());	    
			    
			    found_something = true;
			    break;			
			}		    		    
		    }//file
		}//it_jobdir
		
	    
		if (!found_something)
		{
		    std::stringstream error_strm;
		    error_strm << "ERROR  unable to find anything!!!\n\n";
		    error_message(error_strm,false);
		    continue;
		}	    
	    
		std::cerr << "\n\t\t\tout_RD_core_fname = [" << out_RD_core_fname << "]\n"
			<<  "\t\t\tobserved_fname_old = [" << observed_fname_old << "\n\n";
								
		if (boost::filesystem::exists(observed_fname) == 0)
		{  
		    std::cerr << "rename = [" << observed_fname << "\n";
		    boost::filesystem::rename(observed_fname_old, observed_fname);		    
		}
		
	    }//search

		    
	    
	    
			
	    //get "inferred" observed values.
// 	    std::string observed_fname(out_RD_core_fname);
// 	    observed_fname.append(".data.observed");	  
// 	    
// 	    if (boost::filesystem::exists(observed_fname) == 0)
// 	    {
// 		std::swap<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome.first, it_myev->second.the_diploid_Event_outcome.second);
// 		std::swap<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts.first, it_myev->second.the_diploid_profile_brkpts.second);
// 		
// 		std::sprintf(out_RD_core_fname, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s_%s___%s", 
// 					output_dir.c_str(),
// 					it_myev->second.cc_of_event, result_Ev_it->second.UID,
// 					convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first).c_str(),
// 					convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second).c_str(),
// 					genome_name.c_str());	    
// 		
// 		observed_fname = out_RD_core_fname;
// 		observed_fname.append(".data.observed");
// 				
// 	    }
	    
	    const type_map_uint_to_longdouble map_Ref_Genome_position_to_per_read_posterior_sum(read_map_from_file<uint,longdouble>(observed_fname, false));
	
	    
	    
	    
	
	    
	    //get expected values for different hypotheses...
	    const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble map_diploid_outcome_to_absolute_position_to_expected_RD(
					construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region(
								    result_Ev_it->second.chromos[0],
								    expanded_region_of_interest,
								    it_myev->second.the_diploid_Event_outcome,
								    it_myev->second.the_diploid_profile_brkpts,
								    result_Ev_it->second));    
	    
	    
	    
	    

					
					
// 	    {//save
// 		std::string expected_fname(out_RD_core_fname);
// 		expected_fname.append(".data.expected");
// 		
// 		std::ofstream out_expected_fs(expected_fname);
// 		if (check_fstream_for_goodness(out_expected_fs, expected_fname, "create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights", false))
// 		{
// 		    const type_haploid_outcome__haploid_outcome alternative_dip_state(it_myev->second.the_diploid_Event_outcome);
// 		    std::cerr << "alternative_dip_state = " << convert_haploid_outcome_to_string(alternative_dip_state.first) << ", " << convert_haploid_outcome_to_string(alternative_dip_state.second) << "\n";
// 		    
// 		    const type_haploid_outcome__haploid_outcome null_dip_state(type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None));
// 		    
// 		    
// 		    const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_null_dip_state_expected_map 
// 			= map_diploid_outcome_to_absolute_position_to_expected_RD.find(null_dip_state);
// 		    const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_alternative_dip_state_expected_map 
// 			= map_diploid_outcome_to_absolute_position_to_expected_RD.find(alternative_dip_state);
// 					    
// 		    const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_default_del_expected_map 
// 			= map_diploid_outcome_to_absolute_position_to_expected_RD.find(type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__Del));
// 		    const type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble::const_iterator it_default_dup_expected_map 
// 			= map_diploid_outcome_to_absolute_position_to_expected_RD.find(type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__Dup));   
// 					    
// 		    
// 
// 		    if (it_null_dip_state_expected_map != map_diploid_outcome_to_absolute_position_to_expected_RD.end()
// 			and it_alternative_dip_state_expected_map != map_diploid_outcome_to_absolute_position_to_expected_RD.end()
// 			and it_default_del_expected_map != map_diploid_outcome_to_absolute_position_to_expected_RD.end() 
// 			and it_default_dup_expected_map != map_diploid_outcome_to_absolute_position_to_expected_RD.end()  )
// 		    {
// 			for (type_map_uint_to_longdouble::const_iterator it_null = it_null_dip_state_expected_map->second.begin(),
// 									it_expect = it_alternative_dip_state_expected_map->second.begin(),
// 									it_def_del = it_null_dip_state_expected_map->second.begin(),
// 									it_def_dup = it_null_dip_state_expected_map->second.begin();
// 				it_null != map_diploid_outcome_to_absolute_position_to_expected_RD.at(null_dip_state).end();
// 				++it_null, ++it_expect)
// 			{ 
// 			    out_expected_fs 
// 				<< it_null->first << "\t" 
// 				<< it_null->second << "\t"
// 				<< it_expect->second << "\t"
// 				<< it_def_del->second << "\t"
// 				<< it_def_dup->second << "\n";
// 			}
// 		    }
// 		    else
// 		    {  std::cerr << "\n\tunrecognized alternative dip sates!!!\n\n";  }
// 		    
// 		    out_expected_fs.close();	
// 		}
// 	    }//save
// 
// 		
		
	    std::string nonunique_fname(out_RD_core_fname);
	    nonunique_fname.append(".data.nonunique");
	    type_set_uint  nonunique_positions_in_region_of_interest;
	    if (boost::filesystem::exists(nonunique_fname) == 0)
	    {
			nonunique_positions_in_region_of_interest = get_nonunique_positions_in_interval(
								    load_Non_unique_regions_of_genome(),
								    result_Ev_it->second.chromos[0], expanded_region_of_interest);
	    }
	    else
	    {
			nonunique_positions_in_region_of_interest = read_set_from_file<uint>(nonunique_fname);
	    }


		

		
	    
	    {//RD hypothesis test
	    
		BOOST_Interval affected_region(empty_BI);
		
		if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None))
		{
		    affected_region = result_Ev_it->second.region_between_and_including_the_LCRs_themselves;	    
		}
		else
		{
		    for (uint hap=0; hap<2; ++hap)
		    {
			if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) != hap_outcome__None)
			{
			    const type_BI__BI affected_regions_each_LCR(
						    BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(
							    convert_each_profile_breakpoint_to_absolute_breakpoints(
								    pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
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
			= map_diploid_outcome_to_affected_region_absolute_position_to_expected_RD.find(it_myev->second.the_diploid_Event_outcome);
			
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
		    continue;
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
		
		
		
		
		
// 		if (!test_if_haploid_outcome_is_GeneConversion(it_myev->second.the_diploid_Event_outcome.first) and  !test_if_haploid_outcome_is_GeneConversion(it_myev->second.the_diploid_Event_outcome.second))
// 		{//save to common file
// 		    //not GeneConv
// 		    
// 		    std::string commonfname(output_dir);
// 		    
// 		    if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None))	    
// 			commonfname.append("/all_RD_tests.null");
// 		    else
// 			commonfname.append("/all_RD_tests.NAHR");		
// 					
// 		    std::ofstream outfs(commonfname, std::ios::app);
// 		    
// 		    if (check_fstream_for_goodness(outfs, commonfname, "redo_all_inferred_RD_tests_and_make_ratios", false))
// 		    {
// 			outfs
// 	// 		    << genome_name << "\t"
// 	// 		    << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first) << "\t"
// 	// 		    << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second) << "\t"
// 			    << "\n" << odds_ratio__alt___inferred.toString(10)
// 			    << "\n" << odds_ratio__alt___unique.toString(10)
// 			    
// 			    << "\n" << odds_ratio__def_del___inferred.toString(10)
// 			    << "\n" << odds_ratio__def_del___unique.toString(10)
// 			    
// 			    << "\n" << odds_ratio__def_dup___inferred.toString(10)
// 			    << "\n" << odds_ratio__def_dup___unique.toString(10)
// 			    
// 			    << "\n" << odds_ratio__best_default__inferred.toString(10)
// 			    << "\n" << odds_ratio__best_default__unique.toString(10);
// 		    }//good
// 		    
// 		    outfs.close();
// 		}//save to common file
// 	    
	    
	    }//RD hypothesis test
	    
	}//it_myev
    }//it_gen
            		

}//redo_all_inferred_RD_tests_and_make_ratios








type_set_uint get_nonunique_positions_in_interval
		    (const type_map_uint_to_list_BI &non_uniq_regions_of_genome,
		     const uint &chromo,
		     const BOOST_Interval &expanded_region_of_interest)
{
    type_set_uint  nonunique_positions_in_region_of_interest;
    type_set_uint::iterator it_insert_nonunique = nonunique_positions_in_region_of_interest.begin();
        
    type_map_uint_to_list_BI::const_iterator it_chr_non = non_uniq_regions_of_genome.find(chromo);	
    if (it_chr_non != non_uniq_regions_of_genome.end() )
    {
	for (type_list_BI::const_iterator it_non = it_chr_non->second.begin();
		it_non != it_chr_non->second.end();
		++it_non)
	{
	    const BOOST_Interval nonuniq_intersection(BOOST_intersect(*it_non, expanded_region_of_interest));
	    if (!BOOST_empty(nonuniq_intersection))
	    {
		for (uint pos = nonuniq_intersection.lower(); pos <= nonuniq_intersection.upper(); ++pos)
		{
		    it_insert_nonunique = nonunique_positions_in_region_of_interest.insert(it_insert_nonunique, pos);
		}
	    }//overlap
	}//bi
    }//found chr    
    
    return  nonunique_positions_in_region_of_interest;
    
}//get_nonunique_positions_in_interval


































type_map_BI_to_longdouble bin_by_mapykey
		(const type_map_uint_to_longdouble &some_map,
		const uint &binsize)
{	
	if (some_map.empty())
	{  return type_map_BI_to_longdouble();  }
	
	type_map_BI_to_longdouble binned_map;
	
	
	type_map_uint_to_longdouble::const_iterator it_map = some_map.begin();	
	
	uint starting_window_position = it_map->first;
	uint last_window_position = starting_window_position;
	uint current_size = 1;
	longdouble current_bin_aggregate_sum = it_map->second;
	++it_map;
	
	while (it_map != some_map.end())
	{
		if (current_size < binsize)
		{
			current_bin_aggregate_sum += it_map->second;
			last_window_position = it_map->first;
			++current_size;			
		}
		else
		{
			binned_map.insert(std::pair<BOOST_Interval, longdouble>
									(BOOST_Interval(starting_window_position, last_window_position), current_bin_aggregate_sum));
			
			starting_window_position = it_map->first;
			last_window_position = starting_window_position;
			current_size = 1;
			current_bin_aggregate_sum = it_map->second;	
		}//add, start anew.
		
		++it_map;
	}//it_map
	
// 	//any remaining windows.
// 	if (current_size > 0) 
// 	{
// 		binned_map.insert(std::pair<BOOST_Interval, longdouble>
// 						(BOOST_Interval(starting_window_position, last_window_position), current_bin_aggregate_sum));
// 	}
	
	
	return binned_map;
	
}//bin_by_mapykey
























real  perform_Wilcoxon_signed_rank_test
				(const uint &binsize,
				const type_map_uint_to_longdouble &map_position_to_expected_RD___NULL,
				const type_map_uint_to_longdouble &map_position_to_Observed_RD,
				const Sampled_diploid_Event_data &some_result)
{
	//Wilcoxon signed rank test:  http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
	//
	// Bin according to consecutive non-overlapping windows of length M.  Say M = 50.
	// Windows are therefore independent.
	// We form D(x) = B(x) - N(x),   where x is a window, B(x) is the observed read count in window x, and N(x) is the expected read-depth in window x under the Null distribution.
	// We then perform a one-sided Wilcoxon signed rank test on D.  The null hypothesis is that the median of D is 0.  If we suspect a deletion, then the alternative is one-sided, i.e. median < 0.  If the call was null, then the test is two-sided.  We obtain a p-value from this.
		
	if (map_position_to_expected_RD___NULL.empty() or map_position_to_Observed_RD.empty())
	{  return  -1;  }
	
	
	typedef  std::pair<longdouble, bool>  type_longdouble__bool;
	typedef  std::list<type_longdouble__bool> type_list_longdouble__bool;
	
	//extract only the region in quesiton.
	const type_map_BI_to_longdouble binned_observations(bin_by_mapykey(map_position_to_Observed_RD, binsize));
	const type_map_BI_to_longdouble binned_expected_null(bin_by_mapykey(map_position_to_expected_RD___NULL, binsize));
	
	type_map_BI_to_longdouble binned_difference(binned_observations);
	{
		type_map_BI_to_longdouble::iterator sub_it = binned_difference.begin();
		for (type_map_BI_to_longdouble::const_iterator it_null = binned_expected_null.begin();
			it_null != binned_expected_null.end();
			++it_null)
		{
			sub_it->second -= it_null->second;			
			++sub_it;
		}//null
	}
	
	
	
	type_list_longdouble__bool rankings;
	for (type_map_BI_to_longdouble::const_iterator it_diff = binned_difference.begin();
		it_diff != binned_difference.end();
		++it_diff)
	{			
		rankings.push_back(type_longdouble__bool(std::abs<longdouble>(it_diff->second), (it_diff->second > 0)));			
	}//diff
	
	rankings.sort();
	
	int observed_signed_ranksum = 0;
	uint rank_number = 1;
	for (type_list_longdouble__bool::const_iterator it_rank = rankings.begin();
		it_rank != rankings.end();
		++it_rank)
	{
		if (it_rank->second)
			observed_signed_ranksum += rank_number;			
		else
			observed_signed_ranksum -= rank_number;	
		
		++rank_number;
	}//rank
	
	
	
	
	// CLT approximation:
	// under H_0, observed_signed_ranksum has mean = 0 and variance = n(n + 1)(2n + 1)/6    where   n = number of elements in ranking.
	const uint number_of_elts_in_ranking = rankings.size();
	const double varnormal = (double)(number_of_elts_in_ranking*(number_of_elts_in_ranking+1)*(number_of_elts_in_ranking*2 + 1))/6;
	boost::math::normal_distribution<mpfr_class> my_standard_normal(0,1);
	
	const double CLT_transformed_observed = (double)observed_signed_ranksum/std::sqrt(varnormal);
	
	std::cerr 
		<< "\t\tn = " << rankings.size() << "\n"
		<< "\t\tobserved_signed_ranksum = " << observed_signed_ranksum << "\n"				
		<< "\t\tvarnormal = " << varnormal 
		<< "\t\tCLT_transformed_observed = " << CLT_transformed_observed << "\n";
	
	real tail_prob(0);
	
	if (some_result.the_diploid_Event_outcome.first == hap_outcome__Del
		or some_result.the_diploid_Event_outcome.second == hap_outcome__Del)		
	{
		tail_prob = real(boost::math::cdf<mpfr_class>(my_standard_normal, mpfr_class(CLT_transformed_observed)).get_mpfr_t());
	}
	else if (some_result.the_diploid_Event_outcome.first == hap_outcome__Dup
		or some_result.the_diploid_Event_outcome.second == hap_outcome__Dup)		
	{
		tail_prob = real(boost::math::cdf<mpfr_class>(my_standard_normal, mpfr_class(CLT_transformed_observed*(-1))).get_mpfr_t());
	}
	else
	{
		const double abs_CLT_transform = std::abs(CLT_transformed_observed);
		tail_prob = real(boost::math::cdf<mpfr_class>(my_standard_normal, mpfr_class(abs_CLT_transform*(-1))).get_mpfr_t())*2;
	}
			
	
	return tail_prob;	

}//perform_Wilcoxon_signed_rank_test















uint randomly_choose_bin_of_same_GC_count_and_count_observed
			(const uint &chr_of_current_bin,
			const BOOST_Interval &region_of_current_bin,
			const type_map_uint_to_list_BI &NON_unique_regions_of_genome,
			BamTools::BamReader &my_BAM_reader)
{
	const uint skip_amount_at_end_of_chromos = 10000;
	const uint preupload_size = 10000;
	
	const uint bin_width = BOOST_width_inclusive(region_of_current_bin);
	const uint desired_GC_count = calculate_GC_content_of_sequence(
										get_specific_DNA_sequence_from_Reference_genome(chr_of_current_bin, region_of_current_bin));	
	
	std::cerr << "desired_GC_count = " << desired_GC_count << "\n";
	
	boost::random::mt19937 random_generator;
	boost::random::uniform_int_distribution<> uniform_chromosome_dist(1, 22);		
	
	
	
	uint chr_of_random_bin = uniform_chromosome_dist(random_generator);
	BOOST_Interval region_of_random_bin;			
	
	
	std::cerr << "chr_of_random_bin = " << chr_of_random_bin << "\n";			
	
	const BOOST_Interval acceptable_p_arm(skip_amount_at_end_of_chromos, Event::centromere_coordinates.at(chr_of_random_bin).first - skip_amount_at_end_of_chromos);
	
	const BOOST_Interval acceptable_q_arm(Event::centromere_coordinates.at(chr_of_random_bin).second + skip_amount_at_end_of_chromos,
										Event::chromosome_lengths.at(chr_of_random_bin) - skip_amount_at_end_of_chromos);
	
	boost::random::uniform_int_distribution<> uniform_position_dist(1, Event::chromosome_lengths.at(chr_of_random_bin));
	
	
	
	
	bool good_draw = false;
	uint num_attempts = 0;
	while (!good_draw)
	{				
		bool found_good_GC_content = false;
		while (!found_good_GC_content)
		{
			++num_attempts;
			const uint low_endpoint = uniform_position_dist(random_generator);
			region_of_random_bin.set(low_endpoint, low_endpoint + bin_width - 1);
			
			if (BOOST_subset(region_of_random_bin, acceptable_p_arm)  or  BOOST_subset(region_of_random_bin, acceptable_q_arm))
			{		
				std::string uploaded_big_region;
			
				#pragma omp critical(upload_big_region_for_checking_wilcoxon_test_matched)
				uploaded_big_region = get_specific_DNA_sequence_from_Reference_genome(chr_of_random_bin, BOOST_Interval(low_endpoint, low_endpoint + preupload_size));		
			
				for (uint start = low_endpoint; start < low_endpoint + preupload_size - bin_width; ++start)
				{
					const std::string test_seq(uploaded_big_region.substr(start - low_endpoint, bin_width));
					const uint gc_count__random = calculate_GC_content_of_sequence(test_seq);
					if (gc_count__random == desired_GC_count)
					{
						found_good_GC_content = true;
						region_of_random_bin.set(low_endpoint, low_endpoint + bin_width - 1);
						break;
					}
				}//start
																					
			}
		}//found_good_GC_content
				
		
		std::cerr << "\n\t\tattempting  " << chr_of_random_bin << ": " << region_of_random_bin.lower() << ", " << region_of_random_bin.upper() << ".";
		
		if (BOOST_subset(region_of_random_bin, acceptable_p_arm)  or  BOOST_subset(region_of_random_bin, acceptable_q_arm))
		{
			bool is_unique = true;
			if (NON_unique_regions_of_genome.count(chr_of_random_bin) > 0)
			{
				const type_map_uint_to_list_BI::const_iterator it_chr_nonuniq = NON_unique_regions_of_genome.find(chr_of_random_bin);
				for (type_list_BI::const_iterator it_bi = it_chr_nonuniq->second.begin();
					it_bi != it_chr_nonuniq->second.end();
					++it_bi)
				{
					if (BOOST_overlap(*it_bi, region_of_random_bin))
					{ 
						is_unique = false;
						break;
					}										
				}//bi
			}
			
			if (is_unique)
			{
				if (chr_of_current_bin != chr_of_random_bin
					or !BOOST_overlap(region_of_current_bin, region_of_random_bin))
				{
					good_draw = true;
				}								
			}						
		}//acceptable arms				
	}//while
	
	
	uint number_observed_reads;
	
	#pragma omp critical(wilcoxon_randomly_select_control_set_region_and_count)
	{
		my_BAM_reader.SetRegion(create_BAM_region(chr_of_random_bin, region_of_random_bin, 3));	
		number_observed_reads = count_decent_reads_with_left_endpoints_in_set_region(my_BAM_reader, region_of_random_bin);		
	}
	
	return number_observed_reads;
	
}//randomly_choose_bin_of_same_GC_count_and_count_observed











real  perform_Wilcoxon_signed_rank_test
				(const uint &binsize,
				 const uint &chromo_of_observed,
				const type_map_uint_to_longdouble &map_position_to_Observed_RD,
				const Sampled_diploid_Event_data &some_result,
				const type_map_uint_to_list_BI &NON_unique_regions_of_genome,
				BamTools::BamReader &my_BAM_reader)				
{
	//Wilcoxon signed rank test:  http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
	//
	// Bin according to consecutive non-overlapping windows of length M.  Say M = 50.
	// Windows are therefore independent.
	// We form D(x) = B(x) - N(x),   where x is a window, B(x) is the observed read count in window x, and N(x) is the expected read-depth in window x under the Null distribution.
	// We then perform a one-sided Wilcoxon signed rank test on D.  The null hypothesis is that the median of D is 0.  If we suspect a deletion, then the alternative is one-sided, i.e. median < 0.  If the call was null, then the test is two-sided.  We obtain a p-value from this.
		
	if (map_position_to_Observed_RD.empty())
	{  return  -1;  }
	
	
	typedef  std::pair<longdouble, bool>  type_longdouble__bool;
	typedef  std::list<type_longdouble__bool> type_list_longdouble__bool;
	
	//extract only the region in quesiton.
	const type_map_BI_to_longdouble binned_observations(bin_by_mapykey(map_position_to_Observed_RD, binsize));
	type_map_BI_to_longdouble binned_controls;
	
// 	type_map_uint_to_list_BI NON_unique_regions_of_genome(load_Non_unique_regions_of_genome());
	
	for (type_map_BI_to_longdouble::const_iterator it_bin = binned_observations.begin();
		it_bin != binned_observations.end();
		++it_bin)
	{
		binned_controls[it_bin->first] 
						= randomly_choose_bin_of_same_GC_count_and_count_observed
										(chromo_of_observed,
										it_bin->first,//make_BI(get_endpoints_of_keys_of_map<uint, longdouble>(map_position_to_Observed_RD)),
										NON_unique_regions_of_genome,
										my_BAM_reader);				
	}
	
	
	type_map_BI_to_longdouble binned_difference(binned_observations);
	{
		type_map_BI_to_longdouble::iterator sub_it = binned_difference.begin();
		for (type_map_BI_to_longdouble::const_iterator it_null = binned_controls.begin();
			it_null != binned_controls.end();
			++it_null)
		{
			sub_it->second -= it_null->second;			
			++sub_it;
		}//null
	}
	
	
	
	type_list_longdouble__bool rankings;
	for (type_map_BI_to_longdouble::const_iterator it_diff = binned_difference.begin();
		it_diff != binned_difference.end();
		++it_diff)
	{			
		rankings.push_back(type_longdouble__bool(std::abs<longdouble>(it_diff->second), (it_diff->second > 0)));			
	}//diff
	
	rankings.sort();
	
	int observed_signed_ranksum = 0;
	uint rank_number = 1;
	for (type_list_longdouble__bool::const_iterator it_rank = rankings.begin();
		it_rank != rankings.end();
		++it_rank)
	{
		if (it_rank->second)
			observed_signed_ranksum += rank_number;			
		else
			observed_signed_ranksum -= rank_number;	
		
		++rank_number;
	}//rank
	
	
	
	
	// CLT approximation:
	// under H_0, observed_signed_ranksum has mean = 0 and variance = n(n + 1)(2n + 1)/6    where   n = number of elements in ranking.
	const uint number_of_elts_in_ranking = rankings.size();
	const double varnormal = (double)(number_of_elts_in_ranking*(number_of_elts_in_ranking+1)*(number_of_elts_in_ranking*2 + 1))/6;
	boost::math::normal_distribution<mpfr_class> my_standard_normal(0,1);
	
	const double CLT_transformed_observed = (double)observed_signed_ranksum/std::sqrt(varnormal);
	
	std::cerr 
		<< "\t\tn = " << rankings.size() << "\n"
		<< "\t\tobserved_signed_ranksum = " << observed_signed_ranksum << "\n"				
		<< "\t\tvarnormal = " << varnormal 
		<< "\t\tCLT_transformed_observed = " << CLT_transformed_observed << "\n";
	
	real tail_prob(0);
	
	if (some_result.the_diploid_Event_outcome.first == hap_outcome__Del
		or some_result.the_diploid_Event_outcome.second == hap_outcome__Del)		
	{
		tail_prob = real(boost::math::cdf<mpfr_class>(my_standard_normal, mpfr_class(CLT_transformed_observed)).get_mpfr_t());
	}
	else if (some_result.the_diploid_Event_outcome.first == hap_outcome__Dup
		or some_result.the_diploid_Event_outcome.second == hap_outcome__Dup)		
	{
		tail_prob = real(boost::math::cdf<mpfr_class>(my_standard_normal, mpfr_class(CLT_transformed_observed*(-1))).get_mpfr_t());
	}
	else
	{
		const double abs_CLT_transform = std::abs(CLT_transformed_observed);
		tail_prob = real(boost::math::cdf<mpfr_class>(my_standard_normal, mpfr_class(abs_CLT_transform*(-1))).get_mpfr_t())*2;
	}
			
	
	return tail_prob;	

}//perform_Wilcoxon_signed_rank_test












type_map_uint_to_list_uint get_random_sample_of_bins_from_across_the_genome
								(const uint &bin_width,
								const uint &number_of_samples, //mult of 10
								BamTools::BamReader &my_BAM_reader)	
{

	type_map_uint_to_list_uint sampled_observed_counts_per_GC_content;	

	std::cerr << "load_Non_unique_regions_of_genome...\n";
	const type_map_uint_to_list_BI NON_unique_regions_of_genome(load_Non_unique_regions_of_genome());	
	
	const uint skip_amount_at_end_of_chromos = 10000;
	boost::random::mt19937 random_generator;
	boost::random::uniform_int_distribution<> uniform_chromosome_dist(1, 22);		
	
	for (uint samp=0; samp<number_of_samples; ++samp)
	{
		std::cerr << "\tsamp = " << samp << "\n";
		const uint chr_of_random_bin = uniform_chromosome_dist(random_generator);								
		
		const BOOST_Interval acceptable_p_arm(skip_amount_at_end_of_chromos, Event::centromere_coordinates.at(chr_of_random_bin).first - skip_amount_at_end_of_chromos);
		
		const BOOST_Interval acceptable_q_arm(Event::centromere_coordinates.at(chr_of_random_bin).second + skip_amount_at_end_of_chromos,
											Event::chromosome_lengths.at(chr_of_random_bin) - skip_amount_at_end_of_chromos);
		
		boost::random::uniform_int_distribution<> uniform_position_dist(1, Event::chromosome_lengths.at(chr_of_random_bin));
		
			
		bool good_draw = false;
		uint num_attempts = 0;
		while (!good_draw)
		{				
			++num_attempts;
			const uint low_endpoint = uniform_position_dist(random_generator);
			const BOOST_Interval region_of_random_bin(low_endpoint, low_endpoint + bin_width - 1);		
			
			if (BOOST_subset(region_of_random_bin, acceptable_p_arm)  or  BOOST_subset(region_of_random_bin, acceptable_q_arm))
			{		
				bool is_unique = true;
				if (NON_unique_regions_of_genome.count(chr_of_random_bin) > 0)
				{
					const type_map_uint_to_list_BI::const_iterator it_chr_nonuniq = NON_unique_regions_of_genome.find(chr_of_random_bin);
					for (type_list_BI::const_iterator it_bi = it_chr_nonuniq->second.begin();
						it_bi != it_chr_nonuniq->second.end();
						++it_bi)
					{
						if (BOOST_overlap(*it_bi, region_of_random_bin))
						{ 
							is_unique = false;
							break;
						}										
					}//bi
				}//chr
				
				
				if (is_unique)
				{			
					std::string test_seq;
				
					#pragma omp critical(random_sample_for_controls_for_Wilcoxon_tests)
					test_seq = get_specific_DNA_sequence_from_Reference_genome(chr_of_random_bin, region_of_random_bin);		
				
					const uint gc_count__random = calculate_GC_content_of_sequence(test_seq);		
					
					uint number_observed_reads;
					
					#pragma omp critical(random_sample_for_controls_for_Wilcoxon_tests)
					{
						my_BAM_reader.SetRegion(create_BAM_region(chr_of_random_bin, region_of_random_bin, 3));	
						number_observed_reads = count_decent_reads_with_left_endpoints_in_set_region(my_BAM_reader, region_of_random_bin);		
					}				
					
					good_draw = true;
					sampled_observed_counts_per_GC_content[gc_count__random].push_back(number_observed_reads);				
					
				}//is_unique			
			}//acceptable arms				
		}//while
		
	}//samp
	
	
	return  sampled_observed_counts_per_GC_content;	
	
}//get_random_sample_of_bins_from_across_the_genome


















int randomly_choose_matching_control_from_subsample
			(const uint &desired_GC_count,
			 const type_map_uint_to_list_uint &sampled_control_bins_by_GC_count)
{
	const type_map_uint_to_list_uint::const_iterator it_control_GC = sampled_control_bins_by_GC_count.find(desired_GC_count);
	if (it_control_GC != sampled_control_bins_by_GC_count.end())
	{
		const uint num_subsamples = it_control_GC->second.size();
		assert(num_subsamples > 0);
		
		boost::random::mt19937 random_generator;
		boost::random::uniform_int_distribution<> uniform_sample_dist(1, num_subsamples);
		
		const uint sample_num = uniform_sample_dist(random_generator);
		
		type_list_uint::const_iterator it_sampled_control = it_control_GC->second.begin();
		uint ctr = 1;
		while (ctr < sample_num)
		{
			++ctr;
			++it_sampled_control;			
		}
		
		return *it_sampled_control;
	}
	else
	{
		std::cerr << "\n\nERROR!  unable to find GC count = " << desired_GC_count << "  from subsample!\n\n";
		return -1;
// 		exit(1);
	}
	
	
}//randomly_choose_matching_control_from_subsample
















void  perform_Wilcoxon_signed_rank_test
				(const uint &binsize,
				 const uint &chromo_of_observed,
				const type_map_uint_to_longdouble &map_position_to_Observed_RD,
				const Sampled_diploid_Event_data &some_result,
				 const type_map_uint_to_list_uint &sampled_control_bins_by_GC_count,
				const longdouble &sampled_control_multiplier,
				const type_map_uint_to_list_BI  &nonunique_regions_of_genome,
				Wilcoxon_calculations &observed_wilcoxon_stats___inferred,
				Wilcoxon_calculations &observed_wilcoxon_stats__unique)
{
	//Wilcoxon signed rank test:  http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
	//
	// Bin according to consecutive non-overlapping windows of length M.  Say M = 50.
	// Windows are therefore independent.
	// We form D(x) = B(x) - N(x),   where x is a window, B(x) is the observed read count in window x, and N(x) is the expected read-depth in window x under the Null distribution.
	// We then perform a one-sided Wilcoxon signed rank test on D.  The null hypothesis is that the median of D is 0.  If we suspect a deletion, then the alternative is one-sided, i.e. median < 0.  If the call was null, then the test is two-sided.  We obtain a p-value from this.
	
	const bool verbose = false;
		
	if (map_position_to_Observed_RD.empty())
	{
		return;		
	}
	

	
	const BOOST_Interval region_of_interest(make_BI(get_endpoints_of_keys_of_map<uint,longdouble>(map_position_to_Observed_RD)));
	const std::string region_of_interest_seq(get_specific_DNA_sequence_from_Reference_genome(chromo_of_observed, region_of_interest));				
	
	//extract only the region in question.
	type_map_BI_to_longdouble binned_observations(bin_by_mapykey(map_position_to_Observed_RD, binsize));
// 	if (BOOST_width_inclusive((--binned_observations.end())->first) < binsize)
// 	{  binned_observations.erase(--binned_observations.end());  }
	
	type_map_BI_to_longdouble binned_controls;	
	{
		type_map_BI_to_longdouble::iterator it_bin = binned_observations.begin();
		while (it_bin != binned_observations.end())
		{
			const std::string binned_string(region_of_interest_seq.substr(it_bin->first.lower() - region_of_interest.lower(), binsize));
			const uint num_N = std::count(binned_string.begin(), binned_string.end(), 'N');
			
			if (num_N > 1)
			{
				binned_observations.erase(it_bin++);	
			}
			else
			{
				const uint bin_gc = calculate_GC_content_of_sequence(binned_string);
				const int random_control = randomly_choose_matching_control_from_subsample(bin_gc, sampled_control_bins_by_GC_count);
				
				if (random_control >= 0)
				{
					binned_controls[it_bin->first] = (uint)random_control;
	// 				if (verbose)
	// 				{
	// 					std::cerr << "\tgc=" << bin_gc;
	// 				}
					++it_bin;
				}
				else
				{
					std::cerr << "ERROR! random_control = " << random_control << ",  bin_gc = " << bin_gc << "\n"
							<< "substr = [" << region_of_interest_seq.substr(it_bin->first.lower() - region_of_interest.lower(), binsize) << "].\n"
							<< "substring length = " << region_of_interest_seq.substr(it_bin->first.lower() - region_of_interest.lower(), binsize).size() << "\n";
					
					binned_observations.erase(it_bin++);
				}//error
			}//else
		}//while
	}//binned_controls
	
	
	
	// adjust if you want to test alternative.	
	for (type_map_BI_to_longdouble::iterator it_control = binned_controls.begin();
		it_control != binned_controls.end();
		++it_control)
	{
		it_control->second *= sampled_control_multiplier;
	}
	
	
	
	if (verbose)
	{
		std::stringstream debugstrm;
		debugstrm << "\n\n\nbinned observations and controls:\n";
		for (type_map_BI_to_longdouble::const_iterator it_obs = binned_observations.begin(),
				it_cont = binned_controls.begin();
			it_obs != binned_observations.end();
			++it_obs, ++it_cont)			
		{
			debugstrm << it_obs->first << "\t" << it_obs->second << "\t" << it_cont->second << "\n";						
		}
		
		std::cerr << debugstrm.str() << "\n";
		
// 		print_map_keys_and_values<BOOST_Interval, longdouble>(binned_observations, "binned_observations", NULL, true);
// 		print_map_keys_and_values<BOOST_Interval, longdouble>(binned_controls, "binned_controls", NULL, true);
	}//v
		
		
// 	assert(binned_controls.size() == binned_observations.size());
	
	type_map_BI_to_longdouble binned_difference(binned_observations);
	{
		type_map_BI_to_longdouble::iterator sub_it = binned_difference.begin();
		for (type_map_BI_to_longdouble::const_iterator it_null = binned_controls.begin();
			it_null != binned_controls.end();
			++it_null)
		{
			sub_it->second -= it_null->second;			
			++sub_it;
		}//null
	}
	
	
	
	type_map_BI_to_longdouble binned_difference___unique(binned_difference);
	{		
		type_map_uint_to_list_BI::const_iterator it_non_chr = nonunique_regions_of_genome.find(chromo_of_observed);
		if (it_non_chr != nonunique_regions_of_genome.end())
		{
			type_list_BI nonuniq_in_region_of_interest;					
			for (type_list_BI::const_iterator it_non_bi = it_non_chr->second.begin();
				it_non_bi != it_non_chr->second.end();
				++it_non_bi)
			{
				
				if (BOOST_overlap(*it_non_bi, region_of_interest))
				{
					nonuniq_in_region_of_interest.push_back(*it_non_bi);
				}
			}//it_non_bi
			
			for (type_list_BI::const_iterator it_non_bi = nonuniq_in_region_of_interest.begin();
				it_non_bi != nonuniq_in_region_of_interest.end();
				++it_non_bi)
			{
				type_map_BI_to_longdouble::iterator it_bin_diff = binned_difference___unique.begin();
				while (it_bin_diff != binned_difference___unique.end())
				{
					if (BOOST_overlap(it_bin_diff->first, *it_non_bi))
					{
						binned_difference___unique.erase(it_bin_diff++);
					}	
					else if (it_bin_diff->first.lower() > it_non_bi->upper())
					{  break;  }
					else
					{  ++it_bin_diff;  }
				}//it_bin_diff
			}//it_non_bi			
		}//chr		
	}//unique
	
	
	
	compute_signed_rank_sum(binned_difference, verbose, observed_wilcoxon_stats___inferred.signed_rank_sum_value, 
														observed_wilcoxon_stats___inferred.number_of_elements_in_ranking);

	compute_signed_rank_sum(binned_difference___unique, verbose, observed_wilcoxon_stats__unique.signed_rank_sum_value,
																observed_wilcoxon_stats__unique.number_of_elements_in_ranking);	
		
			
	
	// CLT approximation:
	// under H_0, observed_wilcoxon_stats___inferred.signed_rank_sum_value has mean = 0 and variance = n(n + 1)(2n + 1)/6    where   n = number of elements in ranking.
	const double varnormal__inferred = (double)(observed_wilcoxon_stats___inferred.number_of_elements_in_ranking*(observed_wilcoxon_stats___inferred.number_of_elements_in_ranking+1)*(observed_wilcoxon_stats___inferred.number_of_elements_in_ranking*2 + 1))/6;
	boost::math::normal_distribution<mpfr_class> my_standard_normal___inferred(0,1);
	
	observed_wilcoxon_stats___inferred.variance_of_normal = varnormal__inferred;
	
	const double varnormal__unique = (double)(observed_wilcoxon_stats__unique.number_of_elements_in_ranking*(observed_wilcoxon_stats__unique.number_of_elements_in_ranking+1)*(observed_wilcoxon_stats__unique.number_of_elements_in_ranking*2 + 1))/6;
	boost::math::normal_distribution<mpfr_class> my_standard_normal___unique(0,1);		
	
	observed_wilcoxon_stats__unique.variance_of_normal = varnormal__unique;
	
	const double CLT_transformed_observed__inferred = (double)(observed_wilcoxon_stats___inferred.signed_rank_sum_value/std::sqrt(varnormal__inferred));
	const double CLT_transformed_observed__unique = (double)(observed_wilcoxon_stats__unique.signed_rank_sum_value/std::sqrt(varnormal__unique));
	
		
	observed_wilcoxon_stats___inferred.wilcoxon_likelihood = real(boost::math::pdf<mpfr_class>(my_standard_normal___inferred,
																	   mpfr_class(CLT_transformed_observed__inferred)).get_mpfr_t());
	
	observed_wilcoxon_stats__unique.wilcoxon_likelihood = real(boost::math::pdf<mpfr_class>(my_standard_normal___unique,
																	   mpfr_class(CLT_transformed_observed__unique)).get_mpfr_t());	
	
	
	
	std::cerr 
		<< "\tinferred:\n"
		<< "\t\tn = " << observed_wilcoxon_stats___inferred.number_of_elements_in_ranking << "\n"
		<< "\t\tobserved_signed_ranksum = " << observed_wilcoxon_stats___inferred.signed_rank_sum_value << "\n"				
		<< "\t\tvarnormal = " << varnormal__inferred  << "\n"
		<< "\t\tCLT_transformed_observed = " << CLT_transformed_observed__inferred << "\n"
		
		<< "\tunique:\n"
		<< "\t\tn = " << observed_wilcoxon_stats__unique.number_of_elements_in_ranking << "\n"
		<< "\t\tobserved_signed_ranksum = " << observed_wilcoxon_stats__unique.signed_rank_sum_value << "\n"				
		<< "\t\tvarnormal = " << varnormal__unique << "\n"
		<< "\t\tCLT_transformed_observed = " << CLT_transformed_observed__unique << "\n";		
	

	
	if (some_result.the_diploid_Event_outcome.first == hap_outcome__Del
		or some_result.the_diploid_Event_outcome.second == hap_outcome__Del)		
	{
		if (observed_wilcoxon_stats___inferred.number_of_elements_in_ranking > 0)
		{
			observed_wilcoxon_stats___inferred.wilcoxon_pvalue = real(boost::math::cdf<mpfr_class>(my_standard_normal___inferred, mpfr_class(CLT_transformed_observed__inferred)).get_mpfr_t());
		}
		
		if (observed_wilcoxon_stats__unique.number_of_elements_in_ranking > 0)		
		{  
			observed_wilcoxon_stats__unique.wilcoxon_pvalue = real(boost::math::cdf<mpfr_class>(my_standard_normal___unique, mpfr_class(CLT_transformed_observed__unique)).get_mpfr_t());  			
		}
	}
	else if (some_result.the_diploid_Event_outcome.first == hap_outcome__Dup
		or some_result.the_diploid_Event_outcome.second == hap_outcome__Dup)		
	{
		if (observed_wilcoxon_stats___inferred.number_of_elements_in_ranking > 0)
		{		
			observed_wilcoxon_stats___inferred.wilcoxon_pvalue = real(boost::math::cdf<mpfr_class>(my_standard_normal___inferred, mpfr_class(CLT_transformed_observed__inferred*(-1))).get_mpfr_t());
		}
		
		if (observed_wilcoxon_stats__unique.number_of_elements_in_ranking > 0)
		{
			observed_wilcoxon_stats__unique.wilcoxon_pvalue = real(boost::math::cdf<mpfr_class>(my_standard_normal___unique, mpfr_class(CLT_transformed_observed__unique*(-1))).get_mpfr_t());  			
		}
	}
	else
	{
		if (observed_wilcoxon_stats___inferred.number_of_elements_in_ranking > 0)
		{		
			const double abs_CLT_transform___inferred = std::abs(CLT_transformed_observed__inferred);
			observed_wilcoxon_stats___inferred.wilcoxon_pvalue = real(boost::math::cdf<mpfr_class>(my_standard_normal___inferred, mpfr_class(abs_CLT_transform___inferred*(-1))).get_mpfr_t())*2;
		}
		
		if (observed_wilcoxon_stats__unique.number_of_elements_in_ranking > 0)
		{
			const double abs_CLT_transform__unique = std::abs(CLT_transformed_observed__unique);
			observed_wilcoxon_stats__unique.wilcoxon_pvalue = real(boost::math::cdf<mpfr_class>(my_standard_normal___unique, mpfr_class(abs_CLT_transform__unique*(-1))).get_mpfr_t())*2;
		}
	}
						

}//perform_Wilcoxon_signed_rank_test












void compute_signed_rank_sum
				(const type_map_BI_to_longdouble &binned_difference,
				 const bool &verbose,
				 longdouble &observed_signed_ranksum,
				uint &number_of_elts_in_ranking)
{
	if (binned_difference.size() < 2)
	{
		observed_signed_ranksum = -1;
		number_of_elts_in_ranking = 0;
		return;		
	}
	
	type_list_longdouble__bool rankings;
	for (type_map_BI_to_longdouble::const_iterator it_diff = binned_difference.begin();
		it_diff != binned_difference.end();
		++it_diff)
	{			
		rankings.push_back(type_longdouble__bool(std::abs<longdouble>(it_diff->second), (it_diff->second > 0)));			
	}//diff
	
	
	
	rankings.sort();
	
	{//eliminate 0's
		type_list_longdouble__bool::iterator it_rank = rankings.begin();
		while (it_rank != rankings.end() and std::round(it_rank->first) == 0)  
		{
			it_rank = rankings.erase(it_rank);						
		}						
	}//eliminate 0's
	
	
	
	observed_signed_ranksum = 0;
	{//observed_signed_ranksum
		uint rank_number = 1;
		
		type_list_longdouble__bool::const_iterator it_rank = rankings.begin();
		while (it_rank != rankings.end())
		{			
			const longdouble rounded_difference_magnitude__it_rank = std::round(it_rank->first);
			
			type_list_longdouble__bool::const_iterator it_next = it_rank;
			++it_next;
			
			uint number_of_equivalents = 1;			
			uint equiv_rank_sum = rank_number;
			while (it_next != rankings.end()  and  std::round(it_next->first) == rounded_difference_magnitude__it_rank)
			{
				equiv_rank_sum += (rank_number + number_of_equivalents);
				++number_of_equivalents;
				++it_next;
			}
			
			
			const longdouble average_rank_among_equals = (longdouble)equiv_rank_sum/number_of_equivalents;
			
			for (type_list_longdouble__bool::const_iterator it_dummy = it_rank;
				it_dummy != it_next;
				++it_dummy)
			{
				if (it_dummy->second)  // positive
					observed_signed_ranksum += average_rank_among_equals;
				else
					observed_signed_ranksum -= average_rank_among_equals;				
			}//dummy
			
			
			//increment:
			it_rank = it_next;
			rank_number += number_of_equivalents;			
		}//rankings
		
	}//observed_signed_ranksum	
	
	
	if (verbose)
	{
		std::stringstream debugstrm;		
		debugstrm << "rankings:\n\t";
		for (type_list_longdouble__bool::const_iterator it_rank = rankings.begin();
			it_rank != rankings.end();
			++it_rank)
		{
			debugstrm << "[" << it_rank->first << ", " << it_rank->second << "]\n\t";									
		}
		std::cerr<< debugstrm.str() << "\n";	
	}	
	
	
	number_of_elts_in_ranking = rankings.size();
	
}//compute_signed_rank_sum













int get_uniform_random_draw
			(const int &max_num)
{	
	boost::random::mt19937 random_gen_for_unif_int;
	boost::random::uniform_int_distribution<> uniform_shuffle_dist(0, max_num-1);		
	
	return uniform_shuffle_dist(random_gen_for_unif_int);	
	
}//get_uniform_random_draw

















void  do_signed_rank_test_for_Read_Depth
				(const std::string &postdatatdir,
				 type_map_uint_to_CC &Conn_Comps)
{
	
	const uint binwidth = 100;	
	
	const type_map_uint_to_list_BI nonuniqe_regions_of_genome(load_Non_unique_regions_of_genome());
	
    print_line_of_markers("==");
    std::cerr << " \n\n\n\n\n\ndo_signed_rank_test_for_Read_Depth:\n\n";	

    std::string calls_fname(postdatatdir);
    calls_fname.append("/all_calls");
    
    type_map_string_to_uint_to_Sampled_diploid_Event_data map_genome_to_event_to_sampled_Event(read_Sampled_Events_from_file(calls_fname));
	
	const std::vector<type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator>  loopable_genome_to_sampled_Events(		
							get_loopable_iterators<type_map_string_to_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event));
	
	
	print_set<std::string>(extract_keys_of_map_and_return_as_set<std::string, type_map_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event),
							"genomes to loop over-MPI", NULL, true);	
	
	
	
	// MPI
	for (int genrank = my_MPI_rank; genrank < loopable_genome_to_sampled_Events.size(); genrank += size_MPI_WORLD)
	{
		std::cerr << "\n\n\nrank " << my_MPI_rank << " looking for next genome...\n\n";
    
//     for (type_map_string_to_uint_to_Sampled_diploid_Event_data::reverse_iterator it_gen = map_genome_to_event_to_sampled_Event.rbegin();
// 	    it_gen != map_genome_to_event_to_sampled_Event.rend();
// 	    ++it_gen)
//     {
		
		const type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator it_gen = loopable_genome_to_sampled_Events.at(genrank);
	
		genome_name = it_gen->first;
		
		std::string full_job_dirname;
		{//full_job_dirname
			const boost::filesystem::directory_iterator it_end;  
			for (boost::filesystem::directory_iterator it_jobdir(postdatatdir);
				it_jobdir != it_end;
				++it_jobdir)
			{		
				if (boost::filesystem::is_directory(it_jobdir->path()))
				{
					if (it_jobdir->path().filename().string().find(it_gen->first) != std::string::npos)
					{
						full_job_dirname = boost::filesystem::complete(it_jobdir->path()).string();
						break;			
					}		    		    
				}//dir				
			}//it_jobdir	    	    	    
		}//full_job_dirname
		
		
		std::cerr << "\n\t\tfull_job_dirname = [" << full_job_dirname << "]\n";
		
		output_dir = full_job_dirname;
	
		
		std::string bam_filename;	
		{//full_job_dirname
			const boost::filesystem::directory_iterator it_end;  
			for (boost::filesystem::directory_iterator it_bamdir("/gpfs/scratch/mmparks/SAMs_and_BAMs");
				it_bamdir != it_end;
				++it_bamdir)
			{		
				if (!boost::filesystem::is_directory(it_bamdir->path()))
				{
					if (it_bamdir->path().filename().string().find(it_gen->first) != std::string::npos
					and  it_bamdir->path().filename().string().find(".bas") == std::string::npos
					and  it_bamdir->path().filename().string().find(".bai") == std::string::npos
					and  it_bamdir->path().filename().string().find(".bam") == it_bamdir->path().filename().string().size() - 4)			
					{
						bam_filename = boost::filesystem::complete(it_bamdir->path()).string();
						break;			
					}		    		    
				}//file
			}//it_jobdir	    	    	    
		}//full_job_dirname	
		
		std::cerr << "\t\tbam_filename = [" << bam_filename << "]\n\n";
		
			
		
		assert(!full_job_dirname.empty());
		assert(!bam_filename.empty());			
		
		
		{
			char common_fname[option_size];
			std::sprintf(common_fname, "%s/Wilcoxon_signed_rank_sum_test.ratios.matched", output_dir.c_str());
			boost::filesystem::remove(common_fname);			
		}
		
		
		BamTools::BamReader my_BAM_reader;

		prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population
						(bam_filename,
						std::string(),//gender_and_population_filename,
						my_BAM_reader);		


		read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda(bam_filename, false, NULL);
		std::cerr.flush();
		std::cerr << "done.";
		

	
		
		type_map_uint_to_list_uint sampled_control_bins_by_GC_count;
		
		print_map_to_list<uint,uint>(sampled_control_bins_by_GC_count, "sampled_control_bins_by_GC_count", NULL, true);
		
		{//check sampled controls
			char common_fname[option_size];
			std::sprintf(common_fname, "%s/sampled_controls.%s.%u",
										output_dir.c_str(), genome_name.c_str(), binwidth);
			
			if (boost::filesystem::exists(common_fname))
			{
				std::cerr << "\n\tuploading already-saved sampled controls for genome " << it_gen->first << "...\n";
				std::ifstream infs(common_fname);
				boost::archive::text_iarchive inarch(infs);
				
				inarch >> sampled_control_bins_by_GC_count;
				
				infs.close();				
			}
			
			if (sampled_control_bins_by_GC_count.empty())
			{
				// if stil empty, perform sampling
				std::cerr << "\nget_random_sample_of_bins_from_across_the_genome for genome " << it_gen->first << "...\n";
				sampled_control_bins_by_GC_count = 
									get_random_sample_of_bins_from_across_the_genome
											(binwidth,
											5000, //number_of_samples,
											my_BAM_reader);
				std::cerr  << " DONE!\n";
				
				{//save
					std::ofstream outfs(common_fname);
					boost::archive::text_oarchive outarch(outfs);
					
					outarch << sampled_control_bins_by_GC_count;

					outfs.close();			
				}//save
			}//sample
			
		}//check sampled controls
		

		
		
		
		const std::vector<type_map_uint_to_Sampled_diploid_Event_data::iterator>  myev__loop(
						get_loopable_iterators<type_map_uint_to_Sampled_diploid_Event_data>(it_gen->second));
	
		
		#pragma omp parallel
		{
			#pragma omp for schedule(dynamic,10)
			for (int myev_ctr = 0; myev_ctr < myev__loop.size(); ++myev_ctr)
			{    
				std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() 
						<< "  myev_ctr = " << myev_ctr << " out of " << myev__loop.size() << ".\n";
				
				const type_map_uint_to_Sampled_diploid_Event_data::iterator it_myev = myev__loop.at(myev_ctr);
						
				it_myev->second.print_this_Sampled_diploid_Event_data();
				
				const type_map_uint_to_Event::iterator result_Ev_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID); 
				
				
				if (result_Ev_it->second.chromos[0] == 24  and   gender_of_individual_is_female)
				{
					std::cerr << "Event is on Y chromosome, but individual is female.  Skipping read-depth graph in in \"create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights\".\n";
					continue;
				}
				
				const BOOST_Interval expanded_region_of_interest( 
							safe_subtract_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.lower(), 100000),
							safe_chromo_add_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.upper(), 100000, result_Ev_it->second.chromos[0]));        
							
				//identify all events intersecting this region:
				type_map_uint_to_Event_ptr  all_events_near_expanded_region;
				for (type_map_uint_to_Event::iterator it_ev = Conn_Comps.at(it_myev->second.cc_of_event).events.begin();
					it_ev != Conn_Comps.at(it_myev->second.cc_of_event).events.end();
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
				
				
				
				
				
				
				char out_RD_core_fname[option_size];
				
				std::sprintf(out_RD_core_fname, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s_%s___%s", 
							output_dir.c_str(),
							it_myev->second.cc_of_event, result_Ev_it->second.UID,
							convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first).c_str(),
							convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second).c_str(),
							genome_name.c_str());		    
						
				
		// 	    //if already computed, skip)
		// 	    {//check already computed
		// 		std::string testname(out_RD_core_fname);
		// 		testname.append(".tests");
		// 		
		// 		if (boost::filesystem::exists(testname) != 0)
		// 		{
		// 		    std::ifstream infs(testname);
		// 		    
		// 		    std::string inteststr(option_size*10, '\0');
		// 		    char intestfile[option_size*10];
		// 		    std::sprintf(intestfile, "%s", inteststr.c_str());
		// 		    		    
		// 		    infs.read(intestfile, option_size*10);
		// 		    infs.close();
		// 		    
		// 		    inteststr.assign(intestfile);
		// 		    const uint num_newlines = std::count(inteststr.begin(), inteststr.end(), '\n');		    
		// 		    if (num_newlines >= 31)
		// 		    {
		// 			std::cerr << "\n\tnum_newlines = " << num_newlines << " !!!  continue!\n";
		// 			continue;			
		// 		    }
		// 		}//does not exit
		// 	    }//check already computed
				
				
				
				
				std::string observed_fname(out_RD_core_fname);
				observed_fname.append(".data.observed");	    
				
				if (boost::filesystem::exists(observed_fname) == 0)
				{//search
					std::string searchdir(output_dir);
					searchdir.append("/Read_Depth");
									
					std::stringstream cc_strm;
					cc_strm << "cc_" << it_myev->second.cc_of_event << "_";
					std::stringstream ev_strm;
					ev_strm << "ev_" << it_myev->second.event_UID << "_";
					
					char observed_fname_old[option_size];
					
					bool found_something = false;
					const boost::filesystem::directory_iterator it_end;  
					for (boost::filesystem::directory_iterator it_rdfile(searchdir);
						it_rdfile != it_end;
						++it_rdfile)
					{		
						if (!boost::filesystem::is_directory(it_rdfile->path()))
						{
							if (it_rdfile->path().filename().string().find(cc_strm.str()) != std::string::npos
								and  it_rdfile->path().filename().string().find(ev_strm.str()) != std::string::npos
								and  it_rdfile->path().filename().string().find(genome_name) != std::string::npos
								and it_rdfile->path().filename().string().find("RD_hypothesis_test__cc_") != std::string::npos)
							{
								char scanformat[option_size];
								std::sprintf(scanformat, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%%s", 
											output_dir.c_str(),
											it_myev->second.cc_of_event,
											it_myev->second.event_UID);	
								
								char in_calls_and_genome[option_size];		    
												
				// 			    std::cerr << "\n\t\t\t\ttry name: [" << boost::filesystem::complete(it_rdfile->path()).string().c_str() << "\n\t\t\t\tscanformat = [ " << scanformat << "]\n";
								std::sscanf(boost::filesystem::complete(it_rdfile->path()).string().c_str(), scanformat,
										in_calls_and_genome);
								
								std::string calls_CHOPPED(in_calls_and_genome);
								calls_CHOPPED = calls_CHOPPED.substr(0,calls_CHOPPED.find(genome_name));
								
				// 			    std::cerr << "\n\t\t\t\tcalls_CHOPPED = [" << calls_CHOPPED << "]\n";
								
								std::sprintf(observed_fname_old, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s%s.data.observed", 
									output_dir.c_str(),
									it_myev->second.cc_of_event, result_Ev_it->second.UID,
									calls_CHOPPED.c_str(),
									genome_name.c_str());	    
								
								found_something = true;
								break;			
							}		    		    
						}//file
					}//it_jobdir
					
					
					if (!found_something)
					{
						std::stringstream error_strm;
						error_strm << "ERROR  unable to find anything!!!\n\n";
						error_message(error_strm,false);
						continue;
					}	    
					
					std::cerr << "\n\t\t\tout_RD_core_fname = [" << out_RD_core_fname << "]\n"
						<<  "\t\t\tobserved_fname_old = [" << observed_fname_old << "\n\n";
											
					if (boost::filesystem::exists(observed_fname) == 0)
					{  
						std::cerr << "rename = [" << observed_fname << "\n";
						boost::filesystem::rename(observed_fname_old, observed_fname);		    
					}
					
				}//search

				
				const type_map_uint_to_longdouble map_Ref_Genome_position_to_per_read_posterior_sum(read_map_from_file<uint,longdouble>(observed_fname, false));
			
				
				
				{//RD hypothesis test
				
					BOOST_Interval affected_region(empty_BI);
					
					if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None))
					{
						affected_region = result_Ev_it->second.region_between_and_including_the_LCRs_themselves;	    
					}
					else
					{
						for (uint hap=0; hap<2; ++hap)
						{
							if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) != hap_outcome__None)
							{
								const type_BI__BI affected_regions_each_LCR(
											BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(
												convert_each_profile_breakpoint_to_absolute_breakpoints(
													pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
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
							

					
// 					std::cerr << "\ncalculating Wilcoxon_pvalue__inferred...\n";
					
					



					//null:
					Wilcoxon_calculations observed_wilcoxon_stats___null___inferred;
					Wilcoxon_calculations observed_wilcoxon_stats___null___unique;
					
					perform_Wilcoxon_signed_rank_test
								(binwidth,
								result_Ev_it->second.chromos[0],
								map_affected_Ref_Genome_position_to_per_read_posterior_sum,
								it_myev->second,
								sampled_control_bins_by_GC_count,
								1,
								nonuniqe_regions_of_genome,
								observed_wilcoxon_stats___null___inferred,
								observed_wilcoxon_stats___null___unique);
								
								
					//alternative:
					Wilcoxon_calculations observed_wilcoxon_stats___alternative___inferred;
					Wilcoxon_calculations observed_wilcoxon_stats___alternative___unique;
					
					perform_Wilcoxon_signed_rank_test
								(binwidth,
								result_Ev_it->second.chromos[0],
								map_affected_Ref_Genome_position_to_per_read_posterior_sum,
								it_myev->second,
								sampled_control_bins_by_GC_count,
								0.5,
								nonuniqe_regions_of_genome,
								observed_wilcoxon_stats___alternative___inferred,
								observed_wilcoxon_stats___alternative___unique);								


								
					
					std::cerr.flush();
					std::cerr << " DONE!\n";
					

							

// 					std::cerr << "wilcoxon_pvalue___inferred = " << wilcoxon_pvalue___inferred << "\n"
// 								<< "wilcoxon_pvalue___unique = " << wilcoxon_pvalue___unique << "\n";
					

				
					#pragma omp critical(write_wilcoxon_signed_rank_to_common_file)
					{//common file		
						//write to common filename
						char common_fname[option_size];
						std::sprintf(common_fname, "%s/Wilcoxon_signed_rank_sum_test.ratios.matched", output_dir.c_str());
						
						std::ofstream out_wilcoxon_fs(common_fname, std::ios_base::app);
						
						boost::math::cauchy_distribution<mpfr_class> my_standard_Cauchy_dist(0,1);
						
						out_wilcoxon_fs 
							//info:
							<< genome_name << "\t"
							<< it_myev->second.cc_of_event << "\t"
							<< it_myev->second.event_UID << "\t"
							
							//p-values:
							<< observed_wilcoxon_stats___null___inferred.wilcoxon_pvalue << "\t"
							<< observed_wilcoxon_stats___null___unique.wilcoxon_pvalue << "\t"
							<< observed_wilcoxon_stats___alternative___inferred.wilcoxon_pvalue << "\t"
							<< observed_wilcoxon_stats___alternative___unique.wilcoxon_pvalue << "\t"
							
							//likelihoods:
							<< observed_wilcoxon_stats___null___inferred.wilcoxon_likelihood << "\t"
							<< observed_wilcoxon_stats___null___unique.wilcoxon_likelihood << "\t"
							<< observed_wilcoxon_stats___alternative___inferred.wilcoxon_likelihood << "\t"
							<< observed_wilcoxon_stats___alternative___unique.wilcoxon_likelihood << "\t"
							
							//signed_rank_sums:
							<< observed_wilcoxon_stats___null___inferred.signed_rank_sum_value << "\t"
							<< observed_wilcoxon_stats___null___unique.signed_rank_sum_value << "\t"
							<< observed_wilcoxon_stats___alternative___inferred.signed_rank_sum_value << "\t"
							<< observed_wilcoxon_stats___alternative___unique.signed_rank_sum_value << "\t"
							
							//Likelihood ratios:
							<< observed_wilcoxon_stats___null___inferred.wilcoxon_likelihood / observed_wilcoxon_stats___alternative___inferred.wilcoxon_likelihood << "\t"
							<< observed_wilcoxon_stats___null___unique.wilcoxon_likelihood / observed_wilcoxon_stats___alternative___unique.wilcoxon_likelihood << "\t"							
							
							//ratios:
							<< observed_wilcoxon_stats___null___inferred.signed_rank_sum_value /
								observed_wilcoxon_stats___alternative___inferred.signed_rank_sum_value << "\t"
							<< observed_wilcoxon_stats___null___unique.signed_rank_sum_value /
								observed_wilcoxon_stats___alternative___unique.signed_rank_sum_value << "\t"
								
							//Cauchy CDF (standard Cauchy)
							<< real(boost::math::cdf<mpfr_class>(my_standard_Cauchy_dist, mpfr_class((double)observed_wilcoxon_stats___null___inferred.signed_rank_sum_value) /
								mpfr_class((double)observed_wilcoxon_stats___alternative___inferred.signed_rank_sum_value)).get_mpfr_t()) << "\t"
							<< real(boost::math::cdf<mpfr_class>(my_standard_Cauchy_dist, mpfr_class((double)observed_wilcoxon_stats___null___unique.signed_rank_sum_value) /
								mpfr_class((double)observed_wilcoxon_stats___alternative___unique.signed_rank_sum_value)).get_mpfr_t()) << "\n";
								
							
							
						out_wilcoxon_fs.close();
					}//common file
					
					std::cerr << "\n\twriting_is_DONE.\n";
				
				}//RD hypothesis test
			
			}//it_myev
			
			std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "  done with it_myev loop...\n";
			
		}//parallel
		std::cerr << "\nrank " << my_MPI_rank << " moving to next genome...\n";
		
    }//it_gen            		
	
}//do_signed_rank_test_for_Read_Depth






















void  bernoulli_lyapunov_nonparametric_test
				(const std::string &postdatatdir,
				 type_map_uint_to_CC &Conn_Comps)
{
	
// 	//debug
// 	boost::math::normal_distribution<mpfr_class> my_debug_normal((double)914.169, (double)std::sqrt((longdouble)114677));
// 	std::cerr << "\n\ndebug-normal:  boost::math::cdf<longdouble>(my_debug_normal, 1235.91) = " 
// 			<< real(boost::math::cdf<mpfr_class>(my_debug_normal, (double)1235.91).get_mpfr_t()) <<  "\n\n";
// 			
// 	exit(1);
// 	//debug
	
	
	
	const uint binwidth = 100;	
	
	const type_map_uint_to_list_BI nonuniqe_regions_of_genome(load_Non_unique_regions_of_genome());
	
	print_line_of_markers("==");
	std::cerr << " \n\n\n\n\n\ndo_signed_rank_test_for_Read_Depth:\n\n";	

	std::string calls_fname(postdatatdir);
	calls_fname.append("/all_calls");
	
	type_map_string_to_uint_to_Sampled_diploid_Event_data map_genome_to_event_to_sampled_Event(read_Sampled_Events_from_file(calls_fname));
	
	const std::vector<type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator>  loopable_genome_to_sampled_Events(		
							get_loopable_iterators<type_map_string_to_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event));
	
	
	print_set<std::string>(extract_keys_of_map_and_return_as_set<std::string, type_map_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event),
							"genomes to loop over-MPI", NULL, true);	
	
	
	
	// MPI
	for (int genrank = my_MPI_rank; genrank < loopable_genome_to_sampled_Events.size(); genrank += size_MPI_WORLD)
	{
		std::cerr << "\n\n\nrank " << my_MPI_rank << " looking for next genome...\n\n";

		
		const type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator it_gen = loopable_genome_to_sampled_Events.at(genrank);
	
		genome_name = it_gen->first;
		
		std::string full_job_dirname;
		{//full_job_dirname
			const boost::filesystem::directory_iterator it_end;  
			for (boost::filesystem::directory_iterator it_jobdir(postdatatdir);
				it_jobdir != it_end;
				++it_jobdir)
			{		
				if (boost::filesystem::is_directory(it_jobdir->path()))
				{
					if (it_jobdir->path().filename().string().find(it_gen->first) != std::string::npos)
					{
						full_job_dirname = boost::filesystem::complete(it_jobdir->path()).string();
						break;			
					}		    		    
				}//dir				
			}//it_jobdir	    	    	    
		}//full_job_dirname
		
		
		std::cerr << "\n\n\nrank " << my_MPI_rank << "\tfull_job_dirname = [" << full_job_dirname << "]\n";
		
		output_dir = full_job_dirname;
	
		
		std::string bam_filename;	
		{//full_job_dirname
			const boost::filesystem::directory_iterator it_end;  
			for (boost::filesystem::directory_iterator it_bamdir("/gpfs/scratch/mmparks/SAMs_and_BAMs");
				it_bamdir != it_end;
				++it_bamdir)
			{		
				if (!boost::filesystem::is_directory(it_bamdir->path()))
				{
					if (it_bamdir->path().filename().string().find(it_gen->first) != std::string::npos
					and  it_bamdir->path().filename().string().find(".bas") == std::string::npos
					and  it_bamdir->path().filename().string().find(".bai") == std::string::npos
					and  it_bamdir->path().filename().string().find(".bam") == it_bamdir->path().filename().string().size() - 4)			
					{
						bam_filename = boost::filesystem::complete(it_bamdir->path()).string();
						break;			
					}		    		    
				}//file
			}//it_jobdir	    	    	    
		}//full_job_dirname	
		
		std::cerr << "\n\n\nrank " << my_MPI_rank << "\ttbam_filename = [" << bam_filename << "]\n\n";
		
			
		
		assert(!full_job_dirname.empty());
		assert(!bam_filename.empty());			
		
		
// 		{
// 			char common_fname[option_size];
// 			std::sprintf(common_fname, "%s/Wilcoxon_signed_rank_sum_test.lyap.normalized", output_dir.c_str());
// 			boost::filesystem::remove(common_fname);			
// 		}
		
		
		BamTools::BamReader my_BAM_reader;

		prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population
						(bam_filename,
						std::string(),//gender_and_population_filename,
						my_BAM_reader);		


		read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda(bam_filename, false, NULL);
		std::cerr.flush();
		std::cerr << "\n\n\nrank " << my_MPI_rank << "\tprepare-done.";
		

	
		
		type_map_uint_to_list_uint sampled_control_bins_by_GC_count;
		
		print_map_to_list<uint,uint>(sampled_control_bins_by_GC_count, "sampled_control_bins_by_GC_count", NULL, true);
		
		{//check sampled controls
			char common_fname[option_size];
			std::sprintf(common_fname, "%s/sampled_controls.%s.%u",
										output_dir.c_str(), genome_name.c_str(), binwidth);
			
			if (boost::filesystem::exists(common_fname))
			{
				std::cerr << "\n\n\nrank " << my_MPI_rank << "\tuploading already-saved sampled controls for genome " << it_gen->first << "...\n";
				std::ifstream infs(common_fname);
				boost::archive::text_iarchive inarch(infs);
				
				inarch >> sampled_control_bins_by_GC_count;
				
				infs.close();				
			}
			
			if (sampled_control_bins_by_GC_count.empty())
			{
				// if stil empty, perform sampling
				std::cerr<< "\n\n\nrank " << my_MPI_rank << "\tget_random_sample_of_bins_from_across_the_genome for genome " << it_gen->first << "...\n";
				sampled_control_bins_by_GC_count = 
									get_random_sample_of_bins_from_across_the_genome
											(binwidth,
											5000, //number_of_samples,
											my_BAM_reader);
				std::cerr  << " DONE!\n";
				
				{//save
					std::ofstream outfs(common_fname);
					boost::archive::text_oarchive outarch(outfs);
					
					outarch << sampled_control_bins_by_GC_count;

					outfs.close();			
				}//save
			}//sample
			
		}//check sampled controls
		

		
		
		
		const std::vector<type_map_uint_to_Sampled_diploid_Event_data::iterator>  myev__loop(
						get_loopable_iterators<type_map_uint_to_Sampled_diploid_Event_data>(it_gen->second));
	
		
		
		#pragma omp parallel
		{
			uint thread_work_ctr = 0;  //for debug
			
			#pragma omp for schedule(static)
			for (int myev_ctr = 0; myev_ctr < myev__loop.size(); ++myev_ctr)
			{    				
				std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() 
						<< "  myev_ctr = " << myev_ctr << " out of " << myev__loop.size() << ".\n";
						
						
				std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "-work-counter = " << thread_work_ctr << "\n";
				++thread_work_ctr;
				
				const type_map_uint_to_Sampled_diploid_Event_data::iterator it_myev = myev__loop.at(myev_ctr);
												
				
// // 				//DEBUG!!!
// 				if (it_myev->second.event_UID != 17095 or it_myev->second.associated_genome.compare("NA12716") != 0)
// 					continue;
// // 				//DEBUG!!
						
							
				const type_map_uint_to_Event::iterator result_Ev_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID); 
				
				
				if (result_Ev_it->second.chromos[0] == 24  and   gender_of_individual_is_female)
				{
					std::cerr << "Event is on Y chromosome, but individual is female.  Skipping read-depth graph in in \"create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights\".\n";
					continue;
				}
				
				const BOOST_Interval expanded_region_of_interest( 
							safe_subtract_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.lower(), 100000),
							safe_chromo_add_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.upper(), 100000, result_Ev_it->second.chromos[0]));        
							
				//identify all events intersecting this region:
				type_map_uint_to_Event_ptr  all_events_near_expanded_region;
				for (type_map_uint_to_Event::iterator it_ev = Conn_Comps.at(it_myev->second.cc_of_event).events.begin();
					it_ev != Conn_Comps.at(it_myev->second.cc_of_event).events.end();
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
				
				
				
				
				
				
				char out_RD_core_fname[option_size];
				
				std::sprintf(out_RD_core_fname, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s_%s___%s", 
							output_dir.c_str(),
							it_myev->second.cc_of_event, result_Ev_it->second.UID,
							convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first).c_str(),
							convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second).c_str(),
							genome_name.c_str());		    
						
				
		// 	    //if already computed, skip)
		// 	    {//check already computed
		// 		std::string testname(out_RD_core_fname);
		// 		testname.append(".tests");
		// 		
		// 		if (boost::filesystem::exists(testname) != 0)
		// 		{
		// 		    std::ifstream infs(testname);
		// 		    
		// 		    std::string inteststr(option_size*10, '\0');
		// 		    char intestfile[option_size*10];
		// 		    std::sprintf(intestfile, "%s", inteststr.c_str());
		// 		    		    
		// 		    infs.read(intestfile, option_size*10);
		// 		    infs.close();
		// 		    
		// 		    inteststr.assign(intestfile);
		// 		    const uint num_newlines = std::count(inteststr.begin(), inteststr.end(), '\n');		    
		// 		    if (num_newlines >= 31)
		// 		    {
		// 			std::cerr << "\n\tnum_newlines = " << num_newlines << " !!!  continue!\n";
		// 			continue;			
		// 		    }
		// 		}//does not exit
		// 	    }//check already computed
				
				
				
				
				std::string observed_fname(out_RD_core_fname);
				observed_fname.append(".data.observed");	    
				
				if (boost::filesystem::exists(observed_fname) == 0)
				{//search
					std::string searchdir(output_dir);
					searchdir.append("/Read_Depth");
									
					std::stringstream cc_strm;
					cc_strm << "cc_" << it_myev->second.cc_of_event << "_";
					std::stringstream ev_strm;
					ev_strm << "ev_" << it_myev->second.event_UID << "_";
					
					char observed_fname_old[option_size];
					
					bool found_something = false;
					const boost::filesystem::directory_iterator it_end;  
					for (boost::filesystem::directory_iterator it_rdfile(searchdir);
						it_rdfile != it_end;
						++it_rdfile)
					{		
						if (!boost::filesystem::is_directory(it_rdfile->path()))
						{
							if (it_rdfile->path().filename().string().find(cc_strm.str()) != std::string::npos
								and  it_rdfile->path().filename().string().find(ev_strm.str()) != std::string::npos
								and  it_rdfile->path().filename().string().find(genome_name) != std::string::npos
								and it_rdfile->path().filename().string().find("RD_hypothesis_test__cc_") != std::string::npos)
							{
								char scanformat[option_size];
								std::sprintf(scanformat, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%%s", 
											output_dir.c_str(),
											it_myev->second.cc_of_event,
											it_myev->second.event_UID);	
								
								char in_calls_and_genome[option_size];		    
												
				// 			    std::cerr << "\n\t\t\t\ttry name: [" << boost::filesystem::complete(it_rdfile->path()).string().c_str() << "\n\t\t\t\tscanformat = [ " << scanformat << "]\n";
								std::sscanf(boost::filesystem::complete(it_rdfile->path()).string().c_str(), scanformat,
										in_calls_and_genome);
								
								std::string calls_CHOPPED(in_calls_and_genome);
								calls_CHOPPED = calls_CHOPPED.substr(0,calls_CHOPPED.find(genome_name));
								
				// 			    std::cerr << "\n\t\t\t\tcalls_CHOPPED = [" << calls_CHOPPED << "]\n";
								
								std::sprintf(observed_fname_old, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s%s.data.observed", 
									output_dir.c_str(),
									it_myev->second.cc_of_event, result_Ev_it->second.UID,
									calls_CHOPPED.c_str(),
									genome_name.c_str());	    
								
								found_something = true;
								break;			
							}		    		    
						}//file
					}//it_jobdir
					
					
					if (!found_something)
					{
						std::stringstream error_strm;
						error_strm << "ERROR  unable to find anything!!!\n\n";
						error_message(error_strm,false);
						continue;
					}	    
					
					std::cerr << "\n\t\t\tout_RD_core_fname = [" << out_RD_core_fname << "]\n"
						<<  "\t\t\tobserved_fname_old = [" << observed_fname_old << "\n\n";
											
					if (boost::filesystem::exists(observed_fname) == 0)
					{  
						std::cerr << "rename = [" << observed_fname << "\n";
						boost::filesystem::rename(observed_fname_old, observed_fname);		    
					}
					
				}//search

				
				const type_map_uint_to_longdouble map_Ref_Genome_position_to_per_read_posterior_sum(read_map_from_file<uint,longdouble>(observed_fname, false));
			
				
				
				{//RD hypothesis test
				
					BOOST_Interval affected_region(empty_BI);
					
					if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None))
					{
						affected_region = result_Ev_it->second.region_between_and_including_the_LCRs_themselves;	    
					}
					else
					{
						for (uint hap=0; hap<2; ++hap)
						{
							if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) != hap_outcome__None)
							{
								const type_BI__BI affected_regions_each_LCR(
											BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(
												convert_each_profile_breakpoint_to_absolute_breakpoints(
													pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
													result_Ev_it->second.compressed_map_LCR_to_profile__for_each_LCR)));
											
								const BOOST_Interval affected_hull(affected_regions_each_LCR.first.lower(), affected_regions_each_LCR.second.upper());
																							
								if (BOOST_empty(affected_region))
								{  affected_region = affected_hull;  }
								else
								{  affected_region = BOOST_hull(affected_region, affected_hull);  }
							}//occurs
						}//hap
					}
					
					
					
					std::stringstream  out_ss;
				
					out_ss  << "\nrank " << my_MPI_rank << ", thread " << omp_get_thread_num() 
							<< "\n\tmyev_ctr = " << myev_ctr
							<< "\n\tgenome = " << genome_name
							<< "\n\tcc = " << it_myev->second.cc_of_event
							<< "\n\tevent = " << it_myev->second.event_UID
							<< "\n\taffected_region = " << affected_region 
							<< "\n\tchrome = " <<  result_Ev_it->second.chromos[0]
							<< "\n\tchrom-length = " << Event::chromosome_lengths.at(result_Ev_it->second.chromos[0]) 
							<< "\n\toutcome = " << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first) 
							<< ", " << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second) 
							<< "\n";


					//restrict to affected region:
					const type_map_uint_to_longdouble  map_affected_Ref_Genome_position_to_per_read_posterior_sum(
							get_map_intersection_of_map_keys_and_interval<uint,longdouble>(map_Ref_Genome_position_to_per_read_posterior_sum, affected_region));
							

// 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get affected_sequence__padded\n";
// 					std::cerr << "\ncalculating Wilcoxon_pvalue__inferred...\n";
					const uint pad_amount = 5000;
					const uint padded_max_chromo_sensitive = safe_chromo_add_base_1(affected_region.upper(), pad_amount, result_Ev_it->second.chromos[0]);
					//std::max<uint>(affected_region.upper() + pad_amount, );
					
					const std::string affected_sequence__padded(get_specific_DNA_sequence_from_Reference_genome(result_Ev_it->second.chromos[0],
														BOOST_Interval(affected_region.lower(), padded_max_chromo_sensitive)));
					
					
// 					type_vector_longdouble bernoulli_p_by_position(BOOST_width_inclusive(affected_region), 0);
					type_vector_vector_longdouble poisson_lambda_by_position_and_readgroup___NULL(BOOST_width_inclusive(affected_region),
																    type_vector_longdouble(Readgroup_stats.size(), 0));		
					
					type_vector_uint bernoulli_observations_by_position(BOOST_width_inclusive(affected_region), 0);
					
					boost::math::normal_distribution<mpfr_class> my_standard_normal(0,1);
					
// 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get bernoulli_p_by_position\n";
					//get bernoulli p at each position by summing over all readgroups, using GC content at that position.
					for (uint le = 0; le < BOOST_width_inclusive(affected_region); ++le)
					{
						uint rg_indx = 0;
						for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
							it_rg != Readgroup_stats.end();
							++it_rg)
						{
							const uint gc_content = calculate_GC_content_of_sequence(affected_sequence__padded.substr(le, it_rg->second.median_insert_size));
							
							poisson_lambda_by_position_and_readgroup___NULL.at(le).at(rg_indx)
									= it_rg->second.haploid_fragmentation_rate__per_GC_count.at(gc_content) * 2; // null case is diploid
							
							++rg_indx;
						}
					}//bernoulli p
					
					
					
					type_vector_longdouble bernoulli_p_0_by_position__null_model(BOOST_width_inclusive(affected_region), 0);
					for (uint i=0; i < bernoulli_p_0_by_position__null_model.size(); ++i)
					{
						longdouble total_lambda_over_readgroups = 0;
						for (uint rg=0; rg < poisson_lambda_by_position_and_readgroup___NULL.at(i).size(); ++rg)
						{  total_lambda_over_readgroups += poisson_lambda_by_position_and_readgroup___NULL.at(i).at(rg);  }
						
						bernoulli_p_0_by_position__null_model.at(i) = std::exp(-total_lambda_over_readgroups);
					}//i
					
					
					
					
					
// 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get bernoulli_observations_by_position\n";
					//determine bernoulli observation at each point.
					uint num_positives = 0;
					for (type_map_uint_to_longdouble::const_iterator it_obs = map_affected_Ref_Genome_position_to_per_read_posterior_sum.begin();
						it_obs != map_affected_Ref_Genome_position_to_per_read_posterior_sum.end();
						++it_obs)
					{
						if (it_obs->second > 0.5)
						{
							bernoulli_observations_by_position.at(it_obs->first - affected_region.lower()) = 1;
							++num_positives;
						}
					}//observed
					
					
					
					out_ss << "\n\tnum_obs/num_pos = " << (longdouble)num_positives / map_affected_Ref_Genome_position_to_per_read_posterior_sum.size() << "\n\n";
					
					
					const Lyap_CLT_data stats_null_inferred(calculate_Lyapunov_CLT_sum(bernoulli_p_0_by_position__null_model, 
													   poisson_lambda_by_position_and_readgroup___NULL,
													bernoulli_observations_by_position));
					
					
// 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get poisson_lambda_by_position___alternative\n";
					type_vector_vector_longdouble poisson_lambda_by_position___alternative(poisson_lambda_by_position_and_readgroup___NULL);
					
					if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
					{
						out_ss << "\talternative = diploid dup\n";
						for (uint le=0; le < poisson_lambda_by_position___alternative.size(); ++le)
						{  
							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position___alternative.at(le).size(); ++rg_indx)
							{  poisson_lambda_by_position___alternative.at(le).at(rg_indx) *= 2;  }
						}													
					}//diploid dup
					else if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
					{
						out_ss << "\talternative = diploid del\n";
						for (uint le=0; le < poisson_lambda_by_position___alternative.size(); ++le)
						{  
							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position___alternative.at(le).size(); ++rg_indx)
							{  poisson_lambda_by_position___alternative.at(le).at(rg_indx) /= 10;  }
						}													
					}//diploid del		
					else if (test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
					{
						out_ss << "\talternative = haploid dup\n";
						for (uint le=0; le < poisson_lambda_by_position___alternative.size(); ++le)
						{  
							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position___alternative.at(le).size(); ++rg_indx)
							{  poisson_lambda_by_position___alternative.at(le).at(rg_indx) *= 1.5;  }
						}
					}//haploid dup
					else if (test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
					{
						out_ss << "\talternative = haploid del\n";
						for (uint le=0; le < poisson_lambda_by_position___alternative.size(); ++le)
						{  
							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position___alternative.at(le).size(); ++rg_indx)
							{  poisson_lambda_by_position___alternative.at(le).at(rg_indx) /= 2;  }
						}
					}//haploid del
					
					
					
					
					
					Lyap_CLT_data stats_alt_inferred;
					bool best_alt_is_del = false;
					bool modified_diploid_call_to_haploid_same_call = false;
					
					real Lyap_CLT_tail_alt___inferred(-1);
					real Lyap_CLT_normalized_stat_alt___inferred(0);
					
					
					if (test_if_at_least_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome))
					{		
						stats_alt_inferred = calculate_Lyapunov_CLT_sum(bernoulli_p_0_by_position__null_model,
												poisson_lambda_by_position___alternative, 
												bernoulli_observations_by_position);
						
						Lyap_CLT_normalized_stat_alt___inferred = compute_normalized_statistic(stats_null_inferred, stats_alt_inferred);
						Lyap_CLT_tail_alt___inferred = compute_lower_tail_probability(stats_null_inferred, stats_alt_inferred);
						
						
						out_ss
							<< "\tgiven-alternative:    "
							<< Lyap_CLT_tail_alt___inferred 
							<< "  =  N(  " << stats_null_inferred.sum_log_P_obs << " - " << stats_alt_inferred.sum_log_P_obs
								<< " = " << (stats_null_inferred.sum_log_P_obs - stats_alt_inferred.sum_log_P_obs) 
							<< "  ;  " << stats_null_inferred.sum_mu << " - " << stats_alt_inferred.sum_mu	
								<< " = " << (stats_null_inferred.sum_mu - stats_alt_inferred.sum_mu)
							<< "  ,  " <<  stats_null_inferred.sn_2 << " + " << stats_alt_inferred.sn_2 
								<< " = " << (stats_null_inferred.sn_2 + stats_alt_inferred.sn_2)
							<< "  )\n"
							<< "\t\tnormalized-stat = " << Lyap_CLT_normalized_stat_alt___inferred <<  "   -->  " 
							<< real(boost::math::cdf<mpfr_class>(my_standard_normal, Lyap_CLT_normalized_stat_alt___inferred.toDouble()).get_mpfr_t())
							<< "\n";
					}
					
					
					
					if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome)  
						or  !test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome))
					{//alternatives
					
						type_vector_vector_longdouble poisson_lambda_by_position__alternative__hap_del(poisson_lambda_by_position_and_readgroup___NULL);
						type_vector_vector_longdouble poisson_lambda_by_position__alternative__hap_dup(poisson_lambda_by_position_and_readgroup___NULL);
						for (uint le=0; le < poisson_lambda_by_position_and_readgroup___NULL.size(); ++le)
						{  
							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position_and_readgroup___NULL.at(le).size(); ++rg_indx)
							{ 
								poisson_lambda_by_position__alternative__hap_del.at(le).at(rg_indx) /= 2; 
								poisson_lambda_by_position__alternative__hap_dup.at(le).at(rg_indx) *= 1.5;
							}
						}
						
						const Lyap_CLT_data stats_alt_del_inferred(
									calculate_Lyapunov_CLT_sum(bernoulli_p_0_by_position__null_model, 
												   poisson_lambda_by_position__alternative__hap_del, 
													bernoulli_observations_by_position));
						const Lyap_CLT_data stats_alt_dup_inferred(
									calculate_Lyapunov_CLT_sum(bernoulli_p_0_by_position__null_model,
												   poisson_lambda_by_position__alternative__hap_dup,
												   bernoulli_observations_by_position));
						
						const real Lyap_CLT_tail___del(compute_lower_tail_probability(stats_null_inferred, stats_alt_del_inferred));
						const real Lyap_CLT_normalized_stat___del(compute_normalized_statistic(stats_null_inferred, stats_alt_del_inferred));

						const real Lyap_CLT_tail___dup(compute_lower_tail_probability(stats_null_inferred, stats_alt_dup_inferred));
						const real Lyap_CLT_normalized_stat___dup(compute_normalized_statistic(stats_null_inferred, stats_alt_dup_inferred));					
						
						if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
						{
							if (Lyap_CLT_normalized_stat___del < Lyap_CLT_normalized_stat_alt___inferred)
							{
								out_ss << "\tLyap_CLT_normalized_stat___del = " << Lyap_CLT_normalized_stat___del << " < "
									<< Lyap_CLT_normalized_stat_alt___inferred << " = Lyap_CLT_normalized_stat_alt___inferred\n";
									
								best_alt_is_del = true;
								stats_alt_inferred = stats_alt_del_inferred;
								Lyap_CLT_tail_alt___inferred = Lyap_CLT_tail___del;
								Lyap_CLT_normalized_stat_alt___inferred = Lyap_CLT_normalized_stat___del;
								modified_diploid_call_to_haploid_same_call = true;
								poisson_lambda_by_position___alternative = poisson_lambda_by_position__alternative__hap_del;
							}
						}
						else if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
						{
							if (Lyap_CLT_normalized_stat___dup < Lyap_CLT_normalized_stat_alt___inferred)
							{
								out_ss << "\tLyap_CLT_normalized_stat___dup = " << Lyap_CLT_normalized_stat___dup << " < "
									<< Lyap_CLT_normalized_stat_alt___inferred << " = Lyap_CLT_normalized_stat_alt___inferred\n";
									
								best_alt_is_del = false;
								stats_alt_inferred = stats_alt_dup_inferred;
								Lyap_CLT_tail_alt___inferred = Lyap_CLT_tail___dup;
								Lyap_CLT_normalized_stat_alt___inferred = Lyap_CLT_normalized_stat___dup;
								modified_diploid_call_to_haploid_same_call = true;
								poisson_lambda_by_position___alternative = poisson_lambda_by_position__alternative__hap_dup;
							}
						}
						else // contains neither dup nor del
						{
							if (Lyap_CLT_normalized_stat___del < Lyap_CLT_normalized_stat___dup)
							{
								out_ss << "\tLyap_CLT_normalized_stat___del = " << Lyap_CLT_normalized_stat___del << " < " 
									<< Lyap_CLT_normalized_stat___dup << " = Lyap_CLT_normalized_stat___dup\n";
								best_alt_is_del = true;
								stats_alt_inferred = stats_alt_del_inferred;
								Lyap_CLT_tail_alt___inferred = Lyap_CLT_tail___del;
								Lyap_CLT_normalized_stat_alt___inferred = Lyap_CLT_normalized_stat___del;
								poisson_lambda_by_position___alternative = poisson_lambda_by_position__alternative__hap_del;
							}
							else
							{
								out_ss << "\tLyap_CLT_normalized_stat___dup = " << Lyap_CLT_normalized_stat___dup << " < " << 
									Lyap_CLT_normalized_stat___del << " = Lyap_CLT_normalized_stat___del\n";
								best_alt_is_del = false;
								stats_alt_inferred = stats_alt_dup_inferred;
								Lyap_CLT_tail_alt___inferred = Lyap_CLT_tail___dup;
								Lyap_CLT_normalized_stat_alt___inferred = Lyap_CLT_normalized_stat___dup;
								poisson_lambda_by_position___alternative = poisson_lambda_by_position__alternative__hap_dup;
							}
						}
																		
						
						out_ss
							<< "\tbest alternative:  "
							<< Lyap_CLT_tail_alt___inferred 
							<< "  =  N(  " << stats_null_inferred.sum_log_P_obs << " - " << stats_alt_inferred.sum_log_P_obs
								<< " = " << (stats_null_inferred.sum_log_P_obs - stats_alt_inferred.sum_log_P_obs) 
							<< "  ;  " << stats_null_inferred.sum_mu << " - " << stats_alt_inferred.sum_mu	
								<< " = " << (stats_null_inferred.sum_mu - stats_alt_inferred.sum_mu)
							<< "  ,  " <<  stats_null_inferred.sn_2 << " + " << stats_alt_inferred.sn_2 
								<< " = " << (stats_null_inferred.sn_2 + stats_alt_inferred.sn_2)							
							<< "  )\nmodified = " << modified_diploid_call_to_haploid_same_call << "\n"
							<< "\t\tnormalized-stat = " << Lyap_CLT_normalized_stat_alt___inferred <<  "   -->  " 
							<< real(boost::math::cdf<mpfr_class>(my_standard_normal, Lyap_CLT_normalized_stat_alt___inferred.toDouble()).get_mpfr_t())
							<< "\n";					
					}//find best alt
					
														

// 					std::cerr << "get is_unique\n";	
					type_vector_bool is_unique(BOOST_width_inclusive(affected_region), true);
					{
						const type_map_uint_to_list_BI::const_iterator it_found_chr = nonuniqe_regions_of_genome.find(result_Ev_it->second.chromos[0]); 
						if (it_found_chr != nonuniqe_regions_of_genome.end())
						{
							for (type_list_BI::const_iterator it_bi = it_found_chr->second.begin();
								it_bi != it_found_chr->second.end();
								++it_bi)
							{
								const BOOST_Interval intersection(BOOST_intersect(*it_bi, affected_region));
								if (!BOOST_empty(intersection))
								{
									for (uint pos = intersection.lower(); pos <= intersection.upper(); ++pos)
									{  is_unique.at(pos - affected_region.lower()) = false;  }						
								}														
							}
						}
					}
					
// 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get number_of_unique_positions\n";	
					uint number_of_unique_positions = 0;
					for (type_vector_bool::const_iterator it_uniq = is_unique.begin(); it_uniq != is_unique.end();  ++it_uniq)
					{
						if (*it_uniq)
						{  ++number_of_unique_positions;  }						
					}
					
					
// 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get Lyap_CLT_statistic___unique\n";	


					Lyap_CLT_data stats_null_unique;
					Lyap_CLT_data stats_alt_unique;
					
					real Lyap_CLT_tail___unique(-1);
					real Lyap_CLT_normalized_stat___unique(-1);

					if (number_of_unique_positions > 999)
					{
						stats_null_unique = calculate_Lyapunov_CLT_sum(bernoulli_p_0_by_position__null_model,
											       poisson_lambda_by_position_and_readgroup___NULL,
											       bernoulli_observations_by_position, is_unique);
						stats_alt_unique  = calculate_Lyapunov_CLT_sum(bernoulli_p_0_by_position__null_model,
											       poisson_lambda_by_position___alternative,
											       bernoulli_observations_by_position, is_unique);
						
						Lyap_CLT_tail___unique = compute_lower_tail_probability(stats_null_unique, stats_alt_unique);
						Lyap_CLT_normalized_stat___unique = compute_normalized_statistic(stats_null_unique, stats_alt_unique);
						
						
						out_ss
							<< "\tunique:    "
							<< Lyap_CLT_tail___unique 
							<< "  =  N(  " << stats_null_unique.sum_log_P_obs << " - " << stats_alt_unique.sum_log_P_obs
								<< " = " << (stats_null_unique.sum_log_P_obs - stats_alt_unique.sum_log_P_obs) 
							<< "  ;  " << stats_null_unique.sum_mu << " - " << stats_alt_unique.sum_mu	
								<< " = " << (stats_null_unique.sum_mu - stats_alt_unique.sum_mu)
							<< "  ,  " <<  stats_null_unique.sn_2 << " + " << stats_alt_unique.sn_2 
								<< " = " << (stats_null_unique.sn_2 + stats_alt_unique.sn_2)
							<< "  )\n";
					}

					
					
						
					std::cerr << out_ss.str();
					
				
					#pragma omp critical(write_wilcoxon_signed_rank_to_common_file)
					{//common file		
						//write to common filename
						char common_fname[option_size];
						std::sprintf(common_fname, "%s/Wilcoxon_signed_rank_sum_test.lyap.normalized", output_dir.c_str());
						
						std::ofstream out_wilcoxon_fs(common_fname, std::ios_base::app);
						
						
						out_wilcoxon_fs 
							//info:
							<< genome_name << "\t"
							<< it_myev->second.cc_of_event << "\t"
							<< it_myev->second.event_UID << "\t"
							
							//p-values:
							<< Lyap_CLT_normalized_stat_alt___inferred << "\t"
							<< Lyap_CLT_normalized_stat___unique << "\n";
							
						out_wilcoxon_fs.close();
					}//common file
					
// 					std::cerr  << "\n\trank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "writing_is_DONE.\n";
				
				}//RD hypothesis test			
			}//it_myev
			
			std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "  done with it_myev loop...\n";
			
		}//parallel
		std::cerr << "\nrank " << my_MPI_rank << " moving to next genome...\n";
		
    }//it_gen            		
	
}//bernoulli_lyapunov_nonparametric_test









// type_vector_longdouble  form_P_0_by_position_from_NULL_poisson_lambdas_and_desired_outcome
// 			(const type_vector_longdouble &poisson_lambda_per_position_NULL,
// 			const Event &some_Event,
// 			const type_haploid_outcome__haploid_outcome &desired_diploid_outcome,
// 			const type_uint__uint &breakpoints)
// {
// 	type_vector_longdouble  bernoulli_P_0__by_position_alternative();
// 	
// 	if (test_if_is_diploid_No_Event(desired_diploid_outcome))
// 	{
// 		
// 	}
// 	
// 	
// 	
// }//form_P_fail_by_position_from_NULL_poisson_lambdas_and_desired_outcome










Lyap_CLT_data calculate_Lyapunov_CLT_sum
			(const type_vector_longdouble &bernoulli_p_0_by_position__null_model,
			 const type_vector_vector_longdouble &poisson_lambda_summed_over_readgroup_per_position,
			const type_vector_uint &bernoulli_observations_by_position,
			const type_vector_bool &is_unique)
{
	
// 	if (bernoulli_p_by_position.empty() or bernoulli_observations_by_position.empty())
// 		return -1;
	assert(poisson_lambda_summed_over_readgroup_per_position.size() == bernoulli_observations_by_position.size());
	
	if (!is_unique.empty())
	{  assert(poisson_lambda_summed_over_readgroup_per_position.size() == is_unique.size());  }
	
	Lyap_CLT_data CLT_sums;
	
	uint number_of_good_positions = 0;
	for (uint le=0; le < poisson_lambda_summed_over_readgroup_per_position.size(); ++le)
	{
		if (is_unique.empty()  or  is_unique.at(le))
		{
			++number_of_good_positions;
			
			//add poisson fragmentation rates:
			longdouble total_poisson_frag_rates_across_RGs = 0;
			for (uint rg_indx = 0; rg_indx < poisson_lambda_summed_over_readgroup_per_position.at(le).size(); ++rg_indx)
			{  total_poisson_frag_rates_across_RGs += poisson_lambda_summed_over_readgroup_per_position.at(le).at(rg_indx);  }
			
			const longdouble P_0_model = std::exp(-total_poisson_frag_rates_across_RGs);
			const longdouble P_1_model = 1.00L - P_0_model;
			
			const longdouble P_0_null = bernoulli_p_0_by_position__null_model.at(le);
			const longdouble P_1_null = 1.00L - P_0_null;
			
			
			const longdouble mu_A = P_0_null * std::log(P_0_model) + P_1_null * std::log(P_1_model);
			const longdouble var_A =  P_0_null * std::pow(std::log(P_0_model) - mu_A,2)  +  P_1_null * std::pow(std::log(P_1_model) - mu_A,2);
			
			CLT_sums.sum_mu += mu_A;
			CLT_sums.sn_2 += var_A;
									
			if (bernoulli_observations_by_position.at(le) > 0)
			{  CLT_sums.sum_log_P_obs += std::log(P_1_model);  }
			else
			{  CLT_sums.sum_log_P_obs += std::log(P_0_model);  }
			

			if (isnan(mu_A) or isnan(var_A))
			{
				std::cerr << "\nNAN encountered!\n"
						<< "le = " << le << "\n"
						<< "mu_A = " << mu_A << "\n"
						<< "var_A = " << var_A << "\n";
						
				std::cerr << "\tbernoulli_p_by_position__and_readgroup.at(le).size() = " << poisson_lambda_summed_over_readgroup_per_position.at(le).size() << "\n";
				for (uint rg = 0; rg < poisson_lambda_summed_over_readgroup_per_position.at(le).size(); ++rg)
				{
					std::cerr << "\t\t"
    // 						<<"bernoulli_p_by_position.at(le) = " << bernoulli_p_by_position.at(le) << "\n"
						<< "poisson_lambda_summed_over_readgroup_per_position.at(le).at(rg) = " << poisson_lambda_summed_over_readgroup_per_position.at(le).at(rg) << "\n";
				}
			}
		}
	}//le	
	
	
	if (number_of_good_positions > 0)
	{
		assert(CLT_sums.sn_2 > 0);
		return   CLT_sums;
	}
	else
	{  return CLT_sums;  }
	
	
}//calculate_Lyapunov_CLT_sum








































void sample_control_positions_from_across_genome
					(const uint &number_of_samples,
					 const type_map_uint_to_list_BI &NON_unique_regions_of_genome,
					BamTools::BamReader &my_BAM_reader,
					type_map_uint_to_list_uint &sampled_observed_counts_per_GC_content,
					uint &mean_frag_length)
{	
	type_map_uint_to_uint non_sex_chromosome_lengths(Event::chromosome_lengths);
	non_sex_chromosome_lengths.erase(23);
	non_sex_chromosome_lengths.erase(24);
	
	type_map_uint_to_string map_chromo_to_seq;
	
	std::cerr << "\nloading entire genome into strings...\n";	
	for (uint chr_indx = 1; chr_indx < 23; ++chr_indx)
	{
		map_chromo_to_seq[chr_indx] = get_specific_DNA_sequence_from_Reference_genome(chr_indx, BOOST_Interval(1, non_sex_chromosome_lengths.at(chr_indx)));				
	}//chr_indx
	
	
	std::cerr << "making cumulative genome map to chromosome...\n";
	uint length_of_haploid_nonsex_genome = 0;
	type_map_BI_to_uint map_cumulative_genome_position_to_chromo;
	{
		uint cumulative_index = 0;
		for (type_map_uint_to_uint::const_iterator it_chr = non_sex_chromosome_lengths.begin();
			it_chr != non_sex_chromosome_lengths.end();
			++it_chr)
		{
			map_cumulative_genome_position_to_chromo[BOOST_Interval(cumulative_index + 1, cumulative_index + it_chr->second)] = it_chr->first;
			cumulative_index += it_chr->second;
			length_of_haploid_nonsex_genome += it_chr->second;
		}//chromosome
	}
	
	std::cerr << "length_of_haploid_nonsex_genome = " << length_of_haploid_nonsex_genome << "\n";
	
	
	{//mean_frag_length
		uint weighted_frag_length_total = 0;
		longdouble mean_frag_len__exact = 0;
		
		for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
			it_rg != Readgroup_stats.end();
			++it_rg)
		{
			weighted_frag_length_total += it_rg->second.number_mapped_Proper_Pairs_reads;
		}//rg
		
		for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
			it_rg != Readgroup_stats.end();
			++it_rg)
		{
			mean_frag_len__exact += (longdouble)it_rg->second.median_insert_size * ((longdouble)it_rg->second.number_mapped_Proper_Pairs_reads / weighted_frag_length_total);
		}//rg
		
		mean_frag_length = (uint)mean_frag_len__exact;
	}//mean_frag_length
	std::cerr << "mean_frag_length = " << mean_frag_length << "\n";

	
		
	const uint skip_amount_at_end_of_chromos = std::max<uint>(500, mean_frag_length);
	
	
	
	boost::random::mt19937 random_generator;	
	boost::random::uniform_int_distribution<uint> uniform_genome_dist(1, length_of_haploid_nonsex_genome);
	
	
	#pragma omp parallel for schedule(dynamic,5)
	for (uint samp=0; samp<number_of_samples; ++samp)
	{
		uint chr_of_random_draw;
		uint chromo_position_of_random_draw;
		
		bool draw_is_acceptable = false;
		uint num_attempts = 0;
		while (!draw_is_acceptable)
		{
// 			std::cerr << "\tattempting successful draw...\n";
			const uint cumulative_genome_position_random_draw = uniform_genome_dist(random_generator);
// 			std::cerr << "\t\tcumulative_genome_position_random_draw = " << cumulative_genome_position_random_draw << "\n";
		
			++num_attempts;
			chr_of_random_draw = 0;
			chromo_position_of_random_draw = 0;
			//get chromo and position			
			for (type_map_BI_to_uint::const_iterator it_cumu = map_cumulative_genome_position_to_chromo.begin();
				it_cumu != map_cumulative_genome_position_to_chromo.end();
				++it_cumu)
			{
				if (BOOST_in(cumulative_genome_position_random_draw, it_cumu->first))
				{
					chr_of_random_draw = it_cumu->second;
					chromo_position_of_random_draw = cumulative_genome_position_random_draw - it_cumu->first.lower() + 1;
					break;
				}//inside
			}//cumulative		
			assert(chr_of_random_draw > 0);
			assert(chromo_position_of_random_draw > 0);
			
			const BOOST_Interval acceptable_p_arm(skip_amount_at_end_of_chromos, Event::centromere_coordinates.at(chr_of_random_draw).first - skip_amount_at_end_of_chromos);
			
			const BOOST_Interval acceptable_q_arm(Event::centromere_coordinates.at(chr_of_random_draw).second + skip_amount_at_end_of_chromos,
											non_sex_chromosome_lengths.at(chr_of_random_draw) - skip_amount_at_end_of_chromos);
			
// 			std::cerr << "\t\tcheck_if_position_is_in_NON_unique_region:  " << check_if_position_is_in_NON_unique_region(NON_unique_regions_of_genome, chr_of_random_draw, chromo_position_of_random_draw) << "\n"
// 				<< "\t\tBOOST_in(chromo_position_of_random_draw, acceptable_p_arm) = " << BOOST_in(chromo_position_of_random_draw, acceptable_p_arm) << "\n"
// 				<< "\t\tBOOST_in(chromo_position_of_random_draw, acceptable_q_arm) = " << BOOST_in(chromo_position_of_random_draw, acceptable_q_arm) << "\n\n";
			
			draw_is_acceptable = (!check_if_position_is_in_NON_unique_region(NON_unique_regions_of_genome, chr_of_random_draw, chromo_position_of_random_draw)
						and (BOOST_in(chromo_position_of_random_draw, acceptable_p_arm) or BOOST_in(chromo_position_of_random_draw, acceptable_q_arm)));
			
			assert(num_attempts < 20);
			
		}//unacceptable draw
			
		
		const std::string test_seq(map_chromo_to_seq.at(chr_of_random_draw).substr(chromo_position_of_random_draw-1, mean_frag_length));
		const uint gc_count__random = calculate_GC_content_of_sequence(test_seq);		
		
		
		uint number_observed_reads;
		
		#pragma omp critical(random_sample_matched_control_position_from_good_part_of_genome)
		{
			my_BAM_reader.SetRegion(create_BAM_region(chr_of_random_draw, BOOST_Interval(chromo_position_of_random_draw-1, chromo_position_of_random_draw+1), 2));	
			number_observed_reads = count_decent_reads_with_left_endpoints_in_set_region(my_BAM_reader, BOOST_make_point(chromo_position_of_random_draw));		
		}				
		
		sampled_observed_counts_per_GC_content[gc_count__random].push_back(number_observed_reads);
		
	}//samp		
	
}//sample_control_positions_from_across_genome

























bool check_if_position_is_in_NON_unique_region
	(const type_map_uint_to_list_BI &nonuniq_regions,
	const uint &chromo,
	const uint &pos_on_chromo)
{
	const type_map_uint_to_list_BI::const_iterator it_chr = nonuniq_regions.find(chromo);
	if (it_chr != nonuniq_regions.end())
	{
		for (type_list_BI::const_iterator it_nonuniq = it_chr->second.begin();
			it_nonuniq != it_chr->second.end();
			++it_nonuniq)
		{
			if (BOOST_in(pos_on_chromo, *it_nonuniq))
			{  return true;  }
			else if (pos_on_chromo < it_nonuniq->lower())
			{  return false;  }
		}//nonuniq
	}//chr
	
	return false;
	
}//check_if_position_is_in_NON_unique_region









bool check_if_region_overlaps_NON_unique_region
	(const type_map_uint_to_list_BI &nonuniq_regions,
	const uint &chromo,
	const BOOST_Interval &reg_on_chromo)
{
	const type_map_uint_to_list_BI::const_iterator it_chr = nonuniq_regions.find(chromo);
	if (it_chr != nonuniq_regions.end())
	{
		for (type_list_BI::const_iterator it_nonuniq = it_chr->second.begin();
			it_nonuniq != it_chr->second.end();
			++it_nonuniq)
		{
			if (BOOST_overlap(reg_on_chromo, *it_nonuniq))
			{  return true;  }
			else if (reg_on_chromo.upper() < it_nonuniq->lower())
			{  return false;  }
		}//nonuniq
	}//chr
	
	return false;
	
}//check_if_region_overlaps_NON_unique_region































// void sample_control_positions_for_all_genomes
// 				(const std::string &postdatatdir)
// {		
// 	const type_map_uint_to_list_BI nonuniqe_regions_of_genome(load_Non_unique_regions_of_genome());
// 	
// 	print_line_of_markers("==");
// 	std::cerr << " \n\n\n\n\n\nsample_control_poisitons_for_all_genomes:\n\n";	
// 
// 	std::string calls_fname(postdatatdir);
// 	calls_fname.append("/all_calls");
// 	
// 	type_map_string_to_uint_to_Sampled_diploid_Event_data map_genome_to_event_to_sampled_Event(read_Sampled_Events_from_file(calls_fname));
// 	
// 	const std::vector<type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator>  loopable_genome_to_sampled_Events(		
// 							get_loopable_iterators<type_map_string_to_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event));
// 	
// 	
// 	print_set<std::string>(extract_keys_of_map_and_return_as_set<std::string, type_map_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event),
// 							"genomes to loop over-MPI", NULL, true);	
// 	
// 	
// 	// MPI
// 	for (int genrank = my_MPI_rank; genrank < loopable_genome_to_sampled_Events.size(); genrank += size_MPI_WORLD)
// 	{
// 		std::cerr << "\n\n\nrank " << my_MPI_rank << " looking for next genome...\n\n";    
// 		
// 		const type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator it_gen = loopable_genome_to_sampled_Events.at(genrank);
// 	
// 		genome_name = it_gen->first;
// 		
// 		std::string full_job_dirname;
// 		{//full_job_dirname
// 			const boost::filesystem::directory_iterator it_end;  
// 			for (boost::filesystem::directory_iterator it_jobdir(postdatatdir);
// 				it_jobdir != it_end;
// 				++it_jobdir)
// 			{		
// 				if (boost::filesystem::is_directory(it_jobdir->path()))
// 				{
// 					if (it_jobdir->path().filename().string().find(it_gen->first) != std::string::npos)
// 					{
// 						full_job_dirname = boost::filesystem::complete(it_jobdir->path()).string();
// 						break;			
// 					}		    		    
// 				}//dir				
// 			}//it_jobdir	    	    	    
// 		}//full_job_dirname
// 		
// 		
// 		std::cerr << "\n\n\nrank " << my_MPI_rank << "\tfull_job_dirname = [" << full_job_dirname << "]\n";
// 		
// 		output_dir = full_job_dirname;
// 	
// 		
// 		std::string bam_filename;	
// 		{//full_job_dirname
// 			const boost::filesystem::directory_iterator it_end;  
// 			for (boost::filesystem::directory_iterator it_bamdir("/gpfs/scratch/mmparks/SAMs_and_BAMs");
// 				it_bamdir != it_end;
// 				++it_bamdir)
// 			{		
// 				if (!boost::filesystem::is_directory(it_bamdir->path()))
// 				{
// 					if (it_bamdir->path().filename().string().find(it_gen->first) != std::string::npos
// 					and  it_bamdir->path().filename().string().find(".bas") == std::string::npos
// 					and  it_bamdir->path().filename().string().find(".bai") == std::string::npos
// 					and  it_bamdir->path().filename().string().find(".bam") == it_bamdir->path().filename().string().size() - 4)			
// 					{
// 						bam_filename = boost::filesystem::complete(it_bamdir->path()).string();
// 						break;			
// 					}		    		    
// 				}//file
// 			}//it_jobdir	    	    	    
// 		}//full_job_dirname	
// 		
// 		std::cerr << "\n\n\nrank " << my_MPI_rank << "\ttbam_filename = [" << bam_filename << "]\n\n";
// 		
// 			
// 		
// 		assert(!full_job_dirname.empty());
// 		assert(!bam_filename.empty());			
// 		
// 
// 		
// 		BamTools::BamReader my_BAM_reader;
// 
// 		prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population
// 						(bam_filename,
// 						std::string(),//gender_and_population_filename,
// 						my_BAM_reader);		
// 
// 
// 		read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda(bam_filename, false, NULL);
// 
// 		
// 		uint mean_frag_length;
// 		type_map_uint_to_list_uint sampled_observed_counts_per_GC_content;
// 		const uint number_of_samples = 20000;
// 		
// 		sample_control_positions_from_across_genome
// 					(number_of_samples,
// 					 nonuniqe_regions_of_genome,
// 					my_BAM_reader,
// 					sampled_observed_counts_per_GC_content,
// 					mean_frag_length);
// 					
// 		
// 		{//check sampled controls
// 			char common_fname[option_size];
// 			std::sprintf(common_fname, "%s/sampled_control_positions.%s", output_dir.c_str(), genome_name.c_str());
// 			
// 			std::ofstream outfs(common_fname);
// 			boost::archive::text_oarchive outarch(outfs);
// 
// 			outarch << mean_frag_length;
// 			outarch << sampled_observed_counts_per_GC_content;
// 
// 			outfs.close();
// 		}//check sampled controls
// 		
// 		
// 		std::cerr << "\nrank " << my_MPI_rank << " moving to next genome...\n";
// 		
//     }//it_gen            
//     
//     
//     std::cerr << "\nrank " << my_MPI_rank << "  DONE with \"sample_control_positions_for_all_genomes\"\n\n";
// 	
// }//sample_control_positions_for_all_genomes
// 





































// void  bernoulli_lyapunov_nonparametric_test_with_samples
// 				(const std::string &postdatatdir,
// 				 type_map_uint_to_CC &Conn_Comps)
// {
// 	
// 	const type_map_uint_to_list_BI nonuniqe_regions_of_genome(load_Non_unique_regions_of_genome());
// 	
// 	print_line_of_markers("==");
// 	std::cerr << " \n\n\n\n\n\nbernoulli_lyapunov_nonparametric_test_with_samples:\n\n";	
// 
// 	std::string calls_fname(postdatatdir);
// 	calls_fname.append("/all_calls");
//     
// 	type_map_string_to_uint_to_Sampled_diploid_Event_data map_genome_to_event_to_sampled_Event(read_Sampled_Events_from_file(calls_fname));
// 	
// 	const std::vector<type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator>  loopable_genome_to_sampled_Events(		
// 							get_loopable_iterators<type_map_string_to_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event));
// 	
// 	
// 	print_set<std::string>(extract_keys_of_map_and_return_as_set<std::string, type_map_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event),
// 							"genomes to loop over-MPI", NULL, true);	
// 	
// 	
// 	
// 	// MPI
// 	for (int genrank = my_MPI_rank; genrank < loopable_genome_to_sampled_Events.size(); genrank += size_MPI_WORLD)
// 	{
// 		std::cerr << "\n\n\nrank " << my_MPI_rank << " looking for next genome...\n\n";
//     	
// 		const type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator it_gen = loopable_genome_to_sampled_Events.at(genrank);
// 	
// 		genome_name = it_gen->first;
// 		
// 		std::string full_job_dirname;
// 		{//full_job_dirname
// 			const boost::filesystem::directory_iterator it_end;  
// 			for (boost::filesystem::directory_iterator it_jobdir(postdatatdir);
// 				it_jobdir != it_end;
// 				++it_jobdir)
// 			{		
// 				if (boost::filesystem::is_directory(it_jobdir->path()))
// 				{
// 					if (it_jobdir->path().filename().string().find(it_gen->first) != std::string::npos)
// 					{
// 						full_job_dirname = boost::filesystem::complete(it_jobdir->path()).string();
// 						break;			
// 					}		    		    
// 				}//dir				
// 			}//it_jobdir	    	    	    
// 		}//full_job_dirname
// 		
// 		
// 		std::cerr << "\n\n\nrank " << my_MPI_rank << "\tfull_job_dirname = [" << full_job_dirname << "]\n";
// 		
// 		output_dir = full_job_dirname;
// 	
// 		
// 		std::string bam_filename;	
// 		{//full_job_dirname
// 			const boost::filesystem::directory_iterator it_end;  
// 			for (boost::filesystem::directory_iterator it_bamdir("/gpfs/scratch/mmparks/SAMs_and_BAMs");
// 				it_bamdir != it_end;
// 				++it_bamdir)
// 			{		
// 				if (!boost::filesystem::is_directory(it_bamdir->path()))
// 				{
// 					if (it_bamdir->path().filename().string().find(it_gen->first) != std::string::npos
// 					and  it_bamdir->path().filename().string().find(".bas") == std::string::npos
// 					and  it_bamdir->path().filename().string().find(".bai") == std::string::npos
// 					and  it_bamdir->path().filename().string().find(".bam") == it_bamdir->path().filename().string().size() - 4)			
// 					{
// 						bam_filename = boost::filesystem::complete(it_bamdir->path()).string();
// 						break;			
// 					}		    		    
// 				}//file
// 			}//it_jobdir	    	    	    
// 		}//full_job_dirname	
// 		
// 		std::cerr << "\n\n\nrank " << my_MPI_rank << "\ttbam_filename = [" << bam_filename << "]\n\n";
// 		
// 			
// 		
// 		assert(!full_job_dirname.empty());
// 		assert(!bam_filename.empty());			
// 		
// 		
// 		{
// 			char common_fname[option_size];
// 			std::sprintf(common_fname, "%s/Wilcoxon_signed_rank_sum_test.ratios.matched", output_dir.c_str());
// 			boost::filesystem::remove(common_fname);			
// 		}
// 		
// 		
// 		BamTools::BamReader my_BAM_reader;
// 
// 		prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population
// 						(bam_filename,
// 						std::string(),//gender_and_population_filename,
// 						my_BAM_reader);		
// 
// 
// 		read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda(bam_filename, false, NULL);
// 		std::cerr.flush();
// 		std::cerr << "\n\n\nrank " << my_MPI_rank << "\tprepare-done.";
// 		
// 
// 	
// 		
// 		type_map_uint_to_list_uint sampled_control_bins_by_GC_count;
// 		
// 		print_map_to_list<uint,uint>(sampled_control_bins_by_GC_count, "sampled_control_bins_by_GC_count", NULL, true);
// 		
// 		{//check sampled controls
// 			char common_fname[option_size];
// 			std::sprintf(common_fname, "%s/sampled_control_positions.%s", output_dir.c_str(), genome_name.c_str());
// 			
// 			if (boost::filesystem::exists(common_fname))
// 			{
// 				uint mean_frag_len;
// 				std::cerr << "\n\n\nrank " << my_MPI_rank << "\tuploading already-saved sampled controls for genome " << it_gen->first << "...\n";
// 				std::ifstream infs(common_fname);
// 				boost::archive::text_iarchive inarch(infs);
// 				
// 				inarch >> mean_frag_len;
// 				inarch >> sampled_control_bins_by_GC_count;
// 				
// 				infs.close();
// 			}
// 			
// 			if (sampled_control_bins_by_GC_count.empty())
// 			{
// 				std::cerr << "\nERROR!  sampled_control_bins_by_GC_count is empty for genome " << it_gen->first << "\n\n";
// 			}//sample
// 		}//check sampled controls
// 		
// 
// 		
// 		
// 		
// 		
// 		const std::vector<type_map_uint_to_Sampled_diploid_Event_data::iterator>  myev__loop(
// 						get_loopable_iterators<type_map_uint_to_Sampled_diploid_Event_data>(it_gen->second));
// 	
// 
// 		#pragma omp parallel
// 		{
// 			uint thread_work_ctr = 0;  //for debug
// 			
// 			#pragma omp for schedule(static)
// 			for (int myev_ctr = 0; myev_ctr < myev__loop.size(); ++myev_ctr)
// 			{    				
// 				std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() 
// 						<< "  myev_ctr = " << myev_ctr << " out of " << myev__loop.size() << ".\n";
// 						
// 						
// 				std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "-work-counter = " << thread_work_ctr << "\n";
// 				++thread_work_ctr;
// 				
// 				const type_map_uint_to_Sampled_diploid_Event_data::iterator it_myev = myev__loop.at(myev_ctr);
// 				
// 				
// 				
// 							
// 				const type_map_uint_to_Event::iterator result_Ev_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID); 
// 				
// 				
// 				if (result_Ev_it->second.chromos[0] == 24  and   gender_of_individual_is_female)
// 				{
// 					std::cerr << "Event is on Y chromosome, but individual is female.  Skipping read-depth graph in in \"create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights\".\n";
// 					continue;
// 				}
// 				
// 				const BOOST_Interval expanded_region_of_interest( 
// 							safe_subtract_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.lower(), 100000),
// 							safe_chromo_add_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.upper(), 100000, result_Ev_it->second.chromos[0]));        
// 							
// 				//identify all events intersecting this region:
// 				type_map_uint_to_Event_ptr  all_events_near_expanded_region;
// 				for (type_map_uint_to_Event::iterator it_ev = Conn_Comps.at(it_myev->second.cc_of_event).events.begin();
// 					it_ev != Conn_Comps.at(it_myev->second.cc_of_event).events.end();
// 					++it_ev)
// 				{	    
// 					if (it_ev->second.chromos[0] == result_Ev_it->second.chromos[0]
// 						and  (BOOST_overlap(it_ev->second.LCRs[0],expanded_region_of_interest)
// 							or  BOOST_overlap(it_ev->second.LCRs[1],expanded_region_of_interest)))
// 					{
// 						all_events_near_expanded_region[it_ev->first] = &it_ev->second;
// 					}
// 				}//ev
// 				
// 				print_map_keys<uint, Event*>(all_events_near_expanded_region, "all_events_near_expanded_region", NULL, true);
// 				
// 				
// 				
// 				
// 				
// 				
// 				char out_RD_core_fname[option_size];
// 				
// 				std::sprintf(out_RD_core_fname, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s_%s___%s", 
// 							output_dir.c_str(),
// 							it_myev->second.cc_of_event, result_Ev_it->second.UID,
// 							convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first).c_str(),
// 							convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second).c_str(),
// 							genome_name.c_str());		    
// 						
// 				
// 		// 	    //if already computed, skip)
// 		// 	    {//check already computed
// 		// 		std::string testname(out_RD_core_fname);
// 		// 		testname.append(".tests");
// 		// 		
// 		// 		if (boost::filesystem::exists(testname) != 0)
// 		// 		{
// 		// 		    std::ifstream infs(testname);
// 		// 		    
// 		// 		    std::string inteststr(option_size*10, '\0');
// 		// 		    char intestfile[option_size*10];
// 		// 		    std::sprintf(intestfile, "%s", inteststr.c_str());
// 		// 		    		    
// 		// 		    infs.read(intestfile, option_size*10);
// 		// 		    infs.close();
// 		// 		    
// 		// 		    inteststr.assign(intestfile);
// 		// 		    const uint num_newlines = std::count(inteststr.begin(), inteststr.end(), '\n');		    
// 		// 		    if (num_newlines >= 31)
// 		// 		    {
// 		// 			std::cerr << "\n\tnum_newlines = " << num_newlines << " !!!  continue!\n";
// 		// 			continue;			
// 		// 		    }
// 		// 		}//does not exit
// 		// 	    }//check already computed
// 				
// 				
// 				
// 				
// 				std::string observed_fname(out_RD_core_fname);
// 				observed_fname.append(".data.observed");	    
// 				
// 				if (boost::filesystem::exists(observed_fname) == 0)
// 				{//search
// 					std::string searchdir(output_dir);
// 					searchdir.append("/Read_Depth");
// 									
// 					std::stringstream cc_strm;
// 					cc_strm << "cc_" << it_myev->second.cc_of_event << "_";
// 					std::stringstream ev_strm;
// 					ev_strm << "ev_" << it_myev->second.event_UID << "_";
// 					
// 					char observed_fname_old[option_size];
// 					
// 					bool found_something = false;
// 					const boost::filesystem::directory_iterator it_end;  
// 					for (boost::filesystem::directory_iterator it_rdfile(searchdir);
// 						it_rdfile != it_end;
// 						++it_rdfile)
// 					{		
// 						if (!boost::filesystem::is_directory(it_rdfile->path()))
// 						{
// 							if (it_rdfile->path().filename().string().find(cc_strm.str()) != std::string::npos
// 								and  it_rdfile->path().filename().string().find(ev_strm.str()) != std::string::npos
// 								and  it_rdfile->path().filename().string().find(genome_name) != std::string::npos
// 								and it_rdfile->path().filename().string().find("RD_hypothesis_test__cc_") != std::string::npos)
// 							{
// 								char scanformat[option_size];
// 								std::sprintf(scanformat, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%%s", 
// 											output_dir.c_str(),
// 											it_myev->second.cc_of_event,
// 											it_myev->second.event_UID);	
// 								
// 								char in_calls_and_genome[option_size];		    
// 												
// 				// 			    std::cerr << "\n\t\t\t\ttry name: [" << boost::filesystem::complete(it_rdfile->path()).string().c_str() << "\n\t\t\t\tscanformat = [ " << scanformat << "]\n";
// 								std::sscanf(boost::filesystem::complete(it_rdfile->path()).string().c_str(), scanformat,
// 										in_calls_and_genome);
// 								
// 								std::string calls_CHOPPED(in_calls_and_genome);
// 								calls_CHOPPED = calls_CHOPPED.substr(0,calls_CHOPPED.find(genome_name));
// 								
// 				// 			    std::cerr << "\n\t\t\t\tcalls_CHOPPED = [" << calls_CHOPPED << "]\n";
// 								
// 								std::sprintf(observed_fname_old, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s%s.data.observed", 
// 									output_dir.c_str(),
// 									it_myev->second.cc_of_event, result_Ev_it->second.UID,
// 									calls_CHOPPED.c_str(),
// 									genome_name.c_str());	    
// 								
// 								found_something = true;
// 								break;			
// 							}		    		    
// 						}//file
// 					}//it_jobdir
// 					
// 					
// 					if (!found_something)
// 					{
// 						std::stringstream error_strm;
// 						error_strm << "ERROR  unable to find anything!!!\n\n";
// 						error_message(error_strm,false);
// 						continue;
// 					}	    
// 					
// 					std::cerr << "\n\t\t\tout_RD_core_fname = [" << out_RD_core_fname << "]\n"
// 						<<  "\t\t\tobserved_fname_old = [" << observed_fname_old << "\n\n";
// 											
// 					if (boost::filesystem::exists(observed_fname) == 0)
// 					{  
// 						std::cerr << "rename = [" << observed_fname << "\n";
// 						boost::filesystem::rename(observed_fname_old, observed_fname);		    
// 					}
// 					
// 				}//search
// 
// 				
// 				const type_map_uint_to_longdouble map_Ref_Genome_position_to_per_read_posterior_sum(read_map_from_file<uint,longdouble>(observed_fname, false));
// 			
// 				
// 				
// 				{//RD hypothesis test
// 				
// 					BOOST_Interval affected_region(empty_BI);
// 					
// 					if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None))
// 					{
// 						affected_region = result_Ev_it->second.region_between_and_including_the_LCRs_themselves;	    
// 					}
// 					else
// 					{
// 						for (uint hap=0; hap<2; ++hap)
// 						{
// 							if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) != hap_outcome__None)
// 							{
// 								const type_BI__BI affected_regions_each_LCR(
// 											BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(
// 												convert_each_profile_breakpoint_to_absolute_breakpoints(
// 													pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
// 													result_Ev_it->second.compressed_map_LCR_to_profile__for_each_LCR)));
// 											
// 								const BOOST_Interval affected_hull(affected_regions_each_LCR.first.lower(), affected_regions_each_LCR.second.upper());
// 																							
// 								if (BOOST_empty(affected_region))
// 								{  affected_region = affected_hull;  }
// 								else
// 								{  affected_region = BOOST_hull(affected_region, affected_hull);  }
// 							}//occurs
// 						}//hap
// 					}
// 					
// 					
// 					
// 					std::stringstream  out_ss;
// 				
// 					out_ss  << "\nrank " << my_MPI_rank << ", thread " << omp_get_thread_num() 
// 							<< "\n\tmyev_ctr = " << myev_ctr
// 							<< "\n\tgenome = " << genome_name
// 							<< "\n\tcc = " << it_myev->second.cc_of_event
// 							<< "\n\tevent = " << it_myev->second.event_UID
// 							<< "\n\taffected_region = " << affected_region 
// 							<< "\n\tchrome = " <<  result_Ev_it->second.chromos[0]
// 							<< "\n\tchrom-length = " << Event::chromosome_lengths.at(result_Ev_it->second.chromos[0]) 
// 							<< "\n\toutcome = " << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first) 
// 							<< ", " << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second) 
// 							<< "\n";
// 
// 
// 					//restrict to affected region:
// 					const type_map_uint_to_longdouble  map_affected_Ref_Genome_position_to_per_read_posterior_sum(
// 							get_map_intersection_of_map_keys_and_interval<uint,longdouble>(map_Ref_Genome_position_to_per_read_posterior_sum, affected_region));
// 							
// 
// // 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get affected_sequence__padded\n";
// // 					std::cerr << "\ncalculating Wilcoxon_pvalue__inferred...\n";
// 					const uint pad_amount = 5000;
// 					const uint padded_max_chromo_sensitive = safe_chromo_add_base_1(affected_region.upper(), pad_amount, result_Ev_it->second.chromos[0]);
// 					//std::max<uint>(affected_region.upper() + pad_amount, );
// 					
// 					const std::string affected_sequence__padded(get_specific_DNA_sequence_from_Reference_genome(result_Ev_it->second.chromos[0],
// 														BOOST_Interval(affected_region.lower(), padded_max_chromo_sensitive)));
// 					
// 					
// // 					type_vector_longdouble bernoulli_p_by_position(BOOST_width_inclusive(affected_region), 0);
// 					type_vector_vector_longdouble bernoulli_p_by_position_and_readgroup___NULL(BOOST_width_inclusive(affected_region),
// 																    type_vector_longdouble(Readgroup_stats.size(), 0));
// 					
// 					type_vector_uint bernoulli_observations_by_position(BOOST_width_inclusive(affected_region), 0);
// 					
// // 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get bernoulli_p_by_position\n";
// 					//get bernoulli p at each position by summing over all readgroups, using GC content at that position.
// 					for (uint le = 0; le < BOOST_width_inclusive(affected_region); ++le)
// 					{
// 						uint rg_indx = 0;
// 						for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
// 							it_rg != Readgroup_stats.end();
// 							++it_rg)
// 						{
// 							const uint gc_content = calculate_GC_content_of_sequence(affected_sequence__padded.substr(le, it_rg->second.median_insert_size));
// 							
// 							bernoulli_p_by_position_and_readgroup___NULL.at(le).at(rg_indx)
// 									= it_rg->second.haploid_fragmentation_rate__per_GC_count.at(gc_content) * 2; // null case is diploid
// 							
// 							++rg_indx;
// 						}
// 					}//bernoulli p
// 					
// 					
// // 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get bernoulli_observations_by_position\n";
// 					//determine bernoulli observation at each point.
// 					uint num_positives = 0;
// 					for (type_map_uint_to_longdouble::const_iterator it_obs = map_affected_Ref_Genome_position_to_per_read_posterior_sum.begin();
// 						it_obs != map_affected_Ref_Genome_position_to_per_read_posterior_sum.end();
// 						++it_obs)
// 					{
// 						if (it_obs->second > 0.5)
// 						{
// 							bernoulli_observations_by_position.at(it_obs->first - affected_region.lower()) = 1;
// 							++num_positives;
// 						}
// 					}//observed
// 					
// 					
// 					
// 					out_ss << "\n\tnum_obs/num_pos = " << (longdouble)num_positives / map_affected_Ref_Genome_position_to_per_read_posterior_sum.size() << "\n\n";
// 					
// 					
// 					const Lyap_CLT_data stats_null_inferred(calculate_Lyapunov_CLT_sum(bernoulli_p_by_position_and_readgroup___NULL, bernoulli_observations_by_position));
// 					
// 					
// // 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get poisson_lambda_by_position___alternative\n";
// 					type_vector_vector_longdouble poisson_lambda_by_position___alternative(bernoulli_p_by_position_and_readgroup___NULL);
// 					
// 					if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
// 					{
// 						out_ss << "\talternative = diploid dup\n";
// 						for (uint le=0; le < poisson_lambda_by_position___alternative.size(); ++le)
// 						{  
// 							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position___alternative.at(le).size(); ++rg_indx)
// 							{  poisson_lambda_by_position___alternative.at(le).at(rg_indx) *= 2;  }
// 						}													
// 					}//diploid dup
// 					else if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
// 					{
// 						out_ss << "\talternative = diploid del\n";
// 						for (uint le=0; le < poisson_lambda_by_position___alternative.size(); ++le)
// 						{  
// 							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position___alternative.at(le).size(); ++rg_indx)
// 							{  poisson_lambda_by_position___alternative.at(le).at(rg_indx) /= 10;  }
// 						}													
// 					}//diploid del		
// 					else if (test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
// 					{
// 						out_ss << "\talternative = haploid dup\n";
// 						for (uint le=0; le < poisson_lambda_by_position___alternative.size(); ++le)
// 						{  
// 							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position___alternative.at(le).size(); ++rg_indx)
// 							{  poisson_lambda_by_position___alternative.at(le).at(rg_indx) *= 1.5;  }
// 						}
// 					}//haploid dup
// 					else if (test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
// 					{
// 						out_ss << "\talternative = haploid del\n";
// 						for (uint le=0; le < poisson_lambda_by_position___alternative.size(); ++le)
// 						{  
// 							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position___alternative.at(le).size(); ++rg_indx)
// 							{  poisson_lambda_by_position___alternative.at(le).at(rg_indx) /= 2;  }
// 						}
// 					}//haploid del
// 					
// 					
// 					
// 					
// 					
// 					Lyap_CLT_data stats_alt_inferred;
// 					bool best_alt_is_del = false;
// 					real Lyap_CLT_tail_alt___inferred(-1);
// 					
// 					if (test_if_at_least_one_haploid_is_NAHR(it_myev->second.the_diploid_Event_outcome))
// 					{		
// 						stats_alt_inferred = calculate_Lyapunov_CLT_sum(poisson_lambda_by_position___alternative, bernoulli_observations_by_position);
// 						
// 						boost::math::normal_distribution<mpfr_class> my_normal_Lyap_alt((double)(stats_null_inferred.sum_mu - stats_alt_inferred.sum_mu),
// 														(double)std::sqrt(stats_null_inferred.sn_2 + stats_alt_inferred.sn_2));
// 						
// 						Lyap_CLT_tail_alt___inferred = real(boost::math::cdf<mpfr_class>(my_normal_Lyap_alt,
// 											(double)(stats_null_inferred.sum_log_P_obs - stats_alt_inferred.sum_log_P_obs)).get_mpfr_t());
// 											
// 						out_ss
// 							<< "\tgiven-alternative:    "
// 							<< Lyap_CLT_tail_alt___inferred 
// 							<< "  =  N(  " << stats_null_inferred.sum_log_P_obs << " - " << stats_alt_inferred.sum_log_P_obs
// 								<< " = " << (stats_null_inferred.sum_log_P_obs - stats_alt_inferred.sum_log_P_obs) 
// 							<< "  ;  " << stats_null_inferred.sum_mu << " - " << stats_alt_inferred.sum_mu	
// 								<< " = " << (stats_null_inferred.sum_mu - stats_alt_inferred.sum_mu)
// 							<< "  ,  " <<  stats_null_inferred.sn_2 << " + " << stats_alt_inferred.sn_2 
// 								<< " = " << (stats_null_inferred.sn_2 + stats_alt_inferred.sn_2)
// 							<< "  )\n";
// 											
// 					}
// 					else
// 					{
// 						type_vector_vector_longdouble poisson_lambda_by_position__alternative__hap_del(bernoulli_p_by_position_and_readgroup___NULL);
// 						for (uint le=0; le < poisson_lambda_by_position__alternative__hap_del.size(); ++le)
// 						{  
// 							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position__alternative__hap_del.at(le).size(); ++rg_indx)
// 							{  poisson_lambda_by_position__alternative__hap_del.at(le).at(rg_indx) /= 2;  }
// 						}
// 												
// 						type_vector_vector_longdouble poisson_lambda_by_position__alternative__hap_dup(bernoulli_p_by_position_and_readgroup___NULL);
// 						for (uint le=0; le < poisson_lambda_by_position__alternative__hap_dup.size(); ++le)
// 						{  
// 							for (uint rg_indx = 0; rg_indx < poisson_lambda_by_position__alternative__hap_dup.at(le).size(); ++rg_indx)
// 							{  poisson_lambda_by_position__alternative__hap_dup.at(le).at(rg_indx) *= 1.5;  }
// 						}
// 						
// 						
// 						const Lyap_CLT_data stats_alt_del_inferred(
// 									calculate_Lyapunov_CLT_sum(poisson_lambda_by_position__alternative__hap_del, bernoulli_observations_by_position));
// 						const Lyap_CLT_data stats_alt_dup_inferred(
// 									calculate_Lyapunov_CLT_sum(poisson_lambda_by_position__alternative__hap_dup, bernoulli_observations_by_position));
// 						
// 						
// 						boost::math::normal_distribution<mpfr_class> my_normal_Lyap_del((double)(stats_null_inferred.sum_mu - stats_alt_del_inferred.sum_mu),
// 														 (double)std::sqrt(stats_null_inferred.sn_2 + stats_alt_del_inferred.sn_2));
// 						
// 						boost::math::normal_distribution<mpfr_class> my_normal_Lyap_dup((double)(stats_null_inferred.sum_mu - stats_alt_dup_inferred.sum_mu),
// 														 (double)std::sqrt(stats_null_inferred.sn_2 + stats_alt_dup_inferred.sn_2));
// 						
// 						const real Lyap_CLT_tail___del(boost::math::cdf<mpfr_class>(my_normal_Lyap_del,
// 										(double)(stats_null_inferred.sum_log_P_obs - stats_alt_del_inferred.sum_log_P_obs)).get_mpfr_t());
// 						
// 						const real Lyap_CLT_tail___dup(boost::math::cdf<mpfr_class>(my_normal_Lyap_del,
// 										(double)(stats_null_inferred.sum_log_P_obs - stats_alt_dup_inferred.sum_log_P_obs)).get_mpfr_t());
// 						
// 						if (Lyap_CLT_tail___del < Lyap_CLT_tail___dup)
// 						{
// 							best_alt_is_del = true;
// 							stats_alt_inferred = stats_alt_del_inferred;
// 							Lyap_CLT_tail_alt___inferred = Lyap_CLT_tail___del;
// 						}
// 						else
// 						{
// 							best_alt_is_del = false;
// 							stats_alt_inferred = stats_alt_dup_inferred;
// 							Lyap_CLT_tail_alt___inferred = Lyap_CLT_tail___dup;
// 						}
// 						
// 						out_ss
// 							<< "\t"
// 							<< Lyap_CLT_tail_alt___inferred 
// 							<< "  =  N(  " << stats_null_inferred.sum_log_P_obs << " - " << stats_alt_inferred.sum_log_P_obs
// 								<< " = " << (stats_null_inferred.sum_log_P_obs - stats_alt_inferred.sum_log_P_obs) 
// 							<< "  ;  " << stats_null_inferred.sum_mu << " - " << stats_alt_inferred.sum_mu	
// 								<< " = " << (stats_null_inferred.sum_mu - stats_alt_inferred.sum_mu)
// 							<< "  ,  " <<  stats_null_inferred.sn_2 << " + " << stats_alt_inferred.sn_2 
// 								<< " = " << (stats_null_inferred.sn_2 + stats_alt_inferred.sn_2)
// 							<< "  )\n";
// 						
// 					}//find best alt
// 					
// 														
// 
// // 					std::cerr << "get is_unique\n";	
// 					type_vector_bool is_unique(BOOST_width_inclusive(affected_region), true);
// 					{
// 						const type_map_uint_to_list_BI::const_iterator it_found_chr = nonuniqe_regions_of_genome.find(result_Ev_it->second.chromos[0]); 
// 						if (it_found_chr != nonuniqe_regions_of_genome.end())
// 						{
// 							for (type_list_BI::const_iterator it_bi = it_found_chr->second.begin();
// 								it_bi != it_found_chr->second.end();
// 								++it_bi)
// 							{
// 								const BOOST_Interval intersection(BOOST_intersect(*it_bi, affected_region));
// 								if (!BOOST_empty(intersection))
// 								{
// 									for (uint pos = intersection.lower(); pos <= intersection.upper(); ++pos)
// 									{  is_unique.at(pos - affected_region.lower()) = false;  }						
// 								}														
// 							}
// 						}
// 					}
// 					
// // 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get number_of_unique_positions\n";	
// 					uint number_of_unique_positions = 0;
// 					for (type_vector_bool::const_iterator it_uniq = is_unique.begin(); it_uniq != is_unique.end();  ++it_uniq)
// 					{
// 						if (*it_uniq)
// 						{  ++number_of_unique_positions;  }						
// 					}
// 					
// 					
// // 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get Lyap_CLT_statistic___unique\n";	
// 
// 
// 					Lyap_CLT_data stats_null_unique;
// 					Lyap_CLT_data stats_alt_unique;
// 					real Lyap_CLT_tail___unique(-1);
// 
// 					if (number_of_unique_positions > 100)
// 					{
// 						stats_null_unique = calculate_Lyapunov_CLT_sum(bernoulli_p_by_position_and_readgroup___NULL, bernoulli_observations_by_position, is_unique);
// 						stats_alt_unique  = calculate_Lyapunov_CLT_sum(poisson_lambda_by_position___alternative, bernoulli_observations_by_position, is_unique);
// 						
// 						const longdouble mean_of_diff = stats_null_unique.sum_mu - stats_alt_unique.sum_mu;
// 						const longdouble var_of_diff = stats_null_unique.sn_2 + stats_alt_unique.sn_2;
// 						const longdouble observation = stats_null_unique.sum_log_P_obs - stats_alt_unique.sum_log_P_obs;
// 						
// 						boost::math::normal_distribution<mpfr_class> my_normal_Lyap_uniq((double)mean_of_diff, (double)std::sqrt(var_of_diff));
// 						Lyap_CLT_tail___unique = real(boost::math::cdf<mpfr_class>(my_normal_Lyap_uniq, (double)observation).get_mpfr_t());
// 						
// 						out_ss
// 							<< "\tunique:    "
// 							<< Lyap_CLT_tail___unique 
// 							<< Lyap_CLT_tail_alt___inferred 
// 							<< "  =  N(  " << stats_null_unique.sum_log_P_obs << " - " << stats_alt_unique.sum_log_P_obs
// 								<< " = " << (stats_null_unique.sum_log_P_obs - stats_alt_unique.sum_log_P_obs) 
// 							<< "  ;  " << stats_null_unique.sum_mu << " - " << stats_alt_unique.sum_mu	
// 								<< " = " << (stats_null_unique.sum_mu - stats_alt_unique.sum_mu)
// 							<< "  ,  " <<  stats_null_unique.sn_2 << " + " << stats_alt_unique.sn_2 
// 								<< " = " << (stats_null_unique.sn_2 + stats_alt_unique.sn_2)
// 							<< "  )\n";
// 					}
// 
// 					
// 					
// 						
// 					std::cerr << out_ss.str();
// 					
// 				
// 					#pragma omp critical(write_wilcoxon_signed_rank_to_common_file)
// 					{//common file		
// 						//write to common filename
// 						char common_fname[option_size];
// 						std::sprintf(common_fname, "%s/Wilcoxon_signed_rank_sum_test.lyap.matched", output_dir.c_str());
// 						
// 						std::ofstream out_wilcoxon_fs(common_fname, std::ios_base::app);
// 						
// 						
// 						out_wilcoxon_fs 
// 							//info:
// 							<< genome_name << "\t"
// 							<< it_myev->second.cc_of_event << "\t"
// 							<< it_myev->second.event_UID << "\t"
// 							
// 							//p-values:
// 							<< Lyap_CLT_tail_alt___inferred << "\t"
// 							<< Lyap_CLT_tail___unique << "\n";
// 							
// 						out_wilcoxon_fs.close();
// 					}//common file
// 					
// // 					std::cerr  << "\n\trank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "writing_is_DONE.\n";
// 				
// 				}//RD hypothesis test			
// 			}//it_myev
// 			
// 			std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "  done with it_myev loop...\n";
// 			
// 		}//parallel
// 		std::cerr << "\nrank " << my_MPI_rank << " moving to next genome...\n";
// 		
//     }//it_gen            		
// 	
// }//bernoulli_lyapunov_nonparametric_test_with_samples
// 






// void determine_empirical_bernoulli_distribution_for_null_and_alternative
// 					(const std::string padded_region_sequence,
// 					const BOOST_Interval &affected_region,
// 					const uint &mean_frag_length,
// 					const type_map_uint_to_list_uint &sampled_control_positions_by_GC_content,
// 					const type_haploid_outcome__haploid_outcome &the_diploid_outcome)
// {
// 	
// 	boost::random::mt19937 random_generator;		
// 	boost::random::uniform_int_distribution<> coin_flipper(1,2);	
// 	
// 	for (uint indx = 0; indx < BOOST_width_inclusive(affected_region); ++indx)
// 	{
// 		const uint gc_content = calculate_GC_content_of_sequence(padded_region_sequence.substr(indx, mean_frag_length));
// 		
// 		const type_map_uint_to_list_uint::const_iterator it_gc_sample = sampled_control_positions_by_GC_content.find(gc_content);		
// 		if (it_gc_sample != sampled_control_positions_by_GC_content.end())
// 		{
// 			longdouble empirical_mean__null = 0;
// 			longdouble empirical_mean__del = 0;
// 			longdouble empirical_mean__dup = 0;
// 			longdouble empirical_mean__alt = 0;
// 			uint num_dels_considered = 0;
// 			
// 			for (type_list_uint::const_iterator it_control = it_gc_sample->second.begin(); it_control != it_gc_sample->second.end(); ++it_control)
// 			{  
// 				empirical_mean__null += *it_control;
// 				
// 				if (coin_flipper(random_generator) == 2)
// 				{  empirical_mean__del += *it_control;  }
// 			}//controls
// 			
// 			empirical_mean__null /= it_gc_sample->second.size();
// 			
// 			
// 			
// 			
// 			
// 		
// 		}//gc
// 	}//indx
// 	
// 	
// 	
// 	
// 	
// }//determine_empirical_bernoulli_distribution_for_null_and_alternative
















real compute_normalized_statistic
		(const Lyap_CLT_data &null_data, 
		 const Lyap_CLT_data &alt_data)
{
	const longdouble mu = null_data.sum_mu - alt_data.sum_mu;
	const longdouble sigma = std::sqrt(null_data.sn_2 + alt_data.sn_2);
	const longdouble observed = null_data.sum_log_P_obs - alt_data.sum_log_P_obs;
	
	return  real(observed - mu) / real(sigma);	
	
}//compute_normalized_statistic






real compute_lower_tail_probability
		(const Lyap_CLT_data &null_data, 
		 const Lyap_CLT_data &alt_data)
{
	const longdouble mu = null_data.sum_mu - alt_data.sum_mu;
	const longdouble sigma = std::sqrt(null_data.sn_2 + alt_data.sn_2);
	const longdouble observed = null_data.sum_log_P_obs - alt_data.sum_log_P_obs;
	
	boost::math::normal_distribution<mpfr_class>  my_normal((double)mu, (double)sigma);
	
	return real(boost::math::cdf<mpfr_class>(my_normal, (double)observed).get_mpfr_t());
	
}//compute_lower_tail_probability
















































void  raw_count_ratios
		(const std::string &postdatatdir,
		type_map_uint_to_CC &Conn_Comps)
{
	
// 	//debug
// 	boost::math::normal_distribution<mpfr_class> my_debug_normal((double)914.169, (double)std::sqrt((longdouble)114677));
// 	std::cerr << "\n\ndebug-normal:  boost::math::cdf<longdouble>(my_debug_normal, 1235.91) = " 
// 			<< real(boost::math::cdf<mpfr_class>(my_debug_normal, (double)1235.91).get_mpfr_t()) <<  "\n\n";
// 			
// 	exit(1);
// 	//debug
	
	
	
	const uint binwidth = 100;	
	
	const type_map_uint_to_list_BI nonuniqe_regions_of_genome(load_Non_unique_regions_of_genome());
	
	print_line_of_markers("==");
	std::cerr << " \n\n\n\n\n\ndo_signed_rank_test_for_Read_Depth:\n\n";	

	std::string calls_fname(postdatatdir);
	calls_fname.append("/all_calls");
	
	type_map_string_to_uint_to_Sampled_diploid_Event_data map_genome_to_event_to_sampled_Event(read_Sampled_Events_from_file(calls_fname));
	
	const std::vector<type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator>  loopable_genome_to_sampled_Events(		
							get_loopable_iterators<type_map_string_to_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event));
	
	
	print_set<std::string>(extract_keys_of_map_and_return_as_set<std::string, type_map_uint_to_Sampled_diploid_Event_data>(map_genome_to_event_to_sampled_Event),
							"genomes to loop over-MPI", NULL, true);	
	
	
	
	// MPI
	for (int genrank = my_MPI_rank; genrank < loopable_genome_to_sampled_Events.size(); genrank += size_MPI_WORLD)
	{
		std::cerr << "\n\n\nrank " << my_MPI_rank << " looking for next genome...\n\n";

		
		const type_map_string_to_uint_to_Sampled_diploid_Event_data::iterator it_gen = loopable_genome_to_sampled_Events.at(genrank);
	
		genome_name = it_gen->first;
		std::string dir_gen_name_search_str(it_gen->first);
		dir_gen_name_search_str.append("_"); // change by MATT SIMULATIONS
		std::string bam_gen_name_search_str(it_gen->first);
		bam_gen_name_search_str.append("."); // change by MATT SIMULATIONS
		
		std::string full_job_dirname;
		{//full_job_dirname
			const boost::filesystem::directory_iterator it_end;  
			for (boost::filesystem::directory_iterator it_jobdir(postdatatdir);
				it_jobdir != it_end;
				++it_jobdir)
			{		
				if (boost::filesystem::is_directory(it_jobdir->path()))
				{
					if (it_jobdir->path().filename().string().find(dir_gen_name_search_str) != std::string::npos)
					{
						full_job_dirname = boost::filesystem::complete(it_jobdir->path()).string();
						break;			
					}		    		    
				}//dir				
			}//it_jobdir	    	    	    
		}//full_job_dirname
		
		
		std::cerr << "\n\n\nrank " << my_MPI_rank << "\tfull_job_dirname = [" << full_job_dirname << "]\n";
		
		output_dir = full_job_dirname;
	
		
		std::string bam_filename;	
		{//full_job_dirname
			const boost::filesystem::directory_iterator it_end;  
			for (boost::filesystem::directory_iterator it_bamdir("/users/mmparks/data/mmparks/simulate"); // CHANGE BY MATT SIMULATION //("/gpfs/scratch/mmparks/SAMs_and_BAMs");
				it_bamdir != it_end;
				++it_bamdir)
			{		
				if (!boost::filesystem::is_directory(it_bamdir->path()))
				{
					if (it_bamdir->path().filename().string().find(bam_gen_name_search_str) != std::string::npos
					and  it_bamdir->path().filename().string().find(".bas") == std::string::npos
					and  it_bamdir->path().filename().string().find(".bai") == std::string::npos
					and  it_bamdir->path().filename().string().find(".bam") == it_bamdir->path().filename().string().size() - 4)			
					{
						bam_filename = boost::filesystem::complete(it_bamdir->path()).string();
						break;			
					}		    		    
				}//file
			}//it_jobdir	    	    	    
		}//full_job_dirname	
		
		std::cerr << "\n\n\nrank " << my_MPI_rank << "\ttbam_filename = [" << bam_filename << "]\n\n";
		
			
		
		assert(!full_job_dirname.empty());
		assert(!bam_filename.empty());			
		
		
		{
			char common_fname[option_size];
			std::sprintf(common_fname, "%s/raw_count_ratios", output_dir.c_str());
			boost::filesystem::remove(common_fname);			
		}
		
		
		BamTools::BamReader my_BAM_reader;

		prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population
						(bam_filename,
						std::string(),//gender_and_population_filename,
						my_BAM_reader);		


		read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda(bam_filename, false, NULL);
		std::cerr.flush();
		std::cerr << "\n\n\nrank " << my_MPI_rank << "\tprepare-done.";
		

	
		
		
		const std::vector<type_map_uint_to_Sampled_diploid_Event_data::iterator>  myev__loop(
						get_loopable_iterators<type_map_uint_to_Sampled_diploid_Event_data>(it_gen->second));
	
		
		
		#pragma omp parallel
		{
			uint thread_work_ctr = 0;  //for debug
			
			#pragma omp for schedule(static)
			for (int myev_ctr = 0; myev_ctr < myev__loop.size(); ++myev_ctr)
			{    				
				std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() 
						<< "  myev_ctr = " << myev_ctr << " out of " << myev__loop.size() << ".\n";
						
						
				std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "-work-counter = " << thread_work_ctr << "\n";
				++thread_work_ctr;
				
				const type_map_uint_to_Sampled_diploid_Event_data::iterator it_myev = myev__loop.at(myev_ctr);
												
				
// // 				//DEBUG!!!
// 				if (it_myev->second.event_UID != 17095 or it_myev->second.associated_genome.compare("NA12716") != 0)
// 					continue;
// // 				//DEBUG!!
						
							
				const type_map_uint_to_Event::iterator result_Ev_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID); 
				
				
				if (result_Ev_it->second.chromos[0] == 24  and   gender_of_individual_is_female)
				{
					std::cerr << "Event is on Y chromosome, but individual is female.  Skipping read-depth graph in in \"create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights\".\n";
					continue;
				}
				
				const BOOST_Interval expanded_region_of_interest( 
							safe_subtract_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.lower(), 100000),
							safe_chromo_add_base_1(result_Ev_it->second.region_between_and_including_the_LCRs_themselves.upper(), 100000, result_Ev_it->second.chromos[0]));        
							
				//identify all events intersecting this region:
				type_map_uint_to_Event_ptr  all_events_near_expanded_region;
				for (type_map_uint_to_Event::iterator it_ev = Conn_Comps.at(it_myev->second.cc_of_event).events.begin();
					it_ev != Conn_Comps.at(it_myev->second.cc_of_event).events.end();
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
				
				
				
				
				
				
				char out_RD_core_fname[option_size];
				
				std::sprintf(out_RD_core_fname, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s_%s___%s", 
							output_dir.c_str(),
							it_myev->second.cc_of_event, result_Ev_it->second.UID,
							convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first).c_str(),
							convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second).c_str(),
							genome_name.c_str());		    
						
	
				
				std::string observed_fname(out_RD_core_fname);
				observed_fname.append(".data.observed");	    
				
				if (boost::filesystem::exists(observed_fname) == 0)
				{//search
					std::string searchdir(output_dir);
					searchdir.append("/Read_Depth");
									
					std::stringstream cc_strm;
					cc_strm << "cc_" << it_myev->second.cc_of_event << "_";
					std::stringstream ev_strm;
					ev_strm << "ev_" << it_myev->second.event_UID << "_";
					
					char observed_fname_old[option_size];
					
					bool found_something = false;
					const boost::filesystem::directory_iterator it_end;  
					for (boost::filesystem::directory_iterator it_rdfile(searchdir);
						it_rdfile != it_end;
						++it_rdfile)
					{		
						if (!boost::filesystem::is_directory(it_rdfile->path()))
						{
							if (it_rdfile->path().filename().string().find(cc_strm.str()) != std::string::npos
								and  it_rdfile->path().filename().string().find(ev_strm.str()) != std::string::npos
								and  it_rdfile->path().filename().string().find(genome_name) != std::string::npos
								and it_rdfile->path().filename().string().find("RD_hypothesis_test__cc_") != std::string::npos)
							{
								char scanformat[option_size];
								std::sprintf(scanformat, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%%s", 
											output_dir.c_str(),
											it_myev->second.cc_of_event,
											it_myev->second.event_UID);	
								
								char in_calls_and_genome[option_size];		    
												
				// 			    std::cerr << "\n\t\t\t\ttry name: [" << boost::filesystem::complete(it_rdfile->path()).string().c_str() << "\n\t\t\t\tscanformat = [ " << scanformat << "]\n";
								std::sscanf(boost::filesystem::complete(it_rdfile->path()).string().c_str(), scanformat,
										in_calls_and_genome);
								
								std::string calls_CHOPPED(in_calls_and_genome);
								calls_CHOPPED = calls_CHOPPED.substr(0,calls_CHOPPED.find(genome_name));
								
				// 			    std::cerr << "\n\t\t\t\tcalls_CHOPPED = [" << calls_CHOPPED << "]\n";
								
								std::sprintf(observed_fname_old, "%s/Read_Depth/RD_hypothesis_test__cc_%u__ev_%u___%s%s.data.observed", 
									output_dir.c_str(),
									it_myev->second.cc_of_event, result_Ev_it->second.UID,
									calls_CHOPPED.c_str(),
									genome_name.c_str());	    
								
								found_something = true;
								break;			
							}		    		    
						}//file
					}//it_jobdir
					
					
					if (!found_something)
					{
						std::stringstream error_strm;
						error_strm << "ERROR  unable to find anything!!!\n\n";
						error_message(error_strm,false);
						continue;
					}	    
					
					std::cerr << "\n\t\t\tout_RD_core_fname = [" << out_RD_core_fname << "]\n"
						<<  "\t\t\tobserved_fname_old = [" << observed_fname_old << "\n\n";
											
					if (boost::filesystem::exists(observed_fname) == 0)
					{  
						std::cerr << "rename = [" << observed_fname << "\n";
						boost::filesystem::rename(observed_fname_old, observed_fname);		    
					}
					
				}//search

				
				const type_map_uint_to_longdouble map_Ref_Genome_position_to_per_read_posterior_sum(read_map_from_file<uint,longdouble>(observed_fname, false));
			
				
				
				{//RD hypothesis test
				
					BOOST_Interval affected_region(empty_BI);
					
					if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None))
					{
						affected_region = result_Ev_it->second.region_between_and_including_the_LCRs_themselves;	    
					}
					else
					{
						for (uint hap=0; hap<2; ++hap)
						{
							if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) != hap_outcome__None)
							{
								const type_BI__BI affected_regions_each_LCR(
											BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(
												convert_each_profile_breakpoint_to_absolute_breakpoints(
													pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
													result_Ev_it->second.compressed_map_LCR_to_profile__for_each_LCR)));
											
								const BOOST_Interval affected_hull(affected_regions_each_LCR.first.lower(), affected_regions_each_LCR.second.upper());
																							
								if (BOOST_empty(affected_region))
								{  affected_region = affected_hull;  }
								else
								{  affected_region = BOOST_hull(affected_region, affected_hull);  }
							}//occurs
						}//hap
					}
					
					

					//restrict to affected region:
					const type_map_uint_to_longdouble  map_affected_Ref_Genome_position_to_per_read_posterior_sum(
							get_map_intersection_of_map_keys_and_interval<uint,longdouble>(map_Ref_Genome_position_to_per_read_posterior_sum, affected_region));
								
							

// 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get affected_sequence__padded\n";
// 					std::cerr << "\ncalculating Wilcoxon_pvalue__inferred...\n";
					const uint pad_amount = 5000;
					const uint padded_max_chromo_sensitive = safe_chromo_add_base_1(affected_region.upper(), pad_amount, result_Ev_it->second.chromos[0]);
					//std::max<uint>(affected_region.upper() + pad_amount, );
					
					const std::string affected_sequence__padded(get_specific_DNA_sequence_from_Reference_genome(result_Ev_it->second.chromos[0],
														BOOST_Interval(affected_region.lower(), padded_max_chromo_sensitive)));
					
					
// 					type_vector_longdouble bernoulli_p_by_position(BOOST_width_inclusive(affected_region), 0);
					type_vector_vector_longdouble poisson_lambda_by_position_and_readgroup___NULL(BOOST_width_inclusive(affected_region),
																    type_vector_longdouble(Readgroup_stats.size(), 0));		
					
					type_vector_uint bernoulli_observations_by_position(BOOST_width_inclusive(affected_region), 0);
					
					boost::math::normal_distribution<mpfr_class> my_standard_normal(0,1);
					
// 					std::cerr  << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "get bernoulli_p_by_position\n";
					//get bernoulli p at each position by summing over all readgroups, using GC content at that position.
					for (uint le = 0; le < BOOST_width_inclusive(affected_region); ++le)
					{
						uint rg_indx = 0;
						for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
							it_rg != Readgroup_stats.end();
							++it_rg)
						{
							const uint gc_content = calculate_GC_content_of_sequence(affected_sequence__padded.substr(le, it_rg->second.median_insert_size));
							
							poisson_lambda_by_position_and_readgroup___NULL.at(le).at(rg_indx)
									= it_rg->second.haploid_fragmentation_rate__per_GC_count.at(gc_content) * 2; // null case is diploid
							
							++rg_indx;
						}
					}//bernoulli p
					
					
					
					type_vector_longdouble poisson_lambda_by_position__null_model(BOOST_width_inclusive(affected_region), 0);
					for (uint i=0; i < poisson_lambda_by_position__null_model.size(); ++i)
					{
						longdouble total_lambda_over_readgroups = 0;
						for (uint rg=0; rg < poisson_lambda_by_position_and_readgroup___NULL.at(i).size(); ++rg)
						{  total_lambda_over_readgroups += poisson_lambda_by_position_and_readgroup___NULL.at(i).at(rg);  }
						
						poisson_lambda_by_position__null_model.at(i) = total_lambda_over_readgroups;
					}//i
					
					
					
					const longdouble total_observed_number_of_reads___inferred = sum_over_map_values<uint, longdouble>(map_affected_Ref_Genome_position_to_per_read_posterior_sum);
					const longdouble expected_number_reads___NULL___inferred = sum_over_vector<longdouble>(poisson_lambda_by_position__null_model);
					
					const longdouble ratio_obs_over_exp_null___inferred = total_observed_number_of_reads___inferred / expected_number_reads___NULL___inferred;
					
					
					

// 					std::cerr << "get is_unique\n";	
					type_vector_bool is_unique(BOOST_width_inclusive(affected_region), true);
					{
						const type_map_uint_to_list_BI::const_iterator it_found_chr = nonuniqe_regions_of_genome.find(result_Ev_it->second.chromos[0]); 
						if (it_found_chr != nonuniqe_regions_of_genome.end())
						{
							for (type_list_BI::const_iterator it_bi = it_found_chr->second.begin();
								it_bi != it_found_chr->second.end();
								++it_bi)
							{
								const BOOST_Interval intersection(BOOST_intersect(*it_bi, affected_region));
								if (!BOOST_empty(intersection))
								{
									for (uint pos = intersection.lower(); pos <= intersection.upper(); ++pos)
									{  is_unique.at(pos - affected_region.lower()) = false;  }						
								}														
							}
						}
					}
					
					
					
					longdouble  expected_number_reads___NULL___unique = 0;
					for (uint indx = 0; indx < BOOST_width_inclusive(affected_region); ++indx)
					{ 
						if (is_unique.at(indx))
						{
							expected_number_reads___NULL___unique += poisson_lambda_by_position__null_model.at(indx);
						}
					}
					
					longdouble total_observed_number_of_reads___unique = 0;
					for (type_map_uint_to_longdouble::const_iterator it_pos = map_affected_Ref_Genome_position_to_per_read_posterior_sum.begin();
						it_pos != map_affected_Ref_Genome_position_to_per_read_posterior_sum.end();
						++it_pos)
					{
						if (is_unique.at(it_pos->first - affected_region.lower()))
						{
							total_observed_number_of_reads___unique += it_pos->second;
						}
					}
					
					
					const longdouble ratio_obs_over_exp_null___unique = total_observed_number_of_reads___unique / expected_number_reads___NULL___unique;
					
					
					std::stringstream  out_ss;
				
					out_ss  << "\nrank " << my_MPI_rank << ", thread " << omp_get_thread_num() 
							<< "\n\tmyev_ctr = " << myev_ctr
							<< "\n\tgenome = " << genome_name
							<< "\n\tcc = " << it_myev->second.cc_of_event
							<< "\n\tevent = " << it_myev->second.event_UID
							<< "\n\taffected_region = " << affected_region 
							<< "\n\tchrome = " <<  result_Ev_it->second.chromos[0]
							<< "\n\tchrom-length = " << Event::chromosome_lengths.at(result_Ev_it->second.chromos[0]) 
							<< "\n\toutcome = " << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first) 
							<< ", " << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second) 
							<< "\n\tratio-inferred(obs/null) = " << ratio_obs_over_exp_null___inferred
							<< "\n\tratio-unique(obs/null) = " << ratio_obs_over_exp_null___unique
							<< "\n\n";
							
					std::cerr << out_ss.str();
					
					
				
					#pragma omp critical(write_wilcoxon_signed_rank_to_common_file)
					{//common file		
						//write to common filename
						char common_fname[option_size];
						std::sprintf(common_fname, "%s/raw_count_ratios", output_dir.c_str());
						
						std::ofstream out_wilcoxon_fs(common_fname, std::ios_base::app);
						
						
						out_wilcoxon_fs 
							//info:
							<< genome_name << "\t"
							<< it_myev->second.cc_of_event << "\t"
							<< it_myev->second.event_UID << "\t"
							
							//p-values:
							<< ratio_obs_over_exp_null___inferred << "\t"
							<< ratio_obs_over_exp_null___unique << "\n";
							
						out_wilcoxon_fs.close();
					}//common file
					
// 					std::cerr  << "\n\trank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "writing_is_DONE.\n";
				
				}//RD hypothesis test			
			}//it_myev
			
			std::cerr << "rank " << my_MPI_rank << ", thread " << omp_get_thread_num() << "  done with it_myev loop...\n";
			
		}//parallel
		std::cerr << "\nrank " << my_MPI_rank << " moving to next genome...\n";
		
    }//it_gen            		
	
}//raw_count_ratios
































void determine_cumulative_GC_count_for_each_chromosome_and_save_to_file()
{
	
	#pragma omp parallel for schedule(dynamic,1)
	for (uint chromo = 1; chromo <= 24; ++chromo)
	{
		//upload entire chromosoome:
		const std::string chromo_seq(get_specific_DNA_sequence_from_Reference_genome(chromo, BOOST_Interval(1, Event::chromosome_lengths.at(chromo))));
		
		type_vector_uint cumulative_GC_count_this_chromo(chromo_seq.size(), 0);
		uint running_GC_count = 0;
		for (uint pos = 0; pos < chromo_seq.size(); ++pos)
		{
			if (chromo_seq.at(pos) == 'C'  or  chromo_seq.at(pos) == 'G')
			{
				++running_GC_count;
			}
			
			cumulative_GC_count_this_chromo.at(pos) = running_GC_count;
		}//pos
		
		char outfnamechar[option_size];
		std::sprintf(outfnamechar, "%s/cumulative_GC_count_chromosome_%u", data_directory.c_str(), chromo);
		std::ofstream outfs(outfnamechar);
		boost::archive::text_oarchive outarch(outfs);
		
		outarch << cumulative_GC_count_this_chromo;
		
		outfs.close();
	}//chromo
	
}//determine_cumulative_GC_count_for_each_chromosome_and_save_to_file

































const type_uint__uint convert_random_draw_to_chromo_and_position
		(const uint &whole_genome_position_draw,
		 const type_map_BI_to_uint &map_cumulative_genome_position_to_chromo)
{
	for (type_map_BI_to_uint::const_iterator it_bi = map_cumulative_genome_position_to_chromo.begin();
		it_bi != map_cumulative_genome_position_to_chromo.end();
		++it_bi)
	{
		if (BOOST_in(whole_genome_position_draw, it_bi->first))
		{
			return  type_uint__uint(it_bi->second, whole_genome_position_draw - it_bi->first.lower() + 1);
		}		
	}//bi
	
	std::stringstream error_strm;
	error_strm << "unable to find whole_genome_position_draw = " << whole_genome_position_draw << "\n";
	print_map_keys_and_values<BOOST_Interval, uint>(map_cumulative_genome_position_to_chromo, "map_cumulative_genome_position_to_chromo", &error_strm, true);
	error_message(error_strm, false);
	
	return type_uint__uint(0,0);
	
}//convert_random_draw_to_chromo_and_position















bool check_if_string_contains_N
		(const std::string &some_seq)
{
	for (uint j=0; j < some_seq.length(); ++j)
	{
		if (some_seq.at(j) == 'N')
		{  return true;  }
	}
	
	return false;
	
}//check_if_string_contains_N









void  randomly_sample_contiguous_region_from_control_genome
				(const type_map_uint_to_list_BI &nonunique_regions_of_genome,
				const type_map_uint_to_set_uint &learning_positions_by_chromo,
				const type_map_BI_to_uint &map_cumulative_genome_position_to_chromo,
				const uint &desired_region_length,
				const type_map_uint_to_string &map_chromo_to_seq,
				type_map_uint_to_set_BI &sampled_test_regions)
{	
	
	const uint nonsex_chromo_len = (--map_cumulative_genome_position_to_chromo.end())->first.upper();
	static boost::random::mt19937 random_generator;  // "static" BY MATT for simluation
	boost::uniform_int<uint> uniform_genome_dist(1, nonsex_chromo_len);
	const uint skip_amount_at_end_of_chromos = 1000;
	
	const uint max_num_attempts = 100000;
	uint attempt_ctr = 0;
	while (attempt_ctr < max_num_attempts)
	{
		++attempt_ctr;
		if (attempt_ctr % 500 == 0)
		{  std::cerr << "\t\tattempt_ctr = " << attempt_ctr << "\n";  }
		
		const type_uint__uint chromo_and_pos(convert_random_draw_to_chromo_and_position(uniform_genome_dist(random_generator), map_cumulative_genome_position_to_chromo));
		const type_uint__BI chromo_and_reg(chromo_and_pos.first, BOOST_Interval(chromo_and_pos.second, chromo_and_pos.second + desired_region_length - 1));
	
		const BOOST_Interval acceptable_p_arm(skip_amount_at_end_of_chromos, Event::centromere_coordinates.at(chromo_and_pos.first).first - skip_amount_at_end_of_chromos);
		const BOOST_Interval acceptable_q_arm(Event::centromere_coordinates.at(chromo_and_pos.first).second + skip_amount_at_end_of_chromos,
											Event::chromosome_lengths.at(chromo_and_pos.first) - skip_amount_at_end_of_chromos);
		
		if (!BOOST_subset(chromo_and_reg.second, acceptable_p_arm)  and  !BOOST_subset(chromo_and_reg.second, acceptable_q_arm))
		{  continue;  }		
		
		if (!check_if_region_overlaps_NON_unique_region(nonunique_regions_of_genome, chromo_and_pos.first, chromo_and_pos.second))
		{
// 			if (learning_positions_by_chromo.count(chromo_and_pos.first) == 0 
// 				or  get_intersection_of_set_and_interval(learning_positions_by_chromo.at(chromo_and_pos.first), chromo_and_reg.second).empty())
// 			{
// 				if (sampled_test_regions.count(chromo_and_reg.first) == 0)
// 				{
					if (!check_if_string_contains_N(map_chromo_to_seq.at(chromo_and_reg.first)
										.substr(chromo_and_reg.second.lower()-1,
											BOOST_width_inclusive(chromo_and_reg.second))))
					{
						sampled_test_regions[chromo_and_reg.first].insert(chromo_and_reg.second);
						return;
					}
// 				}
// 				else
// 				{
// 					bool overlaps_something = false;
// 					for (type_set_BI::const_iterator it_bi = sampled_test_regions.at(chromo_and_reg.first).begin();
// 						it_bi != sampled_test_regions.at(chromo_and_reg.first).end();
// 						++it_bi)
// 					{
// 						if (BOOST_overlap(*it_bi, chromo_and_reg.second))
// 						{
// 							overlaps_something = true;
// 							break;
// 						}
// 						else if (chromo_and_reg.second < it_bi->lower())
// 						{  break;  }
// 					}//bi
// 					
// 					if (!overlaps_something)
// 					{
// 						if (!check_if_string_contains_N(map_chromo_to_seq.at(chromo_and_reg.first)
// 											.substr(chromo_and_reg.second.lower()-1,
// 												BOOST_width_inclusive(chromo_and_reg.second))))	
// 						{
// 							sampled_test_regions[chromo_and_pos.first].insert(chromo_and_reg.second);
// 							return;
// 						}
// 					}
// 				}//has chromo
// 			}
		}
	}//attempt_ctr
	
	
	std::stringstream error_strm;
	error_strm << "too many attempts in \"randomly_sample_contiguous_region_from_control_genome\n";
	print_map_keys_and_values<BOOST_Interval, uint>(map_cumulative_genome_position_to_chromo, "map_cumulative_genome_position_to_chromo", &error_strm, true);
	error_message(error_strm, false);
	
}//randomly_sample_contiguous_region_from_control_genome















void randomly_sample_GC_sensitive_positions_from_control_genome
			(const uint &number_positions_for_sample,
			const type_map_uint_to_list_BI &nonunique_regions_of_genome,
			const type_map_uint_to_set_uint &learning_positions_by_chromo,
			const type_map_BI_to_uint &map_cumulative_genome_position_to_chromo,
			const type_map_uint_to_uint &map_frag_lengths_to_desired_GC,
			const type_map_uint_to_string &map_chromo_to_seq,
			type_map_uint_to_set_uint &sampled_test_positions)
{	
	const uint nonsex_chromo_len = (--map_cumulative_genome_position_to_chromo.end())->first.upper();
	boost::random::mt19937 random_generator;
	boost::uniform_int<uint> uniform_genome_dist(1, nonsex_chromo_len);
	const uint skip_amount_at_end_of_chromos = 1000;
	
	uint attempt_ctr = 0;
	uint samp_ctr = 0;
	while (samp_ctr < number_positions_for_sample)
	{
		++attempt_ctr;
		if (attempt_ctr % 100000 == 0)
		{  std::cerr << "\t\t\t\tattempt_ctr = " << attempt_ctr << "\n";  }		
		
// 		if (samp_ctr % 500 == 0)
// 		{  std::cerr << "\t\t\t\tsamp_ctr = " << samp_ctr << "\n";  }
		
		const type_uint__uint chromo_and_pos(convert_random_draw_to_chromo_and_position(uniform_genome_dist(random_generator), map_cumulative_genome_position_to_chromo));
	
		const BOOST_Interval acceptable_p_arm(skip_amount_at_end_of_chromos, Event::centromere_coordinates.at(chromo_and_pos.first).first - skip_amount_at_end_of_chromos);
		const BOOST_Interval acceptable_q_arm(Event::centromere_coordinates.at(chromo_and_pos.first).second + skip_amount_at_end_of_chromos,
											Event::chromosome_lengths.at(chromo_and_pos.first) - skip_amount_at_end_of_chromos);
		
		if (BOOST_in(chromo_and_pos.second, acceptable_p_arm)  or  BOOST_in(chromo_and_pos.second, acceptable_q_arm))
		{ 
			if (!check_if_position_is_in_NON_unique_region(nonunique_regions_of_genome, chromo_and_pos.first, chromo_and_pos.second))
			{
				if (learning_positions_by_chromo.count(chromo_and_pos.first) == 0 
					or  learning_positions_by_chromo.at(chromo_and_pos.first).count(chromo_and_pos.second) == 0)
				{
					if (sampled_test_positions.count(chromo_and_pos.first) == 0  or  sampled_test_positions.at(chromo_and_pos.first).count(chromo_and_pos.second) == 0)
					{
						//check GC distrib
						bool satisfies_all_GC_requirements = true;
						for (type_map_uint_to_uint::const_iterator it_frag = map_frag_lengths_to_desired_GC.begin();
							it_frag != map_frag_lengths_to_desired_GC.end();
							++it_frag)
						{
							const uint lower_pos_GC = map_chromosome_to_cumulative_GC_counts.at(chromo_and_pos.first).at(chromo_and_pos.second - 1);
							const uint upper_pos_GC = map_chromosome_to_cumulative_GC_counts.at(chromo_and_pos.first).at(chromo_and_pos.second + it_frag->first - 1);
							
							if (upper_pos_GC - lower_pos_GC != it_frag->second)
							{
								satisfies_all_GC_requirements = false;
								break;
							}
						}//frag
						
						if (satisfies_all_GC_requirements)
						{
							if (!check_if_string_contains_N(map_chromo_to_seq.at(chromo_and_pos.first)
											.substr(chromo_and_pos.second-1, maximum_average_fragment_length + 500)))
							{
								sampled_test_positions[chromo_and_pos.first].insert(chromo_and_pos.second);
								++samp_ctr;
							}
						}
					}// sampled
				}// learning
			}// unique
		}// p-/q-arms
	}//samp_ctr
	
}//randomly_sample_positions_from_control_genome












void  draw_count_ratio_using_sampling_strategy_A
			(const type_map_uint_to_CC &Conn_Comps,
			BamTools::BamReader &a_BAM_reader)
{
	std::cerr << "\n\ninside \"draw_count_ratio_using_sampling_strategy_A\"...\n\n";
	
	const type_map_uint_to_list_BI nonunique_regions_of_genome(load_Non_unique_regions_of_genome());
	type_map_uint_to_set_uint  learning_positions_by_chromo; //empty
	const type_map_BI_to_uint  map_cumulative_genome_position_to_chromo(make_map_cumulative_genome_position_to_chromosome());
	
	std::cerr << "\nloading entire genome into strings...\n";
	type_map_uint_to_string map_chromo_to_seq;	
	for (uint chr_indx = 1; chr_indx < 23; ++chr_indx)
	{
		map_chromo_to_seq[chr_indx] = get_specific_DNA_sequence_from_Reference_genome(chr_indx, BOOST_Interval(1, Event::chromosome_lengths.at(chr_indx)));				
	}//chr_indx	
	
	
	std::cerr << "\tsampling test regions...\n";
	type_map_uint_to_set_BI sampled_test_regions;
	
	for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{
		std::cerr << "\t\tcc = " << it_cc->first << "\n";
		for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
			it_ev != it_cc->second.events.end();
			++it_ev)
		{
			std::cerr << "\t\t\tev = " << it_ev->first << "\n";
			
			if (it_ev->second.chromos[0] < 23) //autosomal only
			{
				for (uint rep=0; rep < 10; ++rep)
				{
					randomly_sample_contiguous_region_from_control_genome(
						nonunique_regions_of_genome,
						learning_positions_by_chromo,
						map_cumulative_genome_position_to_chromo,
						BOOST_width_inclusive(it_ev->second.region_between_and_including_the_LCRs_themselves),
						map_chromo_to_seq,
						sampled_test_regions);
				}
			}
		}
	}
	
	
	std::cerr << "\n\tmaking count ratios...\n";

	typedef  std::pair<type_uint__BI,longdouble>  type_uint__BI__longdouble;
	typedef  std::list< type_uint__BI__longdouble>   type_list__uint__BI__longdouble;
	type_list__uint__BI__longdouble  count_ratios;
	
	for (type_map_uint_to_set_BI::const_iterator it_chromo = sampled_test_regions.begin();
		it_chromo != sampled_test_regions.end();
		++it_chromo)
	{
		std::cerr << "\t\tchromo = " << it_chromo->first << "\n";
		for (type_set_BI::const_iterator it_bi = it_chromo->second.begin();
			it_bi != it_chromo->second.end();
			++it_bi)
		{
			std::cerr << "\t\t\tregion = " << *it_bi << "     width = " << BOOST_width_inclusive(*it_bi) << "\n";
			
			
			const uint observation_counts = count_decent_reads_with_left_endpoints_in_region__readgroup_specific(
										a_BAM_reader, it_chromo->first, *it_bi);
			std::cerr << "\t\t\t\tcounts = " << observation_counts << "\n";
			
			
			
			longdouble total_expected_frag_rate = 0;
			
			for (uint pos = it_bi->lower(); pos <= it_bi->upper(); ++pos)
			{
				for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
					it_rg != Readgroup_stats.end();
					++it_rg)
				{
					const uint lower_GC = map_chromosome_to_cumulative_GC_counts.at(it_chromo->first).at(pos-1);
					const uint upper_GC = map_chromosome_to_cumulative_GC_counts.at(it_chromo->first).at(pos + it_rg->second.median_insert_size -1);
					total_expected_frag_rate += it_rg->second.haploid_fragmentation_rate__per_GC_count.at(upper_GC - lower_GC);
				}
			}//pos
			
			std::cerr << "\t\t\t\texpected = " << total_expected_frag_rate << "\n";
			
			count_ratios.push_back(type_uint__BI__longdouble(type_uint__BI(it_chromo->first, *it_bi),  (longdouble)observation_counts / (total_expected_frag_rate*2)));
		}//bi
	}//chromo
	
	
	
	
	{//save
		std::cerr << "\tsave...\n";
		std::string savename(scratch_job_output_dir);
		savename.append("/");
		savename.append(genome_name);
		savename.append(".count_ratios_A");
		std::ofstream out_counts_fs(savename);
		
		for (type_list__uint__BI__longdouble::const_iterator it_control = count_ratios.begin();
			it_control != count_ratios.end();
			++it_control)
		{
			out_counts_fs
				<< it_control->second << "\t"
				<< it_control->first.first << "\t"
				<< it_control->first.second.lower() << "\t"
				<< it_control->first.second.upper() << "\t"
				<< genome_name << "\n";
		}//control
		
		out_counts_fs.close();
	}//save
	
	
	
	std::cerr << "\n\nDONE with \"draw_count_ratio_using_sampling_strategy_A\"!\n\n";
	
}//draw_count_ratio_using_sampling_strategy_A



























void  draw_count_ratio_using_sampling_strategy_B
			(const type_map_uint_to_CC &Conn_Comps,
			BamTools::BamReader &a_BAM_reader)
{
	std::cerr << "\n\ninside \"draw_count_ratio_using_sampling_strategy_B\"...\n\n";
	
	const type_map_uint_to_list_BI nonunique_regions_of_genome(load_Non_unique_regions_of_genome());
	type_map_uint_to_set_uint  learning_positions_by_chromo; //empty
	const type_map_BI_to_uint  map_cumulative_genome_position_to_chromo(make_map_cumulative_genome_position_to_chromosome());
	
	std::cerr << "\nloading entire genome into strings...\n";
	type_map_uint_to_string map_chromo_to_seq;	
	for (uint chr_indx = 1; chr_indx < 23; ++chr_indx)
	{
		map_chromo_to_seq[chr_indx] = get_specific_DNA_sequence_from_Reference_genome(chr_indx, BOOST_Interval(1, Event::chromosome_lengths.at(chr_indx)));				
	}//chr_indx		
	
	
	type_list_longdouble count_ratios;
	
	for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{
		std::cerr << "\tcc = " << it_cc->first << "\n";
		for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
			it_ev != it_cc->second.events.end();
			++it_ev)
		{
			std::cerr << "\t\tev = " << it_ev->first << "\n";
			
			if (it_ev->second.chromos[0] < 23) //autosomal only
			{
				type_map_uint_to_set_uint sampled_test_positions;
				
				std::cerr << "\t\t\tsample\n";
				
				//get sample positions
				
				#pragma ompa parallel
				{
					type_map_uint_to_set_uint sampled_test_positions__thread;
					
					for (uint pos = it_ev->second.region_between_and_including_the_LCRs_themselves.lower(); 
						pos <= it_ev->second.region_between_and_including_the_LCRs_themselves.upper(); ++pos)
					{
						type_map_uint_to_uint map_frag_len_to_GC_distrib;
						for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
							it_rg != Readgroup_stats.end();
							++it_rg)
						{
							const uint lower_GC = map_chromosome_to_cumulative_GC_counts.at(it_ev->second.chromos[0]).at(pos-1);
							const uint upper_GC = map_chromosome_to_cumulative_GC_counts.at(it_ev->second.chromos[0]).at(pos + it_rg->second.median_insert_size -1);
							map_frag_len_to_GC_distrib[it_rg->second.median_insert_size] = upper_GC - lower_GC;
						}//rg
						
						randomly_sample_GC_sensitive_positions_from_control_genome
							(1,
							nonunique_regions_of_genome,
							learning_positions_by_chromo,
							map_cumulative_genome_position_to_chromo,
							map_frag_len_to_GC_distrib,
							map_chromo_to_seq,
							sampled_test_positions__thread);
					}//pos	
					
					#pragma omp critical (add_sampled_test_position_from_each_thrread)
					{
						for (type_map_uint_to_set_uint::const_iterator it_sc = sampled_test_positions__thread.begin();
							it_sc != sampled_test_positions__thread.end();
							++it_sc)
						{
							sampled_test_positions.insert(std::pair<uint, type_set_uint>(it_sc->first, type_set_uint()));
							type_map_uint_to_set_uint::iterator it_found_chr = sampled_test_positions.find(it_sc->first);
							type_set_uint::iterator insert_it = it_found_chr->second.begin();
							for (type_set_uint::const_iterator it_sp = it_sc->second.begin();
								it_sp != it_sc->second.end();	
								++it_sp)
							{
								insert_it = it_found_chr->second.insert(insert_it, *it_sp);
							}//sp
						}//chromo
					}//critical
					
				}//parallel
				
				std::cerr << "\t\t\tcount\n";
				
				
				//count reads at these positions
				longdouble total_expected_frag_rate = 0;
				uint total_observation_counts = 0;
				
				for (type_map_uint_to_set_uint::const_iterator it_chromo = sampled_test_positions.begin();
					it_chromo != sampled_test_positions.end();
					++it_chromo)
				{
					for (type_set_uint::const_iterator it_pos = it_chromo->second.begin();
						it_pos != it_chromo->second.end();
						++it_pos)
					{
						total_observation_counts += count_decent_reads_with_left_endpoints_in_region__readgroup_specific(
													a_BAM_reader, it_chromo->first, BOOST_make_point(*it_pos));
						
						for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
							it_rg != Readgroup_stats.end();
							++it_rg)
						{
							const uint lower_GC = map_chromosome_to_cumulative_GC_counts.at(it_chromo->first).at(*it_pos-1);
							const uint upper_GC = map_chromosome_to_cumulative_GC_counts.at(it_chromo->first).at(*it_pos + it_rg->second.median_insert_size -1);
							total_expected_frag_rate += it_rg->second.haploid_fragmentation_rate__per_GC_count.at(upper_GC - lower_GC);
						}//rg
					}//bi
				}//chromo
				
				count_ratios.push_back((longdouble)total_observation_counts / (total_expected_frag_rate*2));
			}//non-sex chromo
		}//ev
	}//cc
	
	std::cerr << "save\n";
	std::string savedir(scratch_job_output_dir);
	savedir.append("/count_ratios_B");	
	write_singlecontainer_to_file<type_list_longdouble>(count_ratios, "count_ratios_B");
	
	std::cerr << "\n\nDONE with \"draw_count_ratio_using_sampling_strategy_B\"!\n\n";
		
}//draw_count_ratio_using_sampling_strategy_B






















type_map_BI_to_uint  make_map_cumulative_genome_position_to_chromosome()
{
	type_map_uint_to_uint non_sex_chromosome_lengths(Event::chromosome_lengths);
	non_sex_chromosome_lengths.erase(23);
	non_sex_chromosome_lengths.erase(24);
	
	std::cerr << "making cumulative genome map to chromosome...\n";
	uint length_of_haploid_nonsex_genome = 0;
	type_map_BI_to_uint map_cumulative_genome_position_to_chromo;
	{
		uint cumulative_index = 0;
		for (type_map_uint_to_uint::const_iterator it_chr = non_sex_chromosome_lengths.begin();
			it_chr != non_sex_chromosome_lengths.end();
			++it_chr)
		{
			map_cumulative_genome_position_to_chromo[BOOST_Interval(cumulative_index + 1, cumulative_index + it_chr->second)] = it_chr->first;
			cumulative_index += it_chr->second;
			length_of_haploid_nonsex_genome += it_chr->second;
		}//chromosome
	}
	
	return  map_cumulative_genome_position_to_chromo;
	
}//make_map_cumulative_genome_position_to_chromosome


















uint count_decent_reads_with_left_endpoints_in_region__readgroup_specific
					(BamTools::BamReader &a_BAM_reader,
					const uint &chromosome,
					const BOOST_Interval &chromo_interval)
{
	uint observation_count = 0;
	
	const bool success_set = a_BAM_reader.SetRegion(map_chromosome_value_to_BAM_Ref_IDs.at(chromosome), chromo_interval.lower() -2,
							map_chromosome_value_to_BAM_Ref_IDs.at(chromosome), chromo_interval.upper() + 1);
		
	BamTools::BamAlignment temp_core_Balgmt;
	while (a_BAM_reader.GetNextAlignmentCore(temp_core_Balgmt))    
	{
		if ( !temp_core_Balgmt.IsFailedQC())
		{
			if (temp_core_Balgmt.IsProperPair())
			{
				if (BOOST_in((uint)temp_core_Balgmt.Position+1, chromo_interval)  //genome is 1-based, bam-files are 0-based
					and temp_core_Balgmt.Position < temp_core_Balgmt.MatePosition)  //i.e. is first mate
				{  
					++observation_count;
				}
			}
			else if (temp_core_Balgmt.IsPaired())    
			{
				if (temp_core_Balgmt.IsMapped() and  BOOST_in((uint)temp_core_Balgmt.Position+1, chromo_interval))
				{  
					++observation_count;
				}
			}
		}//qc
	}//while
	//Recall that a read is counted only if it STARTS at the considered genome position (i.e. left-endpoint of read lands on that position).  Other reads that overlap are ignored. 
	
	return  observation_count;

}//count_decent_reads_with_left_endpoints_in_set_region








void Count_observation_data::absorb_counts
			(const type_map_string_to_uint &some_count_observations)
{
	for (type_map_string_to_uint::const_iterator it_rg = some_count_observations.begin();
		it_rg != some_count_observations.end();
		++it_rg)
	{
		map_readgroup_to_count[it_rg->first] += it_rg->second;
	}

}//absorb_counts











std::string infer_genome_name_from_BAM_filename
				(const std::string &some_bam_filename)
{
        const int last_slash = some_bam_filename.find_last_of('/');
        //const int dot_at_end_of_genome_name = some_bam_filename.find_first_of('.', last_slash); // BY MATT - original
        const int dot_at_end_of_genome_name = some_bam_filename.find_first_of('.', last_slash + 9);  // BY MATT - simulation
        if (last_slash+1 < some_bam_filename.size() 
		and   0 < dot_at_end_of_genome_name
		and   dot_at_end_of_genome_name - last_slash - 1 > 0)
	{
		return   some_bam_filename.substr(last_slash+1, dot_at_end_of_genome_name - last_slash - 1);
	}
	else 
	{
		return std::string();
	}
	
}//infer_genome_name_from_BAM_filename














uint get_CC_of_EV
			(const type_map_uint_to_CC  &Conn_Comps,
			const uint &ev_id)
{
	for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
		 it_cc != Conn_Comps.end();
		++it_cc)
	{
		if (it_cc->second.events.count(ev_id) > 0)
		{  return it_cc->first;  }
	}
	
	return 0;
	
} // get_CC_of_EV











void simulate_NAHR_events_on_Reference_genome
				(const type_map_uint_to_CC  &Conn_Comps,
				const std::string &sim_gen_name,
				const type_vector_uint &allowable_events,
				const uint &num_events_to_simulate_per_genome)
{
	//only do haploid events for ease and for simplicity of benchmarking across algorithm sake.
	std::cerr << "\n\ninside \"simulate_NAHR_events_on_Reference_genome\"...\n";
	
	static boost::random::mt19937 random_generator_sim;	
	boost::random::uniform_int_distribution<> event_dist(1, allowable_events.size());
	
	type_map_uint_to_Sampled_diploid_Event_data sampled_sim_Evs;			
	type_map_uint_to_set_uint samples; // easier bookkeeping.  cc --> ev
	
	uint number_of_successful_events_drawn = 0;
	while (number_of_successful_events_drawn < num_events_to_simulate_per_genome)
	{
		const uint drawn_ev = allowable_events.at(event_dist(random_generator_sim) - 1);
		std::cerr << "\tdrawn_ev = " << drawn_ev << "\n";
		const uint drawn_cc = get_CC_of_EV(Conn_Comps, drawn_ev);					
		std::cerr << "\tdrawn_cc = " << drawn_cc << "\n";
		
		
		//check that this event does not violate an existing exclusivity constraint.  also check that it has not been drawn before.
		const type_set_uint exclusivity_events(extract_keys_of_map_and_return_as_set<uint, Interlocking_entry<Event> >
												(Conn_Comps.at(drawn_cc).events.at(drawn_ev).local_interlocking_Events));
		
		type_set_uint prev_drawn_relevant_evs;
		if (samples.count(drawn_cc) > 0)
		{
			prev_drawn_relevant_evs = samples.at(drawn_cc);
		}			
		
		
		if (!get_set_intersection_of_two_sets<uint>(prev_drawn_relevant_evs, exclusivity_events).empty()
			or  prev_drawn_relevant_evs.count(drawn_ev) > 0)
		{
			continue;
		}
		
		
		std::cerr << "if made it this far...\n";
		//if made it this far, then it is safe to add this event to the sample.
		
		//draw breakpoint.  first draw a position along the proifle length.  then find the nearest VP (to the right).
		const type_map_uint_to_Event::const_iterator it_drawn_ev = Conn_Comps.at(drawn_cc).events.find(drawn_ev);
		boost::random::uniform_int_distribution<> brkpt_dist(1, it_drawn_ev->second.profile_length);
		
		type_set_uint::const_iterator it_nearest_VP = it_drawn_ev->second.variational_positions_of_profile.end();
		while (it_nearest_VP == it_drawn_ev->second.variational_positions_of_profile.end())
		{
			const uint prof_brkpt = brkpt_dist(random_generator_sim);				
			it_nearest_VP = it_drawn_ev->second.variational_positions_of_profile.lower_bound(prof_brkpt);
		}
		
		
		// draw del or dup.
		boost::random::uniform_int_distribution<> deldup_dist(1, 2);
		const uint drawn_del_or_dup = deldup_dist(random_generator_sim);
		
		type_haploid_outcome drawn_outcome = hap_outcome__Del;
		if (drawn_del_or_dup == 2)
		{
			drawn_outcome = hap_outcome__Dup;
		}

		
		//save sampled event			
		samples[drawn_cc].insert(drawn_ev);
		
		const type_map_uint_to_Sampled_diploid_Event_data::iterator new_sample
								= sampled_sim_Evs.insert(std::pair<uint, Sampled_diploid_Event_data>
															( drawn_ev, Sampled_diploid_Event_data())).first;
		new_sample->second.cc_of_event = drawn_cc;
		new_sample->second.event_UID = drawn_ev;
		new_sample->second.the_diploid_Event_outcome.first = drawn_outcome;
		new_sample->second.the_diploid_profile_brkpts.first = BOOST_make_point(*it_nearest_VP);
		
		++number_of_successful_events_drawn;										
	}//draw events
	
	print_map_to_set<uint, uint>(samples, "samples cc  --> evs");
	
	
	//sort
	type_list_Changed_region affected_list;
	for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_samp = sampled_sim_Evs.begin();
		it_samp != sampled_sim_Evs.end();
		++it_samp)
	{		
		const type_map_uint_to_Event::const_iterator it_the_ev = Conn_Comps.at(it_samp->second.cc_of_event).events.find(it_samp->first);
		
		Changed_region affected_region;
		affected_region.chr = it_the_ev->second.chromos[0];
		affected_region.region = convert_profile_breakpoint_to_absolute_breakpoints(it_samp->second.the_diploid_profile_brkpts.first.lower(),
																		it_the_ev->second.compressed_map_LCR_to_profile__for_each_LCR);
		
		affected_region.rearrangement = it_samp->second.the_diploid_Event_outcome.first;
		affected_region.cc = it_samp->second.cc_of_event;
		affected_region.ev = it_samp->first;
		
		affected_list.push_back(affected_region);		
	}//sampled_sim_Evs
	
	affected_list.sort(compare_Changed_regions);
	
	
	std::cerr << "writing rearranged genome.\n";
	
	//now write the rearranged version of the Reference
	std::ofstream outfs(sim_gen_name);
	
	type_list_Changed_region::const_iterator it_change = affected_list.begin();
	
	for (uint chr = 1; chr <= 24; ++chr)
	{
		std::cerr << "\tchr = " << chr << "\n";
		uint first_pos_to_be_written = 1;
		
		while (it_change != affected_list.end()  and  it_change->chr == chr)
		{
			std::cerr << "\t\t" << it_change->region 
					<< ",  dup =?= "  << (it_change->rearrangement == hap_outcome__Dup) 
					<< ", cc = " << it_change->cc
					<< ", ev = " << it_change->ev
					<< "\n";
					
			//write from the last pos to just before the rearrangement.
			const std::string normal_seq(get_specific_DNA_sequence_from_Reference_genome
												(chr, BOOST_Interval(first_pos_to_be_written, it_change->region.lower() - 1)));
			outfs << normal_seq;
			
			//apply rearrangement;
			if (it_change->rearrangement == hap_outcome__Dup)
			{
				BOOST_Interval reg(it_change->region.lower(), it_change->region.upper()-1);
				const std::string seq(get_specific_DNA_sequence_from_Reference_genome(chr, reg));
				std::cerr << "\t\t\t\treg width = " << BOOST_width_inclusive(reg) << ",  length(seq) = " << seq.length() << "\n";
				outfs << seq;
				outfs << seq; //twice, because is Dup.
				
			}

			first_pos_to_be_written = it_change->region.upper();
		
	
			++it_change;
		}
		
		std::cerr << "\t\trest of chromo...\n";
		
		//rest of chromo.
		const std::string normal_seq(get_specific_DNA_sequence_from_Reference_genome
											(chr, BOOST_Interval(first_pos_to_be_written, Event::chromosome_lengths.at(chr))));		
		
		outfs << normal_seq;						
	}//chr	
	
// unnecessary with pIRS.   if you use wgsim, then this is a good idea.	
// 	// Then write another normal copy of the reference.  (make it diploid genome, with haploid events)
// 	for (uint chr = 1; chr <= 24; ++chr)
// 	{
// 		const std::string normal_seq(get_specific_DNA_sequence_from_Reference_genome
// 											(chr, BOOST_Interval(1, Event::chromosome_lengths.at(chr))));		
// 		
// 		outfs << normal_seq;				
// 	}//chr
			
	outfs.close();
	
	
	
	std::cerr << "save truth.\n";
	
	//save truth	
	
	//write raw coords	
	std::string outfname_truth(sim_gen_name);
	outfname_truth.append(".truth.coords");	
	std::ofstream outfs_truth(outfname_truth);
	
	for (type_list_Changed_region::const_iterator it_change = affected_list.begin();
		 it_change != affected_list.end();
		++it_change)
	{
		outfs_truth 
			<< it_change->chr << "\t" 
			<< it_change->region.lower() << "\t"
			<< it_change->region.upper() << "\t"
			<< (it_change->rearrangement == hap_outcome__Dup) << "\n";
	}//change
	
	outfs_truth.close();
	
	
	std::cerr << "save events.\n";
	//write events	
	outfname_truth = sim_gen_name;
	outfname_truth.append(".truth.evs");
	outfs_truth.open(outfname_truth);		
	
	for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_samp = sampled_sim_Evs.begin();
		it_samp != sampled_sim_Evs.end();
		++it_samp)
	{
		outfs_truth
			<< it_samp->second.cc_of_event << "\t"
			<< it_samp->second.event_UID << "\t"
			<< it_samp->second.the_diploid_profile_brkpts.first.lower() << "\t"
			<< (it_samp->second.the_diploid_Event_outcome.first == hap_outcome__Dup) << "\n";		
	}//sampled_sim_Evs
	
	outfs_truth.close();
	
	std::cerr << "DONE with \"simulate_NAHR_events_on_Reference_genome\"\n\n";
	
	
	
			
} // simulate_NAHR_events_on_Reference_genome







bool compare_Changed_regions
				(const Changed_region &change1, Changed_region &change2)
{
	if (change1.chr < change2.chr)
	{  return true;  }
	else if (change1.chr == change2.chr)
	{
		if (change1.region.lower() < change2.region.lower())
			return true;
	}
	
	return false;
	
}//	compare_Changed_regions	
