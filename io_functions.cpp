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
#include <Interlocking_entry.h>
#include <templates.h>
#include <Coords_and_homologous_region.h>
#include <translations_through_alignment.h>
#include <Sparse_map.h>


#include "globals.h"
#include "io_functions.h"
#include "Conn_Comp.h"
#include "Event.h"
#include "Visitation_object.h"

#include "Call_set.h"

#include "other_functions.h"

//#include <boost/filesystem.hpp>







void read_Connected_Component_ALL_Event_data_from_file
                        (Conn_Comp &some_Conn_Comp)
{
    for (type_map_uint_to_Event::iterator it_ev = some_Conn_Comp.events.begin();
	    it_ev != some_Conn_Comp.events.end();
	    ++it_ev)
    {
	read_Event_ALL_data_from_file(some_Conn_Comp, it_ev->second);    
    }//ev
}//read_Connected_Component_ALL_data_from_file



void read_Connected_Component_basic_Event_data_from_file
                        (Conn_Comp &some_Conn_Comp)
{
    for (type_map_uint_to_Event::iterator it_ev = some_Conn_Comp.events.begin();
	    it_ev != some_Conn_Comp.events.end();
	    ++it_ev)
    {
	read_Event_basic_data_from_file(some_Conn_Comp, it_ev->second);    
    }//ev
}//read_Connected_Component_basic_Event_data_from_file


void read_Connected_Component_interlocking_Event_data_from_file
                        (Conn_Comp &some_Conn_Comp)
{
    for (type_map_uint_to_Event::iterator it_ev = some_Conn_Comp.events.begin();
	    it_ev != some_Conn_Comp.events.end();
	    ++it_ev)
    {
	read_Event_interlocking_data_from_file(some_Conn_Comp, it_ev->second);    
    }//ev
}//read_Connected_Component_interlocking_Event_data_from_file


void read_Connected_Component_regions_Event_data_from_file
                        (Conn_Comp &some_Conn_Comp)
{
    for (type_map_uint_to_Event::iterator it_ev = some_Conn_Comp.events.begin();
	    it_ev != some_Conn_Comp.events.end();
	    ++it_ev)
    {
	read_Event_regions_data_from_file(some_Conn_Comp, it_ev->second);    
    }//ev
}//read_Connected_Component_regions_Event_data_from_file






void read_Event_ALL_data_from_file
                        (Conn_Comp &some_Conn_Comp,
                          Event &some_Event)
{
    some_Event.clear_all_data();    
    
    read_Event_basic_data_from_file(some_Conn_Comp, some_Event);   
    read_Event_interlocking_data_from_file(some_Conn_Comp, some_Event);    
    read_Event_regions_data_from_file(some_Conn_Comp, some_Event);
    
} // read_Event_ALL_data_from_file   






void read_Event_basic_data_from_file
                        ( Conn_Comp &some_Conn_Comp,
                          Event &some_Event)
{

    //std::fprintf(stderr, "event_data_filename = %s\n", event_data_filename.c_str());
    Event *const current_ev = &some_Event;
    
    char event_data_filename[option_size];
    std::sprintf(event_data_filename, "%s/event_data/basic/basic_data_UID_%u_UIDDB_%u", 
                                        data_directory.c_str(), some_Event.UID, some_Event.UID_from_database);
  
    
    FILE *event_file = std::fopen(event_data_filename, "r");      
                                    if (event_file == NULL)
                                    {
                                        std::fprintf(stderr, "ERROR: unable to open file [%s] in \"read_Event_basic_data_from_file\".\n\n",
                                                    event_data_filename );
                                        exit(1);
                                    }
                       

    std::fscanf(event_file, "UID: %*u UID_from_database: %*u\nconnected_component: %*u\n");                                
    
    
    uint in_orient_agree, in_use_comp, in_LCR_0_0, in_LCR_0_1, in_LCR_1_0, in_LCR_1_1, in_swap_qs;
    std::fscanf(event_file, "chromos[0]: %u\nLCR_lengths[0]: %u\nLCRs[0]: %u %u\nchromos[1]: %u\nLCR_lengths[1]: %u\nLCRs[1]: %u %u\nrecomb_type: %u\nprofile_length: %u\norientations_agree: %u\nuse_complement_of_subject: %u\nswapped_query_and_subject_for_consistency: %u\n",
                  &(some_Event.chromos[0]),
                  &(some_Event.LCR_lengths[0]),                 
                  &in_LCR_0_0,
                  &in_LCR_0_1,
                  &(some_Event.chromos[1]),
                  &(some_Event.LCR_lengths[1]),
                  &in_LCR_1_0,
                  &in_LCR_1_1,                 
                  &(some_Event.recomb_type),
                  &(some_Event.profile_length),
                  &in_orient_agree,
                  &in_use_comp,
                  &in_swap_qs
                  );    
		  
    std::fclose(event_file);		 
    
    
    
                  
    some_Event.LCRs[0].set(in_LCR_0_0, in_LCR_0_1);  
    some_Event.LCRs[1].set(in_LCR_1_0, in_LCR_1_1);
                
    some_Event.orientations_agree = (bool)in_orient_agree;
    some_Event.use_complement_of_subject = (bool)in_use_comp;
    some_Event.swapped_query_and_subject_for_consistency = (bool)in_swap_qs;
                    
    
    //finally, initialize some values:
    if (some_Event.recomb_type == recomb_class__DupDel  or    some_Event.recomb_type == recomb_class__Inv)
    {
        some_Event.region_between_and_including_the_LCRs_themselves.set(  some_Event.LCRs[0].lower(),  some_Event.LCRs[1].upper()  );
        some_Event.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs.set(  some_Event.LCRs[0].upper() + 1,  some_Event.LCRs[1].lower() - 1  );
    }
    
    
    
    
    
    
    some_Event.potential_outcomes.clear();
    some_Event.potential_outcomes.insert(hap_outcome__None);
    if (some_Event.recomb_type !=  recomb_class__None)
    {
	some_Event.potential_outcomes.insert(hap_outcome__GeneConv_ABA);
	some_Event.potential_outcomes.insert(hap_outcome__GeneConv_BAB);	    
    }
    
    switch (some_Event.recomb_type)
    {
	case recomb_class__None:	    
	    break;
	case recomb_class__Inv:
	    some_Event.potential_outcomes.insert(hap_outcome__Inv);
	    break;
	case recomb_class__DupDel:
	    some_Event.potential_outcomes.insert(hap_outcome__Del);
	    some_Event.potential_outcomes.insert(hap_outcome__Dup);
	    break;	    
	case recomb_class__intrachromosomal_pq_Transloc:
	    some_Event.potential_outcomes.insert(hap_outcome__Transloc_A);
	    some_Event.potential_outcomes.insert(hap_outcome__Transloc_B);
	    break;
	case recomb_class__interchromosomal_Transloc:
	    some_Event.potential_outcomes.insert(hap_outcome__Transloc_A);
	    some_Event.potential_outcomes.insert(hap_outcome__Transloc_B);	    
	    break;	    
	default:
	    std::cerr <<"\n\n\n ERROR!  unidentified some_Event.recomb_type = " << some_Event.recomb_type << "\n\n\n";
	    exit(1);
	    break;
    };
        
    


    
}  //  read_Event_basic_data_from_file 




                   







void read_Event_interlocking_data_from_file
                        ( Conn_Comp &some_Conn_Comp,
                          Event &some_Event)
{

    Event *const current_ev = &some_Event;
    
    char event_data_filename[option_size];    
    std::sprintf(event_data_filename, "%s/event_data/interlocking/interlocking_events_for_UID_%u_UIDDB_%u", data_directory.c_str(), some_Event.UID, some_Event.UID_from_database);
    
    FILE *event_file = std::fopen(event_data_filename, "r");    
                                                    if (event_file == NULL)
                                                    {
                                                        std::fprintf(stderr, "ERROR: unable to open file [%s] in \"read_Event_mapping_data_from_file\".\n\n",
                                                                    event_data_filename );
                                                        exit(1);
                                                    }                                  
        
    uint num_entries_ctr;    
    
       
    
    
    std::fscanf(event_file, "local_interlocking_Events: %u\n", &num_entries_ctr );
    
    current_ev->local_interlocking_Events.clear();    
    type_map_uint_to_Interlock::iterator current_insert = current_ev->local_interlocking_Events.begin();
    
    for (uint j=0; j< num_entries_ctr; ++j)
    {
        uint in_entry_UID;
        std::fscanf(event_file, "\t%u\n", &in_entry_UID );
        
        current_insert = current_ev->local_interlocking_Events
                            .insert(current_insert, type_uint__Interlock(in_entry_UID,
                                                                        Interlocking_entry<Event>( &( some_Conn_Comp.events.at(in_entry_UID) ) )   ));                                                 
    }  
    
    
    
    
    
    
    std::fscanf(event_file, "global_interlocking_Events: %u\n", &num_entries_ctr );
    
    current_ev->global_interlocking_Events.clear();    
    current_insert = current_ev->global_interlocking_Events.begin();
    
    for (uint j=0; j< num_entries_ctr; ++j)
    {
        uint in_entry_UID;
        std::fscanf(event_file, "\t%u\n", &in_entry_UID );
        
        current_insert = current_ev->global_interlocking_Events
                            .insert(current_insert,
                                    type_uint__Interlock(in_entry_UID,
                                                        Interlocking_entry<Event> ( &(some_Conn_Comp.events.at(in_entry_UID) ) )));                                                  
    }      
       
    
    std::fclose(event_file);
    
    
    
}  // read_Event_interlocking_data_from_file

















void read_Event_regions_data_from_file
                        (Conn_Comp &some_Conn_Comp,
                         Event &some_Event)
{
        
    Event *const current_ev = &some_Event;
    
    char event_data_filename[option_size];      
    
//     std::fprintf(stderr,"regions\n");
    std::sprintf(event_data_filename, "%s/event_data/regions/homologous_regions_UID_%u_UIDDB_%u",
                                    data_directory.c_str(), current_ev->UID, current_ev->UID_from_database);
  
    FILE *event_file = std::fopen(event_data_filename, "r");      
                                                if (event_file == NULL)
                                                {
                                                std::fprintf(stderr, "ERROR: unable to open file [%s] in read_Event_regions_data_from_file.\n\n",
                                                            event_data_filename );
                                                exit(1);
                                                }
    
        


    current_ev->variational_positions_of_profile.clear();
        
    uint var_pos_size;
    type_set_uint::iterator insert_it = current_ev->variational_positions_of_profile.begin();    
    
    std::fscanf(event_file, "variational_positions_of_profile: size = %u\n", &var_pos_size );
    for (uint j=0; j<var_pos_size; ++j)                      
    {
        uint in_varpos_val;
        std::fscanf(event_file, "\t\t%u\n", &in_varpos_val );
        insert_it = current_ev->variational_positions_of_profile.insert(insert_it, in_varpos_val);
    }
    
    
    current_ev->uploaded_variational_positions_of_profile = current_ev->variational_positions_of_profile;

 
    current_ev->large_gaps_of_profile.clear();
 
    uint number_of_large_gaps;
    std::fscanf(event_file, "large_gaps_of_profile: size = %u\n", &number_of_large_gaps);
    for (uint j=0; j<number_of_large_gaps; ++j) 
    {
        uint in_lower_BI, in_upper_BI;
        std::fscanf(event_file, "\t\t[%u, %u]\n", &in_lower_BI, &in_upper_BI);
        current_ev->large_gaps_of_profile.insert( BOOST_Interval(in_lower_BI, in_upper_BI) );
    }
  
      
    current_ev->all_NON_large_gap_intervals_of_profile.clear();
//     current_ev->construct_all_NON_large_gap_intervals_of_profile();
    
    
    
  
  
    
    current_ev->compressed_map_LCR_to_profile__for_each_LCR.clear();
  
    current_ev->compressed_map_LCR_to_profile__for_each_LCR.push_back( type_vector_vector_int() );
    current_ev->compressed_map_LCR_to_profile__for_each_LCR.push_back( type_vector_vector_int() );

    for (uint qs=0; qs<2; ++qs)
    {     
        uint compressed_size;  
        current_ev->compressed_map_LCR_to_profile__for_each_LCR.at(qs).push_back( type_vector_int() );
        current_ev->compressed_map_LCR_to_profile__for_each_LCR.at(qs).push_back( type_vector_int() );
        
        std::fscanf(event_file, "compressed_map_LCR_to_profile__for_each_LCR[%*u]: size = %u\n", &compressed_size );
        
        current_ev->compressed_map_LCR_to_profile__for_each_LCR.at(qs).at(0).reserve(  compressed_size );
        current_ev->compressed_map_LCR_to_profile__for_each_LCR.at(qs).at(1).reserve(  compressed_size );
                            
        for ( uint j=0; j < compressed_size; ++j )
        {
            int in_key, in_val;
            std::fscanf(event_file, "\t\t%d --> %d\n", &in_key, &in_val);
            current_ev->compressed_map_LCR_to_profile__for_each_LCR.at(qs).at(0).push_back(in_key);
            current_ev->compressed_map_LCR_to_profile__for_each_LCR.at(qs).at(1).push_back(in_val);                
        }
    }  
  
  
  
    
    current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.clear();
    
    current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.push_back( type_vector_vector_int() );
    current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.push_back( type_vector_vector_int() );
    
    for (uint qs=0; qs<2; ++qs)
    {     
        uint compressed_size;
        current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs).push_back( type_vector_int() );
        current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs).push_back( type_vector_int() );
        
        std::fscanf(event_file, "compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR[%*u]: size = %u\n", &compressed_size );
        
        current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs).at(0).reserve(compressed_size);
        current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs).at(1).reserve(compressed_size);
                            
        for ( uint j=0; j < compressed_size; ++j)
        {
            int in_key, in_val;
            std::fscanf(event_file, "\t\t%d --> %d\n", &in_key, &in_val);
            current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs).at(0).push_back(in_key);
            current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs).at(1).push_back(in_val);                
        }    
    }

  




    
    uint num_entries_ctr;  
    
    
    
    
    current_ev->regions_homologous_to_directly_affected_region.clear();
    
    if (current_ev->recomb_type == recomb_class__DupDel  or    current_ev->recomb_type == recomb_class__Inv)
    {  
                        
        std::fscanf(event_file, "regions_homologous_to_directly_affected_region: %u\n", &num_entries_ctr ); 
        
        for (uint j=0; j< num_entries_ctr; ++j)
        {
            uint in_map_size, in_chrom;
            int in_map_key, in_map_val;
            type_vector_vector_int in_compressed_map;
            in_compressed_map.push_back(type_vector_int() );
            in_compressed_map.push_back(type_vector_int() );
            
            
            std::fscanf(event_file,"\tchromosome_of_homologous_region: %u\n", &in_chrom);                    
            
            uint in_orientation_of_homologous_region_and_profile_agree, in_homologous_region_must_be_complemented_to_get_to_profile;
            
            std::fscanf(event_file, "\torientation_of_homologous_region_and_profile_agree: %u\n",
                                                            &in_orientation_of_homologous_region_and_profile_agree);
            std::fscanf(event_file, "\thomologous_region_must_be_complemented_to_get_to_profile: %u\n",
                                                            &in_homologous_region_must_be_complemented_to_get_to_profile);                                        
            
            std::fscanf(event_file,"\tcompressed_map_homologous_region_to_profile: size = %u\n", &in_map_size );
            
            
            in_compressed_map.at(0).reserve(in_map_size);
            in_compressed_map.at(1).reserve(in_map_size);
                
            for (uint j=0; j<in_map_size; ++j)
            {
                std::fscanf(event_file,"\t\t%d --> %d\n", &in_map_key, &in_map_val);
                in_compressed_map[0].push_back(in_map_key);
                in_compressed_map[1].push_back(in_map_val);      
            }                                        
            

            current_ev->regions_homologous_to_directly_affected_region
                                        .push_back( Coords_and_homologous_region( in_chrom,
                                                                                  in_compressed_map,
                                                                                  (bool)in_orientation_of_homologous_region_and_profile_agree,
                                                                                  (bool)in_homologous_region_must_be_complemented_to_get_to_profile));
        
        } // end for-loop   num_entries_ctr    
                       
    }  //end if  recomb_type <= 2
    
        
                
    std::fclose(event_file);                
    
  
    
    
    
    
    
    current_ev->set_absolute_coordinate_mapping_into_large_gaps();
    
    
    
    
    //coordinates_over_which_to_consider_read_depth
    current_ev->remaining_coordinates_for_consideration.clear();
    type_set_uint::iterator insert_it_coordinates_to_consider = current_ev->remaining_coordinates_for_consideration.begin();
    
    if (current_ev->recomb_type == recomb_class__DupDel  or    current_ev->recomb_type == recomb_class__Inv)    
    {
        for (uint j=current_ev->region_between_and_including_the_LCRs_themselves.lower(); 
                 j<=current_ev->region_between_and_including_the_LCRs_themselves.upper();
             ++j)
            insert_it_coordinates_to_consider = current_ev->remaining_coordinates_for_consideration.insert(insert_it_coordinates_to_consider, j);
    }
    
//     else    
//         for (uint j=0; j<current_ev->profile_length; ++j)
//             insert_it_coordinates_to_consider = current_ev->remaining_coordinates_for_consideration.insert(insert_it_coordinates_to_consider, j);
        
        
        
        
//     // control size of  breakpoint compute:    
//     {//scope        
//         std::fflush(stderr);
//         
//         if (my_MPI_rank == 0)
//             std::fprintf(stderr, "\nEvent %u:\t\t\ttrue varpos size = %u,", current_ev->UID, var_pos_size);
//         
//         type_set_uint restricted_varpos;
//         type_set_uint::iterator insert_restricted_it = restricted_varpos.begin();        
//         insert_restricted_it = restricted_varpos.insert(insert_restricted_it, *current_ev->variational_positions_of_profile.begin() );               
//         
//         for (type_set_uint::iterator it_check_ahead = ++( ++current_ev->variational_positions_of_profile.begin() ),
//                                      it_next = ++current_ev->variational_positions_of_profile.begin();
//                 it_check_ahead != current_ev->variational_positions_of_profile.end();
//                 ++it_check_ahead,
//                 ++it_next)
//             if (*it_check_ahead - *insert_restricted_it > 10)           
//                 insert_restricted_it = restricted_varpos.insert(insert_restricted_it, *it_next);
//                                    
//         insert_restricted_it = restricted_varpos.insert(insert_restricted_it, *--current_ev->variational_positions_of_profile.end() );  //always insert last one.  
//         
//         current_ev->variational_positions_of_profile = restricted_varpos;
//         if (my_MPI_rank == 0)
//             std::fprintf(stderr, "\t\t\t\trestricted varpos size = %u\n\n", restricted_varpos.size() );
//         
//         std::fflush(stderr);
//     }//scope        
//     
//     



    current_ev->my_original_Breakpoint_complex.variational_positions__ALL = current_ev->variational_positions_of_profile;
    
//     
//     
//     for (uint qs=0; qs<2; ++qs)
//         current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs) =                    
//             convert_full_map_to_compressed_map(
//                 get_inverse_of_map(                              // LCR to remaining varpos
//                     get_map_intersection_of_map_keys_and_set<uint, uint>(
//                             get_inverse_of_map(                  // profile varpos to LCR
//                                 convert_compressed_map_to_full_map(
//                                     current_ev->
//                                         compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs) )), // (begin here)
//                                     current_ev->variational_positions_of_profile  ) )   );



// //check the sparsity of variational positions:
//     if (my_MPI_rank == 0 and current_ev->original_variational_positions.size() > 0)
//     {
//         uint number_of_sparse_vp = 0;
//         for (type_set_uint::const_iterator it_vp_below = current_ev->variational_positions_of_profile.begin(),
//                     it_vp = ++current_ev->variational_positions_of_profile.begin(),
//                     it_vp_above = ++(++current_ev->variational_positions_of_profile.begin());                
//                 it_vp_above != current_ev->variational_positions_of_profile.end();
//                 ++it_vp_below, ++it_vp, ++it_vp_above)
// 	{
//             if (*it_vp_above - *it_vp > maximum_average_fragment_length   or   *it_vp - *it_vp_below > maximum_average_fragment_length)
//                 ++number_of_sparse_vp;
// 	}
//             
//             
// //         std::fprintf(stderr, "\nEvent %u:  number of sparse variational positions =   %u\n\n",
// //                     current_ev->UID, number_of_sparse_vp);
//     }
    
    
    
    
//     if (my_MPI_rank == 0  and  !current_ev->large_gaps_of_profile.empty() )
//     {
//         
//         std::stringstream warning_strm;
//         bool there_is_a_warning = false;
//         
//         
//         print_line_of_markers("WARNING! ", &warning_strm);
//         print_line_of_markers("(", &warning_strm);                            
//         
//         if (get_intersection_of_set_and_interval(   current_ev->variational_positions_of_profile,
//                                                     BOOST_Interval(0, current_ev->large_gaps_of_profile.begin()->lower() - 1)).empty() )
//         {
//             warning_strm << "\n\nEVENT UID " << current_ev->UID <<"  has no varpos  below  first large gap!!!!!!!!!!!!!!\n\n\n";        
//             there_is_a_warning = true;                                                  
//         }
//         if (get_intersection_of_set_and_interval(   current_ev->variational_positions_of_profile,
//                                                     BOOST_Interval((--current_ev->large_gaps_of_profile.end())->upper() + 1, current_ev->profile_length)).empty() )
//         {
//             warning_strm << "\n\nEVENT UID " << current_ev->UID <<"  has no varpos  above  last large gap!!!!!!!!!!!!!!\n\n\n";      
//             there_is_a_warning = true;    
//         }    
//         
//         if (current_ev->large_gaps_of_profile.size() > 1)        
//             for (type_set_BI::const_iterator it_gap = current_ev->large_gaps_of_profile.begin(),
//                                             it_next_gap = ++current_ev->large_gaps_of_profile.begin();
//                     it_next_gap != current_ev->large_gaps_of_profile.end();
//                     ++it_gap, ++it_next_gap)
//                 if (get_intersection_of_set_and_interval(   current_ev->variational_positions_of_profile,
//                                                             BOOST_Interval(it_gap->upper()+1, it_next_gap->lower()-1)).empty()  )
//                 {
//                     warning_strm << "\n\nEVENT UID " << current_ev->UID <<"  has no varpos inbetween large gaps "
//                                 << *it_gap  << " and " << *it_next_gap << "  !!!!!!!!!!!!!!\n\n\n";     
//                     there_is_a_warning = true;    
//                 }  
//                 
//         if (there_is_a_warning)
//         {
//             print_line_of_markers(")", &warning_strm);
//             std::fprintf(stderr, "\n\n%s\n\n", warning_strm.str().c_str() );           
//         }                                   
//     }
    
    
    
    {//last nongap point
	const BOOST_Interval profile_endpoints__LCR_0(get_endpoints_of_values_of_compressed_map(current_ev->compressed_map_LCR_to_profile__for_each_LCR.at(0)));
	const BOOST_Interval profile_endpoints__LCR_1(get_endpoints_of_values_of_compressed_map(current_ev->compressed_map_LCR_to_profile__for_each_LCR.at(1)));
	
	current_ev->last_profile_nongap_index = std::min<uint>(profile_endpoints__LCR_0.upper(), profile_endpoints__LCR_1.upper());
	
	current_ev->variational_positions_of_profile.insert(current_ev->last_profile_nongap_index);
	current_ev->my_original_Breakpoint_complex.variational_positions__ALL.insert(current_ev->last_profile_nongap_index);
	
        for (uint qs=0; qs<2; ++qs)    
	{
            current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)                        
                =  get_compressed_map_inverse_of_compressed_map(  //LCR to profile varpos
                        convert_full_map_to_compressed_map(  
                            get_map_intersection_of_compressed_map_keys_and_set(
                                    get_compressed_map_inverse_of_compressed_map(  //profile varpos to LCR
                                            current_ev->compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)),
                                    current_ev->variational_positions_of_profile)));            
	}
    }//last nongap point

        
        

} //  read_Event_mapping_data_from_file





























void read_connected_component_data_from_file_and_create_Event_shells
                            (const char *const &in_connected_components_filename,
                             type_map_uint_to_CC &Conn_Comps)
{
    Conn_Comps.clear();

    FILE *in_connected_components_file = std::fopen(in_connected_components_filename, "r");
        
    if (in_connected_components_file == NULL)
    {
        std::fprintf(stderr, "ERROR: unable to open file [%s] in read_connected_component_data_from_file_and_create_Event_shells.\n\n",
                    in_connected_components_filename );
        exit(1);
    }
    
    
    type_map_uint_to_CC::iterator current_CC_it = Conn_Comps.begin();
    uint in_UID, in_UID_DB, in_CC_ID;
    
    while ( fscanf(in_connected_components_file, "%u  %u  %u\n", &in_UID, &in_UID_DB, &in_CC_ID) != EOF )
    {
        current_CC_it = Conn_Comps.insert(current_CC_it,  type_uint__CC(in_CC_ID, Conn_Comp(in_CC_ID) ));
        
        current_CC_it->second.events.insert( type_uint__Event(in_UID,
                                                              Event(in_UID, in_UID_DB, in_CC_ID, &(current_CC_it->second)    ) ));
    }  
    
    fclose(in_connected_components_file);
  
  
}  // end of   read_connected_component_data_from_file_and_create_Event_shells



















void read_chromosome_data_from_file()
{

    type_map_uint_to_uint::iterator insert_length_it = Event::chromosome_lengths.begin();
    type_map_uint_to__uint__uint::iterator insert_cent_it = Event::centromere_coordinates.begin();
    
    char chromo_filename[option_size];
    std::sprintf(chromo_filename, "%s/chromosome_data", data_directory.c_str());
    
        
    FILE *chromos_file = std::fopen(chromo_filename, "r");  
                        if (chromos_file == NULL)
                        {
                            std::fprintf(stderr, "\n\nERROR:  unable to open file [%s] in read_chromosome_data_from_file\n\n\n",
                                        chromo_filename);
                            exit(1);
                        }  
    
    uint chromo_val = 1;
    uint length, cent_begin, cent_end;    
    
    while (  std::fscanf(chromos_file, "chromosome %*u: length = %u, centromere = [%u, %u]\n",
                            &length,
                            &cent_begin,
                            &cent_end)           !=   EOF) 
    {
                            
        insert_length_it = Event::chromosome_lengths.insert( insert_length_it,
                                                            type_uint__uint(chromo_val, length) );                              
        insert_cent_it = Event::centromere_coordinates.insert( insert_cent_it,
                                                            type_uint___uint__uint(chromo_val, type_uint__uint(cent_begin, cent_end) )  );
        ++chromo_val;                                                              
    }
                                    
    
    std::fclose(chromos_file);

} //  end  of   read_chromosome_data_from_file
























bool read_permissibile_state_vectors_from_file_for_summing_out_a_given_UID_and_tack_on_empty_state_vector
                                                      (Visitation_object &my_Visit_object)
{
    
  
    //First load the UIDs (correctly ordered by ascending UID in the file) to which the state vector positioins correspond:
    char permissible_state_vector_filename[option_size];  
    std::sprintf(permissible_state_vector_filename, "%s/%s/neighbors/UIDs_for_permissible_state_vectors_UID_%u",
                data_directory.c_str(),
                 visitation_schedule_subdirname.c_str(),
                my_Visit_object.event_for_elimination->UID);    
    
    FILE *UIDs_list = std::fopen(permissible_state_vector_filename, "r");    
                if (UIDs_list == NULL)
                {
                    std::fprintf(stderr, "\n\nERROR:  unable to open file [%s] in read_permissibile_state_vectors_from_file_for_summing_out_a_given_UID\n\n\n",
                                permissible_state_vector_filename);
                    return false;
                }
    
  
  
    
    type_set_uint::iterator in_neigh_it = my_Visit_object.neighbors.begin();
    
    uint in_UID;
    while (std::fscanf(UIDs_list, "%u", &in_UID) != EOF)
    {  in_neigh_it = my_Visit_object.neighbors.insert(in_neigh_it, in_UID);  }
            
    std::fclose(UIDs_list);
    

    
    
  
  
  
  
  
  
  
  //Now upload each state vector itself:
  
    //first add the null state vector (no events occur), since it is always possible, but never actually saved to file (since we deal in sparse maps).
    my_Visit_object.permissible_haploid_state_vectors_sparse.push_back( Sparse_map() );
  
  
  
    std::sprintf(permissible_state_vector_filename, "%s/%s/state_vectors/permissible_sparse_state_vectors_UID_%u",
                data_directory.c_str(),
                 visitation_schedule_subdirname.c_str(),
                my_Visit_object.event_for_elimination->UID);    
    
    FILE *state_vector_file = std::fopen(permissible_state_vector_filename, "r");    
                if (state_vector_file == NULL)
                {
                    std::fprintf(stderr, "\n\nERROR:  unable to open file [%s] in read_permissibile_state_vectors_from_file_for_summing_out_a_given_UID\n\n\n",
                                permissible_state_vector_filename);
                    return false;
                }  
  
  
  
  
  //Sparse_map::read_next_N_Sparse_maps_from_filestream(state_vector_file, -1);  //-1 --> N = infinity
    
    
    uint in_sparse_map_size;
    
    while( std::fscanf(state_vector_file, "%u: ", &in_sparse_map_size) != EOF )
    {
        
        my_Visit_object.permissible_haploid_state_vectors_sparse.push_back( Sparse_map() );
        
        Sparse_map *const newest_state_vector = &my_Visit_object.permissible_haploid_state_vectors_sparse.back();        
        type_map_uint_to_uint::iterator current_sparse_it = newest_state_vector->sparse_map_UID_to_value.begin();
        
        uint in_map_UID, in_map_state_value;
        
        for (uint j=0; j<in_sparse_map_size; ++j)
        {
            std::fscanf(state_vector_file, "[%u, %u] ", &in_map_UID, &in_map_state_value);
            current_sparse_it = newest_state_vector->sparse_map_UID_to_value.insert( current_sparse_it,  type_uint__uint(in_map_UID, in_map_state_value)   );
        }
                    
    }  // end while != EOF
            
    
    std::fclose(state_vector_file);
        
    
    return true;
    

}  //end of    read_permissibile_state_vectors_from_file_for_summing_out_a_given_UID_and_tack_on_empty_state_vector






























int  read_total_diploid_compute_size_from_file_for_given_Conn_Comp
                                     (Conn_Comp &my_CC)
{    

    char schedule_filename [option_size];
    std::sprintf(schedule_filename, "%s/%s/schedules/Conn_Comp_%u",
                                    data_directory.c_str(),
                                    visitation_schedule_subdirname.c_str(),
                                    my_CC.Conn_Comp_ID);                                                                       
                                                                
    FILE *schedule_file = std::fopen(schedule_filename, "r");    
                        if (schedule_file == NULL)
                        {
                            std::fprintf(stderr, "\n\nERROR:  unable to open file [%s] in read_visitation_schedule_from_file_for_given_Conn_Comp\n\n\n",
                                        schedule_filename);
                            return -1;
                        }               
  
    
    int in_total_diploid_compute_size;            
    std::fscanf(schedule_file, "visitation_schedule: total_diploid_compute = %d\n", &in_total_diploid_compute_size); 
        
    
    std::fclose( schedule_file );

    return in_total_diploid_compute_size;
    

}  // read_total_diploid_compute_size_from_file_for_given_Conn_Comp








bool read_visitation_schedule_from_file_for_given_Conn_Comp
                                        (Conn_Comp &my_CC)
{

    char schedule_filename [option_size];
    std::sprintf(schedule_filename, "%s/%s/schedules/Conn_Comp_%u",
                                    data_directory.c_str(),
                                    visitation_schedule_subdirname.c_str(),
                                    my_CC.Conn_Comp_ID);                                                                       
                                                                
    FILE *schedule_file = std::fopen(schedule_filename, "r");    
                        if (schedule_file == NULL)
                        {
                            std::fprintf(stderr, "\n\nERROR:  unable to open file [%s] in read_visitation_schedule_from_file_for_given_Conn_Comp\n\n\n",
                                        schedule_filename);
                            return false;
                        }   
        
    
  
    
    uint in_UID, in_diploid_compute_size;            
    std::fscanf(schedule_file, "visitation_schedule: total_diploid_compute = %u\n", &my_CC.total_diploid_compute_size); 
    

    
    
    while (  EOF  !=  std::fscanf(schedule_file, "UID: %u\tdiploid_compute_size = %u\n",
                                                &in_UID,
                                                &in_diploid_compute_size) )
    {                
        
        type_map_uint_to_Event::iterator located_it = my_CC.events.find(in_UID);
                        if ( located_it == my_CC.events.end() )
                        {
                            std::fprintf(stderr, "ERROR: visitation schedule called for UID %u, but such an event could not be found in Connected Component %u\n\n\n",
                                                in_UID, my_CC.Conn_Comp_ID);
                            return false;
                        }
        
        char temporary_PSF_save_filename[option_size];
        std::sprintf(temporary_PSF_save_filename, "%s/Partial_sum_functions/PSF_eliminating_Event_%u", output_dir.c_str(), in_UID);
                                                                                
        my_CC.variable_elimination_schedule.push_back( Visitation_object( &(located_it->second),
                                                                          in_diploid_compute_size,
                                                                          temporary_PSF_save_filename)   );
    }  // end of   while
    
        
    std::fclose(schedule_file); 
    
    
   
    
    return true;
       
    
}  // end of   read_visitation_schedule_from_file_for_given_Conn_Comp























//copied from my  make_LCR_seqs   program
//This version of obtaining DNA assumes that there are no new characters '\n' in the file except for the header e.g. ">char\n"
std::string get_specific_DNA_sequence_from_Reference_genome
                                              (const uint &chromo,
                                               const BOOST_Interval &coords)
{
//Description: Given coordinates ("coords") on a certain chromosome ("chromo") pertaining to the Reference genome, open the relevant .fa file and extract the DNA sequence specified by "coords" to "destination_seq".

    //std::fprintf(stderr, "chromos = %u, coords = [%u, %u]\n" , chromo, coords.lower(), coords.upper() );
        
	if (BOOST_empty(coords))
	{
		std::stringstream error_strm;
		error_strm << "ERROR!  coords = " << coords << " is EMPTY in \"get_specific_DNA_sequence_from_Reference_genome\"!!!\n";
		error_message(error_strm, true);		
	}
    
    
    char Ref_genome_chromo_filename[option_size];  
        
    if (genome_tag == 0)
    {
        if (chromo < 23)
            std::sprintf(Ref_genome_chromo_filename, "%s%u.fa", Ref_genome_fileprefix.c_str(), chromo);
        else if (chromo == 23) //X chromosome
            std::sprintf(Ref_genome_chromo_filename, "%sX.fa", Ref_genome_fileprefix.c_str());    
        else if (chromo == 24) // Y chromosome
            std::sprintf(Ref_genome_chromo_filename, "%sY.fa", Ref_genome_fileprefix.c_str());
		else
		{
			std::stringstream error_strm;
			error_strm << "ERROR!  invalid chromosome in \"get_specific_DNA_sequence_from_Reference_genome\",  chromo = " << chromo << "\n";
			error_message(error_strm, true);
		}
    }
    else  //simulate    
	{  std::sprintf(Ref_genome_chromo_filename, "%s/reference_genome_%u", data_directory.c_str(), genome_tag);  }
    
    
    FILE *Ref_genome_file = std::fopen(Ref_genome_chromo_filename, "r");                        
                        if (Ref_genome_file == NULL)
                        {    
                            std::stringstream error_strm;
                            error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";                            
                            error_strm << "\n\nERROR!   File [" << Ref_genome_chromo_filename  
                                        << "] could not be opened in function \"get_specific_DNA_sequence_from_Reference_genome\".\n\n"  
                                        << "\nchromo = " << chromo
                                        << "\ncoords = " << coords 
                                        << "\nRef_genome_fileprefix.c_str() = " << Ref_genome_fileprefix.c_str()
                                        << "\n\n";
                                        
                            error_message(error_strm, true);
                        }         
  
         
    //assume first line is just a header  line ("e.g. >chr11"), so skip it.
    while(std::fgetc(Ref_genome_file) != '\n')
    {   }


    const uint length_of_desired_sequence = BOOST_width_inclusive(coords);    
    
    
    
    
    
    int pos_ctr = 1;  //recall that Reference genome is 1-based
    while (pos_ctr < coords.lower() )
    {
        const int nucleotide = std::fgetc(Ref_genome_file);   
        if (nucleotide != '\n')       
		{  ++pos_ctr;  }
    }  
        
    
    //now read in desired nucleotides
    pos_ctr = 0;
    std::string desired_sequence(length_of_desired_sequence, '\0');
    while (pos_ctr < length_of_desired_sequence)
    {
        const int nucleotide = std::fgetc(Ref_genome_file);   
        if (nucleotide != '\n')  
        {
            desired_sequence.at(pos_ctr) = (char)nucleotide;
            ++pos_ctr;
        }
    }  
                  
    std::fclose(Ref_genome_file);
 
    
    //make_sequence_all_uppercase:
    for (uint j=0; j<length_of_desired_sequence; ++j)
    {       
        const int char_as_int = (int)desired_sequence.at(j);
        if (char_as_int > 96)    
            desired_sequence.at(j) = (char)(char_as_int - 32);           
    }
     
     
                 
    return desired_sequence;
    
} //end of get_specific_DNA_sequence_from_Ref
















void upload_universal_variational_positions()
{
    
    char varpos_filename[option_size];
        
    for (type_map_uint_to_uint::const_iterator it_chr = Event::chromosome_lengths.begin();
            it_chr != Event::chromosome_lengths.end();
            ++it_chr)
    {                    
        std::sprintf(varpos_filename, "%s/universal_variational_positions_for_chromosome_%u", variational_positions_directory.c_str(), it_chr->first);             
        
        FILE *varpos_file = std::fopen(varpos_filename, "r");        
                    if (varpos_file == NULL)
                    {
                        std::fprintf(stderr, "ERROR: unable to open file [%s] in function \"upload_universal_variational_positions\"\n\n\n",
                                            varpos_filename);
                        exit(1);
                    }
        
        type_map_uint_to_set_uint::iterator insert_it_chr =  universal_variational_positions.insert( universal_variational_positions.begin(),
                                                                                                     type_uint__set_uint(it_chr->first, type_set_uint()) );
        type_set_uint::iterator insert_universal_varpos_it = insert_it_chr->second.begin();        
        
        uint varpos_value;                
        while ( std::fscanf(varpos_file, "%u\n", &varpos_value) != EOF)
            insert_universal_varpos_it = insert_it_chr->second.insert( insert_universal_varpos_it, varpos_value );                      
        
        std::fclose(varpos_file);                              
    }       

}  // upload_universal_variational_positions





















void upload_universal_indel_positions()
{
    
    char indel_filename[option_size];
        
    for (type_map_uint_to_uint::const_iterator it_chr = Event::chromosome_lengths.begin();
            it_chr != Event::chromosome_lengths.end();
            ++it_chr)
    {           
        {//del
            std::sprintf(indel_filename, "%s/universal_del_positions_for_chromosome_%u", variational_positions_directory.c_str(), it_chr->first);             
            
            FILE *del_file = std::fopen(indel_filename, "r");        
                                if (del_file == NULL)
                                {
                                    std::fprintf(stderr, "ERROR: unable to open file [%s] in function \"upload_universal_indel_positions\"\n\n\n",
                                                        indel_filename);
                                    exit(1);
                                }
                                
            type_map_uint_to_set_uint::iterator insert_it_chr =  universal_del_positions.insert(   universal_del_positions.begin(),
                                                                                                    type_uint__set_uint(it_chr->first, type_set_uint() )   );
            type_set_uint::iterator insert_universal_del_it = insert_it_chr->second.begin();        
            
            uint del_lower;          
            while ( std::fscanf(del_file, "%u\n", &del_lower) != EOF  )
                insert_universal_del_it = insert_it_chr->second.insert( insert_universal_del_it, del_lower);
            
            std::fclose(del_file);  
        }//del
        
        
        
        {//ins
            std::sprintf(indel_filename, "%s/universal_ins_positions_for_chromosome_%u", variational_positions_directory.c_str(), it_chr->first);             
            
            FILE *ins_file = std::fopen(indel_filename, "r");        
                                if (ins_file == NULL)
                                {
                                    std::fprintf(stderr, "ERROR: unable to open file [%s] in function \"upload_universal_indel_positions\"\n\n\n",
                                                        indel_filename);
                                    exit(1);
                                }
                                
            
            universal_ins_positions.insert(   type_uint__set_uint(it_chr->first, type_set_uint() )   );
            type_set_uint::iterator insert_universal_ins_it = universal_ins_positions.at(it_chr->first).begin();        
            
            uint ins_lower;          
            while ( std::fscanf(ins_file, "%u\n", &ins_lower) != EOF  )
                insert_universal_ins_it = universal_ins_positions.at(it_chr->first).insert( insert_universal_ins_it, ins_lower);
            
            std::fclose(ins_file);        
        }//ins
    }       

}  // upload_universal_variational_positions













































void read_Validated_Events_from_file
                            (const char *const &validated_Events_filename,
                             type_list_Validated_rearrangement &Validated_Events)
{        
    
    FILE *validated_file = std::fopen(validated_Events_filename, "r");
                        if (validated_file == NULL)
                        {    
                            std::fprintf(stderr, "ERROR: File [%s] could not be opened in function \"read_Validad_Events_from_file\".\n\n", 
                                                validated_Events_filename);
                            exit(1);
                        }   
    
    
    uint chromo_val, brkpt_begin, brkpt_end, in_recomb_type, in_outcome, number_of_related_Events;
    char sample_subject_name[option_size], supporting_clones[option_size], accession[option_size], published_study[option_size];    
    
    while (EOF !=   std::fscanf(validated_file, "%u:[%u..%u]   %u  --->  %u     related_UIDs:  %u   %s   %s   %s  %s\n",
                                                &chromo_val,
                                                &brkpt_begin, &brkpt_end,
                                                &in_recomb_type, &in_outcome,
                                                &number_of_related_Events,
                                                sample_subject_name, supporting_clones, 
                                                accession,   published_study)    )
        if (number_of_related_Events > 0)
        {
            type_set_uint in_related_UIDs;
            uint in_UID;
            for (uint j=0; j<number_of_related_Events; ++j)
            {
                std::fscanf(validated_file, "UID: %u     %*[^\n]\n", &in_UID);
                in_related_UIDs.insert(in_UID);
            }                               
            
            Validated_Events.push_back(   Validated_rearrangement(                        
                                                     chromo_val,
                                                        in_recomb_type,
                                                        in_outcome,
                                                        BOOST_Interval(brkpt_begin, brkpt_end),
                                                        std::string(accession),
                                                        std::string(supporting_clones), 
                                                        std::string(published_study),
                                                        std::string(sample_subject_name),
                                                        in_related_UIDs)     );
        }
      
    
    std::fclose(validated_file);   
    

}  // read_Validad_Events_from_file








































void save_sampled_Connected_Component_to_file
                        (FILE *&in_file,
                         const type_map_uint_to_Sampled_diploid_Event_data &sample)
{
        
    for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_ev = sample.begin();
            it_ev != sample.end();
            ++it_ev)
    {
        std::fprintf(in_file, "\tEvent  %u\n\t\tstates:\t%u\t%u\n\t\t\tbreakpoints:\t[%u,%u]\t[%u,%u]\n",
                                it_ev->second.event_UID,             
                                it_ev->second.the_diploid_Event_outcome.first, it_ev->second.the_diploid_Event_outcome.second,
                                it_ev->second.the_diploid_profile_brkpts.first.lower(), it_ev->second.the_diploid_profile_brkpts.first.upper(),
				it_ev->second.the_diploid_profile_brkpts.second.lower(), it_ev->second.the_diploid_profile_brkpts.second.upper()  );
    }
                                
    std::fprintf(in_file, "\n\n");        
} // save_sampled_Event_to_file








void load_sampled_Event_from_file
                        (FILE *&out_file,
                         type_map_uint_to__uint__uint__2 &sample)
{
        
    uint in_UID;
    type_uint__uint in_states, in_brkpts;
        
    std::fscanf(out_file, "\t\tEvent  %u\n\t\t\t\tstates:          %u    &&    %u\n\t\t\t\t\t\t\t\t\t\tbreakpoints:     %u    &&    %u\n",
                            &in_UID,             
                            &in_states.first,  &in_states.second,
                            &in_brkpts.first,  &in_brkpts.second  );
                            
    sample.first[in_UID] = in_states;
    sample.second[in_UID] = in_brkpts;
    
} // load_sampled_Event_from_file































bool read_calls_from_file
                    (const std::string &calls_filename,
                     type_map_string_to_list_Call_set  &read_calls,
		     const bool &ignore_gene_conversion_calls,
		     const bool &ignore_second_call_in_diploid_calls)
{
    read_calls.clear();
    
    FILE *in_calls_file = std::fopen( calls_filename.c_str(), "r" );
    
    if ( in_calls_file == NULL)
    {
        std::stringstream error_strm;
        print_line_of_markers("ERROR! ", &error_strm);
        print_line_of_markers("(", &error_strm);
        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
        
        error_strm << "\n\nERROR!   in_calls_file == NULL   in  \"read_calls_from_file\",\n\t\t calls_filename =    "  
        <<  calls_filename  << "\n\n";
        
        print_line_of_markers(")", &error_strm);
        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
        return false;
    }   
            
    
    
    char genome_name[option_size];
    char in_pop[option_size];
    char var_outcome[option_size];
    char P_outcome_char[option_size];
    char brkpt_and_chr_begin[option_size];    
    char brkpt_and_chr_end[option_size];
    char P_brkpt_char[option_size];
    uint amount_of_DNA_affected;
    uint cc_of_event;
    uint event_of_var;
    uint in_brkpt_index;
    char in_gene_conv[option_size];
    
    char in_chr_lower[option_size];
    char in_chr_upper[option_size];
    
//     std::fscanf(in_calls_file, "\n");
    
//format:   genome_name   population   recomb_outcome     P(outcome)    chr_lower     lower_brkpt     chr_upper    upper_brkpt     P(brkpt)    amount_affected      CC_ID   EV_ID            
            
    int number_of_incoming_fields = 0;
    while ( true               )
    {        
	long int last_tell;
	last_tell = ftell(in_calls_file);
	number_of_incoming_fields = std::fscanf(  in_calls_file,  "%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%u\t%u\t%u\t%u\t%[^\n]\n",
                                                    genome_name,
						    in_pop,
                                                    var_outcome,
                                                    P_outcome_char,
                                                    brkpt_and_chr_begin,
                                                    brkpt_and_chr_end,
                                                    P_brkpt_char,
                                                    &amount_of_DNA_affected,
                                                    &cc_of_event,
                                                    &event_of_var,
						    &in_brkpt_index,
                                                    in_gene_conv);
	if (EOF == number_of_incoming_fields)
	    break;
	else if (number_of_incoming_fields < 12 )
	{
	    fseek(in_calls_file, last_tell, SEEK_SET);
	    number_of_incoming_fields = std::fscanf(  in_calls_file,  "%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%u\t%u\t%u\t%[^\n]\n",
                                                    genome_name,
                                                    var_outcome,
                                                    P_outcome_char,
                                                    brkpt_and_chr_begin,
                                                    brkpt_and_chr_end,
                                                    P_brkpt_char,
                                                    &amount_of_DNA_affected,
                                                    &cc_of_event,
                                                    &event_of_var,
                                                    in_gene_conv);	    
	    std::sprintf(in_pop,"?");
	    in_brkpt_index = 999999;
	    if (number_of_incoming_fields != 10)
		assert(number_of_incoming_fields == 10);
	}
	
	
	
	
        longdouble P_outcome = -1.00L;
        longdouble P_brkpt = -1.00L;
        
        if (strcasecmp(P_outcome_char, "\"") == 0)        
            P_outcome = read_calls[ genome_name ].back().P_outcome;
        else
            P_outcome = (longdouble)atof(P_outcome_char);
            
        
        if ( strcasecmp(var_outcome, "none") != 0  )
        {
            if (strcasecmp(P_brkpt_char, "\"") == 0)          
                P_brkpt = read_calls[ genome_name ].back().P_brkpt;        
            else        
                P_brkpt = (longdouble)atof(P_brkpt_char);  
        }
       
        
        read_calls[ genome_name ].push_back(  Call_set(  genome_name,
							in_pop,
                                                         var_outcome,
                                                         P_outcome,
                                                         brkpt_and_chr_begin,
                                                         brkpt_and_chr_end,
                                                         P_brkpt,
                                                         cc_of_event,
                                                         event_of_var,
							 in_brkpt_index)     );   
                                                     
                                                                                                                                                                                  
        read_calls[ genome_name ].back().potentially_gene_conversion =     ( strcasecmp(in_gene_conv, "yes") ==  0)   ;
                                                         
    }
    
    
    std::fclose(  in_calls_file  );     
    
    
    
    
    
    
    if (ignore_gene_conversion_calls)
    {
	for (type_map_string_to_list_Call_set::iterator it_genome = read_calls.begin();
		it_genome != read_calls.end();
		++it_genome)
	{	
	    //ignore gene conversion, and make all calls haploid
	    type_list_Call_set::iterator it_call = ++(it_genome->second.begin());
	    while (it_call != it_genome->second.end() )
	    {	    
		if (it_call->potentially_gene_conversion)
		    it_call = it_genome->second.erase(it_call);
		else
		    ++it_call;
	    }
	}        
    } // ignore_gene_conversion_calls
	
	
	
	
    if (ignore_second_call_in_diploid_calls)
    {
	for (type_map_string_to_list_Call_set::iterator it_genome = read_calls.begin();
		it_genome != read_calls.end();
		++it_genome)
	{		
	    // make all calls haploid
	    type_list_Call_set::iterator it_call = ++(it_genome->second.begin());
	    while (it_call != it_genome->second.end() )
	    {	    		
		type_list_Call_set::iterator it_prev = it_call;
		--it_prev;
		
		if (  it_prev->event_of_var == it_call->event_of_var  )
		    it_call = it_genome->second.erase(it_call);
		else
		    ++it_call;		
	    }    			    
	} // genome	
    }        
    
    
    
    
 
    return true;

} // 




























bool read_selected_connected_components_file
                    (const std::string &selection_dirname,
                     type_set_uint &magic_CCs)
{            
    char selected_CCs_this_genome[option_size];
    std::sprintf(selected_CCs_this_genome, "%s/selected_connected_components_for_genome_%s",
                                            selection_dirname.c_str(),   genome_name.c_str()    );
            
    FILE *selection_file = std::fopen(selected_CCs_this_genome, "r");
                        if (selection_file == NULL)
                        {
                            std::fprintf(stderr, "ERROR: unable to open selection_filename [%s] for  \"read_selected_connected_components_file\".\n\n",
                                                selected_CCs_this_genome );
                            return false;
                        }                   
    
    uint in_CC_val;
    while (  EOF  !=  std::fscanf(selection_file, "%u\n", &in_CC_val)  )
        magic_CCs.insert( in_CC_val );        
    
    std::fclose(selection_file);        

    
    return true;

} // read_selected_connected_components_file






























void append_to_skipped_connected_component_file
                (const Conn_Comp  &skipped_CC,
                 const std::string &reason_for_skipping)
{
    char skipped_filename[option_size];
    std::sprintf(skipped_filename, "%s/skipped_connected_components",  output_dir.c_str()  );
    
    std::ofstream out_skip_strm(  skipped_filename,  std::ios::app);
    
    
    
    out_skip_strm << skipped_CC.Conn_Comp_ID  << "\t" << reason_for_skipping << "\n";        
    
    
    
    out_skip_strm.close();    

} // append_to_skipped_connected_component_file

































void read_Validated_rearrangements_from_file__in_tabular_format
                (const std::string in_validations_filename,
                 type_list_Validated_rearrangement &validated_Events)
{
    
    FILE *in_vals_file = std::fopen( in_validations_filename.c_str() , "r" );
                            if ( in_vals_file == NULL  )
                            {
                                std::stringstream error_strm;
                                print_line_of_markers("ERROR! ", &error_strm);
                                print_line_of_markers("(", &error_strm);
                                
                                error_strm << "\n\nERROR!   unable to open file   in  \"read_Validated_rearrangements_from_file__in_tabular_format\" \n\t\t table_filename =    ["  
                                <<  in_validations_filename  << "]\n\n";
                                
                                print_line_of_markers(")", &error_strm);
                                std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
                                exit(1);
                            }    


    //skip header
    std::fscanf(in_vals_file, "%*[^\n]\n");
//     in_vals_strm << "sample\tchromosome\tbreakpoint begin\tbreakpoint end\tstudy\ttype\tvariant_ID\n";


    char in_sample_name[option_size];
    uint in_chr;
    uint in_brkpt_lower;
    uint in_brkpt_upper;
    char in_published_study[option_size];
    char in_outcome[option_size];
    char in_variant_ID[option_size];
    
    while ( EOF  !=  std::fscanf(in_vals_file,
                                "%[^\t]\tchr%u\t%u\t%u\t%[^\t]\t%[^\t]\t%[^\n]\n",
                                in_sample_name,
                                &in_chr,
                                &in_brkpt_lower ,
                                &in_brkpt_upper,
                                in_published_study,
                                in_outcome,
                                 in_variant_ID)  )
    {
   
        
        uint rearrangement_val;
        uint in_recomb_type;
        if (strcasecmp(in_outcome, "inversion") == 0)
        {
            rearrangement_val = 1;
            in_recomb_type = 1;
        }
        else
        {
            in_recomb_type = 2;
            if (strcasecmp(in_outcome, "duplication") == 0)
                rearrangement_val = 2;
            else if (strcasecmp(in_outcome, "deletion") == 0)
                rearrangement_val = 1;
            else
            {
                std::stringstream error_strm;
                print_line_of_markers("ERROR! ", &error_strm);
                print_line_of_markers("(", &error_strm);
                error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                
                error_strm << "\n\nERROR!    unrecognized rearrangement in validated event  in   \"read_Validated_rearrangements_from_file__in_tabular_format\",  in_validations_filename = [" << in_validations_filename << "]\n\t\t  in_outcome = [" << in_outcome << "]\n\n";
                
                error_strm << "\n\nin_sample_name = [" << in_sample_name
                            << "]\nin_chr = [" << in_chr
                            << "]\nin_brkpt_lower = [" << in_brkpt_lower
                            << "]\nin_brkpt_upper = [" << in_brkpt_upper
                            << "]\nin_published_study = [" << in_published_study
                            << "]\nin_outcome = [" << in_outcome
                            << "]\n\n\n\n";
                
                print_line_of_markers(")", &error_strm);
                std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
                exit(1);
            }
        }
        
        validated_Events.push_back(  Validated_rearrangement(   in_chr,
                                                                in_recomb_type,
                                                                rearrangement_val,
                                                                BOOST_Interval(in_brkpt_lower, in_brkpt_upper),
                                                                in_variant_ID,
                                                                "no_supporting_clone_name",
                                                                in_published_study,
                                                                in_sample_name,
                                                                type_set_uint() )        );        
    }



    std::fclose(in_vals_file);

}  //read_Validated_rearrangements_from_file__in_tabular_format





































void  read_gender_and_population_from_file
                (const std::string &gender_and_population_filename)
{

    FILE *gender_file = std::fopen(gender_and_population_filename.c_str(), "r");  
                                    if (gender_file == NULL)
                                    {
                                        std::fprintf(stderr, "ERROR: unable to open file [%s] in \"read_gender_from_file\".\n\n",
                                                    gender_and_population_filename.c_str() );
                                        exit(1);
                                    }
        
    bool found_sample_in_file = false;
    
    
    char in_sample_name[option_size];
    char in_gender[option_size];
    char in_population[option_size];
    
    while (EOF !=  std::fscanf(gender_file, "%[^\t]\t%[^\t]\t%[^\n]", in_population, in_sample_name, in_gender)  )
        if ( strcasecmp(in_sample_name, genome_name.c_str()) == 0)
        {
            found_sample_in_file = true;
            
            population_abbreviation.assign( in_population );            
            
            if (strcasecmp(in_gender, "female") == 0   or   strcasecmp(in_gender, "f") == 0   )
                gender_of_individual_is_female = true;
            else if (strcasecmp(in_gender, "male") == 0   or   strcasecmp(in_gender, "m") == 0   )
                gender_of_individual_is_female = false;  
            else 
            {
                std::cerr << "\n\n\nERROR! unable to parse gender  for  genome name = [" << genome_name 
                            << "]   that was read in from file ["  << gender_and_population_filename.c_str() << "].     the gender that was readr = [" << in_gender << "]\n\n";
                exit(1);                        
            }
            break;                    
        }
        else
            std:fscanf(gender_file, "\n");                                        
        
        
    std::fclose(gender_file);


    if (!found_sample_in_file)
            {
                std::cerr << "\n\n\nERROR!unable to find  genome name = [" << genome_name 
                            << "]   in file ["  << gender_and_population_filename.c_str() << "].\n\n\n\n";
                exit(1);                        
            }  
            
            
    const std::string gender_name =    gender_of_individual_is_female   ?    "female" :  "male";        
       
    if (my_MPI_rank == 0)
        std:: cerr << "\n\n for genome name  [" << genome_name
                    << "]\n\t\tpopulation = [" << population_abbreviation 
                    << "]\n\t\tgender = " << gender_name  << "\n\n\n";
                                

}//read_gender_from_file

































void save_posterior_samples_to_file
                        (const type_vector_map_uint_to_Sampled_diploid_Event_data &samples_from_posterior,
                         const uint &CC_ID,
                         const uint &number_of_samples,
                         const std::string &some_output_dir)
{
    
    std::fprintf(stderr, "\n\n\nsaving samples...\n");
        
    char save_samples_filename[option_size];
    std::sprintf(save_samples_filename, "%s/Samples_from_posterior/samples_from_posterior_CC_%u__number_of_samples_%u",
                                        some_output_dir.c_str(), 
                                        CC_ID,
                                        number_of_samples);
        
    std::ofstream out_fs( save_samples_filename,   std::ios::binary | std::ios::out  );  //  std::ios_base::binary     
                    if ( !out_fs.good()  )
                    {
                        std::stringstream error_strm;
                        print_line_of_markers("ERROR! ", &error_strm);
                        print_line_of_markers("(", &error_strm);
                        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                        
                        error_strm << "\n\nERROR!   not good ofstream   in  \"save samples\"  in sample_from_posterior__and__save____marginal_distributions__and__outcome_Centroid__and_Shannon__Entropies,\n\t\t save_samples_filename =    "  
                        <<  save_samples_filename  << "\n\n";
                        
                        print_line_of_markers(")", &error_strm);
                        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
                    }
                                                        
    boost::archive::binary_oarchive out_archive( out_fs );
    
    out_archive << samples_from_posterior;            
    out_fs.close();                          


} // save_posterior_samples_to_file

















type_vector_map_uint_to_Sampled_diploid_Event_data load_posterior_samples_from_file
                        (const std::string &some_output_dir,
                         const uint &CC_ID,
                         const uint &number_of_samples)
{
    
    std::fprintf(stderr, "\n\n\nload_posterior_samples_from_file...\n");
        
    char load_samples_filename[option_size];
    std::sprintf(load_samples_filename, "%s/Samples_from_posterior/samples_from_posterior_CC_%u__number_of_samples_%u",
                                        some_output_dir.c_str(), 
                                        CC_ID,
                                        number_of_samples);
        
    std::ifstream in_fs( load_samples_filename,   std::ios::binary | std::ios::in  );  //  std::ios_base::binary     
                    if ( !in_fs.good()  )
                    {
                        std::stringstream error_strm;
                        print_line_of_markers("ERROR! ", &error_strm);
                        print_line_of_markers("(", &error_strm);
                        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                        
                        error_strm << "\n\nERROR!   not good ofstream   in  \"load_posterior_samples_from_file\"\n\t\t load_samples_filename =    "  
                        <<  load_samples_filename  << "\n\n";
                        
                        print_line_of_markers(")", &error_strm);
                        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
                    }
                                                        
    boost::archive::binary_iarchive in_archive( in_fs );
    
    
    type_vector_map_uint_to_Sampled_diploid_Event_data    samples_from_posterior;
    
    in_archive >> samples_from_posterior;            
    in_fs.close(); 
    
    
    return  samples_from_posterior;

} // load_posterior_samples_from_file




type_vector_map_uint_to_Sampled_diploid_Event_data load_posterior_samples_from_file
                        (const std::string &the_fname)
{            
    std::ifstream in_fs(the_fname,   std::ios::binary | std::ios::in  );  //  std::ios_base::binary     
                    if (!in_fs.good())
                    {
                        std::stringstream error_strm;                                               
                        error_strm << "\n\nERROR!   not good ifstream   in  \"load_posterior_samples_from_file\"\n\t\t the_fname =    "   <<  the_fname  << "\n\n";
                        error_message(error_strm, true);
                    }
                                                        
    boost::archive::binary_iarchive in_archive(in_fs);
        
    type_vector_map_uint_to_Sampled_diploid_Event_data    samples_from_posterior;
    
    in_archive >> samples_from_posterior;            
    in_fs.close(); 
    
    
    return  samples_from_posterior;

} // load_posterior_samples_from_file
































void save_marginal_posterior_distributions_to_file
                        (const type_map_uint_to_Marginal_Event_posterior_distribution &marginal_distributions,
                         const std::string &some_output_dir)
{    

    std::fprintf(stderr, "\nsave marginals\n");

    char save_marginals_filename[option_size];
    std::sprintf(save_marginals_filename, "%s/marginal_distributions_via_posterior_sample",
                                        some_output_dir.c_str());
    
    FILE *marginals_file = std::fopen(save_marginals_filename, "a");
                                if (marginals_file == NULL)
                                {
                                    std::fprintf(stderr, "ERROR: unable to open file [%s] in \"save marginals\".\n\n",
                                                save_marginals_filename );
                                    return;
                                }        
    
    
    std::fprintf(marginals_file, "\n\n\n");
    
    for (type_map_uint_to_Marginal_Event_posterior_distribution::const_iterator it_ev = marginal_distributions.begin();
            it_ev != marginal_distributions.end();
            ++it_ev)
    {
        std::fprintf(marginals_file, "\nEvent %u:", it_ev->first);
        
        for (type_map__hap_outcome__hap_outcome__to__longdouble::const_iterator it_outcome = it_ev->second.map_diploid_outcomes_to_posterior_probability.begin();
                it_outcome != it_ev->second.map_diploid_outcomes_to_posterior_probability.end();
                ++it_outcome)
        {                    
            std::fprintf(marginals_file, "\n\t%u,%u\t\t%1.20Lf\n",
                                        it_outcome->first.first, it_outcome->first.second, 
                                        it_outcome->second);
                                        
            for (type_map__BI__BI__to__longdouble::const_iterator it_brkpts 
				= it_ev->second.map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability.at(it_outcome->first).begin();
                    it_brkpts !=  it_ev->second.map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability.at(it_outcome->first).end();
                    ++it_brkpts)                             
                std::fprintf(marginals_file, "\t\t\t\t\t\t[%u,%u],[%u.%u]\t\t%1.20Lf\n",
                                            it_brkpts->first.first.lower(),
					    it_brkpts->first.first.upper(),
                                            it_brkpts->first.second.lower(),
					    it_brkpts->first.second.upper(),			     
                                            it_brkpts->second);
        }                               
    }
    
    
    std::fclose(marginals_file);  
        

} // save_marginal_posterior_distributions_to_file
























void save_calls_tables_to_file
                    (const std::stringstream &output_marginal_Positive_calls_strm,
                    const std::stringstream &output_marginal_Negative_calls_strm,
                    const std::stringstream &output_marginal_all_calls_strm,
                     const std::string &some_output_dir)
{

    char table_filename[option_size];
    
    
    //Positive
    {
        std::sprintf(table_filename, "%s/Positive_calls", 
                                    some_output_dir.c_str()  );      
        std::ofstream  positive_ofstrm(   table_filename ,   std::ios::out  |   std::ios::app);     
        
        if ( !positive_ofstrm.good() )
        {
            std::stringstream error_strm;
            print_line_of_markers("ERROR! ", &error_strm);
            print_line_of_markers("(", &error_strm);
            error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
            
            error_strm << "\n\nERROR!   not good ofstream   in  \"table_filename\"  in  save positive calls\n\t\t table_filename =    ["  
            <<  table_filename  << "]\n\n";
            
            print_line_of_markers(")", &error_strm);
            std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
        }
        else
        {
	    positive_ofstrm.precision(10);
            positive_ofstrm << output_marginal_Positive_calls_strm.str();                     
            positive_ofstrm.close();            
        }
    }
    
    
    
    
    //Negative
    {
        std::sprintf(table_filename, "%s/Negative_calls", 
                                    some_output_dir.c_str()  );      
        std::ofstream  negative_ofstrm(   table_filename ,   std::ios::out  |   std::ios::app);     
        
        if ( !negative_ofstrm.good() )
        {
            std::stringstream error_strm;
            print_line_of_markers("ERROR! ", &error_strm);
            print_line_of_markers("(", &error_strm);
            error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
            
            error_strm << "\n\nERROR!   not good ofstream   in  \"table_filename\"  in  save negative calls\n\t\t table_filename =    ["  
            <<  table_filename  << "]\n\n";
            
            print_line_of_markers(")", &error_strm);
            std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
        }
        else
        {
            negative_ofstrm << output_marginal_Negative_calls_strm.str();                     
            negative_ofstrm.close();            
        }
    }    
    
    
    
    
    
    //All
    {
        std::sprintf(table_filename, "%s/All_calls", 
                                    some_output_dir.c_str()  );      
        std::ofstream  all_ofstrm(   table_filename ,   std::ios::out  |   std::ios::app);     
        
        if ( !all_ofstrm.good() )
        {
            std::stringstream error_strm;
            print_line_of_markers("ERROR! ", &error_strm);
            print_line_of_markers("(", &error_strm);
            error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
            
            error_strm << "\n\nERROR!   not good ofstream   in  \"table_filename\"  in  save All calls\n\t\t table_filename =    ["  
            <<  table_filename  << "]\n\n";
            
            print_line_of_markers(")", &error_strm);
            std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
        }
        else
        {
            all_ofstrm << output_marginal_all_calls_strm.str();                     
            all_ofstrm.close();            
        }
    }        
            


} // save_calls_tables








































type_set_uint  read_totally_excluded_Events_from_previous_job_save_dir
			    (const std::string &old_run_data_dir)
{
    
    type_set_uint  excluded_Events;  
    
    
    char exclu_filename[option_size];
    std::sprintf(exclu_filename, "%s/Totally_excluded_Events", old_run_data_dir.c_str() );    
    
    
    FILE *exclu_file = std::fopen(exclu_filename, "r");
    
    
    if ( exclu_file == NULL  )
    {
	std::stringstream error_strm;
	print_line_of_markers("ERROR! ", &error_strm);
	print_line_of_markers("(", &error_strm);
	error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
	
	error_strm << "\n\nERROR!   not good fopen   in  \"read_totally_excluded_Events_from_previous_job_save_dir\",\n\t\t exclu_filename =    ["  
	<<  exclu_filename  << "]\n\n";
	
	print_line_of_markers(")", &error_strm);
	std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );  
    } 
    else
    {    	
	uint in_excluded_event_uid;    
	while (  std::fscanf(exclu_file, "%u", &in_excluded_event_uid)  !=  EOF)
	{
		excluded_Events.insert(in_excluded_event_uid);
		std::fscanf(exclu_file, "\n");    
	}
	
	std::fclose(exclu_file);
    }
            
    
    return excluded_Events;

} // read_totally_excluded_Events_from_previous_run_diagnostic_output_file


























void remove_directories_associated_with_Partial_Sum_functions_from_scratch_space()
{
    if (my_MPI_rank == 0)
    {            
	//delete Partial Sum Function scratch folders
	const boost::filesystem::directory_iterator it_dir_end;
	for (boost::filesystem::directory_iterator it_dir(scratch_job_output_dir);
		it_dir != it_dir_end;
		++it_dir)
	{
	    int dummy_val;
	    const int sscanf_return_val = std::sscanf(   (*it_dir).path().filename().c_str(), "%d", &dummy_val);
	    
	    if (  boost::filesystem::is_directory( (*it_dir).path() )
		and    sscanf_return_val > 0    )
	    {
		boost::filesystem::remove_all(  (*it_dir).path()  );	    
	    }	    
	}
    } // 0

} // remove_directories_associated_with_Partial_Sum_functions_from_scratch_space





void remove_MPI_rank_directories_from_scratch_space()
{
    if (my_MPI_rank == 0)
    {    
	//delete rank folders
	for (int rank = 0; rank < size_MPI_WORLD; ++rank)
	{
	    char rank_path_name[option_size];
	    std::sprintf(rank_path_name, "%s/rank_%d",
					scratch_job_output_dir.c_str(),
					my_MPI_rank  );
	    
	    boost::filesystem::remove_all( boost::filesystem::path( rank_path_name )  );
	}		
    } // rank 0
    
} // remove_MPI_rank_directories_from_scratch_space














void create_scratch_directory_tree()
{
 
    char sub_out_dirs[option_size]; 

    if (my_MPI_rank == 0)
    {
	std::sprintf(sub_out_dirs, "%s",
				scratch_job_output_dir.c_str() );        
	boost::filesystem::create_directory(sub_out_dirs);
    }
    
    world_MPI_boost_ptr->barrier();
    
    

    std::sprintf(sub_out_dirs, "%s/rank_%d",
			    scratch_job_output_dir.c_str(), my_MPI_rank );        
    boost::filesystem::create_directory(sub_out_dirs);    
    
    
    for (uint thread_ctr = 0; thread_ctr < max_number_of_OMP_threads; ++thread_ctr)
    {
	std::sprintf(sub_out_dirs, "%s/rank_%d/thread_%u",
				scratch_job_output_dir.c_str(), 
				my_MPI_rank, thread_ctr);        
	boost::filesystem::create_directory(sub_out_dirs);                
    }

} // create_scratch_directory_tree







void create_output_directory_tree()
{
    
    char sub_out_dirs[option_size];     	
    
    if (my_MPI_rank == 0)
    {
	std::sprintf(sub_out_dirs, "%s", output_dir.c_str()   );        
	boost::filesystem::create_directory(sub_out_dirs);       
	
	std::sprintf(sub_out_dirs, "%s/Partial_sum_functions", output_dir.c_str());
	boost::filesystem::create_directory(sub_out_dirs);
	
	std::sprintf(sub_out_dirs, "%s/Star_alignments", output_dir.c_str());
	boost::filesystem::create_directory(sub_out_dirs);            
	
	std::sprintf(sub_out_dirs, "%s/Samples_from_posterior", output_dir.c_str());
	boost::filesystem::create_directory(sub_out_dirs);     
	        
	std::sprintf(sub_out_dirs, "%s/Likelihoods",  output_dir.c_str()  );     
	boost::filesystem::create_directory(sub_out_dirs);
	
	std::sprintf(sub_out_dirs, "%s/Read_Depth",  output_dir.c_str()  );     
	boost::filesystem::create_directory(sub_out_dirs);	
    }  
    
    
    world_MPI_boost_ptr->barrier();

} // create_output_directory_tree





























void save_excluded_Events_to_file
		    (const type_set_uint &excluded_Events)
{

    if (   !excluded_Events.empty()   )
    {
	
	char excluded_filename[option_size];
	std::sprintf(excluded_filename, "%s/Totally_excluded_Events", output_dir.c_str() );
	
	FILE *exclu_file = std::fopen(excluded_filename, "a");
	if (exclu_file == NULL)
	{
	    std::stringstream error_strm;
	    print_line_of_markers("ERROR! ", &error_strm);
	    print_line_of_markers("(", &error_strm);
	    error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";	    
	    error_strm << "\n\nERROR!  unable to open file ["  <<   excluded_filename  << "] in \"save totally excluded Events to file\".\n\n";	    
	    print_line_of_markers(")", &error_strm);
	    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                                   
	}	
	else
	{	
	    for(type_set_uint::const_iterator it_exclu = excluded_Events.begin();
		    it_exclu != excluded_Events.end();
		    ++it_exclu)	
		std::fprintf(exclu_file, "%u\n", *it_exclu);
	    
	    std::fclose(exclu_file);
	}	 
	
    } //non-empty
    
    
} // save_excluded_Events_to_file
















void save_remaining_original_breakpoints_to_file
		(const type_map_uint_to_Event  &the_events_of_some_Conn_Comp)
{		
    
//     char orig_brkpts_filename[option_size];
//     std::sprintf(orig_brkpts_filename, "%s/remaining_variational_positions_after_heuristic", output_dir.c_str() );
//     
//     std::ofstream orig_brkpts_file(orig_brkpts_filename, std::ios_base::app);
//     
//     if (check_fstream_for_goodness(orig_brkpts_file, orig_brkpts_filename, "save_remaining_original_breakpoints_to_file", false))
//     {	
// 	for (type_map_uint_to_Event::const_iterator it_ev = the_events_of_some_Conn_Comp.begin();
// 		it_ev != the_events_of_some_Conn_Comp.end();
// 		++it_ev)
// 	{
// 	    orig_brkpts_file << it_ev->first << ":\n";
// 	    orig_brkpts_file << it_ev->second.my_original_Breakpoint_complex.variational_positions__ALL.size() << "\t"
// 			    << it_ev->second.my_original_Breakpoint_complex.breakpoints_NAHR__AB.size() << "\t"
// 			    << it_ev->second.my_original_Breakpoint_complex.breakpoints_NAHR__BA.size() << "\t"
// 			    << it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___ABA.size() << "\t"
// 			    << it_ev->second.my_original_Breakpoint_complex.breakpoints_GeneConversion___BAB.size() << "\n";
// 	    
// 	    write_set_to_file<uint>();
// 	    std::fprintf(orig_brkpts_file, "%u\t%u\n", it_ev->second.UID, it_ev->second.original_variational_positions.size() );
// 	    for (type_set_uint::const_iterator it_orig = it_ev->second.original_variational_positions.begin();
// 		    it_orig != it_ev->second.original_variational_positions.end();
// 		    ++it_orig)
// 		std::fprintf(orig_brkpts_file, "%u\t", *it_orig);
// 	    
// 	    std::fprintf(orig_brkpts_file, "\n");	    
// 	}//ev
// 
// 	std::fclose(orig_brkpts_file);
//     }
    
}//save_remaining_original_breakpoints_to_file














// void load_original_breakpoints_after_heuristic_for_all_Events_in_Connected_Component_from_file
// 		    (const std::string &old_run_data_dir,
// 		     type_map_uint_to_Event  &events_of_some_CC)
// {		    
// 
//     char orig_brkpts_filename[option_size];
//     std::sprintf(orig_brkpts_filename, "%s/remaining_variational_positions_after_heuristic", old_run_data_dir.c_str() );
//     
//     FILE *orig_brkpts_file = std::fopen(orig_brkpts_filename, "r");
//     if (orig_brkpts_file == NULL)
//     {
// 	std::stringstream error_strm;
// 	print_line_of_markers("ERROR! ", &error_strm);
// 	print_line_of_markers("(", &error_strm);
// 	error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";	    
// 	error_strm << "\n\nERROR!  unable to open file ["  <<   orig_brkpts_filename  << "] in \"load orig_brkpts from file\".\n\n";	    
// 	print_line_of_markers(")", &error_strm);
// 	std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() ); 
// 	exit(1);
//     }	
//     else
//     {	
// 	uint in_ev_uid;
// 	uint in_num_orig_brkpts;
// 	
// 	while (std::fscanf(orig_brkpts_file, "%u\t%u\n", &in_ev_uid, &in_num_orig_brkpts)  != EOF)
// 	{
// 	    type_set_uint in_orig_brkpts;
// 	    
// 	    for (uint b=0; b < in_num_orig_brkpts; ++b)
// 	    {
// 		uint in_brkpt;
// 		std::fscanf(orig_brkpts_file, "%u\t", &in_brkpt );
// 		in_orig_brkpts.insert(in_brkpt);
// 	    }
// 	    
// 	    std::fscanf(orig_brkpts_file, "\n");
// 	    
// 	    if (events_of_some_CC.count(in_ev_uid) > 0 )
// 	    {
// 		events_of_some_CC.at(in_ev_uid).original_variational_positions = in_orig_brkpts;
// 		events_of_some_CC.at(in_ev_uid).variational_positions_of_profile = in_orig_brkpts;
// 	    }		    		    
// 	}
// 	
// 	std::fclose(orig_brkpts_file);
//     }
//     
//     
// 
// }//load remaining breakppoints    
// 
// 
























void  save_Events_affected_by_Undefined_regions_to_file
		(const type_map_uint_to_Event  &events_of_some_CC,
		 const type_map_uint_to_list_BI &ignored_regions_fo_genome)
{
        
    type_set_uint Events_affected_by_Undefined_regions;
    
    for (type_map_uint_to_Event::const_iterator it_ev = events_of_some_CC.begin();
	    it_ev != events_of_some_CC.end();
	    ++it_ev)
    {
	for (type_map_uint_to_list_BI::const_iterator it_ignore_chr = ignored_regions_fo_genome.begin();
		it_ignore_chr != ignored_regions_fo_genome.end();
		++it_ignore_chr)
	{
	    if (it_ignore_chr->first == it_ev->second.chromos[0])
	    {
		for (type_list_BI::const_iterator it_ignore_bi = it_ignore_chr->second.begin();
			it_ignore_bi != it_ignore_chr->second.end();
			++it_ignore_bi)
		{
		    if (  BOOST_overlap(*it_ignore_bi, it_ev->second.region_between_and_including_the_LCRs_themselves)  )
		    {
			Events_affected_by_Undefined_regions.insert( it_ev->second.UID );
			break;
		    }			    
		}//BI	
	    }//same chromo		    
	}//ignore chr		
    }//ev		
    
        
    
    char undefined_fname[option_size];
    std::sprintf(undefined_fname, "%s/Events_and_Connected_Components_involving_Undefined_regions_of_the_genome",
				    output_dir.c_str());
				    
    FILE *undefined_file = std::fopen(undefined_fname, "a");
    if (undefined_file == NULL)
    {
	std::stringstream error_strm;
	print_line_of_markers("ERROR! ", &error_strm);
	print_line_of_markers("(", &error_strm);
	error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";	    
	error_strm << "\n\nERROR!  unable to open file undefined_fname = ["  <<   undefined_fname
		    << "] in \"save_Events_affected_by_Undefined_regions_to_file\".\n\n";	    
	print_line_of_markers(")", &error_strm);
	std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() ); 		    
    }	  
    else
    {    
	for (type_set_uint::const_iterator it_aff = Events_affected_by_Undefined_regions.begin();
		it_aff != Events_affected_by_Undefined_regions.end();
		++it_aff)
	    std::fprintf(undefined_file, "%u\t%u\n", *it_aff, events_of_some_CC.at(*it_aff).my_Conn_Comp->Conn_Comp_ID);
	
	std::fclose(undefined_file);   
    }


} // save_Events_affected_by_Undefined_regions_to_file



















type_map_uint_to_list_BI load_Non_unique_regions_of_genome()
{
    
    type_map_uint_to_list_BI nonuniq_regions;
    
    
    char nonuniq_fname[option_size];
    std::sprintf(nonuniq_fname, "%s/Non_unique_regions_of_the_genome", data_directory.c_str()  );
    
    FILE *nonuniq_file = std::fopen(nonuniq_fname, "r");
                    if (nonuniq_file == NULL)
                    {
			std::stringstream error_strm;
			print_line_of_markers("ERROR! ", &error_strm);
			print_line_of_markers("(", &error_strm);
			error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";	    
			error_strm << "\n\nERROR!  unable to open file nonuniq_fname = ["  <<   nonuniq_fname
				    << "] in \"load_Non_unique_regions_of_genome\".\n\n";	    
			print_line_of_markers(")", &error_strm);
			std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
			exit(1);
                    }            
    
    
    uint in_chr;
    uint in_size;
    
    while (std::fscanf(nonuniq_file,  "chr %u\tsize %u\n",  &in_chr, &in_size ) != EOF  )
    {
	for (uint bi=0; bi < in_size; ++bi)
	{
	    uint in_low, in_upp;
	    std::fscanf(nonuniq_file, "\t[%u, %u]\n", &in_low, &in_upp  );
	    
	    nonuniq_regions[in_chr].push_back( BOOST_Interval(in_low, in_upp)  );	
	}        
    } // chr,size
        
    
    std::fclose(nonuniq_file);
    
    
    
    return nonuniq_regions;

} // load_Non_unique_regions_of_genome













type_set_uint upload_connected_components_involving_Undefined_regions_of_genome
			(const std::string &the_data_dir)
{
    char undef_fname[option_size];
    std::sprintf(undef_fname, "%s/Connected_Components_involving_Undefined_regions_of_genome", the_data_dir.c_str() );
    
    
    FILE *undef_file = std::fopen(undef_fname, "r");
                                    if (undef_file == NULL)
                                    {
                                        std::fprintf(stderr, "ERROR: unable to open file undef_fname [%s] in \"upload_connected_components_involving_Undefined_regions_of_genome\".\n\n",
                                                    undef_fname );
                                        exit(1);
                                    }
    
    type_set_uint undef_CCs;
    uint in_CC;
    while (std::fscanf(undef_file, "%u", &in_CC) != EOF)
    {
	undef_CCs.insert(in_CC);
	std::fscanf(undef_file, "\n");
    }
        
    std::fclose(undef_file);
    
    
    return undef_CCs;            
    
} // upload_connected_components_involving_Undefined_regions_of_genome















type_set_string   read_genome_list_from_file
		(const std::string  &file_containing_genome_names)
{
    type_set_string genome_names;
    
    FILE *genome_file = std::fopen(file_containing_genome_names.c_str(), "r");
                                    if (genome_file == NULL)
                                    {
                                        std::fprintf(stderr, "ERROR: unable to open file [%s] in \"read_genome_list_from_file\".\n\n",
                                                    file_containing_genome_names.c_str() );
                                        exit(1);
                                    }
    char in_genome_name[option_size];                         
    while (std::fscanf(genome_file, "%s", in_genome_name) != EOF)
    {
	genome_names.insert(in_genome_name);    
    }
    
    std::fclose(genome_file);
                  
    
    return genome_names;

}//read_genome_list_from_file












































void append_diploid_Sampled_Events_to_file
		(const type_map_uint_to_Sampled_diploid_Event_data &sampled_Events,
		const std::string  &outdir)
{
    std::string sampled_fname(outdir);
    sampled_fname.append("/raw_posterior_marginal_MAP_calls");
    
    std::ofstream  outsamples(sampled_fname, std::ios::app);
    if (check_fstream_for_goodness(outsamples, sampled_fname, "append_diploid_Sampled_Events_to_file", false))
    {
	for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_samp = sampled_Events.begin();
		it_samp != sampled_Events.end();
		++it_samp)
	{
	    outsamples 
		    << it_samp->second.associated_genome << "\t"
		    << it_samp->second.cc_of_event << "\t"
		    << it_samp->second.event_UID << "\t"
		    
		    << convert_haploid_outcome_to_string(it_samp->second.the_diploid_Event_outcome.first) << "\t"
		    << it_samp->second.the_diploid_profile_brkpts.first.lower() << "\t"
		    << it_samp->second.the_diploid_profile_brkpts.first.upper() << "\t"
		    
		    << convert_haploid_outcome_to_string(it_samp->second.the_diploid_Event_outcome.second) << "\t"
		    << it_samp->second.the_diploid_profile_brkpts.second.lower() << "\t"
		    << it_samp->second.the_diploid_profile_brkpts.second.upper() << "\t"
		    
		    << it_samp->second.P_diploid_outcome << "\t"
		    << it_samp->second.P_diploid_breakpoint << "\n";	
	}//samp
	
	outsamples.close();
    }       
    
}//append_diploid_Sampled_Events_to_file











type_map_string_to_uint_to_Sampled_diploid_Event_data read_Sampled_Events_from_file
									(const std::string &fname)
{        
    type_map_string_to_uint_to_Sampled_diploid_Event_data  samps_read_in;
    
    std::ifstream  infs(fname);
    {
	std::string ingen;
	
	while (infs >> ingen)
	{		
	    if (!infs.eof()  and  infs.good())
	    {
		Sampled_diploid_Event_data in_samp;
		
		in_samp.associated_genome = ingen;
		infs >> in_samp.cc_of_event;
		infs >> in_samp.event_UID;
		
		for (uint hap=0; hap<2; ++hap)
		{	
		    std::string in_outcome;
		    infs >> in_outcome;
		    pair_at<type_haploid_outcome>(in_samp.the_diploid_Event_outcome,hap) = convert_string_to_haploid_outcome(in_outcome);
		    
		    uint low;
		    uint upp;
		    infs >> low;
		    infs >> upp;
		    pair_at<BOOST_Interval>(in_samp.the_diploid_profile_brkpts,hap) = BOOST_Interval(low,upp);
		}
		
		infs >> in_samp.P_diploid_outcome;
		infs >> in_samp.P_diploid_breakpoint;
		
		if (!infs.eof()  and  infs.good())
		{
		    samps_read_in[in_samp.associated_genome][in_samp.event_UID] = in_samp;	    
		}	   
	    }
	}//while
	
	infs.close();	
	
    }//check_fstream_for_goodness	
	    
    return  samps_read_in;    
    
}//read_Sampled_Events_from_file















void load_cumulative_GC_counts_for_each_chromosome()
{
	for (uint chromo = 1; chromo < 23; ++chromo)
	{
		map_chromosome_to_cumulative_GC_counts[chromo] = type_vector_uint();
		load_cumulative_GC_counts_for_specific_chromosome(chromo, map_chromosome_to_cumulative_GC_counts.at(chromo));
	}//chromo	
	
}//load_cumulative_GC_counts_for_each_chromosome











void load_cumulative_GC_counts_for_specific_chromosome
				(const uint &chromo,
				type_vector_uint &count_vector)
				
{
	count_vector.clear();

	char infnamechar[option_size];
	std::sprintf(infnamechar, "%s/cumulative_GC_count_chromosome_%u", data_directory.c_str(), chromo);
	std::ifstream infs(infnamechar);
	boost::archive::text_iarchive inarch(infs);
		
	inarch >> count_vector;
	
	infs.close();
	
}//load_cumulative_GC_counts_for_specific_chromosome











type_map_string_to_string__bool load_gender_and_population_map
					(const std::string &gender_and_pop_filename)
{
	type_map_string_to_string__bool  map_genome_name_to_pop_and_female;
	
	std::ifstream infs(gender_and_pop_filename);
	std::string in_genome;
	std::string in_pop;
	std::string in_gender;
	
	while (infs >> in_pop)
	{
		infs >> in_genome >> in_gender;
		
		const bool is_female = (strcasecmp(in_gender.c_str(), "female") == 0);
		map_genome_name_to_pop_and_female[in_genome] = type_string__bool(in_pop, is_female);
	}//in_genome
	
	infs.close();
	
	return map_genome_name_to_pop_and_female;
	
}//load_gender_and_population_map