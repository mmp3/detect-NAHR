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

#include <limits.h>

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


// #include <gmpfrxx.h>
#include <boost/math/bindings/mpfr.hpp>



#include <mpreal.h>
#include <mpfr.h>





#include <boost/random.hpp>

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamAux.h>



#include <general_typedefs.h>
#include <SD_Entry_general.h>
#include <Interlocking_entry.h>
#include <Coords_and_homologous_region.h>
#include <translations_through_alignment.h>
#include <Sparse_map.h>
#include <templates.h>
#include <Readgroup_statistics.h>


#include "Conn_Comp.h"
#include "Event.h"
#include "globals.h"
#include "other_functions.h"
#include "Paired_end_read.h"
#include "MiniRegion.h"
#include "io_functions.h"

#include <CUPA.h>



void Event::print_this_entry(std::stringstream *const &some_ss) const
{
    
    std::stringstream output_str;
    
    bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str;  
    
    print_line_of_markers("[[[", ss_ptr);    
    
    (*ss_ptr) << "\n\n\n\n\nprinting EVENT:\n\n";
    print_this_entry_general(ss_ptr);    
    (*ss_ptr) << "connected_component = " << connected_component << "\n";
    print_map_keys<uint, Interlocking_entry<Event> >(local_interlocking_Events, "local_interlocking_Events", ss_ptr);
    print_map_keys<uint, Interlocking_entry<Event> >(global_interlocking_Events, "global_interlocking_Events", ss_ptr);   
    
    (*ss_ptr) << "\nregions_homologous_to_directly_affected_region.size() = " << regions_homologous_to_directly_affected_region.size() << "\n";
   
    print_line_of_markers("]]]", ss_ptr);    
    

    if (!provided_with_valid_sstream)
    {  std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );  }
    
    
    std::cerr.flush();std::fflush(stderr);
      
} // end of print_this_entry








static bool min_bool(const bool &val_A, const bool &val_B)
{ 
    if (  !val_A    or   !val_B  )
        return false; 
    else
        return true;  
}
























void Event::partition_DAR_by_homology___and___create_natural_poisson_intervals___NOT__gender_sensitive()
{
    
    type_map_uint_to_set_uint partition_of_DAR_by_homology;
    
    //we will have a collection of sets.  each set will be homologous coordinates; this is obtained via, intra-homology.    
    //then, we go through ALL homology maps (including thre intra-, again), and color SETS of coordinates.    
    //Then, we merge sets with like colors.
        //This strategy relies on the assumptiuon that if reigon  A maps to region B and C, then there exists some other map form B to C.  i.e. we assume closure wrt homology.
        
//     typedef std::list<type_set_uint>  type_list_set_uint;
    //type_list_map_uint_to_set_uint  // list of collections of chromosomes and coordinates.
    
    
    
    
    //"Vertical" identification:  any bases tat are aligned together should be identified together.
    std::cerr << "\n\t\t\tintra partitioning...\n";
    
    type_list_set_uint intra_partition;

    //initialize to all individuals:
    for (type_set_uint::const_iterator it = remaining_coordinates_for_consideration.begin();
            it != remaining_coordinates_for_consideration.end();
            ++it)                
    {
        intra_partition.push_back( type_set_uint()  );  
        intra_partition.back().insert(*it);
    }
    
//     std::fprintf(stderr, "intra_partition.size   =   %u   =?=   %u   =   remaining_coordinates_for_consideration.size()\n",
//                             intra_partition.size(),  remaining_coordinates_for_consideration.size()  );
    
    
//     std::fprintf(stderr, "\"intra\" homol regions...\n");    
    for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg = regions_homologous_to_directly_affected_region.begin();
            it_homol_reg != regions_homologous_to_directly_affected_region.end();
            ++it_homol_reg)
    {
        if (  it_homol_reg->chromosome_of_homologous_region == chromos[0]
              and  BOOST_overlap( get_endpoints_of_keys_of_compressed_map(it_homol_reg->compressed_map_homologous_region_to_profile),
                                  region_between_and_including_the_LCRs_themselves)    )
        {                     
            const type_map_uint_to_uint full_map_homologous_region_to_DAR(              
                                                get_map_intersection_of_map_values_and_set<uint, uint>(
                                                    get_map_intersection_of_compressed_map_keys_and_set(                                
                                                        it_homol_reg->compressed_map_homologous_region_to_profile, //to DAR
                                                        remaining_coordinates_for_consideration),
                                                    remaining_coordinates_for_consideration)    );                                                                                                                                                               
                                                                                                
            for (type_map_uint_to_uint::const_iterator it_map_homol_to_DAR = full_map_homologous_region_to_DAR.begin();
                    it_map_homol_to_DAR != full_map_homologous_region_to_DAR.end();
                    ++it_map_homol_to_DAR)  
            { //merge:
                    type_list_set_uint::iterator found_domain, found_range;
                
                    for (type_list_set_uint::iterator find_it = intra_partition.begin();
                            find_it != intra_partition.end();
                            ++find_it)
                        if (find_it->count(it_map_homol_to_DAR->first) > 0)        
                        {
                            found_domain = find_it;
                            break;
                        }
                    
                    for (type_list_set_uint::iterator find_it = intra_partition.begin();
                            find_it != intra_partition.end();
                            ++find_it)
                        if (find_it->count(it_map_homol_to_DAR->second) > 0)        
                        {
                            found_range = find_it;
                            break;
                        }      
                        
                    if (found_domain != found_range)    // i.e. not already merged.
                    {
                        found_domain->insert(found_range->begin(), found_range->end());
                        intra_partition.erase(found_range);
                    }
            }  //homol to DAR                    
        }        
    }
    
    
    
    
    
    
    
    
    std::cerr << "\n\t\t\tcolor intra partition...\n";
    
    //NOW you can color.
    
    typedef  std::pair<type_set_uint, uint>   type_set_uint____uint;
    typedef  std::list<type_set_uint____uint>  type_list___set_uint____uint;
    
    type_list___set_uint____uint coloring_of_sets;                                                            
    type_set_uint current_colors_in_use;
    
    type_list___set_uint____uint::iterator insert_coloring_of_coordinates_it = coloring_of_sets.begin();
    type_set_uint::iterator insert_color_in_use = current_colors_in_use.begin();
    
    insert_color_in_use = current_colors_in_use.insert(insert_color_in_use, 0);
       
    //transfer and initialize color:
    while ( !intra_partition.empty() )              
    {
        insert_coloring_of_coordinates_it = coloring_of_sets.insert( insert_coloring_of_coordinates_it,
                                                                     type_set_uint____uint( *(intra_partition.begin()), 0)    ); 
        intra_partition.erase(intra_partition.begin());                                                                    
    }

     
     
     //color by homology

//     std::fprintf(stderr, "homol regions...  \n");
    for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg = regions_homologous_to_directly_affected_region.begin();
            it_homol_reg != regions_homologous_to_directly_affected_region.end();
            ++it_homol_reg)
    {
        const uint newest_color = sum_over_set<uint>(current_colors_in_use) + 1;
                
        type_map_uint_to_uint full_map_DAR_to_homologous_region(                                         
                                    get_map_intersection_of_compressed_map_keys_and_set(        
                                            get_compressed_map_inverse_of_compressed_map(it_homol_reg->compressed_map_homologous_region_to_profile), //DAR to homol
                                            remaining_coordinates_for_consideration)  );
                                     
                                        
        while( !full_map_DAR_to_homologous_region.empty() )     
        {
            bool found_a_set = false;
            
            for (type_list___set_uint____uint::iterator find_it = coloring_of_sets.begin();
                    find_it != coloring_of_sets.end();
                    ++find_it)       
	    {
                if ( find_it->first.count(full_map_DAR_to_homologous_region.begin()->first) > 0 )
                {
                    found_a_set = true;
                    find_it->second += newest_color;
                    insert_color_in_use = current_colors_in_use.insert( insert_color_in_use,  find_it->second  );
                    
                    for (type_set_uint::const_iterator it_exist = find_it->first.begin();
                            it_exist != find_it->first.end();
                            ++it_exist) 
		    {
                        full_map_DAR_to_homologous_region.erase(*it_exist);                        
		    }
                    break;
                }
	    }
                
            if (!found_a_set)
            {
                if (my_MPI_rank == 0)
                {
                    std::stringstream error_str_strm;
                    print_line_of_markers("ERROR! ", &error_str_strm);
                    print_line_of_markers("(", &error_str_strm);          
		    char color_char[option_size];
                    std::sprintf(color_char, "in \"create_natural poisson intervals\"...  unable to place full_map_DAR_to_homologous_region.begin()->first = %u    (iterator in full:  [%u,%u]    in colored sets!!!\n\n\n\n",
                                full_map_DAR_to_homologous_region.begin()->first,
                                full_map_DAR_to_homologous_region.begin()->first, full_map_DAR_to_homologous_region.begin()->second);
                    error_str_strm <<  color_char;           
                         
		    
                    for (type_list___set_uint____uint::iterator find_it = coloring_of_sets.begin();
                            find_it != coloring_of_sets.end();
                            ++find_it)    
                    {                         
                        std::sprintf(color_char, "color %u:", find_it->second);
                        print_set<uint>(find_it->first,color_char, &error_str_strm);
                    }
                                
                    print_line_of_markers(")", &error_str_strm);
                    std::fprintf(stderr, "\n\n\n\n%s\n\n\n\n", error_str_strm.str().c_str()  );
                }
                full_map_DAR_to_homologous_region.erase(full_map_DAR_to_homologous_region.begin());                       
            } // error
            
        }
    } //  regions_homologous_to_directly_affected_region
    
    
    
    
    
    
    
    

    
    std::cerr << "\n\t\t\tHorizontal Identification...\n";
    
    //"Horizontal Identification"  -  based on Events ONLY.
    {//Horizontal
	std::vector< type_list___set_uint____uint::iterator > coloring_of_sets___loop(
					get_loopable_iterators< type_list___set_uint____uint >( coloring_of_sets ));
	const uint number_of_coloring_sets = coloring_of_sets___loop.size();
	
	
	type_set_uint event_UIDs_for_coloring(  extract_keys_of_map_and_return_as_set<uint, Interlocking_entry<Event> >(local_interlocking_Events));    
	event_UIDs_for_coloring.insert(UID);
	
	for (type_set_uint::const_iterator it_uid = event_UIDs_for_coloring.begin();
		it_uid != event_UIDs_for_coloring.end();
		++it_uid)
	{
	    const type_map_uint_to_Event::const_iterator it_the_ev = my_Conn_Comp->events.find(*it_uid);
	    assert(it_the_ev !=  my_Conn_Comp->events.end());
	    for (uint qs = 0; qs < 2;  ++qs)
	    {
		std::cerr << "\n\n\tuid = " << *it_uid << "...";
		std::cerr << "\t\tLOWER...";
		{//LOWER                
		    const uint newest_color = sum_over_set<uint>(current_colors_in_use) + 1;  		    		    
				
		    #pragma omp parallel
		    {
			type_set_uint  thread___new_colors_in_use;                   
					    
			#pragma omp for nowait schedule(dynamic,500)
			for (uint ec=0; ec < number_of_coloring_sets; ++ec)
			{
			    bool has_a_pos_that_is_lower = false;
			    for (type_set_uint::const_iterator it_ec_pos = coloring_of_sets___loop[ec]->first.begin();
				    it_ec_pos != coloring_of_sets___loop[ec]->first.end();
				    ++it_ec_pos)
			    {
				if (*it_ec_pos < it_the_ev->second.LCRs[qs].lower())
				{
				    has_a_pos_that_is_lower = true;
				    break;
				}				
			    }
			    
			    if (has_a_pos_that_is_lower)
			    {
				coloring_of_sets___loop[ec]->second += newest_color;
				thread___new_colors_in_use.insert(  coloring_of_sets___loop[ec]->second  );
			    }			    
			}//ec								
				
			#pragma omp critical (insert_new_colors_this_thread_for_unique_regions_partition_construction__LOWER)                      
			current_colors_in_use.insert(  thread___new_colors_in_use.begin(),  thread___new_colors_in_use.end()  );                        
		    
			mpfr_free_cache();
		    }//parallel                   
		}//LOWER
		
		
		std::cerr << "\t\tUPPER...";
		
		{//UPPER                
		    const uint newest_color = sum_over_set<uint>(current_colors_in_use) + 1;	    
				
		    #pragma omp parallel
		    {
			type_set_uint  thread___new_colors_in_use;                   
					    
			#pragma omp for nowait schedule(dynamic,500)
			for (uint ec=0; ec < number_of_coloring_sets; ++ec)
			{
			    bool has_a_pos_that_is_greater = false;
			    for (type_set_uint::const_iterator it_ec_pos = coloring_of_sets___loop[ec]->first.begin();
				    it_ec_pos != coloring_of_sets___loop[ec]->first.end();
				    ++it_ec_pos)
			    {
				if (*it_ec_pos > it_the_ev->second.LCRs[qs].upper())
				{
				    has_a_pos_that_is_greater = true;
				    break;
				}				
			    }
			    
			    if (has_a_pos_that_is_greater)
			    {
				coloring_of_sets___loop[ec]->second += newest_color;
				thread___new_colors_in_use.insert(  coloring_of_sets___loop[ec]->second  );
			    }			    
			}//ec								
				
			#pragma omp critical (insert_new_colors_this_thread_for_unique_regions_partition_construction__UPPER)                      
			current_colors_in_use.insert(  thread___new_colors_in_use.begin(),  thread___new_colors_in_use.end()  );                        
		    
			mpfr_free_cache();
		    }//parallel                   
		}//UPPER	    
		
		
		
		
		std::cerr << "\t\tINBETWEEN...";
		
		assert( ((qs + 1) % 2) < 2);
		assert(  !it_the_ev->second.compressed_map_LCR_to_profile__for_each_LCR.empty()  );
		
		{//INBETWEEN	    	    
		    type_set_uint  insertions_relative_to_other_LCR;	    
		    {
			
			type_map_uint_to_uint map_profile_to_thiiiis_LCR(	    
							    convert_compressed_map_to_full_map(
							    get_compressed_map_inverse_of_compressed_map( 
								    it_the_ev->second.compressed_map_LCR_to_profile__for_each_LCR.at(qs) )  ));	
								    
			const  type_map_uint_to_uint map_profile_to_OTHER_LCR(
							convert_compressed_map_to_full_map(
							get_compressed_map_inverse_of_compressed_map( 
								it_the_ev->second.compressed_map_LCR_to_profile__for_each_LCR.at( (qs + 1) % 2) )  ));

						
			for (type_map_uint_to_uint::const_iterator it_del = map_profile_to_OTHER_LCR.begin();
				it_del != map_profile_to_OTHER_LCR.end();
				++it_del)
			    map_profile_to_thiiiis_LCR.erase(it_del->first);
			

			insertions_relative_to_other_LCR = extract_keys_of_map_and_return_as_set<uint,uint>( get_inverse_of_map(map_profile_to_thiiiis_LCR) );
		    }
		    
			
		    uint newest_color = sum_over_set<uint>(current_colors_in_use) + 1;
		    
		    #pragma omp parallel for schedule(dynamic,500)
		    for (uint ec=0; ec < number_of_coloring_sets; ++ec)
		    {
			for (type_set_uint::const_iterator it_ins = insertions_relative_to_other_LCR.begin();
				it_ins != insertions_relative_to_other_LCR.end();
				++it_ins)
			{
			    if (  coloring_of_sets___loop[ec]->first.count(*it_ins) > 0  )
			    {
				#pragma omp critical (add_color_to_coloring_of_sets_and_current_colors_in_use_when_looping_over_insertions_of_Events_LCRs_in_partition_DAR)
				{
				    coloring_of_sets___loop[ec]->second += newest_color;
				    current_colors_in_use.insert( coloring_of_sets___loop[ec]->second );
				    newest_color = sum_over_set<uint>(current_colors_in_use) + 1;
				}
				break;			
			    }		    
			}//ins
		    }//ec											    
		}//INBETWEEN	    	    
	    
	    
	    
	    }//qs	        
	}//uid
    }//Horizontal    

    
    
    
    
    
    
    
    
    std::cerr << "\n\t\t\tcreating intervals according to like colors...\n";   
   

    type_map_uint_to_set_uint::iterator natural_insert_it = partition_of_DAR_by_homology.begin();  
    uint map_ctr = 0;  //essentially renames the color.        
        
    for (type_set_uint::const_iterator it_color = current_colors_in_use.begin();
            it_color != current_colors_in_use.end();
            ++it_color)
    {
	
        type_set_uint newest_interval;
                
        type_list___set_uint____uint::iterator sets_it = coloring_of_sets.begin();   //erased as we go        
        while (sets_it != coloring_of_sets.end() )        
	{
            if (sets_it->second == *it_color)
            {
                newest_interval.insert( sets_it->first.begin(), sets_it->first.end() );
                sets_it = coloring_of_sets.erase(sets_it);
            }
            else
                ++sets_it;  
	}                    
           
           
        if ( !newest_interval.empty() )   //necessary to check this, since a color may have AT ONE TIME been in use, but since got "added to"
        {
            natural_insert_it = partition_of_DAR_by_homology.insert( natural_insert_it,  type_uint__set_uint(map_ctr, newest_interval) );                                                                                                      
            ++map_ctr;            
        }        
    }
    
    
    
    
    std::cerr << "\n\t\t\tfinally, create_natural_poisson_intervals___NOT_gender_sensitive...\n";    
    
    create_natural_poisson_intervals___NOT_gender_sensitive( partition_of_DAR_by_homology );
    
    
    partition_of_DAR_by_homology.clear();        
    
} //  partition_DAR_by_homology

































void Event::create_natural_poisson_intervals___NOT_gender_sensitive( const type_map_uint_to_set_uint &partition_of_DAR_by_homology)
{        

    //init to the DAR partition (but with chromosomes).
    natural_poisson_intervals___gender_INsensitive.clear();
    type_map_uint_to_uint_to_uint_to_longdouble::iterator init_it_projection = natural_poisson_intervals___gender_INsensitive.begin();
    for (type_map_uint_to_set_uint::const_iterator it_nat = partition_of_DAR_by_homology.begin();
            it_nat != partition_of_DAR_by_homology.end();
            ++it_nat)
    {
        const type_map_uint_to_uint_to_longdouble::iterator  it_to_chromo 
                            = natural_poisson_intervals___gender_INsensitive[it_nat->first].insert(  std::pair<uint, type_map_uint_to_longdouble>
                                                                                ( chromos[0] , type_map_uint_to_longdouble()  )      ).first;
                                                                
                            //(it_nat->first, type_set_pos__total_GC_correction()   )  );
        type_map_uint_to_longdouble::iterator insert_set_it = it_to_chromo->second.begin();
                                                                 
        for(type_set_uint::const_iterator it_pos = it_nat->second.begin();
                it_pos != it_nat->second.end();
                ++it_pos)
            insert_set_it = it_to_chromo->second.insert(insert_set_it, type_uint__longdouble(*it_pos, 0.00)  );                                            
    }            
    


    //cycle through homologous regions and add coordinates from each position via mapping:
    for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg = regions_homologous_to_directly_affected_region.begin();
            it_homol_reg != regions_homologous_to_directly_affected_region.end();
            ++it_homol_reg)
    {
        
        const type_map_uint_to_uint full_map_homologous_region_to_DAR(
                                        convert_compressed_map_to_full_map(it_homol_reg->compressed_map_homologous_region_to_profile)   );
                                        
        for (type_map_uint_to_uint::const_iterator it_map_homol_to_DAR = full_map_homologous_region_to_DAR.begin();
                it_map_homol_to_DAR != full_map_homologous_region_to_DAR.end();
                ++it_map_homol_to_DAR)  
	{
            for (type_map_uint_to_set_uint::const_iterator it_nat = partition_of_DAR_by_homology.begin();
                    it_nat != partition_of_DAR_by_homology.end();
                    ++it_nat)    
	    {
                if (it_nat->second.count(it_map_homol_to_DAR->second) > 0)
                {
                    natural_poisson_intervals___gender_INsensitive.at(it_nat->first)[it_homol_reg->chromosome_of_homologous_region]
                                                                        .insert( type_uint__longdouble(it_map_homol_to_DAR->first, 0.00L));
                    break;                
                }    
	    }//nat
	}//map
    } //  regions_homologous_to_directly_affected_region
    
       
       
       
       
       
       
    if (my_MPI_rank == 0)  
    {
        std::stringstream out_strm;       
        for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_nat = natural_poisson_intervals___gender_INsensitive.begin();
                it_nat != natural_poisson_intervals___gender_INsensitive.end();
                ++it_nat)      
        {     
            uint size_of_nat = 0;
            for (type_map_uint_to_uint_to_longdouble::const_iterator it_chr = it_nat->second.begin();
                    it_chr != it_nat->second.end();
                    ++it_chr)
                size_of_nat += it_chr->second.size();            
                
            out_strm << "\t\t\t\tcolor:  " << it_nat->first 
                     << "\t\tsize:  " << size_of_nat << "\n";
//                     <<"\t\texample coords: ["  << *newest_interval.begin() 
//                     <<", " << *(--newest_interval.end())
//                     <<"]\n";
        }
       
        std::fprintf(stderr, "\n%s\n", out_strm.str().c_str() );           
    }
       
       
       
       
    
               

}  //   create_natural_poisson_intervals___gender_INsensitive___NOT_gender_sensitive





























void Event::calculate_total_HAPLOID_GC_rate_for_every_base_by_summing_over_all_Readgroup_GC_rates_at_each_position()
{
    
    //by "total_GC_ratio", I mean across each Readgroup  (i.e. the corrective ratio is summed over each Readgroup).    
    
    for(type_map_uint_to_uint_to_uint_to_longdouble::iterator it_nat = natural_poisson_intervals___gender_INsensitive.begin();
            it_nat != natural_poisson_intervals___gender_INsensitive.end();
            ++it_nat)
    {
        for (type_map_uint_to_uint_to_longdouble::iterator it_chr = it_nat->second.begin();
                it_chr != it_nat->second.end();
                ++it_chr)
        {
	    calculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP
		    (uploaded_pieces_of_chromos,
		     it_chr->first,
		     it_chr->second,
		     my_Conn_Comp->Conn_Comp_ID,
		     UID,
		     it_nat->first);
            
        }  //chr
    }//nat                                    

} // calculate_total_HAPLOID_GC_rate_for_every_base_by_summing_over_all_Readgroup_GC_rates_at_each_position






















void Event::calculate_base_diploid_fragmentation_rate___haploid_sensitive___gender_sensitive__and_base_cumulative_haploid_fragmentation_rate()
{

    
    base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive.clear();
    type_map_uint_to_longdouble::iterator insert_base_diploid_counts = base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive.begin();       
    
    for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_nat = natural_poisson_intervals___gender_INsensitive.begin();
            it_nat != natural_poisson_intervals___gender_INsensitive.end();
            ++it_nat)
    {
        longdouble total_HAPLOID_rate_for_this_color___autosomal = 0.00L;
        longdouble total_HAPLOID_rate_for_this_color___XXXX = 0.00L;
        longdouble total_HAPLOID_rate_for_this_color___YYYY = 0.00L;
        
        
        for (type_map_uint_to_uint_to_longdouble::const_iterator it_chr = it_nat->second.begin();
                it_chr != it_nat->second.end();
                ++it_chr)
        {	    
	    //sort for numerical accuracy (avoid underflow during adding)	    
	    type_list_longdouble haploid_rates_per_pos = extract_values_of_map_and_return_as_list<uint,longdouble>(it_chr->second);	    	    
	    haploid_rates_per_pos.sort();
	    
	    const longdouble haploid_rate_this_chr = sum_over_container<type_list_longdouble, longdouble>(haploid_rates_per_pos);
            
            switch (it_chr->first)
            {
                case 23:
                    total_HAPLOID_rate_for_this_color___XXXX  +=  haploid_rate_this_chr;
                    break;
                case 24:
                    total_HAPLOID_rate_for_this_color___YYYY  +=  haploid_rate_this_chr;
                    break;
                default:
                    total_HAPLOID_rate_for_this_color___autosomal  +=  haploid_rate_this_chr;
                    break;                                
            };                                   
            
        }//chr
        
            
            
        longdouble total_DIPLOID_rate_for_this_color = (total_HAPLOID_rate_for_this_color___autosomal * 2.00L);
                
        if (gender_of_individual_is_female)
	{  total_DIPLOID_rate_for_this_color +=  (total_HAPLOID_rate_for_this_color___XXXX * 2.00L);  }
        else
        {
            total_DIPLOID_rate_for_this_color += total_HAPLOID_rate_for_this_color___XXXX;
            total_DIPLOID_rate_for_this_color += total_HAPLOID_rate_for_this_color___YYYY;        
        }
                        
        
        insert_base_diploid_counts = base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive
                                                    .insert( insert_base_diploid_counts,
                                                             type_uint__longdouble(it_nat->first, total_DIPLOID_rate_for_this_color)); 
    }//nat    
    
    
    
    
    
    
    
    //cumulative:
    base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.clear();
    for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_nat = natural_poisson_intervals___gender_INsensitive.begin();
            it_nat != natural_poisson_intervals___gender_INsensitive.end();
            ++it_nat)
    {
        for (type_map_uint_to_uint_to_longdouble::const_iterator it_chr = it_nat->second.begin();
                it_chr != it_nat->second.end();
                ++it_chr)
        {	    
	    type_map_uint_to_longdouble cumulative_HAPLOID_base_rate_per_position(it_chr->second);	    
	    for (type_map_uint_to_longdouble::iterator it_pos = ++cumulative_HAPLOID_base_rate_per_position.begin(),
							it_prevpos = cumulative_HAPLOID_base_rate_per_position.begin();
		    it_pos != cumulative_HAPLOID_base_rate_per_position.end();
		    ++it_pos, ++it_prevpos)
	    {
		it_pos->second += it_prevpos->second;
	    }	    
	    	    	   
	    base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals[it_nat->first][it_chr->first] = cumulative_HAPLOID_base_rate_per_position;	    
	}//chr
    }//nat
                   

}   //  calculate_base_diploid_fragmentation_rate___haploid_sensitive































void Event::calculate_observed_read_depth_for_natural_poisson_intervals
						(BamTools::BamReader &my_BAM_reader)
{            
    
    type_map_uint_to_uint::iterator insert_nat_it = observed_read_depth_for_natural_poisson_intervals.begin();
    
    for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_color = natural_poisson_intervals___gender_INsensitive.begin();
            it_color != natural_poisson_intervals___gender_INsensitive.end();
            ++it_color)
    {
        uint counts_for_this_color = 0;
        
        for (type_map_uint_to_uint_to_longdouble::const_iterator it_chr = it_color->second.begin();
                it_chr != it_color->second.end();
                ++it_chr)               
        {	    
	    const type_list_BI positions_as_intervals(
				convert_set_to_list_of_contiguous_intervals(extract_keys_of_map_and_return_as_set<uint, longdouble>(it_chr->second)));
	    
	    for (type_list_BI::const_iterator it_bi = positions_as_intervals.begin();
		    it_bi != positions_as_intervals.end();
		    ++it_bi)
	    {
		my_BAM_reader.SetRegion(create_BAM_region(it_chr->first, *it_bi, 3));		
		counts_for_this_color += count_decent_reads_with_left_endpoints_in_set_region(my_BAM_reader, *it_bi);		
	    }
        }//chr
        
        insert_nat_it = observed_read_depth_for_natural_poisson_intervals.insert(insert_nat_it, type_uint__uint(it_color->first, counts_for_this_color));
    }//color
    
    
             
    if (my_MPI_rank == 0)                                           
    {
        print_map_keys_and_values<uint, uint>(observed_read_depth_for_natural_poisson_intervals,
					      "observed_read_depth_for_natural_poisson_intervals",
					      NULL, true);
    }
            
} //  end of   calculate_observed_read_depth_for_natural_poisson_intervals___MPI_and_OMP






















void Event::adjust_expected_diploid_fragmentation_rate___breakpoint_specific
                                            (type_map_uint_to_longdouble &expected_diploid_fragmentation_rate_per_natural_poisson_interval__for_these_states,
                                             const type_map_uint_to_uint &proposed_haploid_state_vector,
                                             const type_map_uint_to_BI &proposed_haploid_breakpoints)  const
{
        
    for (type_map_uint_to_uint::const_iterator it_occurring_event = proposed_haploid_state_vector.begin();
            it_occurring_event != proposed_haploid_state_vector.end();
            ++it_occurring_event)
    {
	//if does not affect RD, skip.
	if (it_occurring_event->second == hap_outcome__Inv)
	{  continue;  }
	
	
        const type_map_uint_to_Event::const_iterator  occurring_interlocking_event_it = my_Conn_Comp->events.find(it_occurring_event->first);
	
	//Note that NAHR breakpoints may not always be "points" - in the case that we are not summing over part of an LCR and thus breakpoints are intervals.
	//format:   [brkpt on LCR 0,  brkpt on LCR 1]
	
	const type_BI__BI  absolute_affected_regions(
				convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals
					(proposed_haploid_breakpoints.at(occurring_interlocking_event_it->first),
					occurring_interlocking_event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
					static_cast<type_haploid_outcome>(it_occurring_event->second)));
                  
	if (test_if_haploid_outcome_is_GeneConversion(static_cast<type_haploid_outcome>(it_occurring_event->second)))
	{//GeneConv

	    const BOOST_Interval  affected_region__LCR_0(absolute_affected_regions.first.lower(), absolute_affected_regions.first.upper()-1);
	    const BOOST_Interval  affected_region__LCR_1(absolute_affected_regions.second.lower(), absolute_affected_regions.second.upper()-1);    

	    {//erased	    
		const BOOST_Interval *const erased_tract__ptr 
				    =   (it_occurring_event->second == hap_outcome__GeneConv_ABA)    ?
					    &affected_region__LCR_0   :   &affected_region__LCR_1;
					    
		const uint  erased_qs  =  (it_occurring_event->second == hap_outcome__GeneConv_ABA)    ?   0 : 1;					    
	    
		type_map_uint_to_longdouble::iterator it_color_expected = expected_diploid_fragmentation_rate_per_natural_poisson_interval__for_these_states.begin();
		for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_color = base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.begin();
			it_color != base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.end();
			++it_color)
		{
		    const type_map_uint_to_uint_to_longdouble::const_iterator it_find_chr = it_color->second.find(occurring_interlocking_event_it->second.chromos[erased_qs]);
		    if (it_find_chr != it_color->second.end())
		    {
			const longdouble adjust_value = get_adjustment_value_from_cumulative_map(*erased_tract__ptr, it_find_chr->second);
			it_color_expected->second -= adjust_value;
		    }                                
			    
		    ++it_color_expected;                                
		} // it_color
	    }//erased
	    
	    {//copied	  
		const BOOST_Interval *const copied_tract__ptr 
					=   (it_occurring_event->second == hap_outcome__GeneConv_ABA)    ?
						&affected_region__LCR_1   :  &affected_region__LCR_0;

		const uint  copied_qs  =  (it_occurring_event->second == hap_outcome__GeneConv_ABA)    ?   1 : 0;    	    
	    
		type_map_uint_to_longdouble::iterator it_color_expected = expected_diploid_fragmentation_rate_per_natural_poisson_interval__for_these_states.begin(); 
		for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_color = base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.begin();
			it_color != base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.end();
			++it_color)
		{
		    const type_map_uint_to_uint_to_longdouble::const_iterator it_find_chr = it_color->second.find(occurring_interlocking_event_it->second.chromos[copied_qs]);
		    if ( it_find_chr != it_color->second.end() )
		    {
			const longdouble adjust_value = get_adjustment_value_from_cumulative_map(*copied_tract__ptr, it_find_chr->second);
			it_color_expected->second += adjust_value;
		    }
		    
		    ++it_color_expected;                                
		} // it_color
	    }//copied
	
	}//GeneConv
	else
	{//NAHR - dup or del
	
	    //assert(BOOST_is_a_point(proposed_haploid_breakpoints.at(occurring_interlocking_event_it->first)));  // not necessarily a point, as we may have skipped some "hidden" brkpts that are not being summed out at the moment...
	
	    const BOOST_Interval affected_region(absolute_affected_regions.first.lower(), absolute_affected_regions.first.upper()-1);
		    //yes, you do this for both dups and dels.
		    //carefully draw some pictures and think hard.
		    //recall that you don't even get to this point in the code for inversions (or translocs).               
	    assert( !BOOST_empty(affected_region));

	    type_map_uint_to_longdouble::iterator it_color_expected = expected_diploid_fragmentation_rate_per_natural_poisson_interval__for_these_states.begin();         
	    for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_color = base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.begin();
		    it_color != base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.end();
		    ++it_color)
	    {
		const type_map_uint_to_uint_to_longdouble::const_iterator it_find_chr = it_color->second.find(occurring_interlocking_event_it->second.chromos[0]);
		if ( it_find_chr != it_color->second.end() )
		{
		    const longdouble adjust_value = get_adjustment_value_from_cumulative_map(affected_region, it_find_chr->second);
		    if (it_occurring_event->second == hap_outcome__Del)  //deletion    
		    {
			it_color_expected->second -= adjust_value;
		    }
		    else // duplication
		    {
			it_color_expected->second += adjust_value;
		    }
		}                                
			
		++it_color_expected;                                
	    } // it_color
	}//NAHR
                                     
    }  // occuring event                    
    
} // adjust_expected_diploid_fragmentation_rate___breakpoint_specific















































real Event::calculate_probability_of_observed_read_depth_given_expected_read_depth_over_natural_poisson_intervals
                                    (const type_map_uint_to_longdouble &expected_diploid_fragmentation_rate_per_natural_poisson_interval)  const
{
        
    real RD_probability(1);
    
    const mpfr_class the_number_one(1);
        
    type_map_uint_to_longdouble::const_iterator it_expected_rate = expected_diploid_fragmentation_rate_per_natural_poisson_interval.begin();    
    for (type_map_uint_to_uint::const_iterator it_observed = observed_read_depth_for_natural_poisson_intervals.begin();  
            it_observed != observed_read_depth_for_natural_poisson_intervals.end();
            ++it_observed)
    {                
	//For some strange reason, boost's "negative binomial"  all of the sudden started going haywire when I gave it large values to compute.  The normal approximation seems to work for some reason...
        const mpfr_class expected_lambda((double)it_expected_rate->second + (double)(0.00001L));   
        const mpfr_class observed_value(it_observed->second);		
	
	
// 	if (it_expected_rate->second > 50)
	{//normal approx
	    const mpfr_class mu(expected_lambda);
	    const mpfr_class sig((double)std::max<longdouble>(3,(sqrt(it_expected_rate->second*(1.00 + sqrt(RD_variance_scale___alpha_squared))))));
	    const boost::math::normal_distribution<mpfr_class>  my_normal(mu,sig);    
	    
	    const real result_normal(boost::math::pdf<mpfr_class>(my_normal, observed_value).get_mpfr_t());
	    
	    if (result_normal <= 0.00L)
	    {	    
		std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "   result_normal = " << result_normal.toString(10) 
			<< "  <= 0.00, color = " << it_observed->first << ",   observed = " << it_observed->second << ",   expected = " <<  it_expected_rate->second 
			<< "\nexpected_lambda = " << expected_lambda
			<< "\nobserved_value = " << observed_value
			<< "\nmu = " << mu
			<< "\nsig = " << sig
			<< "\n";
		
		RD_probability  *=  minimum_representable_positive_real_number;
	    }
	    else                
	    {
		RD_probability  *=  result_normal;        
	    }	    
	    
	}//normal approx
// 	else
// 	{//negative binomial
// 	                    
// 	    // CAREFULLY DERIVED VERSION: 	
// 	    const mpfr_class parameter_BOOST_r( the_number_one / (double)RD_variance_scale___alpha_squared);
// 		    //BOOST_r  <==>  my_k  <==> Gamma shape parameter k.                                                 
// 	    
// 	    const mpfr_class parameter_BOOST_p(the_number_one / (the_number_one  +  (expected_lambda*(double)RD_variance_scale___alpha_squared)));
// 	    
// 	    
// 	    const boost::math::negative_binomial_distribution<mpfr_class> my_negative_binomial(parameter_BOOST_r, parameter_BOOST_p);                
// 	    const real result_negbin(boost::math::pdf<mpfr_class>(my_negative_binomial, observed_value).get_mpfr_t());
//         
// 	    if (result_negbin <= 0.00L)
// 	    {	    
// 		const mpfr_class mu((double)it_expected_rate->second);
// 		const mpfr_class sig((double)sqrt(it_expected_rate->second*(1.00 + sqrt(RD_variance_scale___alpha_squared))));
// 		const boost::math::normal_distribution<mpfr_class>  my_normal(mu,sig);
// 		
// 		
// 		std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "   result_pdf = " << result_negbin.toString(10) 
// 			<< "  <= 0.00, color = " << it_observed->first << ",   observed = " << it_observed->second << ",   expected = " <<  it_expected_rate->second 
// 			<< "\nparameter_BOOST_r = " << parameter_BOOST_r
// 			<< "\nparameter_BOOST_p = " << parameter_BOOST_p
// 			<< "\nexpected_lambda = " << expected_lambda
// 			<< "\nobserved_value = " << observed_value
// 			<< "\nboost::math::mean<mpfr_class>(my_negative_binomial) = " << boost::math::mean<mpfr_class>(my_negative_binomial)
// 			<< "\nboost::math::variance<mpfr_class>(my_negative_binomial) = " << boost::math::variance<mpfr_class>(my_negative_binomial)
// 			<< "\nboost::math::pdf<mpfr_class>(my_negative_binomial, observed_value) = " << boost::math::pdf<mpfr_class>(my_negative_binomial, observed_value)
// 			<< "\nboost::math::pdf<mpfr_class>(my_normal, observed_value) = " << boost::math::pdf<mpfr_class>(my_normal, observed_value)
// 			<< "\n";
// 		
// 		RD_probability  *=  minimum_representable_positive_real_number;
// 	    }
// 	    else                
// 	    {
// 		RD_probability  *=  result_negbin;        
// 	    }
// 	}//negative binomial
        
        
        ++it_expected_rate;                
    }//it_observed
    
            
            
            
            


    return RD_probability;
 
}   //  end of    calculate_probability_of_observed_read_depth_given_expected_read_depth_over_natural_poisson_intervals___gender_INsensitive

























































/////////////////////////////////////////////////////////////////////////////////////////                                 
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////                                 
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//PER / MiniRegion stuff:


// We may find only "half-mapped" PERs -  only one of the two mates is mapped.
// So when we encounter a half-PER within "half_frag_length" of a given position, we must be aware that the other mate may map "behind" the already-mapped-mate.
// Thus a MiniRegion should be padded by     "average_fragment_length" + some padding      on either side of the homologous positon.
// But we should only search within one half-fragment length of the homologous positon - the idea is that the PER should be able to determine is this position is a breakpoint or not (hence half on either side).


bool Event::set_up_all_PERs
		(const bool &consider_GeneConversion_breakpoints,
		BamTools::BamReader &my_BAM_reader)
{  
    bool successful_operations = true;    
 
    if (my_MPI_rank == 0)
    {  std::cerr << "\n\t\t\t\t\tremaining_coordinates_for_consideration.size() = " << remaining_coordinates_for_consideration.size() << "\n\n";  }

    
    double time_begin = omp_get_wtime();
       

    
    
    if (my_MPI_rank == 0)  
    {  std::cerr<< "\n\n\t\t create_search_intervals_from_DAR_and_homologous_regs_for_uploading_PERs...\n";  }
    
    {//search and upload
        const type_map_uint_to_list_BI search_intervals(create_search_intervals_from_DAR_and_homologous_regs_for_uploading_PERs
                                                                (maximum_average_fragment_length, search_block_size)  );   
                
        if (my_MPI_rank == 0)  
	{  std::cerr<< "\n\n\t\t upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals - MPI and OMP!!! ...\n";  }                                                        
        
        successful_operations = upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals_____MPI_and_OMP
                                                    (search_intervals, 20, my_BAM_reader);
                                                    
        if (!successful_operations)                                            
	{  return false;  }
                                    
    }//search and upload   
                                
                                 
    
          
    
    
    
    if (my_MPI_rank == 0)  
    {  std::cerr<< "\n\n\t\t for_all_PERs____set_region_sequences___and___orient...\n";  }
    
    
    pre_upload_homologous_regions_from_the_Reference_genome();
    
    for_all_PERs____set_region_sequences___and___orient___MPI_and_OMP();
            
    uploaded_pieces_of_chromos.clear();
  
    
    
    
    
    
    
    if (my_MPI_rank == 0)  
    {  std::fprintf(stderr, "\n\n\nPERs on this profile.size() = %u\n", PERs_on_this_profile.size() );  }
    


    set_up_and_align_PERs_using_GPU____OMP
			(PERs_on_this_profile,
			 this,
			 true,//create_and_consider_MiniEvents,
			 consider_GeneConversion_breakpoints,  
			 true,//sum_over_all_unaffected_LCRs,  //"false" for heuristic
			 false,//special_save_breakpoints_actually_captured_in_alignment); //"true" for heuristic
			    false);//do_PER_orient_formatting
                   
    
            
    
    
   
    
//     world_MPI_boost_ptr->barrier();  //NECESSARY!
    
    
    
//     const std::vector<type_map_string_to_PER::iterator> PERs___loop(get_loopable_iterators<type_map_string_to_PER>(PERs_on_this_profile));
//     
//     //MPI
//     for (int rank=0; rank < size_MPI_WORLD; ++rank) 
//     {
//         if (my_MPI_rank == rank)
//         {
//             type_map_string_to_PER sending_PERs;
//             type_map_string_to_PER::iterator insert_it = sending_PERs.begin();
//             for (int p = my_MPI_rank; p < (int)PERs___loop.size(); p += size_MPI_WORLD)
//             {            
//                 insert_it = sending_PERs.insert(insert_it,  *PERs___loop.at(p)  );
//                 insert_it->second.deep_copy_from_serialize( PERs___loop.at(p)->second, my_Conn_Comp  );
//             }
//             
//             boost::mpi::broadcast<type_map_string_to_PER>( *world_MPI_boost_ptr, sending_PERs, rank );                                            
//         } // sending rank
//         else
//         {
//             type_map_string_to_PER received_PERs;
//             boost::mpi::broadcast<type_map_string_to_PER>( *world_MPI_boost_ptr, received_PERs, rank );
//             
//             type_map_string_to_PER::iterator it_recv = received_PERs.begin();
//             while (  it_recv   !=   received_PERs.end()  )
//             {
//                 const type_map_string_to_PER::iterator found_it = PERs_on_this_profile.find( it_recv->first );                
//                 found_it->second.deep_copy_from_serialize(it_recv->second, my_Conn_Comp);
//                 
//                 received_PERs.erase(it_recv++);   //post-increment necessary.
//             }            
//         
//         } // receiving rank
//     }
    
//     world_MPI_boost_ptr->barrier();  //NECESSARY!
    
    
    std::cerr<< "\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "  done broadcasting!\n";
        
    std::cerr<< "\n\nrank " << my_MPI_rank <<  "   time to set_up_all_PERs    =     "  << omp_get_wtime() - time_begin << "\n\n";  
  

    return successful_operations;
    
} // end of   set_up_all_PERs



























type_map_uint_to_set_uint Event::determine_homologous_variational_positions_on_immediate_neighbors_and_self()  const
{
        
    type_map_uint_to_set_uint varpos_by_neighbor_and_self;    
    const uint number_neighbors = local_interlocking_Events.size() + global_interlocking_Events.size();
    
    
    std::vector<type_map_uint_to_Event::const_iterator>   neighbors;
    neighbors.reserve(number_neighbors);
                {//scope
                    for (type_map_uint_to_Interlock::const_iterator it_neighb = local_interlocking_Events.begin();
                            it_neighb != local_interlocking_Events.end();
                            ++it_neighb)
		    {  neighbors.push_back(   my_Conn_Comp->events.find( it_neighb->second.entry->UID )   );   }
                    
                    for (type_map_uint_to_Interlock::const_iterator it_neighb = global_interlocking_Events.begin();
                            it_neighb != global_interlocking_Events.end();
                            ++it_neighb)
		    {  neighbors.push_back(  my_Conn_Comp->events.find( it_neighb->second.entry->UID )    );   }   
                    
                    //init sets for neighbors.
                    type_map_uint_to_set_uint::iterator insert_init_it = varpos_by_neighbor_and_self.begin();
                    for (uint n=0; n<number_neighbors; ++n)
		    {  insert_init_it = varpos_by_neighbor_and_self.insert(insert_init_it, type_uint__set_uint(neighbors.at(n)->second.UID, type_set_uint()));  }
                }//scope
    
    
    //neighbors
    for (type_list_Coords_and_homologous_region::const_iterator it_homol = regions_homologous_to_directly_affected_region.begin();
            it_homol != regions_homologous_to_directly_affected_region.end();
            ++it_homol)
    {
        const type_set_uint homologous_coordinates(
                                extract_keys_of_compressed_map_and_return_as_set(it_homol->compressed_map_homologous_region_to_profile)  );
                                                                       
        const BOOST_Interval homologous_coordinates_endpoints(*homologous_coordinates.begin(), *(--homologous_coordinates.end()));                       
                                      
        for (uint n=0; n< number_neighbors; ++n)
	{
            for (uint qs=0; qs<2; ++qs)
	    {
                if (neighbors.at(n)->second.chromos[qs] == it_homol->chromosome_of_homologous_region
                     and BOOST_overlap(homologous_coordinates_endpoints, neighbors.at(n)->second.LCRs[qs]))
                {
                    const type_set_uint intersected_coords(get_intersection_of_set_and_interval(homologous_coordinates, neighbors.at(n)->second.LCRs[qs]));
                    if (!intersected_coords.empty())
                    {                    
                        const type_map_uint_to_set_uint::iterator  it_set_of_intersected_varpos = varpos_by_neighbor_and_self.find(neighbors.at(n)->second.UID);
                        
                        const type_map_uint_to_uint map_LCR_to_varpos( convert_compressed_map_to_full_map(
                                                    neighbors.at(n)->second.compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs))   );
                        for (type_set_uint::const_iterator it_intersect = intersected_coords.begin();
                                it_intersect != intersected_coords.end();
                                ++it_intersect)
                        {
                            const type_map_uint_to_uint::const_iterator it_find = map_LCR_to_varpos.find(*it_intersect);
                            if ( it_find != map_LCR_to_varpos.end() )
			    {  it_set_of_intersected_varpos->second.insert(it_find->second);  }                                  
                        }    
                    }//overlap varpos
                }//overlap
	    }//qs
	}//neighbors
                
    }//homol regs       

    
    //self
    varpos_by_neighbor_and_self.insert(type_uint__set_uint(UID, variational_positions_of_profile));
    
    
    
    return  varpos_by_neighbor_and_self; 

}  // determine_homologous_variational_positions_on_immediate_neighbors


















































real Event::calculate_or_get_Read_Depth_for_these_state_vectors_and_breakpoints
                                        (const Sparse_map &states__haploid_0,
                                         const Sparse_map &states__haploid_1,
                                         const type_map_uint_to_BI &breakpoints__haploid_0,
                                         const type_map_uint_to_BI &breakpoints__haploid_1,
					 std::stringstream *const &display_strm__ptr)   const
{           
    
    const type_map_uint_to_uint_2 diploid_states(states__haploid_0.sparse_map_UID_to_value, states__haploid_1.sparse_map_UID_to_value);
    
    
    type_map_uint_to_longdouble expected_diploid_fragmentation_rate_per_poisson_interval(
					    base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive);
					    
    adjust_expected_diploid_fragmentation_rate___breakpoint_specific
                                            (expected_diploid_fragmentation_rate_per_poisson_interval,
                                            states__haploid_0.sparse_map_UID_to_value,  breakpoints__haploid_0);                                               
                                                                                
    adjust_expected_diploid_fragmentation_rate___breakpoint_specific
                                            (expected_diploid_fragmentation_rate_per_poisson_interval,
                                            states__haploid_1.sparse_map_UID_to_value,  breakpoints__haploid_1); 
                                                          
//     //add noise.  
//     //When a Natural Poisson Interval whose expected read-depth lambda is "0" under these diploid event states, the corresponding NegativeBinomial (or Poisson) distribution used to evaluate the likelihood of the observed data given the expected value lambda is undefined.  Obviously, an easy fix would be to add something like "0.000001" to every expedted value lambda so that the distributions are always defined.   But when lambda is supposd to be 0, and the observed count is say, 10, and under even a haploid deletion, the expecdted lambda is say 200, then setting lambda = "0.00001" will make observing "10" extremely unlikely, even though 10 is obviously much closer to 0 (diploid deletion) than to 200 (haploid deletion).  This means that we never call diploid deletions, even when obvious.  We interpret these 10 PERs to be "noise" (possible junk reads, or by chance there is a near-exact tiny repeat that matches that read between the given region and somewhere else (outside of my SD database) in the genome).  To control for the outlying/edge effects, we add the noise term and make it proportional to the base fragmentation reate. 
//     {//noise
// 	type_map_uint_to_longdouble::const_iterator it_base = base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive.begin();
// 	for (type_map_uint_to_longdouble::iterator it_rd = expected_diploid_fragmentation_rate_per_poisson_interval.begin();
// 		it_rd != expected_diploid_fragmentation_rate_per_poisson_interval.end();
// 		++it_rd)    
// 	{    
// 	    if (it_rd->second < 5)
// 	    {
// 		if (it_rd->second < 0)
// 		{  it_rd->second = base_noise_read_depth_percentage_of_normal_Reference_read_depth*(it_base->second/2);  }
// 		else
// 		{  it_rd->second += base_noise_read_depth_percentage_of_normal_Reference_read_depth*(it_base->second/2); }
// 	    }
// 	    ++it_base;
// 	}//rd
//     }//noise
                                         
                                            
                                                     
    const real P_Read_depth(calculate_probability_of_observed_read_depth_given_expected_read_depth_over_natural_poisson_intervals(
                                      expected_diploid_fragmentation_rate_per_poisson_interval));  
				      
    if (display_strm__ptr != NULL)
    {
	print_map_keys_and_values<uint, longdouble>(expected_diploid_fragmentation_rate_per_poisson_interval,
						    "example_expected_diploid_fragmentation_rate_per_poisson_interval",
						    display_strm__ptr,
						    true);                                                   
	
	(*display_strm__ptr) << "\n\t\t\t\t\t\t\t\t\t" << P_Read_depth << "\n\n";            	    
    }//display
                                      

    return P_Read_depth;                                                          

} // calculate_or_get_Read_Depth_for_these_state_vectors_and_breakpoints














































void Event::pre_upload_homologous_regions_from_the_Reference_genome
			(const int &padding_amount_indicator)
{
    const uint padding_amount =  (padding_amount_indicator < 0)   ?   
					(maximum_average_fragment_length*30   + (uint)maximum_frag_length_absolute_dev*15  + 20000)
					:
					(uint)padding_amount_indicator;
    
    
    type_map_uint_to_list_BI map_chromo_to_homol_regs;  
    
    //region of LCRs
    map_chromo_to_homol_regs[chromos[0]].push_back(region_between_and_including_the_LCRs_themselves);        
//     map_chromo_to_homol_regs[ chromos[0] ]
//                         .push_back( 
//                             BOOST_Interval( 
//                                     (uint)std::max<int>(0,
//                                                         ((int)region_between_and_including_the_LCRs_themselves.lower())
//                                                             - (int)maximum_average_fragment_length*5  - maximum_frag_length_absolute_dev*10 
//                                                             - (int)padding_for_region_homologous_to_a_var_pos),
//                                     std::min<uint>( chromosome_lengths.at(chromos[0]),
//                                                     region_between_and_including_the_LCRs_themselves.upper() 
//                                                             + maximum_average_fragment_length*5 + (uint)maximum_frag_length_absolute_dev*10 
//                                                             + padding_for_region_homologous_to_a_var_pos )      ));
        

    //homologous regions                                
    for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg = regions_homologous_to_directly_affected_region.begin();
            it_homol_reg != regions_homologous_to_directly_affected_region.end();
            ++it_homol_reg)                
    {
        map_chromo_to_homol_regs[it_homol_reg->chromosome_of_homologous_region]
                        .push_back( get_endpoints_of_keys_of_compressed_map(it_homol_reg->compressed_map_homologous_region_to_profile)  );             
    }
    
                
                        
                        
    //expand ALL regions to allow for padding:
    for (type_map_uint_to_list_BI::iterator it_chr = map_chromo_to_homol_regs.begin();
            it_chr != map_chromo_to_homol_regs.end();
            ++it_chr)
        for (type_list_BI::iterator it_BI = it_chr->second.begin();
                it_BI != it_chr->second.end();
                ++it_BI)
        {                   
            const uint LB =  safe_subtract_base_1(it_BI->lower(), padding_amount);            	    
            const uint UB =  safe_chromo_add_base_1(it_BI->upper(), padding_amount, it_chr->first);  
	    
            it_BI->set(LB,UB);
        }    


    merge_intervals_on_each_chromosome(map_chromo_to_homol_regs);                             
    
    uploaded_pieces_of_chromos = upload_list_of_chromosome_regions_from_file(map_chromo_to_homol_regs);   

} // pre_upload_homologous_regions_from_the_Reference_genome
















































// typedef std::map<BOOST_Interval, type_list_Balgmts>   type_map_BI_to_Balgmts;
// typedef  std::map<uint, type_map_BI_to_Balgmts>  type_map_chr_to_BI_to_Balgmts;
// 
// 
// type_map_chr_to_BI_to_Balgmts Event::make_intervals_from_subsequence_of_coordinates
//                         (const type_set_uint &genome_coords,
//                          const uint &chromo,
//                          const uint &padding_amount_below)   // half_minimum_average_fragment_length
// {
//     
//     type_map_chr_to_BI_to_Balgmts  map_regions_to_Balgmts;
//     type_map_BI_to_Balgmts::const_iterator it_insert_BI = map_regions_to_Balgmts.insert(map_regions_to_Balgmts.begin(),
//                                                                                                std::pair<uint, type_map_BI_to_Balgmts>
//                                                                                                       (chromo, type_map_BI_to_Balgmts() ) );
//     type_map_BI_to_Balgmts *const regions_on_this_chromo  =  &map_regions_to_Balgmts.at(chromo);
//     
//     uint current_lower_bound = std::max<int>(0, (int)(*genome_coords.begin()) - (int)padding_amount_below);        
//     
//     
//     for (type_set_uint::const_iterator  it_coord  =  ++genome_coords.begin(),
//                                         it_coord_prev = genome_coords.begin();
//             it_coord  !=  genome_coords.end();
//             ++it_coord, ++it_coord_prev)
//         if ( *it_coord_prev + padding_amount_below < *it_coord)  // by adding, we avoid the danger of going below 0  (i.e. if we did   *it_coord  - padding_amount_below).
//         {
//             it_insert_BI = regions_on_this_chromo->insert( it_insert_BI, 
//                                                            std::pair<BOOST_Interval, type_list_Balgmts>
//                                                                 ( BOOST_Interval(current_lower_bound, *it_coord_prev),  type_list_Balgmts()    )   );
//                                                                 
//             current_lower_bound =  *it_coord - padding_amount_below;   //safe to subtract in this case, because we know it is > 0.                                                      
//         }    
// 
// 
//     return  map_regions_to_Balgmts;
//     
// }                        
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
// void Event::preload_all_PERs_in_given_region
//                 ( )
// {
// 
//     type_map_chr_to_BI_to_Balgmts  map_regions_to_Balgmts;
//     
//     type_map_chr_to_BI_to_Balgmts new_map_regions_to_Balgmts(  make_intervals_from_subsequence_of_coordinates(remaining_coordinates_for_consideration,
//                                                                                                                   chromos[0],
//                                                                                                                   half_minimum_average_fragment_length)  );
//     
//     for (type_map_chr_to_BI_to_Balgmts::const_iterator it_new = new_map_regions_to_Balgmts.begin();
//             it_new != new_map_regions_to_Balgmts.end();
//             ++it_new)
//     {
//         type_map_chr_to_BI_to_Balgmts::iterator it_find_chromo = map_regions_to_Balgmts.find(it_new->first);
//         
//         if (it_find_chromo == map_regions_to_Balgmts.end() )
//             map_regions_to_Balgmts.insert( *it_new );
//         else
//         {
//             it_find_chromo->second.
//             
//         
//         
//         }
//             
//     
//     }
//    
//     
//     for (type_set_uint::const_iterator it_pos = remaining_coordinates_for_consideration.begin(),
//                                        it_next_pos = ++remaining_coordinates_for_consideration.begin();
//             it_pos != remaining_coordinates_for_consideration.end();   // CAREFUL!!!!
//             ++it_pos, ++it_next_pos)           
//     
//     ;
// 
// 
// 
// 
// 
// 
// }
























type_map_uint_to_list_BI Event::create_search_intervals_from_DAR_and_homologous_regs_for_uploading_PERs
                                                                (const uint &search_interval_padding_amount,
                                                                 const uint &subdivison_size)
{

    type_map_uint_to_list_BI map_chromo_to_homol_regs;  
    
    
    if ( !remaining_coordinates_for_consideration.empty() )
    {//remaining_DAR
        int last_endpoint__homol_reg =  *remaining_coordinates_for_consideration.begin();                  
        type_map_uint_to_list_BI::iterator it_chr = map_chromo_to_homol_regs.insert(std::pair<uint, type_list_BI>(chromos[0], type_list_BI())).first;
                    
        for (type_set_uint::const_iterator it_rem = remaining_coordinates_for_consideration.begin(),
                                           it_rem_next = ++remaining_coordinates_for_consideration.begin();
                it_rem_next != remaining_coordinates_for_consideration.end();
                ++it_rem, ++it_rem_next)
        {     
            if (  *it_rem_next - *it_rem  >  search_interval_padding_amount )
            {
                it_chr->second.push_back(BOOST_Interval(last_endpoint__homol_reg, *it_rem));                
                last_endpoint__homol_reg = *it_rem_next;
            }
        }   
        
        //last
        it_chr->second.push_back(BOOST_Interval((uint)last_endpoint__homol_reg, *(--remaining_coordinates_for_consideration.end()) ));
    }//remaining DAR
    
    
    

    //homologous regions                                
    for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg = regions_homologous_to_directly_affected_region.begin();
            it_homol_reg != regions_homologous_to_directly_affected_region.end();
            ++it_homol_reg)  
    {
        
        if (it_homol_reg->compressed_map_homologous_region_to_profile.empty()  or  it_homol_reg->compressed_map_homologous_region_to_profile.at(0).empty())
            continue;        
        
        type_map_uint_to_list_BI::iterator it_chr = map_chromo_to_homol_regs.insert( std::pair<uint, type_list_BI>
                                                                                     (it_homol_reg->chromosome_of_homologous_region, type_list_BI())  ).first;
                
                                                                                     
        int last_endpoint__homol_reg = it_homol_reg->compressed_map_homologous_region_to_profile.at(0).at(0);        
        int next_putative_endpoint__homol_reg = last_endpoint__homol_reg;
        
        for (type_vector_int::const_iterator it_hr_pos = it_homol_reg->compressed_map_homologous_region_to_profile.at(0).begin(),
                                             it_DAR_pos = it_homol_reg->compressed_map_homologous_region_to_profile.at(1).begin();
                it_hr_pos != it_homol_reg->compressed_map_homologous_region_to_profile.at(0).end();
                ++it_hr_pos, ++it_DAR_pos)
            if ( *it_hr_pos < 0)                
	    {  next_putative_endpoint__homol_reg  +=  *it_DAR_pos;  }
            else
            {
                if ( *it_hr_pos - next_putative_endpoint__homol_reg  >  (int)search_interval_padding_amount )
                {
                    it_chr->second.push_back(  BOOST_Interval( (uint)last_endpoint__homol_reg,  (uint)next_putative_endpoint__homol_reg)   );

                    last_endpoint__homol_reg = *it_hr_pos;
                    next_putative_endpoint__homol_reg = last_endpoint__homol_reg;                                      
                }
                else
		{  next_putative_endpoint__homol_reg = *it_hr_pos;  }
                
//                 if ( next_putative_endpoint__homol_reg - *it_hr_pos  >  (int)search_interval_padding_amount )
//                 {
//                     it_chr->second.push_back(  BOOST_Interval(   (uint)std::min<int>(*it_hr_pos, last_endpoint__homol_reg),
//                                                                  (uint)std::max<int>(*it_hr_pos, last_endpoint__homol_reg)   ));
//                                                                  
//                     last_endpoint__homol_reg = next_putative_endpoint__homol_reg;                                            
//                     
// //                     it_chr->second.push_back(  BOOST_Interval( (uint)last_endpoint__homol_reg,  (uint)next_putative_endpoint__homol_reg)   );
//                     
// //                     last_endpoint__homol_reg = *it_hr_pos;
// //                     next_putative_endpoint__homol_reg = last_endpoint__homol_reg;                         
//                 }
//                 else
//                     next_putative_endpoint__homol_reg = *it_hr_pos;
              
            }
            
        //last
        it_chr->second.push_back(BOOST_Interval((uint)last_endpoint__homol_reg, (uint)next_putative_endpoint__homol_reg));
    }
      
      
      
      
      
    //expand ALL regions to allow for padding:
    map_chromo_to_homol_regs = pad_all_regions(map_chromo_to_homol_regs, search_interval_padding_amount);    

    merge_intervals_on_each_chromosome(map_chromo_to_homol_regs);    
    
//     subdivide_intervals(map_chromo_to_homol_regs, subdivison_size);
                               
    
    return map_chromo_to_homol_regs;    
    
}//create_search_intervals_from_DAR_and_homologous_regs_for_uploading_PERs
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
bool Event::upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals_____MPI_and_OMP
                                (const type_map_uint_to_list_BI &map_chromo_to_homol_regs,
                                 const uint &MR_padding_amount_for_proper_pair,
				BamTools::BamReader &my_BAM_reader)
{
    
    //we now have the set of regions for which we should upload PERs from!!!
//     PERs_on_this_profile.clear();
        
    uint number_of_searching_regions = 0;
                                    for (type_map_uint_to_list_BI::const_iterator it_init = map_chromo_to_homol_regs.begin();
                                            it_init != map_chromo_to_homol_regs.end();
                                            ++it_init)                                           
				    {  number_of_searching_regions += it_init->second.size();  }
                        
    std::vector<type_uint__BI>  searching_regions;
    searching_regions.reserve(number_of_searching_regions);
                                    {//init
                                        for (type_map_uint_to_list_BI::const_iterator it_init = map_chromo_to_homol_regs.begin();
                                                it_init != map_chromo_to_homol_regs.end();
                                                ++it_init)   
					{
                                            for (type_list_BI::const_iterator it_bi = it_init->second.begin();
                                                    it_bi != it_init->second.end();
                                                    ++it_bi)
					    {  searching_regions.push_back( type_uint__BI(it_init->first, *it_bi)  );  }
					}
                                    }//init        
                
    
    
    
    
    
    const type_map_uint_to_uint map_BAM_Ref_IDs_to_chromosome_value(get_inverse_of_map(map_chromosome_value_to_BAM_Ref_IDs));

    
    {
        uint amount_of_search = 0;
        for (uint ctr=0; ctr < number_of_searching_regions; ++ctr)
	{  amount_of_search += BOOST_width_inclusive( searching_regions.at(ctr).second );  }
        
        std::stringstream diag_ss;
        diag_ss << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num()
                << "   in \"upload_BAM_alginments...\"  has    " << amount_of_search << "  bases of genome to upload from...\n\n";
        std::fprintf(stderr, "\n\n%s\n\n", diag_ss.str().c_str() );       
    }
    
    
    
    
    
    
    
    
    const double time_upload_raw_PERs_begin = omp_get_wtime();
    
    
    std::vector<type_list_PER> found_PERs___by_search_region(number_of_searching_regions, type_list_PER());
    
    //raw upload all PERs, pair them, and give them spawning MiniRegions.
    
    bool success_raw_upload = true;
    
    #pragma omp parallel
    {        
        //MPI and OMP !!!!!!
        #pragma omp for schedule(dynamic,1)        
        for (int region_ctr = my_MPI_rank;  region_ctr < number_of_searching_regions; region_ctr += size_MPI_WORLD)
        {
            type_list_Balgmts found_Balgmts;
            
            #pragma omp critical (upload_one_at_a_time_because_of_read_head)
            {//one upload at a time
                
                const bool success_set_region 
                    = my_BAM_reader.SetRegion(create_BAM_region(searching_regions.at(region_ctr).first, searching_regions.at(region_ctr).second, 0));
                    
                if (!success_set_region)
                {
                    std::stringstream error_strm;
                    print_line_of_markers("ERROR! ", &error_strm);
                    print_line_of_markers("(", &error_strm);
                    error_strm << "success_set_region  failed in  \"gather_all_PERs_that_appropriately_intersect_this_MiniRegion_and__format_them_so_that_they_should_align_directly_to_the_spawning_Event_profile\"\n\n";
                    char tempout[option_size];
                    std::sprintf(tempout, ".SetRegion( map_chromosome_value_to_BAM_Ref_IDs.at(%u) = %u, search_portion_of_region_interval.lower() - 1 = %u\n, map_chromosome_value_to_BAM_Ref_IDs.at(%u) = %u, search_portion_of_region_interval.upper() - 1 = %u",
                                searching_regions.at(region_ctr).first,
                                map_chromosome_value_to_BAM_Ref_IDs.at(searching_regions.at(region_ctr).first),
                                searching_regions.at(region_ctr).second.lower() - 1,
                                searching_regions.at(region_ctr).first,
                                map_chromosome_value_to_BAM_Ref_IDs.at(searching_regions.at(region_ctr).first),
                                searching_regions.at(region_ctr).second.upper() - 1  );
                    error_strm << tempout;       
                                
                    print_line_of_markers(")", &error_strm);
                    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );  
                    
                    #pragma omp critical (raw_upload_fail_set)
                    success_raw_upload = false;
                }                           
                else
                {                
                    BamTools::BamAlignment Balgmt;
                    
                    while (my_BAM_reader.GetNextAlignmentCore(Balgmt))  
                    {
                        if (!Balgmt.IsFailedQC())  
                        {
                            if (Balgmt.IsProperPair() or   ( Balgmt.IsPaired()  and  Balgmt.IsMapped()  and   Balgmt.IsMateMapped()))  
                            {
                                Balgmt.BuildCharData();
                                found_Balgmts.push_back(Balgmt);                    
                            }
                        }
                    }//while
                }//success            
        
            //Now pair up these reads to form Paired_end_reads;
	    const type_uint__uint failed_PERs_and_discarded__counts
			= create_list_of_unformatted_Paired_end_reads_from_list_of_Balgmts(found_Balgmts, my_BAM_reader,
											    found_PERs___by_search_region.at(region_ctr));
            found_Balgmts.clear(); //should already be erased anyway!!  (just for emphasis).
 
        }//one upload at a time

            
            
            //by construction below, the first Miniregion is always the spawning one...
            for (type_list_PER::iterator it_per = found_PERs___by_search_region.at(region_ctr).begin();
                    it_per != found_PERs___by_search_region.at(region_ctr).end();
                    ++it_per)
	    {
                if (it_per->mates.first.IsProperPair())
                {
                    //consistency
                    if (it_per->mates.first.Position  >   it_per->mates.second.Position)
		    {  std::swap<BamTools::BamAlignment>(it_per->mates.first, it_per->mates.second);  }                     
                    
                    const uint chr_val_of_first_mate = map_BAM_Ref_IDs_to_chromosome_value.at(  it_per->mates.first.RefID  );
                    
                    const BOOST_Interval spawning_MR_region_interval(
                        (uint)std::max<int>(1, it_per->mates.first.Position+1 - (int)MR_padding_amount_for_proper_pair ),
                        std::min<uint>(chromosome_lengths.at(chr_val_of_first_mate), 
                                       (uint)it_per->mates.second.Position+1 + (uint)it_per->mates.second.Length  + MR_padding_amount_for_proper_pair )  );
                    
                        
                        
                    bool orientations_agree, must_complement;                    
                    determine_relation_of_region_to_Event_profile(
                                                    orientations_agree,
                                                    must_complement,
                                                    chr_val_of_first_mate,
                                                    spawning_MR_region_interval  );
                        
                    const MiniRegion spawning_MR(
                                            chr_val_of_first_mate,
                                            spawning_MR_region_interval,
                                            orientations_agree,
                                            must_complement);
                            
                    it_per->add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately( spawning_MR );
                
                }//ProperPair
                else  // only Paired not Proper
                {
                    
                    const int MR_additional_padding_NON_proper  
                                    =         it_per->my_Readgroup_ptr->median_insert_size____absolute_standard_deviation < 20           ?
                                                    it_per->my_Readgroup_ptr->median_insert_size____absolute_standard_deviation*10
                                                    :
                                                    (int)ceil(1.5*(float)it_per->my_Readgroup_ptr->median_insert_size____absolute_standard_deviation);
                    
                                                                         
                    {//MR around first mate 
                        const type_map_uint_to_uint::const_iterator it_find_BAM_Ref_ID_to_chro_val__mate_first
                                            = map_BAM_Ref_IDs_to_chromosome_value.find(it_per->mates.first.RefID);
                        
                        assert(  it_find_BAM_Ref_ID_to_chro_val__mate_first != map_BAM_Ref_IDs_to_chromosome_value.end() );
                                              
                        const BOOST_Interval spawning_MR_region_interval(
                                (uint)std::max<int>(1, it_per->mates.first.Position+1
                                                            - it_per->my_Readgroup_ptr->median_insert_size - MR_additional_padding_NON_proper  ),
                                std::min<uint>( chromosome_lengths.at( it_find_BAM_Ref_ID_to_chro_val__mate_first->second),
                                                ((uint)it_per->mates.first.Position+1  + (uint)it_per->mates.first.Length )
                                                        + it_per->my_Readgroup_ptr->median_insert_size + MR_additional_padding_NON_proper)  );
                        
                                                        
                        bool orientation_spawning_MR_and_profile_agree, must_complement_spawning_MR_to_get_to_profile;                            
                        determine_relation_of_region_to_Event_profile(
                                                    orientation_spawning_MR_and_profile_agree,
                                                    must_complement_spawning_MR_to_get_to_profile,
                                                    it_find_BAM_Ref_ID_to_chro_val__mate_first->second,
                                                    spawning_MR_region_interval );                                                                                             
                                                        
                        const MiniRegion  spawning_MR(
                                it_find_BAM_Ref_ID_to_chro_val__mate_first->second,
                                spawning_MR_region_interval,
                                orientation_spawning_MR_and_profile_agree,
                                must_complement_spawning_MR_to_get_to_profile); 
                                
                        it_per->add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately( spawning_MR );                    

                    }//MR around first mate  
                                                                        
                }//not proper pair
	    }
                                               
        }//region ctr
     
        mpfr_free_cache();
    }//omp parallel
        
        
        
        
        
        
        
    searching_regions.clear();
        
        
        
    
        
        
        
        
    {
        std::stringstream  diag_ss;
        diag_ss   << "\n\t\tinside \"upload_BAM_alignments...\".  uploaded raw PERs.  time = " << omp_get_wtime() - time_upload_raw_PERs_begin << "  s.\n"
		  <<"\t\tNow adding miniregions...\n\n";
	std::cerr << diag_ss.str();
    }
        
    const double time_add_MRs_begin = omp_get_wtime();
        
        
        
        
        
        
        
        
    // add miniregions

    type_map_string_to_PER  PERs_this_Event__this_rank; 
    
    
    
//     std::cerr << "\nsize_MPI_WORLD = " << size_MPI_WORLD << "\t\tnumber_of_searching_regions = " << number_of_searching_regions << "\n";
    

    for (int region_ctr = my_MPI_rank;  region_ctr < number_of_searching_regions; region_ctr += size_MPI_WORLD)
    {
        
//         std::cerr << "\nregion_ctr = " << region_ctr << "\n";
        
        const uint number_of_PERS__this_search_region = found_PERs___by_search_region.at(region_ctr).size();
        
        std::vector<type_list_PER::iterator>  found_PERs___by_search_region___loop( number_of_PERS__this_search_region, found_PERs___by_search_region.at(region_ctr).end() );
                                            {//init
                                                uint init_ctr_reg = 0;
                                                for (type_list_PER::iterator it_PER_search_reg = found_PERs___by_search_region.at(region_ctr).begin();
                                                        it_PER_search_reg != found_PERs___by_search_region.at(region_ctr).end();
                                                        ++it_PER_search_reg)
						{  found_PERs___by_search_region___loop.at(init_ctr_reg++) = it_PER_search_reg;  }
                                            }//init   
                                            
//         std::cerr << "\nnumber_of_PERS__this_search_region = " << number_of_PERS__this_search_region << "\n";                                  
        
        #pragma omp parallel
        {
            
            #pragma omp for schedule(dynamic,20)
            for (uint per_ctr = 0; per_ctr < number_of_PERS__this_search_region; ++per_ctr)       
            {               
            
                type_list_BI homol_DAR_regs;   
                
                
                const MiniRegion *const first_spawning_MR_ptr = found_PERs___by_search_region___loop.at(per_ctr)->get_ptr_to_first_MiniRegion();  
                assert(  first_spawning_MR_ptr  != NULL  );
                
                if (    first_spawning_MR_ptr->chromosome_of_region == chromos[0]
                        and   BOOST_overlap( first_spawning_MR_ptr->region_interval, region_between_and_including_the_LCRs_themselves )  )
                {
                    homol_DAR_regs.push_back(first_spawning_MR_ptr->region_interval);                                
                }
                
                
                if (BOOST_empty(first_spawning_MR_ptr->region_interval))
                {
                    std::stringstream error_strm;                    
                    error_strm << "\n\nWARNING!    SOMEHOW  the first_spawning MR  has an empty reigon interval!!\n\n";                    
                    first_spawning_MR_ptr->print_this_MiniRegion( &error_strm );
                    found_PERs___by_search_region___loop.at(per_ctr)->print_this_PER( &error_strm );                    
                    warning_message(error_strm, false);
                }
                
//                 const double time_begin_all = omp_get_wtime();
                
//                 double time_piece = time_begin_all;
                
                {//get homol regs on DAR                                                                      
                    for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg = regions_homologous_to_directly_affected_region.begin();
                            it_homol_reg != regions_homologous_to_directly_affected_region.end();
                            ++it_homol_reg)
		    {
                        if ( first_spawning_MR_ptr->chromosome_of_region == it_homol_reg->chromosome_of_homologous_region
                                and BOOST_overlap( first_spawning_MR_ptr->region_interval,
                                                get_endpoints_of_keys_of_compressed_map(it_homol_reg->compressed_map_homologous_region_to_profile)  )    )
                        {
                            
                            const type_map_uint_to_uint full_map_MR_to_DAR(        
                                        convert_compressed_map_to_full_map(
                                                get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                        it_homol_reg->compressed_map_homologous_region_to_profile,
                                                        first_spawning_MR_ptr->region_interval )
                                                ));                            
//                                         get_map_intersection_of_map_keys_and_interval<uint,uint>(
//                                                     convert_compressed_map_to_full_map(it_homol_reg->compressed_map_homologous_region_to_profile),
//                                                     first_spawning_MR_ptr->region_interval)    );
                            
                            if (  ! full_map_MR_to_DAR.empty()  )
                            {                                
                                uint DAR_lower_endpoint;
                                uint DAR_upper_endpoint;                                                          
                                            {
                                                const type_map_uint_to_uint full_map_DAR_to_MR(  get_inverse_of_map( full_map_MR_to_DAR )   );                                
                                                DAR_lower_endpoint = full_map_DAR_to_MR.begin()->first;
                                                DAR_upper_endpoint = (--full_map_DAR_to_MR.end())->first;                                                                                                                                    
                                            }                    
                                
                                
                                //extend
                                const int diff_from_lower = ((int)full_map_MR_to_DAR.begin()->first) - (int)first_spawning_MR_ptr->region_interval.lower();
                                if (  diff_from_lower > 0  )                        
                                {
                                    if (it_homol_reg->orientation_of_homologous_region_and_profile_agree)
				    {  DAR_lower_endpoint = safe_subtract_base_1(DAR_lower_endpoint, diff_from_lower);  }
                                    else
				    {  DAR_upper_endpoint = safe_chromo_add_base_1(DAR_upper_endpoint, diff_from_lower, chromos[0]);  }
                                }
                                
                                const int diff_from_upper = ((int)first_spawning_MR_ptr->region_interval.upper()) - (int)(--full_map_MR_to_DAR.end())->first;
                                if (  diff_from_upper > 0  )
                                {
                                    if (it_homol_reg->orientation_of_homologous_region_and_profile_agree)
				    {
					DAR_upper_endpoint = safe_chromo_add_base_1(DAR_upper_endpoint, diff_from_upper, chromos[0]);
				    }
                                    else                                                                         
				    {  DAR_lower_endpoint = safe_subtract_base_1(DAR_lower_endpoint, diff_from_upper);  }                                                                                                      
                                }
                            
                                assert(  DAR_lower_endpoint  <  DAR_upper_endpoint );
                                homol_DAR_regs.push_back(   BOOST_Interval(   DAR_lower_endpoint,  DAR_upper_endpoint )   ); 
                            }//not empty
                        }
		    }
                }//get homol_regs on DAR
                
                
//                 std::cerr << "\n\tget homol_regs on DAR:  time = " << omp_get_wtime() - time_piece << ".   size = " << homol_DAR_regs.size() << "\n";
//                 time_piece = omp_get_wtime();
                
                //now go through and add homologous MRs.
                for (type_list_BI::const_iterator it_DAR_reg = homol_DAR_regs.begin(); 
                        it_DAR_reg != homol_DAR_regs.end();
                        ++it_DAR_reg)
                {
                    bool near_or_contained_in_LG = false;
//                     bool near_lg_or_end__below = false;
//                     bool near_lg__or_end__above = false;

                    
                    for (type_set_BI::const_iterator it_lg = large_gaps_of_profile.begin();
                            it_lg != large_gaps_of_profile.end();
                            ++it_lg)
                    {                              
                        const BOOST_Interval padded_lg_profile( (uint)std::max<int>(0, ((int)it_lg->lower()) - 10),
                                                                it_lg->upper() + 10 );
                        
                        for (uint qs = 0; qs < 2 ; ++qs)
                        {
                            
                            const type_map_uint_to_uint map_LCR_to_profile__lg(
                                    convert_compressed_map_to_full_map(
                                            get_compressed_map_inverse_of_compressed_map(
                                                get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                        get_compressed_map_inverse_of_compressed_map(  compressed_map_LCR_to_profile__for_each_LCR.at(qs)  ),
                                                        padded_lg_profile  )
                                                )));                                                        
//                                 get_inverse_of_map(  // LCR to profile
//                                     get_map_intersection_of_map_keys_and_interval<uint,uint>(
//                                         get_inverse_of_map(                                                             // profile to LCR
//                                                 convert_compressed_map_to_full_map(  compressed_map_LCR_to_profile__for_each_LCR.at(qs)  ) ),  //LCR to profile
//                                         padded_lg_profile  )  ));
                                        
                            if (  ! map_LCR_to_profile__lg.empty()  )
                            {
                                const BOOST_Interval padded_lg_on_LCR(  map_LCR_to_profile__lg.begin()->first, (--(map_LCR_to_profile__lg.end()))->first  );                                    
                                
                                if (  BOOST_overlap(padded_lg_on_LCR, *it_DAR_reg)  )
                                {
                                    near_or_contained_in_LG = true;
//                                     if (  abs( ((int)padded_lg_on_LCR.lower())  -  (int)it_DAR_reg->lower() ) < 20  )
//                                         near_lg_or_end__below = true;
//                                     if (  abs( ((int)padded_lg_on_LCR.upper())  -  (int)it_DAR_reg->upper() ) < 20  )
//                                         near_lg__or_end__above = true;                        
                                }                                                                                                
                            }
                        } //qs
                                                                                                                                                                                
                    }// lg                                        
                    


                    for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg = regions_homologous_to_directly_affected_region.begin();
                            it_homol_reg != regions_homologous_to_directly_affected_region.end();
                            ++it_homol_reg)
		    {
                        if ( BOOST_overlap(*it_DAR_reg, get_endpoints_of_values_of_compressed_map(it_homol_reg->compressed_map_homologous_region_to_profile)))
                        {                                                  
                            const type_map_uint_to_uint map_DARreg_to_some_homol_reg(
                                        convert_compressed_map_to_full_map(
                                            get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
                                                get_compressed_map_inverse_of_compressed_map(  it_homol_reg->compressed_map_homologous_region_to_profile  ),
                                                *it_DAR_reg   )
                                            )  );
					    
                            
                            if (!map_DARreg_to_some_homol_reg.empty())
                            {
                                
                                uint homolreg_lower_endpoint
                                            = std::min<uint>(map_DARreg_to_some_homol_reg.begin()->second, (--map_DARreg_to_some_homol_reg.end())->second);
                                uint homolreg_upper_endpoint
                                            = std::max<uint>(map_DARreg_to_some_homol_reg.begin()->second, (--map_DARreg_to_some_homol_reg.end())->second);                                
                                
                                const int diff_from_lower = ((int)map_DARreg_to_some_homol_reg.begin()->first) - (int)it_DAR_reg->lower();
                                if (  diff_from_lower > 0  )                        
                                {                                                                                
                                    if (it_homol_reg->orientation_of_homologous_region_and_profile_agree)
				    {  homolreg_lower_endpoint = safe_subtract_base_1(homolreg_lower_endpoint, diff_from_lower);  }
                                    else
				    {
					homolreg_upper_endpoint = safe_chromo_add_base_1(homolreg_upper_endpoint, diff_from_lower,
											it_homol_reg->chromosome_of_homologous_region);
				    }
                                }
                                
                                const int diff_from_upper = ((int)it_DAR_reg->upper())  -  (int)(--map_DARreg_to_some_homol_reg.end())->first;
                                if (  diff_from_upper > 0  )
                                {                                    
                                    if (it_homol_reg->orientation_of_homologous_region_and_profile_agree)
				    {
                                        homolreg_upper_endpoint 
                                                = safe_chromo_add_base_1(homolreg_upper_endpoint, 
									 diff_from_upper,
									 it_homol_reg->chromosome_of_homologous_region);
				    }
                                    else                                                                         
				    {  homolreg_lower_endpoint = safe_subtract_base_1(homolreg_lower_endpoint, diff_from_upper);  }				                                                                                                     
                                }      
                                                                                                
                                
//                                 if (near_lg_or_end__below  or  near_lg__or_end__above)
                                if (near_or_contained_in_LG)
                                {        
				    homolreg_lower_endpoint
					= safe_subtract_base_1(homolreg_lower_endpoint, 
							    found_PERs___by_search_region___loop.at(per_ctr)->my_Readgroup_ptr->median_insert_size);
						
				    homolreg_upper_endpoint
					= safe_chromo_add_base_1(homolreg_upper_endpoint,
							   found_PERs___by_search_region___loop.at(per_ctr)->my_Readgroup_ptr->median_insert_size,
							   it_homol_reg->chromosome_of_homologous_region);				                                                                    
                                }
                                
                                
                                bool orientations_agree = it_homol_reg->orientation_of_homologous_region_and_profile_agree;
                                bool must_complement = it_homol_reg->homologous_region_must_be_complemented_to_get_to_profile;
                                
                                if (recomb_type == recomb_class__Inv  and  BOOST_overlap(*it_DAR_reg, LCRs[1])  )   
                                {                                       
                                    orientations_agree = !orientations_agree;  //flip
                                    must_complement = !must_complement;        //flip                           
                                }
                                
                                found_PERs___by_search_region___loop.at(per_ctr)->add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately(
                                                MiniRegion(    it_homol_reg->chromosome_of_homologous_region,
                                                            BOOST_Interval(  homolreg_lower_endpoint, homolreg_upper_endpoint ),
                                                            orientations_agree,
                                                            must_complement  )  );
                                                            
                                if (homolreg_upper_endpoint < 1000)                                
                                {
                                    std::stringstream error_strm;                                    
                                    error_strm << "\n\nERROR!   MR is at end of chromo!!!!"  
                                                << "\n\thomolreg_lower_endpoint = "  <<  homolreg_lower_endpoint  
                                                << "\n\thomolreg_upper_endpoint = "  <<  homolreg_upper_endpoint 
                                                << "\n\tit_homol_reg->chromosome_of_homologous_region = " << it_homol_reg->chromosome_of_homologous_region
                                                << "\n\tit_DAR_reg = " << *it_DAR_reg << "\n\n";                                                
                                    print_map_keys_and_values<uint,uint>(map_DARreg_to_some_homol_reg, "map_DARreg_to_some_homol_reg", &error_strm, true); 
                                    
                                    error_message(error_strm, false);                                
                                }
                                                            
                                                            
                                
                            }//not empty                                        
                                                                                                                                                                                                                                                            
                        } // check this homol reg  
		    }
                        
                } // homolreg on DAR
                
//                 std::cerr << "\n\tget all MR loop total:  time = " << omp_get_wtime() - time_piece << ".\n";
//                 time_piece = omp_get_wtime();
                                                
            } // per loop                                            
                                
            mpfr_free_cache();                              
        }//parallel                
        
    
//         std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "  condense from threads...\n";
    
        // NOT PARALLEL !!!!!!!!!
        //transfer to global list.  for this rank.
                        
                        
                        
                        
        //condense for this rank
        
        for (type_list_PER::const_iterator it_found_per = found_PERs___by_search_region.at(region_ctr).begin();
                it_found_per != found_PERs___by_search_region.at(region_ctr).end();
                ++it_found_per)
        {
            const std::pair<type_map_string_to_PER::iterator, bool> result_insert_it 
                            = PERs_this_Event__this_rank.insert( std::pair<std::string, Paired_end_read>(it_found_per->name, *it_found_per)   );
                            
            result_insert_it.first->second.deep_copy_from_serialize( *it_found_per, my_Conn_Comp );
            
            if (!result_insert_it.second) // no new element inserted -  already existed
	    {  result_insert_it.first->second.add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately( *it_found_per );  }
        }                        
        
        
        found_PERs___by_search_region.at(region_ctr).clear();
        
        
        mpfr_free_cache();
        
                                        
    }//region_ctr
    
    
    
    
    found_PERs___by_search_region.clear();
    
    
    
    
    
    

        

//     std::cerr << "\n\nrank " << my_MPI_rank << "      waiting at   all_reduce ...";

    if (size_MPI_WORLD > 1)
    {

        bool world_MPI__success_val;
        boost::mpi::all_reduce<bool>(  *world_MPI_boost_ptr, success_raw_upload, world_MPI__success_val, min_bool  );
        
        if (!world_MPI__success_val)
            return false;
    }
    else if (!success_raw_upload)
    {  return false;  }
        
    
    
    
    
    
    
    
    
    
    
    
    
    
//     std::cerr << "\n\nrank " << my_MPI_rank << "      sending PERs to rank 0 ...";        
    
    
    //condense for MPI world
    
    if (my_MPI_rank > 0 )
    {
        world_MPI_boost_ptr->send<type_map_string_to_PER>( 0, my_MPI_rank, PERs_this_Event__this_rank );
        PERs_this_Event__this_rank.clear();
    }
    else
    {
        //condense  from self  
        {
            type_map_string_to_PER::iterator it_self_PER = PERs_this_Event__this_rank.begin();
            while (  it_self_PER != PERs_this_Event__this_rank.end()  )
            {                        
                const std::pair<type_map_string_to_PER::iterator, bool> result_insert_it 
                                = PERs_on_this_profile.insert( std::pair<std::string, Paired_end_read>(it_self_PER->first, it_self_PER->second)   );
                                
                result_insert_it.first->second.reconnect_pointers_after_uploading_from_a_serialization( my_Conn_Comp );
                
                if (!result_insert_it.second) // no new element inserted -  already existed
                    result_insert_it.first->second.add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately(it_self_PER->second);                                          
            
                PERs_this_Event__this_rank.erase( it_self_PER++ );  //post-increment necessary!!!
            }  
        }        
        PERs_this_Event__this_rank.clear();
                
        
        
        //condense  from other ranks
        for (int rank = 1; rank < size_MPI_WORLD; ++rank )        
        {
            type_map_string_to_PER recvd_PERs_from_rank;            
            world_MPI_boost_ptr->recv<type_map_string_to_PER>(rank, rank, recvd_PERs_from_rank);
            
            
            
            //condense            
            type_map_string_to_PER::iterator it_recv_PER = recvd_PERs_from_rank.begin();
            while (  it_recv_PER != recvd_PERs_from_rank.end()  )
            {                        
                const std::pair<type_map_string_to_PER::iterator, bool> result_insert_it 
                                = PERs_on_this_profile.insert( std::pair<std::string, Paired_end_read>(it_recv_PER->first, it_recv_PER->second)   );
                                
                result_insert_it.first->second.deep_copy_from_serialize(  it_recv_PER->second, my_Conn_Comp  );                
                
                
                if (!result_insert_it.second) // no new element inserted -  already existed
		{  result_insert_it.first->second.add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately(it_recv_PER->second);  }
            
                recvd_PERs_from_rank.erase(it_recv_PER++);  //post-increment necessary!!!
            }                                
        } // rank    
        
    } // 0
    
    
    






    
//     std::cerr << "\n\nrank " << my_MPI_rank << "     at broadcast PERs_on_this_profile ...";    

    boost::mpi::broadcast<type_map_string_to_PER>( *world_MPI_boost_ptr, PERs_on_this_profile, 0);
           
    
    for (type_map_string_to_PER::iterator it_PER = PERs_on_this_profile.begin();
	    it_PER != PERs_on_this_profile.end();
	    ++it_PER)
    {  it_PER->second.reconnect_pointers_after_uploading_from_a_serialization(my_Conn_Comp);  }
    
 
    
    
//     world_MPI_boost_ptr->barrier();
    
    
    std::cerr << "\n\t\ttime to add MRs:  " << omp_get_wtime() - time_add_MRs_begin << "  s.\n";
    
    

    mpfr_free_cache(); 
    
    
    return   success_raw_upload;
    
}//upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals_____MPI_and_OMP
























void Event::for_all_PERs____set_region_sequences___and___orient___MPI_and_OMP()
{      
    
    const std::vector<type_map_string_to_PER::iterator>  found_PERs___loop(get_loopable_iterators<type_map_string_to_PER>(PERs_on_this_profile));    
              
    
    //prepare:
    
    #pragma omp parallel 
    {//parallel
    
        #pragma omp for  schedule(dynamic,20)
        for (uint j = 0; j < found_PERs___loop.size();  ++j) 
        {         
            // make sure everything is suff long
            found_PERs___loop.at(j)->second.extend_MRs_to_minimum_length();
            found_PERs___loop.at(j)->second.set_sequence_for_all_MiniRegions(uploaded_pieces_of_chromos); 
        }
        
        mpfr_free_cache();
    }//parallel
    
            
    
    //align as necessary
    
    set_up_and_align_PERs_using_GPU____OMP
			(PERs_on_this_profile,
			 this,//spawning_Event,
			 false,//create_and_consider_MiniEvents,
			 false,//consider_GeneConversion_breakpoints,  //irrelevant if "create_and_consider_MiniEvents" = false
			 false,//sum_over_all_unaffected_LCRs,  //"false" for heuristic
			 false,//special_save_breakpoints_actually_captured_in_alignment, //"true" for heuristic
			 true);//do_PER_orient_formatting)  // orthogonal to all other options    
    
    
    //finish orienting
    
    #pragma omp parallel 
    {//parallel
    
        #pragma omp for  schedule(dynamic,20)
        for (uint j = 0; j < found_PERs___loop.size();  ++j) 
        {         
	    found_PERs___loop.at(j)->second.orient_and_complement_each_mate_according_to_its_spawning_Event_profile();
        }
        
        mpfr_free_cache();
    }//parallel    
    
      
    mpfr_free_cache();     

}//orient





















bool Event::determine_relation_of_region_to_Event_profile
                                (bool &orientations_of_region_and_Event_profile_agree,
                                 bool &must_complement_region_to_get_to_Event_profile,
                                 const uint &chromo_of_region,
                                 const BOOST_Interval interval_of_region)
{
    //init
    orientations_of_region_and_Event_profile_agree = true;
    must_complement_region_to_get_to_Event_profile = false;
    
    
    if (chromo_of_region == chromos[0]  and  BOOST_overlap(interval_of_region, LCRs[0]) )
    {
        orientations_of_region_and_Event_profile_agree = true;
        must_complement_region_to_get_to_Event_profile = false;
        return true;
    }
    else if (chromo_of_region == chromos[0]  and  BOOST_overlap(interval_of_region, LCRs[1]) )
    {
        orientations_of_region_and_Event_profile_agree =    (recomb_type != recomb_class__Inv);
        must_complement_region_to_get_to_Event_profile =    (recomb_type == recomb_class__Inv);
        return true;
    }      
    else if (chromo_of_region == chromos[0] and  BOOST_overlap(interval_of_region, region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs)    )
    {
        orientations_of_region_and_Event_profile_agree = true;
        must_complement_region_to_get_to_Event_profile = false;
        return true;
    }       
    else
    {
        for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg = regions_homologous_to_directly_affected_region.begin();
                it_homol_reg != regions_homologous_to_directly_affected_region.end();
                ++it_homol_reg)      
            if ( it_homol_reg->chromosome_of_homologous_region ==  chromo_of_region 
                and  BOOST_overlap(interval_of_region, get_endpoints_of_keys_of_compressed_map(it_homol_reg->compressed_map_homologous_region_to_profile)) )
            {
                if ( recomb_type == recomb_class__Inv
                     and   BOOST_overlap( get_endpoints_of_values_of_compressed_map(it_homol_reg->compressed_map_homologous_region_to_profile), LCRs[1])    )
                {
                    orientations_of_region_and_Event_profile_agree = ! it_homol_reg->orientation_of_homologous_region_and_profile_agree;
                    must_complement_region_to_get_to_Event_profile = ! it_homol_reg->homologous_region_must_be_complemented_to_get_to_profile;                
                }
                else
                {
                    orientations_of_region_and_Event_profile_agree = it_homol_reg->orientation_of_homologous_region_and_profile_agree;
                    must_complement_region_to_get_to_Event_profile = it_homol_reg->homologous_region_must_be_complemented_to_get_to_profile;                                        
                }
                
                return true;
            }    
            
        return false;
    }//else
    
            
}  // get_relation_of_region_to_Event_profile





































































void Event::clear_all_data()
{
    global_interlocking_Events.clear();
    local_interlocking_Events.clear();
    remaining_coordinates_for_consideration.clear();
    my_original_Breakpoint_complex.clear_all_data();
    
    natural_poisson_intervals___gender_INsensitive.clear();
    base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive.clear();
    base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.clear();
    observed_read_depth_for_natural_poisson_intervals.clear();
    uploaded_pieces_of_chromos.clear();
    PERs_on_this_profile.clear();
    filter_for_NON_large_gap_of_neighbors_and_self__that_affect_Read_depth.clear();
    RD_affecting_neighbors_self.clear();
    all_NON_large_gap_intervals_of_profile.clear();
    large_gaps_of_profile.clear();                        
    variational_positions_on_LCR_immediately_following_small_insertion_in_other_LCR.clear();
                            
    regions_homologous_to_LCR_profile.clear();
    regions_homologous_to_directly_affected_region.clear();
    variational_positions_of_profile.clear();  
    compressed_map_LCR_to_profile__for_each_LCR.clear();
    compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.clear();   
    
    absolute_coordinates_inside_of_large_gaps_of_profile.first.clear();
    absolute_coordinates_inside_of_large_gaps_of_profile.second.clear();
        

}// clear_all_data








































void Event::upload_or_create_natural_poisson_intervals___gender_INsensitive()
{
    
    natural_poisson_intervals___gender_INsensitive.clear();
    
    
    char nat_pois_filename[option_size];
    std::sprintf(nat_pois_filename, "%s/%s/natural_poisson_intervals/natural_poisson_intervals_for_Event_%u",
                                                    data_directory.c_str(),
                                                    visitation_schedule_subdirname.c_str(),
                                                    UID );
                                                    
    FILE *in_file = std::fopen(nat_pois_filename, "r");         
    if (in_file == NULL)
    {
        std::stringstream warning_strm;
        print_line_of_markers("WARNING! ", &warning_strm);
        print_line_of_markers("(", &warning_strm);
        warning_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
        
        warning_strm << "\n\nWARNING!!!   unable to find/open natural poisson intervals  file for Event " << UID
                     << ",     checked filename = [" << nat_pois_filename << "].\n\n\n";        
                     
        print_line_of_markers(")", &warning_strm);
        std::fprintf(stderr, "\n\n%s\n\n", warning_strm.str().c_str() );
    
        partition_DAR_by_homology___and___create_natural_poisson_intervals___NOT__gender_sensitive();        
    }        
    else
    {

        uint number_of_natural_poisson_colors;
        
        std::fscanf(in_file, "number_of_colors = %u", &number_of_natural_poisson_colors);
        
        
        
        type_map_uint_to_uint_to_uint_to_longdouble::iterator insert_color_it = natural_poisson_intervals___gender_INsensitive.begin();
        
        for (uint color = 0; color < number_of_natural_poisson_colors; ++color)
        {
            insert_color_it = natural_poisson_intervals___gender_INsensitive.insert(  insert_color_it,
                                                                    std::pair<uint, type_map_uint_to_uint_to_longdouble>
                                                                    (color, type_map_uint_to_uint_to_longdouble()  )    );
                                                                    
            uint number_of_chromos_this_color;                                                     
            std::fscanf(in_file, "\ncolor:  %*d      number_chromos =  %u", &number_of_chromos_this_color);
            
            
            
            
            type_map_uint_to_uint_to_longdouble::iterator insert_chr_it = insert_color_it->second.begin();
            
            for (uint chr_ctr = 0; chr_ctr < number_of_chromos_this_color; ++chr_ctr)
            {          
                
                uint compressed_size;
                uint chr_val;
                
                std::fscanf(in_file, "\n\t\tchromo  %u       compressed_size = %u", &chr_val, &compressed_size);
                                        
                
                insert_chr_it = insert_color_it->second.insert( insert_chr_it,
                                                                std::pair<uint, type_map_uint_to_longdouble>
                                                                        (chr_val, type_map_uint_to_longdouble()  )     );
                
                                                                
                                                                        
                                                                        
                                                                        
                                                                        
                type_map_uint_to_longdouble::iterator insert_coord_it =  insert_chr_it->second.begin();      
        
                
                type_vector_vector_int compressed_map_coords;
    
                compressed_map_coords.push_back( type_vector_int() );
                compressed_map_coords.push_back( type_vector_int() );                        
                compressed_map_coords.at(0).reserve( compressed_size );
                compressed_map_coords.at(1).reserve( compressed_size );
                                                            
                for ( uint j=0; j < compressed_size; ++j )
                {
                    int in_key, in_val;
                    std::fscanf(in_file, "\n\t\t\t%d      %d", &in_key, &in_val);
                    compressed_map_coords.at(0).push_back(in_key);
                    compressed_map_coords.at(1).push_back(in_val);                
                }
    
    
                const type_map_uint_to_uint full_coords_as_map(     convert_compressed_map_to_full_map( compressed_map_coords )    );
        
                for (type_map_uint_to_uint::const_iterator it_coord = full_coords_as_map.begin();
                        it_coord != full_coords_as_map.end();
                        ++it_coord)
                    insert_coord_it = insert_chr_it->second.insert(  insert_coord_it,  type_uint__longdouble( it_coord->first,  0.00L)   );                                                                  
            } // chr                                
        }//nat
        
        
        std::fclose( in_file );
                                
                
//         in_fs.close();
        
//         std::cerr << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n\n   DONE readsing nat poisson intervals...\n";
        
                    
        if (my_MPI_rank == 0)  
        {
            std::stringstream out_strm;       
            for (type_map_uint_to_uint_to_uint_to_longdouble::const_iterator it_nat = natural_poisson_intervals___gender_INsensitive.begin();
                    it_nat != natural_poisson_intervals___gender_INsensitive.end();
                    ++it_nat)      
            {     
                uint size_of_nat = 0;
                for (type_map_uint_to_uint_to_longdouble::const_iterator it_chr = it_nat->second.begin();
                        it_chr != it_nat->second.end();
                        ++it_chr)
                    size_of_nat += it_chr->second.size();            
                    
                out_strm << "\t\t\t\tcolor:  " << it_nat->first 
                        << "\t\tsize:  " << size_of_nat << "\n";
    //                     <<"\t\texample coords: ["  << *newest_interval.begin() 
    //                     <<", " << *(--newest_interval.end())
    //                     <<"]\n";
            }
        
            std::fprintf(stderr, "\n%s\n", out_strm.str().c_str() );           
        }        
    
    } // found it.








//     //make gender sensitive:
// 
//     if (gender_of_individual_is_female)  // erase Y
//     {
//         type_map_uint_to_uint_to_uint_to_longdouble::iterator it_nat = natural_poisson_intervals___gender_INsensitive.begin();
//         while ( it_nat != natural_poisson_intervals___gender_INsensitive.end() )
//         {
//             type_map_uint_to_uint_to_longdouble::iterator it_chr = it_nat->second.begin();
//             while (it_chr != it_nat->second.end())
// 	    {
//                 if (it_chr->first == 24)
//                     it_nat->second.erase( it_chr++ );
//                 else
//                     ++it_chr;
// 	    }
//             
//             if (it_nat->second.empty())
//                 natural_poisson_intervals___gender_INsensitive.erase( it_nat++ );
//             else
//                 ++it_nat;                    
//         }             
//     }//female





}  // upload_or_create_natural_poisson_intervals___gender_INsensitive

































void Event::ignore_regions_of_the_genome_from_consideration_by_this_Event
                    (const type_map_uint_to_list_BI  &regions_to_ignore)
{        
        
    
    {//erase from DAR
        const type_map_uint_to_list_BI::const_iterator it_find_chr = regions_to_ignore.find(  chromos[0]  );
        if ( it_find_chr != regions_to_ignore.end() )
        {
            for (type_list_BI::const_iterator it_bi = it_find_chr->second.begin();
                    it_bi != it_find_chr->second.end();
                    ++it_bi)
            {                        
                type_set_uint::iterator it_rem = remaining_coordinates_for_consideration.begin();
                while (it_rem != remaining_coordinates_for_consideration.end())
		{
                    if (  BOOST_in(*it_rem, *it_bi)  )			
                        remaining_coordinates_for_consideration.erase( it_rem++ );
                    else
                        ++it_rem;
		}
            }
            
            
            bool affects_the_LCRs = false;
            for (type_list_BI::const_iterator it_bi = it_find_chr->second.begin();
                    it_bi != it_find_chr->second.end();
                    ++it_bi)   
	    {
                if (  BOOST_overlap(*it_bi, LCRs[0])   or   BOOST_overlap(*it_bi, LCRs[1])  )
                {
                    affects_the_LCRs = true;
                    break;                
                }
	    }
            
            
            
            if (affects_the_LCRs)
            {            
//                for (uint qs=0; qs<2; ++qs)        
//                    compressed_map_LCR_to_profile__for_each_LCR.at(qs)
//                        =  convert_full_map_to_compressed_map(
//                                get_map_intersection_of_compressed_map_keys_and_set(
//                                    compressed_map_LCR_to_profile__for_each_LCR.at(qs),
//                                    remaining_coordinates_for_consideration  ));
                
                for (uint qs=0; qs<2; ++qs) 
		{
                    compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)
                        =  convert_full_map_to_compressed_map(
                                get_map_intersection_of_compressed_map_keys_and_set(
                                    compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs),
                                    remaining_coordinates_for_consideration  ));     
		}
                                                                                                    
                                
                variational_positions_of_profile.clear();  
                
                for (uint qs=0; qs<2; ++qs)
                {
                    const type_set_uint  varpos_that_remain_LCR_qs(
                                            extract_keys_of_compressed_map_and_return_as_set(
                                                get_compressed_map_inverse_of_compressed_map( 
                                                        compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs)   )));
                                        
                    variational_positions_of_profile.insert(  varpos_that_remain_LCR_qs.begin(), varpos_that_remain_LCR_qs.end() );                                
                }
                        

                type_list_Coords_and_homologous_region::iterator it_homol = regions_homologous_to_directly_affected_region.begin();
                while (  it_homol != regions_homologous_to_directly_affected_region.end()  )
                {
                    it_homol->compressed_map_homologous_region_to_profile
                            =  get_compressed_map_inverse_of_compressed_map(  //homol reg to DAR
                                    convert_full_map_to_compressed_map(
                                        get_map_intersection_of_compressed_map_keys_and_set(
                                            get_compressed_map_inverse_of_compressed_map(  it_homol->compressed_map_homologous_region_to_profile  ),  //DAR to homol reg
                                            remaining_coordinates_for_consideration)                ));
                                            
                    if (it_homol->compressed_map_homologous_region_to_profile.at(0).empty())                        
                        it_homol = regions_homologous_to_directly_affected_region.erase(it_homol);
                    else
                        ++it_homol;
                }//regions_homologous_to_directly_affected_region    
            
            }//affects_the_LCRs
            
            
        }//it_find_chr -  chromos[0]
        
    }//erase from DAR
    
    
    
    
    
    
    {//erase from other  homol regs
        type_list_Coords_and_homologous_region::iterator it_homol = regions_homologous_to_directly_affected_region.begin();
        while (  it_homol != regions_homologous_to_directly_affected_region.end()  )
        {
            const type_map_uint_to_list_BI::const_iterator it_find_chr = regions_to_ignore.find(  it_homol->chromosome_of_homologous_region  );
            if (it_find_chr == regions_to_ignore.end() )
                ++it_homol;
            else
            {
                for (type_list_BI::const_iterator it_bi = it_find_chr->second.begin();
                        it_bi != it_find_chr->second.end();
                        ++it_bi)   
                    if (   BOOST_overlap(  *it_bi, get_endpoints_of_keys_of_compressed_map(it_homol->compressed_map_homologous_region_to_profile) )    )
                    {
                        it_homol->compressed_map_homologous_region_to_profile =                         
                                            get_compressed_map_intersection_of_map_keys_and_COMPLEMENT_of_interval__using_compressed_map(
                                                it_homol->compressed_map_homologous_region_to_profile,
                                                *it_bi);                                                                                                                                         
                    }
                    
                if (it_homol->compressed_map_homologous_region_to_profile.at(0).empty())                        
                    it_homol = regions_homologous_to_directly_affected_region.erase(it_homol);
                else
                    ++it_homol;
            }// found chr
        }//regions_homologous_to_directly_affected_region      
    
    }//erase from other  homol regs            



} // ignore_regions_of_the_genome_from_consideration_by_this_Event

























//.first = GC sensitive and information-sensitive.
//.second = only GC sensitive.
type_longdouble__longdouble  Event::calculate_expected_number_of_PERs_at_breakpoint
			    (const BOOST_Interval &profile_breakpoint,
			    const MiniRegion &mini_LCR_0,
			    const MiniRegion &mini_LCR_1,
			     const bool  &AB_ABA__false_____BA_BAB__true)  const
{
    
    const bool is_NAHR = BOOST_is_a_point(profile_breakpoint);
    
    const MiniRegion *MR_LCRs[2];
    const BOOST_Interval* full_LCRs[2];
    type_set_uint hybrid_varpos_indeces__0_1[2];
    BOOST_Interval absolute_brkpt_interval__on_LCR[2];
    type_set_uint absolute_variational_positions__on_LCR[2];    
    BOOST_Interval outside_LCR__0_1[2];    
    uint index_to_get_to_true_index[2] = {0, 1};
	    
    
    MR_LCRs[0] = &mini_LCR_0;
    MR_LCRs[1] = &mini_LCR_1;
    
    full_LCRs[0] = &LCRs[0];
    full_LCRs[1] = &LCRs[1];    
        
            
    {//get some breakpoint info
	BOOST_split_overlapping_intervals_into_lower_region_and_upper_region(
			    convert_each_profile_breakpoint_to_absolute_breakpoints(profile_breakpoint, compressed_map_LCR_to_profile__for_each_LCR),
			    absolute_brkpt_interval__on_LCR);

	for (uint qs=0; qs<2; ++qs)
	{
	    const type_map_uint_to_uint  map_profile_varpos_only_to_mini
			= convert_compressed_map_to_full_map(
				get_compressed_map_inverse_of_compressed_map(
				    get_compressed_map_intersection_of_map_keys_and_interval__using_compressed_map(
					    compressed_map_LCR_to_variational_positions_of_profile__for_each_LCR.at(qs),
					    MR_LCRs[qs]->region_interval)  ));
					    
	    absolute_variational_positions__on_LCR[qs] 
			= extract_values_of_map_and_return_as_set<uint,uint>(map_profile_varpos_only_to_mini);
			
	    absolute_variational_positions__on_LCR[qs].insert(
		    pair_at<type_set_uint>(absolute_coordinates_inside_of_large_gaps_of_profile, qs).lower_bound(MR_LCRs[qs]->region_interval.lower()),
		    pair_at<type_set_uint>(absolute_coordinates_inside_of_large_gaps_of_profile, qs).upper_bound(MR_LCRs[qs]->region_interval.upper()));
	}//qs
    }//get some breakpoint info 
	
      
    
    // recall that "get_variational_indeces_mapped_onto_hybrid_sequence"  and  "get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint" are designed so that you the first input sequence always comes first.  i.e. it does not know about the type of hybrid AB vs BA  (or ABA vs BAB).         
    
    //rearrange, as appropriate:
    if (AB_ABA__false_____BA_BAB__true)
    {
	std::swap<const MiniRegion*>(MR_LCRs[0], MR_LCRs[1]);
	std::swap<const BOOST_Interval*>(full_LCRs[0], full_LCRs[1]);
	std::swap<type_set_uint>(hybrid_varpos_indeces__0_1[0], hybrid_varpos_indeces__0_1[1]);
	std::swap<BOOST_Interval>(absolute_brkpt_interval__on_LCR[0], absolute_brkpt_interval__on_LCR[1]);
	std::swap<type_set_uint>(absolute_variational_positions__on_LCR[0], absolute_variational_positions__on_LCR[1]);
	std::swap<uint>(index_to_get_to_true_index[0], index_to_get_to_true_index[1]);
    }

    

    bool must_be_reversed__0_1[2];
    bool must_complement__0_1[2];    
    
    
    if (!AB_ABA__false_____BA_BAB__true)
    {//AB, ABA
	must_be_reversed__0_1[0] = false;
	must_be_reversed__0_1[1] = !(recomb_type == recomb_class__DupDel);    
	
	must_complement__0_1[0] = false;
	must_complement__0_1[1] = (recomb_type == recomb_class__Inv);	
	
    }//AB, ABA
    else
    {//BA, BAB
	must_be_reversed__0_1[0] = !(recomb_type == recomb_class__DupDel);
	must_be_reversed__0_1[1] = false;
	
	must_complement__0_1[0] = (recomb_type == recomb_class__Inv);
	must_complement__0_1[1] = false;
    }//BA, BAB    
    
    get_variational_indeces_mapped_onto_hybrid_sequence
	    (*MR_LCRs[0],
	    absolute_brkpt_interval__on_LCR[0],    //not to be included!!
	    must_be_reversed__0_1[0],   //always want beginning of A
	    *MR_LCRs[1],
	    absolute_brkpt_interval__on_LCR[1],
	    must_be_reversed__0_1[1],
	    absolute_variational_positions__on_LCR[0],
	    absolute_variational_positions__on_LCR[1],
	    is_NAHR,
	    hybrid_varpos_indeces__0_1[0],
	    hybrid_varpos_indeces__0_1[1]);	
	    

	    
	    
//     std::cerr << "\n\nabsolute_brkpt_interval__on_LCR[0] = " << absolute_brkpt_interval__on_LCR[0] 
// 	    << "\nabsolute_brkpt_interval__on_LCR[1] = " << absolute_brkpt_interval__on_LCR[1] 
// 	    << "\n\n";
	    
//     print_set<uint>(absolute_variational_positions__on_LCR[0], "absolute_variational_positions__on_LCR[0]", NULL, true);
//     print_set<uint>(absolute_variational_positions__on_LCR[1], "absolute_variational_positions__on_LCR[1]", NULL, true);
	    
//     print_set<uint>(hybrid_varpos_indeces__0_1[0], "hybrid_varpos_indeces__0_1[0]", NULL, true);
//     print_set<uint>(hybrid_varpos_indeces__0_1[1], "hybrid_varpos_indeces__0_1[1]", NULL, true);
	        	    	    

    const type_string__BI hybrid_seq_and_brkpt_indeces(
				get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint(
					*MR_LCRs[0],
					absolute_brkpt_interval__on_LCR[0],
					must_complement__0_1[0],
					must_be_reversed__0_1[0],
					*MR_LCRs[1],
					absolute_brkpt_interval__on_LCR[1],
					must_complement__0_1[1],
					must_be_reversed__0_1[1],
					is_NAHR));   		    					

    const type_vector_int absolute_labels_of_hybrid_indeces(
				    label_absolute_coordinates_on_hybrid
					(*MR_LCRs[index_to_get_to_true_index[0]],
					absolute_brkpt_interval__on_LCR[index_to_get_to_true_index[0]],
					*MR_LCRs[index_to_get_to_true_index[1]],
					absolute_brkpt_interval__on_LCR[index_to_get_to_true_index[1]],
					is_NAHR,
					AB_ABA__false_____BA_BAB__true,
					recomb_type));
					
    type_map_uint_to_uint absolute_positions_to_hybrid_indeces;
    type_map_uint_to_uint::iterator insert_it_abs = absolute_positions_to_hybrid_indeces.begin();
    for (uint j=0; j < absolute_labels_of_hybrid_indeces.size(); ++j)
    {
	insert_it_abs = absolute_positions_to_hybrid_indeces.insert(insert_it_abs, type_uint__uint((uint)absolute_labels_of_hybrid_indeces.at(j), j));
    }
	
//     print_map_keys_and_values<uint,uint>(absolute_positions_to_hybrid_indeces, "absolute_positions_to_hybrid_indeces", NULL, true);


	std::cerr << "hybrid_seq_and_brkpt_indeces.second = " << hybrid_seq_and_brkpt_indeces.second << "\n";
	print_set<uint>(hybrid_varpos_indeces__0_1[0], "hybrid_varpos_indeces__0_1[0]");
	print_set<uint>(hybrid_varpos_indeces__0_1[1], "hybrid_varpos_indeces__0_1[1]");
	
	
	
    
    const BOOST_Interval hybrid_indeces_below_LCR_0(
				get_positions_of_hybrid_lying_BELOW_LCR(absolute_positions_to_hybrid_indeces, LCRs[0]));
				
    const BOOST_Interval hybrid_indeces_inbetween_LCRs(
				get_positions_of_hybrid_lying_INBETWEEN_LCRs(absolute_positions_to_hybrid_indeces, 
									     region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs));

    const BOOST_Interval hybrid_indeces_above_LCR_1(
				get_positions_of_hybrid_lying_ABOVE_LCR(absolute_positions_to_hybrid_indeces, LCRs[1], chromos[1]));
    

				
				
				
				
    for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
	    it_rg != Readgroup_stats.end();
	    ++it_rg)
    {  std::cerr << "\n\t" << it_rg->first << "  -  average mate length  = " << it_rg->second.average_mate_length << "\n";  }
		
		       
    longdouble expected_fragmentation_rate___GC_sensitive_and_information_aware__ALL_RG = 0.00L;  // GC and mates must display breakpoint
    longdouble expected_fragmentation_rate___GC_sensitive_only = 0.00L;  // GC only - mates need not display both breakpoints, but just overlap
//     longdouble expected_fragmentation_rate___totally_naive = 0.00L;
    longdouble expected_fragmentation_rate___informative_but_naive_coverage = 0.00L;
        
    for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
	    it_rg != Readgroup_stats.end();
	    ++it_rg)
    {			
	    const longdouble lambda_per_haploid__0_this_group = ((longdouble)it_rg->second.number_mapped_Proper_Pairs_reads / 2) / (length_of_haploid_genome*2); 
	    
	const uint mean = (uint)it_rg->second.median_insert_size;
	const uint sd = (uint)it_rg->second.median_insert_size____absolute_standard_deviation;
	//const BOOST_Interval test_frag_range(safe_subtract_base_1(mean, sd*4), mean + 4*sd);
	const BOOST_Interval test_frag_range(BOOST_make_point(mean));
	
	boost::math::normal_distribution<longdouble> normal_dist__frag(mean, sd);
	
	for (uint frag_len = test_frag_range.lower(); frag_len <= test_frag_range.upper(); ++frag_len)
	{	
		const longdouble cont_correction_prob =  1;  //boost::math::cdf<longdouble>(normal_dist__frag, (longdouble)frag_len + 0.5)
							     //- boost::math::cdf<longdouble>(normal_dist__frag, (longdouble)frag_len - 0.5);		
	    
	    
	    int leftmost_index = std::max<int>(0, (int)hybrid_seq_and_brkpt_indeces.second.lower() - (frag_len-1));	
	    while (leftmost_index < hybrid_seq_and_brkpt_indeces.second.lower())
	    {	    	    
		BOOST_Interval mates[2];
		mates[0].set( (uint)leftmost_index, (uint)leftmost_index + it_rg->second.average_mate_length- 1);
		mates[1].set( (uint)leftmost_index + frag_len - 1 - (it_rg->second.average_mate_length- 1),
			    (uint)leftmost_index + frag_len - 1);
		
			    
		const bool  informs_LCR_0 =  (  !get_intersection_of_set_and_interval(hybrid_varpos_indeces__0_1[0], mates[0]).empty()
						or   !get_intersection_of_set_and_interval(hybrid_varpos_indeces__0_1[0], mates[1]).empty() );
	    
		const bool  informs_LCR_1 =  (  !get_intersection_of_set_and_interval(hybrid_varpos_indeces__0_1[1], mates[0]).empty()
						or   !get_intersection_of_set_and_interval(hybrid_varpos_indeces__0_1[1], mates[1]).empty() );
			    
						
		const bool  intersects_below_LCR_0(BOOST_overlap(mates[0],hybrid_indeces_below_LCR_0)  or  BOOST_overlap(mates[1],hybrid_indeces_below_LCR_0));
		const bool  intersects_between_LCRs(BOOST_overlap(mates[0],hybrid_indeces_inbetween_LCRs)  or  BOOST_overlap(mates[1],hybrid_indeces_inbetween_LCRs));
		const bool  intersects_above_LCR_1(BOOST_overlap(mates[0],hybrid_indeces_above_LCR_1)  or  BOOST_overlap(mates[1],hybrid_indeces_above_LCR_1));
				
		const uint gc_count = calculate_GC_content_of_sequence(hybrid_seq_and_brkpt_indeces.first.substr(leftmost_index, it_rg->second.median_insert_size));
														    
// 		expected_fragmentation_rate___GC_sensitive_only += it_rg->second.haploid_fragmentation_rate__per_GC_count.at(gc_count);
			
// 		std::cerr << "\t" <<  informs_LCR_0 << "  " << informs_LCR_1 << "  " << intersects_below_LCR_0 << "  " << intersects_between_LCRs << "  "  << intersects_above_LCR_1 << "\n";
		
		if ((informs_LCR_0  and  informs_LCR_1)
		    or   (informs_LCR_0  and  (intersects_between_LCRs or intersects_above_LCR_1))
		    or   (informs_LCR_1  and  (intersects_between_LCRs or intersects_below_LCR_0)) )  //then is informative of this breakpoint.
		{							
		    expected_fragmentation_rate___GC_sensitive_and_information_aware__ALL_RG
				+= cont_correction_prob * it_rg->second.haploid_fragmentation_rate__per_GC_count.at(gc_count);
		}
		
		if (BOOST_subset(hybrid_seq_and_brkpt_indeces.second, BOOST_Interval(mates[0].lower(), mates[1].upper())))
		{
			expected_fragmentation_rate___GC_sensitive_only += cont_correction_prob * it_rg->second.haploid_fragmentation_rate__per_GC_count.at(gc_count);
// 			expected_fragmentation_rate___totally_naive += cont_correction_prob * haploid_coverage;
			expected_fragmentation_rate___informative_but_naive_coverage += cont_correction_prob * lambda_per_haploid__0_this_group;
		}
		
		
		++leftmost_index;
	    }//relevant				    
	}//frag_len
	
    }//rg
    
    std::cerr << "internal:  expected_fragmentation_rate___informative_but_naive_coverage = " << expected_fragmentation_rate___informative_but_naive_coverage << "  (without RD_zygous multiplier)\n";
//     std::cerr << "internal:  expected_fragmentation_rate___totally_naive = " << expected_fragmentation_rate___totally_naive << "  (without RD_zygous multiplier)\n";
    
    return  type_longdouble__longdouble(expected_fragmentation_rate___GC_sensitive_and_information_aware__ALL_RG, expected_fragmentation_rate___GC_sensitive_only);

}//calculate_expected_number_of_PERs_at_breakpoint
























void Event::set_absolute_coordinate_mapping_into_large_gaps()
{
    for (uint qs=0; qs<2; ++qs)
    {    
	type_set_uint  profile_coords_in_LGs;
	type_set_uint::iterator insert_it_prof = profile_coords_in_LGs.begin();
	for (type_set_BI::const_iterator it_lg = large_gaps_of_profile.begin();
		it_lg != large_gaps_of_profile.end();
		++it_lg)
	{	
	    for (uint j=it_lg->lower();  j<= it_lg->upper(); ++j)
	    {  insert_it_prof = profile_coords_in_LGs.insert(insert_it_prof, j);  }    
	}
	
	pair_at<type_set_uint>(absolute_coordinates_inside_of_large_gaps_of_profile, qs)
	    =  extract_keys_of_map_and_return_as_set<uint,uint>(
		    get_map_intersection_of_map_keys_and_set<uint,uint>(
			    convert_compressed_map_to_full_map(
					get_compressed_map_inverse_of_compressed_map(compressed_map_LCR_to_profile__for_each_LCR.at(qs))),
			    profile_coords_in_LGs));
    }//qs

}//set_absolute_coordinate_mapping_into_large_gaps













void Sampled_diploid_Event_data::print_this_Sampled_diploid_Event_data
				(std::stringstream *const &some_ss)  const
{
    
    std::stringstream output_str;
    
    bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str;  
    
    print_line_of_markers("[[[", ss_ptr);    
    
    (*ss_ptr) << "\n\n\n\n\nprinting Sampled_diploid_Event_data:\n\n"; 
    (*ss_ptr) 	<< "\t\tcc_of_event = " << cc_of_event
		<< "\n\t\tevent_UID = " << event_UID
		<< "\n\t\tthe_diploid_profile_brkpts = " << the_diploid_profile_brkpts.first << "   ,    " << the_diploid_profile_brkpts.second
		<< "\n\t\tthe_diploid_Event_outcome = " << convert_haploid_outcome_to_string(the_diploid_Event_outcome.first) 
			<< ",  " << convert_haploid_outcome_to_string(the_diploid_Event_outcome.second)
		<< "\n\t\tP_diploid_outcome = " << P_diploid_outcome
		<< "\n\t\tP_diploid_breakpoint = " << P_diploid_breakpoint
		<< "\n\n";
    

    print_line_of_markers("]]]", ss_ptr);    
    

    if (!provided_with_valid_sstream)
    {  std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );  }
    
    std::cerr.flush();std::fflush(stderr);
      
} // end of print_this_entry



























void Breakpoint_complex::print_this_Breakpoint_complex
					(std::stringstream *const &some_ss) const
{
    
    std::stringstream output_str;
    
    bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str;  
      
    
    (*ss_ptr) << "\nBreakpoint complex:";   
    
    print_set<uint>(variational_positions__ALL, "variational_positions__ALL", ss_ptr, false);
    print_set<BOOST_Interval, compare_BI>(breakpoints_NAHR__AB, "breakpoints_NAHR__AB", ss_ptr, false);
    print_set<BOOST_Interval, compare_BI>(breakpoints_NAHR__BA, "breakpoints_NAHR__BA", ss_ptr, false);
    print_set<BOOST_Interval, compare_BI>(breakpoints_GeneConversion___ABA, "breakpoints_GeneConversion___ABA", ss_ptr, false);
    print_set<BOOST_Interval,compare_BI>(breakpoints_GeneConversion___BAB, "breakpoints_GeneConversion___BAB", ss_ptr, false);
    
    (*ss_ptr) << "\n\tsizes:"
	    << "\n\t\tvarpos:  " << variational_positions__ALL.size()
	    << "\n\t\tAB:  " << breakpoints_NAHR__AB.size()
	    << "\n\t\tBA:  " << breakpoints_NAHR__BA.size()
	    << "\n\t\tABA:  " << breakpoints_GeneConversion___ABA.size()
	    << "\n\t\tBAB:  " << breakpoints_GeneConversion___BAB.size()
	    << "\n";        

    if (!provided_with_valid_sstream)
    {  std::cerr << "\n" << ss_ptr->str() << "\n";  }
      
} // end of print_this_entry








void Breakpoint_complex::clear_all_data()
{
    variational_positions__ALL.clear(); 

    breakpoints_NAHR__AB.clear(); 
    breakpoints_NAHR__BA.clear(); 
    
    breakpoints_GeneConversion___ABA.clear(); 
    breakpoints_GeneConversion___BAB.clear(); 

}//clear_all_data













