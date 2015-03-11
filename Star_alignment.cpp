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


#include <mpreal.h>
#include <mpfr.h>



#include <api/BamAlignment.h>



#include <general_typedefs.h>
#include <translations_through_alignment.h>
#include <templates.h>


#include "globals.h"
#include "Star_alignment.h"
#include "PER_bundle_wrapper.h"



static  std::set<char>  initialize_allowable_nucleotides()
{
    std::set<char>  dummy_allowable_nucleotide_values;
    dummy_allowable_nucleotide_values.insert('A');
    dummy_allowable_nucleotide_values.insert('C');
    dummy_allowable_nucleotide_values.insert('G');
    dummy_allowable_nucleotide_values.insert('T');
    dummy_allowable_nucleotide_values.insert('N');
    
    return dummy_allowable_nucleotide_values;
}
static std::set<char>  allowable_nucleotide_values = initialize_allowable_nucleotides();








void Star_alignment::print_this_Star_alignment(std::stringstream *const &some_ss)  const
{
    
    std::stringstream output_str;
    
    const bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str; 

    
    (*ss_ptr) << "\nprinting Star Alignment...\n";
    
    print_line_of_markers("[[", ss_ptr);


    (*ss_ptr) << "\n\nprogressive_star_alignment.size() = " << progressive_star_alignment.size() 
                <<"\nstar_profile_size = " << star_profile_size
                << "\nstar_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
                << "\nstar_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = " << star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX
                << "\nactivated_star_alignment = " << activated_star_alignment 
                << "\n\n";
		
    (*ss_ptr) << "\nLCR_0__LCR_1__XY___YX___XYX___YXY = " << LCR_0__LCR_1__XY___YX___XYX___YXY << "\n";
                
    print_vector<std::string>(names_of_rows, "names_of_rows", ss_ptr, true );
       
    print_vector<int>(progressive_star_alignment.at(1), "progressive_star_alignment.foundation_seq", ss_ptr, true);        
    
    std::cerr << "\n\n\nfailed_alignment = " << failed_alignment << "\n\n";

    print_line_of_markers("]]", ss_ptr);

    if (!provided_with_valid_sstream)
    {
        std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );       
        std::cerr.flush();std::fflush(stderr);
    }
    
    
} // print_this_Star_alignment















void Star_alignment::init_Star_alignment
			    (const type_vector_int &in_foundation_coords,
			    const std::string &in_foundation_seq,
			    const std::string &foundation_seq_name,
			    const BOOST_Interval &foundation_index_brkpts,
			    const bool &is_NAHR,
			    const int &in__LCR_0__LCR_1__XY___YX___XYX___YXY)
{
    clear_all_data();
    
    activated_star_alignment = true;
    star_profile_size = in_foundation_seq.size();
    failed_alignment = false;
    
    LCR_0__LCR_1__XY___YX___XYX___YXY = in__LCR_0__LCR_1__XY___YX___XYX___YXY;
    
    
        
    star_profile_index_of_first_Y_nucleotide_in_hybrid_XY = foundation_index_brkpts.lower();
    if (is_NAHR)
    {  star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = -1;  }
    else
    {  star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = foundation_index_brkpts.upper();  }
	    
	    
    // You create a desired MiniRegion of the genome "foundation", possibly a hybrid.
    // If it is a Gene Converison hybrid, then at least one of "foundation_index_of_first_Y_nucleotide_of_hybrid_XY" 
    //	or "foundation_index_of_first_second_X_nucleotide_in_hybrid_XYX" is in (0, foundation_length)
    // If it is an NAHR hybrid, then "foundation_index_of_first_Y_nucleotide_of_hybrid_XY" is in (0, foundation_length).
    //		In the NAHR case, then "star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX" shouldf be set to "foundation_length".
    if (star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX < 0)
    {  star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = in_foundation_seq.size();  }
	    
    assert(star_profile_index_of_first_Y_nucleotide_in_hybrid_XY <= star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX);
            
    
    std::string coords_name("coordinates");    
    
    progressive_star_alignment.push_back(in_foundation_coords);                                    
    progressive_star_alignment.push_back(convert_string_to_ints(in_foundation_seq));            
            
            
            
    if (   (star_profile_size  !=  (uint)in_foundation_coords.size())
	or   ! ( -1 <= star_profile_index_of_first_Y_nucleotide_in_hybrid_XY  and  star_profile_index_of_first_Y_nucleotide_in_hybrid_XY <= (int)star_profile_size )
	or   ! ( -1 <= star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX  and  star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX <= (int)star_profile_size ))
    {
	std::stringstream error_strm;                
	error_strm << "\n\nERROR!   in  \"Star_alignment\" constructor...\n\n";                
	error_strm << "star_profile_size = " << star_profile_size
		    << "\nstar_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
		    << "\nstar_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = " << star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX
		    << "\n(uint)in_foundation_coords.size() = " << (uint)in_foundation_coords.size()
		    << "\nfoundation_seq_name = [" << foundation_seq_name << "]\n"
		    << "\n(star_profile_size  !=  (uint)in_foundation_coords.size()) = " << (star_profile_size  !=  (uint)in_foundation_coords.size())
		    << "\n(-1 <= star_profile_index_of_first_Y_nucleotide_in_hybrid_XY) = " << (-1 <= star_profile_index_of_first_Y_nucleotide_in_hybrid_XY)
		    << "\n(star_profile_index_of_first_Y_nucleotide_in_hybrid_XY <= (int)star_profile_size) = " << (star_profile_index_of_first_Y_nucleotide_in_hybrid_XY <= (int)star_profile_size)
		    << "\nfoundation_index_brkpts = " << foundation_index_brkpts
		    << "\n\nin_foundation_seq:\n" << in_foundation_seq << "\n\n";
		    
	print_vector<int>(in_foundation_coords, "in_foundation_coords", &error_strm, true);                                
	error_message(error_strm, false);                   	
	
	coords_name = "!!!!!!!!!!!!ERROR-IN-STAR_ALIGNMENT-CONSTRUCTOR!!!!!!!!!!!!!!!!!!";
	
	const int diff_in_sizes = (int)in_foundation_coords.size() - star_profile_size;
	if ( diff_in_sizes > 0)
	{
	    for (int j=0; j < diff_in_sizes; ++j)
	    {  progressive_star_alignment.at(1).push_back(3);  }
	}
	else if (diff_in_sizes < 0)
	{
	    for (int j=0; j < (-1)*diff_in_sizes; ++j)
	    {  progressive_star_alignment.at(1).pop_back();  }                               
	}
	else 
	{
	    star_profile_index_of_first_Y_nucleotide_in_hybrid_XY = (int)star_profile_size;
	    star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = (int)star_profile_size;
	}
	
	failed_alignment = true;	
    }//sanity check
    

            
    //dummy values
    names_of_rows.push_back( coords_name );
    names_of_rows.push_back( foundation_seq_name );    
    
    //dummy values
    vector_PER_info.push_back("coordinates");
    vector_PER_info.push_back("genome sequence");
    
    //dummy
    vector_PER_base_Qualities.push_back(std::string(" "));
    vector_PER_base_Qualities.push_back(std::string(" "));
    
    
}//init_Star_alignment

























//"absolute_coords_hybrid_X_to_pairwise_algmt_profile"  and  "absolute_coords_hybrid_Y_to_pairwise_algmt_profile"  are interchangeable.  In fact,t hey could just be combined into one map, if desired.
bool Star_alignment::add_pairwise_alignment_to_progressive_star_alignment
                                   (const real &P_PER_cond_this_sequence,
                                    const type_string__string &PER_and_hybrid_algmt,
                                    const type_map_uint_to_uint &absolute_coords_hybrid_X_to_pairwise_algmt_profile,
                                    const type_map_uint_to_uint &absolute_coords_hybrid_Y_to_pairwise_algmt_profile,                             
                                    const bool &orientation_of_pairwise_algmt_agrees_with_foundation_sequence,
                                    const bool &must_complement_algmt_to_get_to_foundation_sequence,
				    type_Balgmt__Balgmt  mates_of_PER, //NOT BY REFERENCE
				    const real &sum_over_all_P_locations)
{             
    assert(activated_star_alignment);
    
    if (failed_alignment)
    {  return false;  }
    
    const int pairwise_algmt_size = PER_and_hybrid_algmt.first.size();    
    
    type_string__string prepared_PER_and_hybrid_algmt(PER_and_hybrid_algmt);
    
		
                {//error-check
                    const BOOST_Interval valid_interval(0, pairwise_algmt_size);
                    bool contains_valid_indeces = true;
                    
		    for (type_map_uint_to_uint::const_iterator it_map = absolute_coords_hybrid_X_to_pairwise_algmt_profile.begin();
			    it_map != absolute_coords_hybrid_X_to_pairwise_algmt_profile.end();
			    ++it_map) 
		    {
			if ( !BOOST_in(it_map->second, valid_interval) )
			{
			    contains_valid_indeces = false;
			    break;            
			}
		    }
			
		    for (type_map_uint_to_uint::const_iterator it_map = absolute_coords_hybrid_Y_to_pairwise_algmt_profile.begin();
			    it_map != absolute_coords_hybrid_Y_to_pairwise_algmt_profile.end();
			    ++it_map)    
		    {
			if ( !BOOST_in(it_map->second, valid_interval) )
			{
			    contains_valid_indeces = false;
			    break;            
			}            
		    }
                    
                    
                    if (!contains_valid_indeces)
                    {
                        std::stringstream error_strm;                        
                        error_strm << "\n\nERROR!   STAR-ALIGNMENT-ERROR:   contains_valid_indeces = false!!!!\n\n";                        
                        error_strm << "\nPER_name = " << mates_of_PER.first.Name << "\n\n";
			error_strm << "\nLCR_0__LCR_1__XY___YX___XYX___YXY = " << LCR_0__LCR_1__XY___YX___XYX___YXY << "\n";
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_X_to_pairwise_algmt_profile, "absolute_coords_hybrid_X_to_pairwise_algmt_profile",
                                                            &error_strm, true );
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_Y_to_pairwise_algmt_profile, "absolute_coords_hybrid_Y_to_pairwise_algmt_profile",
                                                            &error_strm, true );                                                                  
                        error_strm << "\n\nPER_and_hybrid_algmt:\nPER:\n" << PER_and_hybrid_algmt.first << "\nHybrid:\n" << PER_and_hybrid_algmt.second << "\n\n";
                        error_strm << "\norientation_of_pairwise_algmt_agrees_with_foundation_sequence = " << orientation_of_pairwise_algmt_agrees_with_foundation_sequence
                                    << "\nmust_complement_algmt_to_get_to_foundation_sequence = " << must_complement_algmt_to_get_to_foundation_sequence
                                    << "\n\n";        
				    
			error_message(error_strm, false);  
			
			failed_alignment = true;
                        return false;
                    }    
                }//error-check
    
            
    
                {//error-check
                    bool contains_valid_coords = true;
                    
                    const type_set_int coords(convert_vector_to_set<int>(progressive_star_alignment.at(0)));
		    
		    type_set_uint not_found_coords__X;
		    type_set_uint not_found_coords__Y;
                    
		    for (type_map_uint_to_uint::const_iterator it_map = absolute_coords_hybrid_X_to_pairwise_algmt_profile.begin();
			    it_map != absolute_coords_hybrid_X_to_pairwise_algmt_profile.end();
			    ++it_map)   
		    {
			if (coords.count((int)it_map->first) == 0)
			{
			    not_found_coords__X.insert(it_map->first);
			    contains_valid_coords = false;
			}
		    }
			
		    for (type_map_uint_to_uint::const_iterator it_map = absolute_coords_hybrid_Y_to_pairwise_algmt_profile.begin();
			    it_map != absolute_coords_hybrid_Y_to_pairwise_algmt_profile.end();
			    ++it_map)    
		    {
			if (coords.count((int)it_map->first) == 0)
			{
			    not_found_coords__Y.insert(it_map->first);
			    contains_valid_coords = false;
			}            
		    }
                    
                    
                    if (!contains_valid_coords)
                    {
                        std::stringstream error_strm;
                        
                        error_strm << "\n\nERROR!  STAR-ALIGNMENT-ERROR:  contains_valid_coords = false!!!!\n\n";                        
                        error_strm << "\nPER_name = " << mates_of_PER.first.Name << "\n\n";
			
			error_strm << "\nLCR_0__LCR_1__XY___YX___XYX___YXY = " << LCR_0__LCR_1__XY___YX___XYX___YXY << "\n";
			
			print_set<uint>(not_found_coords__X, "not_found_coords__X", &error_strm, true);
			print_set<uint>(not_found_coords__Y, "not_found_coords__Y", &error_strm, true);
			
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_X_to_pairwise_algmt_profile, "absolute_coords_hybrid_X_to_pairwise_algmt_profile",
                                                            &error_strm, true );
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_Y_to_pairwise_algmt_profile, "absolute_coords_hybrid_Y_to_pairwise_algmt_profile",
                                                            &error_strm, true );   
                                                            
                        print_vector<int>( progressive_star_alignment.at(0), "star_coordinates__ie_first_row", &error_strm, true);
                                                            
                        error_strm << "\n\nPER_and_hybrid_algmt:\nPER:\n" << PER_and_hybrid_algmt.first << "\nHybrid:\n" << PER_and_hybrid_algmt.second << "\n\n";
                        error_strm << "\norientation_of_pairwise_algmt_agrees_with_foundation_sequence = " << orientation_of_pairwise_algmt_agrees_with_foundation_sequence
                                    << "\nmust_complement_algmt_to_get_to_foundation_sequence = " << must_complement_algmt_to_get_to_foundation_sequence
                                    << "\n\n";             
				    
			error_strm <<"\nmates_of_PER.first.Name = " << mates_of_PER.first.Name << "\n\n";
                        
                        error_message(error_strm, false);   
			
			failed_alignment = true;
                        return false;
                    }    
                }//error-check
    
    
    
    
    
    
    
    
    
    type_vector_int absolute_pairwise_algmt_coords__prepared(pairwise_algmt_size,  -1); // thus, a -1 that remains after filling means ins in PER wrt ref...
    
    try 
    {
    
    {
        const type_map_uint_to_uint   map_pairwise_algmt_profile_to_hybrid_X_absolute(get_inverse_of_map(absolute_coords_hybrid_X_to_pairwise_algmt_profile));
        for (type_map_uint_to_uint::const_iterator it_prof_to_abs_A = map_pairwise_algmt_profile_to_hybrid_X_absolute.begin();
                it_prof_to_abs_A != map_pairwise_algmt_profile_to_hybrid_X_absolute.end();
                ++it_prof_to_abs_A)
	{  absolute_pairwise_algmt_coords__prepared.at(it_prof_to_abs_A->first) = (int)it_prof_to_abs_A->second;  }
    }    
    {
        const type_map_uint_to_uint   map_pairwise_algmt_profile_to_hybrid_Y_absolute (  get_inverse_of_map(  absolute_coords_hybrid_Y_to_pairwise_algmt_profile  )   );
        for (type_map_uint_to_uint::const_iterator it_prof_to_abs_B = map_pairwise_algmt_profile_to_hybrid_Y_absolute.begin();
                it_prof_to_abs_B != map_pairwise_algmt_profile_to_hybrid_Y_absolute.end();
                ++it_prof_to_abs_B)
	{  absolute_pairwise_algmt_coords__prepared.at(it_prof_to_abs_B->first) = (int)it_prof_to_abs_B->second;  }
    }    
    
    }//try
    catch (int nnnnnnn)
    {
                        std::stringstream error_strm;
                        error_strm << "\n\nERROR!   STAR-ALIGNMENT-ERROR:  preparing \"absolute_pairwise_algmt_coords__prepared\"   failed!!!!!!\n\n";
                        error_strm << "\nPER_name = " << mates_of_PER.first.Name << "\n\n";
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_X_to_pairwise_algmt_profile, "absolute_coords_hybrid_X_to_pairwise_algmt_profile",
                                                            &error_strm, true );
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_Y_to_pairwise_algmt_profile, "absolute_coords_hybrid_Y_to_pairwise_algmt_profile",
                                                            &error_strm, true );   
                                                            
                        print_vector<int>( progressive_star_alignment.at(0), "star_coordinates__ie_first_row", &error_strm, true);
                        error_strm << "\n\nPER_and_hybrid_algmt:\nPER:\n" << PER_and_hybrid_algmt.first << "\nHybrid:\n" << PER_and_hybrid_algmt.second << "\n\n";
                        error_strm << "\norientation_of_pairwise_algmt_agrees_with_foundation_sequence = " << orientation_of_pairwise_algmt_agrees_with_foundation_sequence
                                    << "\nmust_complement_algmt_to_get_to_foundation_sequence = " << must_complement_algmt_to_get_to_foundation_sequence
                                    << "\n\n";                                        
                        
                        error_message(error_strm, false);
			
			failed_alignment = true;
                        return false;
    }
     




    
    names_of_rows.push_back(mates_of_PER.first.Name);
    progressive_star_alignment.push_back(type_vector_int(star_profile_size, ' '));  //new PER        
    
    
    //prepare:    
    PER_bundle_wrapper__space::remove_N_tails_from_read_and_quality(mates_of_PER.first.QueryBases, mates_of_PER.first.Qualities);
    PER_bundle_wrapper__space::remove_N_tails_from_read_and_quality(mates_of_PER.second.QueryBases, mates_of_PER.second.Qualities);    
    
    for (uint mate=0; mate<2; ++mate)
    {
	for (uint q=0; q < pair_at<BamTools::BamAlignment>(mates_of_PER, mate).Qualities.size(); ++q)
	{  pair_at<BamTools::BamAlignment>(mates_of_PER, mate).Qualities.at(q) -= ASCII_quality_offset;  }
    }
    
    if (must_complement_algmt_to_get_to_foundation_sequence)
    {
        prepared_PER_and_hybrid_algmt.first = get_complement_of_sequence___UNSAFE(prepared_PER_and_hybrid_algmt.first);
        prepared_PER_and_hybrid_algmt.second = get_complement_of_sequence___UNSAFE(prepared_PER_and_hybrid_algmt.second);    
    }
    
    if (!orientation_of_pairwise_algmt_agrees_with_foundation_sequence)
    {
        prepared_PER_and_hybrid_algmt.first = get_reverse_of_sequence(prepared_PER_and_hybrid_algmt.first);
        prepared_PER_and_hybrid_algmt.second = get_reverse_of_sequence(prepared_PER_and_hybrid_algmt.second);
        
        absolute_pairwise_algmt_coords__prepared = get_reverse_of_vector<int>(absolute_pairwise_algmt_coords__prepared);
	
	std::swap<BamTools::BamAlignment>(mates_of_PER.first, mates_of_PER.second);
	for (uint mate=0; mate<2; ++mate)
	{  
	    pair_at<BamTools::BamAlignment>(mates_of_PER, mate).Qualities = get_reverse_of_sequence(pair_at<BamTools::BamAlignment>(mates_of_PER, mate).Qualities);  
	    pair_at<BamTools::BamAlignment>(mates_of_PER, mate).QueryBases = get_reverse_of_sequence(pair_at<BamTools::BamAlignment>(mates_of_PER, mate).QueryBases);//just for consistency.
	}
    }
    
    

    
    
    {//PER_data	
	const type_map_uint_to_uint map_BAM_RefIDs_to_chromo_val( get_inverse_of_map(map_chromosome_value_to_BAM_Ref_IDs)  );		
	
	std::stringstream PER_info;
	PER_info 
	    << "Likelihood: " << P_PER_cond_this_sequence.toString(10) << "\t" 
	    << "Posterior: " << (P_PER_cond_this_sequence / sum_over_all_P_locations).toString(10) << "\t"
	    << "Proper: " << mates_of_PER.first.IsProperPair() << "\t"
	    << "A_Map_Qual: " << mates_of_PER.first.MapQuality << "\t"
	    << "B_Map_Qual: " << mates_of_PER.second.MapQuality << "\t"
	    << "A_Mapped: " << mates_of_PER.first.IsMapped() << "\t"
	    << "B_Mapped: " << mates_of_PER.first.IsMateMapped() << "\t";
	    
	for (uint mate=0; mate<2; ++mate)
	{
	    int original_chromo_val = -1;	    
	    if (map_BAM_RefIDs_to_chromo_val.count(pair_at<BamTools::BamAlignment>(mates_of_PER, mate).RefID) != 0)
	    {  original_chromo_val = (int)map_BAM_RefIDs_to_chromo_val.at(pair_at<BamTools::BamAlignment>(mates_of_PER, mate).RefID);  }
	
	    PER_info << "chr_" << original_chromo_val << "__" <<  pair_at<BamTools::BamAlignment>(mates_of_PER, mate).Position+1 << "\t";
	}
	
	PER_info << "A_Base_Quals:  ";
	for (uint q=0; q < mates_of_PER.first.Qualities.size(); ++q)	
	{
	    PER_info << (int)mates_of_PER.first.Qualities.at(q) << " ";	
	}
	PER_info << "\t";
	
	PER_info << "B_Base_Quals  ";
	for (int q=mates_of_PER.second.Qualities.size()-1; q >= 0; --q)  //remember that mate B is always aligne din reverse directon.
	{  PER_info << (int)mates_of_PER.second.Qualities.at(q) << " ";  }
	
	vector_PER_info.push_back(PER_info.str());
	
	std::string all_Quals(mates_of_PER.first.Qualities);
	all_Quals.append(get_reverse_of_sequence(mates_of_PER.second.Qualities));
	vector_PER_base_Qualities.push_back(all_Quals);	
    }//PER_data
    
    
    
    
    
    const type_vector_vector_int::iterator newest_PER_on_star_profile = --progressive_star_alignment.end();    
    
            {//check degeneracy
                bool contains_real_num = false;
                for (uint j=0; j<absolute_pairwise_algmt_coords__prepared.size(); ++j)
		{
                    if (absolute_pairwise_algmt_coords__prepared.at(j) != -1)
                    {
                        contains_real_num = true;
                        break;
                    }
		}
                    
                if (!contains_real_num)
                {
                    std::stringstream error_strm;                    
                    error_strm << "\n\nERROR!  STAR-ALIGNMENT-ERROR:   contains_real_num = false!!!!\n\n";
                    error_strm << "\nPER_name = " << mates_of_PER.first.Name << "\n\n";
                    print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_X_to_pairwise_algmt_profile, "absolute_coords_hybrid_X_to_pairwise_algmt_profile",
                                                        &error_strm, true );
                    print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_Y_to_pairwise_algmt_profile, "absolute_coords_hybrid_Y_to_pairwise_algmt_profile",
                                                        &error_strm, true );    
                    error_strm << "\n\nPER_and_hybrid_algmt:\nPER:\n" << PER_and_hybrid_algmt.first << "\nHybrid:\n" << PER_and_hybrid_algmt.second << "\n\n";
                    error_strm << "\norientation_of_pairwise_algmt_agrees_with_foundation_sequence = " << orientation_of_pairwise_algmt_agrees_with_foundation_sequence
                                << "\nmust_complement_algmt_to_get_to_foundation_sequence = " << must_complement_algmt_to_get_to_foundation_sequence
                                << "\n\n";                                     
		    error_message(error_strm, false);
		    
		    failed_alignment = true;
                    return false;
                }
            }//check degeneracy
        
    
    int universal_star_index = 0;

    
    //iterate until you get to the  first aligned index:
    int pairwise_algmt_index__of_first_aligned_pos_of_hybrid = 0;
    while (pairwise_algmt_index__of_first_aligned_pos_of_hybrid < (int)pairwise_algmt_size  and 
	    absolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid) == -1)
    {  ++pairwise_algmt_index__of_first_aligned_pos_of_hybrid;  }
    
                    if (pairwise_algmt_index__of_first_aligned_pos_of_hybrid == (int)pairwise_algmt_size)
                    {
                        std::stringstream error_strm;                          
                        error_strm << "\n\nERROR!    STAR-ALIGNMENT-ERROR:  pairwise_algmt_index__of_first_aligned_pos_of_hybrid ==  pairwise_algmt_size!!!\n\n";                        
                        error_strm << "\npairwise_algmt_index__of_first_aligned_pos_of_hybrid = " << pairwise_algmt_index__of_first_aligned_pos_of_hybrid
                                    << "\nabsolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid) = " << absolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid)
                                    << "\nstar_profile_size = " << star_profile_size << "\n\n";
                           
                        print_vector<int>( absolute_pairwise_algmt_coords__prepared, "absolute_pairwise_algmt_coords__prepared", &error_strm, true);                                                   
                        print_vector<int>( progressive_star_alignment.at(0), "star_coordinates__ie_first_row", &error_strm, true);            
                                                                                                            
                        error_strm << "\nPER_name = " << mates_of_PER.first.Name << "\n\n";
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_X_to_pairwise_algmt_profile, "absolute_coords_hybrid_X_to_pairwise_algmt_profile",
                                                            &error_strm, true );
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_Y_to_pairwise_algmt_profile, "absolute_coords_hybrid_Y_to_pairwise_algmt_profile",
                                                            &error_strm, true );                                                                
                        error_strm << "\n\nPER_and_hybrid_algmt:\nPER:\n" << PER_and_hybrid_algmt.first << "\nHybrid:\n" << PER_and_hybrid_algmt.second << "\n\n";
                        error_strm << "\norientation_of_pairwise_algmt_agrees_with_foundation_sequence = " << orientation_of_pairwise_algmt_agrees_with_foundation_sequence
                                    << "\nmust_complement_algmt_to_get_to_foundation_sequence = " << must_complement_algmt_to_get_to_foundation_sequence
                                    << "\n\n";                                        
                        
                        error_message(error_strm, false);
			
			failed_alignment = true;
                        return false;
                    }       
    
    
    //we now have the first aligned index.
    
    absolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid);
    
    
    while(universal_star_index  <  star_profile_size)
    {        
        if (progressive_star_alignment.at(0).at(universal_star_index) 
		!=   absolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid))
	{  ++universal_star_index;  }
        else
	{  break;  }
    }
           
                    if (universal_star_index == star_profile_size)
                    {
                        std::stringstream error_strm;                        
                        error_strm << "\n\nERROR!    STAR-ALIGNMENT-ERROR:  unable to match first coords!!!\n\n";
                        error_strm << "\npairwise_algmt_index__of_first_aligned_pos_of_hybrid = " << pairwise_algmt_index__of_first_aligned_pos_of_hybrid
                                    << "\nabsolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid) = " << absolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid)
                                    << "\nstar_profile_size = " << star_profile_size << "\n\n";
                        print_vector<int>( absolute_pairwise_algmt_coords__prepared, "absolute_pairwise_algmt_coords__prepared", &error_strm, true);         
                        print_vector<int>( progressive_star_alignment.at(0), "star_coordinates__ie_first_row", &error_strm, true);            
                        
                        error_strm << "\nPER_name = " << mates_of_PER.first.Name << "\n\n";
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_X_to_pairwise_algmt_profile, "absolute_coords_hybrid_X_to_pairwise_algmt_profile",
                                                            &error_strm, true );
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_Y_to_pairwise_algmt_profile, "absolute_coords_hybrid_Y_to_pairwise_algmt_profile",
                                                            &error_strm, true );   
                  
                        error_strm << "\n\nPER_and_hybrid_algmt:\nPER:\n" << PER_and_hybrid_algmt.first << "\nHybrid:\n" << PER_and_hybrid_algmt.second << "\n\n";
                        error_strm << "\norientation_of_pairwise_algmt_agrees_with_foundation_sequence = " << orientation_of_pairwise_algmt_agrees_with_foundation_sequence
                                    << "\nmust_complement_algmt_to_get_to_foundation_sequence = " << must_complement_algmt_to_get_to_foundation_sequence
                                    << "\n\n";                                        
                        
                        error_message(error_strm,false);
			
			failed_alignment = true;
                        return false;
                    }                                    
    
    
    
    
    
    
    
    
    
    if( pairwise_algmt_index__of_first_aligned_pos_of_hybrid > 0)
    {
        //then there were a bunch of insertions in the PER relative to the ref.  so we should insert them before this universal position in the star profile.
        
        int number_of_gaps_in_star_profile_before_beginning_of_pairwise_algmt_on_star_algmt = 0;
                        {//get number_of_gaps_in_star_profile_before_beginning_of_pairwise_algmt_on_star_algmt
                                    int prev_gap_ctr = universal_star_index - 1;
                                    while (   prev_gap_ctr >= 0    and    progressive_star_alignment.at(0).at(prev_gap_ctr) == -1   )
                                    {
                                        ++number_of_gaps_in_star_profile_before_beginning_of_pairwise_algmt_on_star_algmt;
                                        --prev_gap_ctr;        
                                    }
                        }//get number_of_gaps_in_star_profile_before_beginning_of_pairwise_algmt_on_star_algmt
        
        
        
        if (number_of_gaps_in_star_profile_before_beginning_of_pairwise_algmt_on_star_algmt == 0)
        {            
            expand_star_profile_for_inserted_bases_in_a_PER(  universal_star_index,
                                                              newest_PER_on_star_profile,
                                                              prepared_PER_and_hybrid_algmt.first.substr(0, pairwise_algmt_index__of_first_aligned_pos_of_hybrid)  );                                                 
        } // make new gaps
        else  // there are some existing gaps
        {
            int next_index_in_pairwise_algmt_to_be_inserted = 0;
            int next_star_gap_to_be_filled_by_PER_insertion = universal_star_index - number_of_gaps_in_star_profile_before_beginning_of_pairwise_algmt_on_star_algmt;                     
            
            while (        next_index_in_pairwise_algmt_to_be_inserted < pairwise_algmt_index__of_first_aligned_pos_of_hybrid  
                     and   next_star_gap_to_be_filled_by_PER_insertion < universal_star_index   )
            {
                newest_PER_on_star_profile->at(next_star_gap_to_be_filled_by_PER_insertion)
                                                =   (int)prepared_PER_and_hybrid_algmt.first.at(next_index_in_pairwise_algmt_to_be_inserted);
                                                
                ++next_index_in_pairwise_algmt_to_be_inserted;   
                ++next_star_gap_to_be_filled_by_PER_insertion;            
            }                        
            
            
            if ( next_index_in_pairwise_algmt_to_be_inserted < pairwise_algmt_index__of_first_aligned_pos_of_hybrid)//leftover !!!  more ins than existing gaps!
            {
                expand_star_profile_for_inserted_bases_in_a_PER
                        ( universal_star_index,
                        newest_PER_on_star_profile,
                        prepared_PER_and_hybrid_algmt.first.substr( next_index_in_pairwise_algmt_to_be_inserted,
                                                              pairwise_algmt_index__of_first_aligned_pos_of_hybrid - next_index_in_pairwise_algmt_to_be_inserted)  );
            }
                
        } // fill in some gaps
    
    } //ins in beginning of PER
    
    
    
    
    
    
    
    
    
    int pairwise_algmt_ctr = pairwise_algmt_index__of_first_aligned_pos_of_hybrid;
    
    while (pairwise_algmt_ctr < pairwise_algmt_size)
    {               
        //ORDER OF IF-STATEMENTS MATTERS!!!!                  
        if (  absolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_ctr) == progressive_star_alignment.at(0).at(universal_star_index)  ) //as real coords, or as -1
        {
            newest_PER_on_star_profile->at(universal_star_index) = (int)prepared_PER_and_hybrid_algmt.first.at(pairwise_algmt_ctr);
            ++pairwise_algmt_ctr;
            ++universal_star_index;            
        }
        else if (  progressive_star_alignment.at(0).at(universal_star_index) == -1)  //there is a gap in the star alignment, but none in this PER
        {
            while (progressive_star_alignment.at(0).at(universal_star_index) == -1)
            {
                newest_PER_on_star_profile->at(universal_star_index) = (int)'-';
                ++universal_star_index;
            }        
        }
        else    //ins in PER but not yet in star algmt.  expand ENTIRE star alignment.
        {
            expand_star_profile_for_inserted_bases_in_a_PER( universal_star_index,
                                                             newest_PER_on_star_profile,
                                                             prepared_PER_and_hybrid_algmt.first.at(pairwise_algmt_ctr) );            
            ++pairwise_algmt_ctr;           
        }  
        
        
        
                    if (universal_star_index >= star_profile_size)
                    {
                        std::stringstream error_strm;
                        error_strm << "\n\nERROR!    STAR-ALIGNMENT-ERROR:  universal_star_index = " << universal_star_index
                                    << "    >=    " <<  star_profile_size << " = star_profile_size   !!!\n\n";
                                    
                        error_strm << "\n\npairwise_algmt_ctr = " << pairwise_algmt_ctr
                                    << "\npairwise_algmt_size = " << pairwise_algmt_size << "\n\n";
                        
                        error_strm << "\npairwise_algmt_index__of_first_aligned_pos_of_hybrid = " << pairwise_algmt_index__of_first_aligned_pos_of_hybrid
                                    << "\nabsolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid) = " << absolute_pairwise_algmt_coords__prepared.at(pairwise_algmt_index__of_first_aligned_pos_of_hybrid)
                                    << "\nstar_profile_size = " << star_profile_size << "\n\n";
                           
                        print_vector<int>( absolute_pairwise_algmt_coords__prepared, "absolute_pairwise_algmt_coords__prepared", &error_strm, true);               
                                    
                        print_vector<int>( progressive_star_alignment.at(0), "star_coordinates__ie_first_row", &error_strm, true);      
                        error_strm << "\nPER_name = " << mates_of_PER.first.Name << "\n\n";
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_X_to_pairwise_algmt_profile, "absolute_coords_hybrid_X_to_pairwise_algmt_profile",
                                                            &error_strm, true );
                        print_map_keys_and_values<uint,uint>(absolute_coords_hybrid_Y_to_pairwise_algmt_profile, "absolute_coords_hybrid_Y_to_pairwise_algmt_profile",
                                                            &error_strm, true );  
                                                            
                        error_strm << "\n\nPER_and_hybrid_algmt:\nPER:\n" << PER_and_hybrid_algmt.first << "\nHybrid:\n" << PER_and_hybrid_algmt.second << "\n\n";
                        error_strm << "\norientation_of_pairwise_algmt_agrees_with_foundation_sequence = " << orientation_of_pairwise_algmt_agrees_with_foundation_sequence
                                    << "\nmust_complement_algmt_to_get_to_foundation_sequence = " << must_complement_algmt_to_get_to_foundation_sequence
                                    << "\n\n";                                        
                        
                        error_message(error_strm, false);
			
			failed_alignment = true;
                        return false;
                    }               
        
        
    }//while pairwise
    
    
    return true;
    
}  //    add_pairwise_alignment_to_progressive_star_alignment






























void Star_alignment::expand_star_profile_for_inserted_bases_in_a_PER
                                                (int &universal_star_index,
                                                 const type_vector_vector_int::iterator &newest_PER_on_star_profile,
                                                 const char &some_char)
{   
    expand_star_profile_for_inserted_bases_in_a_PER
                                        (universal_star_index,
                                         newest_PER_on_star_profile,
                                         std::string(1, some_char));
                                         
}// expand_star_profile_for_inserted_bases_in_a_PER         




void Star_alignment::expand_star_profile_for_inserted_bases_in_a_PER
                                                (int &universal_star_index,
                                                 const type_vector_vector_int::iterator &newest_PER_on_star_profile,
                                                 const std::string &insertion_string   )
{
    if (failed_alignment)
    {  return;  }
    
    
    const uint length_of_insertion = insertion_string.size();
    const type_vector_int insertion_vector(convert_string_to_ints(insertion_string));
    
    
    progressive_star_alignment.at(0).insert( progressive_star_alignment.at(0).begin() + universal_star_index,
                                            length_of_insertion,
                                            -1  );

    progressive_star_alignment.at(1).insert(   progressive_star_alignment.at(1).begin() + universal_star_index,
                                                length_of_insertion,
                                                (int)'-'  );                                            
                                                                                        
    for (uint seq_in_star = 2; seq_in_star < progressive_star_alignment.size()-1; ++seq_in_star)    //reads
    {
        if (universal_star_index >= star_profile_size)
            progressive_star_alignment.at(seq_in_star).insert(   progressive_star_alignment.at(seq_in_star).begin() + universal_star_index,
                                                                length_of_insertion,
                                                                (int)'-'  );            
        else if ( progressive_star_alignment.at(seq_in_star).at(universal_star_index) ==  (int)' ' )
            progressive_star_alignment.at(seq_in_star).insert(   progressive_star_alignment.at(seq_in_star).begin() + universal_star_index,
                                                                length_of_insertion,
                                                                (int)' '  );
        else if (universal_star_index > 0   and    progressive_star_alignment.at(seq_in_star).at(universal_star_index - 1) == (int)' ' )
            progressive_star_alignment.at(seq_in_star).insert(   progressive_star_alignment.at(seq_in_star).begin() + universal_star_index,
                                                                length_of_insertion,
                                                                (int)' '  );
        else // == 0 or previous base was not a space
            progressive_star_alignment.at(seq_in_star).insert(   progressive_star_alignment.at(seq_in_star).begin() + universal_star_index,
                                                                length_of_insertion,
                                                                (int)'-'  );       
    }
        
        
        

    newest_PER_on_star_profile->insert(  newest_PER_on_star_profile->begin() + universal_star_index,
                                         insertion_vector.begin(),
                                         insertion_vector.end()    );
    
    
    if (universal_star_index  <=   star_profile_index_of_first_Y_nucleotide_in_hybrid_XY)                                        
    {  star_profile_index_of_first_Y_nucleotide_in_hybrid_XY +=  (int)length_of_insertion;  }
    
    if (universal_star_index  <=  star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX)
    {  star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX += (int)length_of_insertion;  }
                                        
    universal_star_index += (int)length_of_insertion;
    star_profile_size +=  (int)length_of_insertion;
    
    
    
    for (uint r=1; r<progressive_star_alignment.size(); ++r)
    {
	if (progressive_star_alignment.at(r).size() != progressive_star_alignment.at(0).size())
	{
	    std::stringstream debug_strm;
	    debug_strm << "progressive_star_alignment.at(" << r << ").size() = " << progressive_star_alignment.at(r).size() 
			<< "  !=  " <<  progressive_star_alignment.at(0).size() << " = progressive_star_alignment.at(0).size() "
			<< "\nthat was in \"expand_star_profile_for_inserted_bases_in_a_PER\", "
			<< "\nuniversal_star_index = " << universal_star_index
			<< "\ninsertion_string = " << insertion_string 
			<< "\n\n";
	    std::cerr << "\n\n" << debug_strm.str() << "\n\n";
	}
    }
    
        
}  //  expand_star_profile_for_inserted_bases_in_a_PER






































bool Star_alignment::create_HTML_output_for_star_alignment
                                                (const int &chromo_of_region_X,
                                                 const int &chromo_of_region_Y,
                                                 const bool &crop_to_relevant_region,
                                                 type_vector_string &star_html_output)
{
    assert(activated_star_alignment);
    
    if (failed_alignment)
    {
	star_html_output.clear();
	star_html_output.push_back("<span style=\"font-size: 36pt\">ERROR IN STAR ALIGNMENT</span>");
	return false;    
    }    
    
    int prepared_crop_amount_below;
    int prepared_crop_amount_above;
  
    const type_vector_vector_int prepared___progressive_star_alignment(
						sort_by_first_alignment_position_in_profile__and__crop(crop_to_relevant_region,
												       prepared_crop_amount_below,
													prepared_crop_amount_above));
    
    const int prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY 
                    =   (star_profile_index_of_first_Y_nucleotide_in_hybrid_XY  >=   star_profile_size - prepared_crop_amount_above) 
							?
                                    prepared___progressive_star_alignment.at(0).size()       
                                    :        star_profile_index_of_first_Y_nucleotide_in_hybrid_XY - prepared_crop_amount_below;
                                
    const int prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX
                    =   (star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX  >=   star_profile_size - prepared_crop_amount_above) 
							?
                                    prepared___progressive_star_alignment.at(0).size()       
                                    :        star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX - prepared_crop_amount_below;    
    
//     std::cerr	
// 	    <<" \n\n\tstar_profile_size = " << star_profile_size
// 	    << "\n\tstar_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
// 	    << "\n\tstar_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = " << star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX
// 	    << "\n\tprepared_crop_amount_above = " << prepared_crop_amount_above
// 	    << "\n\tprepared_crop_amount_below = " << prepared_crop_amount_below
// 	    << "\n\tprepared____star_profile_size = " << prepared___progressive_star_alignment.at(0).size()
// 	    << "\n\tprepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
// 	    << "\n\tprepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = " << prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX
// 	    << "\n\n";
		
    
				    
    
    
    star_html_output = type_vector_string( progressive_star_alignment.size(), std::string());
        
    
    
    std::string begin_ref_name;
    std::string end_ref_name;  
    
    uint longest_name__padded;
    
    {//names
                    int first_coord_val = -1;                                        
                    for (int j=0; j< prepared___progressive_star_alignment.at(0).size(); ++j)
		    {
                        if (prepared___progressive_star_alignment.at(0).at(j) != -1)
                        {
                            first_coord_val = prepared___progressive_star_alignment.at(0).at(j);
                            break;                        
                        }
		    }
                      
                    int last_coord_val = -1;  
                    for (int j=prepared___progressive_star_alignment.at(0).size()-1; j >= 0; --j)
		    {
                        if (prepared___progressive_star_alignment.at(0).at(j) != -1)
                        {
                            last_coord_val = prepared___progressive_star_alignment.at(0).at(j);
                            break;                        
                        } 
		    }
                        
                    
    
                                                                              
                    {                        
                        char chr_and_coords[option_size];
                        std::sprintf(chr_and_coords, "chr %d:  %d", chromo_of_region_X, first_coord_val);
                        begin_ref_name.assign(chr_and_coords);
                    
                        
                        if (LCR_0__LCR_1__XY___YX___XYX___YXY <= 1)
                            std::sprintf(chr_and_coords, "chr %d:  %d", chromo_of_region_X, last_coord_val);
                        else
                            std::sprintf(chr_and_coords, "chr %d:  %d", chromo_of_region_Y, last_coord_val);
                        
                        end_ref_name.assign(chr_and_coords);
                    }
    
    
    
                    
                    longest_name__padded = begin_ref_name.size();
                    
                    for (uint j=2; j < progressive_star_alignment.size(); ++j)    
		    {  longest_name__padded  = std::max<uint>(longest_name__padded, names_of_rows.at(j).size()  );  }
                    
                    longest_name__padded += 5;  // for spacing
                    
                    
                    star_html_output.at(0).append( std::string( longest_name__padded,  ' ') );
                    star_html_output.at(1).append( begin_ref_name );
                    star_html_output.at(1).append(  std::string( longest_name__padded - begin_ref_name.size(),  ' ')  );
                    
                    for (uint row = 2; row < progressive_star_alignment.size(); ++row)
                    {
                        star_html_output.at(row).append(  names_of_rows.at(row)   );
                        star_html_output.at(row).append(  std::string( longest_name__padded - names_of_rows.at(row).size(),  ' ')   );
                    }
    }//names   
        
        
    
    if (   (chromo_of_region_X < 0   and
	    (star_profile_index_of_first_Y_nucleotide_in_hybrid_XY > 0   or  star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX < star_profile_size) )
        or (chromo_of_region_Y < 0   and 
	    (star_profile_index_of_first_Y_nucleotide_in_hybrid_XY < star_profile_size  and  star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX >= 0) ))
    {
        std::stringstream error_strm;
        
        error_strm << "\n\nERROR!     STAR-ALIGNMENT-ERROR:   in   \"create_HTML_output_for_star_alignment\"  for Star Alignment!!!\n\n\n"
                   << "chromo_of_region_X = " << chromo_of_region_X
                   << "\nstar_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
                   << "\nchromo_of_region_Y = " << chromo_of_region_Y
                   << "\nstar_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
                   << "\nstar_profile_size = " << star_profile_size
		    << "\n\n\nprogressive_star_alignment.size() = " << progressive_star_alignment.size()
                    << "\nstar_profile_size = " << star_profile_size
                    << "\nstar_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
                    << "\nstar_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = " << star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX;
		    
        error_strm << "\n\n\n\n\nprogressive alignment:\n\n\n";
        for (uint j=0; j < star_profile_size; ++j)           
	{  error_strm << names_of_rows.at(j) << "\t\t\t" <<  convert_ints_to_string(   progressive_star_alignment.at(j)   ) << "\n\n";  }           
                
        error_message(error_strm, false);     
	
	failed_alignment = true;
	{
	    star_html_output.clear();
	    star_html_output.push_back("<span style=\"font-size: 36pt\">ERROR IN STAR ALIGNMENT</span>");  
	}	
	
        return  false;
    }
            
            
    
    
    
    
    
    
    type_vector_bool variational_positions_of_star_profile(prepared___progressive_star_alignment.at(0).size(), false);        
    {//varpos        
        if (prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY > (int)prepared___progressive_star_alignment.at(0).size())  //sanity check
        {
            std::stringstream error_strm;
            error_strm << "\n\nERROR!     STAR-ALIGNMENT-ERROR:    in  \"create_HTML_output_for_star_alignment\"...\n\n";            
            error_strm << "star_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
                        << "\nstar_profile_size = " << star_profile_size 
                        << "\nprepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY = " << prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY
                        << "\nprepared____star_profile_size = " << prepared___progressive_star_alignment.at(0).size() << "\n\n\n\n";"\n\n\n\n";                                                
            print_this_Star_alignment(&error_strm);            
            error_message(error_strm, false);      
	    
	    failed_alignment = true;
	    {
		star_html_output.clear();
		star_html_output.push_back("<span style=\"font-size: 36pt\">ERROR IN STAR ALIGNMENT</span>");  
	    }		    
            return  false;
        }//sanity check            
        
        
        
    
        if (prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY > 0
	    or   prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX < prepared___progressive_star_alignment.at(0).size())  //chromo A is valid
        {
            const type_map_uint_to_set_uint::const_iterator universal_variational_positions_on_chromo_of_X_it = universal_variational_positions.find( (uint)chromo_of_region_X );
            const type_map_uint_to_set_uint::const_iterator universal_del_positions_on_chromo_of_X_it = universal_del_positions.find( (uint)chromo_of_region_X);        
            const type_map_uint_to_set_uint::const_iterator universal_ins_positions_on_chromo_of_X_it = universal_ins_positions.find( (uint)chromo_of_region_X);              
                
	    //until first breakpoint
            for (int pos = 0; pos < std::min<int>(prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY, prepared___progressive_star_alignment.at(0).size()); ++pos)
	    {
                if ( prepared___progressive_star_alignment.at(0).at(pos) != -1         // a real coordinate
                     and  
                     (     universal_variational_positions_on_chromo_of_X_it->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) ) > 0
                        or universal_del_positions_on_chromo_of_X_it->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) )  >  0
                        or universal_ins_positions_on_chromo_of_X_it->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) )  >  0   )     )  
		{  variational_positions_of_star_profile.at(pos) = true;  }
	    }
	    
            for (int pos = std::max<int>(0,prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX);
		    pos < prepared___progressive_star_alignment.at(0).size();
		    ++pos)
	    {
                if ( prepared___progressive_star_alignment.at(0).at(pos) != -1         // a real coordinate
                     and  
                     (     universal_variational_positions_on_chromo_of_X_it->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) ) > 0
                        or universal_del_positions_on_chromo_of_X_it->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) )  >  0
                        or universal_ins_positions_on_chromo_of_X_it->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) )  >  0   )     )  
		{  variational_positions_of_star_profile.at(pos) = true;  }
	    }	    	    
        }//A
        
        
        if (  prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY < prepared___progressive_star_alignment.at(0).size()
		and  prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX  >= 0)  //chromo B is valid
        {
            
            const type_map_uint_to_set_uint::const_iterator universal_variational_positions_on_chromo_of_Y_ptr = universal_variational_positions.find( (uint)chromo_of_region_Y);  
            const type_map_uint_to_set_uint::const_iterator universal_del_positions_on_chromo_of_Y_ptr = universal_del_positions.find( (uint)chromo_of_region_Y);        
            const type_map_uint_to_set_uint::const_iterator universal_ins_positions_on_chromo_of_Y_ptr = universal_ins_positions.find( (uint)chromo_of_region_Y);              
                
            for (int pos = std::max<int>(0, prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY);
		    pos < std::min<int>(prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX, prepared___progressive_star_alignment.at(0).size());
		    ++pos)
	    {
                if ( prepared___progressive_star_alignment.at(0).at(pos) != -1         // a real coordinate
                     and  
                     (     universal_variational_positions_on_chromo_of_Y_ptr->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) ) > 0
                        or universal_del_positions_on_chromo_of_Y_ptr->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) )  >  0
                        or universal_ins_positions_on_chromo_of_Y_ptr->second.count( (uint)prepared___progressive_star_alignment.at(0).at(pos) )  >  0   )     )  
		{  variational_positions_of_star_profile.at(pos) = true;  }
	    }	    
        }//B
    }//varpos
    
    
    
    
    std::string color_for_X;
    std::string color_for_Y;    
    
    switch (LCR_0__LCR_1__XY___YX___XYX___YXY)
    {
        case 0:        //  0
            color_for_X.assign("Red");
            break;
        case 1:        //  1
            color_for_X.assign("Blue");
            break;
        case 2:        //  AB
            color_for_X.assign("Red");
            color_for_Y.assign("Blue");            
            break;
        case 3:       //  BA
            color_for_X.assign("Blue");
            color_for_Y.assign("Red");            
            break;  
	case 4:  //ABA
            color_for_X.assign("Red");
            color_for_Y.assign("Blue");            
            break;
	case 5:  //BAB
            color_for_X.assign("Blue");
            color_for_Y.assign("Red");            
            break;  	    	    
        default:
            std::stringstream error_strm;
            error_strm << "\n\nERROR!    STAR-ALIGNMENT-ERROR:   unrecognized value    LCR_0__LCR_1__XY___YX___XYX___YXY ="
		    << LCR_0__LCR_1__XY___YX___XYX___YXY
		    << " in  \"Star_alignment::create_HTML_output_for_star_alignment\"\n\n";                        
            print_this_Star_alignment( &error_strm );            
            error_message(error_strm, false);
	    {
		star_html_output.clear();
		star_html_output.push_back("<span style=\"font-size: 36pt\">ERROR IN STAR ALIGNMENT</span>");  
	    }
	    
            return false;
            break;
    };
    
    
    
    const type_vector_vector_int::const_iterator ref_algmt_in_star =  ++prepared___progressive_star_alignment.begin();
    
    for (int universal_star_pos = 0; universal_star_pos < prepared___progressive_star_alignment.at(0).size(); ++universal_star_pos)
    {
        if (universal_star_pos == prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY)
	{  star_html_output.at(0).append( "Breakpoint!");  }
        else if (universal_star_pos == prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX)
	{  star_html_output.at(0).append( "Breakpoint!");  }
        else if (   (universal_star_pos < prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY  
		      or universal_star_pos > prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY + 10)
		and (universal_star_pos < prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX  
		      or universal_star_pos > prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX + 10) )//10 = length of word "Breakpoint!"  after the 'B'
	{  star_html_output.at(0).push_back(' ');  }
        
        if (universal_star_pos >= ref_algmt_in_star->size())
	{
			
	    for (uint r=0; r<progressive_star_alignment.size(); ++r)
	    {
		if (progressive_star_alignment.at(r).size() != progressive_star_alignment.at(0).size())
		{
		    std::stringstream debug_strm;
		    debug_strm << "prepared___progressive_star_alignment.at(" << r << ").size() = " << prepared___progressive_star_alignment.at(r).size();
		    std::cerr << "\n\n" << debug_strm.str() << "\n\n";
		}
	    }
	}
	
        
        if (ref_algmt_in_star->at(universal_star_pos)  ==  '-')
        {//gap in ref                         
            star_html_output.at(1).push_back( '-' );
            
            for (uint row = 2; row < progressive_star_alignment.size(); ++row)
            {
                if ( prepared___progressive_star_alignment.at(row).at(universal_star_pos) !=  '-' 
                     and  prepared___progressive_star_alignment.at(row).at(universal_star_pos)  != ' ')   // 'x' should be impossible.
                {
                    star_html_output.at(row).append( "<span style=\"background-color: " );
                    star_html_output.at(row).append( "Gold"  );
                    star_html_output.at(row).append("\">");            
                    star_html_output.at(row).push_back(prepared___progressive_star_alignment.at(row).at(universal_star_pos));            
                    star_html_output.at(row).append("</span>");                                      
                }
                else  //should be   '-'   or   ' '
                {
                    star_html_output.at(row).push_back(prepared___progressive_star_alignment.at(row).at(universal_star_pos));
                }
            }//rows
        }//gap in ref                        
        else  //nucleotide in ref
        {
            std::string ref_color;
            
	    //color reference/foundation sequence
            if (variational_positions_of_star_profile.at(universal_star_pos))  //is a varpos
            {                
                if (LCR_0__LCR_1__XY___YX___XYX___YXY <= 1) //no event
		{  ref_color  = color_for_X;  }
                else if (LCR_0__LCR_1__XY___YX___XYX___YXY <= 3) // NAHR
		{
		    ref_color =  universal_star_pos  <  prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY   ?
						    color_for_X   :   color_for_Y;
		}
		else
		{//GeneConv
		    if (universal_star_pos  <  prepared____star_profile_index_of_first_Y_nucleotide_in_hybrid_XY)
			ref_color = color_for_X;
		    else if (universal_star_pos  <  prepared____star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX)
			ref_color = color_for_Y;
		    else
			ref_color = color_for_X;
		}//GeneConv
		                         
                star_html_output.at(1).append("<span style=\"background-color: ");
                star_html_output.at(1).append(ref_color);
                star_html_output.at(1).append(";color: White\">");            
                star_html_output.at(1).push_back(ref_algmt_in_star->at(universal_star_pos));            
                star_html_output.at(1).append("</span>");  
            }
            else
	    {  star_html_output.at(1).push_back(ref_algmt_in_star->at(universal_star_pos));  }
            
            //color reads
            for (uint row = 2; row < progressive_star_alignment.size(); ++row)
            {
                if (       prepared___progressive_star_alignment.at(row).at(universal_star_pos)  ==  'x'
                      or   prepared___progressive_star_alignment.at(row).at(universal_star_pos)  ==  ' ')  //not important
                {
                    star_html_output.at(row).push_back(prepared___progressive_star_alignment.at(row).at(universal_star_pos));    
                }
                else   //ref has a vlaue, and PER has a value (maybe '-')
                {                 
                    std::string  PER_nucleotide_color;
		    if (prepared___progressive_star_alignment.at(row).at(universal_star_pos) == ref_algmt_in_star->at(universal_star_pos))
		    {  PER_nucleotide_color = ref_color;  }
		    else
		    {
			const int qual_index =  get_total_PER_mate_index_at_star_index(
								prepared___progressive_star_alignment,
								   (int)row,
								   universal_star_pos);
								   
			if (vector_PER_base_Qualities.at(row).at(qual_index) < 10)
			{  PER_nucleotide_color = "Yellow";  }
			else
			{  PER_nucleotide_color = "Gold";  }
		    }
                                                         
                    if (!PER_nucleotide_color.empty())      // is empty if this is a non-varpos of ref
                    {
                        star_html_output.at(row).append( "<span style=\"background-color: " );
                        star_html_output.at(row).append( PER_nucleotide_color  );
                        if (  !(PER_nucleotide_color.compare("Gold") == 0  or PER_nucleotide_color.compare("Yellow") == 0)  )
			{  star_html_output.at(row).append(";color: White");  }
                        star_html_output.at(row).append("\">");            
                        star_html_output.at(row).push_back(  prepared___progressive_star_alignment.at(row).at(universal_star_pos)  );            
                        star_html_output.at(row).append("</span>");
                    }
                    else                                        
		    {  star_html_output.at(row).push_back(prepared___progressive_star_alignment.at(row).at(universal_star_pos));  }
                }                                
            }//rows            
        } //nucleotide in ref                         
    }// universal_star_pos
    
    
    
    star_html_output.at(1).append(  std::string(5, ' ')  );
    star_html_output.at(1).append(  end_ref_name  );
    
    for (uint row=2; row < progressive_star_alignment.size(); ++row)
    {
        star_html_output.at(row).append(  std::string(5, ' ')  );  
	star_html_output.at(row).append(  vector_PER_info.at(row)   );
    }
    
    
    
    
    return  true;        

}//  create_HTML_output_for_star_alignment

































// void Star_alignment::save_star_alignments_to_file_as_html
//                                 (const uint &Conn_Comp_ID,
//                                  const uint &EV_UID,
//                                  const uint &desired_breakpoint,
//                                  const int &LCR_0__LCR_1__XY___YX___XYX___YXY,   // 0...3
//                                  const int &chromo_MR_A,
//                                  const int &chromo_MR_B,
//                                  const bool &crop_to_relevant_region)       const
// {
// 
//     std::string suffix_type;
//     switch (LCR_0__LCR_1__XY___YX___XYX___YXY)
//     {
//         case 0: 
//             suffix_type = "LCR_A";
//             break;
//         case 1:
//             suffix_type = "LCR_B";
//             break;            
//         case 2:
//             suffix_type = "Hybrid_AB";
//             break;      
//         case 3:
//             suffix_type = "Hybrid_BA";
//             break;  
//         default:
//             std::stringstream error_strm;
//             print_line_of_markers("ERROR! ", &error_strm);
//             print_line_of_markers("(", &error_strm);
//             error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
//             
//             error_strm << "\n\nERROR!    STAR-ALIGNMENT-ERROR:    undefined  value for  LCR_0__LCR_1__XY___YX___XYX___YXY = " << LCR_0__LCR_1__XY___YX___XYX___YXY 
//                         <<  "\n in \"save_star_alignments_to_file_as_html\"\n\n";
//             
//             print_line_of_markers(")", &error_strm);
//             std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
//             return;
//             break;    
//     };
//     
//     
//     
//     
// 
//     
//     
//     char star_algmts_filename[option_size];        
//     std::sprintf(star_algmts_filename, "%s/star_alignments___%s__cc_%u__ev_%u___breakpoint_%u___%s.html",
//                                     output_dir.c_str(),
// 				    genome_name.c_str(),
//                                     Conn_Comp_ID, EV_UID,
//                                     desired_breakpoint,
//                                     suffix_type.c_str()  );       
//                         
//     std::ofstream  save_star_algmts_filestream( star_algmts_filename );        
//                                 if (  !save_star_algmts_filestream.is_open()  )
//                                 {
//                                     std::fprintf(stderr, "\n\n\n\n\nERROR: Command [%s] could not be opened for \"save_star_alignments_to_file_as_html\".\n\n\n\n\n\n", 
//                                                         star_algmts_filename);
//                                     return;
//                                 }        
//     
//     
//     
//     save_star_algmts_filestream << "\n<pre>\n<font size=1>\n\n";
//           
//     save_star_algmts_filestream << "\n<span style=\"font-size: 24pt\">Genome:  " <<  genome_name << "    Connected component  " <<  Conn_Comp_ID
//                                 << ",     Event  " << EV_UID << ",  breakpoint  " << desired_breakpoint
//                                 << "\n\n" << suffix_type
//                                 << "</span>\n\n\n";
//     
//         
//     type_vector_string  html_output;
//     const bool success_create_html =  create_HTML_output_for_star_alignment(   
//                                                                 chromo_MR_A, chromo_MR_B,
//                                                                 LCR_0__LCR_1__XY___YX___XYX___YXY,
//                                                                 crop_to_relevant_region,
//                                                                 html_output    );
//                                                                                                                                 
//     
//     for (uint row=0; row < html_output.size(); ++row)
//         save_star_algmts_filestream << html_output.at(row) << "\n";
//                                 
//     
//     
//     save_star_algmts_filestream << "\n\n</font>\n</pre>\n";
//     save_star_algmts_filestream.close();
//     
// 
//     
//     
//     if (!success_create_html)
//         boost::filesystem::remove(    boost::filesystem::path( star_algmts_filename )    );
//     
//     
// } // save_star_alignments_to_file_as_html



















void Star_alignment::clear_all_data()
{
    
    progressive_star_alignment.clear();
    
    star_profile_size = 0;
    
    names_of_rows.clear();
    
    star_profile_index_of_first_Y_nucleotide_in_hybrid_XY = -1;
    star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX = -1;
    
    vector_PER_info.clear();
    vector_PER_base_Qualities.clear();
    
    mpfr_free_cache();    

}//clear_all_data


        





int Star_alignment::get_total_PER_mate_index_at_star_index
			    (const type_vector_vector_int  &prepared___progressive_star_alignment,
			     const int &row_index_for_PER_in_prepared__progressive_alignment,
			     const int &base_index_in_prepared__progressive_alignment)
{
    
    int number_of_actual_bases_before_desired_position = 0;
    for (uint indx = 0; indx < base_index_in_prepared__progressive_alignment; ++indx)
    {
	if (allowable_nucleotide_values.count(prepared___progressive_star_alignment.at(row_index_for_PER_in_prepared__progressive_alignment).at(indx)) > 0)
	{ ++number_of_actual_bases_before_desired_position;  }	    
    }    
    
    return number_of_actual_bases_before_desired_position;    

}//get_total_PER_mate_index_at_star_index











type_vector_vector_int Star_alignment::sort_by_first_alignment_position_in_profile__and__crop
							(const bool &crop_to_relevant_region,
							 int &prepared_crop_amount_below,
							int &prepared_crop_amount_above)
{
    type_vector_vector_int  prepared___progressive_star_alignment(progressive_star_alignment);
    prepared_crop_amount_below = 0;
    prepared_crop_amount_above = 0;
    
    if (progressive_star_alignment.size() > 2)
    {
	int earliest_start_pos = star_profile_size;
	int latest_end_pos = 0; 
	
	
	type_multimap_uint_to_uint  star_firstpos_to_oldrow;
	
	//determine extemes across al PER alignments:
	for (uint row = 2; row < progressive_star_alignment.size(); ++row)
	{
	    const std::string star_algmt_as_string(convert_ints_to_string(progressive_star_alignment.at(row)));            
	    {
		const int indx_of_first_non_space = star_algmt_as_string.find_first_not_of(' ');            
		if (indx_of_first_non_space != std::string::npos  
		    and  indx_of_first_non_space < star_profile_size  and  indx_of_first_non_space < earliest_start_pos)
		{  earliest_start_pos = indx_of_first_non_space;  }
		
		star_firstpos_to_oldrow.insert(type_uint__uint((uint)indx_of_first_non_space, row));		
	    }
	    {
		const int indx_of_last_non_space =  star_algmt_as_string.find_last_not_of(' ');                        
		if (indx_of_last_non_space != std::string::npos
		    and  indx_of_last_non_space < star_profile_size  and  indx_of_last_non_space > latest_end_pos)
		{    latest_end_pos = indx_of_last_non_space;  }
	    }
	}
	
	
	{//sort
	    type_vector_vector_int  sorted_progressive_algmt(progressive_star_alignment);
	    type_vector_string sorted_names(names_of_rows);
	    type_vector_string  sorted_info(vector_PER_info);
	    type_vector_string  sorted_quals(vector_PER_base_Qualities);
	    
	    uint row_ctr = 2;
	    for (type_multimap_uint_to_uint::const_iterator it_sort = star_firstpos_to_oldrow.begin();
		    it_sort != star_firstpos_to_oldrow.end();
		    ++it_sort)
	    {
		sorted_progressive_algmt.at(row_ctr) = progressive_star_alignment.at(it_sort->second);
		sorted_names.at(row_ctr) = names_of_rows.at(it_sort->second);
		sorted_info.at(row_ctr) = vector_PER_info.at(it_sort->second);	
		sorted_quals.at(row_ctr) = vector_PER_base_Qualities.at(it_sort->second);	
		++row_ctr;
	    }
	     
	     progressive_star_alignment = sorted_progressive_algmt;
	     prepared___progressive_star_alignment = sorted_progressive_algmt;
	     names_of_rows = sorted_names;
	     vector_PER_info = sorted_info;
	     vector_PER_base_Qualities = sorted_quals;
	}//sort
	
	
	//crop
	if (crop_to_relevant_region)
	{        
	    earliest_start_pos = std::max<int>(0, earliest_start_pos - 10);
	    latest_end_pos = std::min<int>(star_profile_size - 1,  latest_end_pos + 10);	    
	    
	    assert(earliest_start_pos <= latest_end_pos);
	    
	    for (uint row = 0; row < prepared___progressive_star_alignment.size(); ++row)
	    {
		prepared___progressive_star_alignment.at(row).assign(prepared___progressive_star_alignment.at(row).begin()  +  earliest_start_pos,
									prepared___progressive_star_alignment.at(row).begin()  +  latest_end_pos + 1);        
	    }                       
	    
	    prepared_crop_amount_below = earliest_start_pos;
	    prepared_crop_amount_above = (star_profile_size - 1) - latest_end_pos;                                                                    
	}
    }    
    
    
    return prepared___progressive_star_alignment;
    
}//sort_by_first_alignment_position_in_profile__and__crop