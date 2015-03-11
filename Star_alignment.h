#ifndef STAR_ALIGNMENT_H_INCLUDED
#define STAR_ALIGNMENT_H_INCLUDED


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
#include <templates.h>
#include <translations_through_alignment.h>


#include "globals.h"








                    
                    
class Star_alignment
{
    public:
        Star_alignment()   :   progressive_star_alignment(),
                               star_profile_size(0),
                               failed_alignment(true),
                               LCR_0__LCR_1__XY___YX___XYX___YXY(-1),
                               activated_star_alignment(false)
        { }
        
        
        
        Star_alignment( const type_vector_int &in_foundation_coords,
                         const std::string &in_foundation_seq,
                         const std::string &foundation_seq_name,
			const BOOST_Interval &foundation_index_brkpts,
			const bool &is_NAHR,
			const int &in__LCR_0__LCR_1__XY___YX___XYX___YXY)	
        {
	    init_Star_alignment
			    (in_foundation_coords,
			    in_foundation_seq,
			    foundation_seq_name,
			    foundation_index_brkpts,
			    is_NAHR,
			    in__LCR_0__LCR_1__XY___YX___XYX___YXY);
        }    
                        
        
        
        
        ~Star_alignment()
        { } 
        
        
        
        
        
        //static variables:
// 	static std::set<char>  allowable_values;
        
        
        //variables:        
                
        type_vector_vector_int  progressive_star_alignment;    
                                    //organization:
                                    //                  0    -     coordinates              (e.g.:    736  737   738   739   740  750    )
                                    //                  1    -     foundation sequence      (e.g.:     A    C     C     G     T    G     )
                                    //                  2    -     some Paired-end read     (e.g.:                C     G     T    x    x    x    A   G    T     )
                                    //                  ...
                        
        
        int star_profile_size; // i.e. cols
        
        
        type_vector_string  names_of_rows;
        
        int star_profile_index_of_first_Y_nucleotide_in_hybrid_XY;
	int star_profile_index_of_first_second_X_nucleotide_in_hybrid_XYX;
        
      
	type_vector_string  vector_PER_info;
	type_vector_string  vector_PER_base_Qualities;
	
	
	bool failed_alignment;
	
	bool activated_star_alignment;
	
	
	int LCR_0__LCR_1__XY___YX___XYX___YXY;
        
        
	
        //functions:
        
        void print_this_Star_alignment(std::stringstream *const &some_ss = NULL)  const;
	

	void init_Star_alignment
			    (const type_vector_int &in_foundation_coords,
			    const std::string &in_foundation_seq,
			    const std::string &foundation_seq_name,
			    const BOOST_Interval &foundation_index_brkpts,
			    const bool &is_NAHR,
			    const int &in__LCR_0__LCR_1__XY___YX___XYX___YXY);	
        
	
	bool add_pairwise_alignment_to_progressive_star_alignment
                                (const real &P_PER_cond_this_sequence,
                                    const type_string__string &PER_and_hybrid_algmt,
                                    const type_map_uint_to_uint &absolute_coords_hybrid_X_to_pairwise_algmt_profile,
                                    const type_map_uint_to_uint &absolute_coords_hybrid_Y_to_pairwise_algmt_profile,                             
                                    const bool &orientation_of_pairwise_algmt_agrees_with_foundation_sequence,
                                    const bool &must_complement_algmt_to_get_to_foundation_sequence,
				    type_Balgmt__Balgmt  mates_of_PER,   //NOT BY REFERENCE
				    const real &sum_over_all_P_locations);
                            
                                    
         bool create_HTML_output_for_star_alignment
                                                (const int &chromo_of_region_X,
                                                 const int &chromo_of_region_Y,
                                                 const bool &crop_to_relevant_region,
                                                 type_vector_string &star_html_output);
						 
	int get_total_PER_mate_index_at_star_index
			    (const type_vector_vector_int  &prepared___progressive_star_alignment,
			     const int &row_index_for_PER_in_prepared__progressive_alignment,
			     const int &base_index_in_prepared__progressive_alignment);						 
                              
        void clear_all_data();                                 
                
        
    private:
        void expand_star_profile_for_inserted_bases_in_a_PER
                                                (int &universal_star_index,
                                                 const type_vector_vector_int::iterator &newest_PER_on_star_profile,
                                                 const char &some_char   );

        void expand_star_profile_for_inserted_bases_in_a_PER
                                                (int &universal_star_index,
                                                 const type_vector_vector_int::iterator &newest_PER_on_star_profile,
                                                 const std::string &insertion_string   );
						
	type_vector_vector_int sort_by_first_alignment_position_in_profile__and__crop
						     (const bool &crop_to_relevant_region,
						      int &prepared_crop_amount_below,
						     int &prepared_crop_amount_above);						     
                                    
};
                    
      
typedef  std::pair<BOOST_Interval, Star_alignment>    type_BI__Star_alignment;      
typedef  std::list<type_BI__Star_alignment>      type_list_BI__Star_alignment;
typedef  std::map<uint, type_list_BI__Star_alignment>    type_map_uint_to_list_BI__Star_alignment;
                    
                    
                    
                    
#endif