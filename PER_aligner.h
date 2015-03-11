#ifndef PER_ALIGNER_H_INCLUDED
#define PER_ALIGNER_H_INCLUDED


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


#include <api/BamAlignment.h>



#include <mpreal.h>
#include <mpfr.h>


#include <general_typedefs.h>

#include "globals.h"


class Paired_end_read;
class MiniRegion;



class PER_aligner
{                
    public:
        PER_aligner()  :  my_PER_ptr(NULL)
        { }            
                
        PER_aligner(const Paired_end_read *const &in_PER_ptr)
                          :  my_PER_ptr(in_PER_ptr)                          
        {  }             
        
        //pseudo-copy constructor
        PER_aligner(const PER_aligner &another_PER_aligner, const Paired_end_read *const new_PER_ptr)
                                    :   my_PER_ptr(new_PER_ptr)                                                                          
        {   }             
                        
        ~PER_aligner()
        { }                            
           
                              

        //variables:
        const Paired_end_read *my_PER_ptr; 
        
        
            
            
        //functions:         
        real calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0 ,
                            type_BI__BI &subset_of_region_that_is_aligned_to_mates__A__B___base_0)    const;
                            
        real calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0)    const;
                            
        real calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt)    const;             
                            
        real calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0)    const;                            
                            
        real calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion)    const;                                    
                            
        	    
	static std::string convert_uint_to_nucleotide_context
					(const uint &context_val);   
        
	static uint convert_nucleotide_to_index
                    (const char &some_nucleotide);  

        static char convert_index_to_nucleotide
                            (const uint &some_index); 
           
					

        //enum
        enum  State_type
        {
            state__NULL = -1,    state__Match = 0,    state__Ins_in_PER  = 1,    state__Del_in_PER = 2  
        };
        
        

    private:                 
                                            
        static uint convert_nucleotide_context_to_uint
                                    ( const std::string &prepared_REFerence,
                                      const int &position___base_1);                                  
                                    			                                                       
                                        
                                                                   
        void  appropriately_orient_mate_and_reference___and___complement_mate_and_determine_strandedness_and_ordering
                           (const BamTools::BamAlignment &the_naive_mate,
                            const bool &is_stored_as_first_naive_mate,
                            
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            
                            const std::string &in_reference, 
                            
                            //output                            
                            std::string &prepared_read,
                            type_vector_int &prepared_read_observed_nucleotide_pmfs,

                            bool &is_prepared_mate_A,
                            
                            std::string &prepared_REFerence,
                            bool &prepared_REFerence_is_complement_of_original_REFerence)    const;                                                                      
                                        
                                      

                                                
        void  do_Viterbi_traceback___and__get_lower_endpoint_on_REFerence
                                (const State_type *const *const  &DP_table___Viterbi_states,
                                
                                const std::string &prepared_REFerence,
                                    
                                const std::string &prepared_read,        
                                 
                                const bool &global__true____gloCAL_false,                                                                  
                                const int &ref_index_of_Viterbi_rightmost_nucleotide,
                                 
                                int &y_0,
                                const bool &return_Viterbi_sequence,
                                type_string__string  &PER_and_hybrid_algmt)   const;                                         
                            
                            
                            
        void  perform_sum_Forward__on__prepared_read_and_prepared_REF_____PRECISE___ie___GLOBAL
                                       (longdouble *const *const  &DP_table___Forward__M,
                                        longdouble *const *const  &DP_table___Forward__D,
                                        longdouble *const *const  &DP_table___Forward__I,

                                        longdouble *const *const  &DP_table___Viterbi_log10probabilities,
                                        State_type *const *const  &DP_table___Viterbi_states,
                                        
                                        const std::string &prepared_REFerence,
//                                         const type_vector_int  &amount_inside_relevant_homopolymers_along_REFerence___base_1,
                                         
                                        const std::string &prepared_read,                                      
                                        const type_vector_int  &prepared_read_observed_base_pmfs,
                                        
                                        const bool &return_Viterbi_alignment)    const;
                                        
                                        
                                        
                                        
        void  perform_sum_Forward__on__prepared_read_and_prepared_REF_____for_Viterbi_approximation___ie__gloCAL
                                       (longdouble *const *const  &DP_table___Viterbi_log10probabilities,
                                        State_type *const *const  &DP_table___Viterbi_states,
                                        
                                        const std::string &prepared_REFerence,
//                                         const type_vector_int  &amount_inside_relevant_homopolymers_along_REFerence___base_1,
                                         
                                        const std::string &prepared_read,                                      
                                        const type_vector_int  &prepared_read_observed_base_pmfs)    const;                                        
                                                   
       
       
        real PRECISELY_align_single_mate_of_PER_to_PRECISE_reference_region
                            (const BamTools::BamAlignment &the_mate,
                             const bool &is_stored_as_first_mate,
                             const std::string &m_reference_PRECISE,
                             const bool &orientation_of_PER_and_MiniRegion_agree,
                             const bool &complement_PER_to_get_to_MiniRegion,
                             
                             type_string__string &PER_and_hybrid_algmt,
                             const bool &return_Viterbi_alignment    )    const;
                                    
                                    
        void get_Viterbi_alignment_and_endpoints_of_reference_in_Viterbi_alignment
                            (const BamTools::BamAlignment &mate_A,//const std::string &in_read_A,
                            const BamTools::BamAlignment &mate_B,//const std::string &in_read_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0,
                            type_BI__BI &subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                            const bool &return_Viterbi_alignment_sequence)    const;

                                           
                 
                                        
        real calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference____IMPL
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0,
                            type_BI__BI &subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                             const bool &return_Viterbi_alignment_sequence)    const;                                        
                                        
                                    
		 
    public:
		 
	longdouble  calculate_P_transition__cond__everything_else
		    (const uint &cycle,
		     const std::string &prepared_reference_sequence,
		     const uint &reference_position__base_1,		     
		     const State_type &state_FROM,
		     const State_type &state_TO  )   const;

		 
	longdouble  calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
		    (const std::string &prepared_read_sequence,
		     const uint &position_along_read___base_1,
		     const type_vector_int &prepared_sequence_base_Qualities,
		     const std::string &prepared_reference_sequence,
		     const uint &reference_position__base_1,		     	     
		     const State_type &state_of_emission)  const;		 


		     
	uint  identify_homopolymer_length_at_given_position
				(const std::string &prepared_REFerence,
				 const uint &position_along_reference___base_1)   const;     
		     

				 
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {              
            //do nothing
        };     
	
	
};//PER_aligner







#endif