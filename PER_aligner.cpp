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

#include <array>

#include <omp.h>

#include <stdexcept>


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

#include <api/BamAlignment.h>


#include <mpreal.h>
#include <mpfr.h>




#include <general_typedefs.h>
#include <templates.h>
#include <translations_through_alignment.h>

#include "PER_aligner.h"
#include "globals.h"
#include "Paired_end_read.h"
#include "MiniRegion.h"

#include "PER_bundle_wrapper.h"


#include "other_functions.h"






#define    max_advance_size   (100000)
#define    min_advance_size   (-100000)



//IMPLEMENTATION NOTES:
//  1) Naming convention:  "mate A": refers to the read which "comes first" wrt the indexing (coordinates) of the refernece.  The strand (forward vs .reverse) of the Reference which generated the mates is irrelevant as far as the designation of  "mate A"  vs.  "mate B"  is concerned.
//
//      Reference coordinates:                  m    m+1    m+2    m+3    m+4    m+5    m+6    m+7    m+8    m+9    m+10    m+11    m+12    m+13
//                                                   --------------------->                                  <--------------------------
//                                                          Mate A                                                      Mate B
//
//
// 2)  Naming convention:   a "prepared mate" has been orientated/complemtned approriately so that it "matches" the input reference sequence.
// 3)  Naming convention:   a "prepared reference" has been reversed if in the context of mate B (since mate B is always generated in the descending direction), and has been appropriately complemented in the context of either mate (according to the strand which generated the mate). 
// 4)  To protect against bugs/coding errors, I have chosen to only use ONE algorithm for forward/backward sums.  Thus, for mate B, which was created in a descending direction (bases elongation began at positon n, then n-1, the n-2, ..., 1), the reference is REVERSED, so that the forward and backward algorithms are the same for mate A as for mate B.  Further, for (appropriately) both mate A and mate B, the reference should be COMPLEMENTED if the mate was generated using the reverse strand 		// Thus, all of the aligning procedure is the same for both mates.  We only have to be careful at the very end - we must be sure to un-reverse and un-complement the refernece sequence and coordinates as appropriate.  This must always be done for mate B, and for mate A, it depends which strand (forward or reverse) generated mate A.





















//static variables:   







uint PER_aligner::convert_nucleotide_to_index
                    (const char &some_nucleotide)
{
    switch(some_nucleotide)
    {
        case 'A':
            return 0;
            break;
        case 'C':
            return 1;
            break;
        case 'G':
            return 2;
            break;
        case 'T':
            return 3;
            break;
	case 'N':
	    return 4;
	    break;
        default:

            std::stringstream error_strm;
            print_line_of_markers("ERROR! ", &error_strm);
            print_line_of_markers("(", &error_strm);
            error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
            
            error_strm << "\n\nERROR!    unable to convert nucleotide [" << some_nucleotide << "]  to an index  in  \"PER_aligner::convert_nucleotide_to_index\"\n\n";
            
            print_line_of_markers(")", &error_strm);
            std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                            
            abort();  exit(1);
	    break;
    };
}  // PER_aligner::convert_nucleotide_to_index








char PER_aligner::convert_index_to_nucleotide
                    (const uint &some_index)
{
    switch(some_index)
    {
        case 0:
            return 'A';
            break;
        case 1:
            return 'C';
            break;
        case 2:
            return 'G';
            break;
        case 3:
            return 'T';
            break;
	case 4:
	    return 'N';
	    break;
        default:
            std::stringstream error_strm;
            print_line_of_markers("ERROR! ", &error_strm);
            print_line_of_markers("(", &error_strm);
            error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
            
            error_strm << "\n\nERROR!    unable to convert index [" << some_index << "]  to a nucleotide  in  \"PER_aligner::convert_index_to_nucleotide\"\n\n";
            
            print_line_of_markers(")", &error_strm);
            std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                            
            abort(); exit(1);
            break;
    };

}  // PER_aligner::convert_index_to_nucleotide








uint PER_aligner::convert_nucleotide_context_to_uint
                                    ( const std::string &prepared_REFerence,
                                      const int &position___base_1)
{
    // you add 1 or 5 or 21 so that all contexts can be converted to an integer and not overlap.  this is easier for indexing
    switch (position___base_1)
    {
        case 0:   //this is the case for the "before first generated base" scenario
            return 0;
            break;            
            
        case 1:
            return convert_nucleotide_to_index( prepared_REFerence.at(0) ) + 1;
            break;            
            
        case 2:
            return  convert_nucleotide_to_index( prepared_REFerence.at(0) )
                    +  (   convert_nucleotide_to_index( prepared_REFerence.at(1) )   <<  2   ) + 5;
            break;            
            
        default:
            return      convert_nucleotide_to_index( prepared_REFerence.at(position___base_1-3) ) 
                        +    (    convert_nucleotide_to_index( prepared_REFerence.at(position___base_1-2) ) << 2   )
                        +    (    convert_nucleotide_to_index( prepared_REFerence.at(position___base_1-1) ) << 4   )  +  21; 
            break;                                                                
    } // switch 
              
};





// Here we apply a "mask"
std::string PER_aligner::convert_uint_to_nucleotide_context
					(const uint &context_val)
{
    
    if (context_val == 0)
    {  return std::string();  }
    else if (context_val < 5)
    {  return std::string(1, convert_index_to_nucleotide(context_val-1));  }
    else if (context_val < 21)
    {
	const uint adjusted_context_val = context_val - 5;
	const uint zero_binary_digit = adjusted_context_val & 3;
	const uint first_binary_digit =   (  adjusted_context_val & (3 << 2)   )   >>  2;
	
	std::string context_nucleotides;
	context_nucleotides.push_back(  convert_index_to_nucleotide(zero_binary_digit)  );
	context_nucleotides.push_back(  convert_index_to_nucleotide(first_binary_digit)  );	    
	
	return context_nucleotides;
    }
    else
    {
	const uint adjusted_context_val = context_val - 21;
	
	const uint zero_binary_digit = adjusted_context_val & 3;
	const uint first_binary_digit =     (  adjusted_context_val & (3 << 2)   )   >>  2;	
	const uint second_binary_digit =    (  adjusted_context_val & (3 << 4)   )   >>  4;
	
	std::string context_nucleotides;
	context_nucleotides.push_back(  convert_index_to_nucleotide(zero_binary_digit)  );
	context_nucleotides.push_back(  convert_index_to_nucleotide(first_binary_digit)  );
	context_nucleotides.push_back(  convert_index_to_nucleotide(second_binary_digit)  );
	
	return context_nucleotides;	    
    }

} // PER_aligner::convert_uint_to_nucleotide_context


















uint  PER_aligner::identify_homopolymer_length_at_given_position
				(const std::string &prepared_REFerence,
				 const uint &position_along_reference___base_1)  const
// 				 const bool &use_expanded_simple_repeat_definition)
{
    //use_expanded_simple_repeat_definition
    if (position_along_reference___base_1 < 3)  //irrelevant - declared not long enough to have a homopolymer bias.
    {  return 1;  }
    
    
    int current_ref_checking_pos___base_1;//check the substring starting at this base_1 pos and of length = "query_seq.size()"
    std::string query_seq;
    {
	const std::string ref_triplet(prepared_REFerence.substr((int)position_along_reference___base_1-1-2, 3));
	
	if (ref_triplet.at(0) == ref_triplet.at(1)  and ref_triplet.at(0) == ref_triplet.at(2))  // ZZZ
	{
	    query_seq = ref_triplet.at(0);
	    current_ref_checking_pos___base_1 = ((int)position_along_reference___base_1);  // inefficient, but convenient.
	}
	else if (ref_triplet.at(0) == ref_triplet.at(2))  // ZYZ
	{
	    query_seq = ref_triplet.substr(1,2);
	    current_ref_checking_pos___base_1 = ((int)position_along_reference___base_1) - 1;
	}
	else // XYZ
	{
	    query_seq = ref_triplet;
	    current_ref_checking_pos___base_1 = ((int)position_along_reference___base_1) - 2;
	}
    }
    
    
    const int query_size = query_seq.size();               
    uint length_of_simple_repeat = 0;
    
    while ( current_ref_checking_pos___base_1 > 0  
	    and
	    query_seq.compare(prepared_REFerence.substr(current_ref_checking_pos___base_1-1, query_size)) == 0)
    {
	current_ref_checking_pos___base_1 -= query_size;
	++length_of_simple_repeat;
    }
    
    
    return   std::max<uint>(1,length_of_simple_repeat);

}//identify_homopolymer_length_at_given_position




























































//Principles:
// 1.  First, the Reference is fixed, and the read is oriented&complemented to match the Reference. (e.g. if the read came from another member of the repeat family to which this region of the Reference genome belongs).
// 2.  We are concnered with how contexts (e.g. triplets) GENERATE the read.  Thus, once the read is oriented&complemented to the Reference, we should complement the Refernece (if necessary) so that it looks like the sequence which generated the read.  Hence, it matters if the read looks like the forward or reverse strand of the Reference (after orienting&complementing the read initially).



//If someday you try to use single-ended reads, this function will have to be modified so that "is_mate_A" is always true (i.e. the prepared_REF is never reversed).  ( I think that is all you need to do .....maybe..... ?????)

void  PER_aligner::appropriately_orient_mate_and_reference___and___complement_mate_and_determine_strandedness_and_ordering
                           (const BamTools::BamAlignment &the_naive_mate,
                            const bool &is_stored_as_first_naive_mate,
                            
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            
                            const std::string &in_reference,   
                            
                            //output                            
                            std::string &prepared_read,
                            type_vector_int &prepared_read_observed_nucleotide_qualities,

                            bool &is_prepared_mate_A,
                            
                            std::string &prepared_REFerence,
                            bool &prepared_REFerence_is_complement_of_original_REFerence)    const
{
// We assume the Illumina model:
//                                                   --------------------->                                  <--------------------------
//                                                          Mate A                                                      Mate B

    
    // Naive assignment:
    prepared_read.assign(the_naive_mate.QueryBases);
    
    
    prepared_read_observed_nucleotide_qualities.resize(the_naive_mate.Qualities.size(), -1);
    for (uint j=0; j< the_naive_mate.Qualities.size(); ++j)	    
    {  					
	prepared_read_observed_nucleotide_qualities[j] = (int)the_naive_mate.Qualities[j] - ASCII_quality_offset;
    }
    
    
    
    

    
    
    bool the_bases_emitted_during_read_generation_looked_like_the_FORWARD_strand___ie_the_reverse_strand_served_as_a_template
			    = ! the_naive_mate.IsReverseStrand();    
	    //as always, the bases of the read were approrpriately complemented by the industrial aligner (e.g. BWA) and so the bases given (in "the_naive_mate") ALWAYS look like the forward strand.  the point of this boolean variable is to determine what the ORIGINAL bases of naive_mate looked like, BEFORE the industrial aligner appropriately complemented them for convenience.
	    
    
    //Naive initialization:
    is_prepared_mate_A = is_stored_as_first_naive_mate;                    
    prepared_REFerence_is_complement_of_original_REFerence = false;
            
    
    if (!orientation_of_PER_and_MiniRegion_agree)
    {
        prepared_read = get_reverse_of_sequence( prepared_read );
        prepared_read_observed_nucleotide_qualities = get_reverse_of_vector<int>( prepared_read_observed_nucleotide_qualities );
	
        is_prepared_mate_A = !is_prepared_mate_A;
    }
    
    

    if ( complement_PER_to_get_to_MiniRegion )
    {
        prepared_read = get_complement_of_sequence(prepared_read);   
        the_bases_emitted_during_read_generation_looked_like_the_FORWARD_strand___ie_the_reverse_strand_served_as_a_template
		=  ! the_bases_emitted_during_read_generation_looked_like_the_FORWARD_strand___ie_the_reverse_strand_served_as_a_template;
    } 
    
    
    
    //Now the read is orented&complemented to match the reference. 
    // We then format the REFERENCE so that all alignments are in the forward direction and all triplet motifs are according to the forward strand of the reference genome.
    
// We assume the Illumina model:
//                                                   --------------------->                                  <--------------------------
//                                                          Mate A                                                      Mate B
//                                                                          REF    
    
    if (is_prepared_mate_A)
    {         
        if (the_bases_emitted_during_read_generation_looked_like_the_FORWARD_strand___ie_the_reverse_strand_served_as_a_template)
	{
            prepared_REFerence.assign(in_reference);
	}
        else
        {
	    prepared_read.assign(get_complement_of_sequence(prepared_read)); // now the read bases really do look like those generated.
	    
            prepared_REFerence.assign(get_complement_of_sequence(in_reference));
            prepared_REFerence_is_complement_of_original_REFerence = true;
        }
    }
    else // mate B - then the read & Reference must be reversed  so that all alignment algorithms are the same (i.e. simple)
    {	
        if (the_bases_emitted_during_read_generation_looked_like_the_FORWARD_strand___ie_the_reverse_strand_served_as_a_template)
	{
	    prepared_read = get_reverse_of_sequence(prepared_read);	    
	    prepared_read_observed_nucleotide_qualities = get_reverse_of_vector<int>(prepared_read_observed_nucleotide_qualities);
	    
	    prepared_REFerence.assign(get_reverse_of_sequence(in_reference));	
	}
        else
        {
	    prepared_read.assign(get_reverse_complement_of_sequence___ie_inversion(prepared_read));
	    prepared_read_observed_nucleotide_qualities = get_reverse_of_vector<int>(prepared_read_observed_nucleotide_qualities);
	    
            prepared_REFerence.assign(get_reverse_complement_of_sequence___ie_inversion(in_reference));            
            prepared_REFerence_is_complement_of_original_REFerence = true;
        }
    }
         
         
    PER_bundle_wrapper__space::remove_N_tails_from_read_and_quality(prepared_read, prepared_read_observed_nucleotide_qualities);    
            

}  // appropriately_orient_mate_and_reference___and___complement_mate_and_determine_strandedness_and_ordering







































    

void  PER_aligner::do_Viterbi_traceback___and__get_lower_endpoint_on_REFerence
                                (const State_type *const *const  &DP_table___Viterbi_states,
                                
                                const std::string &prepared_REFerence,
                                    
                                const std::string &prepared_read,
                                 
                                const bool &global__true____gloCAL_false,                                                                  
                                const int &ref_index_of_Viterbi_rightmost_nucleotide,
                                 
                                int &y_0,
                                const bool &return_Viterbi_sequence,
                                type_string__string  &PER_and_hybrid_algmt)   const
{            
    //Viterbi traceback:      
    y_0 = -1;                    
    
    
    int traceback_indx__REFerence = ref_index_of_Viterbi_rightmost_nucleotide;
    int traceback_indx__read = prepared_read.size();
    
    type_string__string PER_and_hybrid_algmt__backwards; 
    PER_and_hybrid_algmt__backwards.first.reserve( prepared_read.size()+20 );
    PER_and_hybrid_algmt__backwards.second.reserve( prepared_read.size()+20 );
    
    
    
    while (traceback_indx__read > 0  and    traceback_indx__REFerence > 0)
    {
        switch (  DP_table___Viterbi_states[traceback_indx__read][traceback_indx__REFerence]  )
        {
            case state__Match:
            {
                PER_and_hybrid_algmt__backwards.first.push_back(  prepared_read.at(traceback_indx__read-1)  );
                PER_and_hybrid_algmt__backwards.second.push_back( prepared_REFerence.at(traceback_indx__REFerence-1)   );  
                --traceback_indx__read;
                --traceback_indx__REFerence;                        
                break;
            }
            case state__Del_in_PER:
            {
                PER_and_hybrid_algmt__backwards.first.push_back('-');
                PER_and_hybrid_algmt__backwards.second.push_back( prepared_REFerence.at(traceback_indx__REFerence-1)   );  
                --traceback_indx__REFerence;
                break;
            }
            case state__Ins_in_PER:  
            {
                PER_and_hybrid_algmt__backwards.first.push_back(  prepared_read.at(traceback_indx__read-1)  );
                PER_and_hybrid_algmt__backwards.second.push_back('-');
                --traceback_indx__read;                                                                    
                break;
            }
            default:
            {
                std::stringstream error_strm;
                print_line_of_markers("ERROR! ", &error_strm);
                print_line_of_markers("(", &error_strm);
                error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                
                error_strm << "\n\nERROR!   unrecognized traceback code  in \"do_Viterbi_traceback___and__get_lower_endpoint_on_REFerence\"!!!!\n\n";
                
                error_strm << "\nDP_table_A__Viterbi_states[traceback_indx__read][traceback_indx__REFerence]  =  " 
                            << DP_table___Viterbi_states[traceback_indx__read][traceback_indx__REFerence]
                                << "\ntraceback_indx__read = " << traceback_indx__read << "\n"
                                << "\ntraceback_indx__REFerence = " << traceback_indx__REFerence << "\n";
                                
                error_strm << "\nprepared_read = \n" << prepared_read
                            << "\nprepared_REFerence = \n" << prepared_REFerence
                            << "\n\nprepared_REFerence.size() = " << prepared_REFerence.size() << "\n\n";
			    
			    
		error_strm << "\n\nsome DP_table___Viterbi_states...\n\n";
			    
		for (int v=std::max<int>(0,traceback_indx__REFerence-20); v <= std::min<int>(traceback_indx__REFerence+20, (int)prepared_REFerence.size()); ++v)
		{   error_strm << "\tDP_table___Viterbi_states[" << traceback_indx__read << "][" << v << "] = "
				<< DP_table___Viterbi_states[traceback_indx__read][v] << "\n";  }
                            
//                             error_strm << "\nPER name = " << my_PER_ptr->name << "\n\n";  
//                         my_PER_ptr->print_this_PER( &error_strm );     
//                 error_strm << "\n\norientation_of_PER_and_MiniRegion_agree = " << orientation_of_PER_and_MiniRegion_agree << "\n";
//                 error_strm << "\n\ncomplement_PER_to_get_to_MiniRegion = " << complement_PER_to_get_to_MiniRegion << "\n";
                
//                 print_Balgmt(mate_A, &error_strm);
//                 print_Balgmt(mate_B, &error_strm);
                
//                 my_PER_ptr->print_this_PER( &error_strm );
                
//                         print_2d_array<int>(DP_table_A__Viterbi_states,
//                                             prepared_read_A_length+1, 
//                                             ref_length_minus_1 + 1,
//                                             "DP_table_A__Viterbi_states");   
//                                             
//                         print_2d_array<longdouble>(DP_table_A__Viterbi_log10probabilities,
//                                                     prepared_read_A_length+1, 
//                                                     ref_length_minus_1 + 1,
//                                                     "DP_table_A__Viterbi_log10probabilities"); 
                
                print_line_of_markers(")", &error_strm);
                std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
                abort();  exit(1);;
                break; 
            }
        }//switch
    }
        
        
        
        
        
        
    y_0 = traceback_indx__REFerence + 1;  // +1 because the alignment ALWAYS ends in a match
    
    
    
    
    
    if (return_Viterbi_sequence)
    {                
        if (traceback_indx__read > 0)  // then there is still a string of insertions to add
        {
            PER_and_hybrid_algmt__backwards.first.append(  get_reverse_of_sequence( prepared_read.substr(0, traceback_indx__read) )   ); // add backwards, for consistency.
            PER_and_hybrid_algmt__backwards.second.append(  std::string(traceback_indx__read, '-')   );
        }    
        else if (global__true____gloCAL_false  and  traceback_indx__REFerence > 0   )
        {
            PER_and_hybrid_algmt__backwards.first.append(  std::string(traceback_indx__REFerence, '-')   );
            PER_and_hybrid_algmt__backwards.second.append(   get_reverse_of_sequence( prepared_REFerence.substr(0, traceback_indx__REFerence)  )  ); // add backwards, for consistency.       
            y_0 = 1;        
        }            
            
        PER_and_hybrid_algmt.first = get_reverse_of_sequence(PER_and_hybrid_algmt__backwards.first);        
        PER_and_hybrid_algmt.second = get_reverse_of_sequence(PER_and_hybrid_algmt__backwards.second);           
    }// return_Viterbi_sequence
    
   
} //  do_Viterbi_traceback___and__get_endpoints_on_REFerence
























void  PER_aligner::perform_sum_Forward__on__prepared_read_and_prepared_REF_____PRECISE___ie___GLOBAL
                                       (longdouble *const *const  &DP_table___Forward__M,
                                        longdouble *const *const  &DP_table___Forward__D,
                                        longdouble *const *const  &DP_table___Forward__I,

                                        longdouble *const *const  &DP_table___Viterbi_log10probabilities,
                                        State_type *const *const  &DP_table___Viterbi_states,
                                        
                                        const std::string &prepared_REFerence,
                                         
                                        const std::string &prepared_read,                                        
                                        const type_vector_int  &prepared_read_observed_base_Qualities,
                                        
                                        const bool &return_Viterbi_alignment)    const
{
    //assume DP_tables are already allocated correctly.
                                                                            
    //init A    

    //[0][0]
    DP_table___Forward__M[0][0] = 1.00L;
    DP_table___Forward__D[0][0] = 0.00L;
    DP_table___Forward__I[0][0] = 0.00L;
    
    DP_table___Viterbi_log10probabilities[0][0] = 0.00L;
    DP_table___Viterbi_states[0][0] = state__Match;
    
    	
    

    
    //0th column:    
    DP_table___Forward__M[1][0] = 0.00L;
    DP_table___Forward__D[1][0] = 0.00L;
    DP_table___Forward__I[1][0] =  
                                    calculate_P_transition__cond__everything_else
						    (0,
						    prepared_REFerence,
						    0,		     	    
						    state__Match,
						    state__Ins_in_PER)   
				    * calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
					    (prepared_read,
					     1,
					    prepared_read_observed_base_Qualities,
					    prepared_REFerence,
					    0,		     
					    state__Ins_in_PER);	
    					    
    DP_table___Viterbi_log10probabilities[1][0] = log10(  DP_table___Forward__I[1][0]  );
    DP_table___Viterbi_states[1][0] = state__Ins_in_PER;                                                                                 
                            
    //rest of 0-column
    for (int i = 2; i <= prepared_read.size(); ++i)
    {
        DP_table___Forward__M[i][0] =  0.00L;
        DP_table___Forward__D[i][0] =  0.00L;
        
        DP_table___Forward__I[i][0] =    DP_table___Forward__I[i-1][0]
                                    *  calculate_P_transition__cond__everything_else
						    (i-1,
						    prepared_REFerence,
						    0,		     	    
						    state__Ins_in_PER,
						    state__Ins_in_PER)
				    * calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
					    (prepared_read,
					     i,
					    prepared_read_observed_base_Qualities,
					    prepared_REFerence,
					    0,		     
					    state__Ins_in_PER);
                                                            
        DP_table___Viterbi_log10probabilities[i][0] = log10(  DP_table___Forward__I[i][0]  );
        DP_table___Viterbi_states[i][0] = state__Ins_in_PER;                                                                
    }
        



    //0th row:
    DP_table___Forward__M[0][1] = 0.00L;
    DP_table___Forward__D[0][1] =  calculate_P_transition__cond__everything_else
						    (0,
						    prepared_REFerence,
						    0,
						    state__Match,
						    state__Del_in_PER);
    DP_table___Forward__I[0][1] = 0.00L;   
    
    DP_table___Viterbi_log10probabilities[0][1] = log10(  DP_table___Forward__D[0][1]  );
    DP_table___Viterbi_states[0][1] = state__Del_in_PER;      
    
    
    //rest of 0th row
    for (int j = 2 ; j <= prepared_REFerence.size(); ++j)
    {
        DP_table___Forward__M[0][j] =  0.00L;
        DP_table___Forward__I[0][j] =  0.00L;        
        DP_table___Forward__D[0][j] =   DP_table___Forward__D[0][j-1]
                                    *   calculate_P_transition__cond__everything_else
						    (0,
						    prepared_REFerence,
						    j-1,		     	    
						    state__Del_in_PER,
						    state__Del_in_PER);
						    
                                    
        DP_table___Viterbi_log10probabilities[0][j] = log10(  DP_table___Forward__D[0][j]  );
        DP_table___Viterbi_states[0][j] = state__Del_in_PER;                                      
    }
    
    
    
    //rest of DP matrix:    
    //fill                                                    
    for (int j = 1; j <= prepared_REFerence.size(); ++j)    
    {			
        for (int i = 1; i <= prepared_read.size(); ++i)        
        {   	        
	    
            const longdouble marginalized_Match_emission_for_this_position
				    =  calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
					    (prepared_read,
					     i,
					    prepared_read_observed_base_Qualities,
					    prepared_REFerence,
					    j,		     
					    state__Match); 
            {//to Match                                                
                DP_table___Forward__M[i][j] 
                    = (    
                            (    DP_table___Forward__M[i-1][j-1] 
                                    *   calculate_P_transition__cond__everything_else
						    (i-1,
						    prepared_REFerence,
						    j-1,		     	    
						    state__Match,
						    state__Match)
                                    *  marginalized_Match_emission_for_this_position
                            )  // Match to Match
                        +
                            (    DP_table___Forward__I[i-1][j-1] 
                                    *   calculate_P_transition__cond__everything_else
						    (i-1,
						    prepared_REFerence,
						    j-1,		     		    
						    state__Ins_in_PER,
						    state__Match)
                                    *  marginalized_Match_emission_for_this_position
                            )  // Ins to Match                            
                        +
                            (    DP_table___Forward__D[i-1][j-1] 
                                    *   calculate_P_transition__cond__everything_else
						    (i-1,
						    prepared_REFerence,
						    j-1,		     	    
						    state__Del_in_PER,
						    state__Match)
                                    *  marginalized_Match_emission_for_this_position
                            )  // Del to Match                            
                        );                    
            }// to Match
                    
                    
                    
            const longdouble marginalized_Insertion_emission_for_this_position
				    =  calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
					    (prepared_read,
					     i,
					    prepared_read_observed_base_Qualities,
					    prepared_REFerence,
					    j,		     
					    state__Ins_in_PER);                            
            {//to Insertion                                           
                DP_table___Forward__I[i][j] 
                    = (
                            (   DP_table___Forward__M[i-1][j] 
                                    *   calculate_P_transition__cond__everything_else
						    (i-1,
						    prepared_REFerence,
						    j,		     
						    state__Match,
						    state__Ins_in_PER)
                                        *  marginalized_Insertion_emission_for_this_position
                            )
                        +   
                            (   DP_table___Forward__I[i-1][j] 
                                    *   calculate_P_transition__cond__everything_else
						    (i-1,
						    prepared_REFerence,
						    j,		         
						    state__Ins_in_PER,
						    state__Ins_in_PER)
                                        *  marginalized_Insertion_emission_for_this_position                                                        
                            )                                        
                        );                        
            }//to Insertion
          
                    
            
            {// to Deletion                                
                DP_table___Forward__D[i][j] 
                    =  (
                            (    DP_table___Forward__M[i][j-1] 
                                    *   calculate_P_transition__cond__everything_else
						    (i,
						    prepared_REFerence,
						    j-1,		     
						    state__Match,
						    state__Del_in_PER)  ) 
                            +
                            (    DP_table___Forward__D[i][j-1] 
                                    *   calculate_P_transition__cond__everything_else
						    (i,
						    prepared_REFerence,
						    j-1,		       
						    state__Del_in_PER,
						    state__Del_in_PER)  )                                                                                                   
                        );                       
            }// to Deletion 
            
            
            
            
            
            
            
            
           
            
            // Viterbi:
            if (return_Viterbi_alignment)
            {
                
                //init to match
                longdouble best_log10_prob                
                            =  DP_table___Viterbi_log10probabilities[i-1][j-1]
                                   + log10(  calculate_P_transition__cond__everything_else
						    (i-1,
						    prepared_REFerence,
						    j-1,		     
						    DP_table___Viterbi_states[i-1][j-1],
						    state__Match)       )
                                +  log10(  marginalized_Match_emission_for_this_position  );
                                                                                            
                State_type best_state_transition = state__Match;
                                            
                
                
                //ins
                const longdouble log10_P_ins
                            =   DP_table___Viterbi_states[i-1][j] == state__Del_in_PER  ?
                                        LOG10_OF_ZERO
                                        :                                        
                                        (   DP_table___Viterbi_log10probabilities[i-1][j]
					    +  log10(  calculate_P_transition__cond__everything_else
								(i-1,
								prepared_REFerence,
								j,		     
								DP_table___Viterbi_states[i-1][j],
								state__Ins_in_PER)    )
					    +  log10(  marginalized_Insertion_emission_for_this_position  )   );
                                        
                                                
                                                
                //del                        
                const longdouble log10p_del_here
                            =   DP_table___Viterbi_states[i][j-1] == state__Ins_in_PER  ?
                                        LOG10_OF_ZERO
                                        :                                        
                                        (   DP_table___Viterbi_log10probabilities[i][j-1]
						+  log10(  calculate_P_transition__cond__everything_else
								    (i,
								    prepared_REFerence,
								    j-1,		     
								    DP_table___Viterbi_states[i][j-1],
								    state__Del_in_PER)  ));                                                                       
                                

                if (log10_P_ins > best_log10_prob)
                {
                    best_state_transition = state__Ins_in_PER;
                    best_log10_prob = log10_P_ins;                                    
                }                    
                
                
                if (log10p_del_here > best_log10_prob)
                {
                    best_log10_prob = log10p_del_here;
                    best_state_transition = state__Del_in_PER;                
                }
                                                                                                                            
                
                
                //save
                DP_table___Viterbi_states[i][j] = best_state_transition;
                DP_table___Viterbi_log10probabilities[i][j] = best_log10_prob;    
                
            }//Viterbi 
                
     
        }// i      
    } // j
        
       

} // perform_sum_Forward__on__prepared_read_and_prepared_REF_____PRECISE___ie___GLOBAL





























	// [context]   -- triplet (including homopolymer), double, single, or nothing
	// [rep length]  -- only for homopolymer
	// [cycle]	
	// [state]   -- match = 0, insertion = 1	
	

void  PER_aligner::perform_sum_Forward__on__prepared_read_and_prepared_REF_____for_Viterbi_approximation___ie__gloCAL
                                       (longdouble *const *const  &DP_table___Viterbi_log10probabilities,
                                        State_type *const *const  &DP_table___Viterbi_states,
                                        
                                        const std::string &prepared_REFerence,
                                         
                                        const std::string &prepared_read,                                        
                                        const type_vector_int  &prepared_read_observed_base_Qualities)    const
{
    //assume DP_tables are already allocated correctly.
   

    //[0][0]    
    DP_table___Viterbi_log10probabilities[0][0] = 0.00L;
    DP_table___Viterbi_states[0][0] = state__Match;
                     
    
    //Must do [1][0] specially.                                                            
    DP_table___Viterbi_log10probabilities[1][0] =
					log10(   calculate_P_transition__cond__everything_else
									(0,
									prepared_REFerence,
									0,		     
									state__Match,
									state__Ins_in_PER)    )
					+  log10( calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
							(prepared_read,
							1,
							prepared_read_observed_base_Qualities,
							prepared_REFerence,
							0,		     
							state__Ins_in_PER)  );
    DP_table___Viterbi_states[1][0] = state__Ins_in_PER;
								       
								       
    
    //rest of 0-column
    for (int i = 2; i <= prepared_read.size(); ++i)
    {    
        DP_table___Viterbi_log10probabilities[i][0] =  DP_table___Viterbi_log10probabilities[i-1][0]
							+  log10(   calculate_P_transition__cond__everything_else
									(i-1,
									prepared_REFerence,
									0,		         
									state__Ins_in_PER,
									state__Ins_in_PER)	)        
                                                       +  log10(  calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
								    (prepared_read,
								    i,
								    prepared_read_observed_base_Qualities,
								    prepared_REFerence,
								    0,		     
								    state__Ins_in_PER)  );
        DP_table___Viterbi_states[i][0] = state__Ins_in_PER;                                                                
    }
        

               
    
    //rest of 0-row      
    for (int j = 1; j <= prepared_REFerence.size(); ++j)
    {                                    
        DP_table___Viterbi_log10probabilities[0][j] = 0.00L;  // because we allow it to start farther down on the REF.
        DP_table___Viterbi_states[0][j] = state__Match;       // because we allow it to start farther down on the REF.
    }
    

    
    
        
    
    
    //fill                                                    
    for (int j = 1; j <= prepared_REFerence.size(); ++j) 
    {
        for (int i = 1; i <= prepared_read.size(); ++i)        
        {   		    
	    longdouble best_log10_prob;
	    State_type best_state_transition;
						
            
	    {//init to match						    
		const longdouble marginalized_Match_emission_for_this_position
				    =  calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
					    (prepared_read,
					     i,
					    prepared_read_observed_base_Qualities,
					    prepared_REFerence,
					    j,		     
					    state__Match);
		
		best_log10_prob                
			    =  DP_table___Viterbi_log10probabilities[i-1][j-1]
                                    +  log10(   calculate_P_transition__cond__everything_else
						    (i-1,
						    prepared_REFerence,
						    j-1,		        
						    DP_table___Viterbi_states[i-1][j-1],
						    state__Match)	)
				+  log10(  marginalized_Match_emission_for_this_position  );
											    
		best_state_transition = state__Match;		
	    }//init to match
                                        
            
            
            
	    {//ins	    						    						    
		const longdouble marginalized_Insertion_emission_for_this_position
				    =  calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
					    (prepared_read,
					     i,
					    prepared_read_observed_base_Qualities,
					    prepared_REFerence,
					    j,		     
					    state__Ins_in_PER);					    
						
		const longdouble log10_P_ins
			    =   (DP_table___Viterbi_states[i-1][j] == state__Del_in_PER)  ?
					LOG10_OF_ZERO
					:                                        
					(   DP_table___Viterbi_log10probabilities[i-1][j]
						+  log10(   calculate_P_transition__cond__everything_else
								(i-1,
								prepared_REFerence,
								j,		     
								DP_table___Viterbi_states[i-1][j],
								state__Ins_in_PER)	)
						+  log10(  marginalized_Insertion_emission_for_this_position  )    );
					    
		if (log10_P_ins > best_log10_prob)
		{
		    best_state_transition = state__Ins_in_PER;
		    best_log10_prob = log10_P_ins;                                    
		}           					
	    }//ins
                                    
                                            
                                            
            {//del                                                                                             
		const longdouble log10p_del_here
			    =   (DP_table___Viterbi_states[i][j-1] == state__Ins_in_PER)  ?
					LOG10_OF_ZERO
					:                                        
					(   DP_table___Viterbi_log10probabilities[i][j-1]
						+  log10(   calculate_P_transition__cond__everything_else
								(i,
								prepared_REFerence,
								j-1,		     
								DP_table___Viterbi_states[i][j-1],
								state__Del_in_PER)	)	);
								
		if (log10p_del_here > best_log10_prob)
		{
		    best_log10_prob = log10p_del_here;
		    best_state_transition = state__Del_in_PER;                
		}
	    }//del
                                                                                                                                        
            
            //save
            DP_table___Viterbi_states[i][j] = best_state_transition;
            DP_table___Viterbi_log10probabilities[i][j] = best_log10_prob;                                                    
            
        }// i       
    }//j
                

} // perform_sum_Forward__on__prepared_read_and_prepared_REF_____for_Viterbi_approximation___ie__gloCAL










































real PER_aligner::PRECISELY_align_single_mate_of_PER_to_PRECISE_reference_region
                            (const BamTools::BamAlignment &the_mate,
                             const bool &is_stored_as_first_mate,
                             const std::string &m_reference_PRECISE,
                             const bool &orientation_of_PER_and_MiniRegion_agree,
                             const bool &complement_PER_to_get_to_MiniRegion,
                             
                             type_string__string &PER_and_hybrid_algmt,
                             const bool &return_Viterbi_alignment    )    const
{
        
    
    std::string  prepared_read;
    std::string  prepared_REFerence;
    bool  is_mate_A;
    
    type_vector_int prepared_read_observed_base_Qualities;
        
    bool prepared_REFerence_is_complement_of_original_REFerence;    

    appropriately_orient_mate_and_reference___and___complement_mate_and_determine_strandedness_and_ordering
                                                (the_mate,
                                                    is_stored_as_first_mate,                            
                                                    orientation_of_PER_and_MiniRegion_agree,
                                                    complement_PER_to_get_to_MiniRegion,                            
                                                    m_reference_PRECISE, 
                                                    prepared_read,
                                                    prepared_read_observed_base_Qualities,
                                                    is_mate_A,                            
                                                    prepared_REFerence,
                                                 prepared_REFerence_is_complement_of_original_REFerence); 
						 
    const int prepared_read_length = prepared_read.size();
    const int prepared_REFerence_length = m_reference_PRECISE.size();   						 
    
                 

                                                    
    longdouble         **DP_table___Viterbi_log10probabilities;
    State_type         **DP_table___Viterbi_states;
    
    longdouble         **DP_table____M, **DP_table____I, **DP_table____D;
            
                                                                        
    allocate_and_initialize_2d_array_to_value<longdouble>(DP_table____M, 
                                                                        prepared_read_length+1,
                                                                        prepared_REFerence_length + 1,
                                                                        0.00L );  
    allocate_and_initialize_2d_array_to_value<longdouble>(DP_table____I, 
                                                                        prepared_read_length+1,
                                                                        prepared_REFerence_length + 1,
                                                                        0.00L );       
    allocate_and_initialize_2d_array_to_value<longdouble>(DP_table____D, 
                                                                        prepared_read_length+1,
                                                                        prepared_REFerence_length + 1,
                                                                        0.00L );                                                                               
            
                                                                        

    allocate_and_initialize_2d_array_to_value<longdouble>(DP_table___Viterbi_log10probabilities, 
                                                                        prepared_read_length+1,
                                                                        prepared_REFerence_length + 1,
                                                                        LOG10_OF_ZERO );
    allocate_and_initialize_2d_array_to_value<State_type>(DP_table___Viterbi_states,
                                                                        prepared_read_length+1,
                                                                        prepared_REFerence_length + 1,
                                                                        state__NULL);                                                                            

                                                                        
                                                                            
    perform_sum_Forward__on__prepared_read_and_prepared_REF_____PRECISE___ie___GLOBAL
                                    (DP_table____M,
                                    DP_table____D,
                                    DP_table____I,

                                    DP_table___Viterbi_log10probabilities,
                                    DP_table___Viterbi_states,
                                    
                                    prepared_REFerence,
                                        
                                    prepared_read,                                    
                                    prepared_read_observed_base_Qualities,
                                    
                                    return_Viterbi_alignment);                              
        
       
        
        
    if (return_Viterbi_alignment)   
    {
        int dummy_y_0;        
        do_Viterbi_traceback___and__get_lower_endpoint_on_REFerence
                                (DP_table___Viterbi_states,                                
                                prepared_REFerence,                                    
                                prepared_read,                                        
                                true,      
				 prepared_REFerence.size(),
                                dummy_y_0,
                                 true,
                                PER_and_hybrid_algmt); 
                                
        if (prepared_REFerence_is_complement_of_original_REFerence)                        
        {
            PER_and_hybrid_algmt.first = get_complement_of_sequence( PER_and_hybrid_algmt.first );
            PER_and_hybrid_algmt.second = get_complement_of_sequence( PER_and_hybrid_algmt.second );        
        }
        
        if (!is_mate_A)
        {
            PER_and_hybrid_algmt.first = get_reverse_of_sequence( PER_and_hybrid_algmt.first );
            PER_and_hybrid_algmt.second = get_reverse_of_sequence( PER_and_hybrid_algmt.second );        
        }
    }
        

        
        
    
    
    
    real marginal_probability(  DP_table____M[prepared_read_length][prepared_REFerence_length]  );
    marginal_probability    +=  DP_table____I[prepared_read_length][prepared_REFerence_length];
    marginal_probability    +=  DP_table____D[prepared_read_length][prepared_REFerence_length];        
    
    
    //delete                        
    delete_2d_array<longdouble>(DP_table____M, (uint)(prepared_read_length+1));
    delete_2d_array<longdouble>(DP_table____I, (uint)(prepared_read_length+1));
    delete_2d_array<longdouble>(DP_table____D, (uint)(prepared_read_length+1));        
    
    delete_2d_array<longdouble>(DP_table___Viterbi_log10probabilities, (uint)(prepared_read_length+1));
    delete_2d_array<State_type>(DP_table___Viterbi_states, (uint)(prepared_read_length+1));    
    
        
                                        
    return    marginal_probability;                                            
        
                      
}  // PRECISELY_align_single_mate_of_PER_to_PRECISE_reference_region
























































































  
  
  
  
//return probability of Viterbi alignment (NOT MAP, i.e. not posterior);
void PER_aligner::get_Viterbi_alignment_and_endpoints_of_reference_in_Viterbi_alignment
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0,
                            type_BI__BI &subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                            const bool &return_Viterbi_alignment_sequence)    const
{
    
    std::string prepared_read___true_A;
    type_vector_int prepared_read_observed_nucleotide_Quality___true_A;
    bool naive_mate_A__is_prepared_mate___true_A;
    std::string prepared_REFerence___true_A;
    bool prepared_REFerence_is_complement_of_original_REFerence___true_A;

    
    appropriately_orient_mate_and_reference___and___complement_mate_and_determine_strandedness_and_ordering
                    (mate_A,
                     true,
                     orientation_of_PER_and_MiniRegion_agree,
                     complement_PER_to_get_to_MiniRegion,
                     m_reference,
                     prepared_read___true_A,
                     prepared_read_observed_nucleotide_Quality___true_A,
                     naive_mate_A__is_prepared_mate___true_A,
                     prepared_REFerence___true_A,
                     prepared_REFerence_is_complement_of_original_REFerence___true_A);  
		     
    

    std::string prepared_read___true_B;
    type_vector_int prepared_read_observed_nucleotide_Quality___true_B;
    bool naive_mate_B__is_prepared_mate___true_A;
    std::string prepared_REFerence___true_B;
    bool prepared_REFerence_is_complement_of_original_REFerence___true_B;
    
    
    
    appropriately_orient_mate_and_reference___and___complement_mate_and_determine_strandedness_and_ordering
                    (mate_B,
                     false,
                     orientation_of_PER_and_MiniRegion_agree,
                     complement_PER_to_get_to_MiniRegion,
                     m_reference,
                     prepared_read___true_B,
                     prepared_read_observed_nucleotide_Quality___true_B,
                     naive_mate_B__is_prepared_mate___true_A,
                     prepared_REFerence___true_B,
                     prepared_REFerence_is_complement_of_original_REFerence___true_B); 
                     
    if (  !naive_mate_A__is_prepared_mate___true_A  )
    {
        std::swap<std::string>(prepared_read___true_A, 
			       prepared_read___true_B);
        std::swap<type_vector_int>(prepared_read_observed_nucleotide_Quality___true_A,
						 prepared_read_observed_nucleotide_Quality___true_B);
        std::swap<std::string>(prepared_REFerence___true_A,
			       prepared_REFerence___true_B); 
        std::swap<bool>(prepared_REFerence_is_complement_of_original_REFerence___true_A,
			prepared_REFerence_is_complement_of_original_REFerence___true_B);
    }
               
              

            
    // mate A:    
    longdouble         **DP_table_A__Viterbi_log10probabilities;
    State_type         **DP_table_A__Viterbi_states;
               
    allocate_and_initialize_2d_array_to_value<longdouble>(DP_table_A__Viterbi_log10probabilities, 
                                                                        prepared_read___true_A.size()+1,
                                                                        m_reference.size()+1,
                                                                        LOG10_OF_ZERO );
    allocate_and_initialize_2d_array_to_value<State_type>(DP_table_A__Viterbi_states,
                                                                        prepared_read___true_A.size()+1,
                                                                        m_reference.size()+1,
                                                                        state__NULL);  
                                            
                                                                        
    {// Viterbi approx                                                                   
        perform_sum_Forward__on__prepared_read_and_prepared_REF_____for_Viterbi_approximation___ie__gloCAL
                                    (DP_table_A__Viterbi_log10probabilities,
                                     DP_table_A__Viterbi_states,
                                    
                                    prepared_REFerence___true_A,
                                        
                                    prepared_read___true_A,
                                    prepared_read_observed_nucleotide_Quality___true_A);
    }// Viterbi approx                 
                                              
                                              
                                              
                                              
                                              
                                              
                                   
                                   
                                   
                                   
    // mate B:    
    longdouble         **DP_table_B__Viterbi_log10probabilities;
    State_type         **DP_table_B__Viterbi_states;
               
    allocate_and_initialize_2d_array_to_value<longdouble>(DP_table_B__Viterbi_log10probabilities, 
                                                                        prepared_read___true_B.size()+1,
                                                                        m_reference.size()+1,
                                                                        LOG10_OF_ZERO );
    allocate_and_initialize_2d_array_to_value<State_type>(DP_table_B__Viterbi_states,
                                                                        prepared_read___true_B.size()+1,
                                                                        m_reference.size()+1,
                                                                        state__NULL);  
                                            
                                                                        
    {// Viterbi approx                          
        perform_sum_Forward__on__prepared_read_and_prepared_REF_____for_Viterbi_approximation___ie__gloCAL
                                    (DP_table_B__Viterbi_log10probabilities,
                                     DP_table_B__Viterbi_states,
                                    
                                    prepared_REFerence___true_B,
                                        
                                    prepared_read___true_B,
                                    prepared_read_observed_nucleotide_Quality___true_B); 
    }// Viterbi approx                                                 
                                              
                                              
         


            
            
            
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
            
        
        
        
    //  link the DP of mate A to the DP of mate B.
    
    longdouble log10likelihood_Viterbi_algmt = LOG10_OF_ZERO*2;
    int ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A = -1;    
    int ref_index_of_Viterbi_leftmost_nucleotide_of_mate_B___on_same_REF_as_A___NOT_B = -1;   

    
    int length_of_hidden_advance;
    
            

				
    for (int j_REF__base1 = m_reference.size();  j_REF__base1 > 0; --j_REF__base1)
    {   
        // j_REF__base1 is the leftmost Ref index base 1 of the alignment of read B to Ref.
        
        //best jump
        longdouble log10p_best_jump_from_mate_A = LOG10_OF_ZERO*2;
        int best_ref_index_of_rightmost_mate_A_nucleotide__base1 = -1;                                
        {       
            int  lower_bound   =   std::min<int>(   std::max<int>(1, j_REF__base1 - max_advance_size - 1),
                                                          m_reference.size()     );
                                                        
            int  upper_bound  =   std::min<int>(   std::max<int>(1, j_REF__base1 - min_advance_size - 1),
                                                          m_reference.size()     );                                               
            
	    assert(lower_bound < upper_bound);
//             if (upper_bound - prepared_read___true_A.size() - 1 > j_REF__base1 )
//                 upper_bound = j_REF__base1 + prepared_read___true_A.size() - 1;   // so that mate B really is to right of mate A
//             
//             if (lower_bound - prepared_read___true_A.size() - 1 > j_REF__base1 )
//                 lower_bound = j_REF__base1 + prepared_read___true_A.size() - 1;   // so that mate B really is to right of mate A


            for (int k = lower_bound;   k <= upper_bound;  ++k) 
	    {
                if (DP_table_A__Viterbi_log10probabilities[prepared_read___true_A.size()][k]  >  log10p_best_jump_from_mate_A)
                {
                    log10p_best_jump_from_mate_A = DP_table_A__Viterbi_log10probabilities[prepared_read___true_A.size()][k];
                    best_ref_index_of_rightmost_mate_A_nucleotide__base1 = k;                
                }                               
	    }
        }   
                   
                   
                   
        //add best jump to Viterbi alignment                              
        DP_table_B__Viterbi_log10probabilities[prepared_read___true_B.size()][ m_reference.size() -  (j_REF__base1 - 1) ] += log10p_best_jump_from_mate_A; 
                                                    //since mate B always has the REF reversed so that forward/backward operations can use same function.
                
                                                
        
        // get best (Viterbi) algmt
        if (DP_table_B__Viterbi_log10probabilities[prepared_read___true_B.size()][ m_reference.size() -  (j_REF__base1 - 1) ] > log10likelihood_Viterbi_algmt)
        {
            log10likelihood_Viterbi_algmt = DP_table_B__Viterbi_log10probabilities[prepared_read___true_B.size()][ m_reference.size() -  (j_REF__base1 - 1) ];
            ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A = best_ref_index_of_rightmost_mate_A_nucleotide__base1;
            ref_index_of_Viterbi_leftmost_nucleotide_of_mate_B___on_same_REF_as_A___NOT_B = j_REF__base1;
        }
        
    }//end        
                      


    length_of_hidden_advance
		    = ref_index_of_Viterbi_leftmost_nucleotide_of_mate_B___on_same_REF_as_A___NOT_B
			- ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A - 1;


    
  
			
//     for (uint z=0; z< prepared_REFerence___true_A.size(); ++z)
//     {  std:: cerr << "\n\tA[end][" << z << "] = " << DP_table_A__Viterbi_log10probabilities[prepared_read___true_A.size()][z]
// 		    << ",  B[end][" << z << "] = " << DP_table_B__Viterbi_log10probabilities[prepared_read___true_B.size()][z];  }
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
        

        
        
        
    
    //Viterbi traceback:     
    
    int y_a_0 = -1;
    const int y_a_1 = ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A;
        
    { // A traceback        
	type_string__string PER_and_hybrid_algmt___A;
	
	
	
	do_Viterbi_traceback___and__get_lower_endpoint_on_REFerence
			(DP_table_A__Viterbi_states,
			prepared_REFerence___true_A,
			prepared_read___true_A,
			false,
			ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A,
			y_a_0,
			return_Viterbi_alignment_sequence,
			PER_and_hybrid_algmt___A  );
	
			
	if (return_Viterbi_alignment_sequence)
	{  			
	    if (prepared_REFerence_is_complement_of_original_REFerence___true_A)
	    {
		PER_and_hybrid_algmt.first = get_complement_of_sequence___UNSAFE( PER_and_hybrid_algmt___A.first );
		PER_and_hybrid_algmt.second = get_complement_of_sequence___UNSAFE( PER_and_hybrid_algmt___A.second );        
	    }
	    else
	    {
		PER_and_hybrid_algmt = PER_and_hybrid_algmt___A;	    
	    }
	    
	    // append hidden advance.                                
	    if (length_of_hidden_advance > 0)
	    {
		PER_and_hybrid_algmt.first.append( std::string(length_of_hidden_advance, 'x') );        
		PER_and_hybrid_algmt.second.append(   m_reference.substr(ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A + 1 - 1, 
											length_of_hidden_advance)   );                                                        
	    } 
	    
	}     			                    
    } // A traceback
                    
                     
    
    
    
    
    
    
    const int y_b_0 = ref_index_of_Viterbi_leftmost_nucleotide_of_mate_B___on_same_REF_as_A___NOT_B;
    int y_b_1 = -1;
        
    {// B traceback
	type_string__string PER_and_hybrid_algmt___B;
	
	
	do_Viterbi_traceback___and__get_lower_endpoint_on_REFerence
			(DP_table_B__Viterbi_states,
			prepared_REFerence___true_B,
			prepared_read___true_B,
			false,
			m_reference.size() - (ref_index_of_Viterbi_leftmost_nucleotide_of_mate_B___on_same_REF_as_A___NOT_B - 1),
			y_b_1,
			return_Viterbi_alignment_sequence,
			PER_and_hybrid_algmt___B  );                     
			
	// since the REF for mate B is ALWAYS reversed so that all alignments have a "standardized" direction, then we need to convert y_b_1 to the original ref value.
// 	y_b_1 =   y_b_0 + (prepared_read___true_B.size() - 1) - (y_b_1 - 1);
	y_b_1 =   m_reference.size() - (y_b_1 - 1);
			
			
	if (return_Viterbi_alignment_sequence)
	{         
	    if (prepared_REFerence_is_complement_of_original_REFerence___true_B)
	    {
		PER_and_hybrid_algmt___B.first = get_complement_of_sequence___UNSAFE( PER_and_hybrid_algmt___B.first );
		PER_and_hybrid_algmt___B.second = get_complement_of_sequence___UNSAFE( PER_and_hybrid_algmt___B.second );        
	    }
	    
	    //since mate B and its reference are always reversed for ease of forward/backward functions...                 
	    PER_and_hybrid_algmt___B.first = get_reverse_of_sequence( PER_and_hybrid_algmt___B.first );
	    PER_and_hybrid_algmt___B.second = get_reverse_of_sequence( PER_and_hybrid_algmt___B.second );
	    
	    
	    //recall:   
		    // length_of_hidden_advance = ref_index_of_Viterbi_leftmost_nucleotide_of_mate_B___on_same_REF_as_A___NOT_B - ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A - 1;
	    
	    if (length_of_hidden_advance  >=  0)
	    {
		PER_and_hybrid_algmt.first.append( PER_and_hybrid_algmt___B.first );
		PER_and_hybrid_algmt.second.append( PER_and_hybrid_algmt___B.second );
	    }
	    else  // < 0
	    {        
		//thus:       ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A   >=   ref_index_of_Viterbi_leftmost_nucleotide_of_mate_B___on_same_REF_as_A___NOT_B                                                                
		
		const int length_of_algmt__B = PER_and_hybrid_algmt___B.second.size();
		
		
		int algmt_start_adjusted = 0;
		
		//we don't want any part of the alignment before mate_A alignment is done...
		{
		    int ref_index_ctr = ref_index_of_Viterbi_leftmost_nucleotide_of_mate_B___on_same_REF_as_A___NOT_B;
		    while (  algmt_start_adjusted < length_of_algmt__B   and   ref_index_ctr  <=  ref_index_of_Viterbi_rightmost_nucleotide_of_mate_A)
		    {
			if (  PER_and_hybrid_algmt___B.second.at(algmt_start_adjusted) != '-'  )
			    ++ref_index_ctr;
			
			++algmt_start_adjusted;
		    }
		    
		    while ( algmt_start_adjusted  <  length_of_algmt__B    and    PER_and_hybrid_algmt___B.second.at(algmt_start_adjusted) == '-' )
			++algmt_start_adjusted;                        
		}
			
		    
		if (algmt_start_adjusted  <  length_of_algmt__B)     
		{
		    PER_and_hybrid_algmt.first.append( PER_and_hybrid_algmt___B.first.substr( algmt_start_adjusted  ) );
		    PER_and_hybrid_algmt.second.append( PER_and_hybrid_algmt___B.second.substr(  algmt_start_adjusted  ) );                 
		}                                                            
	    }  //  < 0
	    
	    
	    
	} // return_Viterbi_alignment_sequence
			
    
    } // B traceback
    
  
  
  
    
                                          
             
    
    assert( y_a_0 > 0 );
    assert( y_a_1 > 0 );
    assert( y_b_0 > 0 );
    assert( y_b_1 > 0 );
    
    
    
    
    //convert to base 0:
    subset_of_region_that_is_aligned_to_mates__A__B___base_0.first.set(    (uint)y_a_0 - 1,   (uint)y_a_1 - 1     );
                                                  
    subset_of_region_that_is_aligned_to_mates__A__B___base_0.second.set(   (uint)y_b_0 - 1,   (uint)y_b_1 - 1      );   
    
    
//     subset_of_region_that_is_aligned_to_mates__A__B___base_0.first.set(   subset_of_region_that_is_aligned_to_mates__A__B___base_0.first.lower() - 1,
//                                                                           subset_of_region_that_is_aligned_to_mates__A__B___base_0.first.upper() - 1  );
//                                                                           
//     subset_of_region_that_is_aligned_to_mates__A__B___base_0.second.set(   subset_of_region_that_is_aligned_to_mates__A__B___base_0.second.lower() - 1,
//                                                                            subset_of_region_that_is_aligned_to_mates__A__B___base_0.second.upper() - 1  );
                                                                          
                                                                          
    subset_of_region_that_is_aligned___base_0  =  BOOST_hull(   subset_of_region_that_is_aligned_to_mates__A__B___base_0.first,
                                                                subset_of_region_that_is_aligned_to_mates__A__B___base_0.second    );
                                                                
                                                                
    //back to original mateA, mate B
    if ( !orientation_of_PER_and_MiniRegion_agree  )
    {
        std::swap<BOOST_Interval>(  subset_of_region_that_is_aligned_to_mates__A__B___base_0.first,
                                    subset_of_region_that_is_aligned_to_mates__A__B___base_0.second   );
    }
                                                   
                                                   
    if (        BOOST_empty(subset_of_region_that_is_aligned___base_0)  
           or   BOOST_empty(subset_of_region_that_is_aligned_to_mates__A__B___base_0.first)
           or   BOOST_empty(subset_of_region_that_is_aligned_to_mates__A__B___base_0.second)   )    
    {
        std::stringstream error_strm;
        print_line_of_markers("ERROR! ", &error_strm);
        print_line_of_markers("(", &error_strm);
        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
        
        error_strm << "\n\nERROR!   in  PER_aligner___base \"align_PER_to_region\":    subset_of_region_that_is_aligned___base_0  is empty!!!!!\n\n";
         
        
        error_strm << "\n\nprepared_read___true_A:\n" << prepared_read___true_A
                    << "\nprepared_read___true_B:\n" << prepared_read___true_B 
                    << "\n\nm_reference:\n" << m_reference
                    << "\n\n";
                    
        error_strm  << "\n\nThese are base 1:\ny_a_0 = " << y_a_0
                    << "\ny_a_1 = " << y_a_1
                    << "\ny_b_0 = " << y_b_0
                    << "\ny_b_1 = " << y_b_1
                    << "\n\n";
                    
                    
        error_strm << "\nPER_and_hybrid_algmt:\n" << PER_and_hybrid_algmt.first << "\n\n" << PER_and_hybrid_algmt.second << "\n\n";                      
                    
                                        
        print_1d_array<longdouble>(DP_table_B__Viterbi_log10probabilities[1],
                            m_reference.size(),
                            "DP_table_B__Viterbi_log10probabilities[1]");             
                                                            
        
        error_strm << "\nsubset_of_region_that_is_aligned = " << subset_of_region_that_is_aligned___base_0 << "\n";
        
        error_strm << "subset_of_region_that_is_aligned_to_mates__A__B___base_0.first:  "  
                            << subset_of_region_that_is_aligned_to_mates__A__B___base_0.first
                    << "\nsubset_of_region_that_is_aligned_to_mates__A__B.second :   " 
                                    << subset_of_region_that_is_aligned_to_mates__A__B___base_0.second << "\n\n";
                                    
//         my_PER_ptr->print_this_PER( &error_strm );                            
        
        print_line_of_markers(")", &error_strm);
        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
    
    }                                                   
                                                   
                                                   
                            
                            
    
    //delete
    delete_2d_array<longdouble>(DP_table_A__Viterbi_log10probabilities, (uint)(prepared_read___true_A.size()+1));
    delete_2d_array<State_type>(DP_table_A__Viterbi_states, (uint)(prepared_read___true_A.size()+1));        

    
    //delete
    delete_2d_array<longdouble>(DP_table_B__Viterbi_log10probabilities, (uint)(prepared_read___true_B.size()+1) );
    delete_2d_array<State_type>(DP_table_B__Viterbi_states, (uint)(prepared_read___true_B.size()+1) );                                                   
    

}  // get_Viterbi_alignment_and_endpoints_of_reference_in_Viterbi_alignment













































//return probability of Viterbi alignment (NOT MAP, i.e. not posterior);
real PER_aligner::calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference____IMPL
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0,
                            type_BI__BI &subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                             const bool &return_Viterbi_alignment_sequence)    const
{
    get_Viterbi_alignment_and_endpoints_of_reference_in_Viterbi_alignment
                            (  mate_A,
                              mate_B,
                             m_reference,
                              orientation_of_PER_and_MiniRegion_agree,
                              complement_PER_to_get_to_MiniRegion,
                            PER_and_hybrid_algmt,
                            subset_of_region_that_is_aligned___base_0,
                            subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                               return_Viterbi_alignment_sequence);
                          
                            
    type_string__string dummy_algmt_1, dummy_algmt_2;          
    const real marginal_probability_of_First_mate_aligning_to_precise_subregion(
                                PRECISELY_align_single_mate_of_PER_to_PRECISE_reference_region
                                        (mate_A,
                                        true,//is_stored_as_first_mate,
                                        m_reference.substr(  subset_of_region_that_is_aligned_to_mates__A__B___base_0.first.lower(),
                                                             BOOST_width_inclusive( subset_of_region_that_is_aligned_to_mates__A__B___base_0.first )  ),
                                        orientation_of_PER_and_MiniRegion_agree,
                                        complement_PER_to_get_to_MiniRegion,
                                        
                                        dummy_algmt_1,
                                        false)       );
					
    const real marginal_probability_of_Second_mate_aligning_to_precise_subregion(
                                PRECISELY_align_single_mate_of_PER_to_PRECISE_reference_region
                                        (mate_B,
                                        false,//is_stored_as_first_mate,
                                        m_reference.substr(  subset_of_region_that_is_aligned_to_mates__A__B___base_0.second.lower(),
                                                             BOOST_width_inclusive( subset_of_region_that_is_aligned_to_mates__A__B___base_0.second )  ),
                                        orientation_of_PER_and_MiniRegion_agree,
                                        complement_PER_to_get_to_MiniRegion,
                                         
                                        dummy_algmt_2,
                                         false)       );   
               
                                        
    return  (     marginal_probability_of_First_mate_aligning_to_precise_subregion
                * marginal_probability_of_Second_mate_aligning_to_precise_subregion    );                                                                                                                            

}//calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference____IMPL































real PER_aligner::calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0,
                            type_BI__BI &subset_of_region_that_is_aligned_to_mates__A__B___base_0)    const
{

    return  calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference____IMPL
                            (mate_A,
                            mate_B,
                            m_reference,
                            orientation_of_PER_and_MiniRegion_agree,
                            complement_PER_to_get_to_MiniRegion,
                            PER_and_hybrid_algmt,
                            subset_of_region_that_is_aligned___base_0,
                            subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                             true);



} // calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference







real PER_aligner::calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,//const std::string &in_read_A,
                            const BamTools::BamAlignment &mate_B,//const std::string &in_read_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0)    const
{
    type_BI__BI subset_of_region_that_is_aligned_to_mates__A__B___base_0;
    
    return  calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference____IMPL
                            (mate_A,//const std::string &in_read_A,
                            mate_B,//const std::string &in_read_B,
                            m_reference,
                            orientation_of_PER_and_MiniRegion_agree,
                            complement_PER_to_get_to_MiniRegion,
                            PER_and_hybrid_algmt,
                            subset_of_region_that_is_aligned___base_0,
                            subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                             true);


}//calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference





real PER_aligner::calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,//const std::string &in_read_A,
                            const BamTools::BamAlignment &mate_B,//const std::string &in_read_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            type_string__string &PER_and_hybrid_algmt)    const
{
    
    BOOST_Interval subset_of_region_that_is_aligned___base_0;
    type_BI__BI subset_of_region_that_is_aligned_to_mates__A__B___base_0;
    
    return  calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference____IMPL
                            (mate_A,//const std::string &in_read_A,
                            mate_B,//const std::string &in_read_B,
                            m_reference,
                            orientation_of_PER_and_MiniRegion_agree,
                            complement_PER_to_get_to_MiniRegion,
                            PER_and_hybrid_algmt,
                            subset_of_region_that_is_aligned___base_0,
                            subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                             true);    




}//calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference                            
                            
                            
                            
                            
                            
                            
                            
real PER_aligner::calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,//const std::string &in_read_A,
                            const BamTools::BamAlignment &mate_B,//const std::string &in_read_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion)    const    
{

    type_string__string PER_and_hybrid_algmt;
    BOOST_Interval subset_of_region_that_is_aligned___base_0;
    type_BI__BI subset_of_region_that_is_aligned_to_mates__A__B___base_0;
    
    return  calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference____IMPL
                            (mate_A,//const std::string &in_read_A,
                            mate_B,//const std::string &in_read_B,
                            m_reference,
                            orientation_of_PER_and_MiniRegion_agree,
                            complement_PER_to_get_to_MiniRegion,
                            PER_and_hybrid_algmt,
                            subset_of_region_that_is_aligned___base_0,
                            subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                             false);      

}//calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference                            






real  PER_aligner::calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
                            (const BamTools::BamAlignment &mate_A,//const std::string &in_read_A,
                            const BamTools::BamAlignment &mate_B,//const std::string &in_read_B,
                            const std::string &m_reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
                            BOOST_Interval &subset_of_region_that_is_aligned___base_0)    const
{
    
    type_string__string PER_and_hybrid_algmt;
    type_BI__BI subset_of_region_that_is_aligned_to_mates__A__B___base_0;
    
    return  calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference____IMPL
                            (mate_A,//const std::string &in_read_A,
                            mate_B,//const std::string &in_read_B,
                            m_reference,
                            orientation_of_PER_and_MiniRegion_agree,
                            complement_PER_to_get_to_MiniRegion,
                            PER_and_hybrid_algmt,
                            subset_of_region_that_is_aligned___base_0,
                            subset_of_region_that_is_aligned_to_mates__A__B___base_0,
                             false);


}//calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference



































































































































































	












































































longdouble  PER_aligner::calculate_P_transition__cond__everything_else
		    (const uint &cycle,
		     const std::string &prepared_reference_sequence,
		     const uint &reference_position__base_1,		     		    
		     const State_type &state_FROM,
		     const State_type &state_TO  )   const
{
    const uint homopolymer_length = identify_homopolymer_length_at_given_position(prepared_reference_sequence, reference_position__base_1);
    
    return 
	((homopolymer_length < 3)
		*(    (state_TO == 0)*TRANSITIONS___BASE__P_MATCH 
		    + (state_TO > 0)*(   (state_FROM == state__Match)*TRANSITIONS___BASE__P_GAP_OPEN 
					 + (state_FROM != state__Match)*TRANSITIONS___BASE__P_GAP_EXT))
	+ (homopolymer_length >= 3)
		*(    (state_TO == 0)*TRANSITIONS___BASE__P_HOMOPOLYMER_MATCH 
		    + (state_TO > 0)*(   (state_FROM == state__Match)*TRANSITIONS___BASE__P_HOMOPOLYMER_GAP_OPEN 
					 + (state_FROM != state__Match)*TRANSITIONS___BASE__P_HOMOPOLYMER_GAP_EXT)));    
    
    
}//calculate_P_transition__cond__everything_else










longdouble  PER_aligner::calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q
		    (const std::string &prepared_read_sequence,
		     const uint &position_along_read___base_1,
		     const type_vector_int &prepared_sequence_base_Qualities,
		     const std::string &prepared_reference_sequence,
		     const uint &reference_position__base_1,    
		     const State_type &state_of_emission)   const
{    
    if (reference_position__base_1 == 0) // will NOT cause warp divergence
    { 
	return 0.25;
    }
    else
    {//non-trivial
    
	const int quality_value = prepared_sequence_base_Qualities.at(position_along_read___base_1-1);
	
	if (quality_value <= BASE_QUAL_CUTOFF)
	{  
	    return 0.25; 
	}
	else
	{//good qual
	    const bool observed_base_and_ref_base_agree =  (prepared_read_sequence.at(position_along_read___base_1-1) == prepared_reference_sequence.at(reference_position__base_1-1));
	
	    double addl_error_prob  = exp10(((double)quality_value * EMISSIONS__QUAL_DECREASE_FACTOR)/(-10));// a little extra, to be conservative		
	
	    //the following will NOT cause warp divergence:
	    if (reference_position__base_1 == 1)
	    {
		addl_error_prob += EMISSIONS__EARLY_ERR; // beginning 2 bases of read is notoriously bad
		
		addl_error_prob = std::min<longdouble>(addl_error_prob, (EMISSIONS___BASE_P_CORRECT-0.25));
		return   observed_base_and_ref_base_agree   ? 
				    (EMISSIONS___BASE_P_CORRECT - addl_error_prob)   :   (EMISSIONS___BASE_P_ERR + addl_error_prob)/3;
	    }
	    else if (reference_position__base_1 == 2)
	    {
		addl_error_prob += (prepared_reference_sequence.at(reference_position__base_1-2) == 'G' or prepared_reference_sequence.at(reference_position__base_1-2) == 'T')*EMISSIONS__GX_or_TX_ERR; // "GX" and "TX"
		
		addl_error_prob = std::min<longdouble>(addl_error_prob, (EMISSIONS___BASE_P_CORRECT-0.25));
		return   observed_base_and_ref_base_agree   ? 
				    (EMISSIONS___BASE_P_CORRECT - addl_error_prob)   :   (EMISSIONS___BASE_P_ERR + addl_error_prob)/3;
	    }
	    else //triplet
	    {
		const uint homopolymer_length = identify_homopolymer_length_at_given_position(prepared_reference_sequence, reference_position__base_1);
		const uint mate_length = prepared_read_sequence.size();				
		
		addl_error_prob += 
			(prepared_reference_sequence.at(reference_position__base_1-3) == 'G' and prepared_reference_sequence.at(reference_position__base_1-2) == 'G')*EMISSIONS__GGX_ERR  // GGX
			+ (prepared_reference_sequence.at(reference_position__base_1-3) == 'A' and prepared_reference_sequence.at(reference_position__base_1-2) == 'G')*EMISSIONS__AGX_ERR  // AGX
			+ (prepared_reference_sequence.at(reference_position__base_1-3) == prepared_reference_sequence.at(reference_position__base_1-2))*EMISSIONS__STUTTER_ERR  // BBX - stutter
			+ (position_along_read___base_1 > 0.25*mate_length)*EMISSIONS__CYCLE_ERR   //cycle errors
			+ (position_along_read___base_1 > 0.50*mate_length)*EMISSIONS__CYCLE_ERR   //cycle errors
			+ (position_along_read___base_1 > 0.75*mate_length)*EMISSIONS__CYCLE_ERR   //cycle errors
			+ (homopolymer_length >= 3)*EMISSIONS__HOMOPOLYMER_ERR
			+ (homopolymer_length >= 6)*EMISSIONS__HOMOPOLYMER_ERR;  //homopolymer
			
		addl_error_prob = std::min<longdouble>(addl_error_prob, (EMISSIONS___BASE_P_CORRECT-0.25));
		return   observed_base_and_ref_base_agree   ? 
				    (EMISSIONS___BASE_P_CORRECT - addl_error_prob)   :   (EMISSIONS___BASE_P_ERR + addl_error_prob)/3;
	    } //triplets	  
	}//good qual
	
    }//non-trivial
                       
    
} // calculate__P_observed_base_and_quality__cond__everything_else__CONDITIONAL_LIKELIHOOD_ON_Q


