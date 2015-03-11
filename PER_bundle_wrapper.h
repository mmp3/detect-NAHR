#ifndef PER_BUNDLE_WRAPPER_H_INCLUDED
#define PER_BUNDLE_WRAPPER_H_INCLUDED


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


#include <general_typedefs.h>

#include "globals.h"

#include <CUPA.h>




class MiniRegion;
class Paired_end_read;
class MiniEvent;


namespace  PER_bundle_wrapper__space
{//PER_bundle_wrapper__space


    class PER_bundle_wrapper
    {
	public:	
	    PER_bundle_wrapper() : is_hybrid(false)
	    {
		std::cerr << "\n\n\nERROR!  should never default construct \"PER_bundle_wrapper\"!!!\n\n\n";
		exit(1);		
	    }	
	    
	    
	    //orienting
	    PER_bundle_wrapper( const CUPA::type_list_Paired_read_alignment_bundle::const_iterator  &in_it_gpu_result__nohybrid,	    
				Paired_end_read *const  &in_ptr__PER)
							    :  it_gpu_result__AB(in_it_gpu_result__nohybrid),
								ptr__PER(in_ptr__PER),
								is_hybrid(false)
	    {  }
	    	    
	    
	    //no-hybrid constructor
	    PER_bundle_wrapper( const CUPA::type_list_Paired_read_alignment_bundle::const_iterator  &in_it_gpu_result__nohybrid,	    
				Paired_end_read *const  &in_ptr__PER,			    
				const std::map<uint,MiniRegion>::const_iterator &in_mini_it,
				const bool &in_mini_orientation_agrees_with_MiniEvent_profile)
				    :  it_gpu_result__AB(in_it_gpu_result__nohybrid),
					ptr__PER(in_ptr__PER),
					mini_A_it(in_mini_it),
					is_hybrid(false),
					mini_A_orientation_agrees_with_MiniEvent_profile(in_mini_orientation_agrees_with_MiniEvent_profile)
	    {  }	
	    
	    
	    //hybrid-constructor
	    PER_bundle_wrapper( const CUPA::type_list_Paired_read_alignment_bundle::const_iterator  &in_it_gpu_result__AB,
				const CUPA::type_list_Paired_read_alignment_bundle::const_iterator  &in_it_gpu_result__BA,
				
				Paired_end_read *const  &in_ptr__PER,
				
				const std::map<uint,MiniRegion>::const_iterator &in_mini_A_it,
				const std::map<uint,MiniRegion>::const_iterator &in_mini_B_it,   
				
				const std::map<uint, MiniEvent>::iterator  &in_it_ME,
				
				const BOOST_Interval  &in_profile_brkpt_interval,
				const BOOST_Interval  &in_hybrid_breakpoint__AB,
				const BOOST_Interval  &in_hybrid_breakpoint__BA,
				const BOOST_Interval  &in_absolute_brkpt_interval__on_A,
				const BOOST_Interval  &in_absolute_brkpt_interval__on_B,
				const bool &in_is_NAHR,
				
				const bool &in_mini_A_orientation_agrees_with_MiniEvent_profile,
				const bool &in_mini_B_orientation_agrees_with_MiniEvent_profile)
				    :  it_gpu_result__AB(in_it_gpu_result__AB),
					it_gpu_result__BA(in_it_gpu_result__BA),
					ptr__PER(in_ptr__PER),
					mini_A_it(in_mini_A_it),
					mini_B_it(in_mini_B_it),
					MiniEvent_intersecting_minis_AB(in_it_ME),
					is_hybrid(true),
					profile_brkpt_interval(in_profile_brkpt_interval),
					hybrid_AB__breakpoint(in_hybrid_breakpoint__AB),
					hybrid_BA__breakpoint(in_hybrid_breakpoint__BA),
					absolute_brkpt_interval__on_A(in_absolute_brkpt_interval__on_A),
					absolute_brkpt_interval__on_B(in_absolute_brkpt_interval__on_B),
					is_NAHR(in_is_NAHR),
					mini_A_orientation_agrees_with_MiniEvent_profile(in_mini_A_orientation_agrees_with_MiniEvent_profile),
					mini_B_orientation_agrees_with_MiniEvent_profile(in_mini_B_orientation_agrees_with_MiniEvent_profile)
	    {  }

	    
	    
		
	    
	    ~PER_bundle_wrapper()
	    {  }	    	    	    
		
		
	    //functions:	    
		void process_results
				(const bool &special_save_breakpoints_actually_captured_in_alignment)  const;
				
		void print_this_PER_bundle_wrapper
					(std::stringstream *const &some_ss = NULL)  const;
				
				
	    //variables:
	    
		//the most important one:  the GPU result:
		CUPA::type_list_Paired_read_alignment_bundle::const_iterator  it_gpu_result__AB;
		CUPA::type_list_Paired_read_alignment_bundle::const_iterator  it_gpu_result__BA;
		
		//PER
		Paired_end_read* ptr__PER;
		
	    private:
	    
		//MR
		std::map<uint,MiniRegion>::const_iterator  mini_A_it;
		std::map<uint,MiniRegion>::const_iterator  mini_B_it;//meaningless if is "no-hybrid"
		
		//ME:
		std::map<uint, MiniEvent>::iterator  MiniEvent_intersecting_minis_AB;//meaningless if is "no-hybrid"
		
		//breakpoint:
		const bool is_hybrid;
		BOOST_Interval  profile_brkpt_interval;
		BOOST_Interval  hybrid_AB__breakpoint;  //relative to length of hybrid
		BOOST_Interval  hybrid_BA__breakpoint;
		BOOST_Interval  absolute_brkpt_interval__on_A;
		BOOST_Interval  absolute_brkpt_interval__on_B;
		bool is_NAHR;	    	    
		
		//orientation
		bool mini_A_orientation_agrees_with_MiniEvent_profile;
		bool mini_B_orientation_agrees_with_MiniEvent_profile;	    

    };//PER_bundle_wrapper



    typedef   std::list<PER_bundle_wrapper>   type_list_PER_bundle_wrapper;
    typedef   std::map<std::string, type_list_PER_bundle_wrapper>   type_map_string_to_list_PER_bundle_wrapper;
    
    
    
    
    
    
    
    
    
    
    
    //helper functions: 
    void remove_N_tails_from_read_and_quality
				    (std::string &prepared_read_bases,
				     std::string &prepared_read_qualities);   
				    
    void remove_N_tails_from_read_and_quality
				    (std::string &prepared_read_bases,
				     type_vector_int &prepared_read_qualities);				    
    
    
    mpfr::mpreal  get_likelihood_of_Viterbi_alignment
			(const CUPA::type_list_Paired_read_alignment_bundle::const_iterator &it_some_gpu_result);
		    
    
    void  prepare_PER_mates_for_GPU_alignment
                            (const BamTools::BamAlignment &mate_A,
                            const BamTools::BamAlignment &mate_B,
			     const std::string &reference,
                            const bool &orientation_of_PER_and_MiniRegion_agree,
                            const bool &complement_PER_to_get_to_MiniRegion,
			     std::string &prepared_read_A,
			     std::string &prepared_base_Qualities__A,
			     std::string &prepared_read_B,			     
			     std::string &prepared_base_Qualities__B,
			     std::string &prepared_reference);
				
    CUPA::type_list_Paired_read_alignment_bundle::const_iterator format_and_submit_alignment_work_to_GPU_workload
				(CUPA::GPU_alignment_workload &GPU_work,
				const BamTools::BamAlignment &mate_A,
				const BamTools::BamAlignment &mate_B,
				const std::string &reference,
				const bool &orientation_of_PER_and_MiniRegion_agree,
				const bool &complement_PER_to_get_to_MiniRegion);
				
				
    void  process_orienting_of_an_illdefined_PER
		    (type_map_string_to_list_PER_bundle_wrapper::const_iterator it_some_PER);				
				


}//PER_bundle_wrapper__space


#endif