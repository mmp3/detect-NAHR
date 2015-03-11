#ifndef PAIRED_END_READ_H_INCLUDED
#define PAIRED_END_READ_H_INCLUDED


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



#include <api/BamAlignment.h>



#include <general_typedefs.h>
#include <Sparse_map.h>
#include <Readgroup_statistics.h>



#include "globals.h"
#include "PER_aligner.h"
#include "MiniEvent.h"









class Event;
class Conn_Comp;
class Star_alignment;
class MiniRegion;

namespace  PER_bundle_wrapper__space
{
    class PER_bundle_wrapper;
};


namespace  CUPA
{
    class GPU_alignment_workload;    
};




class Paired_end_read
{
  public:
    //constructors:
    Paired_end_read() : name("NOT_NAMED__BY_MATT!!!"),
                        my_Readgroup_ptr(NULL),
                        haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal(0, 0.00L),
                        haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome(0, 0.00L),
                        haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome(0, 0.00L),
                        MiniID_ctr(0),
                        failed_PER_flag(true)
    { }             
    
    Paired_end_read( const BamTools::BamAlignment &some_Balgmt,
                     const BamTools::BamAlignment &another_Balgmt)
                                                  : name(some_Balgmt.Name),
                                                    haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal(0, 0.00L),
                                                    haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome(0, 0.00L),
                                                    haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome(0, 0.00L),
                                                    MiniID_ctr(0),
                                                    my_PER_aligner(this)
                                                    
                                                                    
    {                        
	init_PER(some_Balgmt, another_Balgmt);                       
    }      
    
    
    

    
    

    //destructor
    ~Paired_end_read()
    { }
    



    
    //variables:
	std::string name;    
	const Readgroup_statistics *my_Readgroup_ptr;        
	type_Balgmt__Balgmt mates;

	type_uint__uint spawning_profile_chromo_and_position;
	
	bool failed_PER_flag;

	std::map<uint, MiniEvent>  my_MiniEvents;
	std::map<uint, MiniRegion>  my_MiniRegions;
	
	type_set_uint set_of_MiniRegions_on__XXXX_chromosome;
	type_set_uint set_of_MiniRegions_on__YYYY_chromosome;
	
	
	type_map_uint_to_real map_affected_MiniID__to__P_nohybrid;
	
	type_map_uint_to_BI  map_MiniID__to__nohybrid_almgt_interval_absolute;
	
	
	type_uint__real  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal;
	
	type_uint__real  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome;
	type_uint__real  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome;	
	
	
	
	uint MiniID_ctr;   //this is only for counting "add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions"
	

	PER_aligner my_PER_aligner;  
	
	
    	
	
	
	
	
	
	
	
    //functions:   	    
	
    void print_this_PER(std::stringstream *const &some_ss = NULL)  const;    
    
    void init_PER(
		const BamTools::BamAlignment &some_Balgmt,
		const BamTools::BamAlignment &another_Balgmt);    
		    
    const MiniRegion *const get_ptr_to_first_MiniRegion();
    
    void orient_and_complement_each_mate_according_to_its_spawning_Event_profile();                                                                                      
			    
    void  align_to_every_possible_hybrid_outcome_of_MiniRegions_and_Events_and_record_probabilities_and_relevant_breakpoints
			    (CUPA::GPU_alignment_workload &GPU_work,
			     std::map<std::string, std::list<PER_bundle_wrapper__space::PER_bundle_wrapper> >  &record_of_PER_algmt_jobs,
			     const bool &consider_GeneConversion_breakpoints);			    
			    
    void sum_over_all_unaffected_LCRs();                                          
    
    
    void filter_haploid_state_vector_through_MiniEvents
                        (Sparse_map& filtered_sparse_haploid_state_vector,
                        type_map_uint_to_set_uint &sparse_map_occurring_MiniEvent__to__MiniRegions_on_LCRs,
                        type_map_uint_to_set_uint sparse_map_occurring_MiniEvent__to__MiniRegions_in_Inbetween_region,
                        const Sparse_map &proposed_sparse_haploid_state_vector)    const;
    
                                       
    type_uint__real calculate_or_get__haploid_sum_P_mPER__cond__sparse_haploid_state_vector_and_breakpoints
                                                  (const bool &haploid_0_or_1, 
                                                   const Sparse_map &sparse_haploid_state_vector,
                                                   const type_map_uint_to_BI &sparse_haploid_breakpoints)  const;    
                                                                                                                                                                                                              
            
    void add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately
                                        (const MiniRegion &some_MR);    
                                        
    void add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately
                                        (const Paired_end_read &another_PER);    
                                        
    void add_list_of_MiniRegions_to_this_PERs_list_of_MiniRegions_appropriately
                                                        (const std::map<uint, MiniRegion> &some_MRs);                                        
    
    
    type_map_uint_to_uint_to_set_uint create_MiniEvents_and_identify_gender_of_MRs
                                                                        (const Event *const &spawning_Event);                              
                                        
    void make_partners_out_of_MiniRegions_for_each_MiniEvent
                       (const type_map_uint_to_uint_to_set_uint &MiniEvent__to__MiniRegions_on_LCRs__to__intersected_profile_positions___PRELIMINARY);    
                                        
                                        
                                        
    void clear_all_MiniRegion_sequences__for_memory_considerations();    
    
    
    

    void  HEURISTIC__determine_which_breakpoints_are_strongly_supported_by_this_PER_for_a_given_Event____using_pseudo_posterior
                                  (const uint &UID_of_MiniEvent,
                                   const longdouble &posterior_threshold,
				   type_map_uint_to_uint &possibly_supported_breakpoints__01,
				   type_map_uint_to_uint &possibly_supported_breakpoints__10)  const;                  
                                  
    

    

    bool align_to_every_possible_hybrid_outcome_AND_DISPLAY
                            (const uint &EV_UID,
                             const BOOST_Interval &desired_breakpoint,
			     const type_haploid_outcome event_outcome,			     
                             Star_alignment  &star_alignment_for_LCR_0,
                             Star_alignment  &star_alignment_for_LCR_1,
                             Star_alignment  &star_alignment_for_hybrid_AB,
                             Star_alignment  &star_alignment_for_hybrid_BA);   
			    
	bool align_to_every_possible_hybrid_outcome_AND_DISPLAY_REGARDLESS_OF_NOHYBRIDS
                            (const uint &EV_UID,
                             const BOOST_Interval &desired_breakpoint,
			     const type_haploid_outcome event_outcome,
                             Star_alignment  &star_alignment_for_LCR_0,
                             Star_alignment  &star_alignment_for_LCR_1,
                             Star_alignment  &star_alignment_for_hybrid_AB,
                             Star_alignment  &star_alignment_for_hybrid_BA);			    
                             
    
    

                            

    void calculate_probability_of_PER_given_hybrid_AND_DISPLAY
                        (const type_map_uint_to_uint &hybrid_AB_algmt_coords,
                         const type_string__string &PER_and_hybrid_AB_algmt,
                         const uint &index_of_hybrid_AB_that_is_first_B_nucleotide_IE_the_breakpoint_varpos,
                         const uint &absolute_genome_position_of_A_mapping_to_index_0_of_hybrid,  //Note, index 0 of hybrid is NOT  necessarily the first index in the algmt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         const bool &A_mini_was_reversed,
                         const uint &chromosome_of_mini_A,
                         const uint &absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid,
                         const bool &B_mini_was_reversed,
                         const uint &chromosome_of_mini_B,
                         std::stringstream &output_str_stm,
                         const bool &test_is_AB,
                         const Event *const &occurring_Event);     
			
			
			
			
    void  submit_PER_orienting_work_to_GPU
			    (const MiniRegion &some_MR,
			     CUPA::GPU_alignment_workload &GPU_work,
			    std::map<std::string, std::list<PER_bundle_wrapper__space::PER_bundle_wrapper> >  &record_of_PER_algmt_jobs);
			    
                         
                         
                         
//     void extend_MiniRegions
//                         (const type_map_uint_to_list_BI__string &preloaded_sequences);                         
                                   
                        
                        
    std::map<uint, MiniRegion> get_MRs()  const;
 
    void set_sequence_for_all_MiniRegions
                        (type_map_uint_to_list_BI__string &EV_uploaded_pieces_of_chromos);          
         
                                                
    void extend_MRs_to_minimum_length();  
    
    uint determine_observed_frag_length()  const;
    
    
    
	 
					 
					 
    void adjust_orientation_and_complementairy_of_PER_mates_according_its_spawning_Event_profile();
     
                                       
/////////////////////////////////////////////////////////////////////////////////////////                                 
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
  
//This is a strategy for computing   Sum_mappings___P_mPER_cond_state_vector_and_brkpts.  It involves a lot of pre-computing.
              
// new strategy:
// Suppose we want to sum out event E.  Then we go through the LCR profile of E and look for all hybrid-reads.  We do this by looking at every single read that (may) map to E's profile (this is theoretically necessary) - whether at a variational position or not.
// Suppose we are at index x of E's profile:
//      1. We know which homologous reigons map here (down to the exact coordinate mapping), so we create miniRegions (intervals around homologous coordinates).  we then find all PERs from these miniRegions, and each PER is given this set of MiniRegions.  These PERs and their MiniRegions should be stored centrally (by .Name) - i.e. in the Event being summed out itself.
//      2. We then move to x+1.  Most of the PERs and associated miniRegions will be the same (we can check thanks to the central storage location).  Sometimes, we may find the same PER but a new MiniRegion - this MiniRegion should be added to that PER's already existing list.  Other times we may find a new PER - then we add it to the central storage of PERs with appropriate MiniRegions.
//   Thus, while we check every single position of the LCRs (not just variational positions), most of the informartion found is redundant.

//So in the end we have a unqiue set of PERs, and each PER has a unique set of MiniRegions (unique wrt itself).  For each PER, we take Event_being_summed_out itself and its interlocking events F and cross-check them against each PER's Mini-Regions.  As we do this, then for each PER we create some "MiniEvents"  (seen below) - events that will affect the landscape to which we can map this PER.



//Now we pre-compute some relevant mapping poistions, in preparation for summing later.
//For each PER...
  //1.  We check for any MiniRegions not associated with any MiniEvents.  The mapping of the PER to these regions will never change (these regions never change, under any state vector).  So we align the PERs to each of these Event-less MiniRergions, and compute the haploid sum.  This haploid sum will be added to any future sums for state vectors (it is constant!).
  //2.  We align the PER to every MiniRegion as is - this is covers the cases where the MiniRegion is not affected by a particular state vector.
  //3.  Now we cycle through each MiniEvent and align the PER to each rearragnement possibility (hybrids for each relevant breakpoint) and save these probabilities within that MiniEvent.
  
// From this set of PERs, we can aggregate the interlocking Events and aggregate their relevant breakpoints for each such interlocking Event.  Later, when summing for a state vector and all breakpoint combinations, we can easily "filter" these Events and breakpoints through the PER's set of aggregate (e.g. relevant) Events and breakpoints.  This hopefully will reduce the computations as many of the breakpioint combinations (e.g. those not included in the list of relevant ones), will result in the same Sum for that PER.  If we "save as we go", then we can reuse many of the same results rather than re-compute the same sum multiple times.  (e.g., given a state vector & brkpts, filter it, then compute the sum, thenb return the sum AND SAVE the sum.  then, in the future when we get another state vector & brkpts, we filter it, then check if this filtered state&brks have been computed already!)




      
      
      
      


    //then, we will want to pre-compute all relevant sums.  But we don't a priori know  which cobinations of Events are possible (e.g. state vectors) - this is a complicated question dealingw tih determining the variable elimiation schedule.
    //So a better solution is to  "save as we go".  That is, as we are given a state vector (no breakpoints), we:
    //          1.  Filter the state vector so it only contains events on this PER's MiniEvents list (these are the only relevant events, as far as this PER is concerned).
    //          2.  Check if the  sum_mPER  's have already been computed for the given filtered state vector.  If it has simply retrieve the sum, otherwise...
    //          3.  We then compute sum_mPER for all possible combinations of breakpoints for that filtered state vector, INCLUDING EXTREMAL BREAKPOINTS, and store each individually.  When computing the sum, note that we are also summing over all possible (relevant) breakpoints for the event being summed out, and storing THESE individually too.
    
    



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    public:
        
            
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {                                         
            ar  &  name;
            //Readgroup_statistics *my_Readgroup_ptr           // RECONNECT   
            ar  &  mates;
            ar  &  spawning_profile_chromo_and_position; 
            
            
            ar  &  my_MiniEvents;
            ar  &  my_MiniRegions;
            
            ar  &  set_of_MiniRegions_on__XXXX_chromosome;
            ar  &  set_of_MiniRegions_on__YYYY_chromosome;
            
            ar  &  map_affected_MiniID__to__P_nohybrid;
            
            ar  &  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_autosomal;
            ar  &  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_XXXX_chromosome;
            ar  &  haploid_sum_over_MiniRegions_not_affected_by_any_MiniEvents___only_YYYY_chromosome;
            
            
            
            ar  &  MiniID_ctr;
            ar  &  my_PER_aligner;
	    
	    ar  &  failed_PER_flag;
	    
	    ar  &  map_MiniID__to__nohybrid_almgt_interval_absolute;
        };
        
        
        
        
        
        
        void  deep_copy_from_serialize
                                (const Paired_end_read &another_PER,
                                 const Conn_Comp *const &the_conn_comp_ptr); 
                                 
        void reconnect_pointers_after_uploading_from_a_serialization
                                (const Conn_Comp *const &the_conn_comp_ptr);                                 
                                        
                                
                                
                                
                                
                      
                                
    
    
}; //end class   Paired_end_read    

















struct compare_PERs_by_name
{
  bool operator() (const Paired_end_read &LHS, const Paired_end_read &RHS) const
  {
    return   LHS.name.compare( RHS.name ) < 0;
    //return strcmp(LHS->reads.at(0).Name.c_str(), RHS->reads.at(0).Name.c_str())  <  0;
  }
};    





typedef std::set<Paired_end_read, compare_PERs_by_name>  type_set_PER;
typedef std::list<Paired_end_read>  type_list_PER;
typedef std::map<std::string, Paired_end_read>  type_map_string_to_PER;
typedef std::pair<std::string, Paired_end_read>   type_string__PER;













void  set_up_and_align_PERs_using_GPU____OMP
			(type_map_string_to_PER &PERs_for_alignment,
			 const Event *const spawning_Event,
			 const bool &create_and_consider_MiniEvents,
			const bool &consider_GeneConversion_breakpoints,  //irrelevant if "create_and_consider_MiniEvents" = false
			const bool &sum_over_all_unaffected_LCRs,  //"false" for heuristic
			const bool &special_save_breakpoints_actually_captured_in_alignment, //"true" for heuristic
			const bool &do_PER_orient_formatting);  // orthogonal to all other options








#endif
