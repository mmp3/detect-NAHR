#ifndef EVENT_H_INCLUDED
#define EVENT_H_INCLUDED


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
#include <Coords_and_homologous_region.h>

#include "Paired_end_read.h"
#include "MiniRegion.h"




class Conn_Comp;

class Event;
typedef std::pair<uint, Interlocking_entry<Event> >  type_uint__Interlock;
typedef std::map<uint, Interlocking_entry<Event> >  type_map_uint_to_Interlock;










class  Breakpoint_complex
{
    public:
	Breakpoint_complex()
	{  }
	
	~Breakpoint_complex()
	{  }

	    
	//variables:
	type_set_uint variational_positions__ALL; 
    
	type_set_BI  breakpoints_NAHR__AB;
	type_set_BI  breakpoints_NAHR__BA;
	
	type_set_BI   breakpoints_GeneConversion___ABA;
	type_set_BI   breakpoints_GeneConversion___BAB;
	
		
	//functions:
	void print_this_Breakpoint_complex
			    (std::stringstream *const &some_ss = NULL) const;
			    
			    
	void clear_all_data();
		    
    
};//Breakpoint_complex

typedef  std::map<uint, Breakpoint_complex>  type_map_uint_to_Breakpoint_complex;






























class Event : public SD_Entry_general
{
  public:
    //constructors:
    Event() : SD_Entry_general(), my_Conn_Comp(NULL)
    { }
        
    Event(const uint &given_UID,
          const Conn_Comp *const in_Conn_Comp_ptr)
                                : my_Conn_Comp(in_Conn_Comp_ptr),
                                  SD_Entry_general(given_UID)                                    
    { }
    
    Event(const uint &given_UID,
          const uint &given_UIDDB,
          const uint &in_CC,
          const Conn_Comp *const in_Conn_Comp_ptr)          
                                : SD_Entry_general(given_UID, given_UIDDB, in_CC),
                                  my_Conn_Comp(in_Conn_Comp_ptr)
    { }    
    
 
    //destructor:
    ~Event()
    {  }
    
    
    
    

    
    
    //variables:        
    
    const Conn_Comp *const my_Conn_Comp;
    
    
    type_set_haploid_outcome  potential_outcomes;
    
    
    uint last_profile_nongap_index;
    
                        
    //1.
    type_map_uint_to_Interlock local_interlocking_Events;   //indexed by UID
    // 2.
    type_map_uint_to_Interlock global_interlocking_Events;   //indexed by UID
    
     
    Breakpoint_complex  my_original_Breakpoint_complex;
    
    type_set_uint  uploaded_variational_positions_of_profile;
    
    
    
    type_set_uint remaining_coordinates_for_consideration;
             
    type_map_uint_to_uint_to_uint_to_longdouble   natural_poisson_intervals___gender_INsensitive; // includes homologous regions.  
         // color   chromo    pos      aggregate_fragmentation_rate
         
    type_map_uint_to_uint_to_uint_to_longdouble  base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals;
    
    type_map_uint_to_longdouble base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive;
        
    type_map_uint_to_uint observed_read_depth_for_natural_poisson_intervals;            

    type_map_uint_to_list_BI__string uploaded_pieces_of_chromos;           
    
    type_map_string_to_PER PERs_on_this_profile;    //name to PER
        
    type_map_uint_to_set_BI filter_for_NON_large_gap_of_neighbors_and_self__that_affect_Read_depth;  
    type_set_uint RD_affecting_neighbors_self;    
                  
    
    
    type_set_BI  all_NON_large_gap_intervals_of_profile;
                                
                                        
//     type_map__uint__uint___to___longdouble___map_uint_to_longdouble marginal_probability_distribution;
            
    
    
    type__set_uint__set_uint absolute_coordinates_inside_of_large_gaps_of_profile;
    
    
    
    //functions:   
    void print_this_entry(std::stringstream *const &some_ss = NULL) const;
            
                    
    
    void pre_upload_homologous_regions_from_the_Reference_genome
			(const int &padding_amount_indicator = -1);   

    
    
    
    type_map_uint_to_set_uint determine_homologous_variational_positions_on_immediate_neighbors_and_self() const;
    
    
            //read depth:                                                                                              
    
    void partition_DAR_by_homology___and___create_natural_poisson_intervals___NOT__gender_sensitive();
    void create_natural_poisson_intervals___NOT_gender_sensitive( const type_map_uint_to_set_uint &partition_of_DAR_by_homology);
    
    void calculate_total_HAPLOID_GC_rate_for_every_base_by_summing_over_all_Readgroup_GC_rates_at_each_position();    
    void calculate_base_diploid_fragmentation_rate___haploid_sensitive___gender_sensitive__and_base_cumulative_haploid_fragmentation_rate();
       
    void calculate_observed_read_depth_for_natural_poisson_intervals
						(BamTools::BamReader &my_BAM_reader); 
            
    void adjust_expected_diploid_fragmentation_rate___breakpoint_specific
                                            (type_map_uint_to_longdouble &expected_diploid_fragmentation_rate_per_natural_poisson_interval__for_these_states,
                                             const type_map_uint_to_uint &proposed_haploid_state_vector,
                                             const type_map_uint_to_BI &proposed_haploid_breakpoints)  const;
    
    real calculate_probability_of_observed_read_depth_given_expected_read_depth_over_natural_poisson_intervals
                                    (const type_map_uint_to_longdouble &expected_diploid_fragmentation_rate_per_natural_poisson_interval)  const;
    
    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
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
    
    
    
    //functions:
    bool set_up_all_PERs
		(const bool &consider_GeneConversion_breakpoints,
		 BamTools::BamReader &my_BAM_reader);
    
    type_list_MiniRegion create_MiniRegions_from_position( const uint &position_on_Event_profile_or_DAR )  const;            
    


  
            
    real calculate_or_get_Read_Depth_for_these_state_vectors_and_breakpoints
                                        (const Sparse_map &RD_affecting_states_0,
                                         const Sparse_map &RD_affecting_states_1,
                                         const type_map_uint_to_BI &RD_affecting_haploid_breakpoints_0,
                                         const type_map_uint_to_BI &RD_affecting_haploid_breakpoints_1,
					 std::stringstream *const &display_strm__ptr = NULL)  const;      


					 
    
    
    

    type_map_uint_to_list_BI create_search_intervals_from_DAR_and_homologous_regs_for_uploading_PERs
                                                                (const uint &search_interval_padding_amount,
                                                                 const uint &subdivison_size);   
                                                                
                                                                
    bool upload_BAM_alignments___create_PERs___find_all_MiniRegions_______from_search_intervals_____MPI_and_OMP
                                (const type_map_uint_to_list_BI &map_chromo_to_homol_regs,
                                 const uint &MR_padding_amount_for_proper_pair,
				BamTools::BamReader &my_BAM_reader);
                                 
                                 
    void for_all_PERs____set_region_sequences___and___orient___MPI_and_OMP();                              
    
                     
    bool determine_relation_of_region_to_Event_profile
                                (bool &orientations_of_region_and_Event_profile_agree,
                                 bool &must_complement_region_to_get_to_Event_profile,
                                 const uint &chromo_of_region,
                                 const BOOST_Interval interval_of_region);
                                 
                           
                                 
    void clear_all_data();
    
    
    
    void upload_or_create_natural_poisson_intervals___gender_INsensitive();  
    
    
    
    


    void ignore_regions_of_the_genome_from_consideration_by_this_Event
                                    (const type_map_uint_to_list_BI  &regions_to_ignore);    
				    
				    
    type_longdouble__longdouble  calculate_expected_number_of_PERs_at_breakpoint
			    (const BOOST_Interval &profile_breakpoint,
			    const MiniRegion &mini_LCR_0,
			    const MiniRegion &mini_LCR_1,
			     const bool  &AB_ABA__false_____BA_BAB__true)  const;			    
    
    
    void set_absolute_coordinate_mapping_into_large_gaps();   
       
    
    
                                            
    
}; //end class Event    










struct compare_Event_ptrs_by_UID
{
    bool operator() (const Event *const &Aptr, const Event *const &Bptr)
    {
        return  Aptr->UID  <  Bptr->UID;
    }
};








typedef std::map<const Event*, uint, compare_Event_ptrs_by_UID>   type_map_constEvent_ptr_to_uint;
typedef std::pair<const Event*, uint>  type_constEvent_ptr__uint;

typedef std::map<const Event*,  type_haploid_outcome, compare_Event_ptrs_by_UID>   type_map_constEvent_ptr_to_hapout;
typedef std::pair<const Event*, type_haploid_outcome>  type_constEvent_ptr__hapout;


typedef std::map<uint, Event>  type_map_uint_to_Event;
typedef std::pair<uint, Event>   type_uint__Event;

typedef std::map<uint, Event*>  type_map_uint_to_Event_ptr;
typedef std::vector<Event*>  type_vector_Event_ptr;

typedef std::set<Event*, compare_Event_ptrs_by_UID>   type_set_Event_ptr;

typedef std::pair<uint, const Event*>   type_uint__constEvent_ptr;




struct compare_uint__constEvent_ptr
{
    bool operator() (const type_uint__constEvent_ptr &LHS, const type_uint__constEvent_ptr &RHS)
    {
        if (LHS.first < RHS.first)
            return true;
        else if (RHS.first < LHS.first)
            return false;
        else
            return LHS.second->UID < RHS.second->UID;                    
    }                                
};   












class  Sampled_diploid_Event_data
{
  public:
    Sampled_diploid_Event_data()
				:	cc_of_event(0),
					event_UID(0),
					the_diploid_Event_outcome(hap_outcome__None,hap_outcome__None),
					the_diploid_profile_brkpts(empty_BI, empty_BI),
					P_diploid_outcome(-1),
					P_diploid_breakpoint(-1)
					
    { }
    
    Sampled_diploid_Event_data( const uint &in_CC,
				const uint &in_UID,
				const type_haploid_outcome__haploid_outcome &in_outcome,
				const type_BI__BI &in_brkpts)
						:   cc_of_event(in_CC),
						    event_UID(in_UID),
						    the_diploid_Event_outcome(in_outcome),
						    the_diploid_profile_brkpts(in_brkpts),
						    P_diploid_outcome(-1.00),
						    P_diploid_breakpoint(-1.00)
    { }    
    
    ~Sampled_diploid_Event_data()
    { }
    
    
    
    //variables:    
    uint cc_of_event;
    uint event_UID;
    type_haploid_outcome__haploid_outcome  the_diploid_Event_outcome;
    type_BI__BI  the_diploid_profile_brkpts;
    
    longdouble P_diploid_outcome;
    longdouble P_diploid_breakpoint;
    
    std::string  associated_genome;
    

        
    
   
    
    //functions:
    void print_this_Sampled_diploid_Event_data(std::stringstream *const &some_ss = NULL)  const;
    
    
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {            
	//note that we ignore the pointer for the Eliminating Event. this will be restord within Visitation_Object.
	ar  &  cc_of_event;
	ar  &  event_UID;	
	ar  &  the_diploid_Event_outcome;
	ar  &  the_diploid_profile_brkpts;
	ar  &  P_diploid_outcome;
	ar  &  P_diploid_breakpoint;
	ar  &  associated_genome;	
    };    
    
            
};//Sampled_diploid_Event_datas

                        
typedef  std::map<uint, Sampled_diploid_Event_data>   type_map_uint_to_Sampled_diploid_Event_data;
typedef  std::list<Sampled_diploid_Event_data>  type_list_Sampled_diploid_Event_data;        
typedef  std::vector<type_map_uint_to_Sampled_diploid_Event_data>   type_vector_map_uint_to_Sampled_diploid_Event_data;

typedef  std::map<std::string, type_map_uint_to_Sampled_diploid_Event_data>     type_map_string_to_uint_to_Sampled_diploid_Event_data;
typedef  std::map<std::string, type_list_Sampled_diploid_Event_data>  type_map_string_to_list_Sampled_diploid_Event_data;
            

























typedef  std::map<type_BI__BI, longdouble, compare_BI>   type_map__BI__BI__to__longdouble; 
typedef  std::map<type_haploid_outcome__haploid_outcome, type_map__BI__BI__to__longdouble>   type_map__hap_outcome__hap_outcome__to__BI__BI__to__longdouble;
 
class  Marginal_Event_posterior_distribution
{
  public:
    Marginal_Event_posterior_distribution()
    {  }  
    
    Marginal_Event_posterior_distribution(const uint &in_CC_id,
					  const uint &in_EV_ID)
						:  CC_ID(in_CC_id),
						    ev_UID(in_EV_ID)
    {  }        
    
    ~Marginal_Event_posterior_distribution()
    {  }
    
    
    uint CC_ID;
    uint ev_UID;
    
    type_map__hap_outcome__hap_outcome__to__longdouble  map_diploid_outcomes_to_posterior_probability;
    type_map__hap_outcome__hap_outcome__to__BI__BI__to__longdouble  map_diploid_outcomes_to_diploid_breakpoints_to_posterior_probability;  //event outcome to breakpoints to posterior probability of breakpoints    

};//Marginal_Event_posterior_distribution
 
 
typedef  std::map<uint, Marginal_Event_posterior_distribution>  type_map_uint_to_Marginal_Event_posterior_distribution;
























#endif