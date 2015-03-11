#ifndef PARTIAL_SUM_FUNCTION_H_INCLUDED
#define PARTIAL_SUM_FUNCTION_H_INCLUDED


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

// #define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>


#include <mpreal.h>
#include <mpfr.h>



#include <general_typedefs.h>
#include <Sparse_map.h>

#include "globals.h"



class Conn_Comp;
class Event;

class Partial_sum_function; 
typedef std::list<Partial_sum_function>   type_list_Partial_sum_fcn;
typedef std::list<type_list_Partial_sum_fcn::iterator>   type_list_Partial_sum_fcn_iter;
// typedef std::vector<Partial_sum_function>   type_vector_Partial_sum_fcn;

class  Sampled_diploid_Event_data;
typedef  std::map<uint, Sampled_diploid_Event_data>   type_map_uint_to_Sampled_diploid_Event_data;



class Partial_sum_function
{
    public:	
        Partial_sum_function() : eliminating_Event_UID(-2)
        {  }
                     
        
        Partial_sum_function( const type_set_uint &in_haploid_argument_UIDs, 
			      const type_map_uint_to_set_uint& in_current_relevant_variational_positions_per_Event_in_PSF)  
							: haploid_argument_UIDs(in_haploid_argument_UIDs),
							eliminating_Event_UID(-3),
							current_relevant_variational_positions_per_Event_in_PSF(in_current_relevant_variational_positions_per_Event_in_PSF)
        {  }       
                        
        
        Partial_sum_function( const type_set_uint &in_haploid_argument_UIDs,
                              const uint &in__eliminating_Event_UID)
						    : haploid_argument_UIDs(in_haploid_argument_UIDs),
						    eliminating_Event_UID((int)in__eliminating_Event_UID)
        {  }
                                         
        
        
        ~Partial_sum_function()
        {  }
        

    
        //a partial sum function  has an m-tuple as an argument (whose value must be permissible), and returns a number (a double).
            //for simplicity, a partial sum function may accept a map_uint_to_bool of any size, and it can go through and look for the UIDs that pertain to it, and check that the values of the m-tuple are permissible, and finally return the value of the function.


    //since we have to store data for the diploid case, a partial sum function should accept a 2-tuple (pair) of (decimal) integers and return a  long double.        
       
      
      
      
      
      
      
      
        //variables:	
	
        type_set_uint haploid_argument_UIDs; //does NOT include  "eliminating_Event_UID" (if it is defined)
        
        int eliminating_Event_UID;

	type_map_uint_to_set_uint  current_relevant_variational_positions_per_Event_in_PSF; //does NOT include "eliminating_Event_UID" (if it is defined)
	
	
	
	
	
	
	
        
        //functions:   
	void print_this_PSF__states_only
			    (std::stringstream *const &some_ss = NULL) const;        
        
        void print_this_PSF(std::stringstream *const &some_ss = NULL) const;             
        
        void set_eliminating_Event__and__haploid_argument_UIDs_according_to_relevant_PSFs
                                                (const uint &in_eliminating_Event_UID,
                                                  const type_list_Partial_sum_fcn_iter &list_of_related_PSF,
                                                  const type_set_uint &haploid_UID_arguments_for_new_PSF);      
        
        
        bool check_if_event_is_an_argument_to_this_function(const uint &test_UID) const;
                    //only for neighbors, - will return false if you give this function  UID_of_event_being_summed_out
                    
        type_set_uint retrieve_haploid_argument_UIDs() const;
	
	
	
	void record_additional_relevant_variational_positions_per_Event
				(const type_map_uint_to_set_uint &map_more_relevant_varpos_per_Event);
				
        
	
        
        real get_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints
                                            ( const Sparse_map &oversized_haploid_sparse_state_vector_0,
                                            const Sparse_map &oversized_haploid_sparse_state_vector_1,
                                            const type_map_uint_to_BI &oversized_haploid_breakpoints_0,
                                            const type_map_uint_to_BI &oversized_haploid_breakpoints_1 ) const;			    
    
    
                                            
	void  allocate_diploid_breakpoints_for_given_diploid_states
		    ( const Sparse_map &oversized_haploid_sparse_state_vector_0,
			const Sparse_map &oversized_haploid_sparse_state_vector_1,
		      const type_map_uint_to_vector_BI  &breakpoint_space__hap_0,
		      const type_map_uint_to_vector_BI  &breakpoint_space__hap_1);
        
        
        
        //4. for backtrace:
        
        type_multimap_real_to_uint__uint calculate_and_return_eliminated_Event_diploid_outcome_CDF
                                                (const type_map_uint_to_Sampled_diploid_Event_data &map_estimate_for_this_CC)  const;        
        
                                            
        type_multimap_real_to__BI__BI calculate_and_return_eliminated_Event_diploid_breakpoint_CDF
                                                (const type_map_uint_to_Sampled_diploid_Event_data &map_estimate_for_this_CC)  const;        
        
        
        
                                            
        void save_diploid_partial_sum_for_these_sparse_state_vectors_and_breakpoints
                                        (const Sparse_map &oversized_haploid_sparse_state_vector_0,
                                          const Sparse_map &oversized_haploid_sparse_state_vector_1,
                                          const type_map_uint_to_BI &oversized_haploid_breakpoints_0,
                                          const type_map_uint_to_BI &oversized_haploid_breakpoints_1,
                                          const real &contribution_to_partial_sum_value_S);                        
                                            
                                            
            
        void save_Cartesian_product_of_related_PSFs_wrt_eliminating_Event_and_include_probabilities_of_Breakpoint_and_Event
                                            (const Event  &eliminating_Event,
					     const Sparse_map &full_haploid_sparse_state_vector_0,
                                             const Sparse_map &full_haploid_sparse_state_vector_1,                                           
                                             const type_list_Partial_sum_fcn_iter &list_of_related_PSF,
                                             const Partial_sum_function &thread_specific____new_PSF);    
					    
            
            
                    
	void marginalize_over_eliminating_Event_breakpoints_and_outcomes();                                                
                     
        
        void remove_eliminated_Event_from_storage
				(const uint &eliminated_Event);
                                                
        
                                   
                                        
        void clear_data();           
        
        
        void make_numerically_stable_by_multiplying_every_partial_sum_value_by_the_inverse_of_the_smallest_partial_sum_value();
        
        
        
        
        
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {            
            ar  &  haploid_argument_UIDs;
	    ar  &  eliminating_Event_UID;	    
	    ar  &  current_relevant_variational_positions_per_Event_in_PSF;	    
            ar  &  map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage;
        };
        
                                                    

	
	


        class traceback_data_storage
        {
            public:                 
                traceback_data_storage()  : partial_sum_value_S(0.00L)
                { }
                
                traceback_data_storage( const real &in_partial_sum_value_S)  : partial_sum_value_S(in_partial_sum_value_S)
                { }            
                   
        //  S( (E,b)2 , (E,b)3 )  =    Sum_E1                                             prod_read_depth_factors  *  Sum_b1  prod_P_PER_cond_(E,b)1,2,3
        //  ^ partial_sum_value                 ^ P_diploid_eliminating_Event_outcomes      ^ not stored here           ^hybrid_factor^         
        //                                        (each individual term)                                               (sum and each individual term)
        
        
                real partial_sum_value_S;  //summed over diploid eliminating Event outcomes (and thus, implicitly, summed over diploid eliminating Event brkpts.)
                                            //i.e. summed over "P_D__cond__diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states"  for each eliminating Event state.                                                  
                
                
                // the name "hybrid factor" no longer refers to just hybrids, but to read-depth as well.  I leave it for now, although it may cause confusion.
                class hybrid_factor  // sum over breakpoints in eliminating Event ( prod over PERs ( sum over mappings ( P_mPER_cond_states_and_neighbors)))
                {    
                    public:
                        hybrid_factor()  :  P_D__cond__diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states(0.00L)
                        { }                                              
                        
                        type_map_BI__BI___to__real               
                           map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts;
                                                                        //this is for sampling (traceback).   
                                        // key = diploid breakpoint in eliminiating Event
                                        // value =  read_depth_factor * prod over PERs ( sum over mappings ( P_mPER_cond_states_and_neighbors))
                                
                                
                        real P_D__cond__diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states;
                                                    //i.e. marginalized over breakpoints of eliminating_Event given its outcome.
                                                    
                        friend class boost::serialization::access;
                        template<class Archive>
                        void serialize(Archive &ar, const unsigned int version)
                        {
                            ar  &  P_D__cond__diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states;
                            ar  &  map_eliminating_Event_diploid_brkpts__to__P_D_cond_diploid_neighbor_states_brkpts__and__eliminating_Event_diploid_states_and_brkpts;                            
                        }                                                    
                                                    
                };           
                
                
                typedef std::map<type_uint__uint, hybrid_factor>  type_map_eliminating_Event_diploid_outcomes_to_hybrid_factor;
                type_map_eliminating_Event_diploid_outcomes_to_hybrid_factor map_eliminating_Event_diploid_outcomes_to_hybrid_factor;
        
        
                
                friend class boost::serialization::access;
                template<class Archive>
                void serialize(Archive &ar, const unsigned int version)
                {
                    ar  &  partial_sum_value_S;
                    ar  &  map_eliminating_Event_diploid_outcomes_to_hybrid_factor;                    
                }

                
        };    // end of    class traceback_data_storage    
        
        
        
        
        
        

        
        //private data:

        typedef std::pair<type_map_uint_to_uint, type_map_uint_to_uint> type_map_uint_to_uint__2;
            
        typedef std::map< type_map_uint_to_BI__2, traceback_data_storage, compare_pairs_of_maps< uint, BOOST_Interval, uint, BOOST_Interval, compare_BI, compare_BI >  >
                                    type_map_diploid_neighbors_brkpts_to_traceback_data_storage;

			
	typedef   std::map< type_map_uint_to_uint__2, type_map_diploid_neighbors_brkpts_to_traceback_data_storage, 
			    compare_pairs_of_maps<uint, uint, uint, uint, std::less<uint>, std::less<uint> > >
			type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage;
                                        
								
								

	type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage
			    map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage;
				
				
				
                                                                 
        //private functions:
		
        
                        
        void absorb_piece_of_Partial_Sum_Function_into_one_Partial_Sum_Function
                                    (const Partial_sum_function &piece_of_new_Visitation_PSF); 
    
        void  copy_individual_likelihood_files_to_common_likelihood_directory_for_this_Event
                                            (const uint &given_likelihood_eliminating_Event_UID)  const;
    
        static void  combine_individual_likelihood_files_into_common_likelihood_file_for_this_Event
                                            (const uint &given_likelihood_eliminating_Event_UID);
					    
	void upload_all_likelihood_factors_from_file_for_a_given_Event
				    (const std::string &other_job_run_data_dir,
				    const uint &some_eliminating_Event_UID);
                                                                                                                
                     
				    
				    
	static void absorb_piece_of_diploid_neighbor_brkpts
				(const type_map_diploid_neighbors_brkpts_to_traceback_data_storage &diploid_neighbor_brkpts___source,
				 type_map_diploid_neighbors_brkpts_to_traceback_data_storage &diploid_neighbor_brkpts___destination);				    
				    
				    
				    
	static void deep_insert
		    (type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage &receiver,
			const type_map_diploid_neighbor_states_to_diploid_neighbors_brkpts_to_traceback_data_storage &donor);
				 
				 

	type_map_diploid_neighbors_brkpts_to_traceback_data_storage::const_iterator  
	get_iterator_to_traceback_storage_satisfying_diploid_states_and_generalized_diploid_breakpoints
			    (const Sparse_map &oversized_haploid_sparse_state_vector_0,
				const Sparse_map &oversized_haploid_sparse_state_vector_1,
				const type_map_uint_to_BI &oversized_haploid_breakpoints_0,
				const type_map_uint_to_BI &oversized_haploid_breakpoints_1)  const;
				    
};




// Partial_sum_function:

// typedef std::vector<Partial_sum_function>   type_vector_Partial_sum_fcn;
// typedef std::vector<type_vector_Partial_sum_fcn::iterator>   type_vector_Partial_sum_fcn_iter;








// We always need to collect all Partial sum functions and aggregate and store them for each Event.  This is because the traceback is entirely based on going from one Event to another.  This raises an issue though.  When summing forward, a single Event E may pertain to multiple PSF.  By our assumption of "families" of duplicons, each of these PSF has a DISTINCT set of breakpoints of E (this can be accomplished by an apprpriate  breakpoint "erasing/elimination" strategy for the forward sum).  BUT, each PSF will have a different, possible non-disjoint set of argument events.  Thus, in order to aggregate all PSFs into one PSF for Event E, we need to take the "Cartesian product" over the PSFs.  This is slightly complicated by the fact that some PSFs may contain the same event F, though again, with distinct breakpoints; thus the Cartesian product may contain some "redundancies" in terms of Events (but not breakpoints), which we should be wary of (i.e. the same Event cannot have two breakpoints!).

// We always have a state vector that defines a stae for each of the interlocking Events, i.e. every argument of all related PSFs has been given a state.So for each given neighbor state & EE state,we need to Cartesian product over neighbor breakpoints and EE brkpts.




// We have a state vector iterator.  This is relevant to the EE and all related PSF.  Since each PSF and the current EE are disjoint WRT BREAKPOINTS, then we can iterator over the state vectors and calculcate a new PSF for the EE.  
//This leaves us with a set of old, related PSFs, and a new PSF which took into consideration the reminaing coordiantes of the EE.  Now we need to do the Cartesian product business.   Given a state vector....
//      Every PSF accepts this state vector as an argument (maybe they only consider part of it.) We can't do this naively, since some of the PSF may share an Event (albeit distinct breakpoints), andw e wouldn't want to say the same Event had two differetn breakpoints!  Thus, we simply conglomerate each PSF's map's together.
//      For these fixed states, we need to take a Cartesian product over all relevant breakpoints.





#endif