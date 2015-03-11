#ifndef CONN_COMP_H_INCLUDED
#define CONN_COMP_H_INCLUDED


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


#include "Event.h"
#include "Visitation_object.h"
#include "Partial_sum_function.h"
#include "globals.h"








class Sampled_diploid_Event_data;






class Conn_Comp
{
    //We represent a set of globally connected components as a Conn_Comp, where the Conn_Comp holds an array of Events.  For this, we need nested classes.        
    public:         
        //constructor:
        Conn_Comp()  :  Conn_Comp_ID(0)
        { }
        Conn_Comp(const uint &in_cc)  :  Conn_Comp_ID(in_cc)
        { }
        




        //variables: 
        
        const uint Conn_Comp_ID;   
        
        type_map_uint_to_Event events;
        

        type_list_Visitation_object variable_elimination_schedule;
        uint total_diploid_compute_size;
        
        
        type_list_Partial_sum_fcn partial_sum_functions;
        
        
        
        

        
        //functions:                                 
                        
        int perform_complete_variable_elimination_schedule
			(BamTools::BamReader &my_BAM_reader,
			 const bool &consider_GeneConversion_breakpoints_in_Heuristic,
			 const bool &consider_GeneConversion_breakpoints_in_Full_compute);        
                        
        


                                
                                
                                
        
        
                
        bool HEURISTIC__test_each_specific_variational_position_possibility
                                                    (type_set_uint &excluded_Events,
						     const bool &consider_GeneConversion_breakpoints,
						    BamTools::BamReader &my_BAM_reader);
        
        

	void output_reads_covering_certain_breakpoint_regions                             
					(const Sampled_diploid_Event_data &some_result,
					const bool &do_likelihood_ratio_test_on_alignments,
					BamTools::BamReader &my_BAM_reader,
					const bool &display_ALL_reads_not_just_good_ones = false);
					
					
	void  create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights
				(BamTools::BamReader &my_BAM_reader,
				 const Sampled_diploid_Event_data &some_result,
				 const std::string &outdir);
                             
                                                    
      
                                                    
                                                    
        type_map_uint_to_Sampled_diploid_Event_data get_and_save_MAP_estimate();
        
        
        
        type_map_uint_to_Sampled_diploid_Event_data get_a_single_sample()   const;
	
	
        type_map_uint_to_Marginal_Event_posterior_distribution compute_marginal_probabilities_from_posterior_sample
                                                        (const type_vector_map_uint_to_Sampled_diploid_Event_data &posterior_samples)  const;
                            
        void sample_from_posterior__and__save____marginal_distributions__and__outcome_Centroid__and_Shannon__Entropies
                                                        (const uint &number_of_samples);          
                                                        
                                                        
                                                        
        void simulate_elimination_schedule_and_save_natural_poisson_intervals();                                                        
                                                        
                            
	
	int use_saved_likelihoods_to_perform_complete_variable_elimination_schedule_with_current_prior
				(const std::string &old_run_data_dir,
				 const type_set_uint &excluded_Events_from_old_run );	
				 
// 	void display_PERs_on_full_MAP_space_around_breakpoint
//                             (const uint &EV_UID,
//                              const uint &desired_breakpoint);
	
				
// 	int reorganize_saved_Likelihoods_for_all_connected_components
// 				(const std::string &old_run_data_dir );
                            
   
                              
}; //end class Conn_Comp
 
     
     

// Conn_Comp:
typedef std::map<uint, Conn_Comp>  type_map_uint_to_CC;
typedef std::pair<uint, Conn_Comp>   type_uint__CC; 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

#endif