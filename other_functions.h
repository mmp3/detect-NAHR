#ifndef OTHER_FUNCTIONS_H_INCLUDED
#define OTHER_FUNCTIONS_H_INCLUDED


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
#include <api/BamReader.h>
#include <api/BamAux.h>



#include <general_typedefs.h>





class Conn_Comp;
class Call_set;
class Event;
class MiniRegion;
class Visitation_object;
class Paired_end_read;
class Readgroup_statistics;
class Sampled_diploid_Event_data;
class MiniEvent;
class Breakpoint_complex;

           
class  Marginal_Event_posterior_distribution;
typedef  std::map<uint, Marginal_Event_posterior_distribution>  type_map_uint_to_Marginal_Event_posterior_distribution;           
                        




                        
                        
                        
void remove_degenerate_Events
		    (std::map<uint, Conn_Comp> &Conn_Comps);                        
                        
void remove_undesirable_Connected_Components_from_set_of_Connected_Components
		    (std::map<uint, Conn_Comp> &Conn_Comps,
			const type_set_uint &targeted_CC_IDs,
			const type_set_uint &display_Event_and_breakpoint_range,
			const BOOST_Interval &acceptable_Connected_Component_size_range,
			const BOOST_Interval &acceptable_Connected_Component_compute_range,
			const int &last_completed_CC);                       
                        
                        
		    
		    
//BAM operations:                        
                        
void count_decent_reads_with_left_endpoints_in_set_region
	(BamTools::BamReader &a_BAM_reader,  //region already set
	const BOOST_Interval &absolute_genome_interval,
	type_map_uint_to_uint &counts_per_position);
	
uint count_decent_reads_with_left_endpoints_in_set_region
	(BamTools::BamReader &a_BAM_reader,  //region already set
	const BOOST_Interval &absolute_genome_interval);
	
                                                    
type_uint__uint create_list_of_unformatted_Paired_end_reads_from_list_of_Balgmts
				(std::list<BamTools::BamAlignment> &Balgmts_from_a_Mini_region,
				BamTools::BamReader &my_BAM_reader,
				std::list<Paired_end_read> &created_PERs);
// WARNING: this function erases "Balgmts_and_lower_padded_profile_interval_where_found".  It should not be used upon return!                                


BamTools::BamRegion create_BAM_region
				(const uint &some_chromo_value,
				const BOOST_Interval &some_interval,
				const uint interval_padding_amount);				

void set_map_from_chromosome_values_to_RefID_in_BAM_file
                                    (BamTools::BamReader &a_BAM_reader);                                  

 
		    
void prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population
		    (const std::string &BAM_filename,
		     const std::string &gender_and_population_filename,
		    BamTools::BamReader &my_BAM_reader);
		     
                                  
				    
				    




//hybrid manipulation:                                          
                         
//"A" and "B" are abstract - they do not necesssariyl map as "A" --> "LCR 0",  "B" --> "LCR 1".   (it may be reversed!!!)
//This function always creates hybrids of the form "AB" or "ABA"
//Thus, if you want a "BA" or "BAB" hybrid, then you have to "trick" the function by swapping the coordinates.
//"breakpoints_on_A/B_mini_region_absolute" mark the actual breakpoint on each sequence.  Inclusion of the breakpoint depends, of course on NAHR vs. GeneConv.
type_string__BI  get_hybrid_sequence_and_breakpoint_relative_to_hybrid__by_hybridizing_two_MiniRegions_at_specific_breakpoint
                                    (const MiniRegion &A_mini,
                                     const BOOST_Interval &breakpoints_on_A_mini_region_absolute,  //not to be included!!
                                     const bool &must_complement_mini_A_to_get_to_profile,
                                     const bool &A_mini_must_be_reversed,   //  0 = lower() is boundary,  1= upper() is boundary
                                     const MiniRegion &B_mini,
                                     const BOOST_Interval &breakpoints_on_B_mini_region_absolute, //this one is included!!
                                     const bool &must_complement_mini_B_to_get_to_profile,                                     
                                     const bool &B_mini_must_be_reversed, //0 = lower() is boundary,  1= upper() is boundary
				     const bool &is_NAHR); //else, GeneConv
                         
				    
//"A" and "B" are abstract - they do not necesssariyl map as "A" --> "LCR 0",  "B" --> "LCR 1".   (it may be reversed!!!)
//This function always creates hybrids of the form "AB" or "ABA"
//Thus, if you want a "BA" or "BAB" hybrid, then you have to "trick" the function by swapping the coordinates.
//"breakpoints_on_A/B_mini_region_absolute" mark the actual breakpoint on each sequence.  Inclusion of the breakpoint depends, of course on NAHR vs. GeneConv.
void  get_variational_indeces_mapped_onto_hybrid_sequence
			(const MiniRegion &A_mini,
			const BOOST_Interval &breakpoints_on_A_mini_region_absolute,  //not to be included!!
			const bool &A_mini_must_be_reversed,   //  0 = lower() is boundary,  1= upper() is boundary
			const MiniRegion &B_mini,
			const BOOST_Interval &breakpoints_on_B_mini_region_absolute, //this one is included!!
			const bool &B_mini_must_be_reversed,     //  0 = lower() is boundary,  1= upper() is boundary
			 const type_set_uint &absolute_variational_positions__on_A,
			 const type_set_uint &absolute_variational_positions__on_B,
			 const bool &is_NAHR,
			 type_set_uint &hybrid_varpos_indeces__A,
			 type_set_uint &hybrid_varpos_indeces__B); 


//"0" and "1" stand for the TRUE LCR sides.
type_vector_int  label_absolute_coordinates_on_hybrid
		    (const MiniRegion &mini_LCR_0,
		     const BOOST_Interval &absolute_breakpoint_region__LCR_0,
		     const MiniRegion &mini_LCR_1,
		     const BOOST_Interval &absolute_breakpoint_region__LCR_1,
		     const bool &is_NAHR,  //else, GeneConv
		     const bool  &AB_ABA__false_____BA_BAB__true,
		     const type_recomb_class  &recomb_type_of_Event);
		     			                                                                                               
                    
type_map_uint_to_uint get_relative_hybrid_coords_to_algmt_profile
                        (const type_string__string &PER_and_hybrid_algmt,
                         const BOOST_Interval &endpoints_of_hybrid_on_algmt);                                             
                                                                                          			 
			 
void adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
                        (const type_map_uint_to_uint &hybrid_indeces_to_algmt_profile,
                         const BOOST_Interval &hybrid_breakpoints,//&index_of_hybrid_AB_that_is_first_B_nucleotide_IE_the_breakpoint_varpos,
                         const uint &absolute_genome_position_of_A_mapping_to_index_0_of_hybrid,  //Note, index 0 of hybrid is NOT  necessarily the first index in the algmt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         const bool &A_mini_was_reversed,
                         const uint &absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid,
                         const bool &B_mini_was_reversed,
			 const bool &is_NAHR,
                         type_map_uint_to_uint &absolute_A_coords_to_algmt_profile,
                         type_map_uint_to_uint &absolute_B_coords_to_algmt_profile);

void adjust_map_keys_from_relative_hybrid_coords__to__absolute_genome_coordinates
                        (const type_map_uint_to_uint &hybrid_indeces_to_algmt_profile,
                         const BOOST_Interval &hybrid_breakpoints,//&index_of_hybrid_AB_that_is_first_B_nucleotide_IE_the_breakpoint_varpos,
                         const type_uint__uint &absolute_genome_position_of_A_mapping_to_index_0_of_hybrid___and___to_second_GeneConv_breakpoint_on_hybrid,  //Note, index 0 of hybrid is NOT  necessarily the first index in the algmt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         const bool &A_mini_was_reversed,
                         const uint &absolute_genome_position_of_B_mapping_to_breakpoint_of_hybrid,
                         const bool &B_mini_was_reversed,
			 const bool &is_NAHR,
                         type_map_uint_to_uint &absolute_A_coords_to_algmt_profile,
                         type_map_uint_to_uint &absolute_B_coords_to_algmt_profile);
                         
                         
                         

                         
                         
                         
                         
                         
   
			
                     
		     
		     
type_map_uint_to_list_BI  scan_region_for_unidentified__N__nucleotides_and_merge_for_ignoring
                                        (const type_map_uint_to_set_uint &map_chromo_to_pos,
                                         const int &min_allowable_gap_size_between_consecutive_N_nucleotides,
                                         const int &min_size_for_ignoring);   
                                         		     
		     
		     
bool check_that_relevant_regions_are_well_defined__and__hide_undefined_regions_if_necessary
			(const std::map<uint, Event>::iterator  &event_being_summed_out,
			 std::map<uint, Event> &mutable_events_of_this_CC);			 
					 
void  make_gender_sensitive_iterators_for_state_vectors
					(std::vector< type_vector_Sparse_map::const_iterator >  &iterators_to_haploid_0_ie_XXXX_sensitive_state_vectors,
					 std::vector< type_vector_Sparse_map::const_iterator >  &iterators_to_haploid_1_X_or_Y_sensitive_state_vectors,
					 const std::map<uint, Event> &events_of_some_CC,
					 const std::list< Visitation_object >::const_iterator &it_Visit);					 
					 
void eliminate_GeneConversion_event_outcomes
	    (std::vector<type_vector_Sparse_map::const_iterator>  &iterators_to_haploid_state_vectors);			
					 
					 

	    
	    
	    
	    
	    
	    

	    
//miscellaneous	   
	    

uint safe_subtract_base_1
	(const uint &some_value__base_1,
	 const uint &subtract_this_amount);

uint safe_chromo_add_base_1
	(const uint &some_value__base_1,
	 const uint &add_this_amount,
	 const uint &chromo_val);

uint safe_subtract_base_1
	(const uint &some_value__base_1,
	 const int  &subtract_this_amount);

uint safe_chromo_add_base_1
	(const uint &some_value__base_1,
	 const int &add_this_amount,
	 const uint &chromo_val);
		 	 	    
	    

BOOST_Interval get_padded_region_below_and_above_being_wary_of_chromosome_centromere_and_endpoints
                                  (const uint &chromosome_of_position_to_be_padded,
                                   const uint &position_to_be_padded_below_and_above,
                                   const uint &pad_below_amount,
                                   const uint &pad_above_amount);
                            	    
	    
				  
				  
				  
				  
type_BI__BI  convert_each_profile_breakpoint_to_absolute_breakpoints
				(const BOOST_Interval &profile_brkpts,
				 const type_vector_vector_vector_int &compressed_maps_LCRs_to_profile__of_some_Event);					
					
BOOST_Interval  convert_profile_breakpoint_to_absolute_breakpoints
				(const uint &profile_brkpt_index,
				 const type_vector_vector_vector_int &compressed_maps_LCRs_to_profile__of_some_Event);					
					 

				
type_BI__BI convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals
				(const BOOST_Interval &profile_brkpts,
				 const type_vector_vector_vector_int &compressed_maps_LCRs_to_profile__of_some_Event,
				const type_haploid_outcome &corresponding_hap_outcome);	
				
type_BI__BI convert_absolute_breakpoints_to_outcome_meaningful_intervals
		(const type_BI__BI &absolute_brkpts_each_LCR,
		 const type_haploid_outcome &corresponding_hap_outcome);
		
	    
	    
		
		
		
		
	    
//post-processing
		     
		     
		    // make outcome & breakpoint call using from marginal MAP      
void make_MAP_calls_from_marginal_posterior_distributions_per_Event
                                      (const type_map_uint_to_Marginal_Event_posterior_distribution  &marginal_distributions,
                                       const Conn_Comp &the_Conn_Comp,
                                       std::stringstream &output_marginal_Positive_calls_strm,
                                        std::stringstream &output_marginal_Negative_calls_strm,
                                        std::stringstream &output_marginal_all_calls_strm);   		     
		     
	
				 
void  calculate_GC_sensitive_per_position_HAPLOID_fragmentation_rate_of_Reference_genome_over_all_Readgroups___OMP
		    (type_map_uint_to_list_BI__string &pre_uploaded_regions_across_all_chromosomes,
		     const uint &chromosome_value,
		     type_map_uint_to_longdouble &positions_along_chromosome,
		     const uint &CC_ID,
		     const uint &EV_ID,
		     const uint &color_of_interval);				 
				 
				 
				 
				 
		     

type_map_haploid_outcome__haploid_outcome__to_uint_to_longdouble  construct_pseudo_RD_expected_Read_Depth_counts_for_desired_display_region
		    (const uint &chromo_of_region,
		     const BOOST_Interval &expanded_region_of_interest,
		    const type_haploid_outcome__haploid_outcome  &the_diploid_outcome,
		    const type_BI__BI &profile_breakpoints,
		    Event &the_Event);
		     
		     
		     
		     
		     
uint  calculate_NON_uniqueness_of_region
	    (const uint &chr_of_reg,
	     const BOOST_Interval &the_reg,
	     const type_map_uint_to_list_BI &NON_unique_regions_of_genome);
	     
	     
	     
	    
	    
					 
/////////////////////////////////////////////////////////////////////////////////////////////////////					 
// post-processing stuff						
// The following post-processing functions are for producing stats on Calls collected from runs from different genomes. 
//		i.e.:  first run "gather_calls" script to collect calls across different jobs.  Then run all of these things.
                         

// For comparing against validations                  

			


			 
			 
			 





uint  copy_contents_from_previous_output_directory_with_given_jobID_to_current_output_directory__and__return_last_CC_analyzed
			(const std::string  &previous_jobID);



			
			



	 
void check_informativeness_of_every_potential_NAHR_breakpoint
	    (type_list__uint__uint  &informativeness_checking_parameters,
	     const std::map<uint, Conn_Comp> &Conn_Comps,
	     const std::string &some_outputdir,
	     const std::string &BAM_dir);
	 
	 
	
	     
	     
	     
	     
 
			 
class aCGH_validation_container
{
    public:
	aCGH_validation_container()
	{ }
	
	aCGH_validation_container(
			const std::string &in_probe_intensity_filename,
			const uint &in_chr_of_call,
			const BOOST_Interval &in_brkpt_of_call,
			const uint &in_padding_amount,
			const type_set_string &in_individuals_with_call,
			const std::string &in_acceptable_genomes_filename,
			const std::string &in_save_prefix )
				:	probe_intensity_filename(in_probe_intensity_filename),
					chr_of_call(in_chr_of_call),
					brkpt_of_call(in_brkpt_of_call),
					padding_amount(in_padding_amount),
					individuals_with_call(in_individuals_with_call),
					acceptable_genomes_filename(in_acceptable_genomes_filename),
					save_prefix(in_save_prefix)
	{ }
	
	~aCGH_validation_container()
	{ }	
	
    
	std::string probe_intensity_filename;
	uint chr_of_call;
	BOOST_Interval brkpt_of_call;
	uint padding_amount;
	type_set_string individuals_with_call;
	std::string acceptable_genomes_filename;
	std::string save_prefix;  	
};

typedef  std::list<aCGH_validation_container>   type_list_aCGH_validation_container;


void validate_using_probe_intensity_comparison
	    (const aCGH_validation_container &some_aCGH_data);

		    
		    
uint  create_haploid_breakpoint_space
		(type_map_uint_to_vector_BI &brkpt_space_per_Event,
		 const Sparse_map &haploid_outcomes_per_Event,
                 std::map<uint,Breakpoint_complex> &effective_Breakpoint_complex);


void  set_breakpoints_using_breakpoint_space
	(type_map_uint_to_BI &haploid_breakpoint_instance,
	 const type_map_uint_to_vector_BI &haploid_breakpoint_space,
	 const uint &brkpt_ctr);

	 
	 
	 
uint  calculate_GC_content_of_sequence
	    (const std::string &some_sequence);	 
	 
	 
	
std::map<std::string, Readgroup_statistics>::const_iterator  identify_Readgroup_from_name_of_PER
			(const std::string &name_of_PER);	    
	    
	
			
mpfr::mpreal  calculate_poisson_pdf_or_cdf
			    (const longdouble &observation,
			    const longdouble &poisson_parameter,
			    const longdouble &pseudo_amount,
			    const bool &pdf__false__cdf__true);

mpfr::mpreal  calculate_negative_binomial_pdf_or_cdf__transforming_from_poisson_parameter
			    (const longdouble &observation,
			    const longdouble &poisson_parameter,
			    const longdouble &pseudo_amount,
			    const bool &pdf__false__cdf__true);
			
			    
// mpfr::mpreal  calculate_normal_pdf_or_cdf__as_approximation_to_negative_binomial_that_has_been_transformed_from_a_poisson
// 		(const longdouble &observation,
// 		const longdouble &poisson_parameter,
// 		 const longdouble &pseudo_amount,
// 		 const bool &pdf__false__cdf__true);
			
	
void  subdivide_intervals
	    (type_map_uint_to_list_BI  &map_chromo_to_intervals,
	     const uint &subdivison_size);
	     
	     
type_vector_BI split_interval_into_uniformly_sized_subintervals
		(const BOOST_Interval &big_interval,
		 const uint &number_of_subintervals,
		 const uint &minimum_sublength);   
	
type_map_uint_to_list_BI__string  upload_list_of_chromosome_regions_from_file
				    (const type_map_uint_to_list_BI  &map_chr_to_list_of_regions);	
	     
type_map_uint_to_list_BI  pad_all_regions
		(const type_map_uint_to_list_BI &some_map_chr_to_regions,
		 const uint &padding_amount);
				    
	
		 
		 
		 
		 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
BOOST_Interval  get_positions_of_hybrid_lying_BELOW_LCR
		    (const type_map_uint_to_uint  &absolute_positions_to_hybrid_indeces,
		     const BOOST_Interval  &some_LCR_coordinates);
	
	
BOOST_Interval  get_positions_of_hybrid_lying_ABOVE_LCR
		    (const type_map_uint_to_uint  &absolute_positions_to_hybrid_indeces,
		     const BOOST_Interval  &some_LCR_coordinates,
		     const uint &chromosome);		     
	
BOOST_Interval  get_positions_of_hybrid_lying_INBETWEEN_LCRs
		    (const type_map_uint_to_uint  &absolute_positions_to_hybrid_indeces,
		     const BOOST_Interval  &region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs);	
	
	
	
	
type_list_BI convert_set_to_list_of_contiguous_intervals
			(const type_set_uint &some_set_of_coords);	
	
	


				
longdouble get_adjustment_value_from_cumulative_map
			    (const BOOST_Interval  &affected_region,
			    const type_map_uint_to_longdouble &map_positions_to_cumulative_sum);
			
			    
			    
			    
bool  test_if_outcome_should_display__AB_ABA
			    (const type_haploid_outcome &some_event_outcome);

bool  test_if_outcome_should_display__BA_BAB
			    (const type_haploid_outcome &some_event_outcome);
			
	
bool test_if_diploid_state_vectors_contain_non_identifiable_diploid_outcomes
		(const Sparse_map &state_vector_hap_0,
		const Sparse_map &state_vector_hap_1);			    
			    
		
		
		
		
void set_prior_Event_outcome_probabilities
		    (const mpfr::mpreal &desired_probability_of_NO_outcome,
		    const mpfr::mpreal &fraction_of_all_occurring_amount_that_goes_to_GeneConversion_instead_of_NAHR);
		
		
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    

void   add_sparse_vp_to_form_Geneconversion_intervals
		(const Event &the_Event,
		const uint &maximum_number_of_allowable_breakpoints,
		const type_map_uint_to_uint &first_brkpt_to_count,
		const type_map_uint_to_uint &second_brkpt_to_count,
		type_set_BI  &set_of_GeneConv_breakpoints);
		    
		    
type_multimap_uint_to_uint::const_iterator  find_sparsest_vp_less_or_greater_than_value
						    (const type_multimap_uint_to_uint &multimap_sparseness_to_vp,
						    const uint &some_vp,
						    const bool &less_than__false____greater__true);		    
		    

type_multimap_uint_to_uint  get_sparseness_of_variational_positions
				    (const type_set_uint  &some_set_of_vp);	
				    
				    
void  preliminary_refine_breakpoints_for_inversion_heuristic
		(const Event &the_Event,
		const uint &maximum_number_of_allowable_breakpoints,
		const type_map_uint_to_uint &map_varpos_to_number_of_decent_reads__SAVE,
		type_set_uint  &refined_varpos,
		type_map_uint_to_uint  &refined_map_varpos_to_count);
		
		    
type_set_uint refine_breakpoints_using_heuristic
		(const Event &the_Event,
		const uint &maximum_number_of_allowable_breakpoints,
		const type_map_uint_to_uint &map_varpos_to_number_of_decent_reads__SAVE);
		
		
		
		
		
	
		
template <typename Tkey, typename Tval, typename Tcomparator_key = std::less<Tkey>,  typename Tcomparator_val = std::less<Tval> >	
void  prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained
		(std::set<Tkey, Tcomparator_key> &the_set,
		 const std::multimap<Tval, Tkey, Tcomparator_val> &map_count_to_elt,
		const uint &threshold_number_of_elements,
		std::map<Tkey, Tval, Tcomparator_key>  &some_corresponding_map_to_be_pruned_too)
{
    
    typename std::multimap<Tval,Tkey,Tcomparator_val>::const_iterator it_map_count_to_elt = map_count_to_elt.begin();
    
    while (the_set.size() > threshold_number_of_elements
	    and  it_map_count_to_elt != map_count_to_elt.end())	
    {
	the_set.erase(it_map_count_to_elt->second);
	some_corresponding_map_to_be_pruned_too.erase(it_map_count_to_elt->second);
	++it_map_count_to_elt;
    }
    
}//prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained


template <typename Tkey, typename Tval, typename Tcomparator_key = std::less<Tkey>,  typename Tcomparator_val = std::less<Tval> >	
void  prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained
		(std::set<Tkey, Tcomparator_key> &the_set,
		 const std::multimap<Tval, Tkey, Tcomparator_val> &map_count_to_elt,
		const uint &threshold_number_of_elements)
{
    
    typename std::multimap<Tval,Tkey,Tcomparator_val>::const_iterator it_map_count_to_elt = map_count_to_elt.begin();
    
    while (the_set.size() > threshold_number_of_elements
	    and  it_map_count_to_elt != map_count_to_elt.end())	
    {
	the_set.erase(it_map_count_to_elt->second);
	++it_map_count_to_elt;
    }
    
}//prune_low_count_elements_from_set_using_count_map_until_threshold_is_attained
	
	
	
	
	
	
	
	
type_set_BI  group_breakpoints_according_to_relevant_variational_position_constraints
					    (const type_set_uint &relevant_variational_positions,
					    const type_set_BI &the_space_of_potential_breakpoints);
			    
			    
	
void determine_relevant_NAHR_and_GeneConversion_breakpoints_for_each_interlocking_event_and_self_according_to_relevant_variational_positions
                                                        (const Conn_Comp &the_Conn_Comp,
							 const type_map_uint_to_set_uint &relevant_varpos_per_Event,
							std::map<uint,Breakpoint_complex> &relevant_breakpoints_per_Event);
							
	
							
real  calculate_standard_normal_tail
		(const longdouble &standardized_Z_score);	
		
real  calculate_normal_pdf_or_cdf
		(const longdouble &observation,
		const longdouble &mean,
		 const longdouble &sigma,
		 const bool &pdf__false__cdf__true);		
	
		
		
		
		
		
		
		
void subtract_or_add_map_values_from_another_map__within_affected_region
			(const type_map_uint_to_longdouble &map_pos_to_adjustment_value,
			 const BOOST_Interval &affected_region,
			const bool &subtract__false____add_true,
			type_map_uint_to_longdouble &updated_map,
			type_map_uint_to_int  &contributions_per_position);		
	
			
void  adjust_expected_fragmentation_rate_according_to_Event_outcome_and_breakpoints__on_repeat_and_raw_region
				(const Event &the_Event,
				 const type_haploid_outcome &the_haploid_outcome,				    
				const BOOST_Interval &the_haploid_breakpoints,
				const type_map_uint_to_uint_to_longdouble  &map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
				const type_map_uint_to_longdouble &raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
				type_map_uint_to_longdouble &updated_repeat_aware_expected_fragmentation_rate_per_positon,
				type_map_uint_to_int &contributions_per_position_to_updated_expected_fragmentation_rate);

				
void  raw_add_expected_fragmentation_rate_of_homologous_regions__raw_haploid_fragmentation_rate(const Event &the_Event,
				 const type_map_uint_to_uint_to_longdouble  &map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
				type_map_uint_to_longdouble &updated_repeat_aware_expected_fragmentation_rate_per_positon,
				type_map_uint_to_int &contributions_per_position_to_updated_expected_fragmentation_rate);
				
				
void  adjust_expected_fragmentation_rate_by_homologous_regions_raw_haploid_fragmentation_rate
			    (const Event &the_Event,
			    const type_map_uint_to_uint_to_longdouble  &map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
			    const BOOST_Interval &affected_region_on_chromosome_of_Event,
			    const bool &subtract_false____add_true,
			    type_map_uint_to_longdouble &updated_repeat_aware_expected_fragmentation_rate_per_positon,
			    type_map_uint_to_int &contributions_per_position_to_updated_expected_fragmentation_rate);
			    

void appropriately_repeat_aware_adjust_expected_fragmentation_rate_and_number_of_contributions_according_to_haploid_rate_and_event_outcome
				(const Event &the_Event,
				 const type_haploid_outcome &the_haploid_outcome,
				 const BOOST_Interval &the_haploid_breakpoints,
				const type_map_uint_to_longdouble &raw_haploid_fragmentation_rate_per_position_of_expanded_region_of_interest,
				const BOOST_Interval &expanded_region_of_interest,
				const type_map_uint_to_uint_to_longdouble  &map_chromo_to_homol_reg_to_haploid_raw_fragmentation_rate,
				type_map_uint_to_longdouble &updated_repeat_aware_expected_fragmentation_rate_per_positon,
				type_map_uint_to_int &contributions_per_position_to_updated_expected_fragmentation_rate);			    
				
				
				
type_map_uint_to_longdouble  divide_sum_by_contributions
					(const type_map_uint_to_longdouble &the_map_of_sums,
					const type_map_uint_to_int &the_map_of_number_ofcontributions_to_the_sums);

	
					
					
					
void redo_all_marginal_MAP_calls_using_saved_samples_from_posterior	
		(const std::string &basedir,
		const std::map<uint, Conn_Comp> &Conn_Comps);
		

		
longdouble get_standard_deviation_from_mean_for_RD
				(const longdouble &mean);		
	
				
void read_BAS_file_and_GC_rates_and_set_fragment_lengths_and_lambda
		    (const std::string &BAM_filename,
		    const bool &display_info,
		    BamTools::BamReader *const &my_BAM_reader_ptr);
	
void redo_all_inferred_RD_tests_and_make_ratios    
				(const std::string &postdatatdir,
				 std::map<uint, Conn_Comp> &Conn_Comps);		    
		    
	
type_set_uint get_nonunique_positions_in_interval
		    (const type_map_uint_to_list_BI &non_uniq_regions_of_genome,
		     const uint &chromo,
		     const BOOST_Interval &expanded_region_of_interest);				
			
			
			
			
			
			
typedef std::map<BOOST_Interval, longdouble>  type_map_BI_to_longdouble;

type_map_BI_to_longdouble bin_by_mapykey
		(const type_map_uint_to_longdouble &some_map,
		const uint &binsize);			
	

		
		
real  perform_Wilcoxon_signed_rank_test
				(const uint &binsize,
				const type_map_uint_to_longdouble &map_position_to_expected_RD___NULL,
				const type_map_uint_to_longdouble &map_position_to_Observed_RD,
				const Sampled_diploid_Event_data &some_result);		


void  do_signed_rank_test_for_Read_Depth
				(const std::string &postdatatdir,
				 std::map<uint, Conn_Comp> &Conn_Comps);
				
				

uint randomly_choose_bin_of_same_GC_count_and_count_observed
			(const uint &chr_of_current_bin,
			const BOOST_Interval &region_of_current_bin,
			const type_map_uint_to_list_BI &NON_unique_regions_of_genome,
			BamTools::BamReader &my_BAM_reader);
				
				
	
real  perform_Wilcoxon_signed_rank_test
				(const uint &binsize,
				 const uint &chromo_of_observed,
				const type_map_uint_to_longdouble &map_position_to_Observed_RD,
				const Sampled_diploid_Event_data &some_result,
				BamTools::BamReader &my_BAM_reader);
	
				
				


typedef  std::map<uint, type_list_uint>   type_map_uint_to_list_uint;

type_map_uint_to_list_uint get_random_sample_of_bins_from_across_the_genome
								(const uint &bin_width,
								const uint &number_of_samples, //mult of 10
								BamTools::BamReader &my_BAM_reader);				
				
								
int randomly_choose_matching_control_from_subsample
			(const uint &desired_GC_count,
			 const type_map_uint_to_list_uint &sampled_control_bins);								

			
			
			
			
			
class Wilcoxon_calculations
{
	public:
		Wilcoxon_calculations()
						:  wilcoxon_pvalue(-1),
							wilcoxon_likelihood(-1),
							signed_rank_sum_value(0),
							number_of_elements_in_ranking(0),
							variance_of_normal(0)		
		{  }
		
		~Wilcoxon_calculations()
		{  }
		
		
		// variables:
		real wilcoxon_pvalue;
		real wilcoxon_likelihood;
		longdouble signed_rank_sum_value;
		uint number_of_elements_in_ranking;
		uint variance_of_normal;
	
}; // Wilcoxon_calculations

			
			


			
void  perform_Wilcoxon_signed_rank_test
				(const uint &binsize,
				 const uint &chromo_of_observed,
				const type_map_uint_to_longdouble &map_position_to_Observed_RD,
				const Sampled_diploid_Event_data &some_result,
				 const type_map_uint_to_list_uint &sampled_control_bins_by_GC_count,
				const longdouble &sampled_control_multiplier,
				const type_map_uint_to_list_BI  &nonunique_regions_of_genome,
				Wilcoxon_calculations &observed_wilcoxon_stats___inferred,
				Wilcoxon_calculations &observed_wilcoxon_stats__unique);
				
				

				
				
void compute_signed_rank_sum
				(const type_map_BI_to_longdouble &binned_difference,
				 const bool &verbose,
				 longdouble &observed_signed_ranksum,
				uint &number_of_elts_in_ranking);
				
				
int get_uniform_random_draw
			(const int &max_num);
			
	
// template <typename Telt1, typename Telt2>
// void randomly_shuffle_sorted_list_elements_with_equivalent_first_element_in_pair
// 				(std::list< std::pair<Telt1, Telt2> > &the_list)
// {
// 	std::vector< std::pair<Telt1, Telt2> > the_vector(convert_list_to_vector< std::pair<Telt1, Telt2> >(the_list));
// 	
// 	
// 	typename std::vector< std::pair<Telt1, Telt2> >::iterator it_vec = the_vector.begin();		
// 	while (it_vec != the_vector.end())
// 	{
// 		typename std::vector< std::pair<Telt1, Telt2> >::iterator it_next = it_vec; 
// 		++it_next;
// 		int num_equivs = 1;
// 		
// 		while (it_next != the_vector.end()  and  it_next->first == it_vec->first)
// 		{			
// 			++it_next;									
// 			++num_equivs;
// 		}//same and not end		
// 		
// 		if (num_equivs > 1)
// 		{			
// 			std::random_shuffle(it_vec, it_next, get_uniform_random_draw);
// 			std::random_shuffle(it_vec, it_next, get_uniform_random_draw);			
// 		}		
// 	}//end
// 	
// 	
// 	the_list = convert_vector_to_list< std::pair<Telt1, Telt2> >(the_vector);
// 	
// 	
// }//randomly_shuffle_sorted_list_elements_with_equivalent_first_element_in_pair
			
		
	
void  bernoulli_lyapunov_nonparametric_test
				(const std::string &postdatatdir,
				 std::map<uint, Conn_Comp> &Conn_Comps);	
	
				
		
				
class Lyap_CLT_data
{
	public:
		Lyap_CLT_data()
				: sum_log_P_obs(0),
				sum_mu(0),
				sn_2(0)
		{  }
		
		~Lyap_CLT_data()
		{  }
		
		//variables:
		longdouble sum_log_P_obs;
		longdouble sum_mu;
		longdouble sn_2;
	
}; //Lyap_CLT_data
		
		
		

Lyap_CLT_data calculate_Lyapunov_CLT_sum
			(const type_vector_longdouble &bernoulli_p_0_by_position__null_model,
			 const type_vector_vector_longdouble &poisson_lambda_summed_over_readgroup_per_position,
			const type_vector_uint &bernoulli_observations_by_position,
			const type_vector_bool &is_unique = type_vector_bool());
		
		
	
bool check_if_position_is_in_NON_unique_region
	(const type_map_uint_to_list_BI &nonuniq_regions,
	const uint &chromo,
	const uint &pos_on_chromo);
	
bool check_if_region_overlaps_NON_unique_region
	(const type_map_uint_to_list_BI &nonuniq_regions,
	const uint &chromo,
	const BOOST_Interval &reg_on_chromo);	
			
	
	
// void sample_control_positions_from_across_genome
// 					(const uint &number_of_samples,
// 					 const type_map_uint_to_list_BI &NON_unique_regions_of_genome,
// 					BamTools::BamReader &my_BAM_reader,
// 					type_map_uint_to_list_uint &sampled_observed_counts_per_GC_content,
// 					uint &mean_frag_length);
// 	
// 	
// void sample_control_positions_for_all_genomes
// 				(const std::string &postdatatdir);
					
					
	
real compute_normalized_statistic
		(const Lyap_CLT_data &null_data, 
		 const Lyap_CLT_data &alt_data);	
	
real compute_lower_tail_probability
		(const Lyap_CLT_data &null_data, 
		 const Lyap_CLT_data &alt_data);	
		
	
void  raw_count_ratios
		(const std::string &postdatatdir,
		std::map<uint, Conn_Comp> &Conn_Comps);		
		
	
		
		
void determine_cumulative_GC_count_for_each_chromosome_and_save_to_file();
		

				
const type_uint__uint convert_random_draw_to_chromo_and_position
		(const uint &whole_genome_position_draw,
		 const type_map_BI_to_uint &map_cumulative_genome_position_to_chromo);				
				
	
bool check_if_string_contains_N
		(const std::string &some_seq);		
		
void  randomly_sample_contiguous_region_from_control_genome
				(const type_map_uint_to_list_BI &nonunique_regions_of_genome,
				const type_map_uint_to_set_uint &learning_positions_by_chromo,
				const type_map_BI_to_uint &map_cumulative_genome_position_to_chromo,
				const uint &desired_region_length,
				const type_map_uint_to_string &map_chromo_to_seq,
				type_map_uint_to_set_BI &sampled_test_regions);

void randomly_sample_GC_sensitive_positions_from_control_genome
			(const uint &number_positions_for_sample,
			const type_map_uint_to_list_BI &nonunique_regions_of_genome,
			const type_map_uint_to_set_uint &learning_positions_by_chromo,
			const type_map_BI_to_uint &map_cumulative_genome_position_to_chromo,
			const type_map_uint_to_uint &map_frag_lengths_to_desired_GC,
			const type_map_uint_to_string &map_chromo_to_seq,
			type_map_uint_to_set_uint &sampled_test_positions);				
		
type_map_BI_to_uint  make_map_cumulative_genome_position_to_chromosome();		
		
		
uint count_decent_reads_with_left_endpoints_in_region__readgroup_specific
					(BamTools::BamReader &a_BAM_reader,
					const uint &chromosome,
					const BOOST_Interval &chromo_interval);
	
	
	
class Count_observation_data
{
	public:
		Count_observation_data()
		{  }
		
		Count_observation_data
				(const uint &some_chromosome,
				const uint &some_reg,
				const type_map_string_to_uint &some_counts)
					:  chromosome(some_chromosome),
						region(some_reg),
						map_readgroup_to_count(some_counts)
		{  }		
		
		
		//variables:
		uint chromosome;
		BOOST_Interval region;
		type_map_string_to_uint  map_readgroup_to_count;
		
		//functions:
		void absorb_counts
				(const type_map_string_to_uint &some_count_observations);
	
	
};//Count_observation_data

typedef  std::list<Count_observation_data>  type_list_Count_observation_data;




void  draw_count_ratio_using_sampling_strategy_A
			(const std::map<uint, Conn_Comp> &Conn_Comps,
			BamTools::BamReader &a_BAM_reader);
	
void  draw_count_ratio_using_sampling_strategy_B
			(const std::map<uint, Conn_Comp> &Conn_Comps,
			 BamTools::BamReader &a_BAM_reader);
			
	
std::string infer_genome_name_from_BAM_filename
				(const std::string &some_bam_filename);			
			
			
	
				
				
				
				
				
				
				
				
class Changed_region
{
	public:
		uint chr;
		BOOST_Interval region;
		type_haploid_outcome rearrangement;	
		uint cc;
		uint ev;
};

typedef  std::list<Changed_region>  type_list_Changed_region;


bool compare_Changed_regions
				(const Changed_region &change1, Changed_region &change2);
			
	
				
				

				
				
uint get_CC_of_EV
			(const std::map<uint, Conn_Comp>  &Conn_Comps,
			const uint &ev_id);
			
			

void simulate_NAHR_events_on_Reference_genome
				(const std::map<uint, Conn_Comp>  &Conn_Comps,
				const std::string &sim_gen_name,
				const type_vector_uint &allowable_events,
				const uint &num_events_to_simulate_per_genome);			
				
				
				
				
				
			
#endif