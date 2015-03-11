#ifndef POST_H_INCLUDED
#define POST_H_INCLUDED


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
#include <boost/graph/graph_concepts.hpp>


#include <mpreal.h>
#include <mpfr.h>



#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamAux.h>



#include <general_typedefs.h>

#include "Event.h"



class Conn_Comp;
class Call_set;
class Event;
class MiniRegion;
class Visitation_object;
class Paired_end_read;
class Readgroup_statistics;
class MiniEvent;
class Breakpoint_complex;

class Gene;
class Affected_Gene_catalog;
           
class Wilcoxon_calculations;


typedef  std::vector<type_list_Validated_rearrangement::const_iterator>  type_vector_list_Validated_rearrangement_constiterator;


typedef  std::map<type_haploid_outcome__haploid_outcome, type_list_real>    type_map_diploid_outcome_to_list_real;



	



    
//classes:			    
			    

		 
class Gene
{
  public:    
    Gene()
    {  }
    
    Gene(const std::string &in_name,
	 const std::string &in_ensembl_id,
	 const uint &in_chr,
	 const BOOST_Interval &in_coords,
	 const std::string &in_desc)  
			    :	gene_name(in_name),
				ensembl_id(in_ensembl_id),
				chromosome(in_chr),
				coordinates(in_coords),
				description(in_desc)
    {  }    
    
    ~Gene()
    {  }
    
    
    //variables
    std::string gene_name;    
    std::string ensembl_id;
    
    uint chromosome;
    BOOST_Interval coordinates;
    
    std::string description;
    
    
    bool operator< (const Gene &another_Gene)  const;
    
        
};//Gene
		 		 	
typedef   std::map<std::string, Gene>  type_map_string_to_Gene;
typedef   std::map<BOOST_Interval, Gene, compare_BI>  type_map_BI_to_Gene;
		    
typedef   std::map<uint, type_map_BI_to_Gene>  type_map_uint_to_BI_to_Gene;
	    
typedef  std::set<Gene>  type_set_Gene;			    
			    
typedef  std::pair<type_set_Gene, type_set_Gene>  type_set_Gene__2;


typedef  std::map<std::string, type_set_Gene>  type_map_string_to_set_Gene;








			    
			    
typedef  std::map<std::string, Affected_Gene_catalog>  type_map_string_to_Affected_Gene_catalog;			    

class Affected_Gene_catalog
{
    public:
	Affected_Gene_catalog()
	{  }        
	
	~Affected_Gene_catalog()
	{  }
	
	
	//variables:
	    //"ensembl_id" to Gene
	type_map_string_to_Gene Genes_affected_by_Inversions;
	type_map_string_to_Gene Genes_affected_by_Deletions;
	type_map_string_to_Gene Genes_affected_by_Duplications;
	type_map_string_to_Gene Genes_affected_by_GeneConversion__donated_itself_to_elswhere;
	type_map_string_to_Gene Genes_affected_by_GeneConversion__received_from_somewhere_else;
	
	type_set_uint Events_occurring_as_Inversions;
	type_set_uint Events_occurring_as_Deletions;
	type_set_uint Events_occurring_as_Duplications;
	type_set_uint Events_occurring_as_GeneConversions;
	
	
	//functions:
	type_map_string_to_Gene  get_all_affected_Genes()  const;
	
	void collect_from_individual_catalogs
			    (const type_map_string_to_Affected_Gene_catalog &catalogs_per_individual);
    
    
};//Affected_Gene_catalog



class Affected_Gene_tally
{
    public:
	Affected_Gene_tally()
			: num_invs(0),
			  num_dels(0),
			  num_dups(0),
			  num_GeneConv__donor(0),
			  num_GeneConv__receiver(0)
	{  }
	
	~Affected_Gene_tally()
	{  }
	
	
	//variables:
	std::string affected_gene__ensembl_ID;
	
	uint num_invs;
	uint num_dels;
	uint num_dups;
	uint num_GeneConv__donor;
	uint num_GeneConv__receiver;		
	
	
	void add_call_appropriately
			(const type_haploid_outcome &some_event_outcome);
    
};//Affected_Gene_tally


typedef std::map<std::string, Affected_Gene_tally>  type_map_string_to_Affected_Gene_tally;
			    
			    
			    
			    
			    
			    
			    
			    
			    
typedef   std::pair<type_set_string, type_set_string>  type_set_string__2;
		
		
class Called_diploid_wrapper : public Sampled_diploid_Event_data
{
    public:
	Called_diploid_wrapper() : Sampled_diploid_Event_data(), set_RD_odds(false), Wilcoxon_signed_rank_sum_test__pvalue__inferred(-5), raw_count_ratio__inferred(-71), raw_count_ratio__unique(-71), brkpt_odds(-1,-1)
	{  }
	
	Called_diploid_wrapper(const Sampled_diploid_Event_data &some_sample) : Sampled_diploid_Event_data(some_sample), set_RD_odds(false), Wilcoxon_signed_rank_sum_test__pvalue__inferred(-5), raw_count_ratio__inferred(-71), raw_count_ratio__unique(-71), brkpt_odds(-1,-1)
	{  }	
	
	~Called_diploid_wrapper()
	{  }
	
	
	//variables:
	type_set_string associated_other_calls;
	std::string population_name;
	
	type_set_string__2 affected_genes_per_haploid;
	
	uint number_of_unique_bases_in_region;
	
	bool from_a_CC_with_undefined_regions;
	
	real RD_odds_alt;
	real RD_odds_del;
	real RD_odds_dup;
	real RD_odds_best;
	
	bool set_RD_odds;
	
	real Wilcoxon_signed_rank_sum_test__pvalue__inferred;
	real Wilcoxon_signed_rank_sum_test__pvalue__unique;
	
	Wilcoxon_calculations  wilcoxon_matched_info___null___inferred;
	Wilcoxon_calculations  wilcoxon_matched_info___null___unique;
	Wilcoxon_calculations  wilcoxon_matched_info___alternative___inferred;
	Wilcoxon_calculations  wilcoxon_matched_info___alternative___unique;

	real wilcoxon_likelihood_ratio___inferred;
	real wilcoxon_likelihood_ratio___unique;
	real wilcoxon_statistic_ratio___inferred;
	real wilcoxon_statistic_ratio___unique;
	real wilcoxon_Cauchy_cdf___inferred;
	real wilcoxon_Cauchy_cdf___unique;
	
	real wilcoxon_Lyapunov_CLT___inferred;
	real wilcoxon_Lyapunov_CLT___unique;
	
	real raw_count_ratio__inferred;
	real raw_count_ratio__unique;
	
	longdouble logfdr_inferred;
	
	type_real__real  brkpt_odds;
	
// 	real Wilcoxon_MATCHED_signed_rank_sum_test__pvalue__inferred;
// 	real Wilcoxon_MATCHED_signed_rank_sum_test__pvalue__unique;
    
};//Called_diploid_wrapper

typedef  std::map<uint,Called_diploid_wrapper>  type_map_uint_to_Called_diploid_wrapper;
typedef  std::map<std::string, type_map_uint_to_Called_diploid_wrapper>  type_map_string_to_uint_to_Called_diploid_wrapper;

typedef  std::list<const Called_diploid_wrapper*>  type_list_Called_diploid_wrapper_ptr;
typedef  std::map<std::string, type_list_Called_diploid_wrapper_ptr>  type_map_string_to_list_Called_diploid_wrapper_ptr;
typedef  std::map<std::string, type_map_string_to_list_Called_diploid_wrapper_ptr>  type_map_string_to_string_to_list_Called_diploid_wrapper_ptr;






		
class Empirical_properties
{
    public:
	Empirical_properties()			
	{  }
	
	~Empirical_properties()
	{  }
	
	
	//static variables
	static uint MEPS_size;
	static uint informativeness_fragment_length;
	static uint informativeness_mate_size;
	
	
	//variables:	
	type_list_uint sparsity_histogram_of_breakpoints; // key = sparseness, value = count	
	
	type_list_longdouble  mediating_LCR_overall_percent_identity_histogram;
	type_list_uint  mediating_LCR_length_histogram;		
	type_list_uint distance_between_mediating_LCRs_histogram;
	
	type_list_longdouble best_MEPS_homology_histogram;		

	type_list_uint  affected_region_size_histogram; // NAHR only
	
	type_list_uint  length_of_donor_GeneConversion_tract__histogram;
	type_list_uint  length_of_receiving_GeneConversion_tract__histogram;
	type_list_longdouble  ratioGC_content_of_GeneConversion_tract__histogram;	
	
	type_list_uint  number_of_occurring_individuals__vs__number_of_occurring_Events__histogram;
				    //[key] = number of individuals
				    //[value] = number of Events which occurr in [key] individuals
	
	type_list_longdouble posterior_probabilities_of_MAP_event_outcome__histogram;
	type_list_longdouble posterior_probabilities_of_MAP_event_breakpoints__histogram;
	
	//functions:			    
	void absorb_other_Empirical_properties
		    (const Empirical_properties &other_EP);	
		    
	void save_to_file
		    (const std::string &outdir,
		    const std::string &EP_name)   const;
		    
    
};//Empirical_properties
		
		
		
		
		
		
class Summary_statistics
{
    public:
	Summary_statistics()
			    :  total_DNA_affected_by_NAHR___diploid(0),
			    total_DNA_affected_by_GeneConversion___diploid(0),
			    total_NAHR_calls__diploid(0),
			    total_GeneConversion_calls__diploid(0)    
	{  }
	
	~Summary_statistics()
	{  }
	
	
	
        uint total_DNA_affected_by_NAHR___diploid;	
	uint total_DNA_affected_by_GeneConversion___diploid;		
	
	uint total_NAHR_calls__diploid;
	uint total_GeneConversion_calls__diploid;
	
	
	type_map_haploid_outcome_to_uint  map_haploid_outcome_to_counts;
	type_map_haploid_outcome__haploid_outcome_to_uint  map_DIPLOID_outcome_to_counts;	
	
	type_set_uint connected_components_with_some_NAHR_called_positive__haploid;
	type_set_uint connected_components_with_some_GeneConversion_called_positive__haploid;
	
	
	void add_call_appropriately
			    (const std::map<uint, Conn_Comp>  &the_CCs,
			     const Sampled_diploid_Event_data &call);	

    
};//Summary_statistics
		
	
typedef  std::map<std::string, Summary_statistics>  type_map_string_to_Summary_Statistics;	








class  Theoretical_compute_space
{
    public:
	Theoretical_compute_space()
			    :  number_of_deletions_tested__diploid(0),
			    number_of_duplications_tested__diploid(0),
			    number_of_inversions_tested__diploid(0),
			    number_of_GeneConversion_tested__diploid(0),
			    number_of_Connected_Components_tested(0),
			    number_of_Events_tested_diploid(0),
			    amount_of_haploid_Reference_genome_examined_for_potential_events(0)
	{  }
	
	~Theoretical_compute_space()
	{  }
	
	
	
	//variables:
	
	uint number_of_deletions_tested__diploid;
	uint number_of_duplications_tested__diploid;
	uint number_of_inversions_tested__diploid;
	uint number_of_GeneConversion_tested__diploid;
	
	uint number_of_Connected_Components_tested;
	uint number_of_Events_tested_diploid;
	
	uint amount_of_haploid_Reference_genome_examined_for_potential_events;    
	
	type_set_string  potentially_affected_Genes;
	
    
};//Theoretical_compute_space













	

			 
			 
			 
			 
//helper functions/classes for post-processing:

bool some_string_in_set_is_a_substring_of_given_string
		(const type_set_string &set_of_strings,
		 const std::string &given_string);	 


bool compare_calls_by_posterior_prob_then_brkpt_prob
			(const Call_set &LHS, 
			 const Call_set &RHS);
			 
			 
                         
void  write_line_of__Call_vs_Val__of_table__to_ostream
                (const Call_set *const &ptr_to_call,
                const Validated_rearrangement *const &ptr_to_val,
                  std::ostream &some_ostream);                     
			 
			
		

void associate_Events_with_Validated_rearrangements
			    (const type_map_uint_to_CC  &Conn_Comps,
			     type_list_Validated_rearrangement &validations_from_other_studies);
		
		
			    
			    
			
			    
			    
	
		
type_map_string_to_uint_to_Called_diploid_wrapper wrap_Sampled_diploid_Event_data
							    (const type_map_string_to_uint_to_Sampled_diploid_Event_data &some_samples);
		
		






void post_process
	    (std::map<uint, Conn_Comp> &Conn_Comps,
	     const std::string &postdatadir,      
	    const std::string &savedir,
	     const std::string &genes_filename,
	    const std::string &biological_validations_filename,
		const std::string &gender_and_population_filename);





		    
			
			
	
// type_map_string_to_Gene identify_Genes_affected_by_Event_breakpoints
// 					(const Event &occurring_Event,
// 					const type_haploid_outcome &haploid_outcome,
// 					const BOOST_Interval &haploid_brkpts,
// 					const type_map_uint_to_BI_to_Gene &map_chromo_to_region_to_Gene);
				
				
 
	    
	    
	    
	    
	    
	    
//helper functions.	    			    

void  match_validations_with_Calls
		    (type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		     const type_list_Validated_rearrangement &validations_from_other_studies,
		    type_vector_list_Validated_rearrangement_constiterator   &Vals__correctly_matched,
		    type_vector_list_Validated_rearrangement_constiterator   &Vals__alternatively_matched,
		    type_vector_list_Validated_rearrangement_constiterator   &Vals__unmatched,
		    type_map_string_to_string_to_list_Called_diploid_wrapper_ptr &map_genome_to_valname_to_calls__correctly_matched,
		    type_map_string_to_string_to_list_Called_diploid_wrapper_ptr &map_genome_to_valname_to_calls__alternatively_matched,
		    type_map_string_to_uint_to_Called_diploid_wrapper  &novel_mycalls);
		    
	    
type_map_string_to_Gene identify_Genes_potentially_affected_by_an_Event
					(const Event &occurring_Event,
					const type_map_uint_to_BI_to_Gene &map_chromo_to_region_to_Gene);				
				
					
					
					

void determine_Genes_affected_by_each_Event_outcome_on_each_individual
		(const type_map_uint_to_CC &Conn_Comps,
		 const type_map_uint_to_BI_to_Gene &all_Genes,
		const real &logodds_threshold,
		const int &pvalue_for_RD_test_choice,
		type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		type_map_string_to_Affected_Gene_catalog &affected_Genes_by_individual);

	    


void fill_Empirical_properties
	    (const std::map<uint,Conn_Comp> &Conn_Comps,
	     const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
	    const type_map_uint_to_list_BI__string &map_chromo_to_uploaded_regions,
	    Empirical_properties &EP_for_NAHR,
	    Empirical_properties &EP_for_Gene_Conversion,
	    const bool &fill_compute_intensive_properties_too);	    
		
	    
	    
void compile_Summary_statistics
		    (const type_map_uint_to_CC &Conn_Comps,
		     const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		    const real &logodds_threshold,
			const int &pvalue_for_RD_test_choice,
		    type_map_string_to_Summary_Statistics &summary_stats_per_individual,
		    Summary_statistics &summary_stats_across_ALL_individuals);
	   
	    
	    
Theoretical_compute_space compile_Theoretical_compute_space
	    (const type_map_uint_to_CC &Conn_Comps,
	    const type_map_uint_to_BI_to_Gene &Genes_in_the_human_genome);    
	    
	    
	    
	    
// type_map_diploid_outcome_to_list_real collect_RD_hypothesis_test_p_values_for_FDR
// 						    (const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);	    
	    
	    
	    
	    
	    
//helper sub-functions	

	    

void load_RD_hypothesis_test_pvalues
		    (const std::string &postdatadir,
		    type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);





void restrict_to_allowable_genomes
		(type_list_Validated_rearrangement  &calls_from_other_studies,
		const type_set_string &allowable_genomes_names);
	    
	    
	    


void evaluate_impact_of_Call_on_Gene_appropriately
		    (const Event &the_Event,
		     Called_diploid_wrapper &the_call,
		     const uint &hap,
		    const Gene &the_affected_Gene,
		    type_map_string_to_Affected_Gene_catalog  &affected_Genes_per_individual);									    	
		    
	    
	    
uint calculate_GC_content_of_region
		(const uint &chromosome,
		const BOOST_Interval &region,
		const type_map_uint_to_list_BI__string &map_chromo_to_uploaded_regions);			    
			    
	    	    
	    
	    
	    

		    
	
typedef  std::map<uint, type_map_haploid_outcome_to_uint>  type_map_uint_to_haploid_outcome_to_uint;
		    
	    
	    
	    

longdouble determine_best_MEPS_homology_of_profile_breakpoint
		    (const Event &the_Event,
		    const uint &the_profile_breakpoint,
		    const uint &MEPS_size_to_consider);


longdouble determine_LCR_percent_identity__ignoring_large_gaps
			    (const Event &the_Event,
			    const uint &min_gap_size_to_be_considered_large_gap_for_identity_calculation);		    
		    

	    
   
	
			    
    
			    
			    
			    
			    
void write_line_of_Summary_stats
	    (std::ofstream &outfs,
	    const Theoretical_compute_space &the_theoretical_space,
	    const std::string &name_of_line,
	    const Summary_statistics &some_Summary_statistics,
	    const Affected_Gene_catalog &corresponding_Affected_genes,
	    const uint &average_divider = 1);
			    
			    
			  
	    
	    


			    
			    
			    
			    
void write_raw_NAHR_calls_to_stream
	(const std::string &outfname,
	 const type_map_uint_to_CC &Conn_Comps,
	 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
	const type_haploid_outcome &only_this_haploid_outcome);
	

void write_raw_GeneConversion_calls_to_stream
	(const std::string &outfname,
	 const type_map_uint_to_CC &Conn_Comps,
	 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);	
			    
			    
			    
			    
type_map_uint_to_list_BI__string  upload_affected_regions_sequences
				    (const type_map_uint_to_CC &Conn_Comps,
				     const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);    
			    
			    
				    
void fill_uniqueness_and_undefinedness_of_each_call
			    (const type_map_uint_to_CC &Conn_Comps,
			     type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);
			    
			    
			    
	
std::map<uint, std::map<BOOST_Interval, Gene, compare_BI> > read_Genes_from_file
								    (const std::string &gene_fname);
								    
								    
bool check_that_rearrangement_breakpoints_are_sufficiently_close_to_Event_LCRs
	(const BOOST_Interval &rearrangement_breakpoints,
	 const BOOST_Interval *const &LCRs_of_some_Event,
	 const uint &leniancy);
	 
bool some_string_in_set_is_a_substring_of_given_string
		(const type_set_string &set_of_strings,
		 const std::string &given_string);									    
	
				    


void  remove_calls_in_undefined_regions_of_genome
		(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls);		
		
void  remove_calls_without_RD
		(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls);	
		
		

type_map_string_to_uint_to_Called_diploid_wrapper  restrict_to_NAHR_using_RD_threshold
				    (const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
				    const real &logodds_threshold,
					const int &pvalue_for_RD_test_choice);
					
					
		
			    
//save


void save_Theoretical_space_to_table
			    (const std::string &outdir,
			    const Theoretical_compute_space &the_theoretical_space);    
			    
			    
void save_Summary_statistics_to_table
			    (const std::string &outdir,
			     const std::string &addl_id,
			    const Theoretical_compute_space &the_theoretical_space,
			    const type_map_string_to_Summary_Statistics &summary_stats_per_individual,
			    const Summary_statistics &summary_stats_across_ALL_individuals,
			    const type_map_string_to_Affected_Gene_catalog  &affected_Genes_by_individual);
			    	
			    
			    
void save_raw_positive_calls_to_files
	(const std::string &outdir,
	 const std::string &addl_id,
	 const type_map_uint_to_CC &Conn_Comps,
	const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);
	
	    
void make_and_save_RD_FDR_distributions
	    (const std::string &outdir,
	    const type_map_diploid_outcome_to_list_real &RD_hypothesis_p_value_by_diploid_outcome);    
	
	    
	    
void  save_Affected_Genes
	(const std::string &outdir,
	 const type_map_string_to_Affected_Gene_catalog &affected_Genes_by_individual,
	const type_map_uint_to_BI_to_Gene &all_Genes_by_region);
	    
	
	

	
void save_comparisons_to_tables
		(const std::string &outdir,
		 const real &logodds_threshold,
		const int &pvalue_for_RD_test_choice,
		const type_map_uint_to_CC &Conn_Comps,
		const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		const type_list_Validated_rearrangement &validations_from_other_studies,
		const type_vector_list_Validated_rearrangement_constiterator   &Vals__unmatched,
		const type_map_string_to_string_to_list_Called_diploid_wrapper_ptr &map_genome_to_valname_to_calls__correctly_matched,
		const type_map_string_to_string_to_list_Called_diploid_wrapper_ptr &map_genome_to_valname_to_calls__alternatively_matched,
		const type_map_string_to_uint_to_Called_diploid_wrapper  &novel_mycalls);
	
void save_comparisons_to_tables___imp
		(const std::string &outdir,
		 const std::string &addl_id,
		const real &logodds_threshold,
		const int &pvalue_for_RD_test_choice,
		const type_map_uint_to_CC &Conn_Comps,
		const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		const type_list_Validated_rearrangement &validations_from_other_studies,
		const type_vector_list_Validated_rearrangement_constiterator   &Vals__unmatched,
		const type_map_string_to_string_to_list_Called_diploid_wrapper_ptr &map_genome_to_valname_to_calls__correctly_matched,
		const type_map_string_to_string_to_list_Called_diploid_wrapper_ptr &map_genome_to_valname_to_calls__alternatively_matched,
		const type_map_string_to_uint_to_Called_diploid_wrapper  &novel_mycalls);
		
	
void save_histogram_positive_negative_log_odds
			(const std::string &outdir,
			 const std::string &addl_id,
			const int &pvalue_for_RD_test_choice,
			 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);
	
	
	
	
	
	
	
//temkp
void write_collapse_file
	    (const std::string &outdir,
	     const std::string &nahr_or_geneConv,
	     const type_map_uint_to_haploid_outcome_to_uint &map_Event_to_number_of_occuring_individuals);
	    
	    
void collapse_Event_calls
	    (const std::string &outdir,
	     const std::map<uint, Conn_Comp> &Conn_Comps,
	     const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);
	
	    
	
void  remove_calls_on_sex_chromosomes
		(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls,
		const std::map<uint, Conn_Comp> &Conn_Comps);		
		
		
void read_Wilcoxon_pvalues_from_file
			(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls,
			const std::string &postdatadir);	
		
void read_matched_Wilcoxon_pvalues_from_file
			(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls,
			const std::string &postdatadir);	  
			
	
void save_raw_count_ratio_table
			(const std::string &outdir,
			 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
			const std::map<uint, Conn_Comp> &Conn_Comps,
			const std::string &gender_and_population_filename);
	
void debug_low_pvalues
			(const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);			
			
void identify_wilcoxon_bogus
		(type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);			
			
		
void read_all_raw_count_ratios_from_file
			(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls,
			const std::string &postdatadir);
		

void identify_wilcoxon_bad_calls
		(type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		 const std::map<uint, Conn_Comp> &Conn_Comps);
	
void map_genes_to_events_and_save_to_file
		(const type_map_uint_to_BI_to_Gene &Genes_in_the_human_genome,
		 const std::map<uint, Conn_Comp> &Conn_Comps,
		const std::string &filename);
		
type_map_uint_to_BI_to_Gene  read_cancer_sites_from_file
		    (const std::string &gene_fname);	
		    
type_map_string_to_uint_to_Called_diploid_wrapper load_restricted_calls_from_matlab
		(const std::string &restricted_calls_fname);		    
		    
		    
type_map_uint_to_BI_to_Gene load_psi_DR_pseudogene_database
				(const std::string &psidr_fname);		    
		    
void identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs
		(const type_map_uint_to_BI_to_Gene &Genes_in_the_human_genome,
		 const std::map<uint, Conn_Comp> &Conn_Comps);
	
void save_simple_theoretical_stats
		(const std::map<uint, Conn_Comp> &Conn_Comps);		
		
	
uint calculate_distance_to_nearest_MEPS
			(const Event &the_Event,
			const uint &some_breakpoint);

void remove_undefined_from_potential_space
			(std::map<uint, Conn_Comp> &Conn_Comps,
			const bool &remove_entire_cc_too);	
			
void remove_calls_whose_Event_or_Connected_Component_has_been_removed_or_are_sex_chromo
			    (const std::map<uint, Conn_Comp>  &Conn_Comps,
			     type_map_string_to_uint_to_Called_diploid_wrapper &mycalls);
	
void remove_non_dup_dels
		(std::map<uint, Conn_Comp> &Conn_Comps);			    
			
	
void upload_breakpoints_logodds
			(type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
			const std::string &brkpt_odds_fname);
		
	
			
type_map_uint_to_longdouble calculate_breakpoint_distribution_of_Event_for_MEPS
		(const Event &the_Event,
		const uint &MEPS_length,
		const uint &MEPS_num_mismatches,
		const int &probability_counting_method__0_uniform_by_vp____1_uniform_by_position);			
	
		
void calculate_mean_and_var
		(const type_map_uint_to_longdouble &some_dist,
		longdouble &mean, 
		longdouble &var);		
	
		
longdouble calculate_Lyapunov_CLT_statistic_for_MEPS
		(const std::map<uint, Conn_Comp> &Conn_Comps,
		 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		const longdouble &logfdr_threshold_for_brkpts,
		 const uint &MEPS_length,
		const uint &MEPS_num_mismatches,
		const int &probability_counting_method__0_uniform_by_vp____1_uniform_by_position,
		const type_map_uint_to_uint_to_longdouble &maps_ev_to_brkpt_dists);	
		
void calculate_various_MEPS_test_statistics_and_save
		(const std::map<uint, Conn_Comp> &Conn_Comps,
		 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		const longdouble &logfdr_threshold_for_brkpts);	
		
	
std::map<uint, Event>::const_iterator  find_Event
	(const std::map<uint, Conn_Comp> &Conn_Comps,
	const uint &the_ev);		
		
		
#endif