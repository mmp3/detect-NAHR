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

#include <limits.h>


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

#include <boost/math/distributions/normal.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>


#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>



#include <gmpfrxx.h>
#include <boost/math/bindings/mpfr.hpp>


#include <mpreal.h>
#include <mpfr.h>


#include <boost/random.hpp>


#include <api/BamReader.h>
#include <api/BamAux.h>
#include <api/BamAlignment.h>



#include <general_typedefs.h>
#include <translations_through_alignment.h>
#include <templates.h>

#include "Event.h"
#include "other_functions.h"
#include "globals.h"

#include "Conn_Comp.h"
#include "io_functions.h"
#include "Call_set.h"

#include "post.h"




uint Empirical_properties::MEPS_size = 300;
uint Empirical_properties::informativeness_fragment_length = 400;
uint Empirical_properties::informativeness_mate_size = 100;












void post_process
	    (type_map_uint_to_CC &Conn_Comps,
	     const std::string &postdatadir,      
	    const std::string &savedir,
	     const std::string &genes_filename,
	    const std::string &biological_validations_filename,
	const std::string &gender_and_population_filename)		    
{		
	std::cerr << "\n\n\n\tinside \"post_process\"...\n";
	
	const bool only_summary_true___compare_FDR_calls_false = false;
	
	
	remove_undefined_from_potential_space(Conn_Comps, true);//remove_entire_cc_too	
	remove_non_dup_dels(Conn_Comps);
	
	const type_map_uint_to_BI_to_Gene regular_genes(read_Genes_from_file("/users/mmparks/data/mmparks/biomart-no-pseudo.txt"));
	const type_map_uint_to_BI_to_Gene cancer_genes(read_cancer_sites_from_file("/users/mmparks/data/mmparks/biomart-cancer-sites.txt"));
	const type_map_uint_to_BI_to_Gene pseudo_genes(read_Genes_from_file("/users/mmparks/data/mmparks/biomart-pseudo.txt"));
	
// 	std::cerr << "map...\n";    
	map_genes_to_events_and_save_to_file(regular_genes, Conn_Comps, "regular_genes_mapped.txt");    
	map_genes_to_events_and_save_to_file(cancer_genes, Conn_Comps, "cancer_genes_mapped.txt");  
	map_genes_to_events_and_save_to_file(pseudo_genes, Conn_Comps, "pseudo_genes_mapped.txt");  
		        	
	save_simple_theoretical_stats(Conn_Comps);
	identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs(regular_genes, Conn_Comps);	
	
    
	std::cerr << "DONE map genes!!!\n\n";
    
    

    
    
	const real logodds_threshold("0.01");//-8);
	const int pvalue_for_RD_test_choice = 2; //0 = rd ratio, 1 = Wilcooxn test naive, 2 = Wilcoxon random matched control.
    
    
    //load and prepare:
    std::string allowable_genomes_filename(postdatadir);
    allowable_genomes_filename.append("/allowable_genomes");
    
    const type_set_string allowable_genomes(convert_list_to_set<std::string>(read_list_from_file<std::string>(allowable_genomes_filename)));
    
    type_map_string_to_uint_to_Called_diploid_wrapper mycalls;
    {//mycalls
		std::string mycalls_fname(postdatadir);
		mycalls_fname.append("/all_calls");
		
		mycalls = wrap_Sampled_diploid_Event_data(read_Sampled_Events_from_file(mycalls_fname)); 
		
// 		std::cerr << "\tload RD\n";
// 		load_RD_hypothesis_test_pvalues(postdatadir, mycalls);		
		
// 		std::cerr << "\tfill uniqueness\n";
// 		fill_uniqueness_and_undefinedness_of_each_call(Conn_Comps, mycalls);
		
// 		std::cerr << "\tread wilcoxon...\n";
// 		read_Wilcoxon_pvalues_from_file(mycalls, postdatadir);		
		
// 		std::cerr << "\tread wilcoxon matched...\n";
// 		read_matched_Wilcoxon_pvalues_from_file(mycalls, postdatadir);

		std::cerr << "\tread read_all_raw_count_ratios_from_file...\n";
		read_all_raw_count_ratios_from_file(mycalls, postdatadir);	
		
		upload_breakpoints_logodds(mycalls, "/users/mmparks/data/mmparks/breakpoint_ratios.txt");
		
		remove_calls_whose_Event_or_Connected_Component_has_been_removed_or_are_sex_chromo(Conn_Comps, mycalls);		
		
		
// 		std::cerr << "\tremove_calls_without_RD\n";
// 		remove_calls_without_RD(mycalls);
// 		remove_calls_in_undefined_regions_of_genome(mycalls);
// 		remove_calls_on_sex_chromosomes(mycalls, Conn_Comps);
    }//mycalls
    
// 	calculate_various_MEPS_test_statistics_and_save(Conn_Comps, mycalls, -6);
		
    
    std::cerr << "\tcompute theoretical space\n";
    //Theoretical:
    const Theoretical_compute_space the_theoretical_space(compile_Theoretical_compute_space(Conn_Comps, regular_genes)); 
    
    std::cerr << "\tsave theoretical\n";
    save_Theoretical_space_to_table(savedir, the_theoretical_space);    
    
    
// 	identify_wilcoxon_bad_calls(mycalls, Conn_Comps);
// //     debug_low_pvalues(mycalls);
// // 	identify_wilcoxon_bogus(mycalls);
	std::cerr << "\tsave_raw_count_ratio_table\n";
	save_raw_count_ratio_table(savedir, mycalls, Conn_Comps, gender_and_population_filename);
	
	if (only_summary_true___compare_FDR_calls_false)
	{  return;  } 


	//restrict calls according to fdr
	std::cerr << "\nload fdr from matlab...\n";
	const type_map_string_to_uint_to_Called_diploid_wrapper fdr_restricted_calls(load_restricted_calls_from_matlab("/users/mmparks/data/mmparks/job_outputs/results1/restricted_calls.txt"));   	
	{
		std::cerr << "\n\tfill logfdr and erase...\n";
		type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_gen = mycalls.begin();
		while (it_gen != mycalls.end())
		{
			std::cerr << "\t" << it_gen->first << "\n";
			const type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_rest_gen = fdr_restricted_calls.find(it_gen->first);
			if (it_rest_gen == fdr_restricted_calls.end())
			{
				mycalls.erase(it_gen++);
			}
			else
			{
				type_map_uint_to_Called_diploid_wrapper::iterator it_ev = it_gen->second.begin();
				while (it_ev != it_gen->second.end())
				{
					const type_map_uint_to_Called_diploid_wrapper::const_iterator it_rest_ev = it_rest_gen->second.find(it_ev->first);
					if (it_rest_ev == it_rest_gen->second.end())
					{
						it_gen->second.erase(it_ev++);
					}
					else
					{
						it_ev->second.logfdr_inferred = it_rest_ev->second.logfdr_inferred;
						++it_ev;
					}
				}//ev
				
				if (it_gen->second.empty())
				{
					mycalls.erase(it_gen++);
				}
				else
				{			
					++it_gen;
				}
			}				
		}//it_gen
	}//scope
		
	{//scope
		uint num_calls = 0;
		for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator itgen = mycalls.begin();
			itgen != mycalls.end();
			++itgen)
		{
			num_calls += itgen->second.size();
		}

		std::cerr << "mycalls.size() = " << mycalls.size() << ",  total call count = " << num_calls << "\n\n";
	}//scope
    
    
    
    std::cerr << "\tload validations\n";
    
    type_list_Validated_rearrangement validations_from_other_studies;
    {//validations_from_other_studies
		read_Validated_rearrangements_from_file__in_tabular_format(biological_validations_filename,  validations_from_other_studies); 
		restrict_to_allowable_genomes(validations_from_other_studies, allowable_genomes);   
		associate_Events_with_Validated_rearrangements(Conn_Comps, validations_from_other_studies);
    }//validations_from_other_studies
    
    
//     std::cerr << "\tload Genes\n";
    
//     const type_map_uint_to_BI_to_Gene regular_genes(read_Genes_from_file(genes_filename));            
    
//     std::cerr << "\tcompute theoretical space\n";
    //Theoretical:
//     const Theoretical_compute_space the_theoretical_space(compile_Theoretical_compute_space(Conn_Comps, regular_genes));
    
    
    
    
    //Summary:    
    type_map_string_to_Summary_Statistics summary_stats_per_individual;
    Summary_statistics summary_stats_across_ALL_individuals;
        
    std::cerr << "\tcompute summary stats\n";
    
    compile_Summary_statistics(
		    Conn_Comps,
		     mycalls,
		     mpfr::const_infinity(),
			pvalue_for_RD_test_choice,
		    summary_stats_per_individual,
		    summary_stats_across_ALL_individuals);    
    
    
    
    type_map_string_to_Summary_Statistics summary_stats_per_individual___thresh;
    Summary_statistics summary_stats_across_ALL_individuals___thresh;
        
    std::cerr << "\tcompute summary stats\n";
    
    compile_Summary_statistics(
		    Conn_Comps,
		     mycalls,
		    logodds_threshold,
			pvalue_for_RD_test_choice,
		    summary_stats_per_individual___thresh,
		    summary_stats_across_ALL_individuals___thresh);        
    
    
    
    //comparison
    
    type_vector_list_Validated_rearrangement_constiterator  Vals__correctly_matched;
    type_vector_list_Validated_rearrangement_constiterator  Vals__alternatively_matched;
    type_vector_list_Validated_rearrangement_constiterator  Vals__unmatched;
    type_map_string_to_string_to_list_Called_diploid_wrapper_ptr  map_genome_to_valname_to_calls__correctly_matched;
    type_map_string_to_string_to_list_Called_diploid_wrapper_ptr  map_genome_to_valname_to_calls__alternatively_matched;    
    type_map_string_to_uint_to_Called_diploid_wrapper  novel_mycalls;
           
    std::cerr << "\tcompute matches\n";        
    
    match_validations_with_Calls(
		    mycalls,
		     validations_from_other_studies,
		    Vals__correctly_matched,
		    Vals__alternatively_matched,
		    Vals__unmatched,
		    map_genome_to_valname_to_calls__correctly_matched,
		    map_genome_to_valname_to_calls__alternatively_matched,
		    novel_mycalls);
    
    
    
	
	
    //Genes:          
    type_map_string_to_Affected_Gene_catalog  affected_Genes_by_individual;
			
    std::cerr << "\tcompute affected genes\n";
    determine_Genes_affected_by_each_Event_outcome_on_each_individual(
		Conn_Comps,
		 regular_genes,
		 mpfr::const_infinity(),
		pvalue_for_RD_test_choice,
		mycalls,
		affected_Genes_by_individual);    
    
    
    type_map_string_to_Affected_Gene_catalog  affected_Genes_by_individual___thresh;
			
    std::cerr << "\tcompute affected genes\n";
    determine_Genes_affected_by_each_Event_outcome_on_each_individual(
		Conn_Comps,
		 regular_genes,
		 logodds_threshold,
		 pvalue_for_RD_test_choice,
		mycalls,
		affected_Genes_by_individual___thresh);        

		
		
//     std::cerr << "\tupload affected regions\n";
//     const type_map_uint_to_list_BI__string map_chromo_to_uploaded_regions(upload_affected_regions_sequences(Conn_Comps, mycalls));
//     
// 
    const type_map_uint_to_list_BI__string map_chromo_to_uploaded_regions;    
    
    Empirical_properties EP_for_NAHR;
    Empirical_properties EP_for_Gene_Conversion;
    Empirical_properties EP_for_both;
	
    std::cerr << "\tfill_Empirical_properties\n";
    
    fill_Empirical_properties(
			Conn_Comps,
			mycalls,
			map_chromo_to_uploaded_regions,
			EP_for_NAHR,
			EP_for_Gene_Conversion,
			false); 			    
    
    EP_for_both.absorb_other_Empirical_properties(EP_for_NAHR);
    EP_for_both.absorb_other_Empirical_properties(EP_for_Gene_Conversion);
            
    
    
    //save:    
    
    std::cerr << "\tsave theoretical\n";
    save_Theoretical_space_to_table(savedir, the_theoretical_space);
    
    std::cerr << "\tsave summaries\n";
    save_Summary_statistics_to_table(
			    savedir,
			    "full",
			    the_theoretical_space,
			    summary_stats_per_individual,
			    summary_stats_across_ALL_individuals,
			    affected_Genes_by_individual);
    
    save_Summary_statistics_to_table(
			    savedir,
			    "thresh",
			    the_theoretical_space,
			    summary_stats_per_individual___thresh,
			    summary_stats_across_ALL_individuals___thresh,
			    affected_Genes_by_individual___thresh);
        
    
    
//     std::cerr << "\tsave empiricals\n";
//     EP_for_NAHR.save_to_file(savedir, "NAHR");
//     EP_for_Gene_Conversion.save_to_file(savedir, "GeneConversion");
//     EP_for_both.save_to_file(savedir, "all");
    
    std::cerr << "\tsave raw positives\n";
    save_raw_positive_calls_to_files(savedir, "full", Conn_Comps, mycalls);
    
    const type_map_string_to_uint_to_Called_diploid_wrapper restricted_my_calls(restrict_to_NAHR_using_RD_threshold(mycalls, logodds_threshold, pvalue_for_RD_test_choice));
    save_raw_positive_calls_to_files(savedir, "thresh", Conn_Comps, restricted_my_calls);
	    
    
//     std::cerr << "\tsave FDR\n";
//     make_and_save_RD_FDR_distributions(
// 				savedir,
// 				collect_RD_hypothesis_test_p_values_for_FDR(mycalls));
    
    
    std::cerr << "\tsave Affected\n";
    save_Affected_Genes(
		savedir,
		affected_Genes_by_individual,
		regular_genes);  
    
	
    std::cerr << "\tsave comparisons\n";
    save_comparisons_to_tables(
		savedir,
		logodds_threshold,
		pvalue_for_RD_test_choice,
		Conn_Comps,
		mycalls,
		validations_from_other_studies,
		Vals__unmatched,
		map_genome_to_valname_to_calls__correctly_matched,
		map_genome_to_valname_to_calls__alternatively_matched,
		novel_mycalls);
	
    std::cerr << "\tsave positive-negative histogram\n";
    save_histogram_positive_negative_log_odds(savedir, "full", pvalue_for_RD_test_choice, mycalls);    
    save_histogram_positive_negative_log_odds(savedir, "thresh", pvalue_for_RD_test_choice, restricted_my_calls);    
	
// 	std::cerr << "\tsave_wilcoxon_histogram\n";
// 	save_wilcoxon_histogram(savedir, mycalls);
	
    
    //temp
//     collapse_Event_calls(savedir, Conn_Comps, mycalls);    
		
    
}//post_process





































void  match_validations_with_Calls
		    (type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		     const type_list_Validated_rearrangement &validations_from_other_studies,
		    type_vector_list_Validated_rearrangement_constiterator   &Vals__correctly_matched,
		    type_vector_list_Validated_rearrangement_constiterator   &Vals__alternatively_matched,
		    type_vector_list_Validated_rearrangement_constiterator   &Vals__unmatched,
		    type_map_string_to_string_to_list_Called_diploid_wrapper_ptr &map_genome_to_valname_to_calls__correctly_matched,
		    type_map_string_to_string_to_list_Called_diploid_wrapper_ptr &map_genome_to_valname_to_calls__alternatively_matched,
		    type_map_string_to_uint_to_Called_diploid_wrapper  &novel_mycalls)
{

    for (type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
	for (type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_mygen->second.begin();
		it_myev != it_mygen->second.end();
		++it_myev)
	{
	    bool compared_to_something = false;
	    
	    for (type_list_Validated_rearrangement::const_iterator it_val = validations_from_other_studies.begin();
		    it_val != validations_from_other_studies.end();
		    ++it_val)
	    {
		if (it_mygen->first.compare(it_val->NA_subject) == 0  and  it_val->related_UIDs.count(it_myev->first) > 0)
		{
		    compared_to_something = true;
		    
		    if((uint)it_myev->second.the_diploid_Event_outcome.first == it_val->outcome  or  (uint)it_myev->second.the_diploid_Event_outcome.second == it_val->outcome)
		    {
			it_myev->second.associated_other_calls.insert(it_val->variant_accession_name);
			Vals__correctly_matched.push_back(it_val);
			map_genome_to_valname_to_calls__correctly_matched[it_mygen->first][it_val->variant_accession_name].push_back(&it_myev->second);
		    }
		    else if ((test_if_haploid_outcome_is_GeneConversion(it_myev->second.the_diploid_Event_outcome.first) or test_if_haploid_outcome_is_GeneConversion(it_myev->second.the_diploid_Event_outcome.second))
				and (!test_if_haploid_outcome_is_NAHR(it_myev->second.the_diploid_Event_outcome.first)   and  !test_if_haploid_outcome_is_NAHR(it_myev->second.the_diploid_Event_outcome.second)))    
		    {
			it_myev->second.associated_other_calls.insert(it_val->variant_accession_name);
			Vals__alternatively_matched.push_back(it_val);
			map_genome_to_valname_to_calls__alternatively_matched[it_mygen->first][it_val->variant_accession_name].push_back(&it_myev->second);
		    }
		    else
		    {
			Vals__unmatched.push_back(it_val);
		    }		    		    
		}//same gene and event
	    }//it_val
	    
	    if (!compared_to_something  
		and  (test_if_haploid_outcome_is_NAHR(it_myev->second.the_diploid_Event_outcome.first)  or   test_if_haploid_outcome_is_NAHR(it_myev->second.the_diploid_Event_outcome.second))
		and  (!test_if_haploid_outcome_is_GeneConversion(it_myev->second.the_diploid_Event_outcome.first)  and  !test_if_haploid_outcome_is_GeneConversion(it_myev->second.the_diploid_Event_outcome.second))    )
	    {  novel_mycalls[it_mygen->first][it_myev->first] = it_myev->second;  }
	    
	}//it_myev	
    }//it_mygen    
    
    
    std::cerr << "\n\n\n\t"
		<< "validations_from_other_studies.size() = "<< validations_from_other_studies.size() << "\n\t"
		<< "Vals__correctly_matched.size() = "<< Vals__correctly_matched.size() << "\n\t"
		<< "Vals__alternatively_matched.size() = "<< Vals__alternatively_matched.size() << "\n\t"
		<< "Vals__unmatched.size() = "<< Vals__unmatched.size() << "\n\n\n\n";            
    
}//match_validations_with_Calls







void load_RD_hypothesis_test_pvalues
	    (const std::string &postdatadir,
	     type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
    std::string hypothfname(postdatadir);
    hypothfname.append("/all_RD_tests_by_event_and_genome");
    
    std::ifstream infs(hypothfname);
    
    
    std::string ingen;
    std::string inpop;
    uint incc;
    uint inev;    
    std::string inoddsalt;
    std::string inoddsdel;
    std::string inoddsdup;
    std::string inoddsbest;
    
    while (infs >> ingen)
    {
	infs >> inpop >> incc >> inev >> inoddsalt >> inoddsdel >> inoddsdup >> inoddsbest;
// 	std::cerr << "\tingen=[" << ingen << "]\n\tinev = " << inev << "\n\n";
	
	mycalls.at(ingen).at(inev).RD_odds_alt = real(inoddsalt);
	mycalls.at(ingen).at(inev).RD_odds_del = real(inoddsdel);
	mycalls.at(ingen).at(inev).RD_odds_dup = real(inoddsdup);
	mycalls.at(ingen).at(inev).RD_odds_best = real(inoddsbest);
	
	if (mycalls.at(ingen).at(inev).the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None))
		mycalls.at(ingen).at(inev).RD_odds_best = mpfr::log10(mpfr::max(mycalls.at(ingen).at(inev).RD_odds_del, mycalls.at(ingen).at(inev).RD_odds_dup));
// 		assert(mycalls.at(ingen).at(inev).RD_odds_best == mpfr::log10(mpfr::max(mycalls.at(ingen).at(inev).RD_odds_del, mycalls.at(ingen).at(inev).RD_odds_dup)));
	
	mycalls.at(ingen).at(inev).population_name = inpop;
	
	mycalls.at(ingen).at(inev).set_RD_odds = true;
    }
    
    
    infs.close();
    
    
}//load_RD_hypothesis_test_pvalues	    











type_map_string_to_uint_to_Called_diploid_wrapper wrap_Sampled_diploid_Event_data
							    (const type_map_string_to_uint_to_Sampled_diploid_Event_data &some_samples)
{
    
    type_map_string_to_uint_to_Called_diploid_wrapper wrapped_samples;
    
    for (type_map_string_to_uint_to_Sampled_diploid_Event_data::const_iterator it_gen = some_samples.begin();
	    it_gen != some_samples.end();
	    ++it_gen)
    {
	for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_samp = it_gen->second.begin();
		it_samp != it_gen->second.end();
		++it_samp)
	{
	    wrapped_samples[it_gen->first][it_samp->first] = it_samp->second;
	}	
    }    
    
    
    return  wrapped_samples;    
    
}//wrap_Sampled_diploid_Event_data







void associate_Events_with_Validated_rearrangements
			    (const type_map_uint_to_CC  &Conn_Comps,
			     type_list_Validated_rearrangement &validations_from_other_studies)
{
    
    for (type_list_Validated_rearrangement::iterator it_val = validations_from_other_studies.begin();
	    it_val != validations_from_other_studies.end();
	    ++it_val)
    {
	it_val->related_UIDs.clear();
	
	for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{
	    for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
		    it_ev != it_cc->second.events.end();
		    ++it_ev)
	    {
		if (it_val->chromosome == it_ev->second.chromos[0] and  it_val->recomb_type == it_ev->second.recomb_type)
		{   						
		    if (   (it_ev->second.recomb_type == recomb_class__DupDel  and  (it_val->recomb_type == 2 or  it_val->recomb_type  == 1))
			   or  (it_ev->second.recomb_type == recomb_class__Inv  and  it_val->recomb_type == 3) )
		    {
			const bool brkpts_are_suff_close_to_LCRs
					=  check_that_rearrangement_breakpoints_are_sufficiently_close_to_Event_LCRs
							(it_val->breakpoints,
							it_ev->second.LCRs,
							50);
							
			if (brkpts_are_suff_close_to_LCRs) 
			{
			    it_val->related_UIDs.insert(it_ev->second.UID);
			    break;
			}
		    }
		}
	    }//ev
	}//cc
    }    
    
}//associate_Events_with_Validated_rearrangements




































































bool check_that_rearrangement_breakpoints_are_sufficiently_close_to_Event_LCRs
	(const BOOST_Interval &rearrangement_breakpoints,
	 const BOOST_Interval *const &LCRs_of_some_Event,
	 const uint &leniancy)
{
    	    
    const int dist_to_LCR_0 
	    =   BOOST_in(  rearrangement_breakpoints.lower(),   LCRs_of_some_Event[0]  )    ?
		    0   :    std::min<int>(  abs(  (int)rearrangement_breakpoints.lower() - (int)LCRs_of_some_Event[0].lower() ),
					    abs(  (int)rearrangement_breakpoints.lower() - (int)LCRs_of_some_Event[0].upper() )  );                               
					    
    const int dist_to_LCR_1 
	    =   BOOST_in(  rearrangement_breakpoints.lower(),   LCRs_of_some_Event[1]  )    ?
		    0   :    std::min<int>(  abs(  (int)rearrangement_breakpoints.upper() - (int)LCRs_of_some_Event[1].lower() ),
					    abs(  (int)rearrangement_breakpoints.upper() - (int)LCRs_of_some_Event[1].upper() )  ); 
					    
    if (   (  BOOST_in(rearrangement_breakpoints.lower(), LCRs_of_some_Event[0])   or   dist_to_LCR_0 <= leniancy  )
	    and  
	   (  BOOST_in(rearrangement_breakpoints.upper(), LCRs_of_some_Event[1])   or   dist_to_LCR_1 <= leniancy   )  ) 
    {
	return true;
    }
    else
    {
	return false;
    }

    
} // check_that_rearrangement_breakpoints_are_sufficiently_close_to_Event_LCRs





























bool some_string_in_set_is_a_substring_of_given_string
		(const type_set_string &set_of_strings,
		 const std::string &given_string)
{

    for (type_set_string::const_iterator it_str = set_of_strings.begin();
	    it_str != set_of_strings.end();
	    ++it_str)
    {
	if ( given_string.find( *it_str ) != std::string::npos    )		
	{
	    return true;
	}    
    }
    
    return false;

}//some_string_in_set_is_a_substring_of_given_string




































                

                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                            
                
  


                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    


                    
                    
                    
                    
                    
void restrict_to_allowable_genomes
		(type_list_Validated_rearrangement  &calls_from_other_studies,
		const type_set_string &allowable_genomes_names)
{
         
    type_list_Validated_rearrangement::iterator it_val = calls_from_other_studies.begin();    
    while (it_val != calls_from_other_studies.end())
    {
	if (allowable_genomes_names.count(it_val->NA_subject) == 0)               
	    it_val = calls_from_other_studies.erase(it_val);	    
	else	    
	    ++it_val;	    	    
    }        
    
}







































void Summary_statistics::add_call_appropriately
			    (const type_map_uint_to_CC  &the_CCs,
			     const Sampled_diploid_Event_data &call)
{
    
    for (uint hap=0; hap<2; ++hap)
    {	
	++map_haploid_outcome_to_counts[pair_at<type_haploid_outcome>(call.the_diploid_Event_outcome,hap)];
	
	if (pair_at<type_haploid_outcome>(call.the_diploid_Event_outcome,hap) != hap_outcome__None)
	{
	    const type_BI__BI the_affected_regs(
			convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
							    pair_at<BOOST_Interval>(call.the_diploid_profile_brkpts,hap),
							    the_CCs.at(call.cc_of_event).events.at(call.event_UID).compressed_map_LCR_to_profile__for_each_LCR,
							    pair_at<type_haploid_outcome>(call.the_diploid_Event_outcome,hap)));
	    
	    if (test_if_haploid_outcome_is_NAHR(pair_at<type_haploid_outcome>(call.the_diploid_Event_outcome,hap)))
	    {//NAHR
		connected_components_with_some_NAHR_called_positive__haploid.insert(call.cc_of_event);
		total_DNA_affected_by_NAHR___diploid += BOOST_width_inclusive(the_affected_regs.first);
		++total_NAHR_calls__diploid;
	    }//NAHR
	    else
	    {//GeneConv
		connected_components_with_some_GeneConversion_called_positive__haploid.insert(call.cc_of_event);
		total_DNA_affected_by_GeneConversion___diploid += BOOST_width_inclusive(the_affected_regs.first);    
		++total_GeneConversion_calls__diploid;
	    }//GeneConv
	}//none	
    }//hap
    
    
    ++map_DIPLOID_outcome_to_counts[order_a_diploid_outcome(call.the_diploid_Event_outcome)];        

}//add_call_appropriately












				
				
// type_map_string_to_Gene identify_Genes_affected_by_Event_breakpoints
// 					(const Event &occurring_Event,
// 					const type_haploid_outcome &haploid_outcome,
// 					const BOOST_Interval &haploid_brkpts,
// 					const type_map_uint_to_BI_to_Gene &map_chromo_to_region_to_Gene)
// {
//     type_map_string_to_Gene  affected_Genes;    
//     
//     const type_BI__BI absolute_brkpts(
// 			    convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
// 					haploid_brkpts, occurring_Event.compressed_map_LCR_to_profile__for_each_LCR, haploid_outcome));
//     
//     const type_map_uint_to_BI_to_Gene::const_iterator it_chr = map_chromo_to_region_to_Gene.find(occurring_Event.chromos[0]);
//     
//     if (it_chr != map_chromo_to_region_to_Gene.end())
//     {
// 	for (type_map_BI_to_Gene::const_iterator it_bi = it_chr->second.begin();
// 		it_bi != it_chr->second.end();
// 		++it_bi)
// 	{
// 	    if (BOOST_overlap(absolute_brkpts.first, it_bi->first) or  BOOST_overlap(absolute_brkpts.second, it_bi->first))
// 	    {
// 		affected_Genes[it_bi->second.ensembl_id] = it_bi->second;
// 	    }
// 	}//bi	
//     }//chr
//         
//     return affected_Genes;    
//     
// }//identify_Genes_affected_by_Event_breakpoints











type_map_string_to_Gene identify_Genes_potentially_affected_by_an_Event
					(const Event &occurring_Event,
					const type_map_uint_to_BI_to_Gene &map_chromo_to_region_to_Gene)
{
    type_map_string_to_Gene  potentially_affected_Genes;    
    
    const type_map_uint_to_BI_to_Gene::const_iterator it_chr = map_chromo_to_region_to_Gene.find(occurring_Event.chromos[0]);    
    if (it_chr != map_chromo_to_region_to_Gene.end())
    {	
	for (type_map_BI_to_Gene::const_iterator it_bi = it_chr->second.begin();
		it_bi != it_chr->second.end();
		++it_bi)
	{
	    if (BOOST_overlap(occurring_Event.region_between_and_including_the_LCRs_themselves, it_bi->first))
	    {
		potentially_affected_Genes[it_bi->second.ensembl_id] = it_bi->second;
	    }
	}//bi	
    }//chr
        
    return potentially_affected_Genes;    
    
}//identify_Genes_potentially_affected_by_an_Event











void determine_Genes_affected_by_each_Event_outcome_on_each_individual
		(const type_map_uint_to_CC &Conn_Comps,
		 const type_map_uint_to_BI_to_Gene &all_Genes,
		const real &logodds_threshold,
		const int &pvalue_for_RD_test_choice,
		type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		type_map_string_to_Affected_Gene_catalog &affected_Genes_by_individual)
{
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
		affected_Genes_by_individual.insert(std::pair<std::string, Affected_Gene_catalog>(it_mygen->first, Affected_Gene_catalog()));
		
		for (type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{	    
			bool passes_threshold;
			if (pvalue_for_RD_test_choice == 0)
			{
				passes_threshold = mpfr::log10(it_myev->second.RD_odds_alt) <= logodds_threshold;				
			}
			else if (pvalue_for_RD_test_choice == 1)
			{
				passes_threshold = mpfr::log10(it_myev->second.Wilcoxon_signed_rank_sum_test__pvalue__inferred) <= logodds_threshold;
			}
			else if (pvalue_for_RD_test_choice == 2)
			{
				passes_threshold = mpfr::log10(it_myev->second.wilcoxon_matched_info___null___inferred.wilcoxon_pvalue) <= logodds_threshold;
			}
			
			if (passes_threshold)
			{	    
			const type_map_uint_to_Event::const_iterator the_Event_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID);	        
			
			const type_map_uint_to_BI_to_Gene::const_iterator it_genes_chr = all_Genes.find(the_Event_it->second.chromos[0]);
					
			for (uint hap=0; hap<2; ++hap)
			{//hap mycall		    
				switch(pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap))
				{
				case hap_outcome__Del:
					affected_Genes_by_individual.at(it_mygen->first).Events_occurring_as_Deletions.insert(it_myev->second.event_UID);
					break;
				case hap_outcome__Dup:
					affected_Genes_by_individual.at(it_mygen->first).Events_occurring_as_Duplications.insert(it_myev->second.event_UID);
					break;
				case hap_outcome__Inv:
					affected_Genes_by_individual.at(it_mygen->first).Events_occurring_as_Inversions.insert(it_myev->second.event_UID);
					break;
				case hap_outcome__GeneConv_ABA:
				case hap_outcome__GeneConv_BAB:
					affected_Genes_by_individual.at(it_mygen->first).Events_occurring_as_GeneConversions.insert(it_myev->second.event_UID);
					break;
				default:
					break;
				};		
			}//hap
			
			
			
			if (it_genes_chr != all_Genes.end())
			{    
				for (uint hap=0; hap<2; ++hap)
				{//hap mycall
				if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) == hap_outcome__None)
				{  continue;  }
				
				const type_BI__BI the_hap_absolute_brkpts(
						convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
							pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
							the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
							pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap)));

				for (type_map_BI_to_Gene::const_iterator it_bi_gene = it_genes_chr->second.begin();
					it_bi_gene != it_genes_chr->second.end();
					++it_bi_gene)
				{
					if (BOOST_overlap(the_hap_absolute_brkpts.first, it_bi_gene->second.coordinates)
						or  BOOST_overlap(the_hap_absolute_brkpts.second, it_bi_gene->second.coordinates))
					{
					evaluate_impact_of_Call_on_Gene_appropriately(
						the_Event_it->second,
						it_myev->second,
						hap,
						it_bi_gene->second,
						affected_Genes_by_individual);
					}
				}//it_bi_gene
				}//hap mycall
			}//found chr	
			}//thresh
		}//it_myev	
    }//it_mygen   
    
}//determine_Genes_affected_by_each_Event_outcome_on_each_individual










void evaluate_impact_of_Call_on_Gene_appropriately
		    (const Event &the_Event,
		     Called_diploid_wrapper &the_call,
		     const uint &hap,
		    const Gene &the_affected_Gene,
		    type_map_string_to_Affected_Gene_catalog  &affected_Genes_per_individual)
{
    //assume the Event is occurring and really does affect the gene
    if (pair_at<type_haploid_outcome>(the_call.the_diploid_Event_outcome,hap) != hap_outcome__None)
    {    
	pair_at<type_set_string>(the_call.affected_genes_per_haploid,hap).insert(the_affected_Gene.ensembl_id);        
	
	switch(pair_at<type_haploid_outcome>(the_call.the_diploid_Event_outcome,hap))
	{
	    case hap_outcome__Del:
		affected_Genes_per_individual[the_call.associated_genome].Genes_affected_by_Deletions.insert(std::pair<std::string, Gene>(the_affected_Gene.ensembl_id, the_affected_Gene));
		break;
	    case hap_outcome__Dup:
		affected_Genes_per_individual[the_call.associated_genome].Genes_affected_by_Duplications.insert(std::pair<std::string, Gene>(the_affected_Gene.ensembl_id, the_affected_Gene));
		break;
	    case hap_outcome__Inv:
		affected_Genes_per_individual[the_call.associated_genome].Genes_affected_by_Inversions.insert(std::pair<std::string, Gene>(the_affected_Gene.ensembl_id, the_affected_Gene));
		break;
	    case hap_outcome__GeneConv_ABA:
		if (BOOST_overlap(the_Event.LCRs[0], the_affected_Gene.coordinates))
		{
		    affected_Genes_per_individual[the_call.associated_genome].Genes_affected_by_GeneConversion__received_from_somewhere_else.insert(std::pair<std::string, Gene>(the_affected_Gene.ensembl_id, the_affected_Gene));
		}
		
		if (BOOST_overlap(the_Event.LCRs[1], the_affected_Gene.coordinates))			    
		{
		    affected_Genes_per_individual[the_call.associated_genome].Genes_affected_by_GeneConversion__donated_itself_to_elswhere.insert(std::pair<std::string, Gene>(the_affected_Gene.ensembl_id, the_affected_Gene));	    
		}    
		
		break;
	    case hap_outcome__GeneConv_BAB:
		if (BOOST_overlap(the_Event.LCRs[0], the_affected_Gene.coordinates))
		{
		    affected_Genes_per_individual[the_call.associated_genome].Genes_affected_by_GeneConversion__donated_itself_to_elswhere.insert(std::pair<std::string, Gene>(the_affected_Gene.ensembl_id, the_affected_Gene));
		}
		
		if (BOOST_overlap(the_Event.LCRs[1], the_affected_Gene.coordinates))			    
		{
		    affected_Genes_per_individual[the_call.associated_genome].Genes_affected_by_GeneConversion__received_from_somewhere_else.insert(std::pair<std::string, Gene>(the_affected_Gene.ensembl_id, the_affected_Gene));	    
		}    
		
		break;
		
	    default:
		std::stringstream errostrm;
		errostrm << "ERROR! in \"evaluate_impact_of_Call_on_Gene_appropriately\"  !!! invalid haploid outcome = "
			    << convert_haploid_outcome_to_string(pair_at<type_haploid_outcome>(the_call.the_diploid_Event_outcome,hap)) << "\n\n";
		error_message(errostrm, false);
		break;			    
	};//outcomes
    }//occurs
    
    
}//evaluate_impact_of_Call_on_Gene_appropriately




bool Gene::operator< (const Gene &another_Gene)  const
{
    return  (ensembl_id.compare(another_Gene.ensembl_id) < 0);
}






type_map_string_to_Gene Affected_Gene_catalog::get_all_affected_Genes() const
{
    
    type_map_string_to_Gene all_affected_Genes(Genes_affected_by_Inversions);
    all_affected_Genes.insert(Genes_affected_by_Deletions.begin(), Genes_affected_by_Deletions.end());
    all_affected_Genes.insert(Genes_affected_by_Duplications.begin(), Genes_affected_by_Duplications.end());
    all_affected_Genes.insert(Genes_affected_by_GeneConversion__donated_itself_to_elswhere.begin(), Genes_affected_by_GeneConversion__donated_itself_to_elswhere.end());
    all_affected_Genes.insert(Genes_affected_by_GeneConversion__received_from_somewhere_else.begin(), Genes_affected_by_GeneConversion__received_from_somewhere_else.end());
    
    
    return all_affected_Genes;

}//get_all_affected_Genes









void Affected_Gene_catalog::collect_from_individual_catalogs
				    (const type_map_string_to_Affected_Gene_catalog& catalogs_per_individual)
{
    for (type_map_string_to_Affected_Gene_catalog::const_iterator it_ind = catalogs_per_individual.begin();
	    it_ind != catalogs_per_individual.end();
	    ++it_ind)
    {
	Genes_affected_by_Inversions.insert(it_ind->second.Genes_affected_by_Inversions.begin(), it_ind->second.Genes_affected_by_Inversions.end());
	Genes_affected_by_Deletions.insert(it_ind->second.Genes_affected_by_Deletions.begin(), it_ind->second.Genes_affected_by_Deletions.end());
	Genes_affected_by_Duplications.insert(it_ind->second.Genes_affected_by_Duplications.begin(), it_ind->second.Genes_affected_by_Duplications.end());
	
	Genes_affected_by_GeneConversion__donated_itself_to_elswhere.insert(it_ind->second.Genes_affected_by_GeneConversion__donated_itself_to_elswhere.begin(), it_ind->second.Genes_affected_by_GeneConversion__donated_itself_to_elswhere.end());
	Genes_affected_by_GeneConversion__received_from_somewhere_else.insert(it_ind->second.Genes_affected_by_GeneConversion__received_from_somewhere_else.begin(), it_ind->second.Genes_affected_by_GeneConversion__received_from_somewhere_else.end());	
	
	Events_occurring_as_Deletions.insert(it_ind->second.Events_occurring_as_Deletions.begin(), it_ind->second.Events_occurring_as_Deletions.end());
	Events_occurring_as_Duplications.insert(it_ind->second.Events_occurring_as_Duplications.begin(), it_ind->second.Events_occurring_as_Duplications.end());
	Events_occurring_as_Inversions.insert(it_ind->second.Events_occurring_as_Inversions.begin(), it_ind->second.Events_occurring_as_Inversions.end());
	Events_occurring_as_GeneConversions.insert(it_ind->second.Events_occurring_as_GeneConversions.begin(), it_ind->second.Events_occurring_as_GeneConversions.end());
    }    
    

}//collect_from_individual_catalogs








longdouble determine_best_MEPS_homology_of_profile_breakpoint
		    (const Event &the_Event,
		    const uint &the_profile_breakpoint,
		    const uint &MEPS_size_to_consider)
		    		 
{
    const uint prev_vp =  (the_Event.uploaded_variational_positions_of_profile.find(the_profile_breakpoint) == the_Event.uploaded_variational_positions_of_profile.begin())  ?  0  :   *--the_Event.uploaded_variational_positions_of_profile.find(the_profile_breakpoint);
    
    const uint first_possible_brkpt = (prev_vp == 0)  ?   0  :  prev_vp+1;
    
    longdouble best_MEPs_homol = 0;
    
    for (uint possible_brkpt = first_possible_brkpt;  possible_brkpt <= the_profile_breakpoint;  ++possible_brkpt)
    {
	const BOOST_Interval MEPS_region(possible_brkpt, possible_brkpt + MEPS_size_to_consider - 1);
	const uint additional_non_repeat_amount = (MEPS_region.upper() >= the_Event.profile_length)  ?  MEPS_region.upper() - the_Event.profile_length + 1  :  0;
	
	const longdouble test_MEPS_homol = (longdouble)(get_intersection_of_set_and_interval(the_Event.uploaded_variational_positions_of_profile, MEPS_region).size() + additional_non_repeat_amount) / MEPS_size_to_consider;
	
	if (test_MEPS_homol > best_MEPs_homol)
	{  best_MEPs_homol = test_MEPS_homol;  }
    }//possible_brkpt
    
        
    return  best_MEPs_homol;
    
}//determine_best_MEPS_homology_of_profile_breakpoint





uint calculate_GC_content_of_region
		(const uint &chromosome,
		const BOOST_Interval &region,
		const type_map_uint_to_list_BI__string &map_chromo_to_uploaded_regions)
{
    const type_map_uint_to_list_BI__string::const_iterator it_find_chr = map_chromo_to_uploaded_regions.find(chromosome);    
    if (it_find_chr != map_chromo_to_uploaded_regions.end())
    {
	for (type_list__BI__string::const_iterator it_bi_string = it_find_chr->second.begin();
		it_bi_string != it_find_chr->second.end();
		++it_bi_string)
	{
	    if (BOOST_subset(region, it_bi_string->first))
	    {
		const uint relative_start_pos = region.lower() - it_bi_string->first.lower();
		const std::string desired_region_seq(it_bi_string->second.substr(relative_start_pos, BOOST_width_inclusive(region)));
		
		return  calculate_GC_content_of_sequence(desired_region_seq);		
	    }//subset
	}//bi		
    }//chr
    
    
    //If made it this far - ERROR!
    {//ERROR
	std::stringstream errostrm;
	errostrm << "ERROR!   in  \"calculate_GC_content_of_region\".\n"
		<< "(it_find_chr == map_chromo_to_uploaded_regions.end())  =  "  << (it_find_chr == map_chromo_to_uploaded_regions.end()) << "\n"
		<< "chromosome = " << chromosome << "\n"
		<< "region = " << region << "\n\n";
	error_message(errostrm,true);
    }//ERROR
    
    return 0; //dummy
        
}//calculate_GC_content_of_region







void write_collapse_file
	    (const std::string &outdir,
	     const std::string &nahr_or_geneConv,
	     const type_map_uint_to_haploid_outcome_to_uint &map_Event_to_number_of_occuring_individuals)
{
    std::string outfname(outdir);
    outfname.append("/collapsed_");
    outfname.append(nahr_or_geneConv);
    outfname.append(".txt");
    
    std::ofstream outfs(outfname);
    check_fstream_for_goodness(outfs, outfname, "write_collapse_file", true);
    
    outfs
	<< "Event UID"  << "\t"
	<< "recomb type" << "\t"
	<< "number of individuals" << "\n";
	
	
    for (type_map_uint_to_haploid_outcome_to_uint::const_iterator it_ev = map_Event_to_number_of_occuring_individuals.begin();
	    it_ev != map_Event_to_number_of_occuring_individuals.end();
	    ++it_ev)
    {
	for (type_map_haploid_outcome_to_uint::const_iterator it_out = it_ev->second.begin();
		it_out != it_ev->second.end();
		++it_out)
	{
	    outfs
		<< it_ev->first << "\t"
		<< convert_haploid_outcome_to_string(it_out->first) << "\t"
		<< it_out->second << "\n";	    	    
	}//out	
    }//ev
    
    outfs.close();    
    
    
}//write_collapse_file


void collapse_Event_calls
	    (const std::string &outdir,
	     const type_map_uint_to_CC &Conn_Comps,
	     const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
   
    type_map_uint_to_haploid_outcome_to_uint  map_Event_to_number_of_occuring_individuals___NAHR;
    type_map_uint_to_haploid_outcome_to_uint  map_Event_to_number_of_occuring_individuals___GeneConversion;
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
	for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
		it_myev != it_mygen->second.end();
		++it_myev)
	{	    
	    const type_map_uint_to_Event::const_iterator the_Event_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID);
	    
	    const uint hapmax =  (it_myev->second.the_diploid_Event_outcome.first == it_myev->second.the_diploid_Event_outcome.second)   ?   1  :  2;
	    
	    for (uint hap=0; hap < hapmax; ++hap)
	    {
		const type_haploid_outcome the_hap_outcome = pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap);
		
		if (the_hap_outcome != hap_outcome__None)
		{		    
		    const type_BI__BI affected_regions(
				convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
					pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
					the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
					the_hap_outcome));

		    const type_map_uint_to_uint map_vp_to_sparseness(
							invert_multimap_into_map<uint,uint>(
								get_sparseness_of_variational_positions(the_Event_it->second.uploaded_variational_positions_of_profile)));				    
		    
		    
		    if (test_if_haploid_outcome_is_NAHR(the_hap_outcome))
		    {
			++map_Event_to_number_of_occuring_individuals___NAHR[it_myev->second.event_UID][the_hap_outcome];			
		    }//NAHR
		
		    if (test_if_haploid_outcome_is_GeneConversion(the_hap_outcome))
		    {		
			++map_Event_to_number_of_occuring_individuals___GeneConversion[it_myev->second.event_UID][the_hap_outcome];
		    }//GeneConversion				
		}//non-null outcome
	    }//hap	    	    
	}//it_myev	
    }//it_mygen   
    
    
    
    write_collapse_file(outdir, "NAHR",  map_Event_to_number_of_occuring_individuals___NAHR);
    write_collapse_file(outdir, "GeneConversion",  map_Event_to_number_of_occuring_individuals___GeneConversion);
        
    
}//collapse_Event_calls

















void fill_Empirical_properties
	    (const type_map_uint_to_CC &Conn_Comps,
	     const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
	    const type_map_uint_to_list_BI__string &map_chromo_to_uploaded_regions,
	    Empirical_properties &EP_for_NAHR,
	    Empirical_properties &EP_for_Gene_Conversion,
	    const bool &fill_compute_intensive_properties_too)    
{
    const uint threshold_to_be_LG_for_LCR_identity = 100;
    
    type_map_uint_to_haploid_outcome_to_uint  map_Event_to_number_of_occuring_individuals___NAHR;
    type_map_uint_to_haploid_outcome_to_uint  map_Event_to_number_of_occuring_individuals___GeneConversion;
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
	for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
		it_myev != it_mygen->second.end();
		++it_myev)
	{	    
	    const type_map_uint_to_Event::const_iterator the_Event_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID);
	    
	    for (uint hap=0; hap<2; ++hap)
	    {
		const type_haploid_outcome the_hap_outcome = pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap);
		
		if (the_hap_outcome != hap_outcome__None)
		{		    
		    const type_BI__BI affected_regions(
				convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
					pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
					the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
					the_hap_outcome));

		    const type_map_uint_to_uint map_vp_to_sparseness(
							invert_multimap_into_map<uint,uint>(
								get_sparseness_of_variational_positions(the_Event_it->second.uploaded_variational_positions_of_profile)));				    
		    
		    
		    if (test_if_haploid_outcome_is_NAHR(the_hap_outcome))
		    {							
			{//sparseness		    		        		    		    
			    const type_map_uint_to_uint::const_iterator it_vp_to_sparseness
					= map_vp_to_sparseness.find(pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).lower());
					
			    if (it_vp_to_sparseness != map_vp_to_sparseness.end())
			    {  EP_for_NAHR.sparsity_histogram_of_breakpoints.push_back(it_vp_to_sparseness->second);  }			    			    
			}//sparseness				
			
			
			EP_for_NAHR.mediating_LCR_overall_percent_identity_histogram.push_back(
					    determine_LCR_percent_identity__ignoring_large_gaps(the_Event_it->second, threshold_to_be_LG_for_LCR_identity));
			
			EP_for_NAHR.mediating_LCR_length_histogram.push_back((the_Event_it->second.LCR_lengths[0] + the_Event_it->second.LCR_lengths[1])/2);			
			
			EP_for_NAHR.distance_between_mediating_LCRs_histogram.push_back(BOOST_width_inclusive(the_Event_it->second.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs));
			
			
			EP_for_NAHR.affected_region_size_histogram.push_back(BOOST_width_inclusive(affected_regions.first));
			
			if (it_myev->second.the_diploid_Event_outcome.first != it_myev->second.the_diploid_Event_outcome.second   or   hap == 0)
			{  ++map_Event_to_number_of_occuring_individuals___NAHR[it_myev->second.event_UID][the_hap_outcome];  }
			
			if (fill_compute_intensive_properties_too)
			{
			    EP_for_NAHR.best_MEPS_homology_histogram.push_back(determine_best_MEPS_homology_of_profile_breakpoint(
							    the_Event_it->second,
							    pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).lower(),
							    Empirical_properties::MEPS_size));
			}
			
			EP_for_NAHR.posterior_probabilities_of_MAP_event_outcome__histogram.push_back(it_myev->second.P_diploid_outcome);
			EP_for_NAHR.posterior_probabilities_of_MAP_event_breakpoints__histogram.push_back(it_myev->second.P_diploid_breakpoint);
						
		    }//NAHR
		
		    if (test_if_haploid_outcome_is_GeneConversion(the_hap_outcome))
		    {			
			{//sparseness		    		        		    		    
			    const type_map_uint_to_uint::const_iterator it_vp_to_sparseness
					= map_vp_to_sparseness.find(pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).lower());
					
			    if (it_vp_to_sparseness != map_vp_to_sparseness.end())
			    {  EP_for_Gene_Conversion.sparsity_histogram_of_breakpoints.push_back(it_vp_to_sparseness->second);  }
			}//sparseness	
			
			{//sparseness		    		        		    		    
			    const type_map_uint_to_uint::const_iterator it_vp_to_sparseness
					= map_vp_to_sparseness.find(pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).upper());
					
			    if (it_vp_to_sparseness != map_vp_to_sparseness.end())
			    {  EP_for_Gene_Conversion.sparsity_histogram_of_breakpoints.push_back(it_vp_to_sparseness->second);  }
			}//sparseness			
			
			
			EP_for_Gene_Conversion.mediating_LCR_overall_percent_identity_histogram.push_back(
					    determine_LCR_percent_identity__ignoring_large_gaps(the_Event_it->second, threshold_to_be_LG_for_LCR_identity));
			
			EP_for_Gene_Conversion.mediating_LCR_length_histogram.push_back((the_Event_it->second.LCR_lengths[0] + the_Event_it->second.LCR_lengths[1])/2);			
			
			EP_for_Gene_Conversion.distance_between_mediating_LCRs_histogram.push_back(BOOST_width_inclusive(the_Event_it->second.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs));
			
						
			if (it_myev->second.the_diploid_Event_outcome.first != it_myev->second.the_diploid_Event_outcome.second   or   hap == 0)
			{  ++map_Event_to_number_of_occuring_individuals___GeneConversion[it_myev->second.event_UID][the_hap_outcome];  }
			
			if (fill_compute_intensive_properties_too)
			{
			    EP_for_Gene_Conversion.best_MEPS_homology_histogram.push_back(determine_best_MEPS_homology_of_profile_breakpoint(
							    the_Event_it->second,
							    pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).lower(),
							    Empirical_properties::MEPS_size));
			
			    
			    const BOOST_Interval donor_tract = (the_hap_outcome == hap_outcome__GeneConv_ABA)  ?  affected_regions.second  :  affected_regions.first;
			    
			    const BOOST_Interval receiving_tract = (the_hap_outcome == hap_outcome__GeneConv_ABA)  ?  affected_regions.first  :  affected_regions.second;
			    
			    EP_for_Gene_Conversion.length_of_donor_GeneConversion_tract__histogram.push_back(BOOST_width_inclusive(donor_tract));
			    EP_for_Gene_Conversion.length_of_receiving_GeneConversion_tract__histogram.push_back(BOOST_width_inclusive(receiving_tract));
			
			
			    const uint GC_content_of_donor_tract 
					= calculate_GC_content_of_region(
								the_Event_it->second.chromos[0], 
								donor_tract,
								map_chromo_to_uploaded_regions);
			
			    const uint GC_content_of_receiving_tract 
					= calculate_GC_content_of_region(
								the_Event_it->second.chromos[0], 
								receiving_tract,
								map_chromo_to_uploaded_regions);

					
			    EP_for_Gene_Conversion.ratioGC_content_of_GeneConversion_tract__histogram.push_back(
				    ((longdouble)GC_content_of_donor_tract/BOOST_width_inclusive(donor_tract))
				    /((longdouble)GC_content_of_receiving_tract/BOOST_width_inclusive(receiving_tract)));
			}
			
			EP_for_Gene_Conversion.posterior_probabilities_of_MAP_event_outcome__histogram.push_back(it_myev->second.P_diploid_outcome);
			EP_for_Gene_Conversion.posterior_probabilities_of_MAP_event_breakpoints__histogram.push_back(it_myev->second.P_diploid_breakpoint);
			
		    }//GeneConversion				
		}//non-null outcome
	    }//hap	    	    
	}//it_myev	
    }//it_mygen   
    
    
    
    
    for (type_map_uint_to_haploid_outcome_to_uint::const_iterator it_ev_to_hapout = map_Event_to_number_of_occuring_individuals___NAHR.begin();
	    it_ev_to_hapout != map_Event_to_number_of_occuring_individuals___NAHR.end();
	    ++it_ev_to_hapout)
    {
	for (type_map_haploid_outcome_to_uint::const_iterator it_hapout = it_ev_to_hapout->second.begin();
		it_hapout != it_ev_to_hapout->second.end();
		++it_hapout)
	{  EP_for_NAHR.number_of_occurring_individuals__vs__number_of_occurring_Events__histogram.push_back(it_hapout->second);	 }
    }//it_ev_to_hapout       
    
    
    
    for (type_map_uint_to_haploid_outcome_to_uint::const_iterator it_ev_to_hapout = map_Event_to_number_of_occuring_individuals___GeneConversion.begin();
	    it_ev_to_hapout != map_Event_to_number_of_occuring_individuals___GeneConversion.end();
	    ++it_ev_to_hapout)
    {
	for (type_map_haploid_outcome_to_uint::const_iterator it_hapout = it_ev_to_hapout->second.begin();
		it_hapout != it_ev_to_hapout->second.end();
		++it_hapout)
	{  EP_for_Gene_Conversion.number_of_occurring_individuals__vs__number_of_occurring_Events__histogram.push_back(it_hapout->second);	 }
    }//it_ev_to_hapout         
    
    
}//make_histogram_length_of_GeneConversion_tracts










longdouble determine_LCR_percent_identity__ignoring_large_gaps
			    (const Event &the_Event,
			    const uint &min_gap_size_to_be_considered_large_gap_for_identity_calculation)
{
    uint amount_of_profile_that_is_above_LG_threshold = 0;
    uint amount_of_profile_that_is_LG_but_below_threshold = 0;
    
    for (type_set_BI::const_iterator it_lg = the_Event.large_gaps_of_profile.begin();
	    it_lg != the_Event.large_gaps_of_profile.end();
	    ++it_lg)
    {
	if (BOOST_width_inclusive(*it_lg) >= min_gap_size_to_be_considered_large_gap_for_identity_calculation)
	{  amount_of_profile_that_is_above_LG_threshold += BOOST_width_inclusive(*it_lg);  }
	else
	{  amount_of_profile_that_is_LG_but_below_threshold  +=  BOOST_width_inclusive(*it_lg);  }		
    }//lg
    
    
    const uint denom = the_Event.profile_length - amount_of_profile_that_is_above_LG_threshold;
    const uint numer = the_Event.uploaded_variational_positions_of_profile.size() + amount_of_profile_that_is_LG_but_below_threshold;
    
    return  (1.00L - (longdouble)numer/denom);
    
    
}//determine_LCR_percent_identity__ignoring_large_gaps










void Empirical_properties::absorb_other_Empirical_properties(const Empirical_properties& other_EP)
{
    sparsity_histogram_of_breakpoints.insert(sparsity_histogram_of_breakpoints.end(), 
					     other_EP.sparsity_histogram_of_breakpoints.begin(), other_EP.sparsity_histogram_of_breakpoints.end());
    
    mediating_LCR_length_histogram.insert(mediating_LCR_length_histogram.end(),
					other_EP.mediating_LCR_length_histogram.begin(), other_EP.mediating_LCR_length_histogram.end());
    
    mediating_LCR_overall_percent_identity_histogram.insert(mediating_LCR_overall_percent_identity_histogram.end(),
							    other_EP.mediating_LCR_overall_percent_identity_histogram.begin(), other_EP.mediating_LCR_overall_percent_identity_histogram.end());
    
    distance_between_mediating_LCRs_histogram.insert(distance_between_mediating_LCRs_histogram.end(),
						     other_EP.distance_between_mediating_LCRs_histogram.begin(), other_EP.distance_between_mediating_LCRs_histogram.end());
    
    best_MEPS_homology_histogram.insert(best_MEPS_homology_histogram.end(), 
					other_EP.best_MEPS_homology_histogram.begin(), other_EP.best_MEPS_homology_histogram.end());
    
    affected_region_size_histogram.insert(affected_region_size_histogram.end(), 
					  other_EP.affected_region_size_histogram.begin(), other_EP.affected_region_size_histogram.end());
    
    length_of_donor_GeneConversion_tract__histogram.insert(length_of_donor_GeneConversion_tract__histogram.end(),
							   other_EP.length_of_donor_GeneConversion_tract__histogram.begin(), other_EP.length_of_donor_GeneConversion_tract__histogram.end());
    
    length_of_receiving_GeneConversion_tract__histogram.insert(length_of_receiving_GeneConversion_tract__histogram.end(),
							       other_EP.length_of_receiving_GeneConversion_tract__histogram.begin(), other_EP.length_of_receiving_GeneConversion_tract__histogram.end());
    
    ratioGC_content_of_GeneConversion_tract__histogram.insert(ratioGC_content_of_GeneConversion_tract__histogram.end(),
							      other_EP.ratioGC_content_of_GeneConversion_tract__histogram.begin(), other_EP.ratioGC_content_of_GeneConversion_tract__histogram.end());

    number_of_occurring_individuals__vs__number_of_occurring_Events__histogram
				.insert(number_of_occurring_individuals__vs__number_of_occurring_Events__histogram.end(),
					other_EP.number_of_occurring_individuals__vs__number_of_occurring_Events__histogram.begin(), other_EP.number_of_occurring_individuals__vs__number_of_occurring_Events__histogram.end());
				
    posterior_probabilities_of_MAP_event_outcome__histogram
				.insert(posterior_probabilities_of_MAP_event_outcome__histogram.end(),							  
					other_EP.posterior_probabilities_of_MAP_event_outcome__histogram.begin(), other_EP.posterior_probabilities_of_MAP_event_outcome__histogram.end());				

    posterior_probabilities_of_MAP_event_breakpoints__histogram
				.insert(posterior_probabilities_of_MAP_event_breakpoints__histogram.end(),							  
					other_EP.posterior_probabilities_of_MAP_event_breakpoints__histogram.begin(), other_EP.posterior_probabilities_of_MAP_event_breakpoints__histogram.end());				
				
				
}//absorb_other_Empirical_properties













void compile_Summary_statistics
		    (const type_map_uint_to_CC &Conn_Comps,
		     const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		    const real &logodds_threshold,
			const int &pvalue_for_RD_test_choice,
		    type_map_string_to_Summary_Statistics &summary_stats_per_individual,
		    Summary_statistics &summary_stats_across_ALL_individuals)
{
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
		type_map_string_to_Summary_Statistics::iterator it_summ
				= summary_stats_per_individual.insert(std::pair<std::string, Summary_statistics>(it_mygen->first, Summary_statistics())).first;
		for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{					
			bool passes_threshold;
			if (pvalue_for_RD_test_choice == 0)
			{
				passes_threshold = mpfr::log10(it_myev->second.RD_odds_alt) <= logodds_threshold;				
			}
			else if (pvalue_for_RD_test_choice == 1)
			{
				passes_threshold = mpfr::log10(it_myev->second.Wilcoxon_signed_rank_sum_test__pvalue__inferred) <= logodds_threshold;
			}
			else if (pvalue_for_RD_test_choice == 2)
			{
				passes_threshold = mpfr::log10(it_myev->second.wilcoxon_matched_info___null___inferred.wilcoxon_pvalue) <= logodds_threshold;
			}			
			
			if (passes_threshold)
			{
				it_summ->second.add_call_appropriately(Conn_Comps, it_myev->second);
				summary_stats_across_ALL_individuals.add_call_appropriately(Conn_Comps, it_myev->second);	    
			}
		}//it_myev
    }//it_mygen
    
    
}//compile_Summary_statistics













Theoretical_compute_space compile_Theoretical_compute_space
	    (const type_map_uint_to_CC &Conn_Comps,
	    const type_map_uint_to_BI_to_Gene &Genes_in_the_human_genome)
{
    
    Theoretical_compute_space the_theoretical_space;
    
    for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
	    it_cc != Conn_Comps.end();
	    ++it_cc)
    {
	++the_theoretical_space.number_of_Connected_Components_tested;
	
	for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
		it_ev != it_cc->second.events.end();
		++it_ev)
	{	
	    ++the_theoretical_space.number_of_Events_tested_diploid;
	    
	    if (it_ev->second.recomb_type == recomb_class__DupDel)
	    {
		++the_theoretical_space.number_of_deletions_tested__diploid;
		++the_theoretical_space.number_of_duplications_tested__diploid;
		
	
	    }//dupdel
	    else
	    {
		++the_theoretical_space.number_of_inversions_tested__diploid;				
	    }
	    
	    
	    ++the_theoretical_space.number_of_GeneConversion_tested__diploid;
	    
	    the_theoretical_space.amount_of_haploid_Reference_genome_examined_for_potential_events
			+= BOOST_width_inclusive(it_ev->second.region_between_and_including_the_LCRs_themselves);
						
	    
	    const type_set_string potentially_affected_Genes_this_Event(
			extract_keys_of_map_and_return_as_set<std::string, Gene>(
						identify_Genes_potentially_affected_by_an_Event(it_ev->second, Genes_in_the_human_genome)));
	    
	    the_theoretical_space.potentially_affected_Genes.insert(potentially_affected_Genes_this_Event.begin(), potentially_affected_Genes_this_Event.end());
			
	}//ev		
    }//cc        
    
    
    return the_theoretical_space;
    
    
}//compile_Theoretical_compute_space












void save_Theoretical_space_to_table
			    (const std::string &outdir,
			    const Theoretical_compute_space &the_theoretical_space)
{
    std::string outfname(outdir);
    outfname.append("/theoretical_space.txt");
    
    std::ofstream outfs(outfname);
    if (check_fstream_for_goodness(outfs, outfname, "save_Theoretical_space_to_table", true))
    {
	outfs 
	    << "Inversions\t"
	    << "Deletions\t"
	    << "Duplications\t"
	    << "Gene Conversion\t"
	    << "Potentially affected Genes\t"
	    << "Amount of haploid Reference DNA examined\t"
	    << "Connected Components\t"
	    << "Events\n"    
	    

	    << the_theoretical_space.number_of_inversions_tested__diploid << "\t"
	    << the_theoretical_space.number_of_deletions_tested__diploid << "\t"
	    << the_theoretical_space.number_of_duplications_tested__diploid << "\t"
	    << the_theoretical_space.number_of_GeneConversion_tested__diploid << "\t"
	    << the_theoretical_space.potentially_affected_Genes.size() << "\t"
	    << the_theoretical_space.amount_of_haploid_Reference_genome_examined_for_potential_events << "\t"
	    << the_theoretical_space.number_of_Connected_Components_tested << "\t"
	    << the_theoretical_space.number_of_Events_tested_diploid << "\n";	    
	    
	outfs.close();
    }
    
}//save_Theoretical_space_to_table













void write_line_of_Summary_stats
	    (std::ofstream &outfs,
	    const Theoretical_compute_space &the_theoretical_space,
	    const std::string &name_of_line,
	    const Summary_statistics &some_Summary_statistics,
	    const Affected_Gene_catalog &corresponding_Affected_genes,
	    const uint &average_divider)
{    
    
    type_map_haploid_outcome_to_uint map_hapout_to_counts___Non_NULL(some_Summary_statistics.map_haploid_outcome_to_counts);
    map_hapout_to_counts___Non_NULL.erase(hap_outcome__None);
    const uint number_of_called_NAHR_or_GeneConv = map_hapout_to_counts___Non_NULL.size();
    
    type_map_haploid_outcome_to_uint map_hapout_to_counts___Non_NULL__NAHR(map_hapout_to_counts___Non_NULL);
    map_hapout_to_counts___Non_NULL__NAHR.erase(hap_outcome__GeneConv_ABA);
    map_hapout_to_counts___Non_NULL__NAHR.erase(hap_outcome__GeneConv_BAB);
    const uint number_of_called_NAHR = map_hapout_to_counts___Non_NULL__NAHR.size();
    
    type_map_haploid_outcome_to_uint map_hapout_to_counts___Non_NULL__GeneConv(map_hapout_to_counts___Non_NULL);
    map_hapout_to_counts___Non_NULL__GeneConv.erase(hap_outcome__Inv);
    map_hapout_to_counts___Non_NULL__GeneConv.erase(hap_outcome__Del);    
    map_hapout_to_counts___Non_NULL__GeneConv.erase(hap_outcome__Dup);  
    const uint number_of_called_GeneConv = map_hapout_to_counts___Non_NULL__GeneConv.size();
    
    
    const uint num_potential_NAHR = the_theoretical_space.number_of_deletions_tested__diploid 
				+ the_theoretical_space.number_of_duplications_tested__diploid 
				+ the_theoretical_space.number_of_inversions_tested__diploid;
    
    
    type_set_string genes_affected_by_NAHR;
//     extract_map_keys_and_insert_to_set<std::string, Gene>(corresponding_Affected_genes.Genes_affected_by_Inversions, genes_affected_by_NAHR);
    extract_map_keys_and_insert_to_set<std::string, Gene>(corresponding_Affected_genes.Genes_affected_by_Deletions, genes_affected_by_NAHR);
    extract_map_keys_and_insert_to_set<std::string, Gene>(corresponding_Affected_genes.Genes_affected_by_Duplications, genes_affected_by_NAHR);

    
    type_set_string genes_affected_by_GeneConv;
    extract_map_keys_and_insert_to_set<std::string, Gene>(corresponding_Affected_genes.Genes_affected_by_GeneConversion__donated_itself_to_elswhere, genes_affected_by_GeneConv);
    extract_map_keys_and_insert_to_set<std::string, Gene>(corresponding_Affected_genes.Genes_affected_by_GeneConversion__received_from_somewhere_else, genes_affected_by_GeneConv);

        
    type_set_string genes_affected_by_NAHR_or_GeneConv(genes_affected_by_NAHR);
    genes_affected_by_NAHR_or_GeneConv.insert(genes_affected_by_GeneConv.begin(), genes_affected_by_GeneConv.end());    
    
    const uint num_inv = (some_Summary_statistics.map_haploid_outcome_to_counts.count(hap_outcome__Inv) > 0)  ?    
			    some_Summary_statistics.map_haploid_outcome_to_counts.at(hap_outcome__Inv)  :  0;

    const uint num_del = (some_Summary_statistics.map_haploid_outcome_to_counts.count(hap_outcome__Del) > 0)  ?    
			    some_Summary_statistics.map_haploid_outcome_to_counts.at(hap_outcome__Del)  :  0;
				
    const uint num_dup = (some_Summary_statistics.map_haploid_outcome_to_counts.count(hap_outcome__Dup) > 0)  ?    
			    some_Summary_statistics.map_haploid_outcome_to_counts.at(hap_outcome__Dup)  :  0;
				
    uint num_GeneConv = (some_Summary_statistics.map_haploid_outcome_to_counts.count(hap_outcome__GeneConv_ABA) > 0)  ?    
			    some_Summary_statistics.map_haploid_outcome_to_counts.at(hap_outcome__GeneConv_ABA)  :  0;
				
    if (some_Summary_statistics.map_haploid_outcome_to_counts.count(hap_outcome__GeneConv_BAB) > 0)
    {  some_Summary_statistics.map_haploid_outcome_to_counts.at(hap_outcome__GeneConv_BAB);  }
				
				
				
    const uint number_theoretical_NAHR
			= the_theoretical_space.number_of_inversions_tested__diploid 
			+ the_theoretical_space.number_of_deletions_tested__diploid 
			+ the_theoretical_space.number_of_duplications_tested__diploid;
			
    type_set_uint ccs_with_some_call(some_Summary_statistics.connected_components_with_some_NAHR_called_positive__haploid);
    ccs_with_some_call.insert(some_Summary_statistics.connected_components_with_some_GeneConversion_called_positive__haploid.begin(),
				some_Summary_statistics.connected_components_with_some_GeneConversion_called_positive__haploid.end());    
    
    
    
    outfs 
	<< name_of_line << "\t"
	    
// 	    << "Number of Inversions\t\t"
// 	    << "Number of Deletions\t\t"
// 	    << "Number of Duplications\t\t"  	
// 	    << "Number of Gene Conversions\t"
	
	<< ((longdouble)num_inv/average_divider) << "\t(" << ((longdouble)num_inv/average_divider)/the_theoretical_space.number_of_inversions_tested__diploid << ")\t"
	<< ((longdouble)num_del/average_divider) << "\t(" << ((longdouble)num_del/average_divider)/the_theoretical_space.number_of_deletions_tested__diploid << ")\t"
	<< ((longdouble)num_dup/average_divider) << "\t(" << ((longdouble)num_dup/average_divider)/the_theoretical_space.number_of_duplications_tested__diploid << ")\t"
	<< ((longdouble)num_GeneConv/average_divider) << "\t(" << ((longdouble)num_GeneConv/average_divider)/the_theoretical_space.number_of_GeneConversion_tested__diploid << ")\t"
	
	
// 	    << "Number of genes affected by NAHR\t\t"	    
// 	    << "Number of genes affected by Gene Conversion\t\t"
// 	    << "Number of genes affected by NAHR or Gene Conversion\t\t"  
	
	<< ((longdouble)genes_affected_by_NAHR.size()/average_divider) << "\t(" << ((longdouble)genes_affected_by_NAHR.size()/average_divider)/the_theoretical_space.potentially_affected_Genes.size() << ")\t"
	<< ((longdouble)genes_affected_by_GeneConv.size()/average_divider) << "\t(" << ((longdouble)genes_affected_by_GeneConv.size()/average_divider)/the_theoretical_space.potentially_affected_Genes.size() << ")\t"
	<< ((longdouble)genes_affected_by_NAHR_or_GeneConv.size()/average_divider) << "\t(" << ((longdouble)genes_affected_by_NAHR_or_GeneConv.size()/average_divider)/the_theoretical_space.potentially_affected_Genes.size() << ")\t"
	
	
// 	    << "Total DNA affected by NAHR\t\t"
// 	    << "Total DNA affected by Gene Conversion\t\t"
// 	    << "Total DNA affected by NAHR or Gene Conversion\t\t"    
	
	<< ((longdouble)some_Summary_statistics.total_DNA_affected_by_NAHR___diploid/average_divider) << "\t(" << ((longdouble)some_Summary_statistics.total_DNA_affected_by_NAHR___diploid/average_divider)/length_of_haploid_genome << ")\t"
	
	<< ((longdouble)some_Summary_statistics.total_DNA_affected_by_GeneConversion___diploid/average_divider) << "\t(" << ((longdouble)some_Summary_statistics.total_DNA_affected_by_GeneConversion___diploid/average_divider)/length_of_haploid_genome << ")\t"
	
	<< ((longdouble)(some_Summary_statistics.total_DNA_affected_by_NAHR___diploid + some_Summary_statistics.total_DNA_affected_by_GeneConversion___diploid))/average_divider << "\t(" << ((longdouble)(some_Summary_statistics.total_DNA_affected_by_NAHR___diploid + some_Summary_statistics.total_DNA_affected_by_GeneConversion___diploid)/average_divider)/length_of_haploid_genome << ")\t"
	
	
// 	    << "Number of Events experiencing NAHR\t\t"
// 	    << "Number of Events experiencing Gene Conversion\t\t"
// 	    << "Number of Events experiencing NAHR or Gene Conversion\t\t"
	
	<< ((longdouble)number_of_called_NAHR/average_divider) << "\t(" << (longdouble)number_of_called_NAHR/number_theoretical_NAHR << ")\t"
	<< ((longdouble)number_of_called_GeneConv/average_divider) << "\t(" << (longdouble)number_of_called_GeneConv/the_theoretical_space.number_of_GeneConversion_tested__diploid << ")\t"
	<< ((longdouble)(number_of_called_NAHR + number_of_called_GeneConv)/average_divider) << "\t(" << ((longdouble)(number_of_called_NAHR + number_of_called_GeneConv)/average_divider)/(number_theoretical_NAHR + the_theoretical_space.number_of_GeneConversion_tested__diploid) << ")\t"	    
	
	
// 	    << "Number of Connected Components experiencing NAHR\t\t"
// 	    << "Number of Connected Components experiencing Gene Conversion\t\t"
// 	    << "Number of Connected Components experiencing NAHR or Gene Conversion\t\t\n";

	<< (longdouble)some_Summary_statistics.connected_components_with_some_NAHR_called_positive__haploid.size()/average_divider << "\t(" << ((longdouble)some_Summary_statistics.connected_components_with_some_NAHR_called_positive__haploid.size()/average_divider)/the_theoretical_space.number_of_Connected_Components_tested
	<< (longdouble)some_Summary_statistics.connected_components_with_some_GeneConversion_called_positive__haploid.size()/average_divider << "\t(" << ((longdouble)some_Summary_statistics.connected_components_with_some_GeneConversion_called_positive__haploid.size()/average_divider)/the_theoretical_space.number_of_Connected_Components_tested
	<< (longdouble)ccs_with_some_call.size()/average_divider << "\t(" << ((longdouble)ccs_with_some_call.size()/average_divider)/the_theoretical_space.number_of_Connected_Components_tested << ")\n";
		    
            
}//write_line_of_Summary_stats



















void save_Summary_statistics_to_table
			    (const std::string &outdir,
			     const std::string &addl_id,
			    const Theoretical_compute_space &the_theoretical_space,
			    const type_map_string_to_Summary_Statistics &summary_stats_per_individual,
			    const Summary_statistics &summary_stats_across_ALL_individuals,
			    const type_map_string_to_Affected_Gene_catalog  &affected_Genes_by_individual)
{
    std::string outfname(outdir);
    outfname.append("/summary_statistics.");
    outfname.append(addl_id);
    outfname.append(".txt");
    
    std::ofstream outfs(outfname);
    if (check_fstream_for_goodness(outfs, outfname, "save_Summary_statistics_to_table", true))
    {
	outfs 
	    << "Sample name\t"
	    << "Number of Inversions\t\t"
	    << "Number of Deletions\t\t"
	    << "Number of Duplications\t\t"  	
	    << "Number of Gene Conversions\t\t"
	    
	    << "Number of genes affected by NAHR\t\t"	    
	    << "Number of genes affected by Gene Conversion\t\t"
	    << "Number of genes affected by NAHR or Gene Conversion\t\t"
	    
	    << "Total DNA affected by NAHR\t\t"
	    << "Total DNA affected by Gene Conversion\t\t"
	    << "Total DNA affected by NAHR or Gene Conversion\t\t"
	    
	    << "Number of Events experiencing NAHR\t\t"
	    << "Number of Events experiencing Gene Conversion\t\t"
	    << "Number of Events experiencing NAHR or Gene Conversion\t\t"
	    
	    << "Number of Connected Components experiencing NAHR\t\t"
	    << "Number of Connected Components experiencing Gene Conversion\t\t"
	    << "Number of Connected Components experiencing NAHR or Gene Conversion\t\t\n";
	    
	    
	    type_set_uint affected_Events;
	    uint number_of_affectedgenes_total = 0;
	    
	    //Individual summaries:
	    for (type_map_string_to_Summary_Statistics::const_iterator it_sum = summary_stats_per_individual.begin();
		    it_sum != summary_stats_per_individual.end();
		    ++it_sum)
	    {
		write_line_of_Summary_stats(
				outfs,
				the_theoretical_space,
			       it_sum->first,
			      it_sum->second,
			      affected_Genes_by_individual.at(it_sum->first));
	    }//it_sum
	    
	    //Totals:
	    Affected_Gene_catalog affected_genes_on_ALL_individuals;
	    affected_genes_on_ALL_individuals.collect_from_individual_catalogs(affected_Genes_by_individual);
	    
	    write_line_of_Summary_stats(
			    outfs,
			    the_theoretical_space,
			    "Total",
			    summary_stats_across_ALL_individuals,
			    affected_genes_on_ALL_individuals);
	    
	    //Average:
	    write_line_of_Summary_stats(
			    outfs,
			    the_theoretical_space,
			    "Average",
			    summary_stats_across_ALL_individuals,
			    affected_genes_on_ALL_individuals,
			    summary_stats_per_individual.size());	
	    
	    
	outfs 
	    << "\n\n"
	    << "number of distinct events:\t"
	    << affected_genes_on_ALL_individuals.Events_occurring_as_Inversions.size() << "\t("
	    << (double)affected_genes_on_ALL_individuals.Events_occurring_as_Inversions.size()/the_theoretical_space.number_of_inversions_tested__diploid << ")\t"
	    << affected_genes_on_ALL_individuals.Events_occurring_as_Deletions.size() << "\t("
	    << (double)affected_genes_on_ALL_individuals.Events_occurring_as_Deletions.size()/the_theoretical_space.number_of_deletions_tested__diploid << ")\t"
	    << affected_genes_on_ALL_individuals.Events_occurring_as_Duplications.size() << "\t("
	    << (double)affected_genes_on_ALL_individuals.Events_occurring_as_Duplications.size()/the_theoretical_space.number_of_duplications_tested__diploid << ")\t"
	    << affected_genes_on_ALL_individuals.Events_occurring_as_GeneConversions.size() << "\t("
	    << (double)affected_genes_on_ALL_individuals.Events_occurring_as_GeneConversions.size()/the_theoretical_space.number_of_GeneConversion_tested__diploid << ")\n"
	    
	    <<"number of distinct genes affected:\t"
	    << affected_genes_on_ALL_individuals.Genes_affected_by_Inversions.size() << "\t("
	    << (double)affected_genes_on_ALL_individuals.Genes_affected_by_Inversions.size()/the_theoretical_space.potentially_affected_Genes.size() << ")\t"
	    << affected_genes_on_ALL_individuals.Genes_affected_by_Deletions.size() << "\t("
	    << (double)affected_genes_on_ALL_individuals.Genes_affected_by_Deletions.size()/the_theoretical_space.potentially_affected_Genes.size() << ")\t"
	    << affected_genes_on_ALL_individuals.Genes_affected_by_Duplications.size() << "\t("
	    << (double)affected_genes_on_ALL_individuals.Genes_affected_by_Duplications.size()/the_theoretical_space.potentially_affected_Genes.size() << ")\t"
	    << (affected_genes_on_ALL_individuals.Genes_affected_by_GeneConversion__donated_itself_to_elswhere.size() + affected_genes_on_ALL_individuals.Genes_affected_by_GeneConversion__received_from_somewhere_else.size()) << "\t("
	    << (double)(affected_genes_on_ALL_individuals.Genes_affected_by_GeneConversion__donated_itself_to_elswhere.size() + affected_genes_on_ALL_individuals.Genes_affected_by_GeneConversion__received_from_somewhere_else.size())/the_theoretical_space.potentially_affected_Genes.size() << ")\n";	    
	    
	    
	type_set_uint events_occuring__dels_and_dups(affected_genes_on_ALL_individuals.Events_occurring_as_Deletions);
	type_set_string genes_affected__dels_and_dups(extract_keys_of_map_and_return_as_set<std::string, Gene>( affected_genes_on_ALL_individuals.Genes_affected_by_Deletions));
	
	events_occuring__dels_and_dups.insert(affected_genes_on_ALL_individuals.Events_occurring_as_Duplications.begin(),
					    affected_genes_on_ALL_individuals.Events_occurring_as_Duplications.end());	    
	
	const type_set_string affected_genes__dups(extract_keys_of_map_and_return_as_set<std::string, Gene>(affected_genes_on_ALL_individuals.Genes_affected_by_Duplications));
	genes_affected__dels_and_dups.insert(affected_genes__dups.begin(), affected_genes__dups.end());
	
	outfs
	    << "events combined across dels and dups:\t"
	    << events_occuring__dels_and_dups.size() << "\t"
	    << "(" << (double)events_occuring__dels_and_dups.size() / the_theoretical_space.number_of_deletions_tested__diploid << ")\n"
	    
	    << "affected genes across dels and dups:\t"
	    << genes_affected__dels_and_dups.size() << "\t"
	    << "(" << (double)genes_affected__dels_and_dups.size() / the_theoretical_space.potentially_affected_Genes.size() << ")\n";
    

	    
	    	    
	outfs.close();
    }    
    
    
}//save_Summary_statistics_to_table














void Empirical_properties::save_to_file
		    (const std::string& outdir, 
		     const std::string& EP_name)  const
{
    std::string outrootname(outdir);
    outrootname.append("/empirical__");
    outrootname.append(EP_name);    
    
    std::string outfname = outrootname;
    outfname.append(".breakpoint_sparsity");    
    write_singlecontainer_to_file<type_list_uint>(sparsity_histogram_of_breakpoints, outfname);
    
    outfname = outrootname;
    outfname.append(".LCR_identity");    
    write_singlecontainer_to_file<type_list_longdouble>(mediating_LCR_overall_percent_identity_histogram, outfname);
    
    outfname = outrootname;
    outfname.append(".LCR_length");    
    write_singlecontainer_to_file<type_list_uint>(mediating_LCR_length_histogram, outfname);    
    
    outfname = outrootname;
    outfname.append(".LCR_distance");    
    write_singlecontainer_to_file<type_list_uint>(distance_between_mediating_LCRs_histogram, outfname);      
    
    outfname = outrootname;
    outfname.append(".MEPS");    
    write_singlecontainer_to_file<type_list_longdouble>(best_MEPS_homology_histogram, outfname);        
    
    outfname = outrootname;
    outfname.append(".affected_region");    
    write_singlecontainer_to_file<type_list_uint>(affected_region_size_histogram, outfname);       
    
    outfname = outrootname;
    outfname.append(".donor_tract_length");    
    write_singlecontainer_to_file<type_list_uint>(length_of_donor_GeneConversion_tract__histogram, outfname);
    
    outfname = outrootname;
    outfname.append(".receiver_tract_length");    
    write_singlecontainer_to_file<type_list_uint>(length_of_receiving_GeneConversion_tract__histogram, outfname);
    
    outfname = outrootname;
    outfname.append(".GC_tract_ratio");    
    write_singlecontainer_to_file<type_list_longdouble>(ratioGC_content_of_GeneConversion_tract__histogram, outfname);
    
    outfname = outrootname;
    outfname.append(".event_multiplicity");    
    write_singlecontainer_to_file<type_list_uint>(number_of_occurring_individuals__vs__number_of_occurring_Events__histogram, outfname);
    
    outfname = outrootname;
    outfname.append(".posterior_outcome");     
    write_singlecontainer_to_file<type_list_longdouble>(posterior_probabilities_of_MAP_event_outcome__histogram, outfname);
    
    outfname = outrootname;
    outfname.append(".posterior_breakpoint");     
    write_singlecontainer_to_file<type_list_longdouble>(posterior_probabilities_of_MAP_event_breakpoints__histogram, outfname);    
                  
	
}//save_to_file













// void write_raw_negative_calls_to_stream
// 	(const std::string &outfname,
// 	 const type_map_uint_to_CC &Conn_Comps,
// 	 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
// {
//     std::ofstream outfs(outfname);
//     check_fstream_for_goodness(outfs, outfname, "write_raw_negative_calls_to_stream", true);
//     
//     
//     outfs
// 	<< "Sample name" << "\t"
// 	<< "population" << "\t"
// 	<< "call type" << "\t"
// 	
// 	<< "chromosome" << "\t"
// 	<< "begin" << "\t"
// 	<< "end" << "\t"
// 	<< "affected region size" << "\t"
// 	
// 	<< "matching 1000 Genomes calls" << "\t"
// 	<< "genes affected" << "\t"
// 	
// 	<< "Read-depth hypothesis test p-value" << "\t"
// 	
// 	<< "P(outcome)" << "\t"
// 	<< "P(breakpoint)" << "\t"
// 	<< "connected component" << "\t"
// 	<< "event" << "\t"
// 	<< "profile brkpt" << "\n";
//     
//     
//     for (type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_mygen = mycalls.begin();
// 	    it_mygen != mycalls.end();
// 	    ++it_mygen)
//     {
// 	for (type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_mygen->second.begin();
// 		it_myev != it_mygen->second.end();
// 		++it_myev)
// 	{
// 	    for (uint hap=0; hap<2; ++hap)
// 	    {
// 		if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) == only_this_haploid_outcome)
// 		{		
// 		    const type_map_uint_to_Event::const_iterator the_Event_it =  Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID);
// 		    
// 		    const BOOST_Interval affected_region =  (only_this_haploid_outcome == hap_outcome__None)   ?  
// 								the_Event_it->second.region_between_and_including_the_LCRs_themselves
// 								:
// 								convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
// 									pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
// 									the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
// 									only_this_haploid_outcome).first;
// 		    
// 		    outfs
// 			<< it_mygen->first << "\t"
// 			<< it_myev->second.population_name << "\t"
// 			<< convert_haploid_outcome_to_string(only_this_haploid_outcome) << "\t"
// 			
// 			<< the_Event_it->second.chromos[0] << "\t"
// 			<< affected_region.lower() << "\t"
// 			<< affected_region.upper() << "\t"
// 			<< BOOST_width_inclusive(affected_region) << "\t"
// 			
// 			<< concatenate_singlecontainer_elements_into_string<type_set_string>(it_myev->second.associated_other_calls) << "\t"
// 			<< concatenate_singlecontainer_elements_into_string<type_set_string>(pair_at<type_set_string>(it_myev->second.affected_genes_per_haploid,hap)) << "\t"
// 			
// 			<< it_myev->second.RD_hypoth_tests_p_value<< "\t"
// 			
// 			<< it_myev->second.P_diploid_outcome << "\t"
// 			<< it_myev->second.P_diploid_breakpoint << "\t"
// 			<< it_myev->second.cc_of_event << "\t"
// 			<< it_myev->second.event_UID << "\t"
// 			<< pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).lower() << "\n";
// 		}//outcome filter
// 	    }//hap
// 		
// 	}//ymev
//     }//mygen
//     
// }//write_raw_negative_calls_to_stream
// 








void write_raw_NAHR_calls_to_stream
	(const std::string &outfname,
	 const type_map_uint_to_CC &Conn_Comps,
	 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
	const type_haploid_outcome &only_this_haploid_outcome)
{
    std::ofstream outfs(outfname);
    check_fstream_for_goodness(outfs, outfname, "write_raw_NAHR_calls_to_stream", true);
    
    
    outfs
	<< "Sample name" << "\t"
	<< "population" << "\t"
	<< "call type" << "\t"
	
	<< "chromosome" << "\t"
	<< "begin" << "\t"
	<< "end" << "\t"
	<< "affected region size" << "\t"
	<< "number unique positions in region" << "\t"
	<< "involves undefined regions" << "\t"
	
	<< "matching 1000 Genomes calls" << "\t"
	<< "genes affected" << "\t"
	
	<< "Read-depth hypothesis test log-odds (log[null/alternative])" << "\t"
	
	<< "Wilcoxon signed rank-sum test log-pvalue" << "\t"
	
	<< "P(outcome)" << "\t"
	<< "P(breakpoint)" << "\t"
	<< "connected component" << "\t"
	<< "event" << "\t"
	<< "profile brkpt" << "\n";
    
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
		for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{
			for (uint hap=0; hap<2; ++hap)
			{
				if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) == only_this_haploid_outcome)
				{		
					const type_map_uint_to_Event::const_iterator the_Event_it =  Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID);
					
					const BOOST_Interval affected_region =  (only_this_haploid_outcome == hap_outcome__None)   ?  
										the_Event_it->second.region_between_and_including_the_LCRs_themselves
										:
										convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
											pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
											the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
											only_this_haploid_outcome).first;
					
					outfs
					<< it_mygen->first << "\t"
					<< it_myev->second.population_name << "\t"
					<< convert_haploid_outcome_to_string(only_this_haploid_outcome) << "\t"
					
					<< the_Event_it->second.chromos[0] << "\t"
					<< affected_region.lower() << "\t"
					<< affected_region.upper() << "\t"
					<< BOOST_width_inclusive(affected_region) << "\t"
					<< it_myev->second.number_of_unique_bases_in_region << "\t";
					
					if (it_myev->second.from_a_CC_with_undefined_regions)
					{  outfs << "yes\t";  }
					else
					{  outfs << "no\t";  }
					
					outfs
					<< concatenate_singlecontainer_elements_into_string<type_set_string>(it_myev->second.associated_other_calls) << "\t"
					<< concatenate_singlecontainer_elements_into_string<type_set_string>(pair_at<type_set_string>(it_myev->second.affected_genes_per_haploid,hap)) << "\t";
					
					if (only_this_haploid_outcome == hap_outcome__None)
						outfs << mpfr::log10(mpfr::max(it_myev->second.RD_odds_del, it_myev->second.RD_odds_dup)) << "\t";
					else
						outfs << mpfr::log10(it_myev->second.RD_odds_alt) << "\t";
					
					
					outfs
					<< it_myev->second.Wilcoxon_signed_rank_sum_test__pvalue__inferred << "\t"
					<< it_myev->second.P_diploid_outcome << "\t"
					<< it_myev->second.P_diploid_breakpoint << "\t"
					<< it_myev->second.cc_of_event << "\t"
					<< it_myev->second.event_UID << "\t"
					<< pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).lower() << "\n";
				}//outcome filter
			}//hap
			
		}//ymev
    }//mygen
    
    
    outfs.close();
    
}//write_raw_NAHR_calls_to_stream











void write_raw_GeneConversion_calls_to_stream
	(const std::string &outfname,
	 const type_map_uint_to_CC &Conn_Comps,
	 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
    std::ofstream outfs(outfname);
    check_fstream_for_goodness(outfs, outfname, "write_raw_GeneConversion_calls_to_stream", true);
    
    
    outfs
	<< "Sample name" << "\t"
	<< "population" << "\t"
	<< "call type" << "\t"
	
	<< "donor chromosome" << "\t"
	<< "donor begin" << "\t"
	<< "donor end" << "\t"
	<< "donor size" << "\t"
	
	<< "receiver chromosome" << "\t"
	<< "receiver begin" << "\t"
	<< "receiver end" << "\t"
	<< "receiver size" << "\t"
	
	<< "number unique positions in region" << "\t"
	<< "involves undefined regions" << "\t"	
	
	<< "genes affected" << "\t"
	
// 	<< "Read-depth hypothesis test p-value" << "\t"
	
	<< "P(outcome)" << "\t"
	<< "P(breakpoint)" << "\t"
	<< "connected component" << "\t"
	<< "event" << "\t"
	<< "profile brkpt low" << "\t"
	<< "profile brkpt high" << "\n";
    
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
	for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
		it_myev != it_mygen->second.end();
		++it_myev)
	{
	    for (uint hap=0; hap<2; ++hap)
	    {
		if (test_if_haploid_outcome_is_GeneConversion(pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap)))
		{		
		    const type_map_uint_to_Event::const_iterator the_Event_it =  Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID);
		    
		    const type_BI__BI affected_region =  (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) == hap_outcome__None)   ?  
								type_BI__BI(the_Event_it->second.region_between_and_including_the_LCRs_themselves, the_Event_it->second.region_between_and_including_the_LCRs_themselves)
								:
								convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
									pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
									the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
									pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap));
								
		    const BOOST_Interval donor_region =   pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) == hap_outcome__GeneConv_ABA  ?  
								affected_region.second  :  affected_region.first;    
								
		    const BOOST_Interval receiver_region =   pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap) == hap_outcome__GeneConv_ABA  ?  
								affected_region.first  :  affected_region.second;
								
		    
		    outfs
			<< it_mygen->first << "\t"
			<< it_myev->second.population_name << "\t"
			<< "GeneConversion" << "\t"
			
			<< the_Event_it->second.chromos[0] << "\t"
			<< donor_region.lower() << "\t"
			<< donor_region.upper() << "\t"
			<< BOOST_width_inclusive(donor_region) << "\t"
			
			<< the_Event_it->second.chromos[0] << "\t"
			<< receiver_region.lower() << "\t"
			<< receiver_region.upper() << "\t"
			<< BOOST_width_inclusive(receiver_region) << "\t"
			
			<< it_myev->second.number_of_unique_bases_in_region << "\t";
			
			if (it_myev->second.from_a_CC_with_undefined_regions)
			{  outfs << "yes\t";  }
			else
			{  outfs << "no\t";  }			
			
		    outfs			
			<< concatenate_singlecontainer_elements_into_string<type_set_string>(pair_at<type_set_string>(it_myev->second.affected_genes_per_haploid,hap)) << "\t"
			
// 			<< it_myev->second.RD_hypoth_tests_p_value<< "\t"
			
			<< it_myev->second.P_diploid_outcome << "\t"
			<< it_myev->second.P_diploid_breakpoint << "\t"
			<< it_myev->second.cc_of_event << "\t"
			<< it_myev->second.event_UID << "\t"
			<< pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).lower() << "\t"
			<< pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap).upper() << "\n";
			
		}//outcome filter
	    }//hap
		
	}//ymev
    }//mygen
    
    outfs.close();
    
}//write_raw_GeneConversion_calls_to_stream















void save_raw_positive_calls_to_files
	(const std::string &outdir,
	 const std::string &addl_id,
	 const type_map_uint_to_CC &Conn_Comps,
	const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
    std::string outrootname(output_dir);
    outrootname.append("/My_calls.");
    outrootname.append(addl_id);
    
        
    std::string outfname;
    
    
    outfname = outrootname;
    outfname.append(".none.txt");        
    write_raw_NAHR_calls_to_stream(
	outfname,
	 Conn_Comps,
	 mycalls,
	hap_outcome__None);
    
    outfname = outrootname;
    outfname.append(".inversions.txt");        
    write_raw_NAHR_calls_to_stream(
	outfname,
	 Conn_Comps,
	 mycalls,
	hap_outcome__Inv);        
    
    
    outfname = outrootname;
    outfname.append(".deletions.txt");        
    write_raw_NAHR_calls_to_stream(
	outfname,
	 Conn_Comps,
	 mycalls,
	hap_outcome__Del);
    
    
    outfname = outrootname;
    outfname.append(".duplications.txt");       
    write_raw_NAHR_calls_to_stream(
	outfname,
	 Conn_Comps,
	 mycalls,
	hap_outcome__Dup);    
    
    
    outfname = outrootname;
    outfname.append(".geneconversions.txt");       
    write_raw_GeneConversion_calls_to_stream(
	 outfname,
	 Conn_Comps,
	 mycalls);
    
    
}//save_raw_positive_calls_to_files















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
		const type_map_string_to_uint_to_Called_diploid_wrapper  &novel_mycalls)
{
    save_comparisons_to_tables___imp
		(outdir,
		 "full",
		real(mpfr::const_infinity()),
		 pvalue_for_RD_test_choice,
		Conn_Comps,
		mycalls,
		validations_from_other_studies,
		Vals__unmatched,
		map_genome_to_valname_to_calls__correctly_matched,
		map_genome_to_valname_to_calls__alternatively_matched,
		novel_mycalls);
		
    save_comparisons_to_tables___imp
		(outdir,
		 "thresh",
		logodds_threshold,
		pvalue_for_RD_test_choice,
		Conn_Comps,
		mycalls,
		validations_from_other_studies,
		Vals__unmatched,
		map_genome_to_valname_to_calls__correctly_matched,
		map_genome_to_valname_to_calls__alternatively_matched,
		novel_mycalls);	
		
		
    mpfr_free_cache();
    
}//save_comparisons_to_tables













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
		const type_map_string_to_uint_to_Called_diploid_wrapper  &novel_mycalls)
{
    typedef  std::map<std::string, Validated_rearrangement> type_map_string_to_Validated_rearrangement;
    typedef  std::map<std::string, type_map_string_to_Validated_rearrangement>   type_map_string_to_string_to_Validated_rearrangement;
    
    bool use_threshold = mpfr::isfinite(logodds_threshold);

    uint number_below_theshold__match = 0;
    uint number_below_theshold__UNMATCH = 0;
    uint number_below_theshold__alternative = 0;
   
    
    type_map_string_to_string_to_Validated_rearrangement map_genome_to_callID_to_Validated_rearrangement;
    for (type_list_Validated_rearrangement::const_iterator it_val = validations_from_other_studies.begin();
	    it_val != validations_from_other_studies.end();
	    ++it_val)
    {
	map_genome_to_callID_to_Validated_rearrangement[it_val->NA_subject][it_val->variant_accession_name] = *it_val;	
    }//val
    
    
     const type_map_uint_to_list_BI nonunique_regions(load_Non_unique_regions_of_genome());
    
    
                                                                          
    {//unmatched
	std::string outfname(outdir);
	outfname.append("/validations.");
	outfname.append(addl_id);
	outfname.append(".unmatched.txt");
	
	std::ofstream outfs(outfname);
	check_fstream_for_goodness(outfs, outfname, "save_comparisons_to_tables", true);		
	
	outfs
	    << "Sample name" << "\t"
	    << "publication" << "\t"
	    << "validated call ID" << "\t"	    
	    << "call type" << "\t"
	    << "chromosome" << "\t"
	    << "begin" << "\t"
	    << "end" << "\t"
	    << "number unique positions in region" << "\t"
	    << "related event UIDs" << "\t"
	    << "log odds (log [P null/ P deletion])" << "\t"
		<< "log Wilcoxon signed rank-sum test p-value" << "\t"
		<< "log Wilcoxon MATCHED signed rank-sum test p-value" << "\n";
		
	    
	    
	std::string outhistfname(outdir);
	outhistfname.append("/validations.");
	outhistfname.append(addl_id);
	outhistfname.append(".unmatched.hist");
	
	std::ofstream outhistfs(outhistfname);
	check_fstream_for_goodness(outhistfs, outhistfname, "save_comparisons_to_tables", true);	    
	    
	
    
	for (type_vector_list_Validated_rearrangement_constiterator::const_iterator it_unmatch = Vals__unmatched.begin();
		it_unmatch != Vals__unmatched.end();
		++it_unmatch)
	{	  
	    if (!(*it_unmatch)->related_UIDs.empty())
	    {
			const real logodds(mpfr::log10(mycalls.at((*it_unmatch)->NA_subject).at(*(*it_unmatch)->related_UIDs.begin()).RD_odds_del));
			const real logWilcoxon(mpfr::log10(mycalls.at((*it_unmatch)->NA_subject).at(*(*it_unmatch)->related_UIDs.begin()).Wilcoxon_signed_rank_sum_test__pvalue__inferred));
			const real logWilcoxon_matched(mpfr::log10(mycalls.at((*it_unmatch)->NA_subject).at(*(*it_unmatch)->related_UIDs.begin()).wilcoxon_matched_info___null___inferred.wilcoxon_pvalue));
			bool passes_threshold;
			
			if (pvalue_for_RD_test_choice == 0)
			{
				passes_threshold = logodds <= logodds_threshold;				
			}
			else if (pvalue_for_RD_test_choice == 1)
			{
				passes_threshold = logWilcoxon <= logodds_threshold;
			}			
			else if (pvalue_for_RD_test_choice == 2)
			{
				passes_threshold = logWilcoxon_matched <= logodds_threshold;
			}				
		
		if (passes_threshold or !use_threshold)
		{
	    
		    outfs
			<< (*it_unmatch)->NA_subject << "\t"
			<< (*it_unmatch)->published_study << "\t"
			<< (*it_unmatch)->variant_accession_name << "\t"		
			<< convert_haploid_outcome_to_string(static_cast<type_haploid_outcome>((*it_unmatch)->outcome)) << "\t"		
			<< (*it_unmatch)->chromosome << "\t"
			<< (*it_unmatch)->breakpoints.lower() << "\t"
			<< (*it_unmatch)->breakpoints.upper() << "\t"
			<< BOOST_width_inclusive((*it_unmatch)->breakpoints) - calculate_NON_uniqueness_of_region((*it_unmatch)->chromosome, (*it_unmatch)->breakpoints, nonunique_regions) << "\t"
			<< concatenate_singlecontainer_elements_into_string<type_set_uint>((*it_unmatch)->related_UIDs) << "\t"
			<< logodds << "\t"
			<< logWilcoxon << "\t"
			<< logWilcoxon_matched << "\n";
		 
			if (pvalue_for_RD_test_choice == 0)
				outhistfs << logodds << "\n";
			else if (pvalue_for_RD_test_choice == 1)
				outhistfs << logWilcoxon << "\n";
			else if (pvalue_for_RD_test_choice == 2)
				outhistfs << logWilcoxon_matched << "\n";
		
		    ++number_below_theshold__UNMATCH;
		}//threshold
	    }			    
	}//unmatched
	
	
	outfs << "\n\nnumber of unmatched validations with log-odds below " << logodds_threshold << ":\t" << number_below_theshold__UNMATCH << "\n";
	
	outfs.close();
	outhistfs.close();
    }//unmatched
    
    
    
    
		
    
    
    
    
    {//correctly matched
	std::string outfname(outdir);
	outfname.append("/validations.");
	outfname.append(addl_id);
	outfname.append(".matched.txt");
	
	std::ofstream outfs(outfname);
	check_fstream_for_goodness(outfs, outfname, "save_comparisons_to_tables", true);
	
	outfs
	    << "Sample name" << "\t"
	    << "publication" << "\t"
	    << "validated call ID" << "\t"	    
	    << "call type" << "\t"
	    << "chromosome" << "\t"
	    << "begin" << "\t"
	    << "end" << "\t"
	    
	    << "my call begin" << "\t"
	    << "my call end" << "\t"
	    << "amount overlap" << "\t"
	    << "amount non-overlap" << "\t"
	    << "overlap as percent of hull" << "\t"	  
	    
	    << "my connected component ID" << "\t"
	    << "my event UID" << "\t"
	    
	    << "log odds (log [P null/ P alternative])" << "\t"
		<< "log Wilcoxon" << "\t"
		<< "log Wilcoxon matched" << "\n";
	    
	    
	
	std::string outhistfname(outdir);
	outhistfname.append("/validations.");
	outhistfname.append(addl_id);
	outhistfname.append(".matched.hist");
	
	std::ofstream outhistfs(outhistfname);
	check_fstream_for_goodness(outhistfs, outhistfname, "save_comparisons_to_tables", true);
	    
	    
	for (type_map_string_to_string_to_list_Called_diploid_wrapper_ptr::const_iterator it_val_gen = map_genome_to_valname_to_calls__correctly_matched.begin();
		it_val_gen != map_genome_to_valname_to_calls__correctly_matched.end();
		++it_val_gen)
	{
	    for (type_map_string_to_list_Called_diploid_wrapper_ptr::const_iterator it_match_to_calls = it_val_gen->second.begin();
		    it_match_to_calls != it_val_gen->second.end();
		    ++it_match_to_calls)
	    {
		const type_map_string_to_Validated_rearrangement::const_iterator it_match = map_genome_to_callID_to_Validated_rearrangement.at(it_val_gen->first).find(it_match_to_calls->first);
		
		type_list_Called_diploid_wrapper_ptr::const_iterator it_best_matched_mycall = it_match_to_calls->second.end();
		longdouble best_percent_overlap_vs_hull = -1;
		uint best_hap;
		BOOST_Interval best_affected_region;
		
		for (type_list_Called_diploid_wrapper_ptr::const_iterator it_my_match = it_match_to_calls->second.begin();
			it_my_match != it_match_to_calls->second.end();
			++it_my_match)
		{
		    const type_map_uint_to_Event::const_iterator the_Event_it = Conn_Comps.at((*it_my_match)->cc_of_event).events.find((*it_my_match)->event_UID);
		    
		    for (uint hap=0; hap<2; ++hap)
		    {
			if (pair_at<type_haploid_outcome>((*it_my_match)->the_diploid_Event_outcome,hap) == static_cast<type_haploid_outcome>(it_match->second.outcome))
			{
			    const BOOST_Interval affected_region(		
					convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
						    pair_at<BOOST_Interval>((*it_my_match)->the_diploid_profile_brkpts,hap),
						    the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
						    static_cast<type_haploid_outcome>(it_match->second.outcome)).first);
			    
			    const longdouble percent_overlap = BOOST_width_inclusive(BOOST_intersect(it_match->second.breakpoints, affected_region));
			    
			    if (percent_overlap > best_percent_overlap_vs_hull)
			    {
				it_best_matched_mycall = it_my_match;
				best_percent_overlap_vs_hull = percent_overlap;
				best_hap = hap;
				best_affected_region = affected_region;
			    }
			}//same call
		    }//hap
		    
		}//my matches
		
		    
	
		uint amount_nonoverlap = 0;		
		const BOOST_Interval overlap_hull = BOOST_hull(it_match->second.breakpoints, best_affected_region);
		
		if (!BOOST_subset(it_match->second.breakpoints, best_affected_region))
		{
		    amount_nonoverlap += BOOST_width_inclusive(overlap_hull) - BOOST_width_inclusive(best_affected_region);		    
		}
		
		if (!BOOST_subset(best_affected_region, it_match->second.breakpoints))
		{
		    amount_nonoverlap += BOOST_width_inclusive(overlap_hull) - BOOST_width_inclusive(it_match->second.breakpoints);		    
		}
				
				
		const BOOST_Interval overlap_intersection(BOOST_intersect(best_affected_region, it_match->second.breakpoints));
	    		
		const real logodds(mpfr::log10((*it_best_matched_mycall)->RD_odds_alt));
		const real logWilcoxon(mpfr::log10((*it_best_matched_mycall)->Wilcoxon_signed_rank_sum_test__pvalue__inferred));
		const real logWilcoxon_matched(mpfr::log10((*it_best_matched_mycall)->wilcoxon_matched_info___null___inferred.wilcoxon_pvalue));
		bool passes_threshold;
		
		if (pvalue_for_RD_test_choice == 0)
		{
			passes_threshold = logodds <= logodds_threshold;				
		}
		else if (pvalue_for_RD_test_choice == 1)
		{
			passes_threshold = logWilcoxon <= logodds_threshold;
		}				
		else if (pvalue_for_RD_test_choice == 2)
		{
			passes_threshold = logWilcoxon_matched <= logodds_threshold;
		}
			
		
		if (passes_threshold or !use_threshold)
		{		
		    outfs
			<< it_match->second.NA_subject << "\t"
			<< it_match->second.published_study << "\t"
			<< it_match->second.variant_accession_name << "\t"		
			<< convert_haploid_outcome_to_string(static_cast<type_haploid_outcome>(it_match->second.outcome)) << "\t"		
			<< it_match->second.chromosome << "\t"
			<< it_match->second.breakpoints.lower() << "\t"
			<< it_match->second.breakpoints.upper() << "\t"
			
			<< best_affected_region.lower() << "\t"
			<< best_affected_region.upper() << "\t"
			<< BOOST_width_inclusive(overlap_intersection) << "\t"
			<< amount_nonoverlap << "\t"
			<< (longdouble)BOOST_width_inclusive(overlap_intersection)/BOOST_width_inclusive(overlap_hull) << "\t"		    
			
			<< (*it_best_matched_mycall)->cc_of_event << "\t"
			<< (*it_best_matched_mycall)->event_UID << "\t"
			<< logodds << "\t"
			<< logWilcoxon << "\t"
			<< logWilcoxon_matched << "\n";
		 
			if (pvalue_for_RD_test_choice == 0)
				outhistfs << logodds << "\n";
			else if (pvalue_for_RD_test_choice == 1)
				outhistfs << logWilcoxon << "\n";
			else if (pvalue_for_RD_test_choice == 2)
				outhistfs << logWilcoxon_matched << "\n";
			
			++number_below_theshold__match;
		}//thresh
		    
	    }//it_match	
	    
	}//it_val_gen
	
	outfs.close();
	outhistfs.close();
	    
    }//correctly matched
    
    
    
    
    
    
    
    

    
    {//alternatively matched
	std::string outfname(outdir);
	outfname.append("/validations.");
	outfname.append(addl_id);
	outfname.append(".alternative.txt");
	
	std::ofstream outfs(outfname);
	check_fstream_for_goodness(outfs, outfname, "save_comparisons_to_tables", true);
	
	outfs
	    << "Sample name" << "\t"
	    << "publication" << "\t"
	    << "validated call ID" << "\t"	    
	    << "call type" << "\t"
	    << "chromosome" << "\t"
	    << "begin" << "\t"
	    << "end" << "\t"
	    << "number unique positions in region" << "\t"
	    
	    << "my call donor begin" << "\t"
	    << "my call donor end" << "\t"	    
	    << "my call receiver begin" << "\t"
	    << "my call receiver end" << "\t"
	    	    
	    << "my connected component ID" << "\t"
	    << "my event UID" << "\t"
	    
	    << "log odds (log [P null/ P deletion])" << "\t"
		<< "log Wilcoxon" << "\n";
	    
	    
	    
	std::string outhistfname(outdir);
	outhistfname.append("/validations.");
	outhistfname.append(addl_id);
	outhistfname.append(".alternative.hist");
	
	std::ofstream outhistfs(outhistfname);
	check_fstream_for_goodness(outhistfs, outhistfname, "save_comparisons_to_tables", true);	    
	    
	    
	
	for (type_map_string_to_string_to_list_Called_diploid_wrapper_ptr::const_iterator it_val_gen = map_genome_to_valname_to_calls__alternatively_matched.begin();
		it_val_gen != map_genome_to_valname_to_calls__alternatively_matched.end();
		++it_val_gen)
	{	    
		for (type_map_string_to_list_Called_diploid_wrapper_ptr::const_iterator it_alternative_to_list = it_val_gen->second.begin();
			it_alternative_to_list != it_val_gen->second.end();
			++it_alternative_to_list)
	    {
		const type_map_string_to_Validated_rearrangement::const_iterator it_alternative = map_genome_to_callID_to_Validated_rearrangement.at(it_val_gen->first).find(it_alternative_to_list->first);
		
		type_list_Called_diploid_wrapper_ptr::const_iterator it_best_matched_mycall = it_alternative_to_list->second.end();
		uint best_hap;
		type_BI__BI best_affected_regions(empty_BI__BI);
		type_haploid_outcome best_GeneConv_outcome;
		
		for (type_list_Called_diploid_wrapper_ptr::const_iterator it_my_alt_match = it_alternative_to_list->second.begin();
			it_my_alt_match != it_alternative_to_list->second.end();
			++it_my_alt_match)
		{
		    const type_map_uint_to_Event::const_iterator the_Event_it = Conn_Comps.at((*it_my_alt_match)->cc_of_event).events.find((*it_my_alt_match)->event_UID);
		    
		    for (uint hap=0; hap<2; ++hap)
		    {
			if (test_if_haploid_outcome_is_GeneConversion(pair_at<type_haploid_outcome>((*it_my_alt_match)->the_diploid_Event_outcome,hap)))
			{
			    const type_BI__BI affected_regions(
					convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
						    pair_at<BOOST_Interval>((*it_my_alt_match)->the_diploid_profile_brkpts,hap),
						    the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
						    pair_at<type_haploid_outcome>((*it_my_alt_match)->the_diploid_Event_outcome,hap)));
			    
			    if (   BOOST_empty(best_affected_regions.first)  or  BOOST_empty(best_affected_regions.second)
				or BOOST_in(it_alternative->second.breakpoints.lower(), affected_regions.first)  or  BOOST_in(it_alternative->second.breakpoints.upper(), affected_regions.first)
				or BOOST_in(it_alternative->second.breakpoints.lower(), affected_regions.second)  or  BOOST_in(it_alternative->second.breakpoints.upper(), affected_regions.second))
			    {
				it_best_matched_mycall = it_my_alt_match;
				best_hap = hap;
				best_affected_regions = affected_regions;
				best_GeneConv_outcome = pair_at<type_haploid_outcome>((*it_my_alt_match)->the_diploid_Event_outcome,hap);
				
				if (it_alternative->second.breakpoints.lower() == affected_regions.first.lower() or it_alternative->second.breakpoints.lower() == affected_regions.first.upper()
				    or  it_alternative->second.breakpoints.lower() == affected_regions.second.lower() or it_alternative->second.breakpoints.lower() == affected_regions.second.upper()
				    or  it_alternative->second.breakpoints.upper() == affected_regions.first.lower() or it_alternative->second.breakpoints.upper() == affected_regions.first.upper()
				    or  it_alternative->second.breakpoints.upper() == affected_regions.second.lower() or it_alternative->second.breakpoints.upper() == affected_regions.second.upper())
				{  break;  }
			    }
			}//same call
		    }//hap
		}//my matches
		
	    
	    
		const BOOST_Interval donor_region = (best_GeneConv_outcome == hap_outcome__GeneConv_ABA)  ?  best_affected_regions.second  :  best_affected_regions.first;
		const BOOST_Interval receiver_region = (best_GeneConv_outcome == hap_outcome__GeneConv_ABA)  ?  best_affected_regions.first  :  best_affected_regions.second;
			
		const real logodds(mpfr::log10((*it_best_matched_mycall)->RD_odds_del));
		const real logWilcoxon(mpfr::log10((*it_best_matched_mycall)->Wilcoxon_signed_rank_sum_test__pvalue__inferred));
		const real logWilcoxon_matched(mpfr::log10((*it_best_matched_mycall)->wilcoxon_matched_info___null___inferred.wilcoxon_pvalue));
		bool passes_threshold;
		
		if (pvalue_for_RD_test_choice == 0)
		{
			passes_threshold = logodds <= logodds_threshold;				
		}
		else if (pvalue_for_RD_test_choice == 1)
		{
			passes_threshold = logWilcoxon <= logodds_threshold;
		}				
		else if (pvalue_for_RD_test_choice == 2)
		{
			passes_threshold = logWilcoxon_matched <= logodds_threshold;
		}			
		
		if (passes_threshold or !use_threshold)				
		{
		    outfs
			<< it_alternative->second.NA_subject << "\t"
			<< it_alternative->second.published_study << "\t"
			<< it_alternative->second.variant_accession_name << "\t"		
			<< convert_haploid_outcome_to_string(static_cast<type_haploid_outcome>(it_alternative->second.outcome)) << "\t"		
			<< it_alternative->second.chromosome << "\t"
			<< it_alternative->second.breakpoints.lower() << "\t"
			<< it_alternative->second.breakpoints.upper() << "\t"
			<< BOOST_width_inclusive(it_alternative->second.breakpoints) - calculate_NON_uniqueness_of_region(it_alternative->second.chromosome, it_alternative->second.breakpoints, nonunique_regions) << "\t"
			
			<< donor_region.lower() << "\t"
			<< donor_region.upper() << "\t"
			<< receiver_region.lower() << "\t"
			<< receiver_region.upper() << "\t"		
					
			<< (*it_best_matched_mycall)->cc_of_event << "\t"
			<< (*it_best_matched_mycall)->event_UID << "\t"
			
			<< logodds << "\t"
			<< logWilcoxon << "\t"
			<< logWilcoxon_matched << "\n";
					
			if (pvalue_for_RD_test_choice == 0 )
				outhistfs << logodds << "\n";
			else if (pvalue_for_RD_test_choice == 1)
				outhistfs << logWilcoxon << "\n";
			else if (pvalue_for_RD_test_choice == 2)
				outhistfs << logWilcoxon_matched << "\n";
		    
		    ++number_below_theshold__alternative;
		}		
		    
	    }//it_alternative_to_list
	}//it_val_gen
	
	outfs << "\n\nnumber of unmatched validations with log-odds below " << logodds_threshold << ":\t" << number_below_theshold__alternative << "\n";	
	    
	outfs.close();
	outhistfs.close();
	    
    }//alternatively matched
    
    
    
    
    
    {//novel_mycalls
	std::string outhistfname(outdir);
	outhistfname.append("/validations.");
	outhistfname.append(addl_id);
	outhistfname.append(".novel.hist");
	
	std::ofstream outhistfs(outhistfname);
	check_fstream_for_goodness(outhistfs, outhistfname, "save_comparisons_to_tables", true);
	
	
	for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = novel_mycalls.begin();
		it_mygen != novel_mycalls.end();
		++it_mygen)
	{
	    for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
		    it_myev != it_mygen->second.end();
		    ++it_myev)
	    {	
			const real logodds(mpfr::log10(it_myev->second.RD_odds_alt));
			const real logWilcoxon(mpfr::log10(it_myev->second.Wilcoxon_signed_rank_sum_test__pvalue__inferred));
			const real logWilcoxon_matched(mpfr::log10(it_myev->second.wilcoxon_matched_info___null___inferred.wilcoxon_pvalue));
			bool passes_threshold;
			
			if (pvalue_for_RD_test_choice == 0)
			{
				passes_threshold = logodds <= logodds_threshold;				
			}
			else if (pvalue_for_RD_test_choice == 1)
			{
				passes_threshold = logWilcoxon <= logodds_threshold;
			}				
			else if (pvalue_for_RD_test_choice == 2)
			{
				passes_threshold = logWilcoxon_matched <= logodds_threshold;
			}					
		
		
			if (passes_threshold or !use_threshold)
			{ 
				if (pvalue_for_RD_test_choice == 0)
					outhistfs << logodds << "\n";
				else if (pvalue_for_RD_test_choice == 1)
					outhistfs << logWilcoxon << "\n";
				else if (pvalue_for_RD_test_choice == 2)
					outhistfs << logWilcoxon_matched << "\n";
			}
		
	    }//ev
	}//gen
	
	outhistfs.close();
	
	
// 	const type_map_string_to_uint_to_Called_diploid_wrapper restricted_novel_mycalls(novel_mycalls, logodds_threshold);
	
    
    }//novel
    
    
    
    
    
    
    {//summary
	std::string outfname(outdir);
	outfname.append("/validations.");
	outfname.append(addl_id);
	outfname.append(".summary.txt");
    
	std::ofstream outfs(outfname);
	check_fstream_for_goodness(outfs, outfname, "save_comparisons_to_tables___imp", true);
	
	const uint total_num = number_below_theshold__match + number_below_theshold__UNMATCH + number_below_theshold__alternative;
	
	outfs
	    << "matched" << "\t"
	    << "completely unmatched" << "\t"
	    << "umatched but called gene conversion" << "\t"
	    << "Total" << "\t"
	    << "sensitivity" << "\n"
	    
	    << number_below_theshold__match << "\t"
	    << number_below_theshold__UNMATCH << "\t"
	    << number_below_theshold__alternative << "\t"
	    << total_num << "\t"
	    << (longdouble)number_below_theshold__match/total_num << "\n";
	    
	outfs.close();
	
    }//summary
    
}//save_comparisons_to_tables___imp















// type_map_diploid_outcome_to_list_real collect_RD_hypothesis_test_p_values_for_FDR
// 						    (const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
// {
//     
//     type_map_diploid_outcome_to_list_real RD_hypothesis_p_value_by_diploid_outcome;
//     
//     
//     for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
// 	    it_mygen != mycalls.end();
// 	    ++it_mygen)
//     {
// 	for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
// 		it_myev != it_mygen->second.end();
// 		++it_myev)
// 	{
// 	    const type_haploid_outcome__haploid_outcome ordered_diploid_outcome(order_a_diploid_outcome(it_myev->second.the_diploid_Event_outcome));
// 	    
// 	    RD_hypothesis_p_value_by_diploid_outcome[ordered_diploid_outcome].push_back(it_myev->second.RD_hypoth_tests_p_value);	    	  
// 	}//myev
//     }//mygen    
//     
//     
//     return RD_hypothesis_p_value_by_diploid_outcome;
//         
// }//collect_RD_hypothesis_test_p_values_for_FDR


















void make_and_save_RD_FDR_distributions
	    (const std::string &outdir,
	    const type_map_diploid_outcome_to_list_real &RD_hypothesis_p_value_by_diploid_outcome)
    
{
    std::string outrootname(outdir);
    outrootname.append("/RD_pvalue_distribution");
    
    
    std::string outfname;
    
    outfname = outrootname;
    outfname.append(".null.txt");
    write_singlecontainer_to_file<type_list_real>(RD_hypothesis_p_value_by_diploid_outcome.at(type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None)), outfname);
    
    outfname = outrootname;
    outfname.append(".nahr.txt");
    {
	std::ofstream outfs(outfname);
	check_fstream_for_goodness(outfs, outfname, "make_and_save_RD_FDR_distributions", true);
	
	for (type_map_diploid_outcome_to_list_real::const_iterator it_dipout = RD_hypothesis_p_value_by_diploid_outcome.begin();
		it_dipout != RD_hypothesis_p_value_by_diploid_outcome.end();
		++it_dipout)
	{
	    if (  (!test_if_haploid_outcome_is_GeneConversion(it_dipout->first.first)  and  !test_if_haploid_outcome_is_GeneConversion(it_dipout->first.second))
		  and (it_dipout->first != type_haploid_outcome__haploid_outcome(hap_outcome__None,hap_outcome__None)))
	    {
		write_singlecontainer_to_file<type_list_real>(it_dipout->second, outfname, &outfs);    
	    }
	}
    }        
    
    
}//make_and_save_RD_FDR_distributions


















type_map_uint_to_list_BI__string  upload_affected_regions_sequences
				    (const type_map_uint_to_CC &Conn_Comps,
				     const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
    
    type_map_uint_to_list_BI  map_chromo_to_affected_BI;    
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
	for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
		it_myev != it_mygen->second.end();
		++it_myev)
	{	    
	    const type_map_uint_to_Event::const_iterator the_Event_it = Conn_Comps.at(it_myev->second.cc_of_event).events.find(it_myev->second.event_UID);
	    
	    for (uint hap=0; hap<2; ++hap)
	    {
		const type_haploid_outcome the_hap_outcome = pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome,hap);
		
		if (the_hap_outcome != hap_outcome__None)
		{		    
		    const type_BI__BI affected_regions(
				convert_profile_breakpoint_to_outcome_meaningful_absolute_intervals(
					pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts,hap),
					the_Event_it->second.compressed_map_LCR_to_profile__for_each_LCR,
					the_hap_outcome));
		    
		    map_chromo_to_affected_BI[the_Event_it->second.chromos[0]].push_back(affected_regions.first);
		    
		    if (test_if_haploid_outcome_is_GeneConversion(the_hap_outcome))
		    {  map_chromo_to_affected_BI[the_Event_it->second.chromos[0]].push_back(affected_regions.second);  }

		}//non-null
	    }//hap
	}//myev
    }//mygen
    
    
    
    map_chromo_to_affected_BI = pad_all_regions(map_chromo_to_affected_BI, 100);
    
    merge_intervals_on_each_chromosome(map_chromo_to_affected_BI);        
    
    
    return  upload_list_of_chromosome_regions_from_file(map_chromo_to_affected_BI);    
    
}//upload_affected_regions_sequences


















void fill_uniqueness_and_undefinedness_of_each_call
			    (const type_map_uint_to_CC &Conn_Comps,
			     type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
    type_map_uint_to_uint map_Event_to_number_of_unique_bases;
    
    const type_map_uint_to_list_BI nonuniq_regions_of_genome(load_Non_unique_regions_of_genome());
    
    for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
	    it_cc != Conn_Comps.end();
	    ++it_cc)
    {
	for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
		it_ev != it_cc->second.events.end();
		++it_ev)
	{
	    map_Event_to_number_of_unique_bases[it_ev->second.UID]
			= BOOST_width_inclusive(it_ev->second.region_between_and_including_the_LCRs_themselves)
				-  calculate_NON_uniqueness_of_region(
						    it_ev->second.chromos[0], 
						    it_ev->second.region_between_and_including_the_LCRs_themselves,
						    nonuniq_regions_of_genome);	
	}//ev
    }//cc
    
    
    const type_set_uint undefined_CCs(upload_connected_components_involving_Undefined_regions_of_genome(data_directory));
    
    
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_mygen = mycalls.begin();
	    it_mygen != mycalls.end();
	    ++it_mygen)
    {
	for (type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_mygen->second.begin();
		it_myev != it_mygen->second.end();
		++it_myev)
	{
	    it_myev->second.number_of_unique_bases_in_region = map_Event_to_number_of_unique_bases.at(it_myev->second.event_UID);
	    
	    it_myev->second.from_a_CC_with_undefined_regions = (undefined_CCs.count(it_myev->second.cc_of_event) > 0);
	    
	}//myev
    }//mygen
    
    
}//fill_uniqueness_of_each_Event


















void save_Affected_Genes
	(const std::string &outdir,
	 const type_map_string_to_Affected_Gene_catalog &affected_Genes_by_individual,
	const type_map_uint_to_BI_to_Gene &all_Genes_by_region)
{    
    
    type_map_string_to_Gene  all_Genes;
    for (type_map_uint_to_BI_to_Gene::const_iterator it_chr = all_Genes_by_region.begin();
	    it_chr != all_Genes_by_region.end();
	    ++it_chr)
    {
	for (type_map_BI_to_Gene::const_iterator it_bi = it_chr->second.begin();
		it_bi != it_chr->second.end();
		++it_bi)
	{
	    all_Genes[it_bi->second.ensembl_id] = it_bi->second;	    	    
	}//bi
    }//chr
    
    
    
    
    typedef std::map<std::string,type_map_haploid_outcome_to_uint>  type_map_string_to_haploid_outcome_to_uint;
    
    type_map_string_to_haploid_outcome_to_uint  map_gene_to_outcome_to_count;
    
    
    //convention:  "hap_outcome__GeneConv_ABA"  <--->  donated itself to elsewhere    
    
    for (type_map_string_to_Affected_Gene_catalog::const_iterator it_sample = affected_Genes_by_individual.begin();
	    it_sample != affected_Genes_by_individual.end();
	    ++it_sample)
    {
	for (type_map_string_to_Gene::const_iterator it_gene = it_sample->second.Genes_affected_by_Inversions.begin();
		it_gene != it_sample->second.Genes_affected_by_Inversions.end();
		++it_gene)
	{  ++map_gene_to_outcome_to_count[it_gene->first][hap_outcome__Inv];  }
	
	for (type_map_string_to_Gene::const_iterator it_gene = it_sample->second.Genes_affected_by_Deletions.begin();
		it_gene != it_sample->second.Genes_affected_by_Deletions.end();
		++it_gene)
	{  ++map_gene_to_outcome_to_count[it_gene->first][hap_outcome__Del];  }
	
	for (type_map_string_to_Gene::const_iterator it_gene = it_sample->second.Genes_affected_by_Duplications.begin();
		it_gene != it_sample->second.Genes_affected_by_Duplications.end();
		++it_gene)
	{  ++map_gene_to_outcome_to_count[it_gene->first][hap_outcome__Dup];  }
	
	for (type_map_string_to_Gene::const_iterator it_gene = it_sample->second.Genes_affected_by_GeneConversion__donated_itself_to_elswhere.begin();
		it_gene != it_sample->second.Genes_affected_by_GeneConversion__donated_itself_to_elswhere.end();
		++it_gene)
	{  ++map_gene_to_outcome_to_count[it_gene->first][hap_outcome__GeneConv_ABA];  }
	
	for (type_map_string_to_Gene::const_iterator it_gene = it_sample->second.Genes_affected_by_GeneConversion__received_from_somewhere_else.begin();
		it_gene != it_sample->second.Genes_affected_by_GeneConversion__received_from_somewhere_else.end();
		++it_gene)
	{  ++map_gene_to_outcome_to_count[it_gene->first][hap_outcome__GeneConv_BAB];  }	
    }//it_sample
    
    
    
    std::string outfname(outdir);
    outfname.append("/affected_gene_tallies.txt");
    
    std::ofstream outfs(outfname);
    check_fstream_for_goodness(outfs, outfname, "save_Affected_Genes", true);
    
    
    outfs
	<< "ensembl ID" << "\t"
	<< "gene name" << "\t"
	
	<< "number inversions" << "\t"
	<< "number deletions" << "\t"
	<< "number duplications" << "\t"
	<< "number donated to elsewhere in Gene Conversion" << "\t"
	<< "number received from elsewhere in Gene Conversion" << "\t"
	
	<< "chromosome" << "\t"
	<< "begin" << "\t"
	<< "end" << "\t"
	<< "description" << "\n";
	
    
    
    for (type_map_string_to_haploid_outcome_to_uint::const_iterator it_gene = map_gene_to_outcome_to_count.begin();
	    it_gene != map_gene_to_outcome_to_count.end();
	    ++it_gene)
    {
	outfs
	    << all_Genes.at(it_gene->first).ensembl_id << "\t"
	    << all_Genes.at(it_gene->first).gene_name << "\t";
	    
	    if (it_gene->second.count(hap_outcome__Inv) > 0)
		outfs << it_gene->second.at(hap_outcome__Inv) << "\t";
	    else
		outfs << "0\t";
	    
	    if (it_gene->second.count(hap_outcome__Del) > 0)
		outfs << it_gene->second.at(hap_outcome__Del) << "\t";
	    else
		outfs << "0\t";
	    
	    if (it_gene->second.count(hap_outcome__Dup) > 0)
		outfs << it_gene->second.at(hap_outcome__Dup) << "\t";
	    else
		outfs << "0\t";
	    
	    if (it_gene->second.count(hap_outcome__GeneConv_ABA) > 0)
		outfs << it_gene->second.at(hap_outcome__GeneConv_ABA) << "\t";
	    else
		outfs << "0\t";
	    
	    if (it_gene->second.count(hap_outcome__GeneConv_BAB) > 0)
		outfs << it_gene->second.at(hap_outcome__GeneConv_BAB) << "\t";
	    else
		outfs << "0\t";
	    
	outfs
	    << all_Genes.at(it_gene->first).chromosome << "\t"
	    << all_Genes.at(it_gene->first).coordinates.lower() << "\t"
	    << all_Genes.at(it_gene->first).coordinates.upper() << "\t"
	    << all_Genes.at(it_gene->first).description << "\n";
			
    }//it_gene

    
    outfs.close();
    
    
}//save_Affected_Genes


    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    

type_map_uint_to_BI_to_Gene  read_Genes_from_file
		    (const std::string &gene_fname)
{
    type_map_uint_to_BI_to_Gene  genes_read_in;
    
    FILE *gfile = std::fopen(gene_fname.c_str(), "r");
    if (check_fstream_for_goodness(gfile, gene_fname, "read_Genes_from_file", false))
    {
	std::fscanf(gfile, "%[^\n]\n");  //skip header 
	char in_whole_line[option_size];
	
// 	while (EOF !=  std::fscanf(gfile, "%[^\t]\t%[^\t]\t%[^\t]\t%u\t%u\t%[^\n]",
// 					in_ensembl, in_name, in_chr, &in_start, &in_end, in_desc))
	while (EOF != std::fscanf(gfile, "%[^\n]", in_whole_line))
	{
		type_vector_string parsed_fields;
		const std::string whole_line_str(in_whole_line);
		boost::split(parsed_fields, whole_line_str, boost::is_any_of("\t"), boost::algorithm::token_compress_on);
// 		assert(parsed_fields.size() >= 5);
		
		std::string description("\"\"");
		if (parsed_fields.size() > 5)
		{  description = parsed_fields.at(5);  }
		
		const int chr_val = convert_chromosome_string_to_number(parsed_fields.at(2));
		const uint begin_coord = std::atoi(parsed_fields.at(3).c_str());
		const uint end_coord = std::atoi(parsed_fields.at(4).c_str());
		
		if (chr_val > 0 and chr_val <= 24)
		{
			genes_read_in[(uint)chr_val][BOOST_Interval(begin_coord, end_coord)]
				= Gene(parsed_fields.at(1), parsed_fields.at(0),
					(uint)chr_val,
					BOOST_Interval(begin_coord, end_coord),
					description);		
		}
		
	    std::fscanf(gfile, "\n");
	}
	
	std::fclose(gfile);
    }    
    
    return  genes_read_in;            

}//read_Genes_from_file

    
    
    
    


    
    
    
    
    
    
    
    
void  remove_calls_in_undefined_regions_of_genome
		(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls)
{    
        
    type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_gen = mycalls.begin();
    while(it_gen != mycalls.end())
    {
		type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_gen->second.begin();
		while(it_myev != it_gen->second.end())
		{
			if (it_myev->second.from_a_CC_with_undefined_regions)
			{
		// 		mycalls__undefined[it_gen->first][it_myev->first] = it_myev->second;
				it_gen->second.erase(it_myev++);
			}
			else
			{
				++it_myev;
			}
		}//it_myev
		
		
		if (it_gen->second.empty())
		{  mycalls.erase(it_gen++);  }
		else
		{  ++it_gen;  }
		
    }//it_gen              
    
}//remove_calls_in_undefined_regions_of_genome








void  remove_calls_without_RD
		(type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{    
        
    type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_gen = mycalls.begin();
    while(it_gen != mycalls.end())
    {
		type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_gen->second.begin();
		while(it_myev != it_gen->second.end())
		{
			if (!it_myev->second.set_RD_odds)
			{
				it_gen->second.erase(it_myev++);
			}
			else
			{
				++it_myev;
			}
		}//it_myev
		
		
		if (it_gen->second.empty())
		{  mycalls.erase(it_gen++);  }
		else
		{  ++it_gen;  }
	
    }//it_gen
             
    
}//remove_calls_without_RD





























void save_histogram_positive_negative_log_odds
			(const std::string &outdir,
			 const std::string &addl_id,
			const int &pvalue_for_RD_test_choice,
			 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
// 	const std::string testid = pvalue_for_RD_test_choice  ?  ".ratio" : ".wilcoxon";	
	
    std::string outfname_pos(outdir);
    outfname_pos.append("/mycalls.");
    outfname_pos.append(addl_id);
// 	outfname_pos.append(testid);	
    outfname_pos.append(".rd.positive.hist");    
    std::ofstream outfs_positive(outfname_pos);
    check_fstream_for_goodness(outfs_positive, outfname_pos, "save_comparisons_to_tables", true);		
    		
    std::string outfname_neg(outdir);
    outfname_neg.append("/mycalls.");
    outfname_neg.append(addl_id);
// 	outfname_neg.append(testid);
    outfname_neg.append(".rd.negative.hist");        
    std::ofstream outfs_negative(outfname_neg);
    check_fstream_for_goodness(outfs_negative, outfname_neg, "save_comparisons_to_tables", true);		


    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
			it_mygen != mycalls.end();
			++it_mygen)
    {
		for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{
			if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None))
			{
				if (pvalue_for_RD_test_choice == 0)
					outfs_negative << mpfr::log10(it_myev->second.RD_odds_best) << "\n";
				else if (pvalue_for_RD_test_choice == 1)
					outfs_negative << (it_myev->second.Wilcoxon_signed_rank_sum_test__pvalue__inferred) << "\n";
				else if (pvalue_for_RD_test_choice == 2)
					outfs_negative << (it_myev->second.wilcoxon_matched_info___null___inferred.wilcoxon_pvalue) << "\n";
// 				std::cerr
// 					<< it_mygen->first 
// 					<< "  " << it_myev->second.cc_of_event 
// 					<< "  " << it_myev->second.event_UID 
// 					<< "  " << mpfr::log10(it_myev->second.RD_odds_best)
// 					<< "\n";
			}
			else if (!test_if_haploid_outcome_is_GeneConversion(it_myev->second.the_diploid_Event_outcome.first) 
				and !test_if_haploid_outcome_is_GeneConversion(it_myev->second.the_diploid_Event_outcome.second))
			{
				if (pvalue_for_RD_test_choice == 0)
					outfs_positive << mpfr::log10(it_myev->second.RD_odds_alt) << "\n";		
				else if (pvalue_for_RD_test_choice == 1)
					outfs_positive << (it_myev->second.Wilcoxon_signed_rank_sum_test__pvalue__inferred) << "\n";
				else if (pvalue_for_RD_test_choice == 2)
					outfs_positive << (it_myev->second.wilcoxon_matched_info___null___inferred.wilcoxon_pvalue) << "\n";				
			}
		}//myev
    }//mygen
    
    
    outfs_positive.close();
    outfs_negative.close();   
    
}//save_histogram_positive_negative_log_odds














type_map_string_to_uint_to_Called_diploid_wrapper  restrict_to_NAHR_using_RD_threshold
				    (const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
				    const real &logodds_threshold,
					const int &pvalue_for_RD_test_choice)
					
{
    type_map_string_to_uint_to_Called_diploid_wrapper  restricted_calls;
    
    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
			it_mygen != mycalls.end();
			++it_mygen)
    {
		for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
				it_myev != it_mygen->second.end();
				++it_myev)
		{
			bool passes_threshold;
			if (pvalue_for_RD_test_choice == 0)
			{
				passes_threshold = mpfr::log10(it_myev->second.RD_odds_alt) <= logodds_threshold;				
			}
			else if (pvalue_for_RD_test_choice == 1)
			{
				passes_threshold = mpfr::log10(it_myev->second.Wilcoxon_signed_rank_sum_test__pvalue__inferred) <= logodds_threshold;
			}
			else if (pvalue_for_RD_test_choice == 2)
			{
				passes_threshold = mpfr::log10(it_myev->second.wilcoxon_matched_info___null___inferred.wilcoxon_pvalue) <= logodds_threshold;
			}			
			
			if (  (test_if_haploid_outcome_is_NAHR(it_myev->second.the_diploid_Event_outcome.first) 
				or test_if_haploid_outcome_is_NAHR(it_myev->second.the_diploid_Event_outcome.second))
				and  passes_threshold)
			{
				restricted_calls[it_mygen->first][it_myev->first] = it_myev->second;
			}
			else if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None))
			{
				restricted_calls[it_mygen->first][it_myev->first] = it_myev->second;
			}		
		}//myev
    }//mygen
    
    
    return restricted_calls;        
    
}//restrict_to_NAHR_using_RD_threshold

















void  remove_calls_on_sex_chromosomes
		(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls,
		const type_map_uint_to_CC &Conn_Comps)
{   
        
    type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_gen = mycalls.begin();
    while(it_gen != mycalls.end())
    {
		type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_gen->second.begin();
		while(it_myev != it_gen->second.end())
		{
			if (Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).chromos[0] >= 23)
			{
				it_gen->second.erase(it_myev++);				
			}//sex chromo
			else
			{
				++it_myev;
			}
		}//it_myev
		
		
		if (it_gen->second.empty())
		{  mycalls.erase(it_gen++);  }
		else
		{  ++it_gen;  }
	
    }//it_gen        
    
}//remove_calls_in_undefined_regions_of_genome












void read_Wilcoxon_pvalues_from_file
			(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls,
			const std::string &postdatadir)
{
	std::string wilcoxon_fname(postdatadir);
	wilcoxon_fname.append("/all_Wilcoxon_signed_rank_sum_test");
	std::ifstream infs(wilcoxon_fname);
	
	std::string in_genome_name;
	uint in_cc;
	uint in_ev;
	real in_wilcoxon_pvalue__inferred;
	real in_wilcoxon_pvalue__unique;
	
	while (infs >> in_genome_name)
	{
		infs >> in_cc >> in_ev >> in_wilcoxon_pvalue__inferred >> in_wilcoxon_pvalue__unique;
// 		std::cerr
// 			<< "line read in: "
// 			<< in_genome_name << "\t" 
// 			<< in_cc << "\t"
// 			<< in_ev << "\t"
// 			<< in_wilcoxon_pvalue__inferred << "\t" 
// 			<< in_wilcoxon_pvalue__unique << "\n";
			
// 		std::cerr
// 			<< "mycalls.count(in_genome_name) = " << mycalls.count(in_genome_name)
// 			<< "\nmycalls.at(in_genome_name).count(in_ev) = " << mycalls.at(in_genome_name).count(in_ev)
// 			<< "\n";
		
		
		mycalls.at(in_genome_name).at(in_ev).Wilcoxon_signed_rank_sum_test__pvalue__inferred = in_wilcoxon_pvalue__inferred;			
		mycalls.at(in_genome_name).at(in_ev).Wilcoxon_signed_rank_sum_test__pvalue__unique = in_wilcoxon_pvalue__unique;		
	}			

	infs.close();		

	
}//read_Wilcoxon_pvalues_from_file












void read_matched_Wilcoxon_pvalues_from_file
			(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls,
			const std::string &postdatadir)
{
	std::string wilcoxon_fname(postdatadir);
	wilcoxon_fname.append("/all_Wilcoxon_signed_rank_sum_test.ratios.matched");
	std::ifstream infs(wilcoxon_fname);
	
	std::string in_genome_name;
	uint in_cc;
	uint in_ev;
	real in_wilcoxon_matched_pvalue__inferred;
	real in_wilcoxon_matched_pvalue__unique;
	
	while (infs >> in_genome_name)
	{
		infs
			//info:
			>> in_cc >> in_ev;
			
// 			std::cerr << "\n\tin_genome_name = [" << in_genome_name << "], in_cc = [" << in_cc << "], in_ev = [" << in_ev << "]\n";
			
		infs
			//p-values:
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___null___inferred.wilcoxon_pvalue
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___null___unique.wilcoxon_pvalue
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___alternative___inferred.wilcoxon_pvalue
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___alternative___unique.wilcoxon_pvalue
			
			//likelihoods:
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___null___inferred.wilcoxon_likelihood
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___null___unique.wilcoxon_likelihood
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___alternative___inferred.wilcoxon_likelihood
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___alternative___unique.wilcoxon_likelihood
			
			//signed_rank_sums:
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___null___inferred.signed_rank_sum_value
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___null___unique.signed_rank_sum_value
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___alternative___inferred.signed_rank_sum_value
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_matched_info___alternative___unique.signed_rank_sum_value
			
			//Likelihood ratios:
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_likelihood_ratio___inferred
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_likelihood_ratio___unique
			
			//ratios:
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_statistic_ratio___inferred
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_statistic_ratio___unique
			
			//Cauchy CDF (standard Cauchy)
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_Cauchy_cdf___inferred
			>> mycalls.at(in_genome_name).at(in_ev).wilcoxon_Cauchy_cdf___unique;
	}			

	infs.close();		

	
}//read_matched_Wilcoxon_pvalues_from_file











void identify_wilcoxon_bad_calls
		(type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		const type_map_uint_to_CC &Conn_Comps)		
{
	
	type_map_uint_to_set_string map_problem_event_to_genomes__null;
	type_map_uint_to_set_string map_problem_event_to_genomes__positive;
	
	const uint min_length_req = 3500;
	
	std::cerr << "\n\n\nproblematic wilcoxon NULL calls:\n\n";
	
	type_map_uint_to_uint  map_event_to_width;
	
	
	
	for (type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_mygen = mycalls.begin();
			it_mygen != mycalls.end();
			++it_mygen)
	{
		for (type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{	
			
			map_event_to_width[it_myev->first] = BOOST_width_inclusive(Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).region_between_and_including_the_LCRs_themselves);
			
			if (test_if_is_diploid_No_Event(it_myev->second.the_diploid_Event_outcome))
			{
				const bool long_enough = BOOST_width_inclusive(Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).region_between_and_including_the_LCRs_themselves)  >  min_length_req;
				
				if (it_myev->second.wilcoxon_Lyapunov_CLT___inferred < -2  and  long_enough)
				{
					std::cerr	
						<< "\tgenome:" << it_myev->second.associated_genome << "\n"
						<< "\tcc:" << it_myev->second.cc_of_event << "\n"
						<< "\tev:" << it_myev->second.event_UID << "\n"
						<< "\tinferred:" << it_myev->second.wilcoxon_Lyapunov_CLT___inferred << "\n"
						<< "\tunique:" << it_myev->second.wilcoxon_Lyapunov_CLT___unique << "\n\n";
						
					std::string  genome_and_call(it_myev->second.associated_genome);
					genome_and_call.append("\t");
					genome_and_call.append(convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first));
					genome_and_call.append(", ");
					genome_and_call.append(convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second));
						
					map_problem_event_to_genomes__null[it_myev->second.event_UID].insert(genome_and_call);
				}
			}	
		}//mycall
	}//mygen	
	
// 	print_map_to_set<uint, std::string>(map_problem_event_to_genomes__null, "map_problem_event_to_genomes__null", NULL, true);	
	{
		std::stringstream outss;
		outss << "\n\nmap_problem_event_to_genomes__null:\n";
		
		for (type_map_uint_to_set_string::const_iterator it_map = map_problem_event_to_genomes__null.begin();
			it_map != map_problem_event_to_genomes__null.end();
			++it_map)
		{
			outss << "\tEvent:  " << it_map->first << ",      width:  " << map_event_to_width.at(it_map->first) << "\n";
			
			for (type_set_string::const_iterator it_gen = it_map->second.begin();
				it_gen != it_map->second.end();
				++it_gen)
			{
				outss << "\t" << *it_gen << "\n";
			}
			
			outss << "\n";
		}
		
		std::cerr << "\n" << outss.str() << "\n";
	}
		
		
	
	
	std::cerr << "\n\n\nproblematic wilcoxon positive calls:\n\n";
	
	
	uint num_good_hap_dels = 0;
	uint num_good_dip_dels = 0;
	uint num_bad_hap_dels = 0;
	uint num_bad_dip_dels = 0;	
	
	uint num_good_hap_dups = 0;
	uint num_good_dip_dups = 0;
	uint num_bad_hap_dups = 0;
	uint num_bad_dip_dups = 0;		
	
	for (type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_mygen = mycalls.begin();
			it_mygen != mycalls.end();
			++it_mygen)
	{
		for (type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{	
			const bool long_enough = BOOST_width_inclusive(Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).region_between_and_including_the_LCRs_themselves)  >  min_length_req;			
			
			if (test_if_at_least_one_haploid_is_NAHR(it_myev->second.the_diploid_Event_outcome)
				and !test_if_at_least_one_haploid_is_Gene_Conversion(it_myev->second.the_diploid_Event_outcome)
				and  long_enough)
			{

				
				if (it_myev->second.wilcoxon_Lyapunov_CLT___inferred > -1)
				{
					std::cerr	
						<< "\tgenome:" << it_myev->second.associated_genome << "\n"
						<< "\tcc:" << it_myev->second.cc_of_event << "\n"
						<< "\tev:" << it_myev->second.event_UID << "\n"
						<< "\tcall: " << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first)
							<< ", " << convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second) << "\n"
						<< "\tinferred:" << it_myev->second.wilcoxon_Lyapunov_CLT___inferred << "\n"
						<< "\tunique:" << it_myev->second.wilcoxon_Lyapunov_CLT___unique << "\n\n";
						
					std::string  genome_and_call(it_myev->second.associated_genome);
					genome_and_call.append("\t");
					genome_and_call.append(convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.first));
					genome_and_call.append(", ");
					genome_and_call.append(convert_haploid_outcome_to_string(it_myev->second.the_diploid_Event_outcome.second));
						
					map_problem_event_to_genomes__positive[it_myev->second.event_UID].insert(genome_and_call);		
					
					if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
						++num_bad_dip_dels;
					else if (test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
						++num_bad_hap_dels;
					else if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
						++num_bad_dip_dups;
					else if (test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
						++num_bad_hap_dups;					
					
				}
				else //good positive
				{
					if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
						++num_good_dip_dels;
					else if (test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Del))
						++num_good_hap_dels;
					else if (test_if_is_diploid_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
						++num_good_dip_dups;
					else if (test_if_exactly_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome, hap_outcome__Dup))
						++num_good_hap_dups;	
				}
			}						
		}//mycall
	}//mygen
	
// 	print_map_to_set<uint, std::string>(map_problem_event_to_genomes__positive, "map_problem_event_to_genomes__positive", NULL, true);
	{
		std::stringstream outss;
		outss << "\n\nmap_problem_event_to_genomes__positive:\n";
		
		for (type_map_uint_to_set_string::const_iterator it_map = map_problem_event_to_genomes__positive.begin();
			it_map != map_problem_event_to_genomes__positive.end();
			++it_map)
		{
			outss << "\tEvent:  " << it_map->first << ",      width:  " << map_event_to_width.at(it_map->first) << "\n";
			
			for (type_set_string::const_iterator it_gen = it_map->second.begin();
				it_gen != it_map->second.end();
				++it_gen)
			{
				outss << "\t" << *it_gen << "\n";
			}
			
			outss << "\n";
		}
		
		std::cerr << "\n" << outss.str() << "\n";
	}
	
	std::cerr << "\n"
		<< "num_good_hap_dels = " << num_good_hap_dels << "\n"
		<< "num_good_dip_dels = " << num_good_dip_dels << "\n"
		<< "num_bad_hap_dels = " << num_bad_hap_dels << "\n"
		<< "num_bad_dip_dels = " << num_bad_dip_dels << "\n"
		
		<< "num_good_hap_dups = " << num_good_hap_dups << "\n"
		<< "num_good_dip_dups = " << num_good_dip_dups << "\n"
		<< "num_bad_hap_dups = " << num_bad_hap_dups << "\n"
		<< "num_bad_dip_dups = " << num_bad_dip_dups << "\n";
	
	
}//identify_wilcoxon_bad_calls
















void identify_wilcoxon_bogus
		(type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
	
    for (type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_mygen = mycalls.begin();
			it_mygen != mycalls.end();
			++it_mygen)
    {
		for (type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{
			if (it_myev->second.wilcoxon_matched_info___null___inferred.wilcoxon_pvalue == -1
				or  it_myev->second.wilcoxon_matched_info___alternative___inferred.wilcoxon_pvalue == -1)
			{
				it_myev->second.wilcoxon_likelihood_ratio___inferred = -1;
				it_myev->second.wilcoxon_statistic_ratio___inferred = -1;
				it_myev->second.wilcoxon_Cauchy_cdf___inferred = -1;
			}
			
			if (it_myev->second.wilcoxon_matched_info___null___unique.wilcoxon_pvalue == -1
				or  it_myev->second.wilcoxon_matched_info___alternative___unique.wilcoxon_pvalue == -1)
			{
				it_myev->second.wilcoxon_likelihood_ratio___unique = -1;
				it_myev->second.wilcoxon_statistic_ratio___unique = -1;				
				it_myev->second.wilcoxon_Cauchy_cdf___unique = -1;
			}			
		}//mycall
	}//mygen
	
}//identify_wilcoxon_bogus










void save_raw_count_ratio_table
			(const std::string &outdir,
			 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
			const type_map_uint_to_CC &Conn_Comps,
			const std::string &gender_and_population_filename)			
{
	std::cerr << "\n\t\tmycalls.size = " << mycalls.size() << "\n";
	
	std::string outfname_combined(outdir);
	outfname_combined.append("/raw_count_ratios.combined.hist");
	std::ofstream outfs_combined(outfname_combined);
	check_fstream_for_goodness(outfs_combined, outfname_combined, "save_comparisons_to_tables", true);		 
		
	std::string outfname_pos(outdir);
	outfname_pos.append("/raw_count_ratios.positive.hist");
	std::ofstream outfs_pos(outfname_pos);
	check_fstream_for_goodness(outfs_pos, outfname_pos, "save_comparisons_to_tables", true);	

	std::string outfname_neg(outdir);
	outfname_neg.append("/raw_count_ratios.negative.hist");
	std::ofstream outfs_neg(outfname_neg);
	check_fstream_for_goodness(outfs_neg, outfname_neg, "save_comparisons_to_tables", true);			
	
	std::cerr << "loading gedner,,,,\n";
	const type_map_string_to_string__bool map_genome_to_pop_and_female(load_gender_and_population_map(gender_and_population_filename));
	std::cerr << "printing gedner,,,,\n";
// 	print_map_keys_and_pair_values<std::string, std::string, bool>(map_genome_to_pop_and_female, "map_genome_to_pop_and_female");
	std::cerr << "looping,,,,\n";
	
	std::cerr << "\tmycalls.size = " << mycalls.size() << "\n";
	for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
				it_mygen != mycalls.end();
				++it_mygen)
	{
		std::cerr << "\t\tgen.size = " << it_mygen->second.size() << "\n";
		for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{
			// 1. inferred
			// 2. unique
			// 3. cc
			// 4. ev
			// 5. call 0
			// 6. call 1
			// 7. genome
// 			std::cerr << "\ttrying:  " << it_myev->second.associated_genome << ":  [";
// 			std::cerr << map_genome_to_pop_and_female.at(it_myev->second.associated_genome).first << "]\n";
// 			
			uint better_brkpt_indx = 0;
			real min_odds_brk(mpfr::const_infinity());
			
			for (uint hap=0; hap<2; ++hap)
			{				
				if (pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome, hap) == hap_outcome__Del
					or  pair_at<type_haploid_outcome>(it_myev->second.the_diploid_Event_outcome, hap) == hap_outcome__Dup)
				{
					if (pair_at<real>(it_myev->second.brkpt_odds, hap) != -1  and  pair_at<real>(it_myev->second.brkpt_odds, hap) < min_odds_brk)
					{
						min_odds_brk = pair_at<real>(it_myev->second.brkpt_odds, hap);
						better_brkpt_indx = (int)hap;
					}
				}
			}//hap			
			
			if (min_odds_brk == mpfr::const_infinity())
			{  min_odds_brk = -1;  }
			
// 			if (it_myev->second.brkpt_odds.first == -1)
// 			{  
// 				min_odds_brk = it_myev->second.brkpt_odds.second;  
// 				better_brkpt_indx = 1;
// 			}
// 			else if (it_myev->second.brkpt_odds.second == -1)
// 			{  
// 				min_odds_brk = it_myev->second.brkpt_odds.first;  
// 				better_brkpt_indx = 0;
// 			}
// 			else 
// 			{
// 				min_odds_brk = mpfr::min(it_myev->second.brkpt_odds.first, it_myev->second.brkpt_odds.second); 
// 				if (it_myev->second.brkpt_odds.first < it_myev->second.brkpt_odds.second)
// 				{  better_brkpt_indx = 0;  }
// 				else
// 				{  better_brkpt_indx = 1;  }
// 			}
			
			BOOST_Interval brkpts(empty_BI);
			if (min_odds_brk != -1)
			{
				brkpts = convert_profile_breakpoint_to_absolute_breakpoints(pair_at<BOOST_Interval>(it_myev->second.the_diploid_profile_brkpts, better_brkpt_indx).lower(), Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).compressed_map_LCR_to_profile__for_each_LCR);   
			}
			
			outfs_combined 
				<< it_myev->second.raw_count_ratio__inferred << "\t"
				<< it_myev->second.raw_count_ratio__unique << "\t"
				<< it_myev->second.cc_of_event << "\t"
				<< it_myev->second.event_UID << "\t"
				<< min_odds_brk << "\t"
				<< Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).chromos[0] << "\t"
				<< brkpts.lower() << "\t"
				<< brkpts.upper() << "\t"
				
				<< (uint)it_myev->second.the_diploid_Event_outcome.first << "\t"
				<< (uint)it_myev->second.the_diploid_Event_outcome.second << "\t"
				<< BOOST_width_inclusive(Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).region_between_and_including_the_LCRs_themselves) << "\t"
				<< it_myev->second.associated_genome << "\t"
//				<< map_genome_to_pop_and_female.at(it_myev->second.associated_genome).first << "\n";  COMMENDTED B Y MATT SIMULATION
				<< "pop\n";
				
				
			
			if (!test_if_at_least_one_haploid_is_del_or_dup(it_myev->second.the_diploid_Event_outcome))
			{
				outfs_neg 
					<< it_myev->second.raw_count_ratio__inferred << "\t"
					<< it_myev->second.raw_count_ratio__unique << "\t"
					<< it_myev->second.cc_of_event << "\t"
					<< it_myev->second.event_UID << "\t"
					<< min_odds_brk << "\t"
					<< brkpts.lower() << "\t"
					<< brkpts.upper() << "\t"					
					
					<< (uint)it_myev->second.the_diploid_Event_outcome.first << "\t"
					<< (uint)it_myev->second.the_diploid_Event_outcome.second << "\t"
					<< BOOST_width_inclusive(Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).region_between_and_including_the_LCRs_themselves) << "\t"
					<< it_myev->second.associated_genome << "\t"
//				<< map_genome_to_pop_and_female.at(it_myev->second.associated_genome).first << "\n";  COMMENDTED B Y MATT SIMULATION
				<< "pop\n";
			}
			else
			{
				outfs_pos 
					<< it_myev->second.raw_count_ratio__inferred << "\t"
					<< it_myev->second.raw_count_ratio__unique << "\t"
					<< it_myev->second.cc_of_event << "\t"
					<< it_myev->second.event_UID << "\t"
					<< min_odds_brk << "\t"
					<< brkpts.lower() << "\t"
					<< brkpts.upper() << "\t"					
					
					<< (uint)it_myev->second.the_diploid_Event_outcome.first << "\t"
					<< (uint)it_myev->second.the_diploid_Event_outcome.second << "\t"
					<< BOOST_width_inclusive(Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).region_between_and_including_the_LCRs_themselves) << "\t"
					<< it_myev->second.associated_genome << "\t"
//				<< map_genome_to_pop_and_female.at(it_myev->second.associated_genome).first << "\n";  COMMENDTED B Y MATT SIMULATION
				<< "pop\n";
			}
		}//myev
	}//mygen
	
	
	outfs_combined.close();
	outfs_pos.close();   
	outfs_neg.close();   
	
	std::cerr << "\nDONE with \"save_raw_count_ratio_table\"!\n\n";
    
}//save_raw_count_ratio_table












void debug_low_pvalues
			(const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
	
	typedef std::map<uint, std::set<std::string> >  type_map_uint_to_set_string;
	type_map_uint_to_set_string low_pvalue_NULL_events;
	
    for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
			it_mygen != mycalls.end();
			++it_mygen)
    {		
		std::cerr << it_mygen->first << "\n";
		
		for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{
			if (it_myev->second.the_diploid_Event_outcome == type_haploid_outcome__haploid_outcome(hap_outcome__None, hap_outcome__None))
			{
				if (it_myev->second.wilcoxon_matched_info___null___unique.wilcoxon_pvalue > -1
					and it_myev->second.wilcoxon_matched_info___null___unique.wilcoxon_pvalue < 0.0001)
				{
					std::cerr << "\t" << it_myev->first << "\t" << it_myev->second.wilcoxon_matched_info___null___unique.wilcoxon_pvalue << "\n";	
					low_pvalue_NULL_events[it_myev->first].insert(it_mygen->first);
				}
			}//no event

		}//myev
    }//mygen	
	
	
	
	std::cerr << "\n\nlow p-value events:\n";	
	for (type_map_uint_to_set_string::const_iterator it_lowev = low_pvalue_NULL_events.begin();
		it_lowev != low_pvalue_NULL_events.end();
		++it_lowev)
	{
		std::cerr << "\t" << it_lowev->first << "\t";
		for (type_set_string::const_iterator it_gen = it_lowev->second.begin();
			it_gen != it_lowev->second.end();
			++it_gen)
		{
			std::cerr << *it_gen << ",  ";
		}
		
		std::cerr << "\n";
	}
	
	
}//debug_low_pvalues
































void read_all_raw_count_ratios_from_file
			(type_map_string_to_uint_to_Called_diploid_wrapper  &mycalls,
			const std::string &postdatadir)
{
	std::string wilcoxon_fname(postdatadir);
	wilcoxon_fname.append("/all_raw_count_ratios");
	std::ifstream infs(wilcoxon_fname);
	
	std::string in_genome_name;
	uint in_cc;
	uint in_ev;
	real in_wilcoxon_matched_pvalue__inferred;
	real in_wilcoxon_matched_pvalue__unique;
	
	while (infs >> in_genome_name)
	{
		infs	>> in_cc >> in_ev;
			
		infs	>> mycalls.at(in_genome_name).at(in_ev).raw_count_ratio__inferred
			>> mycalls.at(in_genome_name).at(in_ev).raw_count_ratio__unique;
	}			

	infs.close();		

	
}//read_all_raw_count_ratios_from_file









void map_genes_to_events_and_save_to_file
		(const type_map_uint_to_BI_to_Gene &Genes_in_the_human_genome,
		 const type_map_uint_to_CC &Conn_Comps,
		const std::string &filename)
{	
	std::cerr << "Genes_in_the_human_genome.size() = " << Genes_in_the_human_genome.size() << "\n";
	
	std::string outfname("/users/mmparks/data/mmparks/");
	outfname.append(filename);
	
	std::ofstream outfs(outfname);
	
	if (check_fstream_for_goodness(outfs, outfname, "map_genes_to_events_and_save_to_file", true))
	{		
		for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
			it_cc != Conn_Comps.end();
			++it_cc)
		{
			for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
				it_ev != it_cc->second.events.end();
				++it_ev)
			{
				const type_map_string_to_Gene potentially_affected_Genes(identify_Genes_potentially_affected_by_an_Event(it_ev->second, Genes_in_the_human_genome));
				
				for (type_map_string_to_Gene::const_iterator it_gene = potentially_affected_Genes.begin();
					it_gene != potentially_affected_Genes.end();
					++it_gene)
				{
					outfs
						<< it_ev->second.connected_component << "\t"
						<<  it_ev->second.UID << "\t" 
						<< it_ev->second.chromos[0] << "\t"
						<< it_ev->second.LCRs[0].lower() << "\t"
						<< it_ev->second.LCRs[0].upper() << "\t"
						<< it_ev->second.LCRs[1].lower() << "\t"
						<< it_ev->second.LCRs[1].upper() << "\t"
						<< it_gene->second.coordinates.lower() << "\t"
						<< it_gene->second.coordinates.upper() << "\t"
						<< it_gene->second.ensembl_id << "\t"
						<< it_gene->second.gene_name << "\t\""
						<< it_gene->second.description << "\"\n";
				}//gene
				
		
			}//ev
		}//cc
		
		outfs.close();
	}// good ofstream
	
}//map_genes_to_events_and_save_to_file







void map_genes_to_events_and_save_to_file
		(const type_map_uint_to_BI_to_Gene &Genes_in_the_human_genome,
		 const type_map_uint_to_BI_to_Gene &cancer_genes,
		 const type_map_uint_to_CC &Conn_Comps)
{	
	std::cerr << "Genes_in_the_human_genome.size() = " << Genes_in_the_human_genome.size() << "\n";
	std::cerr << "cancer_genes.size() = " << cancer_genes.size() << "\n";
	
	std::ofstream outfs("/users/mmparks/data/mmparks/out_preprocessing_0.25Mb/potentially_affected_Genes_by_Event");
	std::ofstream outfs_cancer("/users/mmparks/data/mmparks/out_preprocessing_0.25Mb/potentially_affected_cancer_Genes_by_Event");
	
	if (check_fstream_for_goodness(outfs, "/users/mmparks/data/mmparks/out_preprocessing_0.25Mb/potentially_affected_Genes_by_Event",
				   "map_genes_to_events_and_save_to_file", false))
	{
		check_fstream_for_goodness(outfs_cancer, "/users/mmparks/data/mmparks/out_preprocessing_0.25Mb/potentially_affected_cancer_Genes_by_Event",
				   "map_genes_to_events_and_save_to_file", false);
		
		for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
			it_cc != Conn_Comps.end();
			++it_cc)
		{
			for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
				it_ev != it_cc->second.events.end();
				++it_ev)
			{
				const type_map_string_to_Gene potentially_affected_Genes(identify_Genes_potentially_affected_by_an_Event(it_ev->second, Genes_in_the_human_genome));
				
				for (type_map_string_to_Gene::const_iterator it_gene = potentially_affected_Genes.begin();
					it_gene != potentially_affected_Genes.end();
					++it_gene)
				{
					outfs
						<< it_ev->second.connected_component << "\t"
						<<  it_ev->second.UID << "\t" 
						<< it_ev->second.chromos[0] << "\t"
						<< it_ev->second.LCRs[0].lower() << "\t"
						<< it_ev->second.LCRs[0].upper() << "\t"
						<< it_ev->second.LCRs[1].lower() << "\t"
						<< it_ev->second.LCRs[1].upper() << "\t"
						<< it_gene->second.coordinates.lower() << "\t"
						<< it_gene->second.coordinates.upper() << "\t"
						<< it_gene->second.ensembl_id << "\t"
						<< it_gene->second.gene_name << "\t\""
						<< it_gene->second.description << "\"\n";
				}//gene
				
				
				const type_map_string_to_Gene potentially_affected_CANCER_Genes(identify_Genes_potentially_affected_by_an_Event(it_ev->second, cancer_genes));
				
				for (type_map_string_to_Gene::const_iterator it_gene = potentially_affected_CANCER_Genes.begin();
					it_gene != potentially_affected_CANCER_Genes.end();
					++it_gene)
				{
					outfs_cancer 
						<< it_ev->second.connected_component << "\t"
						<<  it_ev->second.UID << "\t" 
						<< it_ev->second.chromos[0] << "\t"
						<< it_ev->second.LCRs[0].lower() << "\t"
						<< it_ev->second.LCRs[0].upper() << "\t"
						<< it_ev->second.LCRs[1].lower() << "\t"
						<< it_ev->second.LCRs[1].upper() << "\t"
						<< it_gene->second.coordinates.lower() << "\t"
						<< it_gene->second.coordinates.upper() << "\t"
						<< it_gene->second.ensembl_id << "\t"
						<< it_gene->second.gene_name << "\t\""
						<< it_gene->second.description << "\"\n";
				}//gene				
				
				

// 				if (!potentially_affected_Genes.empty())
// 				{
// 					const type_set_string affected_set(extract_keys_of_map_and_return_as_set<std::string, Gene>(potentially_affected_Genes));
// 					const std::string affected_str(concatenate_singlecontainer_elements_into_string<type_set_string>(affected_set));
// 					
// 					outfs <<  it_ev->second.connected_component << "\t" <<  it_ev->second.UID << "\t" <<  affected_str << "\n";
// 				}//not empty
			}//ev
		}//cc
		
		outfs.close();
		outfs_cancer.close();
	}// good ofstream
	
}//map_genes_to_events_and_save_to_file
























type_map_uint_to_BI_to_Gene  read_cancer_sites_from_file
		    (const std::string &gene_fname)
{
    type_map_uint_to_BI_to_Gene  genes_read_in;
    
    FILE *gfile = std::fopen(gene_fname.c_str(), "r");
        
    if (check_fstream_for_goodness(gfile, gene_fname, "read_Genes_from_file", false))
    {
	std::fscanf(gfile, "%[^\n]\n");  //skip header
	char in_whole_line[option_size];
	
// 	while (EOF !=  std::fscanf(gfile, "%[^\t]\t%[^\t]\t%[^\t]\t%u\t%u\t%[^\n]",
// 					in_ensembl, in_name, in_chr, &in_start, &in_end, in_desc))
	while (EOF != std::fscanf(gfile, "%[^\n]", in_whole_line))
	{
		type_vector_string parsed_fields;
		const std::string whole_line_str(in_whole_line);
		boost::split(parsed_fields, whole_line_str, boost::is_any_of("\t"), boost::algorithm::token_compress_on);
		
		if (parsed_fields.size() == 4)
		{
			const std::string description(parsed_fields.at(3));
			
			char in_chr[option_size];
			uint in_low;
			uint in_upp;
			std::sscanf(parsed_fields.at(2).c_str(), "%[^:]:%u-%u", in_chr, &in_low, &in_upp);
			
			const int chr_val = convert_chromosome_string_to_number(in_chr);
			
			if (chr_val > 0 and chr_val <= 24)
			{
				genes_read_in[(uint)chr_val][BOOST_Interval(in_low, in_upp)]
					= Gene(parsed_fields.at(1), parsed_fields.at(0),
						(uint)chr_val,
						BOOST_Interval(in_low, in_upp),
						description);		
			}
		}
		
		std::fscanf(gfile, "\n");
	}
	
	std::fclose(gfile);
    }    
    
    return  genes_read_in;            

}//read_cancer_sites_from_file

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
type_map_string_to_uint_to_Called_diploid_wrapper load_restricted_calls_from_matlab
		(const std::string &restricted_calls_fname)
{
	type_map_string_to_uint_to_Called_diploid_wrapper  in_calls;
	
	std::ifstream infs(restricted_calls_fname);
// 	check_fstream_for_goodness(infs, restricted_calls_fname, "load_restricted_calls_from_matlab", true);
	
	std::string  in_genome;
	while (infs >> in_genome)
	{
		Called_diploid_wrapper incall;
		incall.associated_genome = in_genome;
				
		infs >> incall.cc_of_event >> incall.event_UID >> incall.logfdr_inferred;
						
		in_calls[in_genome][incall.event_UID] = incall;
	}//while
	

	return in_calls;	
	
}//load_restricted_calls_from_matlab



























void identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs
		(const type_map_uint_to_BI_to_Gene &Genes_in_the_human_genome,
		 const type_map_uint_to_CC &Conn_Comps)
{	
	std::cerr << "identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs...\n";
	
	const type_map_uint_to_BI_to_Gene psi_DR(load_psi_DR_pseudogene_database("/users/mmparks/data/mmparks/psiDR.v0.txt"));		
	
	std::cerr << "Genes_in_the_human_genome.size() = " << Genes_in_the_human_genome.size() << "\n";
	std::cerr << "psi_DR.size() = " << psi_DR.size() << "\n";			
	
	
	char fname[option_size];
	
	std::sprintf(fname, "%s/pseudo_genes_on_LCRs.txt", scratch_job_output_dir.c_str());
	std::ofstream outfs_psiDR_LCR(fname);
	check_fstream_for_goodness(outfs_psiDR_LCR, fname, "identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs", true);
	
	std::sprintf(fname, "%s/Ensembl_genes_on_LCRs.txt", scratch_job_output_dir.c_str());
	std::ofstream outfs_Ensembl_LCR(fname);
	check_fstream_for_goodness(outfs_Ensembl_LCR, fname, "identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs", true);
	
	std::sprintf(fname, "%s/pseudo_genes_inbetween_LCRs.txt", scratch_job_output_dir.c_str());
	std::ofstream outfs_psiDR_inbetween(fname);
	check_fstream_for_goodness(outfs_psiDR_inbetween, fname, "identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs", true);
	
	std::sprintf(fname, "%s/Ensembl_genes_inbetween_LCRs.txt", scratch_job_output_dir.c_str());
	std::ofstream outfs_Ensembl_inbetween(fname);
	check_fstream_for_goodness(outfs_Ensembl_inbetween, fname, "identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs", true);
	
	
	
	for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{
		for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
			it_ev != it_cc->second.events.end();
			++it_ev)
		{
			{//Ensembl genes
				const type_map_string_to_Gene potentially_affected_Genes(identify_Genes_potentially_affected_by_an_Event(it_ev->second, Genes_in_the_human_genome));
				
				for (type_map_string_to_Gene::const_iterator it_gene = potentially_affected_Genes.begin();
					it_gene != potentially_affected_Genes.end();
					++it_gene)
				{
					if (BOOST_subset(it_gene->second.coordinates, it_ev->second.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs))
					{
						outfs_Ensembl_inbetween 
							<< it_ev->second.connected_component << "\t"
							<< it_ev->second.UID << "\t"
							<< it_ev->second.chromos[0] << "\t"
							<< it_ev->second.LCRs[0].lower() << "\t"
							<< it_ev->second.LCRs[0].upper() << "\t"
							<< it_ev->second.LCRs[1].lower() << "\t"
							<< it_ev->second.LCRs[1].upper() << "\t"
							<< it_gene->second.coordinates.lower() << "\t"
							<< it_gene->second.coordinates.upper() << "\t"							
							<< it_gene->second.gene_name << "\t"
							<< it_gene->second.ensembl_id << "\n";
					}				
					else if (BOOST_overlap(it_gene->second.coordinates, it_ev->second.region_between_and_including_the_LCRs_themselves))
					{
						outfs_Ensembl_LCR 
							<< it_ev->second.connected_component << "\t"
							<< it_ev->second.UID << "\t"
							<< it_ev->second.chromos[0] << "\t"
							<< it_ev->second.LCRs[0].lower() << "\t"
							<< it_ev->second.LCRs[0].upper() << "\t"
							<< it_ev->second.LCRs[1].lower() << "\t"
							<< it_ev->second.LCRs[1].upper() << "\t"
							<< it_gene->second.coordinates.lower() << "\t"
							<< it_gene->second.coordinates.upper() << "\t"								
							<< it_gene->second.gene_name << "\t"
							<< it_gene->second.ensembl_id << "\n";
					}
				}//gene
			}//Ensembl genes
			
			
			{//pseudogenes
				const type_map_string_to_Gene potentially_affected_Genes(identify_Genes_potentially_affected_by_an_Event(it_ev->second, psi_DR));
				
				for (type_map_string_to_Gene::const_iterator it_gene = potentially_affected_Genes.begin();
					it_gene != potentially_affected_Genes.end();
					++it_gene)
				{
					if (BOOST_subset(it_gene->second.coordinates, it_ev->second.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs))
					{
						outfs_psiDR_inbetween 
							<< it_ev->second.connected_component << "\t"
							<< it_ev->second.UID << "\t"
							<< it_ev->second.chromos[0] << "\t"
							<< it_ev->second.LCRs[0].lower() << "\t"
							<< it_ev->second.LCRs[0].upper() << "\t"
							<< it_ev->second.LCRs[1].lower() << "\t"
							<< it_ev->second.LCRs[1].upper() << "\t"
							<< it_gene->second.coordinates.lower() << "\t"
							<< it_gene->second.coordinates.upper() << "\t"								
							<< it_gene->second.ensembl_id << "\t"
							<< it_gene->second.gene_name << "\n";
					}				
					else if (BOOST_overlap(it_gene->second.coordinates, it_ev->second.region_between_and_including_the_LCRs_themselves))
					{
						outfs_psiDR_LCR 
							<< it_ev->second.connected_component << "\t"
							<< it_ev->second.UID << "\t"
							<< it_ev->second.chromos[0] << "\t"
							<< it_ev->second.LCRs[0].lower() << "\t"
							<< it_ev->second.LCRs[0].upper() << "\t"
							<< it_ev->second.LCRs[1].lower() << "\t"
							<< it_ev->second.LCRs[1].upper() << "\t"
							<< it_gene->second.coordinates.lower() << "\t"
							<< it_gene->second.coordinates.upper() << "\t"								
							<< it_gene->second.ensembl_id << "\t"
							<< it_gene->second.gene_name << "\n";
					}
				}//gene
			}//pseudogenes
		}//ev
	}//cc
	
	
	outfs_Ensembl_inbetween.close();
	outfs_Ensembl_LCR.close();
	outfs_psiDR_inbetween.close();
	outfs_psiDR_LCR.close();	
	
	std::cerr << "DONE with identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs!\n";
	
}//identify_Events_whose_LCRs_are_gene_pseudo_gene_pairs












type_map_uint_to_BI_to_Gene load_psi_DR_pseudogene_database
				(const std::string &psidr_fname)
{
	//assume its weird header has been removed.
	
	std::cerr << "\"load_psi_DR_pseudogene_database\"...\n";
	
// 	type_map_string_to_Gene psiDR_DB;
	type_map_uint_to_BI_to_Gene psiDR_DB_map;
		
	char in_pseudo_id[option_size];
	char in_chromo[option_size];
	uint in_begin;
	uint in_end;
	char in_parent_ensembl[option_size];
	char in_parent_name[option_size];
	
	FILE *file_psidr = std::fopen(psidr_fname.c_str(), "r");
	
	check_fstream_for_goodness(file_psidr, psidr_fname, "load_psi_DR_pseudogene_database", true);
	
	//skip header
	std::fscanf(file_psidr, "%*[^\n]\n");
	
	while (EOF != std::fscanf(file_psidr, "%[^\t]\t%*[^\t]\tchr%[^\t]\t%*[^\t]\t%u\t%u\t%[^\t]\t%*[^\t]\t%[^\t]\t%*[^\n]",
						in_pseudo_id, in_chromo, &in_begin, &in_end, in_parent_ensembl, in_parent_name))
	{
		uint chr_num;
		if (strcasecmp(in_chromo, "X") == 0)
		{  chr_num = 23;  }
		else if (strcasecmp(in_chromo, "Y") == 0)
		{  chr_num = 24;  }
		else
		{  chr_num = atoi(in_chromo);  }
		
		const BOOST_Interval coords(in_begin, in_end);
		
// 		psiDR_DB[in_pseudo_id] = Gene(in_pseudo_id, in_parent_name, in_chromo, BOOST_Interval(in_begin, in_end), in_parent_ensembl);
		psiDR_DB_map[chr_num][coords] =  Gene(in_pseudo_id, in_parent_name, chr_num, coords, in_parent_ensembl);
// 		std::cerr << "\t\tin_pseudo_id = [" << in_pseudo_id << "],  chr [" << chr_num << "]: " << coords << "\n";

		std::fscanf(file_psidr, "\n");
	}//EOF
	
	std:fclose(file_psidr);
	
	std::cerr << "DONE with \"load_psi_DR_pseudogene_database\"!\n";
	
	
	return psiDR_DB_map;
	
}//load_psi_DR_pseudogene_database








void save_simple_theoretical_stats
		(const type_map_uint_to_CC &Conn_Comps)
{
	std::cerr << "\nsave_simple_theoretical_stats...\n";
	
	std::string outfname(scratch_job_output_dir);
	outfname.append("/basic_Event_stats.txt");
	
	std::ofstream  outfs(outfname.c_str());
	check_fstream_for_goodness(outfs, outfname, "save_simple_theoretical_stats", true);
	
	for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{
		for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
			it_ev != it_cc->second.events.end();
			++it_ev)
		{
			if (it_ev->second.recomb_type != recomb_class__DupDel)
			{  continue;  }
			
			uint dist_from_centromere;
			uint dist_from_telomere;
			if (it_ev->second.arms[0])
			{
				dist_from_telomere = it_ev->second.region_between_and_including_the_LCRs_themselves.lower();
				dist_from_centromere = Event::centromere_coordinates.at(it_ev->second.chromos[0]).first - it_ev->second.region_between_and_including_the_LCRs_themselves.upper();
			}
			else
			{
				dist_from_telomere = Event::chromosome_lengths.at(it_ev->second.chromos[0]) - it_ev->second.region_between_and_including_the_LCRs_themselves.upper();
				dist_from_centromere = it_ev->second.region_between_and_including_the_LCRs_themselves.lower() - Event::centromere_coordinates.at(it_ev->second.chromos[0]).second;
			}			
			
			outfs 
				<< it_ev->second.UID << "\t"
				<< BOOST_width_inclusive(it_ev->second.region_between_and_including_the_LCRs_themselves) << "\t"
				<< BOOST_width_inclusive(it_ev->second.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs) << "\t"
				<< it_ev->second.profile_length << "\t"
				<< std::max<uint>(it_ev->second.LCR_lengths[0], it_ev->second.LCR_lengths[1]) << "\t"
				<< (longdouble)it_ev->second.variational_positions_of_profile.size() / it_ev->second.profile_length << "\t"
				<< dist_from_centromere << "\t"
				<< dist_from_telomere << "\t"
				<< BOOST_width_inclusive(it_ev->second.region_inbetween_LCRs_but_NOT_INCLUDING_the_LCRs) + BOOST_width_inclusive(it_ev->second.LCRs[0]) << "\n";
		}//ev
	}//cc
	
	
	outfs.close();	
	
}//save_simple_theoretical_stats
































void remove_undefined_from_potential_space
			(type_map_uint_to_CC &Conn_Comps,
			const bool &remove_entire_cc_too)
{
	std::string fname(data_directory);
	fname.append("/Events_and_Connected_Components_involving_Undefined_regions_of_the_genome");
	std::ifstream infs(fname);
	assert(infs.good());
	
	type_set_uint affected_ccs;
	type_set_uint affected_evs;
	
	uint in_ev;
	while (infs >> in_ev)
	{
		
		uint in_cc;
		infs >> in_cc;	
		
		affected_ccs.insert(in_cc);
		const type_map_uint_to_CC::iterator it_cc = Conn_Comps.find(in_cc);
		if (remove_entire_cc_too)
		{			
			if (it_cc != Conn_Comps.end())
			{
				const type_set_uint affected_evs_this_cc(extract_keys_of_map_and_return_as_set<uint, Event>(it_cc->second.events));
				affected_evs.insert(affected_evs_this_cc.begin(), affected_evs_this_cc.end());
			}			
			Conn_Comps.erase(in_cc);
		}
		else
		{
			if (it_cc != Conn_Comps.end())
			{
				affected_evs.insert(in_ev);
				it_cc->second.events.erase(in_ev);
			}
		}
	}//while
	
	infs.close();
	
	
	std::cerr << "\n\n\tnum affected ccs: " << affected_ccs.size() << "\n"
		<< "\n\tnum affected evs:  " << affected_evs.size() << "\n\n";
		
	
}//remove_undefined_from_potential_space
























void remove_calls_whose_Event_or_Connected_Component_has_been_removed_or_are_sex_chromo
			    (const type_map_uint_to_CC &Conn_Comps,
			     type_map_string_to_uint_to_Called_diploid_wrapper &mycalls)
{
        
    type_map_string_to_uint_to_Called_diploid_wrapper::iterator it_mygen = mycalls.begin();
    while(it_mygen != mycalls.end())	   
    {
	type_map_uint_to_Called_diploid_wrapper::iterator it_myev = it_mygen->second.begin();
	while (it_myev != it_mygen->second.end())
	{
		if (Conn_Comps.count(it_myev->second.cc_of_event) == 0 or Conn_Comps.at(it_myev->second.cc_of_event).events.count(it_myev->second.event_UID) == 0
			or Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID).chromos[0] > 22)
		{  
			it_mygen->second.erase(it_myev++);
		}
		else
		{  ++it_myev;  }
	}//myev
	
	if (it_mygen->second.empty())
	{  mycalls.erase(it_mygen++);  }
	else
	{  ++it_mygen;  }
    }//mygen
    
    
}//remove_calls_whose_Event_or_Connected_Component_has_been_removed_or_are_sex_chromo






void remove_non_dup_dels
		(type_map_uint_to_CC &Conn_Comps)
{
	type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
	while (it_cc != Conn_Comps.end())
	{
		type_map_uint_to_Event::iterator it_ev = it_cc->second.events.begin();
		while (it_ev != it_cc->second.events.end())
		{
			if (it_ev->second.recomb_type != recomb_class__DupDel)
			{
				it_cc->second.events.erase(it_ev++);
			}
			else
			{  ++it_ev;  }
		}//it_ev
		
		if (it_cc->second.events.empty())
		{
			Conn_Comps.erase(it_cc++);
		}
		else
		{  ++it_cc;  }
	}//it_cc	
		
}//remove_non_dup_dels















void upload_breakpoints_logodds
			(type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
			const std::string &brkpt_odds_fname)			
{
	std::cerr << "inside \"upload_breakpoints_logodds\"  with  brkpt_odds_fname = [" << brkpt_odds_fname << "]...\n";
	std::ifstream infs(brkpt_odds_fname);
	
	std::string ingen;
	while (infs >> ingen)
	{
		uint inev;
		uint inb0;
		uint inb1;
		real inbodds;
		
		infs >>  inev >> inb0 >> inb1 >> inbodds;
		
		if (mycalls.count(ingen) > 0  and   mycalls.at(ingen).count(inev) > 0)
		{  
			if (BOOST_equal(mycalls.at(ingen).at(inev).the_diploid_profile_brkpts.first, BOOST_Interval(inb0, inb1)))	
			{  mycalls.at(ingen).at(inev).brkpt_odds.first = inbodds;  }
			else
			{  mycalls.at(ingen).at(inev).brkpt_odds.second = inbodds;  }
		}
	}
	
	
	infs.close();		
	
	std::cerr << "DONE with \"brkpt_odds_fname\".";
		
}//upload_breakpoints_logodds

























uint calculate_distance_to_nearest_MEPS
			(const Event &the_Event,
			const uint &some_breakpoint,
			const uint &MEPS_length,
			const uint &MEPS_num_mismatches)
{	
	// distance from breakpoint
	for (int d = 0; d < the_Event.profile_length; ++d)
	{
		#pragma loop unroll(2)
		for (int sign_mult = -1; sign_mult <= 1; sign_mult+=2)
		{
			// form MEPs
			const int edge = (int)some_breakpoint + d*sign_mult;
			const int other_edge = edge + sign_mult*MEPS_length - sign_mult*1;
			
			if (edge < 0 or other_edge < 0 or edge >= the_Event.profile_length  or other_edge >= the_Event.profile_length)
			{  continue; }
			
			const BOOST_Interval hypothetical_MEPS((uint)std::min<int>(edge, other_edge), (uint)std::max<int>(edge, other_edge));
			
			
			//check identity			
			const uint num_mismatches = get_intersection_of_set_and_interval(the_Event.variational_positions_of_profile, hypothetical_MEPS).size();
			
			if (num_mismatches <= MEPS_num_mismatches)
			{
				//check large gaps
				bool intersects_some_large_gap = false;
				for (type_set_BI::const_iterator it_lg = the_Event.large_gaps_of_profile.begin();
					it_lg != the_Event.large_gaps_of_profile.end();
					++it_lg)
				{
					if (BOOST_overlap(*it_lg, hypothetical_MEPS))
					{
						intersects_some_large_gap = true;
						break;
					}
				}//lg
				
				if (!intersects_some_large_gap)
				{
					return d;
				}
			}//mismatches
		}//sign		
	}//d
	
	
	return the_Event.profile_length;
	
	
// 	for (uint p=0; p<the_Event.profile_length; ++p)
// 	{
// 		const BOOST_Interval hypothetical_MEPS(p, p + MEPS_length - 1);
// 		
// 		const uint num_mismatches = get_intersection_of_set_and_interval(the_Event.variational_positions_of_profile, hypothetical_MEPS).size();
// 		
//  		//if (get_intersection_of_set_and_interval(the_Event.variational_positions_of_profile, hypothetical_MEPS).empty())
// 		if (num_mismatches <= MEPS_num_mismatches)
// 		{
// 			//check large gaps
// 			bool intersects_some_large_gap = false;
// 			for (type_set_BI::const_iterator it_lg = the_Event.large_gaps_of_profile.begin();
// 				it_lg != the_Event.large_gaps_of_profile.end();
// 				++it_lg)
// 			{
// 				if (BOOST_overlap(*it_lg, hypothetical_MEPS))
// 				{
// 					intersects_some_large_gap = true;
// 					break;
// 				}
// 			}//lg
// 			
// 			if (!intersects_some_large_gap)
// 			{
// 				const uint observed_min_dist = (uint)std::min<int>(   std::abs(hypothetical_MEPS.lower() - some_breakpoint),
// 										      std::abs(hypothetical_MEPS.upper() - some_breakpoint));
// 				
// 				if (observed_min_dist < min_distance)
// 				{  min_distance = observed_min_dist;  }
// 			}
// 		}
// 	}//varpos
// 	
// 	
// 	return min_distance;	
	
}//calculate_distance_to_nearest_MEPS








type_map_uint_to_longdouble calculate_breakpoint_distribution_of_Event_for_MEPS
		(const Event &the_Event,
		const uint &MEPS_length,
		const uint &MEPS_num_mismatches,
		const int &probability_counting_method__0_uniform_by_vp____1_uniform_by_position)
{
	type_map_uint_to_longdouble  brkpt_dist;
	
	const std::vector<type_set_uint::const_iterator> vec_vp(get_loopable_constiterators<type_set_uint>(the_Event.variational_positions_of_profile));
	
	const double time_begin = omp_get_wtime();
	std::cerr << "\t\t\t\t\tcalculate_breakpoint_distribution_of_Event_for_MEPS    " << the_Event.UID << "...";
	#pragma omp parallel
	{//parallel
		type_map_uint_to_longdouble  brkpt_dist___thread;
		#pragma omp for schedule(dynamic,3)
		for (uint j=0; j<vec_vp.size(); ++j)
		{
			const uint dist_to_MEPS = calculate_distance_to_nearest_MEPS(the_Event, *vec_vp[j], MEPS_length, MEPS_num_mismatches);
			
			if (probability_counting_method__0_uniform_by_vp____1_uniform_by_position == 0)
			{
				brkpt_dist___thread[dist_to_MEPS] += 1.00L/the_Event.variational_positions_of_profile.size();			
			}//by vp
			else
			{
				if (j > 0)//the_Event.variational_positions_of_profile.begin())
				{
					const uint prev_vp = *vec_vp[j-1];
					
					brkpt_dist___thread[dist_to_MEPS] +=   1.00L*(*vec_vp[j] - prev_vp)/the_Event.profile_length; //(*--the_Event.variational_positions_of_profile.end());
				}//not beginning
				else
				{
					brkpt_dist___thread[dist_to_MEPS] +=   1.00L*(*vec_vp[j])/the_Event.profile_length;
				}//first
			}
		}//vec_brkpts_dist	
		
		#pragma omp critical(combine_MEPS_results_for_distances_for_Event)
		{//critical
			add_map_to_another_direct<uint, longdouble>(brkpt_dist, brkpt_dist___thread);
		}//critical
	}//parallel
	
	std::cerr << "  DONE!   " <<  omp_get_wtime() - time_begin << " s\n";
	
	
	return  brkpt_dist;
	
}//calculate_breakpoint_distribution_of_Event_for_MEPS














longdouble calculate_Lyapunov_CLT_statistic_for_MEPS
		(const type_map_uint_to_CC &Conn_Comps,
		 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		const longdouble &logfdr_threshold_for_brkpts,
		 const uint &MEPS_length,
		const uint &MEPS_num_mismatches,
		const int &probability_counting_method__0_uniform_by_vp____1_uniform_by_position,
		const type_map_uint_to_uint_to_longdouble &maps_ev_to_brkpt_dists)
{
	
	longdouble Lyap_CLT_variances = 0;
	longdouble Lyap_CLT_diffs = 0;

	for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
		it_mygen != mycalls.end();
		++it_mygen)
	{
		for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
			it_myev != it_mygen->second.end();
			++it_myev)
		{
			uint good_brkpt_value = 0;
			if (it_myev->second.brkpt_odds.first > 0  and  mpfr::log10(it_myev->second.brkpt_odds.first) <= logfdr_threshold_for_brkpts)
			{  good_brkpt_value = it_myev->second.the_diploid_profile_brkpts.first.lower();  }
			else if (it_myev->second.brkpt_odds.second > 0  and  mpfr::log10(it_myev->second.brkpt_odds.second) <= logfdr_threshold_for_brkpts)
			{  good_brkpt_value = it_myev->second.the_diploid_profile_brkpts.second.lower();  }
			
			if (good_brkpt_value > 0)
			{
				const uint observed_distance = calculate_distance_to_nearest_MEPS(
									Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID),
									good_brkpt_value, MEPS_length, MEPS_num_mismatches);			
				
				
// 				const type_map_uint_to_longdouble brkpt_dist(
// 						calculate_breakpoint_distribution_of_Event_for_MEPS(
// 								Conn_Comps.at(it_myev->second.cc_of_event).events.at(it_myev->second.event_UID),
// 								MEPS_length, MEPS_num_mismatches,
// 								probability_counting_method__0_uniform_by_vp____1_uniform_by_position));
				
				
				longdouble mean, var;
				calculate_mean_and_var(maps_ev_to_brkpt_dists.at(it_myev->second.event_UID), mean, var);
				
				Lyap_CLT_diffs += (longdouble)observed_distance - mean;
				Lyap_CLT_variances += var;		
				
			}//bodds
		}//ev
	}//gen
	
	std::cerr << "Lyap_CLT_diffs = " << Lyap_CLT_diffs << ",  denom = " << std::sqrt(Lyap_CLT_variances) << "  ";
	
	
	return   Lyap_CLT_diffs / std::sqrt(Lyap_CLT_variances);		
	
}//calculate_Lyapunov_CLT_statistic_for_MEPS








void calculate_mean_and_var
		(const type_map_uint_to_longdouble &some_dist,
		longdouble &mean, 
		longdouble &var)
		
{
	mean = 0;
	var = 0;
	
	for (type_map_uint_to_longdouble::const_iterator it_elt = some_dist.begin();
		it_elt != some_dist.end();
		++it_elt)
	{
		mean += it_elt->first * it_elt->second;
	}//elt
	
	for (type_map_uint_to_longdouble::const_iterator it_elt = some_dist.begin();
		it_elt != some_dist.end();
		++it_elt)
	{
		const longdouble diff = (it_elt->first - mean);
		var += diff*diff;
	}//elt	
	
	var = std::sqrt(var);
	
	
}//calculate_mean_and_sd







void calculate_various_MEPS_test_statistics_and_save
		(const type_map_uint_to_CC &Conn_Comps,
		 const type_map_string_to_uint_to_Called_diploid_wrapper &mycalls,
		const longdouble &logfdr_threshold_for_brkpts)
{
	std::cerr << "calculate_various_MEPS_test_statistics_and_save...\n";
	
	std::string outfname(scratch_job_output_dir);
	outfname.append("/MEPS_test_statistics.txt");
	
	std::ofstream outfs(outfname);
	check_fstream_for_goodness(outfs, outfname, "calculate_various_MEPS_test_statistics_and_save", true);
	
	boost::math::normal_distribution<longdouble> my_standard_normal(0,1);
	
	for (int prob_method = 0; prob_method <= 1; ++prob_method)
	{
		std::cerr << "\tprob_method = " << prob_method << "...\n";
		if (prob_method == 0)
		{  outfs << "uniform over variational positions:\n";  }
		else
		{  outfs << "uniform over genome positions:\n";  }
		
		
		type_map_uint_to_uint_to_longdouble maps_ev_to_brkpt_dists;
		
		//get evs:
		for (type_map_string_to_uint_to_Called_diploid_wrapper::const_iterator it_mygen = mycalls.begin();
			it_mygen != mycalls.end();
			++it_mygen)
		{
			for (type_map_uint_to_Called_diploid_wrapper::const_iterator it_myev = it_mygen->second.begin();
				it_myev != it_mygen->second.end();
				++it_myev)
			{	
				maps_ev_to_brkpt_dists[it_myev->second.event_UID];
			}
		}
		
		
		const std::vector<type_map_uint_to_uint_to_longdouble::iterator> ev_to_brkpt_dist_its = get_loopable_iterators<type_map_uint_to_uint_to_longdouble>(maps_ev_to_brkpt_dists);
		
		std::map<uint, const Event *> all_Event_ptrs;
		for (uint j=0; j<ev_to_brkpt_dist_its.size(); ++j)	
		{  all_Event_ptrs[ev_to_brkpt_dist_its[j]->first] = &find_Event(Conn_Comps, ev_to_brkpt_dist_its.at(j)->first)->second;  }
		
		
		for (uint MEPS_length = 500; MEPS_length <= 500; MEPS_length += 50)
		{
			std::cerr << "\t\tMEPS_length = " << MEPS_length << "...\n";
			for (uint MEPS_num_mismatches = 0; MEPS_num_mismatches <= 2; ++MEPS_num_mismatches)
			{
				std::cerr << "\t\t\tcalculing evs    MEPS_num_mismatches = " << MEPS_num_mismatches << "...\n";
				
// 				#pragma omp parallel for schedule(dynamic,1)
				for (uint j=0; j<ev_to_brkpt_dist_its.size(); ++j)				
// 				for (type_map_uint_to_uint_to_longdouble::iterator it_ev = maps_ev_to_brkpt_dists.begin();
// 					it_ev != maps_ev_to_brkpt_dists.end();
// 					++it_ev)
				{
					ev_to_brkpt_dist_its.at(j)->second = calculate_breakpoint_distribution_of_Event_for_MEPS(
								*all_Event_ptrs.at(ev_to_brkpt_dist_its.at(j)->first), MEPS_length, MEPS_num_mismatches, prob_method);
				}			
				
				std::cerr << "\t\t\tcalculing Lyap...\n";
				const longdouble Lyap_statistic 
						= calculate_Lyapunov_CLT_statistic_for_MEPS(
								Conn_Comps, mycalls, logfdr_threshold_for_brkpts,
								MEPS_length, MEPS_num_mismatches, prob_method, maps_ev_to_brkpt_dists);
				
				const longdouble p_value = boost::math::cdf<longdouble>(my_standard_normal, Lyap_statistic);
						
				outfs << Lyap_statistic << "  & \t";
				
				
			}//MEPS_num_mismatches
			
			outfs << "\\\\\n";
		}//MEPS_length
		
		outfs << "\n\n\n";
	}//method
	

	outfs.close();
	
	std::cerr << "\nDONE calculate_various_MEPS_test_statistics_and_save!\n\n";
	
	
}//calculate_various_MEPS_test_statistics_and_save





type_map_uint_to_Event::const_iterator find_Event
	(const type_map_uint_to_CC &Conn_Comps,
	const uint &the_ev)
{
	for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{
		if (it_cc->second.events.count(the_ev) > 0)
			return it_cc->second.events.find(the_ev);	
	}
	
	
	std::cerr << "DID NOT FIND even = " << the_ev << "!!!!!\n\n";
	exit(1);
}//find_CC_with_EV




// void cal
// 	(const Event &the_Event,
// 	const uint &MEPS_length)
// {
// 	type_vector_int num_mistmaches_this_MEPS;
// 	
// 	for (uint j=0; j < the_Event.profile_length)
// 	
// 	
// 	
// }//
