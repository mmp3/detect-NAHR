#ifndef IO_FUNCTIONS_H_INCLUDED
#define IO_FUNCTIONS_H_INCLUDED


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



#include <SD_Entry_general.h>
#include <general_typedefs.h>



#include "globals.h"





class  Call_set;
class  Conn_Comp;
class  Event;
class  Visitation_object;
class  Sampled_diploid_Event_data;
class  Gene;


typedef  std::map<uint, Sampled_diploid_Event_data>   type_map_uint_to_Sampled_diploid_Event_data;
typedef  std::list<Sampled_diploid_Event_data>  type_list_Sampled_diploid_Event_data;        
typedef  std::vector<type_map_uint_to_Sampled_diploid_Event_data>   type_vector_map_uint_to_Sampled_diploid_Event_data;


class  Marginal_Event_posterior_distribution;
typedef  std::map<uint, Marginal_Event_posterior_distribution>  type_map_uint_to_Marginal_Event_posterior_distribution;


void read_connected_component_data_from_file_and_create_Event_shells
                            (const char *const &in_connected_components_filename,
                             std::map<uint, Conn_Comp> &Conn_Comps);

                             
                             
void read_Event_ALL_data_from_file
                        ( Conn_Comp &some_Conn_Comp,
                          Event &some_Event);                             
                             
void read_Event_basic_data_from_file
                        ( Conn_Comp &some_Conn_Comp,
                          Event &some_Event); 
                          
void read_Event_interlocking_data_from_file
                        ( Conn_Comp &some_Conn_Comp,
                          Event &some_Event);
                          
void read_Event_regions_data_from_file
                        ( Conn_Comp &some_Conn_Comp,
                          Event &some_Event);                          
                             					

void read_Connected_Component_ALL_Event_data_from_file
                        (Conn_Comp &some_Conn_Comp);

void read_Connected_Component_basic_Event_data_from_file
                        (Conn_Comp &some_Conn_Comp);

void read_Connected_Component_interlocking_Event_data_from_file
                        (Conn_Comp &some_Conn_Comp);

void read_Connected_Component_regions_Event_data_from_file
                        (Conn_Comp &some_Conn_Comp);
				
               
                                
                                
                                
                              
void read_chromosome_data_from_file();


bool read_permissibile_state_vectors_from_file_for_summing_out_a_given_UID_and_tack_on_empty_state_vector
                                                      (Visitation_object &my_Visit_object);
               
                                                      
                                                      
int  read_total_diploid_compute_size_from_file_for_given_Conn_Comp
                                     (Conn_Comp &my_CC);                                                      
                                                      
bool read_visitation_schedule_from_file_for_given_Conn_Comp (Conn_Comp &my_CC);


std::string get_specific_DNA_sequence_from_Reference_genome(const uint &chromo,
                                               const BOOST_Interval &coords);

void upload_universal_variational_positions();

                                               
// void read_bas_file();      


void read_Validated_Events_from_file
                            (const char *const &validated_Events_filename,
                             type_list_Validated_rearrangement &Validated_Events);



void upload_universal_indel_positions();



void save_sampled_Connected_Component_to_file
                        (FILE *&out_file,
                         const  type_map_uint_to_Sampled_diploid_Event_data  &sample);
                         
void load_sampled_Event_from_file
                        (FILE *&out_file,
                         type_map_uint_to__uint__uint__2 &sample);                         


bool read_calls_from_file
                    (const std::string &calls_filename,
                     std::map<std::string, std::list<Call_set> >  &read_calls,
		     const bool &ignore_gene_conversion_calls,
		     const bool &ignore_second_call_in_diploid_calls  );

                     
bool read_selected_connected_components_file
                    (const std::string &selection_dirname,
                     type_set_uint &magic_CCs);                   
                     
                     
                     
void append_to_skipped_connected_component_file
                (const Conn_Comp  &skipped_CC,
                 const std::string &reason_for_skipping);                     
                     

void read_Validated_rearrangements_from_file__in_tabular_format
                (const std::string in_validations_filename,
                 type_list_Validated_rearrangement &validated_Events);
                     



void  read_gender_and_population_from_file
                (const std::string &gender_and_population_filename);          
                 
                     

                
void save_posterior_samples_to_file
                        (const type_vector_map_uint_to_Sampled_diploid_Event_data &samples_from_posterior,
                         const uint &CC_ID,                         
                         const uint &number_of_samples,
                         const std::string &some_output_dir);
                         
type_vector_map_uint_to_Sampled_diploid_Event_data  load_posterior_samples_from_file
                        (const std::string &some_output_dir,
                         const uint &CC_ID,
                         const uint &number_of_samples);     
			
			
type_vector_map_uint_to_Sampled_diploid_Event_data load_posterior_samples_from_file
                        (const std::string &the_fname);			
                
                

void save_marginal_posterior_distributions_to_file
                        (const type_map_uint_to_Marginal_Event_posterior_distribution &marginal_distributions,
                         const std::string &some_output_dir);                     
                
                        
                        

void save_calls_tables_to_file
                    (const std::stringstream &output_marginal_Positive_calls_strm,
                    const std::stringstream &output_marginal_Negative_calls_strm,
                    const std::stringstream &output_marginal_all_calls_strm,
                    const std::string &some_output_dir);
                    
                    
		    
type_set_uint  read_totally_excluded_Events_from_previous_job_save_dir
			    (const std::string &old_run_data_dir);            
                

void remove_MPI_rank_directories_from_scratch_space();			    
void remove_directories_associated_with_Partial_Sum_functions_from_scratch_space();

void create_scratch_directory_tree();
void create_output_directory_tree();
			    
	

void save_excluded_Events_to_file
		    (const type_set_uint &excluded_Events);

void save_remaining_original_breakpoints_to_file
		(const std::map<uint, Event>  &the_events_of_some_Conn_Comp);
		
// void  load_original_breakpoints_after_heuristic_for_all_Events_in_Connected_Component_from_file 
// 		    (const std::string &old_run_data_dir,
// 		     std::map<uint, Event>  &events_of_some_CC);
		     
	
void  save_Events_affected_by_Undefined_regions_to_file
		(const std::map<uint, Event>  &events_of_some_CC,
		 const type_map_uint_to_list_BI &ignored_regions_fo_genome);		     
		     

type_map_uint_to_list_BI load_Non_unique_regions_of_genome();		 
		 
	
type_set_uint upload_connected_components_involving_Undefined_regions_of_genome
			(const std::string &the_data_dir);

	

type_set_string   read_genome_list_from_file
		(const std::string  &file_containing_genome_names);			
			
	



template <typename TYPEkey, typename TYPEval>
void  write_map_to_file
	    (const std::string &out_fname,
	     const std::map<TYPEkey, TYPEval>  &the_map)
{
    std::ofstream  out_fs(out_fname.c_str());
                    if ( !out_fs.good()  )
                    {
                        std::stringstream error_strm;
                        print_line_of_markers("ERROR! ", &error_strm);
                        print_line_of_markers("(", &error_strm);
                        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                        
                        error_strm << "\n\nERROR!   not good ofstream   in  \"write_map_to_file\",\n\t\t out_fname =    "  
                        <<  out_fname  << "\n\n";
                        
                        print_line_of_markers(")", &error_strm);
                        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() ); 
			return;
                    }    
    
    for (typename std::map<TYPEkey, TYPEval>::const_iterator it_m = the_map.begin();
	    it_m != the_map.end();
	    ++it_m)
    {
	out_fs << it_m->first  << "\t" << it_m->second << "\n";    
    }
    
    
    out_fs.close();

}//write_map_uint_to_longdouble_to_file		







template <typename TYPEkey>
void  write_set_to_file
	    (const std::string &out_fname,
	     const std::set<TYPEkey>  &the_set,
	    const std::string &delimiter = std::string("\t"))
{
    
    std::ofstream  out_fs(out_fname.c_str());
    if (check_fstream_for_goodness(out_fs, out_fname, "write_set_to_file", false))
    {
	for (typename std::set<TYPEkey>::const_iterator it_m = the_set.begin();
		it_m != the_set.end();
		++it_m)
	{
	    out_fs << *it_m << delimiter;    
	}	
	
	out_fs.close();
    }//good

}//write_set_to_file





	
	
	
template <typename TYPEkey, typename TYPEval>
std::map<TYPEkey, TYPEval>   read_map_from_file
				(const std::string &in_fname,
				 const bool &file_contains_a_header)
{
    std::map<TYPEkey, TYPEval>  inmap;    
    
    std::ifstream in_fs(in_fname.c_str());
                    if ( !in_fs.good()  )
                    {
                        std::stringstream error_strm;
                        print_line_of_markers("ERROR! ", &error_strm);
                        print_line_of_markers("(", &error_strm);
                        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                        
                        error_strm << "\n\nERROR!   not good ifstream   in  \"read_map_from_file\",\n\t\t in_fname =    "  
                        <<  in_fname  << "\n\n";
                        
                        print_line_of_markers(")", &error_strm);
                        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() ); 
			exit(1);
                    }        
    
    
    if (file_contains_a_header)
    {
	std::cerr << "\n\treading header...";
	char dummy_header[option_size];
	in_fs.getline(dummy_header, option_size-1, '\n');
	std::cerr << "DONE WITH HEADER!!!\n\n";
    }
    
    TYPEkey inkey;
    TYPEval inval;                
    while (in_fs >> inkey)
    {
	if (in_fs >> inval)
	{  inmap[inkey] = inval;  }
	else
	{  break;  }
    }//while
    
    
    
    in_fs.close();
    
    return inmap;    

}//read_map_from_file
	
	
	
	
	
	


		    
		    
		    
void append_diploid_Sampled_Events_to_file
		(const type_map_uint_to_Sampled_diploid_Event_data &sampled_Events,
		const std::string  &outdir);		    
	
	
std::map<std::string, std::map<uint, Sampled_diploid_Event_data> >  read_Sampled_Events_from_file
									    (const std::string &fname);		
		
		
									    
template <typename Tcontainer> 
void write_singlecontainer_to_file
				    (const Tcontainer &some_container,
				     const std::string &full_fname = std::string(""),
				    std::ofstream *const &some_ofstream_ptr = NULL)
{
    std::ofstream outfs;       
    std::ofstream *const outfs_ptr = (some_ofstream_ptr == NULL)  ?    &outfs  :  some_ofstream_ptr;
    
    if (some_ofstream_ptr == NULL)
    {  outfs_ptr->open(full_fname);  }

    if (check_fstream_for_goodness(*outfs_ptr, full_fname, "write_singlecontainer_to_file", false))
    {  
	for (typename Tcontainer::const_iterator it = some_container.begin();
		it != some_container.end();
		++it)
	{
	    (*outfs_ptr) << *it << "\n";   	    
	}//it
    }         
    
    if (some_ofstream_ptr == NULL)
    {  outfs_ptr->close();  }
    
}//write_singlecontainer_to_file









template <typename Telt>
std::list<Telt> read_list_from_file
			(const std::string &filename)
{
    std::list<Telt> list_from_file;
    
    std::ifstream infs(filename);
    
    Telt inelt;
    while (infs >> inelt)
    {
	list_from_file.push_back(inelt);		
    }
    
    infs.close();
    
    
    return list_from_file;    
    
}//read_singlecontainer_from_file






template <typename Telt>
std::set<Telt> read_set_from_file
			(const std::string &filename)
{
    std::set<Telt> set_from_file;
    typename std::set<Telt>::iterator insert_it = set_from_file.begin();
    
    std::ifstream infs(filename);
    
    Telt inelt;
    while (infs >> inelt)
    {
	insert_it = set_from_file.insert(insert_it, inelt);
    }
    
    infs.close();
    
    
    return set_from_file;    
    
}//read_singlecontainer_from_file






void load_cumulative_GC_counts_for_each_chromosome();

void load_cumulative_GC_counts_for_specific_chromosome
				(const uint &chromo,
				type_vector_uint &count_vector);



type_map_string_to_string__bool load_gender_and_population_map
					(const std::string &gender_and_pop_filename);
									    
	
#endif		