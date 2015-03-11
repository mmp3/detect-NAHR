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



// #include <openmpi/mpi.h>
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
// #include "/usr/include/gmp-i386.h"
// #include "/usr/include/gmpxx.h"
// #include "/usr/include/gmp.h"


#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamAux.h>


#include <general_typedefs.h>
#include <templates.h>
#include <Readgroup_statistics.h>


#include "Conn_Comp.h"
#include "Event.h"
#include "io_functions.h"
#include "globals.h"
#include "Visitation_object.h"
#include "other_functions.h"

#include "PER_aligner.h"
// #include "Call_set.h"
#include <boost/graph/graph_concepts.hpp>

#include "post.h"


int main(int argc, char** argv)
{    
                
    // setup MPI 

    
    boost::mpi::environment my_boost_mpi_env(argc, argv);
    boost::mpi::communicator world_MPI_boost;        
    world_MPI_boost_ptr = &world_MPI_boost;
    
    
    my_MPI_rank = world_MPI_boost.rank();
    size_MPI_WORLD = world_MPI_boost.size(); 
    
    if (my_MPI_rank == 0)
    {  std::cerr <<"\n\nrank " << my_MPI_rank << "   size_MPI_WORLD = " << size_MPI_WORLD << "\n\n";  }
    
    
    
    if (size_MPI_WORLD > 1)
    {
        if (my_MPI_rank == 0)
            std::cerr << "\n\n\n\ntesting MPI:  rank 0 will send the value [71] to rank 1...\n\n";
        
        
        if (my_MPI_rank == 0)
        {
            world_MPI_boost_ptr->send<int>(1, 0, 71);
            std::cerr << "\n\t rank 0 has sent value [71] to rank 1...\n";
        }
        
        if (my_MPI_rank == 1)
        {
            int recv_val = -1;
            world_MPI_boost_ptr->recv<int>(0, 0, recv_val);
            std::cerr << "\n\t rank 1 has received value [" << recv_val << "] from rank 1.\n";
        }
        
        std::cerr << "\n\ndone with MPI testing.\n";        
    }  //test
    
    
    
    
    
    
    
    
    
    
    // setup OMP
    
    max_number_of_OMP_threads = (uint)omp_get_max_threads();
    std::fprintf(stderr, "\n\nrank %d:  max_number_of_OMP_threads = %u\n\n", my_MPI_rank, max_number_of_OMP_threads);
    
    std::fflush(stderr);
    world_MPI_boost_ptr->barrier();    
    
    
    
    

  
  
// IMPORTANT NOTE:  "IsReverseStrand" doesn't mean that the nucleotide sequence given actually maps to the reverse strand.  The nucleotide sequence displayed always maps to the forward strand (at least, it does with WGSIM).  I suppose this is just extra stoed information regarding the read characteristics.  Thus, we don't ever need to consider "IsReverseStrand:", since it doesn't affect hte nucleotide sequence that we are actually given.  i.e. never complement or reverse.  This is manifsted in the function of class  Paired_end_read  called   "orient_and_complement_each_mate_according_to_its_spawning_Event_profile".
  
  
// CAREFUL WHEN YOU ARE UPDATING DURING THE PERFORMV VARIABLE ELIMINATION SCHEDULE.  we need to be careful of which profile positions are completely deleted, vs. which ones are only deleted for read-depth or PERs.  also check various member variables of Event.  for example, in Event::aggregate_all_relevant_breakpoints_for_each_interlocking_event_and_self_according_to_relevant_brkpts_for_PERs, we draw from   variational_positions.  Thus, as we go down the eliminaition schedulee, we need to be careful to eliminate positions from "variational_positions" as well.
  
  
// Suppose Event E is a dup/del; let's think about just the deletion case for now.  Suppose F is some other event (maybe a translocation).

//  Event E:      [     LCR_0     ]-----------------------------------------------------------[     LCR_1     ]

//  Event F:                                                         [ LCR_0 ]
//                                                                   [ LCR_1 ]     
//
//   Suppose that we sum out Event E first, and then event F at some later step.
//
//      NOTE, importantly, that when we sum out Event E, we will consider the Read Depth for the paired-end reads (PERs) that align to the deleted region of E - that is, ALL of the PERs that align to event F's profile - but for these PERs, we will NOT consider their hybriud probabilities (alignments/mappings).  Thus, if we permanently erase these positons from consideration, then when we get to eliminating Event F, there will be no PERs to consider for its hybrid mappings!
//      Thus we must be careful in how we keep track of "eliminated" genome/profile positions and "eliminated" PERs!
//      We see from the above schematic example that the read-depth may be evaluated in advance for certain positions ("in advance" meaning that the red depth for Event F is considered before eliminating Event F - actually during Event E).  Thus, the coordinates for which we consider read depth should be smaller than the coordinates over which we consider hybrid PERs.     
  

    
//     {
// 	const mpfr_class the_number_one(1);
//         const mpfr_class expected_lambda(2453.749602);   
// 	const uint observed = 2956;
//         const mpfr_class observed_value(observed);		
//         
//         // CAREFULLY DERIVED VERSION: 	
//         const mpfr_class parameter_BOOST_r( the_number_one / (double)RD_variance_scale___alpha_squared);
//                 //BOOST_r  <==>  my_k  <==> Gamma shape parameter k.                                                 
//         
//         const mpfr_class parameter_BOOST_p(the_number_one / (the_number_one  +  (expected_lambda*(double)RD_variance_scale___alpha_squared)));
//         
//         
//         const boost::math::negative_binomial_distribution<mpfr_class> my_negative_binomial(parameter_BOOST_r, parameter_BOOST_p);                
//         const mpfr_class result_pdf(boost::math::pdf<mpfr_class>(my_negative_binomial, observed_value));	
// 	
// 	        
// 
//     std::cerr <<"\n\nresult_pdf = " << result_pdf
// 	    << "\n\ndouble check:\n"
// 	    << "\nparameter_BOOST_r = " << parameter_BOOST_r
// 	    << "\nparameter_BOOST_p = " << parameter_BOOST_p
// 	    << "\nobserved_value = " << observed_value
// 	    << "\nboost::math::mean<mpfr_class>(my_negative_binomial) = " << boost::math::mean<mpfr_class>(my_negative_binomial)
// 	    << "\nboost::math::variance<mpfr_class>(my_negative_binomial) = " << boost::math::variance<mpfr_class>(my_negative_binomial)
// 	    << "\n";   
// }
//     return 0; 
    
    
    
    
    
    
    
    

    // setup REAL    
        
    mpfr::mpreal::set_emax( mpfr::mpreal::get_emax_max() );
    mpfr::mpreal::set_emin( mpfr::mpreal::get_emin_min() );
    mpfr::mpreal::set_default_prec(8192);
    
    minimum_representable_positive_real_number = 0.75L;
    minimum_representable_positive_real_number.set_exp( (-1)*(int)(10000000000LL) );  //10 billion // mpfr::minval();
    for (uint j=0; j < 6; ++j)
        minimum_representable_positive_real_number *= minimum_representable_positive_real_number;
    //we need a number that is smaller than any "read-depth" outcome we are likely to see, and smaller than any product-over-reads that we are likely to see.  But we don't want it SO small (i.e. = mpfr::minval()) that it will cause underflow (i.e. "nan" values) when multiplied by other numbers and propagated through the program.
    
    
    //default values:    
    if (my_MPI_rank == 0)
    {
        std:: cerr << "\n\n\n\tMPFR_PREC_MAX = " << MPFR_PREC_MAX 
                    << "\nmpfr::mpreal::get_default_prec() =  " << mpfr::mpreal::get_default_prec()
                    << "\n\tmpfr::mpreal::get_emax_max() = " << mpfr::mpreal::get_emax_max()
                   << "\n\tmpfr::mpreal::get_emin_min() = " << mpfr::mpreal::get_emin_min()
                   << "\n\tmpfr::mpreal::get_emax() = " << mpfr::mpreal::get_emax()
                   << "\n\tmpfr::mpreal::get_emin() = " << mpfr::mpreal::get_emin()
                   << "\n\tmpfr::minval() = minimum_representable_positive_real_number = " << minimum_representable_positive_real_number
                   << "\n\n\n\n";
    }
                   
    std::cerr<< "" <<  ""  << "\n\n";
    std::fflush(stderr);
    world_MPI_boost_ptr->barrier();
    
    
    
    
    std::cerr.precision(10);    
    
    
    
    set_prior_Event_outcome_probabilities(0.9999999L, 0.999L);
    
    // variables
    
    std::string BAM_filename;    
    
    
    BamTools::BamReader  my_BAM_reader;
    
    
    std::string display_MAP_dir;        
        
    type_list_Sampled_diploid_Event_data  display_Event_and_breakpoint_range;
    bool display_ALL_reads_not_just_good_ones = false;
    
    BOOST_Interval  acceptable_Connected_Component_size_range(empty_BI);
    BOOST_Interval  acceptable_Connected_Component_compute_range(empty_BI);            
    type_set_uint targeted_CC_IDs;
    
    uint number_of_samples_from_posterior = 0;
        
    int number_of_simulated_reads = 0;                      
    
    int last_completed_CC = -1;            
    
    bool used_debug_view = false;
    
    bool  only_compute_and_save_natural_poisson_intervals = false;        
    
    bool  do_post_processing_stats_on_calls = false;
    std::string  dir_of_calls;
        
    std::string validations_dirname;                  
    
    std::string gender_and_population_filename;       
        
    bool copy_data_from_most_recent_job_of_same_name = false;     
    std::string copy_data_from_job_with_jobID;
    
    bool did_a_compare_to_vals = false;
        
    std::string use_previously_computed_likelihoods___old_data_dir; 
       
    type_list_int  list_of_prior_probability_exponents_for_computing;        
    
    type_list_Sampled_diploid_Event_data  pseudo_RD_hypothesis_test_work;   

    type_list_uint__uint__uint  unique_regions_for_PER_displays;                

    type_list__uint__uint  informativeness_checking_parameters;    
    
    type_list_aCGH_validation_container   acgh_data_checks;
    
    bool consider_GeneConversion_breakpoints_in_Heuristic = false;
    bool consider_GeneConversion_in_Full_compute = true;
    
    
    
    std::string Genes_in_human_genome__filename;
    std::string biological_validations_filename;
    std::string postdatadir;
    
    std::string redo_brkpts_dir;
    std::string redo_RD_dir;
    
    bool do_raw_count_ratios__strategy_A = false;
	bool do_raw_count_ratios__strategy_B = false;
	
	
	bool make_simulated_genomes = false;
	uint num_sim_genomes = 0;
	uint num_ev_per_sim_genome = 0;
	std::string simulated_genome_base_fname;
	std::string allowable_events_fname;
	
    
    //  OPTIONS
    
    
    char *option = new char[10*option_size];      
    
    for (int arg_indx=0; arg_indx<argc; ++arg_indx)
    {
        char *temp_input = new char[10*option_size];
        
        unsigned int arg_indx_char = 0;  
        while(argv[arg_indx][arg_indx_char] != '='
            && argv[arg_indx][arg_indx_char] != '\0' //'\0' = null character for terminating strings            
            && arg_indx_char < (option_size-1) ) // leave room for a '\0' at end.
        {
            option[arg_indx_char] = argv[arg_indx][arg_indx_char];
            ++arg_indx_char;    
        }
    
        option[arg_indx_char] = '\0'; //null-terminate for safety.
        
        if (argv[arg_indx][arg_indx_char] != '\0')
            ++arg_indx_char;    
        
        
        
        
        
        if (strcmp(option,"-d") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            data_directory.assign(temp_input);
            
            if ( data_directory.empty() )
            {
                std::printf("ERROR: data_directory not specificed.\n\n");
                return -1;
            }     
            else
                if (my_MPI_rank == 0)
                    std::fprintf(stderr, "\n\tdata_directory set to [%s]\n", data_directory.c_str() );
        } //end if "-data_directory"    
        
        
        
        
        
        else if (strcmp(option,"-R") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            Ref_genome_fileprefix.assign( temp_input );
            if (Ref_genome_fileprefix.empty())
            {
                std::printf("ERROR: Reference genome directory not specificed.\n\n");
                //show_how_to_use(argv[0]);
                return -1;
            }      
            else
                if (my_MPI_rank == 0)
                    std::fprintf(stderr, "\n\tReference genome file prefix set to [%s]\n", Ref_genome_fileprefix.c_str() );            
        } //end if "-Reference genome"     
                

        
        
        
        else if (strcmp(option,"-b") == 0)
        {
            char tmp_in_filename_bam[option_size];
            strcpy(tmp_in_filename_bam, &argv[arg_indx][arg_indx_char]);
            
            BAM_filename.assign(tmp_in_filename_bam);
            
            if ( BAM_filename.empty() )
            {
                std::fprintf(stderr, "ERROR: BAM file not specificed.\n\n");
                //show_how_to_use(argv[0]);
                return -1;
            }     
            else
                if (my_MPI_rank == 0)                                             
                    std::fprintf(stderr, "\n\tBAM file set to [%s]\n", BAM_filename.c_str() );
        } //end if "-in_BAM_filename" 
    
    
    
    
    
//         else if (strcmp(option,"-min") == 0)
//         {
//             strcpy(temp_input, &argv[arg_indx][arg_indx_char]);            
//             min_CC_size = atoi(temp_input);   
//             if (my_MPI_rank == 0)
//                 std::fprintf(stderr, "\n\tset min_CC_size = %u\n", min_CC_size);             
//         
//         } //end if "-min"     
//         
//         
//         
//         else if (strcmp(option,"-max") == 0)
//         {
//             strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
//             
//             max_CC_size = atoi(temp_input); 
//             if (my_MPI_rank == 0)
//                 std::fprintf(stderr, "\n\tset max_CC_size = %u\n", max_CC_size);                 
//         
//         } //end if "-max"     
        
        
        
        else if (strcmp(option,"-CC_size") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
                        
            uint in_min_CC_size, in_max_CC_size;
            
            if (  EOF  !=   std::sscanf(temp_input, "%u,%u",  &in_min_CC_size, &in_max_CC_size )     )            
                acceptable_Connected_Component_size_range.set(  in_min_CC_size, in_max_CC_size );
            else
            {
                std::cerr << "\n\n\n\nFailed to read    \"acceptable_Connected_Component_size_range\"  correctly !!!!!! \n\n\n\n";
                exit(1);
            }            
            
            if (my_MPI_rank == 0)
                std::cerr <<  "\n\tset acceptable_Connected_Component_size_range = " << acceptable_Connected_Component_size_range << "\n";        
        
        } //end if "-CC_size"  
        
        
        
        
        else if (strcmp(option,"-compute_size") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
                        
            uint in_min_compute_size, in_max_compute_size;
            
            if (  EOF  !=   std::sscanf(temp_input, "%u,%u",  &in_min_compute_size, &in_max_compute_size )     )            
                acceptable_Connected_Component_compute_range.set(  in_min_compute_size, in_max_compute_size );
            else
            {
                std::cerr << "\n\n\n\nFailed to read    \"acceptable_Connected_Component_compute_range\"  correctly !!!!!! \n\n\n\n";
                exit(1);
            }            
            
            if (my_MPI_rank == 0)
                std::cerr <<  "\n\tset acceptable_Connected_Component_compute_range = " << acceptable_Connected_Component_compute_range << "\n";        
        
        } //end if "-compute_size"         
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-magic") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
	    
	    std::string in_targeted_str(  &argv[arg_indx][arg_indx_char]  );
	    
	    while ( !  in_targeted_str.empty()  )
	    {
		const uint firstcomma = std::min<uint>(in_targeted_str.find(','), in_targeted_str.size() );
		const uint desired_CC = (uint)atoi(in_targeted_str.substr(0, firstcomma).c_str());
		
		targeted_CC_IDs.insert(  desired_CC    );
                if (my_MPI_rank == 0)
                    std::fprintf(stderr, "added target CC ID = %u\n",   desired_CC  );				
		
		if (firstcomma == in_targeted_str.size() or  firstcomma == in_targeted_str.size()-1)	
		    in_targeted_str.clear();
		else
		    in_targeted_str = in_targeted_str.substr(firstcomma+1);		
	    }
            
//             if (strlen(temp_input) != 0)
//             {
//             
//                 const uint some_target_CC_ID = (uint)atoi(temp_input);   
//                 targeted_CC_IDs.insert(some_target_CC_ID);
//                 
//                 if (my_MPI_rank == 0)
//                     std::fprintf(stderr, "added some_target_CC_ID = %u\n", some_target_CC_ID);
//             }
        
        } //end if "-magic"    
    
    
    
        
        
        
        else if (strcmp(option,"-tag") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            genome_tag = (uint)atoi(temp_input);   
            
            if (my_MPI_rank == 0)
                std::fprintf(stderr, "set tag = %u\n", genome_tag);
        
        } //end if "-tag"    
        
        

        else if (strcmp(option,"-n") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            number_of_simulated_reads = (int)atoi(temp_input); 
            
            if (my_MPI_rank == 0)
                fprintf(stderr, "set number_of_simulated_reads = %d\n", number_of_simulated_reads);
        
        } //end if "-number_of_simulated_reads"    
        
        
        
        else if (strcmp(option,"-w") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            maximum_average_fragment_length = (uint)atoi(temp_input);
            
            if (my_MPI_rank == 0)
                fprintf(stderr, "set maximum_average_fragment_length = %u\n", maximum_average_fragment_length);
        
        } //end if "-average_fragment_length"
        
        
        
        
        
        
        else if (strcmp(option,"-l") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            length_of_haploid_genome = (uint)atoi(temp_input);   
            
            if (my_MPI_rank == 0)
                fprintf(stderr, "set length_of_haploid_genome = %u\n", length_of_haploid_genome);
        
        } //end if "-average_fragment_length"   
    
    
    
    
    
    
        else if (strcmp(option,"-e") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
                      
            
//             PER_aligner::likelihood_P_single_nucleotide_error  = (longdouble)atof(temp_input);  
//             PER_aligner::likelihood_1_minus_P_single_nucleotide_error  =  1.00L - PER_aligner::likelihood_P_single_nucleotide_error;
//             
//             PER_aligner::log10_likelihood_P_single_nucleotide_error  =  log10(PER_aligner::likelihood_P_single_nucleotide_error);
//             PER_aligner::log10_likelihood_1_minus_P_single_nucleotide_error = log10(PER_aligner::likelihood_1_minus_P_single_nucleotide_error);
//             
//             if (my_MPI_rank == 0)
//                 fprintf(stderr, "set likelihood error rate = %1.10Lf\n", PER_aligner::likelihood_P_single_nucleotide_error);
            
        } //end if "-average_fragment_length"       
        
        
        
        
        
        
        else if (strcmp(option,"-p") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
	    const real p_no_event__new = std::string(temp_input);    
	    set_prior_Event_outcome_probabilities(p_no_event__new, 0.999L);
	    
        } //end if "-average_fragment_length"      
        
        
        
        
        
        
        else if (strcmp(option,"-display") == 0)
        {                       
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
	    uint in_CC;
            uint in_EV;
	    uint in_profile_brkpt_low__hap0;
	    uint in_profile_brkpt_upp__hap0;	    
	    char in_outcome__hap0[option_size];
	    uint in_profile_brkpt_low__hap1;
	    uint in_profile_brkpt_upp__hap1;	    
	    char in_outcome__hap1[option_size];
	    
            std::sscanf(temp_input, "%u,%u,[%u,%u],%[^,],[%u,%u],%s", 
				    &in_CC, &in_EV,
				    &in_profile_brkpt_low__hap0, &in_profile_brkpt_upp__hap0, in_outcome__hap0,
				    &in_profile_brkpt_low__hap1, &in_profile_brkpt_upp__hap1, in_outcome__hap1); 
			
    
            display_Event_and_breakpoint_range.push_back(
		    Sampled_diploid_Event_data(in_CC, in_EV,
					    type_haploid_outcome__haploid_outcome(convert_string_to_haploid_outcome(in_outcome__hap0),
										  convert_string_to_haploid_outcome(in_outcome__hap1)),
					    type_BI__BI(BOOST_Interval(in_profile_brkpt_low__hap0,in_profile_brkpt_upp__hap0),
							BOOST_Interval(in_profile_brkpt_low__hap1,in_profile_brkpt_upp__hap1))
					    )); 

	    if (my_MPI_rank == 0)
	    {
		std::cerr << "\n\twill display custom result:\n";
		display_Event_and_breakpoint_range.back().print_this_Sampled_diploid_Event_data();	    
	    }    	    
          
        } //end if "-display"             
        
                
                
	
        else if (strcmp(option,"-display_ALL") == 0)
        {                       
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
	    uint in_CC;
            uint in_EV;
	    uint in_profile_brkpt_low__hap0;
	    uint in_profile_brkpt_upp__hap0;	    
	    char in_outcome__hap0[option_size];
	    uint in_profile_brkpt_low__hap1;
	    uint in_profile_brkpt_upp__hap1;	    
	    char in_outcome__hap1[option_size];
	    
            std::sscanf(temp_input, "%u,%u,[%u,%u],%[^,],[%u,%u],%s", 
				    &in_CC, &in_EV,
				    &in_profile_brkpt_low__hap0, &in_profile_brkpt_upp__hap0, in_outcome__hap0,
				    &in_profile_brkpt_low__hap1, &in_profile_brkpt_upp__hap1, in_outcome__hap1); 
			
    
	    display_ALL_reads_not_just_good_ones = true;
            display_Event_and_breakpoint_range.push_back(
		    Sampled_diploid_Event_data(in_CC, in_EV,
					    type_haploid_outcome__haploid_outcome(convert_string_to_haploid_outcome(in_outcome__hap0),
										  convert_string_to_haploid_outcome(in_outcome__hap1)),
					    type_BI__BI(BOOST_Interval(in_profile_brkpt_low__hap0,in_profile_brkpt_upp__hap0),
							BOOST_Interval(in_profile_brkpt_low__hap1,in_profile_brkpt_upp__hap1))
					    )); 

	    if (my_MPI_rank == 0)
	    {
		std::cerr << "\n\twill display custom result:\n";
		display_Event_and_breakpoint_range.back().print_this_Sampled_diploid_Event_data();	    
	    }    	    
          
        } //end if "-display"    	
	
	
	
	
	
	
	
	
	
	
	
	
                
                
                
                
                
        else if (strcmp(option,"-s") == 0)
        {            
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            scratch_job_output_dir.assign( temp_input );
            if (scratch_job_output_dir.empty())
            {
                std::fprintf(stderr, "ERROR: scratch_job_output_dir not specificed.\n\n");
                //show_how_to_use(argv[0]);
                return -1;
            }     
            else
                if (my_MPI_rank == 0)
                    std::fprintf(stderr, "\n\tscratch_job_output_dir directory set to [%s]\n", scratch_job_output_dir.c_str() );
        } //end if "-s"         
        
        
        

        
        
        
        else if (strcmp(option,"-samples") == 0)
        {
            
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            number_of_samples_from_posterior  =  (uint)atoi(temp_input);
            
            
            if (strlen(temp_input) == 0)
            {
                std::fprintf(stderr, "ERROR: \"number_of_samples_from_posterior\" not specificed.\n\n");
                exit(1);
            }
            else
                if (my_MPI_rank == 0)
                    std::fprintf(stderr, "\n\tset  number_of_samples_from_posterior  =  %u\n", number_of_samples_from_posterior);
                
                
        } //end if "-samples"            
        
        
        

        
        
        
//         else if (strcmp(option,"-v") == 0)
//         {
//             
//             strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
//             validated_Events_filename.assign( temp_input );
//             if (validated_Events_filename.empty())
//             {
//                 std::fprintf(stderr, "ERROR: validated_Events_filename not specificed.\n\n");
//                 return -1;
//             }     
//             else
//                 if (my_MPI_rank == 0)
//                     std::fprintf(stderr, "\n\tvalidated_Events_filename set to [%s]\n", validated_Events_filename.c_str() );
//         } //end if "-v"
        
        
        



        else if (strcmp(option,"-alpha") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            const longdouble scale_parameter__alpha = (longdouble)atof(temp_input);
            RD_variance_scale___alpha_squared = scale_parameter__alpha*scale_parameter__alpha;
            
            if (my_MPI_rank == 0)            
                std::cerr<< "set   RD_variance_scale___alpha_squared  =  "  <<  RD_variance_scale___alpha_squared <<  "  =   (" << scale_parameter__alpha << ")^2\n";
            
            if (my_MPI_rank == 0   and   RD_variance_scale___alpha_squared >= 1.00L )
            {
                std::stringstream error_strm;
                print_line_of_markers("ERROR! ", &error_strm);
                print_line_of_markers("(", &error_strm);
                error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                
                error_strm << "\n\nERROR!    RD_variance_scale___alpha_squared   =   "  <<  RD_variance_scale___alpha_squared  << "  >=  1.00 !!!!\n\n";
                
                print_line_of_markers(")", &error_strm);
                std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
                exit(1);            
            }
            
        } //end if "-average_fragment_length"      
        
        
        
        
        
        
        
        else if (strcmp(option,"-hPt") == 0) // heuristic posterior threshhold
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            heuristic_posterior_PROBABILITY_threshold  =  (longdouble)atof(temp_input);
            
            if (my_MPI_rank == 0)            
	    {  std::cerr<< "set   heuristic_posterior_PROBABILITY_threshold  =  "  <<  heuristic_posterior_PROBABILITY_threshold <<  "\n";  }
                        
        } //end if  "heuristic_posterior_PROBABILITY_threshold"
        
        
        
        
        else if (strcmp(option,"-hCt") == 0) // heuristic posterior threshhold
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            heuristic_posterior_COUNT_threshold  =  (uint)atoi(temp_input);
            
            if (my_MPI_rank == 0)            
	    {  std::cerr<< "set   heuristic_posterior_COUNT_threshold  =  "  <<  heuristic_posterior_COUNT_threshold <<  "\n";  }
                        
        } //end if  "heuristic_posterior_COUNT_threshold"        
        
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-o") == 0)
        {
            char temp_out_dir[option_size];
            strcpy(temp_out_dir, &argv[arg_indx][arg_indx_char]);
            
            output_dir.assign(temp_out_dir);
            
            if ( output_dir.empty() )
            {
                std::cerr << "ERROR:  output_dir not specificed. temp_out_dir   =   [" << temp_out_dir << "]\n\n";
                //show_how_to_use(argv[0]);
                return -1;
            }      
            else
                if (my_MPI_rank == 0)
                    std::fprintf(stderr, "\n\toutput_dir set to [%s]\n", output_dir.c_str() );            
        } //end if "-o"             
        
        
        
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-display_MAP") == 0)
        {            
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            display_MAP_dir.assign( temp_input );
            if (display_MAP_dir.empty())
            {
                std::fprintf(stderr, "ERROR: display_MAP_dir not specificed.\n\n");
                //show_how_to_use(argv[0]);
                return -1;
            }     
            else
                if (my_MPI_rank == 0)
                    std::fprintf(stderr, "\n\tdisplay_MAP_dir set to [%s]\n", display_MAP_dir.c_str() );
        } //end if "-s"          
        
        
        
        
        
        
        else if (strcmp(option,"-keep") == 0) // heuristic posterior threshhold
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            keep_sparse_variational_positions_for_dupdels =    ( atoi(temp_input) == 1 );
            
            if (my_MPI_rank == 0)            
                std::cerr<< "set   keep_sparse_variational_positions_for_dupdels  =  "  <<  keep_sparse_variational_positions_for_dupdels <<  "\n";
                        
        } //end if  "heuristic_posterior_COUNT_threshold"          
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-vs") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            visitation_schedule_subdirname.assign( temp_input );
            if (visitation_schedule_subdirname.empty())
            {
                std::fprintf(stderr, "ERROR: visitation_schedule_subdirname not specificed.\n\n");
                return -1;
            }     
            else
                if (my_MPI_rank == 0)
                    std::fprintf(stderr, "\n\tvisitation_schedule_subdirname set to [%s].\n", visitation_schedule_subdirname.c_str());
        } //end if "-vs"         
        
        
                
        
        
        
        else if ( strcmp(option, "-hM") == 0 ) // heuristic max
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            maximum_number_of_breakpoints_to_consider_after_heuristic  =  (uint)atoi(temp_input);
	    
	    if (maximum_number_of_breakpoints_to_consider_after_heuristic < 2)
	    {
		std::cerr << "\n\nmaximum_number_of_breakpoints_to_consider_after_heuristic = " << maximum_number_of_breakpoints_to_consider_after_heuristic << "   <  2  !!! ILLEGAL because of Gene Conversion requirement!!!\n\n\n";
	    }
            
            if (my_MPI_rank == 0)            
	    {
                std::cerr  <<  "set   maximum_number_of_breakpoints_to_consider_after_heuristic  =  "  
                           <<   maximum_number_of_breakpoints_to_consider_after_heuristic <<  "\n";
	    }
                        
        } //end if  "heuristic max"               
        
        
        
        
        
        else if ( strcmp(option, "-after") == 0 ) // heuristic max
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            last_completed_CC  =  atoi(temp_input);
            
            if (my_MPI_rank == 0)            
                std::cerr  <<  "set   last_completed_CC  =  "  
                           <<   last_completed_CC <<  "\n";
                           
//             assert( last_completed_CC >= 0);
                        
        } //end if  "heuristic max"            
        
        
        
        
        
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-nat_poisson") == 0)
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            if (  atoi(temp_input)  == 1  )
                only_compute_and_save_natural_poisson_intervals = true;
            
            std::cerr<< "\n\nonly_compute_and_save_natural_poisson_intervals  =  " << only_compute_and_save_natural_poisson_intervals << " !!!!!\n\n\n\n";            
            
        } //end if "-average_fragment_length"              
        
        
        
        
        
        
        
        else if (strcmp(option,"-phred") == 0)
        {
	    std::sscanf(&argv[arg_indx][arg_indx_char], "%u", &ASCII_quality_offset);
            
            if (my_MPI_rank == 0)
	    {  std::cerr<< "\n\nASCII_quality_offset  =  " << ASCII_quality_offset << ".\n\n\n";  }
            
        } //end if "-average_fragment_length"               
        
        
        
        
        
        else if ( strcmp(option, "-block") == 0 ) // heuristic max
        {
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            search_block_size  =  (uint)atoi(temp_input);
            
            if (my_MPI_rank == 0)            
                std::cerr  <<  "set   search_block_size  =  "  
                           <<   search_block_size <<  "\n";
                                                   
        } //end if  "heuristic max"            
                
        
        
        
        
        
        
        
        
        
        
        
        else if ( strcmp(option, "-postprocess") == 0 ) // heuristic max
        {
	    std::cerr << "\n\treading post-processing filenames...\n";
	    
	    char in_genes_fname[option_size];
	    char in_vals_fname[option_size];
	    char in_postdatadir[option_size];
	    
	    std::sscanf(&argv[arg_indx][arg_indx_char], "%[^,],%[^,],%s", in_genes_fname, in_vals_fname, in_postdatadir);
	    
	    Genes_in_human_genome__filename.assign(in_genes_fname);
	    biological_validations_filename.assign(in_vals_fname);
	    postdatadir.assign(in_postdatadir);
	    
	    std::cerr << "\t\tGenes_in_human_genome__filename  = [ " << Genes_in_human_genome__filename << "]\n";
	    std::cerr << "\t\t  = [ " << biological_validations_filename << "]\n";
	    
	    do_post_processing_stats_on_calls = true;
                                                   
        } //end if  "postprocess"                    
        
        
        
       
       
        
        
        
        
        
        
        
        else if (strcmp(option,"-selections_dir") == 0)
        {            
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
         
            validations_dirname.assign( temp_input );
            
            if (validations_dirname.empty())        
            {
                if (my_MPI_rank == 0)  
                    std::fprintf(stderr, "ERROR:  failed to read selected connected components from file in option \"-selections_file\".\n\n");                            
                exit(1);
            }      
            else
                if (my_MPI_rank == 0)             
                    std::fprintf(stderr, "\n\tsuccessfully read selected connected components from file [%s]\n",  validations_dirname.c_str() );                  
        } //end if "selections_dir"                 
        
        
        
        
        
        
        
        else if (strcmp(option,"-varpos") == 0)
        {            
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
         
            if (strlen(temp_input) == 0)
            {
                std::cerr << "\n\n\nERROR!   failed to provide varpos directory!!!\n\n\n";
                exit(1);            
            }
            else
            {   
                variational_positions_directory.assign( temp_input );
                
                if (my_MPI_rank == 0)            
                    std::fprintf(stderr, "\n\tvariational_positions_directory set to [%s]\n",  variational_positions_directory.c_str() );      
            }
        } //end if "-s"                 
        
                
        
        
        
        else if (strcmp(option,"-G") == 0)
        {            
            genome_name.assign(  &argv[arg_indx][arg_indx_char]  );
            
         
            if ( !genome_name.empty()   and  my_MPI_rank == 0)            
                    std::fprintf(stderr, "\n\tgenome_name set to [%s]\n",  genome_name.c_str() );                  
        } //end if "-s"                 
        
                        
        
        
        
  
        
        
        
        
        
        else if (strcmp(option,"-gender_and_pop") == 0)
        {            
            gender_and_population_filename.assign( &argv[arg_indx][arg_indx_char] );
            
            if (gender_and_population_filename.empty())
            {
                std::cerr << "\n\n\nERROR!   failed to provide gender and population file!!!!!!\n\n\n";
                exit(1);            
            }
            else if (my_MPI_rank == 0)                           
                std::cerr << "\n\nset gender_and_population_filename = [" << gender_and_population_filename << "]\n\n";                       
                                      
        } //end if "-gender"                        
        
        
        
        
        
        
        
        else if (strcmp(option,"-copy_from_most_recent_of_same_job_name") == 0) // heuristic posterior threshhold
        {            
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
            
            if (strlen(temp_input) == 1  and   ( atoi(temp_input) == 1 )    )
                copy_data_from_most_recent_job_of_same_name = true;            
            
            if (my_MPI_rank == 0)            
                std::cerr<< "set   copy_data_from_most_recent_job_of_same_name  =  "  <<  copy_data_from_most_recent_job_of_same_name <<  "\n";
                        
        } //end if  "copy_from_most_recent_of_same_job_name"          
                
        
        else if (strcmp(option,"-copy_from_specific_job_ID") == 0) // heuristic posterior threshhold
        {            
	    copy_data_from_job_with_jobID.assign(&argv[arg_indx][arg_indx_char]);
	    copy_data_from_most_recent_job_of_same_name = true;      
	    	    
            if (my_MPI_rank == 0)            
                std::cerr<< "set   copy_data_from_job_with_jobID  =  "  <<  copy_data_from_job_with_jobID <<  "\n";	
	                                              
        } //end if  "copy_from_most_recent_of_same_job_name"           
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-reformat_calls") == 0)
        {            
            uint lower_job_ID_range;
            uint upper_job_ID_range;
            char temp_job_outs_dir[option_size];
            char temp_data_dir[option_size];
            
            std::sscanf( &argv[arg_indx][arg_indx_char],  "%[^,],%[^,],%u,%u",
                                                    temp_data_dir,
                                                    temp_job_outs_dir,
                                                    &lower_job_ID_range,
                                                    &upper_job_ID_range  ); 
                                                                                       
	    data_directory.assign(temp_data_dir);
	    
	    
            const BOOST_Interval job_ID_number_range(lower_job_ID_range, upper_job_ID_range); 
            std::cerr << "\n\n\n\nreformat calls job ID range:  " << job_ID_number_range << "\n\n";
            
            
            type_map_uint_to_CC Conn_Comps_for_reformat;
            {//load CCs
            
		char in_connected_components_filename[option_size];
                std::sprintf(in_connected_components_filename, "%s/UID_and_connected_component", data_directory.c_str() );                   
                
                
                read_connected_component_data_from_file_and_create_Event_shells(in_connected_components_filename, Conn_Comps_for_reformat);
                
                for (type_map_uint_to_CC::iterator it_cc = Conn_Comps_for_reformat.begin();
                        it_cc != Conn_Comps_for_reformat.end();
                        ++it_cc)
                    for( type_map_uint_to_Event::iterator it_ev = it_cc->second.events.begin();
                            it_ev != it_cc->second.events.end();
                            ++it_ev)
                        read_Event_basic_data_from_file(it_cc->second, it_ev->second);           
            }//load CCs
            
            
            
            
            const boost::filesystem::directory_iterator it_end;            
            for (boost::filesystem::directory_iterator  it_jobdir(temp_job_outs_dir);
                    it_jobdir != it_end;
                    ++it_jobdir)                    
            {
		if (  !  boost::filesystem::is_directory(  (*it_jobdir).path()  )  )
		    continue;
		
		output_dir.assign(  (*it_jobdir).path().c_str()  );
//                 const std::string job_name_str( (*it_jobdir).path().c_str() );
                const uint pos_of_dot_o = output_dir.find(".o");
//                 const uint pos_of_dot_mgt = job_name_str.find(".mgt");
                
                if (pos_of_dot_o == std::string::npos )//  or  pos_of_dot_mgt == std::string::npos )
                    continue;
                 
                const uint job_id_val = (uint)atoi(   output_dir.substr(pos_of_dot_o+2).c_str()    );
                //(uint)atoi(   job_name_str.substr(pos_of_dot_o+2, pos_of_dot_mgt - pos_of_dot_o - 2 ).c_str()    );
                
                if (  !BOOST_in(job_id_val, job_ID_number_range)  )
                    continue;  
                
                
                std::cerr << "\n\n\nreformatting for job name = [" << output_dir << "]...\n\n";
                
                {//remove
                    boost::filesystem::remove(  std::string(output_dir).append("/Positive_calls")  );
                    boost::filesystem::remove(  std::string(output_dir).append("/Negative_calls")  );
                    boost::filesystem::remove(  std::string(output_dir).append("/All_calls")  );
		    boost::filesystem::remove(  std::string(output_dir).append("/marginal_distributions_via_posterior_sample")  );		    
                }
                
                
                {
                    char temp_genome_name[option_size];
                    const uint pos_last_slash = output_dir.find_last_of("/");
                    std::sscanf(  output_dir.substr(pos_last_slash+1).c_str(), "nf_%[^_]_%*s", temp_genome_name   );
                    genome_name.assign( temp_genome_name );
                    std::cerr << "\n\n\tdetected genome_name = [" << genome_name << "]\n";
                }
                
                
                
                
                                
                for (boost::filesystem::directory_iterator  it_samples(  std::string(output_dir).append("/Samples_from_posterior")   );
                        it_samples != it_end;
                        ++it_samples) 
                {
                    uint in_CC_id;
                    uint in_num_samples;                                        
                    
                    int sscanfreturnval 
			    =  std::sscanf(   (*it_samples).path().filename().c_str(),   "samples_from_posterior_CC_%u__number_of_samples_%u", 
											    &in_CC_id, &in_num_samples);
                    
                    std::cerr << "\n(*it_samples).path().c_str() = [" << (*it_samples).path().c_str() << "]\n\t\t\tin_CC_id = " << in_CC_id << ",   in_num_samples = " << in_num_samples << "\n";
                    
                    if (sscanfreturnval == EOF)
                        std::cerr << "\n\n\n\n\nsscanfreturnval = EOF!!!!\n\n\n";
                    
                    
                    
                    for( type_map_uint_to_Event::iterator it_ev = Conn_Comps_for_reformat.at(in_CC_id).events.begin();
                            it_ev != Conn_Comps_for_reformat.at(in_CC_id).events.end();
                            ++it_ev)
                        read_Event_ALL_data_from_file( Conn_Comps_for_reformat.at(in_CC_id), it_ev->second );
                    
                    
                    const  type_vector_map_uint_to_Sampled_diploid_Event_data  samples_from_posterior(
                                            load_posterior_samples_from_file(output_dir, in_CC_id, in_num_samples)     );
                    
				
					    
		    //Marginals
                    const type_map_uint_to_Marginal_Event_posterior_distribution marginal_distributions
                                            ( Conn_Comps_for_reformat.at(in_CC_id).compute_marginal_probabilities_from_posterior_sample(samples_from_posterior)  );      
					    					    
		    save_marginal_posterior_distributions_to_file(  marginal_distributions, output_dir );
					    
                                                
                                            
		    // marginal MAP calls
                    std::stringstream output_marginal_Positive_calls_strm;
                    std::stringstream output_marginal_Negative_calls_strm;
                    std::stringstream output_marginal_all_calls_strm;
		    
		    output_marginal_Positive_calls_strm.precision(10);
		    output_marginal_Negative_calls_strm.precision(10);
		    output_marginal_all_calls_strm.precision(10);
        
                    make_MAP_calls_from_marginal_posterior_distributions_per_Event
                                                (marginal_distributions,
                                                Conn_Comps_for_reformat.at(in_CC_id),
                                                output_marginal_Positive_calls_strm,
                                                output_marginal_Negative_calls_strm,
                                                output_marginal_all_calls_strm   );                            
                                 
                    save_calls_tables_to_file
                        (output_marginal_Positive_calls_strm,
                            output_marginal_Negative_calls_strm,
                            output_marginal_all_calls_strm,
                         output_dir); 
                            
                            
                            
                    for( type_map_uint_to_Event::iterator it_ev = Conn_Comps_for_reformat.at(in_CC_id).events.begin();
                            it_ev != Conn_Comps_for_reformat.at(in_CC_id).events.end();
                            ++it_ev)       
                        it_ev->second.clear_all_data();
                            
                        
                }//it_samples
                                               
            }//it_jobdir
            
            
            std::cerr << "\n\n\nreformatting complete!!!\n\n\n";
            return 0;                                                                                                                                                         
                                                    
        } //end if  "reformat_calls"          
        
        
        
        
        
        
        
        else if (strcmp(option,"-use_saved_likelihoods") == 0)
        {         
	    
	    std::cerr << "\n\t\tusing saved likelihoods...\n\n";
	    
	    std::string input_dirname__and_prior_values(  &argv[arg_indx][arg_indx_char]  );
	    
	    int begin_exp;
	    int exp_step;
	    int last_exp;
	    
	    char in_old_data_dir[option_size];            	    
            std::sscanf( input_dirname__and_prior_values.c_str() ,  "%[^,],%d,%d,%d",   
								    in_old_data_dir,
								   &begin_exp,  &exp_step,   &last_exp ); 
	    use_previously_computed_likelihoods___old_data_dir.assign(  in_old_data_dir  );
	    
	    
	    for (int exp_val = begin_exp; exp_val >= last_exp;  exp_val += exp_step)
		list_of_prior_probability_exponents_for_computing.push_back( exp_val );	    
	    
	    if (my_MPI_rank == 0)
		print_list<int>(list_of_prior_probability_exponents_for_computing, "list_of_prior_probability_exponents_for_computing", NULL, true);
	                                                  	    	    
	    
	} // use_saved_likelihoods
	
	
	
	
	
	
        else if (strcmp(option,"-pseudo_RD_display") == 0)
        {         	    
	    std::cerr << "\n\t\tpseudo_RD_display  option...\n\n";
	        
	    uint in_cc;
	    uint in_ev_uid;
	    char in_outcome_hap_0[option_size];
	    char in_outcome_hap_1[option_size];
	    uint in_profile_brkpt_00, in_profile_brkpt_01;
	    uint in_profile_brkpt_10, in_profile_brkpt_11;
	    
            std::sscanf( &argv[arg_indx][arg_indx_char],  "%u,%u,%[^,],[%u,%u],%[^,],[%u,%u]",
						    &in_cc,
                                                    &in_ev_uid,
                                                    in_outcome_hap_0,
                                                    &in_profile_brkpt_00, &in_profile_brkpt_01,
						    in_outcome_hap_1,
						    &in_profile_brkpt_10, &in_profile_brkpt_11); 	    
        	    
	    pseudo_RD_hypothesis_test_work.push_back(
		    Sampled_diploid_Event_data(
			    in_cc, in_ev_uid,
			    type_haploid_outcome__haploid_outcome(convert_string_to_haploid_outcome(in_outcome_hap_0),
								  convert_string_to_haploid_outcome(in_outcome_hap_1)),
			    type_BI__BI(BOOST_Interval(in_profile_brkpt_00,in_profile_brkpt_01),
					BOOST_Interval(in_profile_brkpt_10,in_profile_brkpt_11) )));
	    

	    targeted_CC_IDs.insert(in_cc);
	    
	    
	    
	    std::cerr << "\n\n\t\t\tread in:"
			<< "\n\t\t\t\tin_cc = " << in_cc
			<< "\n\t\t\t\tin_ev_uid = " << in_ev_uid
			<< "\n\t\t\t\tin_outcome_hap = [" << in_outcome_hap_0 << "], [" << in_outcome_hap_1 << "]"
			<< "\n\t\t\t\tin_profile_brkpt = " << in_profile_brkpt_00 << ", " << in_profile_brkpt_01 
							<< ", " << in_profile_brkpt_10 << ", " << in_profile_brkpt_11
			<< "\n\n";
  
	} // pseudo_RD_display
	
	
	
	
	
	
	
	
	
        else if (strcmp(option,"-unique_alignments") == 0)
        {         	    
	    std::cerr << "\n\t\tshowing alignments of PERs in \"unique\" region...\n\n";
	        	 	    
	    uint in_chro;
	    uint in_pos;
	    uint in_padd;
	                
            std::sscanf( &argv[arg_indx][arg_indx_char],  "%u,%u,%u",
							    &in_chro,
							    &in_pos,						    
							    &in_padd); 
						    							    
	    std::cerr << "\n\n\t\t\tread in:"
			<< "\n\t\t\t\tchromo_of_region = " << in_chro
			<< "\n\t\t\t\tpos_of_interest = " << in_pos
			<< "\n\t\t\t\tin_padd = " << in_padd
			<< "\n\n";						
	    
	    unique_regions_for_PER_displays.push_back(   type_uint__uint__uint(type_uint__uint(in_chro, in_pos), in_padd)   );
			

	    
// 	    return 0;	    			    
	} // unique_alignments
     
     
     
        
        
        
        
        
	
        else if (strcmp(option,"-informativeness") == 0)
        {         	    
	    std::cerr << "\n\t\tdetermining \"informativeness\" of all breakpoints...\n\n";	     

	    uint in_frag_start, in_frag_jump, in_frag_end;
	    uint in_mate_start, in_mate_jump, in_mate_end;
	    std::sscanf( &argv[arg_indx][arg_indx_char], "%u:%u:%u,%u:%u:%u",
								    &in_frag_start,
								    &in_frag_jump,
								    &in_frag_end,
								    &in_mate_start,
								    &in_mate_jump,
								    &in_mate_end   );
								    
	    std::cerr << "\n\n\n\nWill determine informativeness of paired-end read combinations...\n\n\n\n";
	    for (uint frag=in_frag_start; frag <= in_frag_end; frag +=in_frag_jump)
	    {
		for (uint mate=in_mate_start; mate <= in_mate_end; mate +=in_mate_jump)
		{
		    informativeness_checking_parameters.push_back(type_uint__uint(frag, mate));
		}//mate
	    }//frag	    
	} // informativeness        
        
        
        
        
        
        
        
        else if (strcmp(option,"-aCGH") == 0)
        {         	    
	    std::cerr << "\n\t\taCGH...\n\n";	     
	     	    	  
	    char outdir[option_size];
	    char in_probe_intensity_fname[option_size];
	    uint in_chr_char;
	    uint in_brkpt_begin;
	    uint in_brkpt_end;
	    uint in_padding;
	    char in_acceptable_genomes_filename[option_size];
// 	    char in_save_prefix[option_size];
	    char in_individuals_with_call[option_size];

	    std::sscanf( &argv[arg_indx][arg_indx_char], "%[^,],%[^,],%[^,],%u:[%u,%u],%u,%s",
			 outdir,
			 in_probe_intensity_fname,
			 in_acceptable_genomes_filename,
			 &in_chr_char,
			 &in_brkpt_begin, &in_brkpt_end,
			 &in_padding,			 
			 in_individuals_with_call);
			 			 
	    type_set_string individuals_with_call_set;	    	    
	    const std::string in_calls_as_str(in_individuals_with_call);
	    
	    std::cerr << "\t\tsplitting...";
	    boost::split(individuals_with_call_set, in_calls_as_str, boost::is_any_of(","), boost::algorithm::token_compress_on);
	    
	    std::cerr << "done.\n";
	    
	    print_set<std::string>(individuals_with_call_set, "individuals_with_call_set", NULL, true);
	    
	    	    
	    acgh_data_checks.push_back( aCGH_validation_container(
							in_probe_intensity_fname,
							in_chr_char,
							BOOST_Interval(in_brkpt_begin, in_brkpt_end),
							in_padding,
							individuals_with_call_set,
							in_acceptable_genomes_filename,
							outdir ));
	}//acgh                
        
        
        
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-GeneConv_in_heuristic") == 0)
        {         	    
	    
	    int in_val;
	    std::sscanf(&argv[arg_indx][arg_indx_char], "%d", &in_val);
	    
	    if (in_val == 1)
	    {	    
		std::cerr << "\n\t\tWILL compute Gene Converison breakpoints in Heuristic.\nWARNING - compute time will increase dramatically.\n\n\n";
		consider_GeneConversion_breakpoints_in_Heuristic = true;
	    }
	    else if (in_val == 0)
	    {	    
		std::cerr << "\n\t\tIGNORING Gene Converison breakpoints in Heuristic.\n\n\n";
		consider_GeneConversion_breakpoints_in_Heuristic = false;
	    }		
	    else 
	    {
		std::cerr << "\n\nunrecognized option = [" << &argv[arg_indx][arg_indx_char] << "] in option \"GeneConv_in_heuristic\".\n\n\n";		
	    }	    
	} // GeneConv_in_heuristic           
        
        
        
        
        else if (strcmp(option,"-GeneConv_in_sample_space") == 0)
        {         	    
	    
	    int in_val;
	    std::sscanf(&argv[arg_indx][arg_indx_char], "%d", &in_val);
	    
	    if (in_val == 1)
	    {	    
		std::cerr << "\n\t\tWILL compute Gene Converison breakpoints in Sample space.\nWARNING - compute time will increase!!!\n\n\n";
		consider_GeneConversion_in_Full_compute = true;
	    }
	    else if (in_val == 0)
	    {	    
		std::cerr << "\n\t\tIGNORING Gene Converison breakpoints in Sample Space.\n\n\n";
		consider_GeneConversion_in_Full_compute = false;
	    }		
	    else 
	    {
		std::cerr << "\n\nunrecognized option = [" << &argv[arg_indx][arg_indx_char] << "] in option \"GeneConv_in_sample_space\".\n\n\n";		
	    }	    	    	    
	} // GeneConv_in_sample_space            
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-redo_calls") == 0)
        {         	    
	    
	    redo_brkpts_dir.assign(&argv[arg_indx][arg_indx_char]);
	    
	    std::cerr << "\n\n-redo_calls:  [" << redo_brkpts_dir << "]\n\n";
		} // redo_calls             
        
        
        
        else if (strcmp(option,"-redo_RD") == 0)
        {         	    	    
	    redo_RD_dir.assign(&argv[arg_indx][arg_indx_char]);
	    
	    std::cerr << "\n\n-redo_RD:  [" << redo_RD_dir << "]\n\n";
		} // redo_RD          
        
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-raw_count_ratio") == 0)
        {         	 
		if (strcasecmp(&argv[arg_indx][arg_indx_char],"A") == 0)
		{  do_raw_count_ratios__strategy_A = true;  }
		else if (strcasecmp(&argv[arg_indx][arg_indx_char],"B") == 0)
		{  do_raw_count_ratios__strategy_B = true;  }
			    
	    std::cerr << "\n\n-raw_count_ratio:  A = " << do_raw_count_ratios__strategy_A << ",   B = " << do_raw_count_ratios__strategy_B << "\n\n";
		} // GeneConv_in_sample_space                
        
        
        
        
        
        
        
        else if (strcmp(option,"-simulate") == 0)
        {         	 
	    
			char in_num_sim_gen[option_size];
			char in_base_sim_fname[option_size];
			char in_num_ev_per_gen[option_size];
			char in_allowable_fname[option_size];
			
			std::sscanf(&argv[arg_indx][arg_indx_char], "%[^,],%[^,],%[^,],%s",
												in_num_sim_gen, in_num_ev_per_gen, in_base_sim_fname, in_allowable_fname);	
			
			num_sim_genomes = atoi(in_num_sim_gen);
			num_ev_per_sim_genome = atoi(in_num_ev_per_gen);
			simulated_genome_base_fname = in_base_sim_fname;
			allowable_events_fname = in_allowable_fname;
			
			std::cerr << "\n\n\tnum_sim_genomes = " << num_sim_genomes 
				<< "\n\tnum_ev_per_sim_genome = " << num_ev_per_sim_genome
				<< "\n\tsimulated_genome_base_fname = [" << simulated_genome_base_fname << "]"
				<< "\n\tallowable_events_fname = [" <<allowable_events_fname << "]\n";
				
			make_simulated_genomes = true;
		} // GeneConv_in_sample_space                   
        
        
        
        
        
        
        
        
        
        
        
        
        else if (strcmp(option,"-debug_view_PSF") == 0)
        {            
     
            strcpy(temp_input, &argv[arg_indx][arg_indx_char]);                        
            
            char temp_view[option_size];
            int is_a_dir, is_a_PSF_UUID;
            std::sscanf(temp_input, "%d,%d,%s", &is_a_dir, &is_a_PSF_UUID,  temp_view);      
            
            std::fprintf(stderr, "\n\tdebug_view [%s],     is_a_dir = %d     is_a_PSF_UUID = %d\n", temp_view, is_a_dir,  is_a_PSF_UUID); 
                        
//             Partial_sum_function::debug_view( std::string(temp_view), (bool)is_a_dir, (bool)is_a_PSF_UUID );
            
            
            
            used_debug_view = true;
        } //end if "-Reference genome"             
        
        
        
        
        
        
//         ///DEBUG
//         else if (strcmp(option,"-test_PER_aligner") == 0)
//         {                       
//             strcpy(temp_input, &argv[arg_indx][arg_indx_char]);
//             
//             char ref_seq_char[10*option_size];
//             char mate_A_char[10*option_size];
//             char mate_B_char[10*option_size];
// 
//             char mate_A_qual[10*option_size];
//             char mate_B_qual[10*option_size];            
//             
//             std::sscanf(temp_input, "%[^,],%[^,],%[^~]~%[^,],%s", ref_seq_char, mate_A_char, mate_A_qual, mate_B_char, mate_B_qual );                                             
//             
//             PER_aligner my_PER_aligner (  -15,   800,    NULL);
//                                 
//             std::string ref(ref_seq_char);
// //             std::string mate_A(mate_A_char);
// //             std::string mate_B(mate_B_char);
//             
//             BamTools::BamAlignment mate_A, mate_B;
//             
//             mate_A.QueryBases = std::string(mate_A_char);                        
//             mate_B.QueryBases = std::string(mate_B_char);            
//             
//             mate_A.SetIsReverseStrand(false);
//             mate_B.SetIsReverseStrand( ! mate_A.IsReverseStrand() ); 
//             
//             
//             
//             mate_A.Qualities = std::string(mate_A_qual);  //std::string(  mate_A.QueryBases.size(),  '5' );
//             mate_B.Qualities = std::string(mate_B_qual); //std::string(  mate_B.QueryBases.size(),  '5' );
//             
// //             mate_A.Qualities = ";DD>>@:BD>A?>BD>43>#6=:B::4,886=*4\"3";  //"6.9:;0==:=<>17>8?=<=?<@48@<9A>>=@>;:"; // "?$C@@BBBDCCB=:?$>69<3E@?;@F?7C>E=D>=$B4@:=>BB5DCC?D=A8DF?E=CDDGEE@EEGG?DH9EB";
// //             mate_B.Qualities = "C9B>>DB@A;@@@AA<A@?AAED@DBDA?;ADADC@"; //"@DDDCBA?DDDDDBCAD?CE9:?;ADAAE@1D<<*;";  //"=@B=A?ADBD6BEABC:FFD@BACEEBBDGBB?@:FEEFD;;4FDAECC?DB@D??@CBBAC?BA";
// //             
// //             mate_A.Qualities = mate_A.Qualities.substr(0, mate_A.QueryBases.size());
// //             mate_B.Qualities = mate_B.Qualities.substr(0, mate_B.QueryBases.size());
// //             mate_A.Qualities = std::string(mate_A.QueryBases.size(), '-');
// //             mate_B.Qualities = std::string(mate_B.QueryBases.size(), '-');                        
//             ASCII_quality_offset = 33;
//             
//             
//             double time_begin = omp_get_wtime();
//             
//             type_string__string  PER_and_hybrid_algmt;
//             BOOST_Interval  subset_of_region_that_is_aligned;
//             type_BI__BI   subset_of_region_that_is_aligned_to_mates__A__B;
//             
//             const real gen_prob( my_PER_aligner.calculate_marginal_probability_that_PER_was_generated_from_some_subregion_of_reference
//                                 (  mate_A,// std::string in_read_A,
//                                   mate_B,// std::string in_read_B,
//                                   ref,
//                                   true,
//                                   false,
//                                  PER_and_hybrid_algmt,
//                                  subset_of_region_that_is_aligned,
//                                  subset_of_region_that_is_aligned_to_mates__A__B)       );
//                               
//                                  
//             std::cerr << "\n\ntime for alignment SEED:  " <<    omp_get_wtime() - time_begin << "\n";                   
//                                 
// //                                 std::fprintf(stderr, "\n made it !\n");                    
//             std::stringstream out_ss;
//             
//             out_ss << "\nref        size " << ref.size() <<":\n" << ref << "\n\n";
//             out_ss << "\nmate_A        size " << mate_A.QueryBases.size() <<":\n" << mate_A.QueryBases << "\n\n";
//             out_ss << "\nmate_B        size " << mate_B.QueryBases.size() <<":\n" << mate_B.QueryBases << "\n\n";
//             out_ss << "\n\nalignment:\n" << PER_and_hybrid_algmt.first << "\n" << PER_and_hybrid_algmt.second << "\n\n";
//             
//             out_ss<< "\n\nprobability = " << gen_prob.toString(15);
//             out_ss<< "\n\nsubset_of_region_that_is_aligned = " << subset_of_region_that_is_aligned
//                     << "\nsubset_of_region_that_is_aligned_to_mates__A__B = " << subset_of_region_that_is_aligned_to_mates__A__B.first
//                     << "\n\t\t\t" << subset_of_region_that_is_aligned_to_mates__A__B << "\n\n";
//             
//         
//             
//             std::fprintf(stderr, "\n\n%s\n\n", out_ss.str().c_str() );            
//             
//             return 0;               
//         } //end if "-test_PER_aligner"     
//         
//         
//         
//         
//         
        
        
        
        else if (strcmp(option,"-test_linker_and_packages") == 0)
        {         	
	    std::cerr << "\n\n\n\n\nTESTING EXECUTABLE!!!!\n\n\n";
	    
	    std::cerr << "\n\ttesting mpreal:";
	    {
		real dummy_val(71.75);
		std::cerr << "\n\t\t\t" << dummy_val.toString(10);
		std::cerr << "\n\t\t\t" << (dummy_val/real(31.6)).toString(15);
	    }//real
	    
	    std::cerr << "\n\ttesting gmpfrxx:";
	    {
		mpfr_class dummy_val(71.75);
		dummy_val.get_mpfr_t();
	    }//gmpfrxx
	    
	    std::cerr << "\n\ttesting boost statistics:";
	    {
		mpfr_class parameter_BOOST_r(23);
		mpfr_class parameter_BOOST_p(0.65);
		boost::math::negative_binomial_distribution<mpfr_class> my_negative_binomial(parameter_BOOST_r, parameter_BOOST_p);
                
		const mpfr_class observed_value(13);
		const real result_pdf( boost::math::pdf<mpfr_class>(my_negative_binomial, observed_value).get_mpfr_t()  );     		
		std::cerr << "\n\t\tmy_negative_binomial = " << result_pdf.toString(15);		
	    }//boost stats
	    
	    std::cerr << "\n\ttesting BAMtools:";
	    {
		BamTools::BamAlignment balgmt;
		BamTools::BamReader myreader;
	    }//BAMtools
	    	    
	    std::cerr << "\n\n\nDONE WITH EXECUTABLE TESTING.\n\n\n";
	    return 0;
	} // test_linker_and_packages           
        
	
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                    
   
   
        delete[] temp_input;
        
    } // end of for-loop: arg_indx 

    delete[] option;  

    // END OF OPTIONS !!!!
                           
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if (did_a_compare_to_vals)
    {  return 0;  }
    
    if (used_debug_view)
    {  return 0;  }
    
    
    
    
    
    
    
    
    
    

    if (my_MPI_rank == 0)
    {  std::fprintf(stderr, "\n\nread_chromosome_data_from_file\n\n" );      }    
    
    read_chromosome_data_from_file();    
    
    
    
    
    if (length_of_haploid_genome == 0)   
    {
        for (type_map_uint_to_uint::const_iterator it_chro = Event::chromosome_lengths.begin();
                it_chro != Event::chromosome_lengths.end();
                ++it_chro)
	{  length_of_haploid_genome += it_chro->second;    }
    }

    if (my_MPI_rank == 0)
    {  std::fprintf(stderr, "\n\nlength_of_haploid_genome = %u\n\n", length_of_haploid_genome);  }    
    
    
    
    
    
    


    
    
    if (!acgh_data_checks.empty())
    {
	std::cerr << "\n\n\nperforming  acgh_data_checks...\n\n\n";
	
	print_map_keys_and_values<uint,uint>(Event::chromosome_lengths, "Event::chromosome_lengths");
	
	for (type_list_aCGH_validation_container::iterator it_acgh = acgh_data_checks.begin();
		it_acgh != acgh_data_checks.end();
		++it_acgh)
	{
	    char addl_suffix[option_size];
	    std::sprintf(addl_suffix, "/aCGH_validation__chr_%u___%u__%u",
			 it_acgh->chr_of_call,
			 it_acgh->brkpt_of_call.lower(),
			 it_acgh->brkpt_of_call.upper() );	    
	    
	    it_acgh->save_prefix.append(addl_suffix);
	    	    	
	    validate_using_probe_intensity_comparison( *it_acgh );
	}    
	
	std::cerr << "\n\n\nDONE with  acgh_data_checks.\nExiting program...\n\n\n";
	return 0;    
    }//acgh
    
    
    
    

    

    

    

                    
    type_map_uint_to_CC Conn_Comps;
    if (!data_directory.empty())
    {//scope    
	char in_connected_components_filename[option_size];
	std::sprintf(in_connected_components_filename, "%s/UID_and_connected_component", data_directory.c_str());   
	
	read_connected_component_data_from_file_and_create_Event_shells
						    (in_connected_components_filename,
						    Conn_Comps);                                                

	std::fprintf(stderr, "DEBUG: read_connected_component_data_from_file_and_create_Event_shells ...done!\n\n");
	    

	std::fprintf(stderr, "DEBUG: eliminating undesirable Conn Comps.\n");   
	
	type_map_uint_to_uint counts;
	for (type_map_uint_to_CC::const_iterator it_cc = Conn_Comps.begin(); it_cc != Conn_Comps.end(); ++it_cc)
	{  counts[it_cc->second.events.size()]++;  }
	
	print_map_keys_and_values<uint,uint>(counts, "num CCs with x events", NULL, true);
	
	
	
	{
	    type_set_uint desired_Events;
	    for (type_list_Sampled_diploid_Event_data::const_iterator it_disp = display_Event_and_breakpoint_range.begin();
		    it_disp != display_Event_and_breakpoint_range.end();
		    ++it_disp)
	    {  desired_Events.insert(it_disp->event_UID);  }
	
	    remove_undesirable_Connected_Components_from_set_of_Connected_Components
				(Conn_Comps,
				targeted_CC_IDs,
				desired_Events,
				acceptable_Connected_Component_size_range,
				acceptable_Connected_Component_compute_range,
				last_completed_CC);
	}
    
			    
	//degenerate CCs:    
	for (type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{
	    for (type_map_uint_to_Event::iterator it_ev = it_cc->second.events.begin();
		    it_ev != it_cc->second.events.end();
		    ++it_ev)          
	    {  read_Event_basic_data_from_file(it_cc->second, it_ev->second);  }
	}//cc			 
	
	remove_degenerate_Events(Conn_Comps);
                							    
    }//scope                                                            
    
    
    
    
    
    
//     std::cerr << "\n\nConnected Component\tEvent\tchromosome\tLCR_A_begin\tLCR_A_end\tLCR_B_begin\tLCR_B_end\trecomb_type\n\n";
//     std::ofstream  outfs("basic_event_data.txt");
//         for (type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
//                 it_cc != Conn_Comps.end();
//                 ++it_cc)
//         {            
//             
//             for (type_map_uint_to_Event::iterator it_ev = it_cc->second.events.begin();
//                     it_ev != it_cc->second.events.end();
//                     ++it_ev)
//             {
//                 read_Event_ALL_data_from_file( it_cc->second, it_ev->second );
//                 outfs << it_cc->first 
//                             << "\t" << it_ev->first 
//                             << "\t" << it_ev->second.chromos[0]
//                             << "\t" << it_ev->second.LCRs[0].lower() << "\t" << it_ev->second.LCRs[0].upper()
//                             << "\t" << it_ev->second.LCRs[1].lower() << "\t" << it_ev->second.LCRs[1].upper()
//                             << "\t";
// //                 if (it_ev->second.recomb_type == 1)
// //                     std::cerr << "inv";
// //                 else
// //                     std::cerr << "dup/del";
//                 
//                 outfs << "\n";  
//                 
// //                 it_ev->second.clear_all_data();
//             }
//             
// //             std::cerr <<"\n\n";
//         }
//     
//     
//     outfs.close();
//     exit(1);






        
        
        
    
    
    
	if (do_raw_count_ratios__strategy_A or  do_raw_count_ratios__strategy_B)
	{
		if (my_MPI_rank == 0)
		{  boost::filesystem::create_directory(scratch_job_output_dir);  }
		std::cerr << "\nload_cumulative_GC_counts_for_each_chromosome...\n";
		load_cumulative_GC_counts_for_each_chromosome();
	}    
    
	if (do_raw_count_ratios__strategy_A)
	{
// 		std::ofstream outfs("/gpfs/scratch/mmparks/SAMs_and_BAMs/coverage_stats.txt");
		
		//loop through bam files
		const boost::filesystem::directory_iterator it_end;  
		for (boost::filesystem::directory_iterator it_bamdir("/users/mmparks/data/mmparks/simulate");  // simulation change
			it_bamdir != it_end;
			++it_bamdir)
		{		
			if (!boost::filesystem::is_directory(it_bamdir->path()))
			{
				if (it_bamdir->path().filename().string().find(".bas") == std::string::npos
					and  it_bamdir->path().filename().string().find(".bai") == std::string::npos
					and  it_bamdir->path().filename().string().find(".bam") == it_bamdir->path().filename().string().size() - 4
					and  it_bamdir->path().filename().string().find("TCGA") == std::string::npos
					and  longdouble(boost::filesystem::file_size(it_bamdir->path()))/1073741824  >= 20)  // GB
				{
					BAM_filename = boost::filesystem::complete(it_bamdir->path()).string();
					genome_name = infer_genome_name_from_BAM_filename(BAM_filename);
					std::cerr << "\n\n\nFor genome_name = [" << genome_name << "]\n\n";
					
					prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population(BAM_filename, gender_and_population_filename, my_BAM_reader);
					
// 					uint total_number_proper_pairs = 0;
// 					longdouble total_coverage = 0;
// 					
// 					for (type_map_string_to_Readgroup_stats::const_iterator it_rg = Readgroup_stats.begin();
// 						it_rg != Readgroup_stats.end();
// 						++it_rg)
// 					{
// 						total_number_proper_pairs += it_rg->second.number_mapped_Proper_Pairs_reads/2;
// 						total_coverage +=  ((longdouble)it_rg->second.number_mapped_Proper_Pairs_reads/2) * it_rg->second.median_insert_size / length_of_haploid_genome;
// 					}
// 					
// 					outfs << genome_name << "\t" << total_number_proper_pairs << "\t" << total_coverage << "\n";
					
					
 					draw_count_ratio_using_sampling_strategy_A(Conn_Comps, my_BAM_reader);
				}		    		    
			}//file
		}//it_jobdir
		
// 		outfs.close();
	}
	
	if (do_raw_count_ratios__strategy_B)
	{
// 		draw_count_ratio_using_sampling_strategy_B(Conn_Comps, my_BAM_reader);
	}
	
	if (do_raw_count_ratios__strategy_A or  do_raw_count_ratios__strategy_B)
	{
		return 0;
	}
	
    
    
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
	if (make_simulated_genomes)
	{
		std::cerr << "\nmake_simulated_genomes...\n";
		std::cerr << "Conn_Comps.size = " << Conn_Comps.size() << "\n";
		for (type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
			it_cc != Conn_Comps.end();
			++it_cc)
		{  read_Connected_Component_ALL_Event_data_from_file(it_cc->second);  }   		
		
		type_vector_uint allowable_events;			
		{
			std::ifstream infs(allowable_events_fname);
			uint in_allow_ev;
			
			while (infs >> in_allow_ev)
			{  allowable_events.push_back(in_allow_ev);  }
			
			infs.close();
		}
		
		print_vector<uint>(allowable_events, "allowable_events", NULL, true);

		for (uint sim_ix = 0; sim_ix < num_sim_genomes; ++sim_ix)
		{
			std::cerr <<  "\n\n\n\nsim_ix = " << sim_ix << "\n";
			char sim_gen_name[option_size];
			std::sprintf(sim_gen_name, "%s.%u.fa", simulated_genome_base_fname.c_str(), sim_ix);		
			
			simulate_NAHR_events_on_Reference_genome(Conn_Comps, sim_gen_name, allowable_events, num_ev_per_sim_genome);
		}
		
		std::cerr << "DONE creating all simulation genomes!\n";
		exit(1);
		
		
	}//make_simulated_genomes
        
        
        
    
    
    
    
    
    
    
    
    if (!redo_brkpts_dir.empty())
    {
	for (type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{  read_Connected_Component_ALL_Event_data_from_file(it_cc->second);  }
	
	std::cerr <<"\n\n\n\nreformatting calls...\n";
	redo_all_marginal_MAP_calls_using_saved_samples_from_posterior(redo_brkpts_dir, Conn_Comps);
	
	std::cerr << "\n\nDONE reformatting calls!\n";
	return 0;	
    }//redo_brkpts_dir
    
    
    
    
    
    if (!redo_RD_dir.empty())
    {
		for (type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
			it_cc != Conn_Comps.end();
			++it_cc)
		{  read_Connected_Component_ALL_Event_data_from_file(it_cc->second);  }
		
		std::cerr <<"\n\n\n\nredoing RD...\n";
	// 	redo_all_inferred_RD_tests_and_make_ratios(redo_RD_dir, Conn_Comps);
	// 	do_signed_rank_test_for_Read_Depth(redo_RD_dir, Conn_Comps);
// 		bernoulli_lyapunov_nonparametric_test(redo_RD_dir, Conn_Comps);
 		raw_count_ratios(redo_RD_dir, Conn_Comps);
// 		sample_control_positions_for_all_genomes(redo_RD_dir);
		
//		determine_cumulative_GC_count_for_each_chromosome_and_save_to_file();
		
		std::cerr << "\nrank " << my_MPI_rank << "  DONE redoing RD!\n";
		
		std::cerr << "\nrank " << my_MPI_rank << " waiting at barrier\n";
		world_MPI_boost_ptr->barrier();
		return 0;	
    }//redo_RD_dir    
    
    
  
    
    
    
    
    

    if (  genome_name.empty()  and !only_compute_and_save_natural_poisson_intervals   and  informativeness_checking_parameters.empty()   and !do_post_processing_stats_on_calls)
    {
        if (my_MPI_rank == 0)
	{  std::cerr << "\n\n\ninferring genome name from input BAM file...\n";  }
        
        genome_name = infer_genome_name_from_BAM_filename(BAM_filename);
        
        const int last_slash = BAM_filename.find_last_of('/');
        const int dot_at_end_of_genome_name = BAM_filename.find_first_of( '.', last_slash); 
        if (genome_name.empty())
        {
		std::cerr << "unable to infer genome name rom this BAM_filename:  [" << BAM_filename << "]   !!!\n\n\n\n";
		exit(1);        
        }
        
        if (my_MPI_rank == 0)
	{  std::cerr << "\t\tinferring genome name to be:  " << genome_name << "\n\n\n\n";  }
    }
    
    
    
    
    
    
    
    
    
    
//     {//print some stats
//     
//     
//         
//         
//         type_map_uint_to_CC Conn_Comps;
//         {//scope    
//                     char in_connected_components_filename[option_size];
//                     std::sprintf(in_connected_components_filename, "%s/UID_and_connected_component", data_directory);   
//                     
//                     read_connected_component_data_from_file_and_create_Event_shells
//                                                                 (in_connected_components_filename,
//                                                                 Conn_Comps);                                                                                
//         }//scope
//         
//         if (my_MPI_rank == 0)
//             std::fprintf(stderr, "DEBUG: read_connected_component_data_from_file_and_create_Event_shells ...done!\n\n");    
//         
//         
//         
//         for (type_map_uint_to_CC::iterator cc_it = Conn_Comps.begin();
//                 cc_it != Conn_Comps.end();
//                 ++cc_it)
//         //  upload Event Data:
//             for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
//                     it_ev != cc_it->second.events.end();
//                     ++it_ev)          
//                 read_Event_ALL_data_from_file(cc_it->second, it_ev->second);
//             
//                         
//         if (my_MPI_rank == 0)
//             std::fprintf(stderr, "Successfully read in all event data.");         
// 
//     
//     
//         uint max_CC_size = 0;
//     
//          for (type_map_uint_to_CC::const_iterator it_cc =  Conn_Comps.begin();
//                 it_cc != Conn_Comps.end();
//                 ++it_cc)
//             if (it_cc->second.events.size() > max_CC_size)
//                 max_CC_size = it_cc->second.events.size();
//             
//             
//         std::stringstream out_ss;    
//             
//         for (uint cc_size = 1; cc_size <= max_CC_size; ++cc_size)        
//             for (type_map_uint_to_CC::const_iterator it_cc =  Conn_Comps.begin();
//                     it_cc != Conn_Comps.end();
//                     ++it_cc)    
//                 if ( it_cc->second.events.size() == cc_size )
//                 {                    
//                     out_ss << "\n\n\nCC " << it_cc->first << "\tsize = " << it_cc->second.events.size() << "   types:   ";
//                     
//                     for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
//                             it_ev != it_cc->second.events.end();
//                             ++it_ev)   
//                         out_ss << it_ev->second.recomb_type << ", ";
//                     
//                     
//                     out_ss << "\tchromos:   ";
//                     
//                     type_set_uint chromos_with_homol;                    
//                     for (type_map_uint_to_Event::const_iterator it_ev = it_cc->second.events.begin();
//                             it_ev != it_cc->second.events.end();
//                             ++it_ev)
//                         for (type_list_Coords_and_homologous_region::const_iterator it_homol_reg 
//                                                 = it_ev->second.regions_homologous_to_directly_affected_region.begin();
//                                 it_homol_reg != it_ev->second.regions_homologous_to_directly_affected_region.end();
//                                 ++it_homol_reg)
//                             chromos_with_homol.insert( it_homol_reg->chromosome_of_homologous_region );
//                         
//                     for (type_set_uint::const_iterator it_chr = chromos_with_homol.begin();
//                             it_chr != chromos_with_homol.end();
//                             ++it_chr)    
//                         out_ss << *it_chr << ", ";                                        
//                     
//                 }
//            
//            
//         std::fprintf(stderr, "\n\n%s\n\n", out_ss.str().c_str() );
//         exit(1);
//     
//     }//print some stats    
    
    
    
    


    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
    
    
    
    
    
    
    
    
    if (!validations_dirname.empty())
    {                                                                    
        const bool success_read_selection = read_selected_connected_components_file(  validations_dirname,  targeted_CC_IDs  );    
        
        if ( !success_read_selection )
        {
            std::cerr << "\n\n\nERROR!  failed to read selected connected components from file [" << validations_dirname << "] !!!\n\n\n";
            exit(1);                    
        }
        else if (my_MPI_rank == 0) 
        {
            std::cerr << "\n\nsuccessfully read selected connected components from file.\n"
                        << "\t\trunning on connedted components...\n\n";
                        
            print_set<uint>(targeted_CC_IDs, "targeted_CC_IDs", NULL, true);
            std::cerr << "\n\n\n";
        }
    }    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if ( !do_post_processing_stats_on_calls 
	and  !only_compute_and_save_natural_poisson_intervals 
	and  informativeness_checking_parameters.empty() )
    {
	std::cerr << "\npreparing BAM reads...\n";
	prepare_BAM_files_for_IO_operations__and__load_BAS_data___and__read_gender_and_population(BAM_filename, gender_and_population_filename, my_BAM_reader);                                   
    }
    
    
    
    

    
    
    
    
    
    


//     if (   do_post_processing_stats_on_calls
// 	    and !only_compute_and_save_natural_poisson_intervals
// 	   and  informativeness_checking_parameters.empty()  )
    {

        {//create directories    
            std::fprintf(stderr, "creating directories via boost:::filesystem():\n");                                	    
	           
            if ( output_dir.empty() )  // was not set in options
            {
                char temp_out_dir[option_size];
                std::sprintf(temp_out_dir,"%s/output",
                                        scratch_job_output_dir.c_str());
                output_dir.assign(temp_out_dir);
            }     
                        
            
            create_scratch_directory_tree();            
            create_output_directory_tree();
                
        }//create directories      
    


        world_MPI_boost_ptr->barrier();    
    
    } //      !do_post_processing_stats_on_calls  or  !only_compute_and_save_natural_poisson_intervals
   
    
    
    
    
    

    
    
    if (!unique_regions_for_PER_displays.empty())
    {
	
	std::cerr << "\n\n\n display_MAP_PER_space_in_some_unique_region_of_the_genome ALL !!!...\n\n";
	
	upload_universal_variational_positions();
	upload_universal_indel_positions();
	
	
	std::cerr << "\n\nBROKEN!!\n\n\n"; 
	exit(1);
// 	for (type_list_uint__uint__uint::const_iterator it_uniq = unique_regions_for_PER_displays.begin();
// 		it_uniq != unique_regions_for_PER_displays.end();
// 		++it_uniq)	    
// 	    display_MAP_PER_space_in_some_unique_region_of_the_genome(it_uniq->first.first, it_uniq->first.second, it_uniq->second);
    
	std::cerr << "\n\n\n  DONE WITH ALL  display_MAP_PER_space_in_some_unique_region_of_the_genome!!!!!!\n\n";
	return 0;
    }
    
    
    
    
    
    
    
    
    if (!pseudo_RD_hypothesis_test_work.empty())
    {
// 	std::cerr << "\n\n\npseudo_RD_exams !!!!\n\n\n";
// 	type_map_uint_to_CC all_CCs_for_repeat_seq;
// 		{
// 		    char in_connected_components_filename[option_size];
// 		    std::sprintf(in_connected_components_filename, "%s/UID_and_connected_component", data_directory.c_str());                   
// 		    
// 		    std::cerr << "\n\nreading CCs....\n\n";		   
// 		    read_connected_component_data_from_file_and_create_Event_shells(in_connected_components_filename, all_CCs_for_repeat_seq);
// 		}
	    
	for (type_list_Sampled_diploid_Event_data::const_iterator it_rd = pseudo_RD_hypothesis_test_work.begin();
		it_rd != pseudo_RD_hypothesis_test_work.end();
		++it_rd)
	{	    
	    read_Connected_Component_ALL_Event_data_from_file(Conn_Comps.at(it_rd->cc_of_event));
	    
	    Conn_Comps.at(it_rd->cc_of_event)
		    .create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights(my_BAM_reader, *it_rd, output_dir);	        
	}//it_rd
	
	
	std::cerr << "\n\n\nsuccessfully completed   all  pseudo_RD_hypothesis_test_work!!!!!\n\n\n";
	return 0;
    }// !pseudo_RD_hypothesis_test_work.empty()
    
    
    
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    

    
    
    
    
    
    
    
    if (!informativeness_checking_parameters.empty())
    {
	std::cerr<< "\n\n\ndetermining informativeness_checking_parameters !!!\n\n\t\t uploading Event data...\n\n";
	char tempoutdir[option_size];
	std::sprintf(tempoutdir, "%s", output_dir.c_str()   );        
	boost::filesystem::create_directory(tempoutdir);
	
	
	//  upload Event Data:
        for (type_map_uint_to_CC::iterator cc_it = Conn_Comps.begin();
		cc_it != Conn_Comps.end();
		++cc_it)
        {                        
            for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
                    it_ev != cc_it->second.events.end();
                    ++it_ev)              
	    {
                read_Event_ALL_data_from_file(cc_it->second, it_ev->second);  	  
		
		it_ev->second.local_interlocking_Events.clear();
		it_ev->second.global_interlocking_Events.clear();
		it_ev->second.my_original_Breakpoint_complex.clear_all_data();
		it_ev->second.remaining_coordinates_for_consideration.clear();
		it_ev->second.natural_poisson_intervals___gender_INsensitive.clear();
		it_ev->second.base_diploid_fragmentation_rate_for_natural_poisson_intervals___gender_sensitive.clear();
		it_ev->second.base_CUMULATIVE_HAPLOID_fragmentation_rate_for_natural_poisson_intervals.clear();
		it_ev->second.observed_read_depth_for_natural_poisson_intervals.clear();
		it_ev->second.uploaded_pieces_of_chromos.clear();
		it_ev->second.PERs_on_this_profile.clear();
		it_ev->second.filter_for_NON_large_gap_of_neighbors_and_self__that_affect_Read_depth.clear();
		it_ev->second.RD_affecting_neighbors_self.clear();			
		it_ev->second.regions_homologous_to_directly_affected_region.clear();
		it_ev->second.regions_homologous_to_LCR_profile.clear();
	    }//ev
	}//cc
	
	std::cerr<< "\n\n\t\tcheck_informativeness_of_every_potential_NAHR_breakpoint informativeness\n";
	remove_degenerate_Events(Conn_Comps);
	
	check_informativeness_of_every_potential_NAHR_breakpoint
		(informativeness_checking_parameters,
		Conn_Comps,
		output_dir,
		BAM_filename.substr(0, BAM_filename.find_last_of('/'))  );
	     
	     
	std::cerr << "\n\n\nDONE checking informativeness!!!!  ending program.\n\n";
	return 0;
    }


	
	





    if (do_post_processing_stats_on_calls)
    {
	for (type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
		it_cc != Conn_Comps.end();
		++it_cc)
	{  read_Connected_Component_ALL_Event_data_from_file(it_cc->second);  }     
			
	
	
	post_process
	    (Conn_Comps,
	     postdatadir,
	     output_dir,     
	     Genes_in_human_genome__filename,
	     biological_validations_filename,
		gender_and_population_filename);	
	
	std::cerr << "\n\nDONE with \"post_process\".  successfully terminating program.\n";	
        return 0;                
    } // do_post_processing_stats_on_calls






    






    if (copy_data_from_most_recent_job_of_same_name)
    {
	if ( !copy_data_from_job_with_jobID.empty()  )
	{	   
	    uint last_cc_analyzed;
	    
	    if (my_MPI_rank == 0)
	    {
		last_cc_analyzed = copy_contents_from_previous_output_directory_with_given_jobID_to_current_output_directory__and__return_last_CC_analyzed(copy_data_from_job_with_jobID);
	    }
	    
	    boost::mpi::broadcast<uint>(*world_MPI_boost_ptr, last_cc_analyzed, 0);
		
	    type_map_uint_to_CC::iterator it_cc = Conn_Comps.begin();
	    while (it_cc->first <= last_cc_analyzed)
	    {  Conn_Comps.erase(it_cc++);  }
	
	    if (my_MPI_rank == 0)
	    {  std::cerr << "\n\nafter removing previously analyzed CC's, number of remaining Conn Comps = " << Conn_Comps.size() << "\n\n";  }	
	}
	else
	{
	
	}	
    
    }//copy_data_from_most_recent_job_of_same_name















    if (only_compute_and_save_natural_poisson_intervals)
    {
        
        std::cerr << "\n\n\n\n\n\n simulating variable elimination schedule and saving Natural Poisson Intervals!!!!\n\n\n\n\n";
    
        type_map_uint_to_CC::iterator cc_it = Conn_Comps.begin();
        while ( cc_it != Conn_Comps.end() )
        {
            
            //  upload Event Data:
            for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
                    it_ev != cc_it->second.events.end();
                    ++it_ev)          
	    {  read_Event_ALL_data_from_file(cc_it->second, it_ev->second);  }
            
                        
            if (my_MPI_rank == 0)
	    {  std::cerr << "\nSuccessfully read in all event data.\n\nreading_visitation_schedule_from_file_for_given_Conn_Comp...\n";  }
           
            //upload visitation schedules            
            const bool success_read_visitation_schedule = read_visitation_schedule_from_file_for_given_Conn_Comp(cc_it->second);
            if (success_read_visitation_schedule)
            {                              
                    //print: 
                    if (my_MPI_rank == 0)
                    {
                        std::fprintf(stderr, "read_visitation_schedule_from_file_for_given_Conn_Comp completed successfully.\n\n");                                  
                        
                        std::fprintf(stderr, "CC_ID = %u\n", cc_it->first);
                        for (type_list_Visitation_object::const_iterator v_it =  cc_it->second.variable_elimination_schedule.begin();
                                v_it != cc_it->second.variable_elimination_schedule.end();
                                ++v_it) 
			{
                            std::fprintf(stderr, "v_it:  event UID = %u,  diploid_compute_size = %u\n",
                                                    v_it->event_for_elimination->UID,
                                                    v_it->diploid_compute_size);    
			}
                    }                                                              
    
            
                const double time_elimination_schedule_begin = omp_get_wtime();
                cc_it->second.simulate_elimination_schedule_and_save_natural_poisson_intervals();  
                                if (my_MPI_rank == 0) 
                                    std::cerr<< "\n\nsimulate_elimination_schedule_and_save_natural_poisson_intervals  for CC:     "  << cc_it->first                                            
                                            << "    total_time    =     "  << omp_get_wtime() - time_elimination_schedule_begin << "\n\n";
            }             
            
            Conn_Comps.erase(cc_it++);
        } // cc_it
        
        
        std::cerr << "\n\n\n\n\n\n\n\n\n\nDone   simulating all variable elimination schedules, and saved all natural poisson intervals!!!\n\n\n\n\n\n";
        
	world_MPI_boost_ptr->barrier();
	remove_directories_associated_with_Partial_Sum_functions_from_scratch_space();
	remove_MPI_rank_directories_from_scratch_space();
        return 0;
                    
    } // only_compute_and_save_natural_poisson_intervals


















































//     if (  !display_MAP_dir.empty()  )
//     {
//         if (my_MPI_rank == 0)
//         {                      
//             
//             upload_universal_variational_positions();
//             upload_universal_indel_positions();
//                         
//             std::fprintf(stderr, "\n\n\nrank 0 displaying pseudo validation for MAP estimates from file [%s].\n\n\n",
//                                     display_MAP_dir.c_str()  );
//             
//             FILE *MAP_file = std::fopen(display_MAP_dir.c_str(), "r");
//                                     if (MAP_file == NULL)
//                                     {
//                                         std::fprintf(stderr, "\n\n\n\n\nERROR: unable to open file [%s] in \"displaying pseudo validation for MAP estimates\".\n\n",
//                                                     display_MAP_dir.c_str() );
//                                         exit(1);
//                                     }            
//                                     
//             
//             uint in_CC_ID;
//             while (  std::fscanf(MAP_file, "\n\nConn_Comp:   %u\n", &in_CC_ID)   != EOF    )
//             {                  
//                 std::fprintf(stderr, "\n\t\tdisplaying for CC = %u\n", in_CC_ID);
//                 Conn_Comp *const the_CC = &Conn_Comps.at(in_CC_ID);
//                 
//                 for (type_map_uint_to_Event::iterator it_ev = the_CC->events.begin();
//                         it_ev != the_CC->events.end();
//                         ++it_ev)          
//                     read_Event_ALL_data_from_file( *the_CC, it_ev->second);                      
//                 
//                 
//                 
//                 std::fprintf(stderr, "\nloading information for Conn Comp   %u\n", in_CC_ID);                
//                 
//                 type_map_uint_to__uint__uint__2 in_MAP_sample;
//                 
//                 for (uint ev_ctr = 0; ev_ctr < the_CC->events.size();  ++ev_ctr)
//                     load_sampled_Event_from_file(MAP_file, in_MAP_sample); 
//                 
//                 std::fscanf(MAP_file, "\t\tP_MAP  =  %*[^\n]\n\n");
//                 
//                 
//                 for (type_map_uint_to__uint__uint::const_iterator it_ev_state = in_MAP_sample.first.begin(),
//                                                                 it_ev_brkpt = in_MAP_sample.second.begin();
//                         it_ev_state != in_MAP_sample.first.end();
//                         ++it_ev_state, ++it_ev_brkpt)
//                 {
//                     if (it_ev_state->second.first != 0)
//                         the_CC->output_reads_covering_certain_breakpoint_regions(it_ev_state->first, 
//                                                                                 it_ev_brkpt->second.first);
//                     if (it_ev_state->second.second != 0)
//                         the_CC->output_reads_covering_certain_breakpoint_regions(it_ev_state->first, 
//                                                                                 it_ev_brkpt->second.second);                                                                                 
//                 }
//                 
//                 std::fprintf(stderr, "\n\t\tDONE displaying for CC = %u\n", in_CC_ID);
//             } 
//                     
//             std::fclose(MAP_file);  
//             
//             std::fprintf(stderr, "\n\n\nrank 0 DONE with pseudo-validation of MAP estimates.\n\n\n");
//         }  // my_MPI_rank == 0        
//             
//         
//         world_MPI_boost_ptr->barrier();
// 	remove_directories_associated_with_Partial_Sum_functions_from_scratch_space();
// 	remove_MPI_rank_directories_from_scratch_space();
//         return 0;               
// 	
// 	
//     } // !display_MAP_dir.empty()



















               
               
                
            
               
               
               
             
  

  
  
  
  
    
    
    world_MPI_boost_ptr->barrier();
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
    if (!display_Event_and_breakpoint_range.empty())
    {                    
	if (my_MPI_rank == 0)                                                                             
	{
	    std::cerr << "\n\n\n\ndisplay_Event_and_breakpoint_range  instead of performing var elim sched...\n\n\n";
	    
	    upload_universal_variational_positions();
	    upload_universal_indel_positions();            
		    
	    for (type_list_Sampled_diploid_Event_data::const_iterator it_disp = display_Event_and_breakpoint_range.begin();
		    it_disp != display_Event_and_breakpoint_range.end();
		    ++it_disp)
	    {
		it_disp->print_this_Sampled_diploid_Event_data();
		
    // 	    read_Event_ALL_data_from_file(Conn_Comps.at(it_disp->cc_of_event), Conn_Comps.at(it_disp->cc_of_event).events.at(it_disp->event_UID));
    // 	    Conn_Comps.at(it_disp->cc_of_event).events.at(it_disp->event_UID).print_this_entry();
		
		Conn_Comps.at(it_disp->cc_of_event).output_reads_covering_certain_breakpoint_regions(*it_disp, true, my_BAM_reader, display_ALL_reads_not_just_good_ones);
		
		print_line_of_markers("++");
	    }                
	    //if (my_MPI_rank == 0)
	    std::fprintf(stderr, "\n\n\n\nDONE display_Event_and_breakpoint_range.\n\n\n");     
	}//0
	
	world_MPI_boost_ptr->barrier();
// 	remove_directories_associated_with_Partial_Sum_functions_from_scratch_space();
// 	remove_MPI_rank_directories_from_scratch_space();
	return 0;
    }
  
      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if (!use_previously_computed_likelihoods___old_data_dir.empty())
    {
	// Now  eliminate variables according to schedule...  
        if (my_MPI_rank == 0)
            std::fprintf(stderr, "\n\n\n\n\nuse saved likelihoods to do variable eliminaiton for diffeent priors!!!...\n\n\n\n\n");
	
	if (my_MPI_rank == 0)
	{//erase unneeded directories
	    char sub_out_dirs[option_size];
	    std::sprintf(sub_out_dirs, "%s/Partial_sum_functions", output_dir.c_str());
	    boost::filesystem::remove_all(sub_out_dirs);
	    
	    std::sprintf(sub_out_dirs, "%s/Star_alignments", output_dir.c_str());
	    boost::filesystem::remove_all(sub_out_dirs);            
	    
	    std::sprintf(sub_out_dirs, "%s/Raw_alignments", output_dir.c_str());
	    boost::filesystem::remove_all(sub_out_dirs); 
	    
	    std::sprintf(sub_out_dirs, "%s/Samples_from_posterior", output_dir.c_str());
	    boost::filesystem::remove_all(sub_out_dirs);     
			
	    std::sprintf(sub_out_dirs, "%s/Likelihoods",  output_dir.c_str()  );     
	    boost::filesystem::remove_all(sub_out_dirs); 
				    
	    //scratch
	    for (int rank=0; rank < size_MPI_WORLD; ++rank)
	    {
		std::sprintf(sub_out_dirs, "%s/rank_%d",
					scratch_job_output_dir.c_str(), rank );        
		boost::filesystem::remove_all(sub_out_dirs);  
	    }		    
	}//erase unneeded directories
	
	
	const type_set_uint excluded_Events_from_old_run(
				  read_totally_excluded_Events_from_previous_job_save_dir(use_previously_computed_likelihoods___old_data_dir)  );  	
	
	
	
		
	const std::string SAVED_output_dir( output_dir );
	const std::string SAVED_scratch_job_output_dir( scratch_job_output_dir );
	
	
        world_MPI_boost_ptr->barrier();	                	
	for (type_list_int::const_iterator it_prior_exp = list_of_prior_probability_exponents_for_computing.begin();
		it_prior_exp != list_of_prior_probability_exponents_for_computing.end();
		++it_prior_exp)
	{
	    P_E_absent_cond_theta = 1.00L;
	    P_E_absent_cond_theta -= mpfr::exp10( real(*it_prior_exp) );
	    
	    P_E_inv_cond_theta = mpfr::exp10( real(*it_prior_exp) );
	    
	    P_E_dup_or_del_cond_theta =  (P_E_inv_cond_theta / 2.00L);
	    	    
//             P_E_absent_cond_theta =   1.00L - (longdouble)exp10( *it_prior_exp );
//             P_E_inv_cond_theta = (1.00L - P_E_absent_cond_theta);
//             P_E_dup_or_del_cond_theta = (P_E_inv_cond_theta / 2.00L);
            
		    if (my_MPI_rank == 0)
		    {
			std::cerr<< "\n\nit_prior_exp = " << *it_prior_exp << "\n";
			std::cerr<< "set    P_E_absent_cond_theta  =  "  <<  P_E_absent_cond_theta <<  "\n";
			std::cerr<< "\t\t\t\tset    P_E_inv_cond_theta  =  "  <<  P_E_inv_cond_theta <<  "\n";
			std::cerr<< "\t\t\t\tset    P_E_dup_or_del_cond_theta  =  "  <<  P_E_dup_or_del_cond_theta <<  "\n\n";
		    }
	    //dirs for this prior
	    if (my_MPI_rank == 0)
	    {
		char suffix_prior[option_size];
		std::sprintf(suffix_prior, "%s/prior_%d", SAVED_output_dir.c_str(), (-1)*(*it_prior_exp));		
		output_dir.assign(suffix_prior);
		boost::filesystem::create_directory(suffix_prior);
		
		std::sprintf(suffix_prior, "%s/prior_%d", SAVED_scratch_job_output_dir.c_str(), (-1)*(*it_prior_exp)  );
		scratch_job_output_dir.assign(suffix_prior);
		boost::filesystem::create_directory(suffix_prior);
	    }
	    world_MPI_boost_ptr->barrier(); // for creating directories
	    	    	    
	    
	    create_scratch_directory_tree();
	    create_output_directory_tree();
	    
	    world_MPI_boost_ptr->barrier(); // for creating directories
	    
	    
	    
	    Conn_Comps.clear();
	    {//scope    
			char in_connected_components_filename[option_size];
			std::sprintf(in_connected_components_filename, "%s/UID_and_connected_component", data_directory.c_str());   
			
			read_connected_component_data_from_file_and_create_Event_shells
								    (in_connected_components_filename,
								    Conn_Comps);                                                                                
	    }//scope	       
	    
	    {		
		type_set_uint desired_Events;
		for (type_list_Sampled_diploid_Event_data::const_iterator it_disp = display_Event_and_breakpoint_range.begin();
			it_disp != display_Event_and_breakpoint_range.end();
			++it_disp)
		{  desired_Events.insert(it_disp->event_UID);  }		
	    
		remove_undesirable_Connected_Components_from_set_of_Connected_Components
				    (Conn_Comps,
				    targeted_CC_IDs,
				    desired_Events,
				    acceptable_Connected_Component_size_range,
				    acceptable_Connected_Component_compute_range,
				    last_completed_CC);
	    }
				
		             
	    type_map_uint_to_CC::iterator cc_it = Conn_Comps.begin();
	    while ( cc_it != Conn_Comps.end() )
	    {	    	    
		
		if (my_MPI_rank == 0)
		{
		    std::cerr << "\n\n\n\nnumber of remaining Conn Comps:   " << Conn_Comps.size() << "\n\n\n";
		    std::cerr << "\n\treading in event data for connected component  " << cc_it->first  << ",  size = " << cc_it->second.events.size() << "\n";
		}		
		//  upload Event Data:
		for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
			it_ev != cc_it->second.events.end();
			++it_ev)              
		    read_Event_ALL_data_from_file(cc_it->second, it_ev->second);            
		
		
			    
		if (my_MPI_rank == 0)
		    std::fprintf(stderr, "Successfully read in all event data.\n\nreading_visitation_schedule_from_file_for_given_Conn_Comp...");                 

	    
		//upload visitation schedules            
		const bool successfully_read_visitation_schedule = read_visitation_schedule_from_file_for_given_Conn_Comp(cc_it->second);
	    
		if (!successfully_read_visitation_schedule)
		{        
		    if (my_MPI_rank == 0)
			append_to_skipped_connected_component_file( cc_it->second,  "failed to read in visitation schedule." );
		    Conn_Comps.erase(cc_it++);
		    continue;            
		}
				//print: 
				if (my_MPI_rank == 0)
				{
				    std::fprintf(stderr, "read_visitation_schedule_from_file_for_given_Conn_Comp completed successfully.\n\n");                                  
				    
				    std::fprintf(stderr, "CC_ID = %u\n", cc_it->first);
				    for (type_list_Visitation_object::const_iterator v_it =  cc_it->second.variable_elimination_schedule.begin();
					    v_it != cc_it->second.variable_elimination_schedule.end();
					    ++v_it)    
					std::fprintf(stderr, "v_it:  event UID = %u,  diploid_compute_size = %u\n",
								v_it->event_for_elimination->UID,
								v_it->diploid_compute_size);    
				}            
		
	
	    

		
				const double time_elimination_schedule_begin = omp_get_wtime();
				
		const int flag_for_variable_elimination
// 			    = cc_it->second.reorganize_saved_Likelihoods_for_all_connected_components
			    = cc_it->second.use_saved_likelihoods_to_perform_complete_variable_elimination_schedule_with_current_prior
					    (use_previously_computed_likelihoods___old_data_dir,
					     excluded_Events_from_old_run  );
					    
				if (my_MPI_rank == 0) 
				{
				    std::cerr<< "\n\nperform_complete_variable_elimination_schedule for CC:     "  << cc_it->first
						<<  "\n\t flag_for_variable_elimination =  "  << flag_for_variable_elimination
					    << "    total_time    =     "  << omp_get_wtime() - time_elimination_schedule_begin << "\n\n";
				}
		    

				    
		mpfr_free_cache();                             
		std::cerr.flush();std::fflush(stderr);std::cerr.flush();
		world_MPI_boost_ptr->barrier();
		

		if (flag_for_variable_elimination == 0)
		{          	                                    
		    //re-load, so that we can get display...
		    for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
			    it_ev != cc_it->second.events.end();
			    ++it_ev)     
		    {
			it_ev->second.clear_all_data();
			read_Event_ALL_data_from_file(cc_it->second, it_ev->second);                    
		    }
		    		    		    
	
		    if (my_MPI_rank == 0)
		    {  std::cerr << "\n\n\n....get_and_save_MAP_estimate....\n\n\n";  }
		    
		    const type_map_uint_to_Sampled_diploid_Event_data pseudo_Validation_Events(cc_it->second.get_and_save_MAP_estimate());
		    
		    mpfr_free_cache();
		    world_MPI_boost_ptr->barrier();
		    

		    if (my_MPI_rank == 0   and   number_of_samples_from_posterior > 0) 
		    {
			cc_it->second.sample_from_posterior__and__save____marginal_distributions__and__outcome_Centroid__and_Shannon__Entropies
												    (number_of_samples_from_posterior);
		    }
																	    
		    
		    mpfr_free_cache();                                                                                                                    
		    std::cerr.flush();std::fflush(stderr);std::cerr.flush();
		    world_MPI_boost_ptr->barrier();
		    
		}  // successful completion of perform_variable_elimination
		else if (my_MPI_rank == 0)   // failure.
		{		    			    		 
		    std::stringstream error_strm;
		    print_line_of_markers("ERROR! ", &error_strm);
		    print_line_of_markers("(", &error_strm);
		    error_strm << "\n\nERROR!   perform_complete_variable_elimination_schedule  failed  with  fail code:    " 
				<<  flag_for_variable_elimination  << "\n\n";
		    
		    print_line_of_markers(")", &error_strm);
		    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );                		                                  
		}		
	    
	    
	    
		//erase:
	    
		for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
			it_ev != cc_it->second.events.end();
			++it_ev)                       
		    it_ev->second.clear_all_data();  
		    
		cc_it->second.partial_sum_functions.clear();
		cc_it->second.variable_elimination_schedule.clear();		    
			
		    	    		             
		world_MPI_boost_ptr->barrier();					
		Conn_Comps.erase(cc_it++);  //post-increment necessary!            
	    }    // all cc_it  
	    
	    
	    remove_directories_associated_with_Partial_Sum_functions_from_scratch_space();
	    remove_MPI_rank_directories_from_scratch_space();
	    
// 	    std::cerr << "\n\n\nDONE WITH 1 PRIOR!!!!!\n\n\n";
// 	    return 0;
	    	    	
	} // prior	
	
	
	
	output_dir.assign( SAVED_output_dir );
	scratch_job_output_dir.assign( SAVED_scratch_job_output_dir );
	
	world_MPI_boost_ptr->barrier();
	remove_directories_associated_with_Partial_Sum_functions_from_scratch_space();
	remove_MPI_rank_directories_from_scratch_space();
	
	if (my_MPI_rank == 0)
	    std::cerr << "\n\n\n\n DONE WITH  \"use_saved_likelihoods_to_perform_complete_variable_elimination_schedule_with_current_prior\" !!!!\n\n\n\n";
	
	return 0;
        
    } // !use_previously_computed_likelihoods___old_data_dir.empty()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


  
  
  
  
      
    
    
    
    //if We made it this far, then we perform the entire variable elimination schedule!!!!!
    {                                      
	// Now  eliminate variables according to schedule...  
        if (my_MPI_rank == 0)
	{  std::fprintf(stderr, "\n\n\n\n\nperform_complete_variable_elimination_schedules!!!...\n\n\n\n\n");  }
                
        
        
        
        world_MPI_boost_ptr->barrier();        
        
        
        while (!Conn_Comps.empty())
        {	   
	    const type_map_uint_to_CC::iterator cc_it = --Conn_Comps.end();
	    
            if (my_MPI_rank == 0)
            {
                std::cerr << "\n\n\n\nnumber of remaining Conn Comps:   " << Conn_Comps.size() << "\n\n\n";
                std::cerr << "\n\treading in event data for connected component  " << cc_it->first  << ",  size = " << cc_it->second.events.size() << "\n";
            }
            
            //  upload Event Data:
            for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
                    it_ev != cc_it->second.events.end();
                    ++it_ev)              
	    {  read_Event_ALL_data_from_file(cc_it->second, it_ev->second);  }          
            
            
                        
            if (my_MPI_rank == 0)
	    {  std::fprintf(stderr, "Successfully read in all event data.\n\nreading_visitation_schedule_from_file_for_given_Conn_Comp...");  }

           
            //upload visitation schedules            
            const bool successfully_read_visitation_schedule = read_visitation_schedule_from_file_for_given_Conn_Comp(cc_it->second);
        
            if (!successfully_read_visitation_schedule)
            {        
                if (my_MPI_rank == 0)
		{  append_to_skipped_connected_component_file( cc_it->second,  "failed to read in visitation schedule." );  }
                Conn_Comps.erase(cc_it);
                continue;            
            }
                            //print: 
                            if (my_MPI_rank == 0)
                            {
                                std::fprintf(stderr, "read_visitation_schedule_from_file_for_given_Conn_Comp completed successfully.\n\n");                                  
                                
                                std::fprintf(stderr, "CC_ID = %u\n", cc_it->first);
                                for (type_list_Visitation_object::const_iterator v_it =  cc_it->second.variable_elimination_schedule.begin();
                                        v_it != cc_it->second.variable_elimination_schedule.end();
                                        ++v_it)    
				{
                                    std::fprintf(stderr, "v_it:  event UID = %u,  diploid_compute_size = %u\n",
                                                            v_it->event_for_elimination->UID,
                                                            v_it->diploid_compute_size);    
				}
                            }
            
    
          
            if (cc_it->second.total_diploid_compute_size > 500000)
            {
                if (my_MPI_rank == 0)
                {
                    std::stringstream reason_strm;
                    reason_strm << "total connected component diploid compute size = " << cc_it->second.total_diploid_compute_size << ", which is deemed too large.";
                    append_to_skipped_connected_component_file( cc_it->second, reason_strm.str() ); 
                }
                Conn_Comps.erase(cc_it);
                continue;                      
            }
            
            
            
	    const double time_elimination_schedule_begin = omp_get_wtime();
			    			    
	    const int flag_for_variable_elimination 
			    = cc_it->second.perform_complete_variable_elimination_schedule(
				    my_BAM_reader,
				    consider_GeneConversion_breakpoints_in_Heuristic,
				    consider_GeneConversion_in_Full_compute);
	    
	    
                            if (my_MPI_rank == 0) 
                                std::cerr<< "\n\nperform_complete_variable_elimination_schedule for CC:     "  << cc_it->first
                                            <<  "\n\t flag_for_variable_elimination =  "  << flag_for_variable_elimination
                                         << "    total_time    =     "  << omp_get_wtime() - time_elimination_schedule_begin << "\n\n";
                

                                
            mpfr_free_cache();                             
            world_MPI_boost_ptr->barrier();
            
            if (flag_for_variable_elimination != 0)  // failure.
            {
                
                for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
                        it_ev != cc_it->second.events.end();
                        ++it_ev)                       
		{  it_ev->second.clear_all_data();  }
                
                if (my_MPI_rank == 0)    
                {
                    std::stringstream error_strm;
                    error_strm << "\n\nERROR!   perform_complete_variable_elimination_schedule  failed  with  fail code:    " 
                                <<  flag_for_variable_elimination  << "\n\n";
                    
                    error_message(error_strm,false);        
                }                                
            }
            else if (use_previously_computed_likelihoods___old_data_dir.empty()) //success
            {          		
		mpfr_free_cache();
		
		upload_universal_variational_positions();
		upload_universal_indel_positions();		
		
		//re-load, so that we can get display...
		read_Connected_Component_ALL_Event_data_from_file(cc_it->second);
		
		if (my_MPI_rank == 0)
		{  std::cerr << "\n\n\n....get_and_save_MAP_estimate....\n\n\n";  }
		
		const type_map_uint_to_Sampled_diploid_Event_data pseudo_Validation_Events(cc_it->second.get_and_save_MAP_estimate());
		
		mpfr_free_cache();
		world_MPI_boost_ptr->barrier();		
		
		for (type_map_uint_to_Sampled_diploid_Event_data::const_iterator it_result = pseudo_Validation_Events.begin();
			it_result != pseudo_Validation_Events.end();
			++it_result)
		{
		    //alignments:		    
		    cc_it->second.output_reads_covering_certain_breakpoint_regions(
					    it_result->second, true, my_BAM_reader);
		    
		    //RD:
		    read_Event_ALL_data_from_file(
					    cc_it->second, cc_it->second.events.at(it_result->second.event_UID));
		    
		    cc_it->second.create_inferred_Read_Depth_graph_using_per_read_posterior_distributions_for_weights(
					    my_BAM_reader, it_result->second, output_dir);
		}
		
		
		
													
		
		universal_variational_positions.clear(); 
		universal_del_positions.clear();
		universal_ins_positions.clear();       
		
		
		if (my_MPI_rank == 0   and   number_of_samples_from_posterior > 0) 
		{
		    read_Connected_Component_ALL_Event_data_from_file(cc_it->second);
		    
		    cc_it->second.sample_from_posterior__and__save____marginal_distributions__and__outcome_Centroid__and_Shannon__Entropies(
												number_of_samples_from_posterior);
		}
				
												
		for (type_map_uint_to_Event::iterator it_ev = cc_it->second.events.begin();
			it_ev != cc_it->second.events.end();
			++it_ev)                       
		{  it_ev->second.clear_all_data();  }                
			    
		mpfr_free_cache();                                          
		
                world_MPI_boost_ptr->barrier();                        
                
            }  // successful completion of perform_variable_elimination
            
            
	    cc_it->second.partial_sum_functions.clear();
	    cc_it->second.variable_elimination_schedule.clear();            
            
            
            
            
            
            if (my_MPI_rank == 0)             
            {       
                std::string last_cc_filename( scratch_job_output_dir );
                last_cc_filename.append( "/" );
                last_cc_filename.append( "last_cc_analyzed" );
                
                std::ofstream last_cc_file( last_cc_filename );
                if (!last_cc_file.good())        
                    std::fprintf(stderr, "\n\n\n\n\nERROR: Command [%s] could not be opened for 'save_last_cc_file'.\n\n", 
                                        last_cc_filename.c_str() );        
                else
                {
                    last_cc_file << cc_it->first;
                    last_cc_file.close();
                }            
            }              
            world_MPI_boost_ptr->barrier();
            
                      
            Conn_Comps.erase(cc_it);  //post-increment necessary!            
        }//cc
    }           
  
  
  


  
  
  
  
  
  
  
  
  
  //assuming our connected components and their eventsd have been created and initialized with loaded data...
  //Then within each connected component, we want first (say), for each event to collection its reads and calculate the probability of each breakpoint (i.e. calculate the likelihood of hybrid read data given breakpoints).
//   for (unsigned int cc=0; cc<num_global_conn_comps; cc++)
//     Conn_Comps[cc].calculate_P_Dhybrid_for_each_event();
  
  
  

    world_MPI_boost_ptr->barrier();
    remove_directories_associated_with_Partial_Sum_functions_from_scratch_space();
    remove_MPI_rank_directories_from_scratch_space();
    
                                                                                
    
    
    std::cerr << "\n\n\n\n\n\n\n\n\n\t\t\t\t\t\t\t\t\tEnd of main.\n\t\t\t\t\t\t\tThe sun is shining\n\t\t\t\t\t\t\tBirds are singing.\n\t\t\t\t\t\t\t   Happiness\n\n\n";
    return 0;
} // end of main
