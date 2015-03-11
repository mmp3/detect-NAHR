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



#include <api/BamReader.h>



#include <general_typedefs.h>
#include <Readgroup_statistics.h>

// 
#include "globals.h"


bool is_simulation = true;



type_map_uint_to_vector_uint  map_chromosome_to_cumulative_GC_counts;


type_map_uint_to_uint map_chromosome_value_to_BAM_Ref_IDs;

boost::mpi::communicator  *world_MPI_boost_ptr;


type_map_uint_to_set_uint universal_variational_positions;
type_map_uint_to_set_uint universal_del_positions;
type_map_uint_to_set_uint universal_ins_positions;



uint genome_tag = 0; //for debugging



longdouble heuristic_posterior_PROBABILITY_threshold = 0.95L;
uint heuristic_posterior_COUNT_threshold = 3;


        
int maximum_frag_length_absolute_dev = 1;        

// int PER_proper_pair_deviation_amount = 40;


uint padding_for_region_homologous_to_a_var_pos = 100;

uint max_number_of_OMP_threads;




uint search_block_size = 5000;






uint ASCII_quality_offset = 33;


std::string  data_directory;

std::string  Ref_genome_fileprefix;
std::string  scratch_job_output_dir;
std::string output_dir;



std::string  display_MAP_dir;


bool keep_sparse_variational_positions_for_dupdels = false;


std::string  visitation_schedule_subdirname(  "visitation_schedules"  );

longdouble RD_variance_scale___alpha_squared = (0.03L*0.03L);
longdouble lambda_0 = 0.00L;
longdouble haploid_coverage = 0.00L;

// type_map_string_to_longdouble coverage_by_RG;

longdouble base_noise_read_depth_percentage_of_normal_Reference_read_depth = 0.0005;



uint length_of_haploid_genome = 0;
uint maximum_average_fragment_length = 0;


uint minimum_natural_poisson_interval_size = 40;




// #define  PRIOR_FRACTION_FOR_ONE_TYPE_OF_GENE_CONVERSION__VS__NAHR    (0.9999L)

real P_E_absent_cond_theta;//(0.999999L);
real P_E_GeneConv_onetype_cond_theta;//(((real(1.00L) - P_E_absent_cond_theta)*PRIOR_FRACTION_FOR_ONE_TYPE_OF_GENE_CONVERSION__VS__NAHR)/2);
real P_E_inv_cond_theta;//((real(1.00L) - (P_E_absent_cond_theta + P_E_GeneConv_onetype_cond_theta*2)));
real P_E_dup_or_del_cond_theta;//((P_E_inv_cond_theta*0.99L)/2);






type_map_string_to_Readgroup_stats Readgroup_stats;






uint maximum_number_of_breakpoints_to_consider_after_heuristic = 4;


std::string variational_positions_directory;

std::string  genome_name;
std::string  population_abbreviation;


bool gender_of_individual_is_female = false;







real minimum_representable_positive_real_number;
























// void print_Balgmt( const BamTools::BamAlignment &some_algmt, std::stringstream *const &some_ss)
// {
//     
//     std::stringstream output_str;
//     
//     const bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
//     std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str; 
// 
//     print_line_of_markers("[[", ss_ptr);
// 
//     (*ss_ptr) << "\nprinting    BamTools::BamAlignment:\n";    
//     
//     (*ss_ptr) << "\n";
//     (*ss_ptr) << "Name: " << some_algmt.Name.c_str()  << "\n";
//     (*ss_ptr) << "\t Length = " <<  some_algmt.Length       << "\n";
//     (*ss_ptr) << "\t QueryBases =   " <<  some_algmt.QueryBases.c_str()  << "\n";
//     (*ss_ptr) << "\t AlignedBases = " <<  some_algmt.AlignedBases.c_str()      << "\n";
//     (*ss_ptr) << "\t Qualities = " <<     some_algmt.Qualities.c_str()   << "\n";
//     (*ss_ptr) << "\t TagData = " <<     some_algmt.TagData.c_str()   << "\n";
//     (*ss_ptr) << "\t RefID = " <<      some_algmt.RefID   << "\n";
//     (*ss_ptr) << "\t Position = " <<      some_algmt.Position   << "\n";
//     (*ss_ptr) << "\t Bin = " <<       some_algmt.Bin   << "\n";
//     (*ss_ptr) << "\t MapQuality = " <<      some_algmt.MapQuality   << "\n";
//     (*ss_ptr) << "\t AlignmentFlag = " <<      some_algmt.AlignmentFlag   << "\n";
//     (*ss_ptr) << "\t MateRefID = " <<      some_algmt.MateRefID   << "\n";
//     (*ss_ptr) << "\t MatePosition = " <<      some_algmt.MatePosition   << "\n";
//     (*ss_ptr) << "\t InsertSize = " <<     some_algmt.InsertSize   << "\n";
//     (*ss_ptr) << "\t Filename = " <<     some_algmt.Filename.c_str()   << "\n";
// 
//     (*ss_ptr) << "\t 0x0001: .IsPaired = " <<     some_algmt.IsPaired()   << "\n";
//     (*ss_ptr) << "\t 0x0002: .IsProperPair = " <<     some_algmt.IsProperPair()   << "\n";
//     (*ss_ptr) << "\t 0x0004: .IsMapped = " <<     some_algmt.IsMapped()   << "\n";
//     (*ss_ptr) << "\t 0x0008: .IsMateMapped = " <<     some_algmt.IsMateMapped()   << "\n";
//     (*ss_ptr) << "\t 0x0010: .IsReverseStrand = " <<     some_algmt.IsReverseStrand()   << "\n";
//     (*ss_ptr) << "\t 0x0020: .IsMateReverseStrand = " <<     some_algmt.IsMateReverseStrand()   << "\n";
//     (*ss_ptr) << "\t 0x0040: .IsFirstMate = " <<     some_algmt.IsFirstMate()   << "\n";
//     (*ss_ptr) << "\t 0x0080: .IsSecondMate = " <<     some_algmt.IsSecondMate()   << "\n";
//     (*ss_ptr) << "\t 0x0100: .IsPrimaryAlignment = " <<    some_algmt.IsPrimaryAlignment()   << "\n";
//     (*ss_ptr) << "\t 0x0200: .IsFailedQC = " <<     some_algmt.IsFailedQC()   << "\n";
//     (*ss_ptr) << "\t 0x0400: .IsDuplicate = " <<     some_algmt.IsDuplicate()   << "\n";
// 
//     print_line_of_markers("]]", ss_ptr);
// 
//     if (!provided_with_valid_sstream)
//     {
//         std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );       
//         std::cerr.flush();std::fflush(stderr);
//     }
//   
// 
// } //end print_algmt


void print_Balgmt( const BamTools::BamAlignment &some_algmt,
                   std::stringstream *const &some_ss)
{
    
    std::stringstream output_str;
    
    bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str;    
    
    
    (*ss_ptr) << "\n\n"
    << "Name: " << some_algmt.Name << "\n"
    << "\t Length = " << some_algmt.Length << "\n"
    << "\t QueryBases =   " << some_algmt.QueryBases << "\n"
    << "\t AlignedBases = " << some_algmt.AlignedBases << "\n"
    << "\t Qualities = " << some_algmt.Qualities << "\n"
//     << "\t TagData = " << some_algmt.TagData << "\n"
    << "\t RefID = " << some_algmt.RefID << "\n"
    << "\t Position = " << some_algmt.Position << "\n"
    << "\t Bin = " << some_algmt.Bin << "\n"
    << "\t MapQuality = " << some_algmt.MapQuality << "\n"
    << "\t AlignmentFlag = " << some_algmt.AlignmentFlag << "\n"
    << "\t MateRefID = " << some_algmt.MateRefID << "\n"
    << "\t MatePosition = " << some_algmt.MatePosition << "\n"
    << "\t InsertSize = " << some_algmt.InsertSize << "\n"
//     << "\t Filename = " << some_algmt.Filename << "\n"

    << "\t 0x0001: .IsPaired = " << (int)some_algmt.IsPaired() << "\n"
    << "\t 0x0002: .IsProperPair = " << (int)some_algmt.IsProperPair() << "\n"
    << "\t 0x0004: .IsMapped = " << (int)some_algmt.IsMapped() << "\n"
    << "\t 0x0008: .IsMateMapped = " << (int)some_algmt.IsMateMapped() << "\n"
    << "\t 0x0010: .IsReverseStrand = " << (int)some_algmt.IsReverseStrand() << "\n"
    << "\t 0x0020: .IsMateReverseStrand = " << (int)some_algmt.IsMateReverseStrand() << "\n"
    << "\t 0x0040: .IsFirstMate = " << (int)some_algmt.IsFirstMate() << "\n"
    << "\t 0x0080: .IsSecondMate = " << (int)some_algmt.IsSecondMate() << "\n"
    << "\t 0x0100: .IsPrimaryAlignment = " << (int)some_algmt.IsPrimaryAlignment() << "\n"
    << "\t 0x0200: .IsFailedQC = " << (int)some_algmt.IsFailedQC() << "\n"
    << "\t 0x0400: .IsDuplicate = " << (int)some_algmt.IsDuplicate() << "\n"

    << "\n\n";
    
    if (!provided_with_valid_sstream)
        std::fprintf(stderr, "\n%s\n", (*ss_ptr).str().c_str() );    

} //end print_Balgmt





