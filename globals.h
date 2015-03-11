#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED


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
#include <api/BamAlignment.h>


#include <general_typedefs.h>
#include <Readgroup_statistics.h>
#include <Sparse_map.h>


extern  bool is_simulation;


extern   type_map_uint_to_vector_uint  map_chromosome_to_cumulative_GC_counts;

typedef  mpfr::mpreal  real;

typedef  std::pair<type_uint__uint, uint>  type_uint__uint__uint;
typedef  std::list<type_uint__uint__uint>  type_list_uint__uint__uint;

typedef  std::pair<type_uint__uint, type_longdouble__longdouble>   type_uint2_longdouble2;
typedef  std::list<type_uint2_longdouble2>   type_list__uint2_longdouble2;    




extern real minimum_representable_positive_real_number;

//parallel
extern boost::mpi::communicator  *world_MPI_boost_ptr;
extern uint max_number_of_OMP_threads;



extern uint ASCII_quality_offset;


//genome
extern type_map_uint_to_uint map_chromosome_value_to_BAM_Ref_IDs;

extern type_map_uint_to_set_uint universal_variational_positions;
extern type_map_uint_to_set_uint universal_del_positions;
extern type_map_uint_to_set_uint universal_ins_positions;



extern uint length_of_haploid_genome;

extern int maximum_frag_length_absolute_dev;

extern uint minimum_natural_poisson_interval_size;



extern uint search_block_size;







extern longdouble RD_variance_scale___alpha_squared;
extern longdouble lambda_0;
extern longdouble haploid_coverage;

// extern  type_map_string_to_longdouble coverage_by_RG;

extern longdouble base_noise_read_depth_percentage_of_normal_Reference_read_depth;

extern longdouble heuristic_posterior_PROBABILITY_threshold;
extern uint heuristic_posterior_COUNT_threshold;




extern uint maximum_average_fragment_length;


//probability:
extern real P_E_absent_cond_theta;
extern real P_E_GeneConv_onetype_cond_theta;
extern real P_E_inv_cond_theta;
extern real P_E_dup_or_del_cond_theta;





//debug:
extern uint genome_tag; //for debugging




extern uint padding_for_region_homologous_to_a_var_pos;









// I/O
extern std::string  data_directory; 

extern std::string  Ref_genome_fileprefix;
extern std::string  scratch_job_output_dir;
extern std::string output_dir;

extern std::string  display_MAP_dir;


extern bool keep_sparse_variational_positions_for_dupdels;



extern std::string  visitation_schedule_subdirname;





extern uint maximum_number_of_breakpoints_to_consider_after_heuristic;



extern  std::string  variational_positions_directory;

extern  std::string  genome_name;
extern  std::string  population_abbreviation;

extern  bool gender_of_individual_is_female;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// typedefs:



typedef std::list<real>  type_list_real;

typedef std::pair<uint, real>   type_uint__real;
typedef std::pair< type_uint__uint,  real >  type_uint__uint___real;

typedef std::pair<real, real>    type_real__real;
typedef  std::pair<type_int__int, type_real__real>   type_int__int___real__real;
typedef std::pair<uint, type_real__real>    type_uint___real__real;
// typedef std::pair<type_uint__uint2, real>   type_uint__uint2___real;
typedef std::map<uint, real>   type_map_uint_to_real; //for partial sum functions
typedef  std::map<uint, type_real__real>   type_map_uint_to___real__real;
typedef std::map<type_uint__uint, real>   type_map_uint__uint_to_real;
// typedef std::map<type_uint__uint2, real>  type_map_uint__uint2___to___real;

typedef std::map<type_BI__BI, real, compare_BI> type_map_BI__BI___to__real; 
                                    
                                                            //eliminating event breakpoints to  hybrid factor value

typedef  std::multimap<real, type_uint__uint>    type_multimap_real_to__uint__uint;
typedef  std::multimap<real, type_BI__BI>  type_multimap_real_to__BI__BI;



typedef std::multimap<real, type_uint__uint>   type_multimap_real_to_uint__uint;




typedef   std::map<type_uint__uint, longdouble>   type_map__uint__uint__to__longdouble;
typedef   std::pair<longdouble, type_map__uint__uint__to__longdouble>  type__longdouble___map__uint__uint__to__longdouble;
typedef   std::map<uint, type__longdouble___map__uint__uint__to__longdouble>    type_map_uint__to___longdouble___map__uint__uint__to__longdouble;
                                                                                    //entropies.


                            //marginal probability of outcome  ,  map breakpoint to marginal probability conditioned on outcome
typedef  std::map<type_uint__uint,type__longdouble___map__uint__uint__to__longdouble>    type_map__uint__uint___to___longdouble___map__uint__uint__to__longdouble;
                            // map  outcome  to      marginal probability of outcome  ,  map breakpoint to marginal probability conditioned on outcome

typedef  std::map<uint, type_map__uint__uint___to___longdouble___map__uint__uint__to__longdouble> 
                                    type_map_uint__to__map__uint__uint___to___longdouble___map__uint__uint__to__longdouble;

                                    

typedef   std::vector<real>   type_vector_real;


typedef  std::map<uint, type_map_uint_to_real>  type_map_uint_to_uint_to_real;



typedef   std::map<BOOST_Interval, type_real__real, compare_BI>   type_map_BI_to__real__real;



typedef  std::pair<BOOST_Interval, type_real__real>   type_BI__real__real;


// BAMTOOLS:
typedef std::list<BamTools::BamAlignment>  type_list_Balgmts;
typedef std::pair<BamTools::BamAlignment, BamTools::BamAlignment>  type_Balgmt__Balgmt;




typedef std::pair<BamTools::BamAlignment, BOOST_Interval>   type_Balgmt__BOOST_Interval;
typedef std::vector<type_Balgmt__BOOST_Interval>   type_vector___Balgmt__BOOST_Interval;







typedef  std::pair<type_map_uint_to_real,type_map_uint_to_real>   type_map_uint_to_real__2;



// class Empirical_PER_data
// {
//     public:
//         Empirical_PER_data()    :   observed_number_of_paired_reads(0),
//                                     observed_average_fragment_length(0.00L),
//                                     observed_average_mate_length(0)
// //                                     observed_average_advance_size(0),
// //                                     empirical_advance_size_standard_deviation(0)
//         { }
//         
//         Empirical_PER_data(     const uint &in_observed_number_of_paired_reads,
//                                 const longdouble &in_observed_average_fragment_length,
//                                 const longdouble &in_observed_average_mate_length   )
// //                                 const uint &in_observed_average_advance_size,
// //                                 const uint &in_empirical_advance_size_standard_deviation)
//                                         :   observed_number_of_paired_reads(in_observed_number_of_paired_reads),
//                                             observed_average_fragment_length(in_observed_average_fragment_length),
//                                             observed_average_mate_length(in_observed_average_mate_length)
// //                                             observed_average_advance_size(in_observed_average_advance_size),
// //                                             empirical_advance_size_standard_deviation(in_empirical_advance_size_standard_deviation)
//         { }
//         
//         ~Empirical_PER_data()
//         { }
//         
//         
//         //data:        
//         const uint observed_number_of_paired_reads;
//         
//         const longdouble observed_average_fragment_length;
//     
//         const longdouble observed_average_mate_length;
//         
// //         const uint observed_average_advance_size;        
// //         const uint empirical_advance_size_standard_deviation;
// };






























namespace boost {
    namespace serialization {

        
        
        
        //mpfr::mpreal
        template<class Archive>
        void save(Archive & ar, const real &some_real, const uint version)
        {                            
            const std::string out_str( some_real.toString() );
            ar <<  out_str;
        }       
        template<class Archive>
        void load(Archive & ar, real &some_real, const uint version)
        {              
            std::string in_str_val;
            ar >> in_str_val;            
            some_real = in_str_val.c_str();
        }       
        
        template<class Archive>
        inline void serialize( Archive & ar, real &some_real, const unsigned int file_version)
        {
            split_free(ar, some_real, file_version); 
        }
        
        
        
        
        
        //BOOST_Interval
        template<class Archive>
        void save(Archive & ar, const BOOST_Interval &some_BI, const unsigned int version)
        {
            ar << some_BI.lower();
            ar << some_BI.upper();                   
        }
        
        template<class Archive>
        void load(Archive & ar, BOOST_Interval &some_BI, const unsigned int version)
        {
            uint in_low, in_upp;
            ar >> in_low;
            ar >> in_upp;
            some_BI.set(in_low, in_upp);            
        }
        
        template<class Archive>
        inline void serialize( Archive & ar, BOOST_Interval &some_BI, const unsigned int file_version)
        {
            split_free(ar, some_BI, file_version); 
        }          
        
        
        
        
        
        //Sparse_map        
        template<class Archive>
        void serialize(Archive & ar, Sparse_map & smap, const unsigned int version)
        {
            ar  &  smap.sparse_map_UID_to_value;
        }
        
        
        
        
        
        
        
        
        
        //BammTools::Bamalignment
        template<class Archive>
        void save(Archive & ar, const BamTools::BamAlignment &some_Balgmt, const unsigned int version)     
        {                                    
            ar  <<   some_Balgmt.Name;               // read name
            ar  <<       some_Balgmt.Length;             // length of query sequence
            ar  <<   some_Balgmt.QueryBases;         // 'original' sequence (as reported from sequencing machine)
            ar  <<   some_Balgmt.AlignedBases;       // 'aligned' sequence (includes any indels, padding, clipping)
            ar  <<   some_Balgmt.Qualities;          // FASTQ qualities (ASCII characters, not numeric values)
            ar  <<   some_Balgmt.TagData;            // tag data (use provided methods to query/modify)
            ar  <<       some_Balgmt.RefID;              // ID number for reference sequence
            ar  <<       some_Balgmt.Position;           // position (0-based) where alignment starts
            ar  <<      some_Balgmt.Bin;                // BAM (standard) index bin number for this alignment
            ar  <<      some_Balgmt.MapQuality;         // mapping quality score
            ar  <<      some_Balgmt.AlignmentFlag;      // alignment bit-flag (use provided methods to query/modify)
//             ar  <<   some_Balgmt.CigarData; // CIGAR operations for this alignment
            ar  <<       some_Balgmt.MateRefID;          // ID number for reference sequence where alignment's mate was aligned
            ar  <<       some_Balgmt.MatePosition;       // position (0-based) where alignment's mate starts
            ar  <<       some_Balgmt.InsertSize;         // mate-pair insert size
//             ar  <<   some_Balgmt.Filename;           // name of BAM file which this alignment comes from            
                    
            const bool some_Balgmt__IsDuplicate  =  some_Balgmt.IsDuplicate();
            const bool some_Balgmt__IsFailedQC  =  some_Balgmt.IsFailedQC();
            const bool some_Balgmt__IsFirstMate  =  some_Balgmt.IsFirstMate();
            const bool some_Balgmt__IsMapped  =  some_Balgmt.IsMapped();
            const bool some_Balgmt__IsMateMapped  =  some_Balgmt.IsMateMapped();
            const bool some_Balgmt__IsMateReverseStrand  =  some_Balgmt.IsMateReverseStrand();
            const bool some_Balgmt__IsPaired  =  some_Balgmt.IsPaired();
            const bool some_Balgmt__IsPrimaryAlignment  =  some_Balgmt.IsPrimaryAlignment();            
            const bool some_Balgmt__IsProperPair  =  some_Balgmt.IsProperPair();
            const bool some_Balgmt__IsReverseStrand  =  some_Balgmt.IsReverseStrand();
            const bool some_Balgmt__IsSecondMate  =  some_Balgmt.IsSecondMate();                    
                    
            ar  <<   some_Balgmt__IsDuplicate;         // returns true if this read is a PCR duplicate
            ar  <<   some_Balgmt__IsFailedQC;          // returns true if this read failed quality control
            ar  <<   some_Balgmt__IsFirstMate;         // returns true if alignment is first mate on read
            ar  <<   some_Balgmt__IsMapped;            // returns true if alignment is mapped
            ar  <<   some_Balgmt__IsMateMapped;        // returns true if alignment's mate is mapped
            ar  <<   some_Balgmt__IsMateReverseStrand; // returns true if alignment's mate mapped to reverse strand
            ar  <<   some_Balgmt__IsPaired;            // returns true if alignment part of paired-end read
            ar  <<   some_Balgmt__IsPrimaryAlignment;  // returns true if reported position is primary alignment
            ar  <<   some_Balgmt__IsProperPair;        // returns true if alignment is part of read that satisfied paired-end resolution
            ar  <<   some_Balgmt__IsReverseStrand;     // returns true if alignment mapped to reverse strand
            ar  <<   some_Balgmt__IsSecondMate;        // returns true if alignment is second mate on read      
            
            //it didn't work doing it directly.  I don't know why...
//             ar  <<   (some_Balgmt.IsDuplicate());         // returns true if this read is a PCR duplicate
//             ar  <<   (some_Balgmt.IsFailedQC());          // returns true if this read failed quality control
//             ar  <<   (some_Balgmt.IsFirstMate());         // returns true if alignment is first mate on read
//             ar  <<   (some_Balgmt.IsMapped());            // returns true if alignment is mapped
//             ar  <<   (some_Balgmt.IsMateMapped());        // returns true if alignment's mate is mapped
//             ar  <<   (some_Balgmt.IsMateReverseStrand()); // returns true if alignment's mate mapped to reverse strand
//             ar  <<   (some_Balgmt.IsPaired());            // returns true if alignment part of paired-end read
//             ar  <<   (some_Balgmt.IsPrimaryAlignment());  // returns true if reported position is primary alignment
//             ar  <<   (some_Balgmt.IsProperPair());        // returns true if alignment is part of read that satisfied paired-end resolution
//             ar  <<   (some_Balgmt.IsReverseStrand());     // returns true if alignment mapped to reverse strand
//             ar  <<   (some_Balgmt.IsSecondMate());        // returns true if alignment is second mate on read
        }//save
        
        
        template<class Archive>
        void load(Archive & ar, BamTools::BamAlignment &some_Balgmt, const uint version)      
        {
            ar  >>   some_Balgmt.Name;               // read name
            ar  >>       some_Balgmt.Length;             // length of query sequence
            ar  >>   some_Balgmt.QueryBases;         // 'original' sequence (as reported from sequencing machine)
            ar  >>   some_Balgmt.AlignedBases;       // 'aligned' sequence (includes any indels, padding, clipping)
            ar  >>   some_Balgmt.Qualities;          // FASTQ qualities (ASCII characters, not numeric values)
            ar  >>   some_Balgmt.TagData;            // tag data (use provided methods to query/modify)
            ar  >>       some_Balgmt.RefID;              // ID number for reference sequence
            ar  >>       some_Balgmt.Position;           // position (0-based) where alignment starts
            ar  >>      some_Balgmt.Bin;                // BAM (standard) index bin number for this alignment
            ar  >>      some_Balgmt.MapQuality;         // mapping quality score
            ar  >>      some_Balgmt.AlignmentFlag;      // alignment bit-flag (use provided methods to query/modify)
//             ar  >>   some_Balgmt.CigarData; // CIGAR operations for this alignment
            ar  >>       some_Balgmt.MateRefID;          // ID number for reference sequence where alignment's mate was aligned
            ar  >>       some_Balgmt.MatePosition;       // position (0-based) where alignment's mate starts
            ar  >>       some_Balgmt.InsertSize;         // mate-pair insert size
//             ar  >>   some_Balgmt.Filename;           // name of BAM file which this alignment comes from                             
            
            bool some_Balgmt__IsDuplicate;
            bool some_Balgmt__IsFailedQC;
            bool some_Balgmt__IsFirstMate;
            bool some_Balgmt__IsMapped;
            bool some_Balgmt__IsMateMapped;
            bool some_Balgmt__IsMateReverseStrand;
            bool some_Balgmt__IsPaired;
            bool some_Balgmt__IsPrimaryAlignment;            
            bool some_Balgmt__IsProperPair;
            bool some_Balgmt__IsReverseStrand;
            bool some_Balgmt__IsSecondMate;
            
            ar  >>   some_Balgmt__IsDuplicate;         // returns true if this read is a PCR duplicate
            ar  >>   some_Balgmt__IsFailedQC;          // returns true if this read failed quality control
            ar  >>   some_Balgmt__IsFirstMate;         // returns true if alignment is first mate on read
            ar  >>   some_Balgmt__IsMapped;            // returns true if alignment is mapped
            ar  >>   some_Balgmt__IsMateMapped;        // returns true if alignment's mate is mapped
            ar  >>   some_Balgmt__IsMateReverseStrand; // returns true if alignment's mate mapped to reverse strand
            ar  >>   some_Balgmt__IsPaired;            // returns true if alignment part of paired-end read
            ar  >>   some_Balgmt__IsPrimaryAlignment;  // returns true if reported position is primary alignment
            ar  >>   some_Balgmt__IsProperPair;        // returns true if alignment is part of read that satisfied paired-end resolution
            ar  >>   some_Balgmt__IsReverseStrand;     // returns true if alignment mapped to reverse strand
            ar  >>   some_Balgmt__IsSecondMate;        // returns true if alignment is second mate on read          
            
            some_Balgmt.SetIsDuplicate(  some_Balgmt__IsDuplicate  );         // sets value of "PCR duplicate" flag
            some_Balgmt.SetIsFailedQC(  some_Balgmt__IsFailedQC  );          // sets value of "failed quality control" flag
            some_Balgmt.SetIsFirstMate(  some_Balgmt__IsFirstMate  );         // sets value of "alignment is first mate" flag
            some_Balgmt.SetIsMapped(  some_Balgmt__IsMapped  );            // sets value of "alignment is mapped" flag
            some_Balgmt.SetIsMateMapped(  some_Balgmt__IsMateMapped  );        // sets value of "alignment's mate is mapped" flag
            some_Balgmt.SetIsMateReverseStrand(  some_Balgmt__IsMateReverseStrand  ); // sets value of "alignment's mate mapped to reverse strand" flag
            some_Balgmt.SetIsPaired(  some_Balgmt__IsPaired  );            // sets value of "alignment part of paired-end read" flag
            some_Balgmt.SetIsPrimaryAlignment(  some_Balgmt__IsPrimaryAlignment  );  // sets value of "position is primary alignment" flag
            some_Balgmt.SetIsProperPair(  some_Balgmt__IsProperPair  );        // sets value of "alignment is part of read that satisfied paired-end resolution" flag
            some_Balgmt.SetIsReverseStrand(  some_Balgmt__IsReverseStrand  );     // sets value of "alignment mapped to reverse strand" flag
            some_Balgmt.SetIsSecondMate(  some_Balgmt__IsSecondMate  );        // sets value of "alignment is second mate on read" flag                                                     
        }//load 
        
        
        template<class Archive>
        inline void serialize( Archive & ar, BamTools::BamAlignment &some_Balgmt, const unsigned int file_version)
        {
            split_free(ar, some_Balgmt, file_version); 
        }             
        
       
       
       
       
       
                

    } // namespace serialization
} // namespace boost








extern type_map_string_to_Readgroup_stats Readgroup_stats;











void print_Balgmt( const BamTools::BamAlignment &some_algmt, std::stringstream *const &some_ss = NULL);





void error_message
	    (const std::stringstream &some_message,
	     const bool &crash_program);
	     
void error_message
	    (const std::string &some_message,
	     const bool &crash_program);	     

void warning_message
	    (const std::stringstream &some_message,
	     const bool &crash_program);



#endif