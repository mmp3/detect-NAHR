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

#include <general_typedefs.h>
#include <SD_Entry_general.h>
#include <templates.h>
#include <translations_through_alignment.h>

#include "globals.h"
#include "MiniRegion.h"
#include "io_functions.h"
#include "other_functions.h"
#include "Event.h"









void MiniRegion::print_this_MiniRegion(std::stringstream *const &some_ss) const
{

    std::cerr.flush();std::fflush(stderr);
    
    std::stringstream output_str;
    
    const bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str;  
    
    
     
    (*ss_ptr) << "\n\nprinting MiniRegion...\n\t\t\t\tMiniID = " << MiniID
              << "\n\t\t\t\tchromosome_of_region " << chromosome_of_region 
              << "\n\t\t\t\tregion_interval = ["<<  region_interval.lower()  <<", "<< region_interval.upper() << "]"
              <<"\n\t\t\t\torientation_of_MiniRegion_and_spawning_Event_profile_agree = "<<  orientation_of_MiniRegion_and_spawning_Event_profile_agree
              <<"\n\t\t\t\tMiniRegion_must_be_complemented_to_get_to_spawning_Event_profile = "<< MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile
              <<"\n\t\t\t\tregion_sequence:\n" << region_sequence   <<"\n";     
               
    
    
    if (!provided_with_valid_sstream)
    {  std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );          }
    std::cerr.flush();std::fflush(stderr);

}  // end   of   print_this_Mini













void MiniRegion::upload_my_sequence_naive()  //why is it naive?
{
  
    region_sequence = get_specific_DNA_sequence_from_Reference_genome(chromosome_of_region, region_interval);                                                   

}  // end of   upload_my_sequence_naive












bool MiniRegion::set_region_sequence_from_preloaded_sequence
                                    (const type_map_uint_to_list_BI__string &preloaded_sequences) 
{
    const type_map_uint_to_list_BI__string::const_iterator  it_uploaded_chromo   =   preloaded_sequences.find( chromosome_of_region );
  
    for(type_list__BI__string::const_iterator it_preload = it_uploaded_chromo->second.begin();
            it_preload != it_uploaded_chromo->second.end();
            ++it_preload)
    {
        if ( BOOST_subset(region_interval, it_preload->first) )
        {            
            region_sequence = it_preload->second.substr(region_interval.lower() - it_preload->first.lower(), BOOST_width_inclusive(region_interval)  );        
            
            if ( region_sequence.size()  != BOOST_width_inclusive(region_interval)  )
            {
                std::stringstream error_strm;
                print_line_of_markers("ERROR! ", &error_strm);
                print_line_of_markers("(", &error_strm);
                error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                
                error_strm << "\n\nERROR!\n\tregion_sequence.size()  = " << region_sequence.size()  
                            << "\n\tBOOST_width_inclusive(region_interval)   "  << BOOST_width_inclusive(region_interval)  
                            << "\n(after normal upload).\n\n";
                
                print_this_MiniRegion(  &error_strm  );            
                            
                print_line_of_markers(")", &error_strm);
                std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );            
            
            }
            
            return true;            
        }     
    }
        
        
        
        
    std::stringstream error_strm;    
    print_line_of_markers("ERROR! ", &error_strm);
    print_line_of_markers("(", &error_strm);
    error_strm << "\nERROR!   unable to find MiniRegion sequence in pre-uploaded strings!!!\n\n";
    error_strm << "it_uploaded_chromo->second.size() = " << it_uploaded_chromo->second.size() << "\n\n";
    for (type_map_uint_to_list_BI__string::const_iterator it_chr = preloaded_sequences.begin();
            it_chr != preloaded_sequences.end();
            ++it_chr)             
    {
        error_strm << "\n\n\tchr " << it_chr->first << ":\t\t\n";
        for(type_list__BI__string::const_iterator it_preload = it_chr->second.begin();
                it_preload != it_chr->second.end();
                ++it_preload)
            error_strm << "BI:  " << it_preload->first << "\t\tstring size =\t" << it_preload->second.size()  << "\n";
    }
                 
    print_this_MiniRegion(&error_strm);
    print_line_of_markers(")", &error_strm);   
    std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );
    
    upload_my_sequence_naive();
    
    
    if ( region_sequence.size()  != BOOST_width_inclusive(region_interval)  )
    {
        std::stringstream error_strm;
        print_line_of_markers("ERROR! ", &error_strm);
        print_line_of_markers("(", &error_strm);
        error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
        
        error_strm << "\n\nERROR!\n\tregion_sequence.size()  = " << region_sequence.size()  
                    << "\n\tBOOST_width_inclusive(region_interval)   "  << BOOST_width_inclusive(region_interval) 
                    << "\n\n(after  \"upload_my_sequence_naive).\"\n\n";
        
        print_this_MiniRegion(  &error_strm  );            
                    
        print_line_of_markers(")", &error_strm);
        std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
    }    
    
    return false;


}  // end of   upload_my_sequence_naive





