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
#include <boost/filesystem.hpp>


#include <mpreal.h>
#include <mpfr.h>



#include <api/BamAlignment.h>



#include <general_typedefs.h>
#include <translations_through_alignment.h>
#include <templates.h>


#include "globals.h"
#include "Call_set.h"










bool  Call_set::check_overlap_with_another_Call_set_by_chromo_and_brkpts  
                            (const Call_set &RHS)    const
{ 
    if (chromo_of_call != RHS.chromo_of_call)
        return false;
    else 
        return BOOST_overlap(brkpts, RHS.brkpts);   

} // compare_two_Call_sets_by_chromo_and_brkpts












void  Call_set::print_this_Call( std::stringstream *const &some_ss)  const
{
    
    std::stringstream output_str;
    
    bool provided_with_valid_sstream   =   (some_ss != NULL) ;    
    std::stringstream *const ss_ptr  =   provided_with_valid_sstream   ?      some_ss    :    &output_str;  
    
    print_line_of_markers("[[[", ss_ptr);    
    
    (*ss_ptr) << "\n\n\n\n\nprinting CALL:\n\n";
    
    (*ss_ptr) << "genome_name = [" << genome_name << "]\n\n";
    (*ss_ptr) << "var_outcome = [" << var_outcome << "]\n\n";
    (*ss_ptr) << "chromo_of_call = [" << chromo_of_call << "]\n\n";
    (*ss_ptr) << "brkpts = [" << brkpts << "]\n\n";
    (*ss_ptr) << "P_outcome = [" << P_outcome << "]\n\n";
    (*ss_ptr) << "P_brkpt = [" << P_brkpt << "]\n\n";
    (*ss_ptr) << "cc_of_event = [" << cc_of_event << "]\n\n";
    (*ss_ptr) << "event_of_var = [" << event_of_var << "]\n\n";    
    (*ss_ptr) << "potentially_gene_conversion = [" << potentially_gene_conversion << "]\n\n";          
   
    print_line_of_markers("]]]", ss_ptr);    
    

    if (!provided_with_valid_sstream)
        std::fprintf(stderr, "\n\n\n%s\n\n\n", (*ss_ptr).str().c_str() );    
        
    std::cerr.flush();std::fflush(stderr);

} // print_this_Call


