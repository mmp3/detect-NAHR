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

#include "globals.h"
#include "Visitation_object.h"
#include "Event.h"








void Visitation_object::save_PSF_to_archive()
{        
    // #pragma omp single !!!!!!!!!!!!!!!!!!!!!!
    
    std::ofstream out_fs( temporary_PSF_serialization_name.c_str(),   std::ios::binary | std::ios::out  );  //  std::ios_base::binary     
                        if ( !out_fs.good()  )
                        {
                            std::stringstream error_strm;
                            print_line_of_markers("ERROR! ", &error_strm);
                            print_line_of_markers("(", &error_strm);
                            error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                            
                            error_strm << "\n\nERROR!   not good ofstream   in  \"save_PSF_to_archive\",\n\t\t temporary_PSF_serialization_name =    "  
                            <<  temporary_PSF_serialization_name  << "\n\n";
                            
                            print_line_of_markers(")", &error_strm);
                            std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
                        }    
    
        
    boost::archive::binary_oarchive out_archive( out_fs );
    
    out_archive << my_Forward_partial_sum_function;    
    
    
    out_fs.close();

}  //  save_PSF_to_archive






void Visitation_object::load_PSF_from_archive()
{

    my_Forward_partial_sum_function.clear_data();    
    
    
    // open the archive
    std::ifstream in_fs( temporary_PSF_serialization_name.c_str(),   std::ios::binary | std::ios::in );  //  std::ios_base::binary    
                            if ( !in_fs.good()  )
                            {
                                std::stringstream error_strm;
                                print_line_of_markers("ERROR! ", &error_strm);
                                print_line_of_markers("(", &error_strm);
                                error_strm << "\n\nrank " << my_MPI_rank << "  thread " << omp_get_thread_num() << "\n";
                                
                                error_strm << "\n\nERROR!   not good ifstream   in  \"load_PSF_from_archive\",\n\t\t temporary_PSF_serialization_name =    "  
                                <<  temporary_PSF_serialization_name  << "\n\n";
                                
                                print_line_of_markers(")", &error_strm);
                                std::fprintf(stderr, "\n\n%s\n\n", error_strm.str().c_str() );    
                            }
    
    
    boost::archive::binary_iarchive in_archive( in_fs );

    // restore  from the archive
    
    in_archive >> my_Forward_partial_sum_function;
    
//     //recall that we ignored the Event pointer when saving.  so it must be reset now.
//     my_Forward_partial_sum_function.eliminating_Event_UID = event_for_elimination->UID;
    
    
    
    in_fs.close();

}  //  load_PSF_from_archive






















void  Visitation_object::refine_state_vectors_according_to_excluded_Events
				(const type_set_uint &excluded_Events)
{
                        
    type_vector_Sparse_map refined_state_vectors;
    const uint naive_number_of_state_vectors = permissible_haploid_state_vectors_sparse.size();
    
    
    for (uint j=0; j<naive_number_of_state_vectors; ++j)
    {
	bool conflicts = false;
	for (type_set_uint::const_iterator it_exc = excluded_Events.begin();
		it_exc != excluded_Events.end();
		++it_exc) 
	{
	    if (permissible_haploid_state_vectors_sparse.at(j).get_state_for_UID(*it_exc) != 0)
	    {
		conflicts = true;
		break;
	    }
	}
	
	if (!conflicts)
	    refined_state_vectors.push_back(  permissible_haploid_state_vectors_sparse.at(j)  );            
    }
    
    

    if (my_MPI_rank == 0)
    {
	std::fprintf(stderr, "\n\n\t\tnumber of (haploid)state vectors excluded: %u",
			    permissible_haploid_state_vectors_sparse.size() - refined_state_vectors.size()   );
    }
			    
			    
    permissible_haploid_state_vectors_sparse = refined_state_vectors;
    diploid_compute_size = (uint)std::pow( permissible_haploid_state_vectors_sparse.size(), 2 ); 

} // refine_state_vectors_according_to_excluded_Events			