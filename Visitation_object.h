#ifndef VISITATION_OBJECT_H_INCLUDED
#define VISITATION_OBJECT_H_INCLUDED


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
#include <Sparse_map.h>

#include "Partial_sum_function.h"



class Event;



class Visitation_object
{
  public:
    //constructors:
    Visitation_object() : event_for_elimination(NULL), diploid_compute_size(0)
    { }
    
    Visitation_object( Event *const &in_event_for_elimination,
                       const uint &in_diploid_compute_size,
                       const char *const &in_temporary_PSF_save_filename)
                                                            : event_for_elimination(in_event_for_elimination),
                                                              diploid_compute_size((int)in_diploid_compute_size),
                                                              temporary_PSF_serialization_name(in_temporary_PSF_save_filename)
    { }    
    
    
    //copy constructor
    Visitation_object( const Visitation_object &other_obj)
                                : event_for_elimination(other_obj.event_for_elimination),
                                  diploid_compute_size(other_obj.diploid_compute_size),
                                  permissible_haploid_state_vectors_sparse(other_obj.permissible_haploid_state_vectors_sparse),
                                  neighbors(other_obj.neighbors),
                                  my_Forward_partial_sum_function(other_obj.my_Forward_partial_sum_function),
                                  temporary_PSF_serialization_name(other_obj.temporary_PSF_serialization_name)                                                                    
    { }
    
    
    //destructors
    ~Visitation_object()
    { }
    
    
        
    //variables:
    
    Event *const event_for_elimination;
    int diploid_compute_size;
    
    type_vector_Sparse_map permissible_haploid_state_vectors_sparse;
    //type_list_Sparse_map permissible_haploid_state_vectors_sparse;
    
    type_set_uint neighbors;  //does not include   event_for_elimination , of course.
    
    
    
    Partial_sum_function my_Forward_partial_sum_function;  //saved for backtrace (sampling).
    const std::string temporary_PSF_serialization_name;
    
    Partial_sum_function  my_Backward_partial_sum_function;
    
    
    //a parital sum function  has an m-tuple as an argument (whose value must be permissible), and returns a number (a double).
        //for simplicity, a partial sum function may except a map_uint_to_bool of any size, and it can go through and look for the UIDs that pertain to it, and check that the values of the m-tuple are permissible, and finally return the value of the function.
    
    //As each event is summed out, a new partial sum function is formed in the next Visitation object
  
  //We should have one set of partial_sum_functions.  As each event is summed out,a  new partial sum function is formed.  When summing over an event, if a certain partial sum function is used (that event is prewsent in its argument), then that partial sum function should be eliminated after it is used.
  
    //serialization:
    void save_PSF_to_archive();
    void load_PSF_from_archive();
    
    
    
    void refine_state_vectors_according_to_excluded_Events
				(const type_set_uint &excluded_Events);
  
  
};





//typedef std::vector<Visitation_object>   type_vector_Visitation_object;
typedef std::list<Visitation_object>   type_list_Visitation_object;



#endif