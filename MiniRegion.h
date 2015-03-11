#ifndef MINIREGION_H_INCLUDED
#define MINIREGION_H_INCLUDED


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



#include "globals.h"




class Paired_end_read;
class Event;

struct compare_PERs_by_name;
typedef std::list<Paired_end_read>  type_list_PER;






class MiniRegion
{
  public:
    //constructors:
    MiniRegion()
    { }       
  
//     MiniRegion(const uint &in_chromosome_of_region,
//                const BOOST_Interval &in_region_interval,
//                const uint &in_absolute_position_homologous_to_spawning_Event_profile_position,
//                const type_uint__uint &in_spawning_Event_chromo_and_profileDAR_position,
//                const bool &in_orientation_of_MiniRegion_and_spawning_Event_profile_agree,
//                const bool &in_MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile)
//                 :  chromosome_of_region(in_chromosome_of_region),
//                     region_interval(in_region_interval),
//                     absolute_position_homologous_to_spawning_Event_profile_position(in_absolute_position_homologous_to_spawning_Event_profile_position),
//                     spawning_Event_chromo_and_profileDAR_position(in_spawning_Event_chromo_and_profileDAR_position),
//                     orientation_of_MiniRegion_and_spawning_Event_profile_agree(in_orientation_of_MiniRegion_and_spawning_Event_profile_agree),
//                     MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile(in_MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile),
//                     region_interval__for_testing(in_region_interval)
//     { }

    MiniRegion(const uint &in_chromosome_of_region,
               const BOOST_Interval &in_region_interval,
               const bool &in_orientation_of_MiniRegion_and_spawning_Event_profile_agree,
               const bool &in_MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile)
                :  chromosome_of_region(in_chromosome_of_region),
                    region_interval(in_region_interval),
                    orientation_of_MiniRegion_and_spawning_Event_profile_agree(in_orientation_of_MiniRegion_and_spawning_Event_profile_agree),
                    MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile(in_MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile),
                    MiniID(0)
    { }
    
    
    //copy constructor - unnecessary.  the Implicit copy constructor will work since there are NO POINTERS.    

    //destructor
    ~MiniRegion()
    { }
    
   
    
    //variables:
    uint MiniID;

    uint chromosome_of_region;
    BOOST_Interval region_interval; //closed interval, of course
//     uint absolute_position_homologous_to_spawning_Event_profile_position;
    
//     type_uint__uint spawning_Event_chromo_and_profileDAR_position;
   
    std::string region_sequence;
       
    
    bool orientation_of_MiniRegion_and_spawning_Event_profile_agree;
    bool MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile;
      //when you finally gather a list of Paired-end Reads (PERs) from all MiniRegions spawned by a Variational position, you will need to align all of these PERs to all of the Mini regions.  To do this, you will need to know how to maniupalte the reads of these PERs - reverse and or complement them.  The best way to do this is to record the PERs s.t. they correctly align to the profile of the Variational Position that spawned these MiniRegions.  Then, when aligning a PER to a mini-region, you simply ask the Mini region how it relates to the profile that spawned it.        
    
    
    //functions:   
    
    void print_this_MiniRegion(std::stringstream *const &some_ss = NULL) const;
    
    
    void upload_my_sequence_naive();  //why is it naive?

    
                         
    type_list_PER gather_all_PERs_that_appropriately_intersect_this_MiniRegion_and__format_them_so_that_they_should_align_directly_to_the_spawning_Event_profile
                        (std::map<std::string, Paired_end_read> &already_found_PERs,
                         const std::list<MiniRegion> &current_list_of_homologous_MiniRegions,
                         BamTools::BamReader *const &a_BAM_reader,
                         const Event *const  &spawning_Event) const;
                         
                         
    bool set_region_sequence_from_preloaded_sequence
                                    (const type_map_uint_to_list_BI__string &preloaded_sequences);                         
    
    
    
                                    
                                 
                                    
                                    
                             
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {            
        ar  &  MiniID;
        
        ar  &  chromosome_of_region;
        ar  &  region_interval;
//         ar  &  absolute_position_homologous_to_spawning_Event_profile_position;
        
//         ar  &  spawning_Event_chromo_and_profileDAR_position;
        
        ar  &  region_sequence;
        
        ar  &  orientation_of_MiniRegion_and_spawning_Event_profile_agree;
        ar  &  MiniRegion_must_be_complemented_to_get_to_spawning_Event_profile;
        
//         ar  &  region_interval__for_testing;
    };                                    
                                
          
                                    
                                    
                              
    
    
}; //end class   MiniRegion    
















struct compare_MiniRegions_by_region
{
    bool operator() (const MiniRegion &LHS, const MiniRegion &RHS) const
    {
        if (LHS.chromosome_of_region > RHS.chromosome_of_region)
            return false;
        else if (LHS.chromosome_of_region < RHS.chromosome_of_region)
            return true;
        else   //  ==    
            return  LHS.region_interval.lower() < RHS.region_interval.lower();
    }
};


struct compare_MiniRegion_ptrs_by_MiniID
{
    bool operator() (const MiniRegion *const &LHS, const MiniRegion *const &RHS) const
    {
        return  LHS->MiniID < RHS->MiniID;
    }
};    




typedef std::pair<MiniRegion, MiniRegion>   type_MiniRegion__MiniRegion;
typedef std::set<MiniRegion, compare_MiniRegions_by_region>  type_set_MiniRegion;
typedef std::set<MiniRegion*, compare_MiniRegion_ptrs_by_MiniID>  type_set_MiniRegion_ptr;



typedef std::map<BOOST_Interval, MiniRegion, compare_BI>
                                                type_map_BOOST_Interval__to__MiniRegion;
typedef std::pair<BOOST_Interval, MiniRegion>  type_BOOST_Interval__MiniRegion;


typedef std::map<uint, MiniRegion>   type_map_uint_to_MiniRegion;
typedef std::pair<uint, MiniRegion>  type_uint__MiniRegion;

typedef std::list<MiniRegion>   type_list_MiniRegion;


#endif
