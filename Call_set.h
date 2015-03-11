#ifndef CALL_SET_H_INCLUDED
#define CALL_SET_H_INCLUDED


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








                    
                    
class Call_set
{
    public:
        Call_set()  :  potentially_gene_conversion(false)
        {  }
        
        
        Call_set (  const std::string  &in_genome_name,
		    const std::string  &in_pop,
                   const std::string  &in_var_outcome,
                   const longdouble  &in_P_outcome,
                   const std::string &in_chr_and_brkpt_begin,
                   const std::string &in_chr_and_brkpt_end,
                   const longdouble &in_P_brkpt,
                   const uint &in_cc,
                   const uint &in_ev,
		    const uint &in_brkpt_index)
                                            :    genome_name(in_genome_name),
						population_of_genome(in_pop),
                                                var_outcome(in_var_outcome),
                                                P_outcome(in_P_outcome),
                                                P_brkpt(in_P_brkpt),
                                                cc_of_event(in_cc),
                                                event_of_var(in_ev),
                                                potentially_gene_conversion(false),
                                                brkpt_index(brkpt_index)
        {
            char *chr_char = new char[option_size];
	    
            uint brkpt_lower, brkpt_upper;
            std::sscanf(  in_chr_and_brkpt_begin.c_str(), "chr %[^:]: %u",
                                               chr_char,
                                               &brkpt_lower);
                       
            std::sscanf(  in_chr_and_brkpt_end.c_str(), "chr %*[^:]: %u",
                                               &brkpt_upper);                                               
            
            if (strcasecmp(chr_char, "X") == 0)                                   
                chromo_of_call = 23;            
            else if (strcasecmp(chr_char, "Y") == 0)                 
                chromo_of_call = 24; 
            else
                chromo_of_call = (uint)atoi(chr_char);
                                   
            brkpts.set(  brkpt_lower,  brkpt_upper );     
	    
	    delete[] chr_char;
        }
        
        
        
        //variables:        
        
        std::string genome_name;
	std::string population_of_genome;
        std::string var_outcome;
        
        uint chromo_of_call;
        BOOST_Interval brkpts;
	uint brkpt_index;
        
        longdouble P_outcome;        
        longdouble P_brkpt;
        
        uint cc_of_event;
        uint event_of_var;
            
        
        bool potentially_gene_conversion;
        
        
        
        
        
        //functions
        
        bool operator<(const Call_set &RHS) const
        {                 
            if (genome_name.compare(RHS.genome_name) < 0)
                return true;
            else if (genome_name.compare(RHS.genome_name) > 0)
                return false;
            
            else if (chromo_of_call < RHS.chromo_of_call)
                return true;
            else if (chromo_of_call > RHS.chromo_of_call)
                return false;
            else
            {
                compare_BI  my_compare;
                if ( my_compare(brkpts, RHS.brkpts) )
                    return true;
                else if ( my_compare(RHS.brkpts, brkpts) )
                    return false;                            
                
                else if (var_outcome.compare(RHS.var_outcome) < 0)
                    return true;
                else if (var_outcome.compare(RHS.var_outcome) > 0)
                    return false;
                
                else if (cc_of_event < RHS.cc_of_event)
                    return true;                
                else if (cc_of_event > RHS.cc_of_event)
                    return false; 
                
                else if (event_of_var < RHS.event_of_var)
                    return true;                
                else if (event_of_var > RHS.event_of_var)
                    return false;                   
                
                else if (P_outcome < RHS.P_outcome)
                    return true;                
                else if (P_outcome > RHS.P_outcome)
                    return false;                    
                
                else if (P_brkpt < RHS.P_brkpt)
                    return true;                
                else if (P_brkpt > RHS.P_brkpt)
                    return false;                  
                
                else if ( potentially_gene_conversion  and  !RHS.potentially_gene_conversion )
                    return false;
                else if ( !potentially_gene_conversion  and  RHS.potentially_gene_conversion )
                    return true;  
                
                else
                    return false;
            }
            
                
            
            
        }//operator<
                    
        
        

        bool  check_overlap_with_another_Call_set_by_chromo_and_brkpts
                                            (const Call_set &RHS)            const;    
					    
					    
	void print_this_Call( std::stringstream *const &some_ss = NULL)  const;
        
        
        
        
};
                    
      
      
      
typedef   std::list<Call_set>  type_list_Call_set;      
      
typedef  std::map<std::string, type_list_Call_set>   type_map_string_to_list_Call_set;
                    
                    
                    
                    
#endif