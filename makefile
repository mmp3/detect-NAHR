MPICXX=mpic++
#CXX = g++
#CXX = gfilt /cand:L
CXXFLAGS=   -O3   -fopenmp   -std=c++0x  $(INCLUDES)  -Werror=return-type		\

INCLUDES=  -I$(GENERAL)    -I$(BOOST_INCLUDE)    -I$(BAMTOOLS)/include   -I$(MPFR)/include    -I$(MPFRCXX)	-I$(MPI_INCLUDE)  -I$(GMP)/include	-I$(GMPFRXX)  -I$(CUPA)	\

LIBS=  $(GMPFRXX)/libgmpfrxx.a	\
		-L$(BOOST_LIB) -ldl  -Wl,-rpath,$(BOOST_LIB)	       -lboost_serialization   -lboost_mpi		-lboost_filesystem	-lboost_system	\
		-L$(MPFR)/lib -ldl -Wl,-rpath,$(MPFR)/lib    -lmpfr						\
		-L$(GMP)/lib -ldl -Wl,-rpath,$(GMP)/lib      -lgmpxx    -lgmp					\
		-L$(BAMTOOLS)/lib -ldl -Wl,-rpath,$(BAMTOOLS)/lib      -lbamtools				\
		-L$(MPI_LIB) -ldl  -Wl,-rpath,$(MPI_LIB)	$(MPI_LIB_NAME_FLAG)				\
		-Wl,-rpath,${LD_RUN_PATH} 	-Wl,-rpath,$(LD_LIBRARY_PATH)	\
		-Wl,-rpath,${MPI_LIB}

#		-L$(BAMTOOLS)/lib -ldl -Wl,-rpath,$(BAMTOOLS)/lib      -lbamtools						\
#Recall:
# -L     gives the path to the libraries you want to link to.
# -l     is the library you want to link to.  MAKE will search through the directories provided via the -L option

#  -std=c++0x  

LIBCUDA=-I$(CUDA_ROOT_DIR)/include   -L${CUDA_ROOT_DIR}/lib64  -Wl,-rpath,$(CUDA__ROOT_DIR)   -lcudart  -lc -lm
OBJCUPA=  $(CUPA)/CUPA.o  $(CUPA)/CUPA.a  $(CUPA)/dlink.o 

MPI_LIB=$(MVAPICH2_DIR)/lib
MPI_INCLUDE=$(MVAPICH2_DIR)/include
MPI_LIB_NAME_FLAG=-lmpichcxx  -lmpich

#MPI_LIB=$(OPENMPI_DIR)/lib
#MPI_INCLUDE=$(OPENMPI_DIR)/include
#MPI_LIB_NAME_FLAG=-lmpi_cxx


BOOST_INCLUDE=${BOOST_ROOT}
BOOST_LIB=${BOOST_ROOT}

#BOOST_INCLUDE=/gpfs/runtime/opt/boost/1.52.0/src/boost_1_52_0
#BOOST_LIB=/gpfs/runtime/opt/boost/1.52.0/lib
#BOOST_INCLUDE=/gpfs/runtime/opt/boost/1.52.0-0x/src/boost_1_52_0
#BOOST_LIB=/gpfs/runtime/opt/boost/1.52.0-0x/lib
#BOOST_INCLUDE=$(BOOST)/include
#BOOST_LIB=$(BOOST/lib





vpath templates.%				$(GENERAL)
vpath Coords_and_homologous_region.%		$(GENERAL)
vpath translations_through_alignment.%		$(GENERAL)
vpath SD_Entry_general.%			$(GENERAL)
vpath Interlocking_entry.%			$(GENERAL)
vpath Sparse_map.%				$(GENERAL)
vpath general_typedefs.%			$(GENERAL)
vpath Readgroup_statistics.%			$(GENERAL)






# 			/gpfs/main/research/compbio/gpfs/home/mmparks/MPFRC++/mpreal.o \
# 			/libmpfr.a \
# 			/libgmp.a \


#		Breakpoint_iterator.o

OBJECTS =	main.o			\
		globals.o		\
		other_functions.o	\
		io_functions.o		\
		MiniRegion.o		\
		MiniEvent.o 		\
		Conn_Comp.o		\
		Event.o			\
		Partial_sum_function.o	\
		Paired_end_read.o	\
		PER_aligner.o		\
		Visitation_object.o	\
		Star_alignment.o	\
		Call_set.o		\
		PER_bundle_wrapper.o	\
		post.o


OBJECTS_GENERAL =	$(GENERAL)/Coords_and_homologous_region.o		\
			$(GENERAL)/translations_through_alignment.o		\
			$(GENERAL)/SD_Entry_general.o			\
			$(GENERAL)/Sparse_map.o				\
			$(GENERAL)/general_typedefs.o			\
			$(GENERAL)/Readgroup_statistics.o




# The special rule .PHONY tells Make which targets are not files. This avoids conflict with files of the same name, and improves performance.
# If a phony target is included as a prerequisite for another target, it will be run every time that other target is required. Phony targets are never up-to-date.
#     from http://linuxdevcenter.com/pub/a/linux/2002/01/31/make_intro.html?page=2
.PHONY : make_general




all : make_general   $(OBJECTS)  
	$(MPICXX)	$(CXXFLAGS) -o detect_NAHR 	$(OBJECTS)   $(OBJECTS_GENERAL)  $(OBJCUPA) $(LIBS)  $(LIBCUDA)	
	@echo "detect_NAHR executable successfully created."

#	$(MPFRCXX)/mpreal.o											\




make_general :
	$(MAKE) -C $(GENERAL)
	#  "$(MAKE)"  is a variable always defined by the "make" utility.  "-C" indicates that there will be a directroy change, to "$(GENERAL)"
	#               from    http://stackoverflow.com/questions/3494999/subdirectories-and-makefiles





main.o :	main.cpp \
		general_typedefs.h    templates.h	Readgroup_statistics.h		 \
		Conn_Comp.h  Event.h  io_functions.h  globals.h  Visitation_object.h    other_functions.h  MiniEvent.h  post.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@



globals.o :	globals.cpp   globals.h \
		general_typedefs.h		Readgroup_statistics.h	Sparse_map.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@


Partial_sum_function.o :	Partial_sum_function.cpp Partial_sum_function.h \
				general_typedefs.h   Sparse_map.h    templates.h   translations_through_alignment.h   \
				globals.h     Event.h  other_functions.h     Conn_Comp.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@


other_functions.o :	other_functions.cpp   other_functions.h \
			general_typedefs.h  translations_through_alignment.h \
			globals.h   MiniRegion.h   Event.h    Paired_end_read.h		\
			Conn_Comp.h	io_functions.h	Call_set.h   Visitation_object.h	Star_alignment.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@


io_functions.o :	io_functions.cpp io_functions.h \
			SD_Entry_general.h   general_typedefs.h  templates.h  Interlocking_entry.h  Coords_and_homologous_region.h  translations_through_alignment.h  Sparse_map.h \
			globals.h  Conn_Comp.h  Visitation_object.h  Event.h	Readgroup_statistics.h	Call_set.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

Conn_Comp.o :	Conn_Comp.cpp  Conn_Comp.h    \
		general_typedefs.h  SD_Entry_general.h  Interlocking_entry.h  templates.h  Coords_and_homologous_region.h  Sparse_map.h   translations_through_alignment.h  Readgroup_statistics.h  \
		globals.h  Event.h  Partial_sum_function.h  other_functions.h  io_functions.h  Visitation_object.h	\
		   Paired_end_read.h   MiniRegion.h     Star_alignment.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@


Event.o :	Event.cpp  Event.h \
		general_typedefs.h  SD_Entry_general.h  Interlocking_entry.h  Coords_and_homologous_region.h  translations_through_alignment.h  Sparse_map.h  templates.h    Readgroup_statistics.h  \
		globals.h  Conn_Comp.h  other_functions.h  Paired_end_read.h  MiniRegion.h       io_functions.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

	



Paired_end_read.o :	Paired_end_read.cpp  Paired_end_read.h \
			Sparse_map.h  general_typedefs.h  SD_Entry_general.h  Interlocking_entry.h  templates.h  translations_through_alignment.h 	Readgroup_statistics.h   \
			globals.h  other_functions.h  Event.h  MiniRegion.h   PER_aligner.h    io_functions.h    Conn_Comp.h	Star_alignment.h   MiniEvent.h   PER_bundle_wrapper.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@




MiniRegion.o :	MiniRegion.cpp   MiniRegion.h \
		general_typedefs.h  SD_Entry_general.h  templates.h  translations_through_alignment.h   \
		globals.h  Paired_end_read.h   io_functions.h  other_functions.h  Event.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

MiniEvent.o :	MiniEvent.cpp  \
		general_typedefs.h  SD_Entry_general.h  templates.h   translations_through_alignment.h   \
		globals.h  Paired_end_read.h  other_functions.h  Event.h  MiniRegion.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@


PER_aligner.o :	PER_aligner.cpp   PER_aligner.h \
		general_typedefs.h  templates.h   translations_through_alignment.h  \
		other_functions.h     Paired_end_read.h   globals.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@



Visitation_object.o : 	Visitation_object.cpp    Visitation_object.h    \
			general_typedefs.h     Sparse_map.h	\
			globals.h    Event.h   Partial_sum_function.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@


Star_alignment.o :	Star_alignment.cpp     Star_alignment.h		\
			general_typedefs.h  templates.h   translations_through_alignment.h  \
			globals.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@


Call_set.o :	Call_set.cpp	Call_set.h	\
		general_typedefs.h		\
		globals.h
	$(CXX) $(CXXFLAGS) -c $< -o $@



PER_bundle_wrapper :	PER_bundle_wrapper.cpp    PER_bundle_wrapper.h	\
				general_typedefs.h	Event.h		other_functions.h	globals.h	MiniRegion.h	Paired_end_read.h	MiniEvent.h	\
				CUPA.h
	$(CXX) $(CXXFLAGS) -c $< -o $@



post.o :	post.cpp  post.h	\
		general_typedefs.h    translations_through_alignment.h      templates.h		\
		Event.h		other_functions.h	globals.h	Conn_Comp.h	io_functions.h	Call_set.h	
	$(CXX) $(CXXFLAGS) -c $< -o $@






clean :
	$(MAKE) -C $(GENERAL) clean
	rm detect_NAHR $(OBJECTS)
	echo "detect_NAHR object files cleaned out."













# # CXX = gfilt
# CXX = g++
# CXXFLAGS = -g -O3 -I$(GENERAL) -I$(BOOST_INCLUDE) -I$(BAMTOOLS)/include -L$(BAMTOOLS)/lib -ldl -Wl,-rpath,$(BAMTOOLS)/lib -lbamtools
# # -std=gnu++0x    <---  if you wantred to use std::unordered_map  (i.e. a hash table
# 
# BOOST_INCLUDE = /gpfs/main/research/compbio/gpfs/home/mmparks/boost/include
# GENERAL = /gpfs/main/research/compbio/gpfs/home/mmparks/NAHR_finder/general_code
# 
# BAM_INCLUDE = /gpfs/main/research/compbio/gpfs/home/mmparks/BAMtools/include
# BAM_LIB = /gpfs/main/research/compbio/gpfs/home/mmparks/BAMtools/lib
# 
# 
# vpath templates.% $(GENERAL)
# vpath Collection_of_intervals.% $(GENERAL)
# vpath Coords_and_homologous_region.% $(GENERAL)
# vpath translations_through_alignment.% $(GENERAL)
# vpath SD_Entry_general.% $(GENERAL)
# vpath Interlocking_entry.% $(GENERAL)
# vpath Sparse_map.% $(GENERAL)
# vpath general_typedefs.% $(GENERAL)
# 
# 
# 
# SOURCES = 	main.cpp\
# 		globals.cpp other_functions.cpp io_functions.cpp\
# 		Conn_Comp.cpp Event.cpp\
# 		Partial_sum_function.cpp\
# 		Paired_end_read.cpp Breakpoint_iterator.cpp\
# 		MiniRegion.cpp MiniEvent.cpp\
# 		$(GENERAL)/templates.cpp $(GENERAL)/Collection_of_intervals.cpp\
# 		$(GENERAL)/Coords_and_homologous_region.cpp $(GENERAL)/translations_through_alignment.cpp\
# 		$(GENERAL)/SD_Entry_general.cpp $(GENERAL)/Sparse_map.cpp
# 
# HEADERS = 	globals.h other_functions.h io_functions.h\
# 		Conn_Comp.h Event.h\
# 		Partial_sum_function.h Visitation_object.h\
# 		Paired_end_read.h  \
# 		MiniRegion.h\
# 		Point_or_interval.h\
# 		templates.h Collection_of_intervals.h\
# 		Coords_and_homologous_region.h translations_through_alignment.h\
# 		SD_Entry_general.h general_typedefs.h\
# 		Interlocking_entry.h Sparse_map.h
# 
# OBJECTS = $(SOURCES:.cpp=.o)
# 
# 
# all : $(OBJECTS)
# 	$(MPICXX) $(CXXFLAGS) $(OBJECTS) -o detect_NAHR
# 	@echo "detect_NAHR executable created."
# 
# %.o : %.cpp $(HEADERS)
# 	$(MPICXX) $(CXXFLAGS) -c $< -o $@
# 
# 
# clean :
# 	rm detect_NAHR $(OBJECTS)
# 	echo "detect_NAHR object files cleaned out."














# CXX = g++
# CXXFLAGS = -g -O3 -I$(BAMTOOLS)/include -L$(BAMTOOLS)/lib -ldl -Wl,-rpath,$(BAMTOOLS)/lib -lbamtools
# OBJECTS = main.o Event.o Conn_Comp.o other_functions.o io_functions.o templates.o
# # globals.o
# BAM_INCLUDE = /gpfs/main/research/compbio/gpfs/home/mmparks/BAMtools/include
# BAM_LIB = /gpfs/main/research/compbio/gpfs/home/mmparks/BAMtools/lib
# 
# NAHR_model : $(OBJECTS)
# 	$(MPICXX) $(CXXFLAGS) $(OBJECTS) -o NAHR_model
# 	@echo "NAHR_model executable created."
# 
# 
# main.o : main.cpp Conn_Comp_and_Event.h other_functions.h templates.h
# 	$(MPICXX) $(CXXFLAGS) -c main.cpp -o $@
# 
# Event.o : Event.cpp Conn_Comp_and_Event.h templates.h
# 	$(MPICXX) $(CXXFLAGS) -c Event.cpp -o $@
# 
# Conn_Comp.o : Conn_Comp.cpp Conn_Comp_and_Event.h templates.h
# 	$(MPICXX) $(CXXFLAGS) -c Conn_Comp.cpp -o $@
# 
# other_functions.o : other_functions.cpp other_functions.h Conn_Comp_and_Event.h
# 	$(MPICXX) $(CXXFLAGS) -c other_functions.cpp -o $@
# 
# io_functions.o : io_functions.cpp Conn_Comp_and_Event.h
# 	$(MPICXX) $(CXXFLAGS) -c io_functions.cpp -o $@
# 
# templates.o : templates.cpp templates.h
# 	$(MPICXX) $(CXXFLAGS) -c templates.cpp -o $@
# 
# # globals.o : global.h global.cpp
# # 	$(MPICXX) $(CXXFLAGS) -c global.cpp -o $@
# 
# clean :
# 	rm NAHR_model $(OBJECTS)
