CXX = mpicxx 
CC  = mipcc 
LD  = mpicxx 

# CXX = charmc -language ampi
# CC  = charmc -language ampi
# LD  = charmc -language ampi

.SUFFIXES: .o .f90 .cpp .cu .s

OFLAGS = -O4 -g -m64 -Wall -msse3 -funroll-all-loops -frounding-math
# OFLAGS = -O0 -g -m64 -Wall -msse3 -frounding-math
CFLAGS = $(OFLAGS) -I$(CUDA_TK)/include  -D__$(MACOSX)__  -I/opt/local/include -I$(HOME)/usr/include 
CFLAGS += -I/$(BOOST_INCLUDEDIR)
# CFLAGS += -D_CH2SIG_
CXXFLAGS = $(CFLAGS) 


OFLAGS += -D_FPESIG_ENABLE_ 
OFLAGS += -DCGAL_NDEBUG -DCGAL_HAS_NO_THREADS

PROG = main
VLIB = fvmhd3d

LIBOBJS = distribute_data.o \
					distribute_work.o \
					myMPI.o \
					build_mesh_active.o \
					build_mesh.o \
					refinement.o \
					set_problem.o \
					exchange_data.o \
					sort_local_data.o \
					compute_update.o \
					compute_pvel.o \
					compute_reconstruction.o \
					compute_fluid_grad.o \
					compute_timestep.o \
					iterate.o \
					dump_binary.o read_binary.o

CUOBJS =

LIBS = -L/opt/local/lib -L$(HOME)/usr/lib -lCGAL  
# LIBS += -lefence


all: svnver $(VLIB) $(PROG)

$(PROG): $(PROG).o
	$(LD) $(CXXFLAGS) $(LIBS) $^ -o $@ -L./ -l$(VLIB)

svnver:
	echo "#define VERSION \"0.7-rev`svnversion -n .`\"" > version.h
#	svnversion -n . >> version.h 
#	echo "\"" >> version.h

$(VLIB): $(LIBOBJS) 
	/bin/rm -rf lib$@.a
	ar qv lib$@.a $^
	ranlib lib$@.a



.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

.cpp.s: 
	$(CXX) $(CXXFLAGS) -S $< -o $@

clean:
	/bin/rm -rf *.o version.h
	/bin/rm -rf $(PROG) lib$(VLIB).a

$(PROG).o:  $(VLIB)
$(LIBOBJS): fvmhd3d.h mytimer.h myMPI.h 
set_problem.o: problems/mhd_turbulence.cpp \
	problems/mhd_turbulence_sph.cpp \
	problems/acoustic_waves.cpp problems/alfven_waves.cpp \
	problems/cloud3d.cpp problems/cloud3d_1.cpp \
	problems/tori_collapse.cpp \
	problems/tori_collapse_circ.cpp \
	problems/mri_disk3d.cpp problems/mri_disk3d_1.cpp \
	problems/capture3d.cpp problems/capture3d_1.cpp
distribute_data.o: distributenew.h

