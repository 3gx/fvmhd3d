OPTS =
# OPTS += -DCMK_OPTIMIZE=1
# OPTS  += -module completion
OPTS += -module GreedyCommLB
# OPTS += -module CkCache
CHARMC = $(CHARMDIR)/bin/charmc $(OPTS) 

CCFLAGS = -DCGAL_NDEBUG -DCGAL_HAS_NO_THREADS -g -O4 -Wall -Wstrict-aliasing=0
CCFLAGS += -m64 -msse3 -funroll-all-loops 
# CCFLAGS += -frounding-math 
CCFLAGS += -ftree-vectorize
CCFLAGS += -ffast-math
CCFLAGS += -I$(HOME)/usr/include -I/$(BOOST_INCLUDEDIR)
LDFLAGS = -L$(HOME)/usr/lib -lCGAL


TARGET = fvmhd3d

default: all
all: $(TARGET)

OBJECTS = fvmhd3d.o System.o \
					loadBalancer.o \
					globalDomains.o \
					moveParticles.o \
					sort_local_data.o \
					globalMesh.o \
					localMesh.o \
					localRefDeref.o \
					computeFluidUpdate.o \
					computeFlux.o \
					computePredictor.o \
					computePvel.o \
					computeReconstruction.o \
					computeTimestep.o \
					slopeLimiter.o \
					Iterate.o \
					ImportFluidData.o \
					IO.o \
					Problem.o 


$(TARGET): $(OBJECTS)
	$(CHARMC) -o $(TARGET) $(LDFLAGS) $(OBJECTS)

.cpp.o:
	$(CHARMC) -c $< -o $@ $(CCFLAGS)

globalMesh.o: globalMesh.cpp
	$(CHARMC) -c $< -o $@ $(CCFLAGS) -frounding-math

localMesh.o: localMesh.cpp
	$(CHARMC) -c $< -o $@ $(CCFLAGS) -frounding-math


.cpp.s:
	$(CHARMC) -S $< -o $@ $(CCFLAGS)

%.decl.h %.def.h :  %.ci
	$(CHARMC) $< 

clean:
	rm -f $(TARGET) *.decl.h *.def.h charmrun $(OBJECTS)

clean_all:
	rm -f $(TARGET) *.decl.h *.def.h charmrun $(OBJECTS) *~

$(TARGET).o: $(OBJECS)
$(OBJECTS): fvmhd3d.ci fvmhd3d.h MeshPoint.h Fluid.h Scheduler.h vector3.h Boundary.h bOctree.h
bOctree.h: vector3.h Boundary.h memory_pool.h
fvmhd3d.h: fvmhd3d.decl.h fvmhd3d.def.h
System.o: distributenew.h
Problem.o: Problem.cpp \
	mri_disk.cpp \
	capture3d.cpp \
	mri_cyl.cpp \
	acoustic.cpp \
	alfven.cpp \
	blast.cpp

