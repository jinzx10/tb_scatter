CXX = g++
CC = g++
#CXXFLAGS = -std=c++14 -O2 -Wall -DARMA_DONT_USE_OPENMP
#LDLIBS = -larmadillo
CXXFLAGS = -std=c++14 -O3 -fopenmp -Wall
LDLIBS = -larmadillo -fopenmp
# implicit compile rule: $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c
# implicit linking rule: $(CC) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) (P.S. LOADLIBES is deprecated)

VPATH = ..

#LDFLAGS = -DNDEBUG
#MAIN = ../test
#OBJS = $(patsubst %.cpp, $(ODIR)/%.o, $(wildcard *.cpp))

# modelconst and main file are private for mono- and diatomic
# the rest are shared
COMMON_OBJS 	= modelconst.o ../mathconst.o ../auxmath.o \
			../GroupOfAtoms.o ../BravLat.o ../LatSupCell.o ../Mol.o ../SupCell.o \
			../BrilZone.o ../LatNeighbor.o ../CplSiteNeighbor.o ../CplOrbNeighbor.o \
			../LatHop.o ../CplHop.o ../Hop.o ../MolOnSite.o \
			../Hamiltonian.o ../LatPot.o ../CplPot.o ../MolPot.o ../Pot.o \
			../LatBandStruc.o ../System.o
EFLD_OBJS 	= ../ElecStruc.o ../Data.o ../EFLD.o
BCME_OBJS 	= ../Diabat.o ../CMEData.o ../BCME.o
ALL_OBJS 	= $(COMMON_OBJS) $(EFLD_OBJS) $(BCME_OBJS) efld_main.o bcme_main.o test.o

efld_main	: efld_main.o $(EFLD_OBJS) $(COMMON_OBJS)
bcme_main 	: bcme_main.o $(BCME_OBJS) $(COMMON_OBJS)
test 		: test.o $(COMMON_OBJS) $(EFLD_OBJS) #$(BCME_OBJS)

JOIN_TREE 	= join.h join.tpp aux.h aux.tpp
ALL_AUX_TREE 	= auxmath.h auxmath.tpp $(JOIN_TREE) mathconst.h
LAT_SUP_TREE 	= LatSupCell.h GroupOfAtoms.h BravLat.h
SUP_CELL_TREE 	= SupCell.h Mol.h $(LAT_SUP_TREE)
ALL_NEIGHBOR 	= LatNeighbor.h CplSiteNeighbor.h CplOrbNeighbor.h
ALL_HOP 	= Hop.h LatHop.h CplHop.h
HOP_TREE 	= $(ALL_HOP) $(ALL_NEIGHBOR) $(SUP_CELL_TREE)
POT_TREE 	= LatPot.h LatNeighbor.h CplSiteNeighbor.h MolPot.h CplPot.h $(SUP_CELL_TREE)
SYSTEM_TREE 	= System.h BrilZone.h Hamiltonian.h MolOnSite.h $(POT_TREE) $(HOP_TREE) \
		$(ALL_AUX_TREE) 

modelconst.o 	: modelconst.h mathconst.h
mathconst.o 	: mathconst.h
auxmath.o 	: $(ALL_AUX_TREE)

GroupOfAtoms.o 	: GroupOfAtoms.h
BravLat.o  	: BravLat.h mathconst.h
LatSupCell.o 	: $(LAT_SUP_TREE) $(ALL_AUX_TREE)
Mol.o 		: Mol.h GroupOfAtoms.h
SupCell.o	: $(SUP_CELL_TREE) $(JOIN_TREE)
BrilZone.o	: BrilZone.h BravLat.h $(ALL_AUX_TREE)

LatNeighbor.o 		: LatNeighbor.h mathconst.h $(ALL_AUX_TREE)
CplSiteNeighbor.o	: CplSiteNeighbor.h $(ALL_AUX_TREE)
CplOrbNeighbor.o	: CplOrbNeighbor.h $(ALL_AUX_TREE)

LatHop.o 	: LatHop.h LatNeighbor.h $(LAT_SUP_TREE) mathconst.h 
CplHop.o 	: CplHop.h CplOrbNeighbor.h $(SUP_CELL_TREE) mathconst.h
Hop.o 		: $(HOP_TREE)

MolOnSite.o 	: MolOnSite.h $(SUP_CELL_TREE) $(ALL_AUX_TREE)
Hamiltonian.o 	: Hamiltonian.h $(HOP_TREE) MolOnSite.h $(JOIN_TREE)
LatBandStruc.o 	: LatBandStruc.h BrilZone.h BravLat.h MolOnSite.h $(SUP_CELL_TREE) $(HOP_TREE)

LatPot.o 	: LatPot.h LatNeighbor.h $(LAT_SUP_TREE)
MolPot.o	: MolPot.h $(SUP_CELL_TREE)
CplPot.o 	: CplPot.h CplSiteNeighbor.h $(SUP_CELL_TREE)
Pot.o 		: Pot.h $(POT_TREE) $(SUP_CELL_TREE)


System.o 	: $(SYSTEM_TREE)

ElecStruc.o 	: ElecStruc.h $(SYSTEM_TREE) $(ALL_AUX_TREE)
Data.o 		: Data.h ElecStruc.h $(SYSTEM_TREE) $(JOIN_TREE)
EFLD.o 		: EFLD.h Data.h ElecStruc.h $(SYSTEM_TREE) $(JOIN_TREE)

Diabat.o 	: Diabat.h $(SYSTEM_TREE) $(ALL_AUX_TREE) 
CMEData.o 	: CMEData.h Diabat.h $(SYSTEM_TREE) $(JOIN_TREE) mathconst.h
BCME.o 		: BCME.h $(SYSTEM_TREE) Diabat.h CMEData.h $(ALL_AUX_TREE)

efld_main.o 	: EFLD.h $(SYSTEM_TREE) Data.h ElecStruc.h $(ALL_AUX_TREE) efld_all.h
bcme_main.o 	: BCME.h $(SYSTEM_TREE) CMEData.h Diabat.h $(ALL_AUX_TREE) bcme_all.h
#test.o 	: EFLD.h Data.h ElecStruc.h BCME.h Diabat.h CMEData.h \
#		$(SYSTEM_TREE) $(ALL_AUX_TREE)
test.o 		: EFLD.h Data.h ElecStruc.h $(SYSTEM_TREE) $(ALL_AUX_TREE)
    			
#$(ODIR)/%.o: $(SDIR)/%.cpp $(SDIR)/*.h $(SDIR)/*.tpp
#	$(CC) -c $< -o $@ $(CFLAGS) $(LIBS)

.PHONY: clean
clean:
	$(RM) $(ALL_OBJS) test efld_main bcme_main
