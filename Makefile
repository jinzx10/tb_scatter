CXX = mpic++
CC = mpic++
CXXFLAGS = -std=c++14 -O3 -Wall -DARMA_DONT_USE_OPENMP
LDLIBS = -larmadillo
# implicit compile rule: $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c
# implicit linking rule: $(CC) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) (P.S. LOADLIBES is deprecated)

VPATH = src

#OBJS = $(patsubst %.cpp, $(ODIR)/%.o, $(wildcard *.cpp))

OBJS 	= Atoms.o BravLat.o LatSupCell.o SupCell.o BrilZone.o CplHop.o LatPot.o LatTB.o \
	  LatNeighbor.o CplSiteNeighbor.o CplOrbNeighbor.o modelconst.o mathconst.o auxmath.o
DEPS 	= System.h System.tpp modelconst.h SupCell.h BravLat.h BrilZone.h LatSupCell.h Atoms.h \
	  LatTB.h CplHop.h LatPot.h LatNeighbor.h CplOrbNeighbor.h CplSiteNeighbor.h \
	  $(ALL_AUX_TREE)

NO_Au_efld	: main_NO_Au_efld.o NO_Au_Diabats.o $(OBJS)
	mkdir -p obj
	$(CXX) $^ $(LDLIBS) -o $@
	mv *.o obj
NO_Au_cme 	: main_NO_Au_cme.o NO_Au_Diabats.o $(OBJS)
	mkdir -p obj
	$(CXX) $^ $(LDLIBS) -o $@
	mv *.o obj
NO_Au_test 	: test_NO_Au.o NO_Au_Diabats.o $(OBJS)
	mkdir -p obj
	$(CXX) $^ $(LDLIBS) -o $@
	mv *.o obj

JOIN_TREE 	= join.h join.tpp aux.h aux.tpp
ALL_AUX_TREE 	= auxmath.h auxmath.tpp $(JOIN_TREE) mathconst.h

modelconst.o 	: modelconst.h mathconst.h
mathconst.o 	: mathconst.h
auxmath.o 	: $(ALL_AUX_TREE)

Atoms.o 	: Atoms.h $(JOIN_TREE) #
BravLat.o  	: BravLat.h $(JOIN_TREE) mathconst.h #
LatSupCell.o 	: LatSupCell.h Atoms.h BravLat.h $(ALL_AUX_TREE) #
SupCell.o	: SupCell.h LatSupCell.h Atoms.h BravLat.h $(JOIN_TREE) mathconst.h #
BrilZone.o	: BrilZone.h BravLat.h $(ALL_AUX_TREE) #

LatNeighbor.o 		: LatNeighbor.h mathconst.h $(ALL_AUX_TREE)
CplSiteNeighbor.o	: CplSiteNeighbor.h $(ALL_AUX_TREE)
CplOrbNeighbor.o	: CplOrbNeighbor.h $(ALL_AUX_TREE)

LatTB.o 	: LatTB.h LatNeighbor.h $(JOIN_TREE) mathconst.h
CplHop.o 	: CplHop.h CplOrbNeighbor.h mathconst.h
LatPot.o 	: LatPot.h LatNeighbor.h
NO_Au_Diabats.o : NO_Au_Diabats.h CplSiteNeighbor.h $(ALL_AUX_TREE)

main_NO_Au_efld.o 	: EFLD.h EFLD.tpp common.h NO_Au_Diabats.h $(DEPS)
main_NO_Au_cme.o 	: CME.h CME.tpp common.h NO_Au_Diabats.h $(DEPS)
test_NO_Au.o 		: CME.h CME.tpp EFLD.h EFLD.tpp common.h NO_Au_Diabats.h $(DEPS)
    			

.PHONY: clean
clean:
	rm -rf obj $(OBJS) NO_Au_efld NO_Au_cme NO_Au_test
