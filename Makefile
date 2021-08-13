BINPATH	:=	bin
VPATH	:=	gen:utils:src:$(BINPATH)

# Compiler settings
COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags)
LINKFLAGS	:=	$(shell root-config --libs)
ADDLINKFLAGS	:=	-lHammerTools -lHammerBase -lHammerCore -lFormFactors -lAmplitudes -lRates
VALLINKFLAGS	:=	-lff_dstaunu


###########
# General #
###########

tools: PrintMCDecay ReweightRDX

.PHONY: clean
clean:
	@rm -rf ./bin/*
	@rm -rf ./gen/*


#########
# Tools #
#########

PrintMCDecay: PrintMCDecay.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) -lEG


###############
# Sample plot #
###############

sample-plots: \
	rdx-run1-Bd2DstMuNu_q2_true.png


####################
# Generic patterns #
####################

# Reweighters with HAMMER
%: %.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

# Weight ntuples
rdx-run1-%-reweighted.root: samples/rdx-run1-%.root ReweightRDX
	./bin/ReweightRDX $< gen/$@ TupleB0/DecayTree

# True q2 plots
%_q2_true.png: gen/%-reweighted.root Makefile
	@echo "Generating $@..."
	@plotbr -n $</TupleB0/DecayTree -o ./gen/$@ \
		-b q2_true q2_true -XD 0 12 --bins 15 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-YL "Normalized" \
		-l Original Reweighted \
		--cuts "true" "flag_ham_ok" \
		--weights "None" "w_ff" \
		--normalize --debug
