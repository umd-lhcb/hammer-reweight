BINPATH	:=	bin
VPATH	:=	utils:src:$(BINPATH)

# Compiler settings
COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags) -Iinclude
LINKFLAGS	:=	$(shell root-config --libs)
ADDLINKFLAGS	:=	-lHammerTools -lHammerBase -lHammerCore -lFormFactors -lAmplitudes -lRates
VALLINKFLAGS	:=	-lff_calc


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

ValidateRDX: ValidateRDX.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS) $(VALLINKFLAGS)


#########
# Plots #
#########

sample-plots: \
	gen/rdx-run1-Bd2DstMuNu_q2_true.png \
	gen/rdx-run1-Bd2DstTauNu_q2_true.png \
	gen/rdx-run2-Bd2DstMuNu_q2_true.png \
	gen/rdx-run2-Bd2DstTauNu_q2_true.png \
	gen/rdx-run2-Bd2DststTauNu_q2_true.png

validation-plots: gen/rdx-run2-validation.root
	@echo "Generating B -> D* validation plot..."
	@plotbr -n $</tree_BDst -o gen/rdx-run2-validation-B02Dst.png \
		-b "w_ff-w_ff_calc" \
		--bins 5 \
		-XL "FF weights" \
		-l "HAMMER - theory" \
		--cuts "ham_ok & ff_calc_ok" \
		--debug
	@plotbr -n $</tree_BDst -o gen/rdx-run2-validation-B02Dst_q2_true.png \
		-b q2_true q2_true q2_true q2_true -XD 3.2 11.0 --bins 20 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-l Generated "HAMMER OK" "HAMMER wt" "Theory wt" \
		--cuts "None" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" \
		--weights "None" "None" "w_ff" "w_ff_calc" \
		--debug
	@echo "Generating B -> D validation plot..."
	@plotbr -n $</tree_BD -o gen/rdx-run2-validation-B2D0.png \
		-b "w_ff-w_ff_calc" \
		--bins 5 \
		-XL "FF weights" \
		-l "HAMMER - theory" \
		--cuts "ham_ok & ff_calc_ok" \
		--debug
	@plotbr -n $</tree_BD -o gen/rdx-run2-validation-B2D0_q2_true.png \
		-b q2_true q2_true q2_true q2_true -XD 3.2 11.0 --bins 20 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-l Generated "HAMMER OK" "HAMMER wt" "Theory wt" \
		--cuts "None" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" \
		--weights "None" "None" "w_ff" "w_ff_calc" \
		--debug

rdx-run2-ntuples: \
	gen/rdx-run2-Bd2DstMuNu-reweighted.root \
	gen/rdx-run2-Bd2DststTauNu-reweighted.root \
	gen/rdx-run2-Bd2DstTauNu-reweighted.root \


####################
# Generic patterns #
####################

.SECONDARY:

# Reweighters with HAMMER
%: %.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

# Weight ntuples
gen/rdx-run1-%-reweighted.root: samples/rdx-run1-%.root ReweightRDX
	./bin/ReweightRDX $< $@ TupleB0/DecayTree run1

gen/rdx-run2-%-reweighted.root: samples/rdx-run2-%.root ReweightRDX
	./bin/ReweightRDX $< $@ TupleB0/DecayTree run2
	./bin/ReweightRDX $< $@ TupleBminus/DecayTree run2

# Validation ntuples
gen/rdx-run2-validation.root: ValidateRDX
	./bin/ValidateRDX $@

# True q2 plots
gen/%_q2_true.png: gen/%-reweighted.root
	@echo "Generating $@..."
	@plotbr -n $</TupleB0/DecayTree -o $@ \
		-b q2_true q2_true -XD 0 12 --bins 15 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-YL "Normalized" \
		-l Original Reweighted \
		--cuts "ham_ok" "ham_ok" \
		--weights "None" "w_ff" \
		--normalize --debug
