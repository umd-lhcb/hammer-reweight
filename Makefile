BINPATH	:=	bin
VPATH	:=	utils:src:$(BINPATH)

# Compiler settings
COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags) -Iinclude
LINKFLAGS	:=	$(shell root-config --libs)
ADDLINKFLAGS	:=	-lHammerTools -lHammerBase -lHammerCore -lFormFactors -lAmplitudes -lRates -lEG
VALLINKFLAGS	:=	-lff_calc


###########
# General #
###########

exe: PrintMCDecay #ReweightRDX ReweightRDXDebug

.PHONY: clean
clean:
	@rm -rf ./bin/*
	@rm -rf ./gen/*

install:
	pip install -r ./requirements.txt


#########
# Tools #
#########

PrintMCDecay: PrintMCDecay.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) -lEG

ValidateRDX: ValidateRDX.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS) $(VALLINKFLAGS)

ReweightRDXDebug: ReweightRDX.cpp
	$(COMPILER) $(CXXFLAGS) -DDEBUG_CLI -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)


#########
# Plots #
#########

.PHONY: sample-plots validation-plots rdx-run2-ntuples

sample-plots: \
	gen/rdx-run1-Bd2DstMuNu_q2_true.png \
	gen/rdx-run1-Bd2DstTauNu_q2_true.png \
	gen/rdx-run2-Bd2DstMuNu_q2_true.png \
	gen/rdx-run2-Bd2DstTauNu_q2_true.png \
	gen/rdx-run2-Bd2DststTauNu_q2_true.png \
	gen/rdx-run2-Bd2DststTauNu_q2_true_D_1.png \
	gen/rdx-run2-Bd2DststTauNu_q2_true_Dst_2.png \
	gen/rdx-run2-Bd2DststTauNu_q2_true_Dst_0.png \
	gen/rdx-run2-Bd2DststTauNu_q2_true_Dst_1.png

validation-plots: gen/rdx-run2-validation.root
	@echo "Generating B -> D* validation plot..."
	@plotbr -n $</tree_BDst -o gen/rdx-run2-validation-B02Dst.png \
		-b "wff-wff_calc" \
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
		--weights "None" "None" "wff" "wff_calc" \
		--debug
	@echo "Generating B -> D validation plot..."
	@plotbr -n $</tree_BD -o gen/rdx-run2-validation-B2D0.png \
		-b "wff-wff_calc" \
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
		--weights "None" "None" "wff" "wff_calc" \
		--debug

# rdx-run2-ntuples: \
#     gen/rdx-run2-Bd2DstMuNu-reweighted.root \
#     gen/rdx-run2-Bd2DststTauNu-reweighted.root \
#     gen/rdx-run2-Bd2DstTauNu-reweighted.root \
#     gen/rdx-run2-Bd2D0DX_MuNu-reweighted.root

rdx-run2-ntuples: \
    gen/rdx-run2-Bd2DstMuNu-reweighted.root


####################
# Generic patterns #
####################

.SECONDARY:

# Reweighters with HAMMER
%: %.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

# Weight ntuples
gen/rdx-run1-%-reweighted.root: samples/rdx-run1-%.root ReweightRDX
	./bin/ReweightRDX $< $@ -r run1

gen/rdx-run2-%-reweighted.root: samples/rdx-run2-%.root ReweightRDXDebug
	./bin/ReweightRDXDebug $< $@ -r run2 | tee gen/$(basename $(notdir $@)).log

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
		--weights "None" "wff" \
		--normalize --debug

# True q2 plots for D**
gen/rdx-run2-Bd2DststTauNu_q2_true.png: gen/rdx-run2-Bd2DststTauNu-reweighted.root
	@echo "Generating $@..."
	@plotbr -n $</TupleB0/DecayTree -o $@ \
		-b q2_true q2_true -XD 2.5 8.5 --bins 15 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-YL "Normalized" \
		-l Original Reweighted \
		--cuts "ham_ok" "ham_ok" \
		--weights "None" "wff" \
		--normalize --debug

gen/rdx-run2-Bd2DststTauNu_q2_true_D_1.png: gen/rdx-run2-Bd2DststTauNu-reweighted.root
	@echo "Generating $@..."
	@plotbr -n $</TupleB0/DecayTree -o $@ \
		-b q2_true q2_true -XD 2.5 8.5 --bins 15 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-YL "Normalized" \
		-l Original Reweighted \
		--cuts "ham_ok & abs(d_meson1_true_id) == 10413" "ham_ok & abs(d_meson1_true_id) == 10413" \
		--weights "None" "wff" \
		--normalize --debug

gen/rdx-run2-Bd2DststTauNu_q2_true_Dst_1.png: gen/rdx-run2-Bd2DststTauNu-reweighted.root
	@echo "Generating $@..."
	@plotbr -n $</TupleB0/DecayTree -o $@ \
		-b q2_true q2_true -XD 2.5 8.5 --bins 15 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-YL "Normalized" \
		-l Original Reweighted \
		--cuts "ham_ok & abs(d_meson1_true_id) == 20413" "ham_ok & abs(d_meson1_true_id) == 20413" \
		--weights "None" "wff" \
		--normalize --debug

gen/rdx-run2-Bd2DststTauNu_q2_true_Dst_0.png: gen/rdx-run2-Bd2DststTauNu-reweighted.root
	@echo "Generating $@..."
	@plotbr -n $</TupleB0/DecayTree -o $@ \
		-b q2_true q2_true -XD 2.5 8.5 --bins 15 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-YL "Normalized" \
		-l Original Reweighted \
		--cuts "ham_ok & (abs(d_meson1_true_id) == 10421 | abs(d_meson1_true_id) == 10411)" "ham_ok & (abs(d_meson1_true_id) == 10421 | abs(d_meson1_true_id) == 10411)" \
		--weights "None" "wff" \
		--normalize --debug

gen/rdx-run2-Bd2DststTauNu_q2_true_Dst_2.png: gen/rdx-run2-Bd2DststTauNu-reweighted.root
	@echo "Generating $@..."
	@plotbr -n $</TupleB0/DecayTree -o $@ \
		-b q2_true q2_true -XD 2.5 8.5 --bins 15 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-YL "Normalized" \
		-l Original Reweighted \
		--cuts "ham_ok & abs(d_meson1_true_id) == 415" "ham_ok & abs(d_meson1_true_id) == 415" \
		--weights "None" "wff" \
		--normalize --debug
