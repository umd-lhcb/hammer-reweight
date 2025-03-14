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

exe: PrintMCDecay ReweightRDX ReweightRDXDebug ReweightRDXDefault ReweightRDXDst10Sig ReweightRDXDstNoCorr ReweightRDXDstNoCorr10Sig ReweightRDXDstRun1 ReweightRDXRemoveRescale

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

ReweightRDX: ReweightRDX.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

ReweightRDXDefault: ReweightRDXDefault.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

ReweightRDXDebug: ReweightRDX.cpp
	$(COMPILER) $(CXXFLAGS) -DDEBUG_CLI -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

ReweightRDXDst10Sig: ReweightRDXDst10Sig.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

ReweightRDXDstNoCorr: ReweightRDXDstNoCorr.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

ReweightRDXDstNoCorr10Sig: ReweightRDXDstNoCorr10Sig.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

ReweightRDXDstRun1: ReweightRDXDstRun1.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

ReweightRDXRemoveRescale: ReweightRDXRemoveRescale.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

ff-params-RDX:
	./utils/gen_ham_params.py ./spec/rdx-run2.yml

ff-params-RDX-no-rescale:
	./utils/gen_ham_params_no_rescale.py ./spec/rdx-run2.yml


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
		--bins 10 \
		-XL "FF weights" \
		-l "HAMMER - theory" \
		--cuts "ham_ok & ff_calc_ok" \
		--debug
	@plotbr -n $</tree_BDst -o gen/rdx-run2-validation-B02Dst_q2_true.png \
		-b q2_true q2_true q2_true q2_true q2_true q2_true q2_true \
		-XD 3.2 11.0 --bins 20 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-l "Theory CLN" "HAM CLN" "HAM BGL (default)" "HAM BGL (lattice)" "HAM BGL +$$\sigma$$ (shift, lattice)" \
		--cuts "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" \
		--weights "wff_calc" "wff" "wff_bgl_n3" "wff_bgl" "wff_bgl_var_p_shift" \
		--colors cornflowerblue black darkgoldenrod crimson limegreen purple deeppink \
		--debug
	@echo "Generating B -> D validation plot..."
	@plotbr -n $</tree_BD -o gen/rdx-run2-validation-B2D0.png \
		-b "wff-wff_calc" \
		--bins 10 \
		-XL "FF weights" \
		-l "HAMMER - theory" \
		--cuts "ham_ok & ff_calc_ok" \
		--debug
	@plotbr -n $</tree_BD -o gen/rdx-run2-validation-B2D0_q2_true.png \
		-b q2_true q2_true q2_true q2_true q2_true q2_true q2_true \
		-XD 3.2 11.0 --bins 20 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-l "Theory CLN" "HAM CLN" "HAM BGL (\$$N=3$$)" "HAM BGL (\$$N=2$$)" "HAM BGL +\$$\\sigma$$ (\$$N=2$$)" "HAM BGL +$$\\sigma$$ (shift, \$$N=2$$)" \
		--cuts "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" \
		--weights "wff_calc" "wff" "wff_bgl_n3" "wff_bgl" "wff_bgl_var_p" "wff_bgl_var_p_shift" \
		--colors cornflowerblue black darkgoldenrod crimson limegreen purple deeppink \
		--debug
	@plotbr -n $</tree_BD -o gen/rdx-run2-validation-B2D0_q2_true_no_shift.png \
		-b q2_true q2_true q2_true q2_true q2_true q2_true q2_true \
		-XD 3.2 11.0 --bins 20 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-l "Theory CLN" "HAM CLN" "HAM BGL (\$$N=3$$)" "HAM BGL (\$$N=2$$)" "HAM BGL +\$$\\sigma$$ (\$$N=2$$)" \
		--cuts "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" \
		--weights "wff_calc" "wff" "wff_bgl_n3" "wff_bgl" "wff_bgl_var_p" "wff_bgl_var_p_shift" \
		--colors cornflowerblue black darkgoldenrod crimson limegreen purple deeppink \
		--debug
	@plotbr -n $</tree_BD -o gen/rdx-run2-validation-B2D0_q2_true.pdf \
		-b q2_true q2_true q2_true q2_true q2_true q2_true q2_true \
		-XD 3.2 11.0 --bins 20 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-l "Theory CLN" "HAM CLN" "HAM BGL (\$$N=3$$)" "HAM BGL (\$$N=2$$)" "HAM BGL +\$$\\sigma$$ (\$$N=2$$)" "HAM BGL +$$\\sigma$$ (shift, \$$N=2$$)" \
		--cuts "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" \
		--weights "wff_calc" "wff" "wff_bgl_n3" "wff_bgl" "wff_bgl_var_p" "wff_bgl_var_p_shift" \
		--colors cornflowerblue black darkgoldenrod crimson limegreen purple deeppink \
		--debug
	@plotbr -n $</tree_BD -o gen/rdx-run2-validation-B2D0_q2_true_no_shift.pdf \
		-b q2_true q2_true q2_true q2_true q2_true q2_true q2_true \
		-XD 3.2 11.0 --bins 20 \
		-XL "True \$$q^2$$ [GeV\$$^2$$]" \
		-l "Theory CLN" "HAM CLN" "HAM BGL (\$$N=3$$)" "HAM BGL (\$$N=2$$)" "HAM BGL +\$$\\sigma$$ (\$$N=2$$)" \
		--cuts "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" "ham_ok & ff_calc_ok" \
		--weights "wff_calc" "wff" "wff_bgl_n3" "wff_bgl" "wff_bgl_var_p" "wff_bgl_var_p_shift" \
		--colors cornflowerblue black darkgoldenrod crimson limegreen purple deeppink \
		--debug

rdx-run2-ntuples: \
	gen/rdx-run2-Bd2Dst0MuNu-sim09k-reweighted.root \
	gen/rdx-run2-Bd2DstMuNu-reweighted.root \
	gen/rdx-run2-Bd2DststTauNu-reweighted.root \
	gen/rdx-run2-Bd2DstTauNu-reweighted.root \
	gen/rdx-run2-Bd2D0DX_MuNu-reweighted.root \
	gen/rdx-run2-Bd2DststMuNu_D0_cocktail-sim09k-reweighted.root

rdx-run2-ntuples-latest: \
	gen/rdx-run2-Bd2DststMuNu_D0_cocktail-sim09k-reweighted.root


####################
# Generic patterns #
####################

.SECONDARY:

# Reweighters with HAMMER
%: %.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)

# Weight ntuples
gen/rdx-run1-%-reweighted.root: samples/rdx-run1-%.root ReweightRDX
	$(word 2, $^) $< $@ -r run1

gen/rdx-run2-%-reweighted.root: samples/rdx-run2-%.root ReweightRDX
	$(word 2, $^) $< $@ | tee gen/$(basename $(notdir $@)).log

# Validation ntuples
gen/rdx-run2-validation.root: ValidateRDX
	$< $@

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
