// Author: Yipeng Sun
// Last Change: Thu Sep 16, 2021 at 06:51 PM +0200

#include <algorithm>
#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <math.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <Hammer/Hammer.hh>
#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>
#include <Hammer/Process.hh>

#include "utils_general.h"
#include "utils_ham.h"

using namespace std;

///////////////////
// Configurables //
///////////////////

//#define SILENT
//#define FORCE_MOMENTUM_CONSERVATION_LEPTONIC
//#define FORCE_MOMENTUM_CONSERVATION_HADRONIC
#define RADIATIVE_CORRECTION

typedef map<vector<Int_t>, unsigned long> DecayFreq;

// clang-format off
auto B_MESON = map<TString, TString>{
  {"TupleBminus/DecayTree", "b"},
  {"TupleB0/DecayTree", "b0"}
};

const auto LEGAL_B_MESON_IDS = vector<Int_t>{511, 521};

void set_input_ff(Hammer::Hammer& ham, TString run) {
  if (run == "run1") {
    ham.setFFInputScheme({
      {"BD", "ISGW2"},  // 12573010
      {"BD*", "ISGW2"}  // 11574010
      // for run 1, B -> D*MuNu is modelled w/ HQET, which is not implemented in HAMMER
    });

  } else if (run == "run2") {
    ham.setFFInputScheme({
      {"BD", "CLN_1"},
      {"BD*", "CLN_2"},
      {"BD**0*", "ISGW2"},
      {"BD**1", "ISGW2"},
      {"BD**1*", "ISGW2"},
      {"BD**2*", "ISGW2"}
    });

    // 12573001, 12573012
    // HQET2(hqetrho2, hqetv1_1, indelta): 1.131 1.035 0.38
    ham.setOptions("BtoDCLN_1: {RhoSq: 1.131, Delta: 0.38, G1: 1.035}");  // HQET3
    // 11574011, 11574021
    // HQET2(hqetrho2, hqetha1_1, hqetr1_1, hqetr2_1, hqetr0_1): 1.122 0.908 1.270 0.852 1.15
    ham.setOptions("BtoD*CLN_2: {RhoSq: 1.122, F1: 0.908, R1: 1.270, R2: 0.852, R0: 1.15}");  // HQET2
    // 11874440
    // ISGW2, which has no configurable parameter
  }
}

void set_output_ff(Hammer::Hammer& ham) {
  ham.addFFScheme("OutputFF", {
      {"BD", "BGL"},
      {"BD*", "BGL"},
      {"BD**0*", "BLR"},
      {"BD**1", "BLR"},
      {"BD**1*", "BLR"},
      {"BD**2*", "BLR"}
  });
}
// clang-format on

void set_decays(Hammer::Hammer& ham) {
  ham.includeDecay("BDTauNu");
  ham.includeDecay("BDMuNu");

  ham.includeDecay("BD*TauNu");
  ham.includeDecay("BD*MuNu");

  ham.includeDecay("BD**0*TauNu");
  ham.includeDecay("BD**1TauNu");
  ham.includeDecay("BD**1*TauNu");
  ham.includeDecay("BD**2*TauNu");

  ham.includeDecay("BD**0*MuNu");
  ham.includeDecay("BD**1MuNu");
  ham.includeDecay("BD**1*MuNu");
  ham.includeDecay("BD**2*MuNu");
}

const Double_t SOFT_PHOTON_THRESH = 0.1;

/////////////////
// Corrections //
/////////////////

string photon_correction(int ref_mom_id, vector<Int_t> photon_mom_id,
                         vector<Hammer::Particle> photon_p,
                         Hammer::Process& proc, Hammer::ParticleIndices& idx,
                         TDatabasePDG* db) {
  stringstream buffer;
  for (auto i = 0; i < photon_p.size(); i++) {
    if (photon_mom_id[i] == ref_mom_id) {
      buffer << "  Adding photon: " << print_p(photon_p[i].p())
             << " to particle " << get_particle_name(ref_mom_id, db) << endl;
      idx.push_back(proc.addParticle(photon_p[i]));
    }
  }

  return buffer.str();
}

Bool_t is_soft_photon(Double_t pe) {
  if (TMath::Abs(pe) > SOFT_PHOTON_THRESH) return false;
  return true;
}

/////////////////
// Reweighting //
/////////////////

typedef pair<unsigned long, unsigned long> RwRate;

RwRate reweight(TFile* input_ntp, TFile* output_ntp, TString tree,
                Hammer::Hammer& ham) {
  auto db = new TDatabasePDG();

  if (B_MESON.find(tree) == B_MESON.end()) return RwRate{0, 0};
  TString b_meson = B_MESON[tree];

  // Input /////////////////////////////////////////////////////////////////////
  TTreeReader reader(tree, input_ntp);

  // B meson truth info
  TTreeReaderValue<Double_t> b_px(reader, b_meson + "_TRUEP_X");
  TTreeReaderValue<Double_t> b_py(reader, b_meson + "_TRUEP_Y");
  TTreeReaderValue<Double_t> b_pz(reader, b_meson + "_TRUEP_Z");
  TTreeReaderValue<Double_t> b_pe(reader, b_meson + "_TRUEP_E");
  TTreeReaderValue<Int_t>    b_id(reader, b_meson + "_TRUEID");
  TTreeReaderValue<Double_t> q2(reader, b_meson + "_True_Q2");
  TTreeReaderValue<Bool_t>   is_tau(reader, b_meson + "_True_IsTauDecay");

  // Muon
  TTreeReaderValue<Double_t> mu_px(reader, b_meson + "_TrueMu_PX");
  TTreeReaderValue<Double_t> mu_py(reader, b_meson + "_TrueMu_PY");
  TTreeReaderValue<Double_t> mu_pz(reader, b_meson + "_TrueMu_PZ");
  TTreeReaderValue<Double_t> mu_pe(reader, b_meson + "_TrueMu_PE");
  TTreeReaderValue<Int_t>    mu_id(reader, "mu_TRUEID");

  // Tau
  TTreeReaderValue<Double_t> tau_px(reader, b_meson + "_TrueTau_PX");
  TTreeReaderValue<Double_t> tau_py(reader, b_meson + "_TrueTau_PY");
  TTreeReaderValue<Double_t> tau_pz(reader, b_meson + "_TrueTau_PZ");
  TTreeReaderValue<Double_t> tau_pe(reader, b_meson + "_TrueTau_PE");

  // The neutrino that is directly associated with the B -> D weak decay
  // Could be a Mu or Tau anti-neutrino.
  TTreeReaderValue<Double_t> anu_px(reader, b_meson + "_TrueNeutrino_PX");
  TTreeReaderValue<Double_t> anu_py(reader, b_meson + "_TrueNeutrino_PY");
  TTreeReaderValue<Double_t> anu_pz(reader, b_meson + "_TrueNeutrino_PZ");
  TTreeReaderValue<Double_t> anu_pe(reader, b_meson + "_TrueNeutrino_PE");

  // Tau neutrino and Mu anti-neutrino associated with the secondary tau decay
  TTreeReaderValue<Double_t> nu_tau_px(reader, b_meson + "_TrueTauNuTau_PX");
  TTreeReaderValue<Double_t> nu_tau_py(reader, b_meson + "_TrueTauNuTau_PY");
  TTreeReaderValue<Double_t> nu_tau_pz(reader, b_meson + "_TrueTauNuTau_PZ");
  TTreeReaderValue<Double_t> nu_tau_pe(reader, b_meson + "_TrueTauNuTau_PE");

  TTreeReaderValue<Double_t> anu_mu_px(reader, b_meson + "_TrueTauNuMu_PX");
  TTreeReaderValue<Double_t> anu_mu_py(reader, b_meson + "_TrueTauNuMu_PY");
  TTreeReaderValue<Double_t> anu_mu_pz(reader, b_meson + "_TrueTauNuMu_PZ");
  TTreeReaderValue<Double_t> anu_mu_pe(reader, b_meson + "_TrueTauNuMu_PE");

  // First D meson and its daughters
  TTreeReaderValue<Double_t> d_idx0_px(reader, b_meson + "_TrueHadron_D0_PX");
  TTreeReaderValue<Double_t> d_idx0_py(reader, b_meson + "_TrueHadron_D0_PY");
  TTreeReaderValue<Double_t> d_idx0_pz(reader, b_meson + "_TrueHadron_D0_PZ");
  TTreeReaderValue<Double_t> d_idx0_pe(reader, b_meson + "_TrueHadron_D0_PE");
  TTreeReaderValue<Int_t>    d_idx0_id(reader, b_meson + "_TrueHadron_D0_ID");

  TTreeReaderValue<Double_t> d_idx0_gd0_px(reader,
                                           b_meson + "_TrueHadron_D0_GD0_PX");
  TTreeReaderValue<Double_t> d_idx0_gd0_py(reader,
                                           b_meson + "_TrueHadron_D0_GD0_PY");
  TTreeReaderValue<Double_t> d_idx0_gd0_pz(reader,
                                           b_meson + "_TrueHadron_D0_GD0_PZ");
  TTreeReaderValue<Double_t> d_idx0_gd0_pe(reader,
                                           b_meson + "_TrueHadron_D0_GD0_PE");
  TTreeReaderValue<Int_t>    d_idx0_gd0_id(reader,
                                        b_meson + "_TrueHadron_D0_GD0_ID");

  TTreeReaderValue<Double_t> d_idx0_gd1_px(reader,
                                           b_meson + "_TrueHadron_D0_GD1_PX");
  TTreeReaderValue<Double_t> d_idx0_gd1_py(reader,
                                           b_meson + "_TrueHadron_D0_GD1_PY");
  TTreeReaderValue<Double_t> d_idx0_gd1_pz(reader,
                                           b_meson + "_TrueHadron_D0_GD1_PZ");
  TTreeReaderValue<Double_t> d_idx0_gd1_pe(reader,
                                           b_meson + "_TrueHadron_D0_GD1_PE");
  TTreeReaderValue<Int_t>    d_idx0_gd1_id(reader,
                                        b_meson + "_TrueHadron_D0_GD1_ID");

  TTreeReaderValue<Double_t> d_idx0_gd2_px(reader,
                                           b_meson + "_TrueHadron_D0_GD2_PX");
  TTreeReaderValue<Double_t> d_idx0_gd2_py(reader,
                                           b_meson + "_TrueHadron_D0_GD2_PY");
  TTreeReaderValue<Double_t> d_idx0_gd2_pz(reader,
                                           b_meson + "_TrueHadron_D0_GD2_PZ");
  TTreeReaderValue<Double_t> d_idx0_gd2_pe(reader,
                                           b_meson + "_TrueHadron_D0_GD2_PE");
  TTreeReaderValue<Int_t>    d_idx0_gd2_id(reader,
                                        b_meson + "_TrueHadron_D0_GD2_ID");

  // Second D meson and its daughters
  TTreeReaderValue<Double_t> d_idx1_px(reader, b_meson + "_TrueHadron_D1_PX");
  TTreeReaderValue<Double_t> d_idx1_py(reader, b_meson + "_TrueHadron_D1_PY");
  TTreeReaderValue<Double_t> d_idx1_pz(reader, b_meson + "_TrueHadron_D1_PZ");
  TTreeReaderValue<Double_t> d_idx1_pe(reader, b_meson + "_TrueHadron_D1_PE");
  TTreeReaderValue<Int_t>    d_idx1_id(reader, b_meson + "_TrueHadron_D1_ID");

  TTreeReaderValue<Double_t> d_idx1_gd0_px(reader,
                                           b_meson + "_TrueHadron_D1_GD0_PX");
  TTreeReaderValue<Double_t> d_idx1_gd0_py(reader,
                                           b_meson + "_TrueHadron_D1_GD0_PY");
  TTreeReaderValue<Double_t> d_idx1_gd0_pz(reader,
                                           b_meson + "_TrueHadron_D1_GD0_PZ");
  TTreeReaderValue<Double_t> d_idx1_gd0_pe(reader,
                                           b_meson + "_TrueHadron_D1_GD0_PE");
  TTreeReaderValue<Int_t>    d_idx1_gd0_id(reader,
                                        b_meson + "_TrueHadron_D1_GD0_ID");

  TTreeReaderValue<Double_t> d_idx1_gd1_px(reader,
                                           b_meson + "_TrueHadron_D1_GD1_PX");
  TTreeReaderValue<Double_t> d_idx1_gd1_py(reader,
                                           b_meson + "_TrueHadron_D1_GD1_PY");
  TTreeReaderValue<Double_t> d_idx1_gd1_pz(reader,
                                           b_meson + "_TrueHadron_D1_GD1_PZ");
  TTreeReaderValue<Double_t> d_idx1_gd1_pe(reader,
                                           b_meson + "_TrueHadron_D1_GD1_PE");
  TTreeReaderValue<Int_t>    d_idx1_gd1_id(reader,
                                        b_meson + "_TrueHadron_D1_GD1_ID");

  TTreeReaderValue<Double_t> d_idx1_gd2_px(reader,
                                           b_meson + "_TrueHadron_D1_GD2_PX");
  TTreeReaderValue<Double_t> d_idx1_gd2_py(reader,
                                           b_meson + "_TrueHadron_D1_GD2_PY");
  TTreeReaderValue<Double_t> d_idx1_gd2_pz(reader,
                                           b_meson + "_TrueHadron_D1_GD2_PZ");
  TTreeReaderValue<Double_t> d_idx1_gd2_pe(reader,
                                           b_meson + "_TrueHadron_D1_GD2_PE");
  TTreeReaderValue<Int_t>    d_idx1_gd2_id(reader,
                                        b_meson + "_TrueHadron_D1_GD2_ID");

  // Third D meson and its daughters
  TTreeReaderValue<Double_t> d_idx2_px(reader, b_meson + "_TrueHadron_D2_PX");
  TTreeReaderValue<Double_t> d_idx2_py(reader, b_meson + "_TrueHadron_D2_PY");
  TTreeReaderValue<Double_t> d_idx2_pz(reader, b_meson + "_TrueHadron_D2_PZ");
  TTreeReaderValue<Double_t> d_idx2_pe(reader, b_meson + "_TrueHadron_D2_PE");
  TTreeReaderValue<Int_t>    d_idx2_id(reader, b_meson + "_TrueHadron_D2_ID");

  TTreeReaderValue<Double_t> d_idx2_gd0_px(reader,
                                           b_meson + "_TrueHadron_D2_GD0_PX");
  TTreeReaderValue<Double_t> d_idx2_gd0_py(reader,
                                           b_meson + "_TrueHadron_D2_GD0_PY");
  TTreeReaderValue<Double_t> d_idx2_gd0_pz(reader,
                                           b_meson + "_TrueHadron_D2_GD0_PZ");
  TTreeReaderValue<Double_t> d_idx2_gd0_pe(reader,
                                           b_meson + "_TrueHadron_D2_GD0_PE");
  TTreeReaderValue<Int_t>    d_idx2_gd0_id(reader,
                                        b_meson + "_TrueHadron_D2_GD0_ID");

  TTreeReaderValue<Double_t> d_idx2_gd1_px(reader,
                                           b_meson + "_TrueHadron_D2_GD1_PX");
  TTreeReaderValue<Double_t> d_idx2_gd1_py(reader,
                                           b_meson + "_TrueHadron_D2_GD1_PY");
  TTreeReaderValue<Double_t> d_idx2_gd1_pz(reader,
                                           b_meson + "_TrueHadron_D2_GD1_PZ");
  TTreeReaderValue<Double_t> d_idx2_gd1_pe(reader,
                                           b_meson + "_TrueHadron_D2_GD1_PE");
  TTreeReaderValue<Int_t>    d_idx2_gd1_id(reader,
                                        b_meson + "_TrueHadron_D2_GD1_ID");

  TTreeReaderValue<Double_t> d_idx2_gd2_px(reader,
                                           b_meson + "_TrueHadron_D2_GD2_PX");
  TTreeReaderValue<Double_t> d_idx2_gd2_py(reader,
                                           b_meson + "_TrueHadron_D2_GD2_PY");
  TTreeReaderValue<Double_t> d_idx2_gd2_pz(reader,
                                           b_meson + "_TrueHadron_D2_GD2_PZ");
  TTreeReaderValue<Double_t> d_idx2_gd2_pe(reader,
                                           b_meson + "_TrueHadron_D2_GD2_PE");
  TTreeReaderValue<Int_t>    d_idx2_gd2_id(reader,
                                        b_meson + "_TrueHadron_D2_GD2_ID");

#ifdef RADIATIVE_CORRECTION
  TTreeReaderValue<Int_t> photon_arr_size(
      reader, b_meson + "_MCTrue_gamma_ArrayLength");
  TTreeReaderArray<Float_t> photon_arr_pe(reader, b_meson + "_MCTrue_gamma_E");
  TTreeReaderArray<Float_t> photon_arr_px(reader, b_meson + "_MCTrue_gamma_PX");
  TTreeReaderArray<Float_t> photon_arr_py(reader, b_meson + "_MCTrue_gamma_PY");
  TTreeReaderArray<Float_t> photon_arr_pz(reader, b_meson + "_MCTrue_gamma_PY");
  TTreeReaderArray<Float_t> photon_arr_mom_id(
      reader, b_meson + "_MCTrue_gamma_mother_ID");
#endif

  // Event ID
  TTreeReaderValue<ULong64_t> eventNumber(reader, "eventNumber");
  TTreeReaderValue<UInt_t>    runNumber(reader, "runNumber");

  // Output ////////////////////////////////////////////////////////////////////
  // Recreate the same folder structure in output
  auto tree_dir  = dirname(string(tree));
  auto tree_name = basename(string(tree));

  output_ntp->cd();
  if (tree_dir != TString("")) {
    output_ntp->mkdir(tree_dir);
    output_ntp->cd(tree_dir);
  }

  // Output branches
  auto output_tree = new TTree(tree_name, tree_name);

  ULong64_t eventNumber_out;
  output_tree->Branch("eventNumber", &eventNumber_out);
  UInt_t runNumber_out;
  output_tree->Branch("runNumber", &runNumber_out);

  Double_t w_ff_out;
  output_tree->Branch("w_ff", &w_ff_out);

  Bool_t ham_ok;
  output_tree->Branch("ham_ok", &ham_ok);
  Bool_t ham_tm_ok;
  output_tree->Branch("ham_tm_ok", &ham_tm_ok);

  Double_t q2_true_out;
  output_tree->Branch("q2_true", &q2_true_out);

  unsigned long num_of_evt        = 0l;
  unsigned long num_of_evt_ham_ok = 0l;
  while (reader.Next()) {
    ham_ok      = false;
    ham_tm_ok   = false;
    w_ff_out    = 1.;
    q2_true_out = *q2 / 1000 / 1000;

    double q2_min = 100 * 100;
    if (*is_tau) q2_min = 1700 * 1700;

    // Check if we have a legal B meson and q2 is large enough to produce a Mu
    if (find_in(LEGAL_B_MESON_IDS, TMath::Abs(*b_id)) && *q2 > q2_min &&
        TMath::Abs(*mu_id) == 13) {
      // Check if we have a legal D meson
      // clang-format off
      auto D_cands =
          PartIdMap{{"D0", *d_idx0_id}, {"D1", *d_idx1_id}, {"D2", *d_idx2_id}};
      auto D_mom = PartMomMap{
          {"D0_PX", *d_idx0_px},
          {"D0_PY", *d_idx0_py},
          {"D0_PZ", *d_idx0_pz},
          {"D0_PE", *d_idx0_pe},

          {"D1_PX", *d_idx1_px},
          {"D1_PY", *d_idx1_py},
          {"D1_PZ", *d_idx1_pz},
          {"D1_PE", *d_idx1_pe},

          {"D2_PX", *d_idx2_px},
          {"D2_PY", *d_idx2_py},
          {"D2_PZ", *d_idx2_pz},
          {"D2_PE", *d_idx2_pe}
      };
      // clang-format on
      auto [D_ok, D_lbl] = is_D_meson(D_cands);

      if (D_ok) {
        ham_tm_ok = true;  // now the naive truth-matching is considered OK

        // clang-format off
        auto D_daughter_id = PartIdMap{
          {"D0_GD0", *d_idx0_gd0_id},
          {"D0_GD1", *d_idx0_gd1_id},
          {"D0_GD2", *d_idx0_gd2_id},

          {"D1_GD0", *d_idx1_gd0_id},
          {"D1_GD1", *d_idx1_gd1_id},
          {"D1_GD2", *d_idx1_gd2_id},

          {"D2_GD0", *d_idx2_gd0_id},
          {"D2_GD1", *d_idx2_gd1_id},
          {"D2_GD2", *d_idx2_gd2_id}
        };

        auto D_daughter_mom = PartMomMap{
          {"D0_GD0_PX", *d_idx0_gd0_px},
          {"D0_GD0_PY", *d_idx0_gd0_py},
          {"D0_GD0_PZ", *d_idx0_gd0_pz},
          {"D0_GD0_PE", *d_idx0_gd0_pe},

          {"D0_GD1_PX", *d_idx0_gd1_px},
          {"D0_GD1_PY", *d_idx0_gd1_py},
          {"D0_GD1_PZ", *d_idx0_gd1_pz},
          {"D0_GD1_PE", *d_idx0_gd1_pe},

          {"D0_GD2_PX", *d_idx0_gd2_px},
          {"D0_GD2_PY", *d_idx0_gd2_py},
          {"D0_GD2_PZ", *d_idx0_gd2_pz},
          {"D0_GD2_PE", *d_idx0_gd2_pe},

          {"D1_GD0_PX", *d_idx1_gd0_px},
          {"D1_GD0_PY", *d_idx1_gd0_py},
          {"D1_GD0_PZ", *d_idx1_gd0_pz},
          {"D1_GD0_PE", *d_idx1_gd0_pe},

          {"D1_GD1_PX", *d_idx1_gd1_px},
          {"D1_GD1_PY", *d_idx1_gd1_py},
          {"D1_GD1_PZ", *d_idx1_gd1_pz},
          {"D1_GD1_PE", *d_idx1_gd1_pe},

          {"D1_GD2_PX", *d_idx1_gd2_px},
          {"D1_GD2_PY", *d_idx1_gd2_py},
          {"D1_GD2_PZ", *d_idx1_gd2_pz},
          {"D1_GD2_PE", *d_idx1_gd2_pe},

          {"D2_GD0_PX", *d_idx2_gd0_px},
          {"D2_GD0_PY", *d_idx2_gd0_py},
          {"D2_GD0_PZ", *d_idx2_gd0_pz},
          {"D2_GD0_PE", *d_idx2_gd0_pe},

          {"D2_GD1_PX", *d_idx2_gd1_px},
          {"D2_GD1_PY", *d_idx2_gd1_py},
          {"D2_GD1_PZ", *d_idx2_gd1_pz},
          {"D2_GD1_PE", *d_idx2_gd1_pe},

          {"D2_GD2_PX", *d_idx2_gd2_px},
          {"D2_GD2_PY", *d_idx2_gd2_py},
          {"D2_GD2_PZ", *d_idx2_gd2_pz},
          {"D2_GD2_PE", *d_idx2_gd2_pe}
        };
        // clang-format on

        // If so, locate the right D meson and loop over all its daughters and
        // add them all
        vector<Hammer::Particle> part_D_daughters{};

        for (auto suffix : vector<TString>{"_GD0", "_GD1", "_GD2"}) {
          auto part_name = D_lbl + suffix;
          auto part_id   = D_daughter_id[part_name];
          if (part_id != 0) {
            part_D_daughters.push_back(
                particle(D_daughter_mom[part_name + "_PE"],
                         D_daughter_mom[part_name + "_PX"],
                         D_daughter_mom[part_name + "_PY"],
                         D_daughter_mom[part_name + "_PZ"], part_id));
          }
        }

        // HAMMER process
        Hammer::Process proc;

        auto d_id       = D_cands[D_lbl];
        auto b_id_fixed = B_id_fix(*b_id, d_id);
        auto part_B     = particle(*b_pe, *b_px, *b_py, *b_pz, b_id_fixed);
        auto part_D =
            particle(D_mom[D_lbl + "_PE"], D_mom[D_lbl + "_PX"],
                     D_mom[D_lbl + "_PY"], D_mom[D_lbl + "_PZ"], d_id);

        Bool_t is_Dst = TMath::Abs(d_id) == 413 || TMath::Abs(d_id) == 423;

        Hammer::Particle part_L, part_NuL, part_TauNuTau, part_TauNuMu, part_Mu;
        Int_t            part_L_id = 0;
        if (*is_tau) {
          part_L_id = Tau_id(*mu_id);
          part_L    = particle(*tau_pe, *tau_px, *tau_py, *tau_pz, part_L_id);
        } else {
          part_L_id = Mu_id(*mu_id);
          part_L    = particle(*mu_pe, *mu_px, *mu_py, *mu_pz, part_L_id);
        }

        part_Mu  = particle(*mu_pe, *mu_px, *mu_py, *mu_pz, Mu_id(*mu_id));
        part_NuL = particle(*anu_pe, *anu_px, *anu_py, *anu_pz,
                            Nu_id(*mu_id, *is_tau));

        // FIXME: From Julian: The nu_mu and nu_tau kinematics are exchanged
        //        below, to correct for the bug in the LHCb EvtGen models
        // part_TauNuMu  = particle(*nu_tau_pe, *nu_tau_px, *nu_tau_py,
        // *nu_tau_pz, Tau_NuMu_id(*mu_id));
        // part_TauNuTau = particle(*anu_mu_pe, *anu_mu_px, *anu_mu_py,
        // *anu_mu_pz, Tau_NuTau_id(*mu_id));
        part_TauNuTau = particle(*nu_tau_pe, *nu_tau_px, *nu_tau_py, *nu_tau_pz,
                                 Tau_NuTau_id(*mu_id));
        part_TauNuMu  = particle(*anu_mu_pe, *anu_mu_px, *anu_mu_py, *anu_mu_pz,
                                Tau_NuMu_id(*mu_id));

        auto part_B_idx = proc.addParticle(part_B);
        auto part_D_idx = proc.addParticle(part_D);

#ifdef FORCE_MOMENTUM_CONSERVATION_LEPTONIC
        part_L.setMomentum(part_B.p() - part_D.p() - part_NuL.p());

        if (*is_tau)
          part_Mu.setMomentum(part_L.p() - part_TauNuMu.p() -
                              part_TauNuTau.p());
#endif

#ifdef FORCE_MOMENTUM_CONSERVATION_HADRONIC
        Hammer::FourMomentum known_mom{0, 0, 0, 0};
        for (auto idx = 0; idx < part_D_daughters.size() - 1; idx++) {
          known_mom += part_D_daughters[idx].p();
        }

        auto bal_mom = part_D.p() - known_mom;
        part_D_daughters[part_D_daughters.size() - 1].setMomentum(bal_mom);
#endif

#ifdef RADIATIVE_CORRECTION
        vector<Int_t>            vec_photon_mom_id{};
        vector<Hammer::Particle> vec_photon_p{};

        for (auto i = 0; i < *photon_arr_size; i++) {
          Float_t g_pe = photon_arr_pe[i];
          Float_t g_px = photon_arr_px[i];
          Float_t g_py = photon_arr_py[i];
          Float_t g_pz = photon_arr_pz[i];

          if (!is_soft_photon(g_pe)) {
            vec_photon_mom_id.push_back(
                static_cast<Int_t>(photon_arr_mom_id[i]));
            vec_photon_p.emplace_back(particle(g_pe, g_px, g_py, g_pz, 22));
          }
        }

        string buf_rad_corr;
#endif

        // Always add the primary leptons
        auto part_L_idx   = proc.addParticle(part_L);
        auto part_NuL_idx = proc.addParticle(part_NuL);

        // Particle indices, also decay tree definition
        Hammer::ParticleIndices part_B_daughters_idx{part_D_idx, part_L_idx,
                                                     part_NuL_idx};
#ifdef RADIATIVE_CORRECTION
        buf_rad_corr +=
            photon_correction(b_id_fixed, vec_photon_mom_id, vec_photon_p, proc,
                              part_B_daughters_idx, db);
#endif
        proc.addVertex(part_B_idx, part_B_daughters_idx);

        // Only add the D daughters if the D meson is a D* (0/+)
        Hammer::ParticleIndices part_D_daughters_idx{};
        if (is_Dst) {
          for (const auto part : part_D_daughters) {
            if (part.pdgId() == 22) {
              // Don't add soft photons
              if (!is_soft_photon(part.p().E())) {
                part_D_daughters_idx.push_back(proc.addParticle(part));
              }
            } else
              part_D_daughters_idx.push_back(proc.addParticle(part));
          }
#ifdef RADIATIVE_CORRECTION
          buf_rad_corr +=
              photon_correction(d_id, vec_photon_mom_id, vec_photon_p, proc,
                                part_D_daughters_idx, db);
#endif
          proc.addVertex(part_D_idx, part_D_daughters_idx);
        }

        if (*is_tau) {
          auto part_Mu_idx       = proc.addParticle(part_Mu);
          auto part_TauNuMu_idx  = proc.addParticle(part_TauNuMu);
          auto part_TauNuTau_idx = proc.addParticle(part_TauNuTau);

          Hammer::ParticleIndices part_L_daughters_idx{
              part_Mu_idx, part_TauNuMu_idx, part_TauNuTau_idx};
#ifdef RADIATIVE_CORRECTION
          buf_rad_corr +=
              photon_correction(part_L_id, vec_photon_mom_id, vec_photon_p,
                                proc, part_L_daughters_idx, db);
#endif
          proc.addVertex(part_L_idx, part_L_daughters_idx);
        }

        // Make sure invariant mass is non-negative
        vector<Hammer::Particle> part_vec{part_B,       part_D,  part_L,
                                          part_NuL,     part_Mu, part_TauNuMu,
                                          part_TauNuTau};
        vector<Bool_t>           part_m2_ok{};
        for (const auto p : part_vec) {
          // if (p.p().mass2() >= 0)  // This sometimes gives us negative
          // numbers for neutrinos due to float-point arithmetics
          if (p.p().mass() >= 0)  // HAMMER automatically apply fabs for very
                                  // small negative floats
            part_m2_ok.push_back(true);
          else
            part_m2_ok.push_back(false);
        }
        Bool_t all_parts_ok = find(part_m2_ok.begin(), part_m2_ok.end(),
                                   false) == part_m2_ok.end();

        Bool_t is_bad_cand = false;

        if (!all_parts_ok) {
          cout << "WARN: Bad kinematics for candidate: " << num_of_evt << endl;
          is_bad_cand = true;
        }

        if (!is_bad_cand) {
          ham.initEvent();

          Int_t proc_id = 0;
          try {
            proc_id = ham.addProcess(proc);
          } catch (const exception& e) {
            cout << "WARN: HAMMER doesn't add process properly: " << num_of_evt
                 << endl;
            cout << e.what() << endl;
            is_bad_cand = true;
          }

          if (proc_id != 0 && !is_bad_cand) {
            Double_t ff_out = 1.;
            try {
              ham.processEvent();
              ff_out = ham.getWeight("OutputFF");
            } catch (const exception& e) {
              cout << "WARN: HAMMER doesn't like candidate for reweighting: "
                   << num_of_evt << endl;
              cout << e.what() << endl;
              is_bad_cand = true;
            }

            if (!isnan(ff_out) && !isinf(ff_out) && !is_bad_cand) {
              w_ff_out = ff_out;
              num_of_evt_ham_ok += 1;
              ham_ok = true;
            }
          }
        }

        if (is_bad_cand) {
#ifndef SILENT
          // Print debug info
          cout << "========" << endl;
          cout << "True q2 (GeV): " << q2_true_out << endl;
          cout << "B meson ID: " << get_particle_name(b_id_fixed, db) << endl;
          cout << "D meson ID: " << get_particle_name(d_id, db) << endl;
          cout << "D daughter 0 ID: "
               << get_particle_name(D_daughter_id[D_lbl + "_GD0"], db) << endl;
          cout << "D daughter 1 ID: "
               << get_particle_name(D_daughter_id[D_lbl + "_GD1"], db) << endl;
          cout << "D daughter 2 ID: "
               << get_particle_name(D_daughter_id[D_lbl + "_GD2"], db) << endl;

          cout << "Is tau decay: " << *is_tau << endl;
          cout << "Current candidate index: " << num_of_evt << endl;

          // More detailed debug messages
          cout << "B meson 4-mom: " << print_p(part_B.p()) << endl;
          cout << "B meson inv.m: " << part_B.p().mass() << endl;
          cout << "D meson 4-mom: " << print_p(part_D.p()) << endl;
          cout << "D meson inv.m: " << part_D.p().mass() << endl;
          for (auto idx = 0; idx < part_D_daughters.size(); idx++) {
            cout << "D daughter idx "s + idx + " 4-mom: "
                 << print_p(part_D_daughters[idx].p()) << endl;
            cout << "D daughter idx "s + idx + " inv.m: "
                 << part_D_daughters[idx].p().mass() << endl;
          }

          if (*is_tau) {
            cout << "Tau 4-mom: " << print_p(part_L.p()) << endl;
            cout << "Tau inv.m: " << part_L.p().mass() << endl;
            cout << "anti-TauNu 4-mom: " << print_p(part_NuL.p()) << endl;
            cout << "anti-TauNu inv.m: " << part_NuL.p().mass() << endl;
            cout << "TauNu 4-mom: " << print_p(part_TauNuTau.p()) << endl;
            cout << "TauNu inv.m: " << part_TauNuTau.p().mass() << endl;
            cout << "Mu 4-mom: " << print_p(part_Mu.p()) << endl;
            cout << "Mu inv.m: " << part_Mu.p().mass() << endl;
            cout << "anti-MuNu 4-mom: " << print_p(part_TauNuMu.p()) << endl;
            cout << "anti-MuNu inv.mass: " << part_TauNuMu.p().mass() << endl;
          } else {
            cout << "Mu 4-mom: " << print_p(part_L.p()) << endl;
            cout << "Mu inv.m: " << part_L.p().mass() << endl;
            cout << "anti-MuNu 4-mom: " << print_p(part_NuL.p()) << endl;
            cout << "anti-MuNu inv.m: " << part_NuL.p().mass() << endl;
          }

          cout << "B meson HAMMER ID: " << part_B_idx << endl;
          cout << "D meson HAMMER ID: " << part_D_idx << endl;
          if (is_Dst) {
            cout << "D daughter HAMMER IDs:";
            for (auto idx : part_D_daughters_idx) cout << "  " << idx;
            cout << endl;
          }
          cout << "Lepton HAMMER ID: " << part_L_idx << endl;
          cout << "Lepton neutrino HAMMER ID: " << part_NuL_idx << endl;

#ifdef FORCE_MOMENTUM_CONSERVATION_HADRONIC
          cout << "Hadronic part known momentum: " << print_p(known_mom)
               << endl;
#endif

#ifdef RADIATIVE_CORRECTION
          cout << "Radiative photon correction info:" << endl;
          cout << buf_rad_corr;
#endif

          cout << "Check if invariant mass of each particle is non-negative: "
               << endl;
          for (auto v : part_m2_ok) cout << "  " << v;
          cout << endl;
          cout << "All particles have OK kinematics: " << all_parts_ok << endl;
          cout << "========" << endl;
#endif
        }
      }
    }

    eventNumber_out = *eventNumber;
    runNumber_out   = *runNumber;

    output_tree->Fill();
    num_of_evt += 1;
  }

  output_ntp->Write();

  // Cleanup ///////////////////////////////////////////////////////////////////
  cout << "Cleanups" << endl;
  delete output_tree;
  delete db;

  return RwRate{num_of_evt, num_of_evt_ham_ok};
}

//////////
// Main //
//////////

int main(int, char** argv) {
  auto    input_ntp  = new TFile(argv[1], "read");
  auto    output_ntp = new TFile(argv[2], "update");
  TString tree       = argv[3];
  TString run        = argv[4];

  Hammer::Hammer ham{};

  set_decays(ham);
  set_input_ff(ham, run);
  set_output_ff(ham);

  ham.setUnits("MeV");
  ham.setOptions("ProcessCalc: {CheckForNaNs: true}");
  ham.initRun();

  auto rate = reweight(input_ntp, output_ntp, tree, ham);

  cout << "Total number of candidates: " << get<0>(rate) << endl;
  cout << "Hammer reweighted candidates: " << get<1>(rate) << endl;
  cout << "Reweighted fraction: "
       << static_cast<float>(get<1>(rate)) / static_cast<float>(get<0>(rate))
       << endl;

  delete input_ntp;
  delete output_ntp;
}
