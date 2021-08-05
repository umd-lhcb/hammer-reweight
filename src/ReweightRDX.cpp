// Author: Yipeng Sun
// Last Change: Thu Aug 05, 2021 at 07:45 PM +0200

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <Hammer/Hammer.hh>
#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>
#include <Hammer/Process.hh>
#include <Hammer/Tools/HammerRoot.hh>

using namespace std;

///////////////////
// Configurables //
///////////////////

#define SILENT
#define FORCE_MOMENTUM_CONSERVATION

typedef map<vector<Int_t>, unsigned long> DecayFreq;

const Double_t Q2_MIN = 100. * 100.;

// clang-format off
auto B_MESON = map<TString, TString>{
  {"TupleBminus/DecayTree", "b"},
  {"TupleB0/DecayTree", "b0"}
};

const auto LEGAL_B_MESON_IDS = vector<Int_t>{511, 521};

void set_input_ff(Hammer::Hammer& ham) {
  ham.setFFInputScheme({
    {"BD*", "ISGW2"},
    {"BD", "ISGW2"}
  });
}

void set_output_ff(Hammer::Hammer& ham) {
  ham.addFFScheme("OutputFF", {
      {"BD*", "BGL"},
      {"BD", "BGL"}
  });
}
// clang-format on

void set_decays(Hammer::Hammer& ham) {
  ham.includeDecay("BD*TauNu");
  ham.includeDecay("BD*MuNu");

  ham.includeDecay("BDTauNu");
  ham.includeDecay("BDMuNu");
}

/////////////////////
// General helpers //
/////////////////////

template <template <typename, typename> class Iterable, typename T,
          typename Allocator>
Bool_t find_in(Iterable<T, Allocator> iter, T elem) {
  if (find(iter.begin(), iter.end(), elem) != iter.end()) return true;
  return false;
}

vector<string> split(const string& s, char delim) {
  stringstream   ss(s);
  string         item;
  vector<string> elems;

  while (getline(ss, item, delim)) {
    elems.push_back(move(item));
  }

  return elems;
}

TString dirname(string s) {
  auto splitted = split(s, '/');
  if (splitted.size() == 1) return TString("");

  TString dir = "";
  for (int idx = 0; idx < splitted.size() - 1; idx++) {
    dir += splitted[idx];
    if (idx < splitted.size() - 2) dir += "/";
  }

  return dir;
}

TString basename(string s) { return TString(split(s, '/').back()); }

Int_t digit_is(Int_t num, Int_t digit, Int_t base = 10) {
  Int_t fac = num / TMath::Power(base, digit - 1);
  return fac % 10;
}

TString print_p(const Hammer::FourMomentum& p) {
  char tmp[80];
  sprintf(tmp, "%.2f, %.2f, %.2f, %.2f", p.E(), p.px(), p.py(), p.pz());
  return TString(tmp);
}

/////////////////////////
// Particle ID helpers //
/////////////////////////

typedef pair<Bool_t, TString>  DMesonPack;
typedef map<TString, Int_t>    PartIdMap;
typedef map<TString, Double_t> PartMomMap;

DMesonPack is_D_meson(const PartIdMap parts) {
  for (const auto [key, val] : parts) {
    auto id = TMath::Abs(val);
    if (digit_is(id, 3) == 4) return DMesonPack{true, key};
  }

  return DMesonPack{false, "none"};
}

Bool_t is_hadron(const Int_t id) {
  if (TMath::Abs(id) > 100) return true;
  return false;
}

Int_t B_id_fix(const Int_t B_id, const Int_t D_id) {
  if (B_id * D_id > 0) return -B_id;
  return B_id;
}

int Mu_id(Int_t mu_id, Int_t true_id = 13) {
  Int_t sign = mu_id / TMath::Abs(mu_id);
  return sign * true_id;
}

int Tau_id(Int_t mu_id) { return Mu_id(mu_id, 15); }

int Nu_id(Int_t mu_id, Bool_t is_tau) {
  if (is_tau) return Mu_id(mu_id, -16);
  return Mu_id(mu_id, -14);
}

int Tau_NuMu_id(Int_t mu_id) { return Mu_id(mu_id, -14); }

int Tau_NuTau_id(Int_t mu_id) { return Mu_id(mu_id, 16); }

////////////////////////////
// HAMMER-related helpers //
////////////////////////////

auto particle(Double_t pe, Double_t px, Double_t py, Double_t pz, Int_t pid) {
  auto four_mom = Hammer::FourMomentum(pe, px, py, pz);
  auto part_id  = static_cast<Hammer::PdgId>(pid);

  return Hammer::Particle(four_mom, part_id);
}

/////////////////
// Reweighting //
/////////////////

typedef pair<unsigned long, unsigned long> RwRate;

RwRate reweight(TFile* input_ntp, TFile* output_ntp, TString tree) {
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
  output_tree->Branch("flag_ham_ok", &ham_ok);

  // Setup HAMMER //////////////////////////////////////////////////////////////
  Hammer::Hammer ham{};

  set_decays(ham);
  set_input_ff(ham);
  set_output_ff(ham);

  ham.setUnits("MeV");
  ham.initRun();

  unsigned long num_of_evt        = 0l;
  unsigned long num_of_evt_ham_ok = 0l;
  while (reader.Next()) {
    ham_ok   = false;
    w_ff_out = 1.;

    // Check if we have a legal B meson and q2 is large enough to produce a Mu
    if (find_in(LEGAL_B_MESON_IDS, TMath::Abs(*b_id)) && *q2 > Q2_MIN &&
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
          if (is_hadron(part_id)) {  // Remove photons otherwise HAMMER get
                                     // stuck for unknown reason
            part_D_daughters.push_back(
                particle(D_daughter_mom[part_name + "_PE"],
                         D_daughter_mom[part_name + "_PX"],
                         D_daughter_mom[part_name + "_PY"],
                         D_daughter_mom[part_name + "_PZ"], part_id));
          }
        }

        // HAMMER process
        Hammer::Process proc;

        auto b_id_fixed = B_id_fix(*b_id, D_cands[D_lbl]);
        auto part_B     = particle(*b_pe, *b_px, *b_py, *b_pz, b_id_fixed);
        auto part_D     = particle(D_mom[D_lbl + "_PX"], D_mom[D_lbl + "_PY"],
                               D_mom[D_lbl + "_PZ"], D_mom[D_lbl + "_PE"],
                               D_cands[D_lbl]);
        Hammer::Particle part_L, part_NuL, part_TauNuTau, part_TauNuMu, part_Mu;
        if (*is_tau)
          part_L = particle(*tau_pe, *tau_px, *tau_py, *tau_pz, Tau_id(*mu_id));
        else
          part_L = particle(*mu_pe, *mu_px, *mu_py, *mu_pz, Mu_id(*mu_id));

        part_Mu       = particle(*mu_pe, *mu_px, *mu_py, *mu_pz, Mu_id(*mu_id));
        part_NuL      = particle(*anu_pe, *anu_px, *anu_py, *anu_pz,
                            Nu_id(*mu_id, *is_tau));
        part_TauNuTau = particle(*nu_tau_pe, *nu_tau_px, *nu_tau_py, *nu_tau_pz,
                                 Tau_NuTau_id(*mu_id));
        part_TauNuMu  = particle(*anu_mu_pe, *anu_mu_px, *anu_mu_py, *anu_mu_pz,
                                Tau_NuMu_id(*mu_id));

        auto part_B_idx = proc.addParticle(part_B);
        auto part_D_idx = proc.addParticle(part_D);

        Hammer::ParticleIndices part_D_daughters_idx{};
        for (const auto part : part_D_daughters) {
          part_D_daughters_idx.push_back(proc.addParticle(part));
        }

#ifdef FORCE_MOMENTUM_CONSERVATION
        // Leptonic part
        part_L.setMomentum(part_B.p() - part_D.p() - part_NuL.p());

        if (*is_tau)
          part_Mu.setMomentum(part_L.p() - part_TauNuMu.p() -
                              part_TauNuTau.p());

        // Hadronic part
        Hammer::FourMomentum known_mom{0, 0, 0, 0};
        for (auto idx = 0; idx < part_D_daughters.size() - 1; idx++) {
          known_mom += part_D_daughters[idx].p();
        }

        part_D_daughters[part_D_daughters.size() - 1].setMomentum(part_D.p() -
                                                                  known_mom);
#endif

        // Always add the primary leptons
        auto part_L_idx   = proc.addParticle(part_L);
        auto part_NuL_idx = proc.addParticle(part_NuL);

        proc.addVertex(part_B_idx, {part_D_idx, part_L_idx, part_NuL_idx});
        proc.addVertex(part_D_idx, part_D_daughters_idx);

        if (*is_tau) {
          auto part_Mu_idx       = proc.addParticle(part_Mu);
          auto part_TauNuMu_idx  = proc.addParticle(part_TauNuMu);
          auto part_TauNuTau_idx = proc.addParticle(part_TauNuTau);

          proc.addVertex(part_L_idx,
                         {part_Mu_idx, part_TauNuMu_idx, part_TauNuTau_idx});
        }

        ham.initEvent();
        auto proc_id = ham.addProcess(proc);

#ifndef SILENT
        // Print debug info for first possibly legal candidate
        cout << "========" << endl;
        cout << "B meson ID: " << b_id_fixed << endl;
        cout << "D meson ID: " << D_cands[D_lbl] << endl;
        cout << "Is tau decay: " << *is_tau << endl;
        cout << "HAMMER process ID: " << proc_id << endl;
        cout << "Current candidate index: " << num_of_evt << endl;

        // More detailed debug messages
        cout << "B meson 4-mom: " << print_p(part_B.p()) << endl;
        cout << "D meson 4-mom: " << print_p(part_D.p()) << endl;
        for (auto idx = 0; idx < part_D_daughters.size(); idx++) {
          cout << "D daughter idx "s + idx + " 4-mom: "
               << print_p(part_D_daughters[idx].p()) << endl;
        }

        if (*is_tau) {
          cout << "Tau 4-mom: " << print_p(part_L.p()) << endl;
          cout << "anti-TauNu 4-mom: " << print_p(part_NuL.p()) << endl;
          cout << "TauNu 4-mom: " << print_p(part_TauNuTau.p()) << endl;
          cout << "Mu 4-mom: " << print_p(part_Mu.p()) << endl;
          cout << "anti-MuNu 4-mom: " << print_p(part_TauNuMu.p()) << endl;
        } else {
          cout << "Mu 4-mom: " << print_p(part_L.p()) << endl;
          cout << "anti-MuNu 4-mom: " << print_p(part_NuL.p()) << endl;
        }

        cout << "B meson HAMMER ID: " << part_B_idx << endl;
        cout << "D meson HAMMER ID: " << part_D_idx << endl;
        cout << "D daughter HAMMER IDs:";
        for (auto idx : part_D_daughters_idx) cout << "  " << idx;
        cout << endl;
        cout << "Lepton HAMMER ID: " << part_L_idx << endl;
        cout << "Lepton neutrino HAMMER ID: " << part_NuL_idx << endl;

        cout << "Hadronic part known momentum: " << print_p(known_mom) << endl;
#endif

        if (proc_id != 0) {
          ham_ok = true;
          ham.processEvent();
          w_ff_out = ham.getWeight("OutputFF");
          num_of_evt_ham_ok += 1;
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

  return RwRate{num_of_evt, num_of_evt_ham_ok};
}

//////////
// Main //
//////////

int main(int, char** argv) {
  auto    input_ntp  = new TFile(argv[1], "read");
  auto    output_ntp = new TFile(argv[2], "recreate");
  TString tree       = argv[3];

  auto rate = reweight(input_ntp, output_ntp, tree);

  cout << "Total number of candidates: " << get<0>(rate) << endl;
  cout << "Hammer reweighted candidates: " << get<1>(rate) << endl;
  cout << "Reweighted fraction: "
       << static_cast<float>(get<1>(rate)) / static_cast<float>(get<0>(rate))
       << endl;

  delete input_ntp;
  delete output_ntp;
}
