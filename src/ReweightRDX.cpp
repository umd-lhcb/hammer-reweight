// Author: Yipeng Sun
// Last Change: Tue Aug 03, 2021 at 08:41 PM +0200

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

typedef map<vector<Int_t>, unsigned long> DecayFreq;

// clang-format off
auto B_MESON = map<TString, TString>{
  {"TupleBminus/DecayTree", "b"},
  {"TupleB0/DecayTree", "b0"}
};

const auto LEGAL_B_MESON_IDS = vector<Int_t>{511, 521};
// clang-format on

const Double_t Q2_MIN = 100. * 100.;

void set_input_ff(Hammer::Hammer& ham, string decay_mode) {
  ;
  ;
}

void set_output_ff(Hammer::Hammer& ham, string decay_mode) {
  ;
  ;
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

  // Muon
  TTreeReaderValue<Double_t> mu_px(reader, b_meson + "_TrueMu_PX");
  TTreeReaderValue<Double_t> mu_py(reader, b_meson + "_TrueMu_PY");
  TTreeReaderValue<Double_t> mu_pz(reader, b_meson + "_TrueMu_PZ");
  TTreeReaderValue<Double_t> mu_pe(reader, b_meson + "_TrueMu_PE");

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

  Bool_t ham_ok;
  output_tree->Branch("flag_ham_ok", &ham_ok);

  ULong64_t eventNumber_out;
  output_tree->Branch("eventNumber", &eventNumber_out);
  UInt_t runNumber_out;
  output_tree->Branch("runNumber", &runNumber_out);

  // Setup HAMMER //////////////////////////////////////////////////////////////
  Hammer::Hammer   ham{};
  Hammer::IOBuffer ham_buf;

  ham.setUnits("MeV");
  ham.initRun();

  unsigned long num_of_evt           = 0l;
  unsigned long num_of_evt_w_b_meson = 0l;
  while (reader.Next()) {
    ham_ok = false;

    // Check if we have a legal B meson and q2 is large enough to produce a Mu
    if (find_in(LEGAL_B_MESON_IDS, TMath::Abs(*b_id)) && *q2 > Q2_MIN) {
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
        // add them if they are hadrons
        vector<Hammer::Particle> d_daughters{};

        for (auto suffix : vector<TString>{"_GD0", "_GD1", "_GD2"}) {
          auto part_name = D_lbl + suffix;
          auto part_id   = D_daughter_id[part_name];
          if (is_hadron(part_id)) {
            d_daughters.push_back(particle(D_daughter_mom[part_name + "PE"],
                                           D_daughter_mom[part_name + "PX"],
                                           D_daughter_mom[part_name + "PY"],
                                           D_daughter_mom[part_name + "PZ"],
                                           part_id));
          }
        }

        // HAMMER process
        Hammer::Process proc;

        auto part_B = particle(*b_pe, *b_px, *b_py, *b_pz,
                               B_id_fix(*b_id, D_cands[D_lbl]));
        auto part_D = particle(D_mom[D_lbl + "_PX"], D_mom[D_lbl + "_PY"],
                               D_mom[D_lbl + "_PZ"], D_mom[D_lbl + "_PE"],
                               D_cands[D_lbl]);

        num_of_evt_w_b_meson += 1;
      }
    }

    num_of_evt += 1;
  }

  // Cleanup ///////////////////////////////////////////////////////////////////
  delete output_tree;

  return RwRate{num_of_evt, num_of_evt_w_b_meson};
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
  cout << "Truth-matched candidates: " << get<1>(rate) << endl;
  cout << "Truth-matched fraction: "
       << static_cast<float>(get<0>(rate)) / static_cast<float>(get<1>(rate))
       << endl;

  delete input_ntp;
  delete output_ntp;
}
