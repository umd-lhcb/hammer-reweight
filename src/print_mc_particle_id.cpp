// Author: Yipeng Sun
// Last Change: Wed Jul 28, 2021 at 05:33 PM +0200

#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>

using namespace std;

typedef map<vector<Int_t>, unsigned long> DecayFreq;

// clang-format off
auto B_MESON = map<TString, TString>{
  {"TupleBminus/DecayTree", "b"},
  {"TupleB0/DecayTree", "b0"}
};

auto DECAY_NAMES = vector<string_view>{
  "B meson ID: ",
  "Mother ID: ",
  "First daughter ID: ",
  "  First G-daughter ID: ",
  "  Second G-daughter ID: ",
  "  Third G-daughter ID: ",
  "Second daughter meson ID: ",
  "  First G-daughter ID: ",
  "  Second G-daughter ID: ",
  "  Third G-daughter ID: ",
  "Third daughter meson ID: ",
  "  First G-daughter ID: ",
  "  Second G-daughter ID: ",
  "  Third G-daughter ID: "
};
// clang-format on

string get_particle_name(Int_t id, TDatabasePDG* db) {
  if (!id) return "None"s;

  auto abs_id = TMath::Abs(id);
  char buf[50];
  sprintf(buf, " (%d)", abs_id);
  auto str_id   = string(buf);
  auto particle = db->GetParticle(abs_id);

  if (particle != nullptr) return string(particle->GetName()) + buf;
  return "Unknown"s + buf;
}

void print_decay_freq(DecayFreq freq) {
  auto db = new TDatabasePDG();

  for (auto const& [key, val] : freq) {
    cout << "======" << endl;
    cout << "The following decay has " << val << " candidates." << endl;
    for (auto idx = 0; idx < key.size(); idx++) {
      if (key[idx]) {
        cout << DECAY_NAMES[idx] << get_particle_name(key[idx], db) << endl;
      }
    }
  }

  delete db;
}

DecayFreq print_id(TFile* input_file, TString tree, int modulo = 40) {
  TTreeReader reader(tree, input_file);
  auto        freq = DecayFreq{};

  if (B_MESON.find(tree) == B_MESON.end()) return freq;
  TString b_meson = B_MESON[tree];

  // B meson truth info
  TTreeReaderValue<Double_t> b_px(reader, b_meson + "_TRUEP_X");
  TTreeReaderValue<Double_t> b_py(reader, b_meson + "_TRUEP_Y");
  TTreeReaderValue<Double_t> b_pz(reader, b_meson + "_TRUEP_Z");
  TTreeReaderValue<Double_t> b_pe(reader, b_meson + "_TRUEP_E");
  TTreeReaderValue<Int_t>    b_id(reader, b_meson + "_TRUEID");
  TTreeReaderValue<Int_t>    mother_id(reader, b_meson + "_TrueHadron_M_ID");

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

  unsigned long counter = 0;
  while (reader.Next()) {
    if (!(counter % modulo)) {
      cout << "======" << endl;
      cout << "Muon E: " << *mu_pe << endl;
      cout << "Primary neutrino E: " << *anu_pe << endl;
      cout << "Tau E: " << *tau_pe << endl;
      cout << "Secondary Tau neutrino E: " << *nu_tau_pe << endl;
      cout << "Secondary Mu neutrino E: " << *anu_mu_pe << endl;
    }

    auto key = vector<Int_t>{};
    key.push_back(*b_id);
    key.push_back(*mother_id);
    key.push_back(*d_idx0_id);
    key.push_back(*d_idx0_gd0_id);
    key.push_back(*d_idx0_gd1_id);
    key.push_back(*d_idx0_gd2_id);
    key.push_back(*d_idx1_id);
    key.push_back(*d_idx1_gd0_id);
    key.push_back(*d_idx1_gd1_id);
    key.push_back(*d_idx1_gd2_id);
    key.push_back(*d_idx2_id);
    key.push_back(*d_idx2_gd0_id);
    key.push_back(*d_idx2_gd1_id);
    key.push_back(*d_idx2_gd2_id);

    if (freq.find(key) == freq.end())
      freq[key] = 1l;
    else
      freq[key] += 1;

    counter += 1;
  }

  return freq;
}

int main(int, char** argv) {
  auto    ntp       = new TFile(argv[1], "read");
  TString tree_name = argv[2];
  auto    db        = new TDatabasePDG();

  auto freq = print_id(ntp, tree_name);
  print_decay_freq(freq);

  delete ntp;
  delete db;
}
