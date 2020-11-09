// Author: Yipeng Sun
// License: GPLv2
// Description: Validation of FF reweighting from ISGW2 -> CLN
// Last Change: Mon Nov 09, 2020 at 01:23 AM +0100

#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <ff_dstaunu.hpp>

using namespace std;

enum BMeson { Charged = 1, Neutral = 0 };

enum FFType { ISGW2 = 0, CLN = 1 };

// in GeV
const Double_t m_B0    = 5.2792;
const Double_t m_Bplus = 5.2792;
const Double_t m_Tau   = 1.7768;

// q2 distributions with a particular FF parameterization
TH1D q2_histo(BMeson b_type, FFType ff_type, Double_t m_lep, const char* name,
              const char* title, Int_t nbinsx, Double_t xlow, Double_t xup,
              Option_t* scale_opt = "width") {
  auto histo   = TH1D(name, title, nbinsx, xlow, xup);
  auto ff_calc = BToDstaunu{};

  ff_calc.SetMasses(b_type);

  for (auto bin = 1; bin <= nbinsx; bin++) {
    auto     q2      = histo.GetBinCenter(bin);
    Double_t q2_dist = ff_calc.Compute(q2, ff_type, m_lep);

    histo.SetBinContent(bin, q2_dist);
  }

  histo.Scale(1 / histo.Integral(scale_opt));
  cout << "Histogram " << title << " has been normalized to "
       << histo.Integral(scale_opt) << endl;

  return histo;
}

TH1D fill_histo(TTree* tree, const char* branch, const char* name,
                const char* title, Double_t nbinsx, Double_t xlow,
                Double_t xup) {
  auto histo = TH1D(name, title, nbinsx, xlow, xup);

  Double_t val;
  tree->SetBranchAddress(branch, &val);

  for (Long64_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    histo.Fill(val);
  }

  return histo;
}

TH1D fill_histo(TTree* tree, const char* branch, const char* weight,
                const char* name, const char* title, Double_t nbinsx,
                Double_t xlow, Double_t xup) {
  auto histo = TH1D(name, title, nbinsx, xlow, xup);

  Double_t val, val_w;
  tree->SetBranchAddress(branch, &val);
  tree->SetBranchAddress(weight, &val_w);

  for (Long64_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    histo.Fill(val, val_w);
  }

  return histo;
}

template <class T>
void debug_histo(T histo, Option_t* scale_opt = "") {
  cout << "Histogram " << histo->GetName() << "has an integral of "
       << histo->Integral(scale_opt) << endl;
  cout << "Maximum of " << histo->GetName() << ": " << histo->GetMaximum()
       << endl;
}

template <class T>
Double_t rescale_histos(T histo1, T histo2) {
  auto ratio = histo2->GetMaximum() / histo1->GetMaximum();
  histo1->Scale(ratio);
  return ratio;
}

int main(int, char** argv) {
  TFile* data_file   = new TFile(argv[1], "read");
  TFile* weight_file = new TFile(argv[2], "read");
  string output_dir  = argv[3];

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TTree* data_tree   = data_file->Get<TTree>("dst_iso");
  TTree* weight_tree = weight_file->Get<TTree>("mc_dst_tau_ff_w");

  weight_tree->BuildIndex("runNumber", "eventNumber");
  data_tree->AddFriend(weight_tree);

  // Reference CLN
  auto histo_ref_cln_B0ToDstTauNu =
      q2_histo(BMeson::Neutral, FFType::CLN, m_Tau, "CLN", "Reference CLN", 200,
               2.5, 12);

  histo_ref_cln_B0ToDstTauNu.SetLineWidth(2);
  histo_ref_cln_B0ToDstTauNu.SetLineColor(kRed);

  // Reference ISGW2
  auto histo_ref_isgw2_B0ToDstTauNu =
      q2_histo(BMeson::Neutral, FFType::ISGW2, m_Tau, "ISGW2",
               "Reference ISGW2", 200, 2.5, 12);

  histo_ref_isgw2_B0ToDstTauNu.SetLineWidth(2);
  histo_ref_isgw2_B0ToDstTauNu.SetLineColor(kBlue);

  // Original ISGW2
  auto histo_orig =
      fill_histo(data_tree, "q2_true", "q2_orig", "q2 original", 70, 2.5, 12);
  histo_orig.Scale(1 / histo_orig.Integral("width"));
  debug_histo(&histo_orig, "width");
  // auto scale_ratio = rescale_histos(&histo_orig,
  // &histo_ref_isgw2_B0ToDstTauNu);

  // We want to align the maximum of the real data points with its reference
  // distribution
  histo_orig.Scale(histo_ref_isgw2_B0ToDstTauNu.GetMaximum() /
                   histo_orig.GetMaximum());

  histo_orig.SetLineWidth(4);
  histo_orig.SetLineColor(kGreen);

  // Reweighted CLN
  auto histo_reweighted =
      fill_histo(data_tree, "q2_true", "w_ff", "q2_reweighted", "q2 reweighted",
                 70, 2.5, 12);
  histo_reweighted.Scale(1 / histo_reweighted.Integral("width"));
  debug_histo(&histo_reweighted, "width");
  // histo_reweighted.Scale(scale_ratio);

  histo_reweighted.SetLineWidth(4);
  histo_reweighted.SetLineColor(kOrange);

  // Plot
  auto canvas = new TCanvas("canvas", "FF validation", 4000, 3000);
  histo_ref_cln_B0ToDstTauNu.Draw("hist C");
  histo_ref_isgw2_B0ToDstTauNu.Draw("same hist C");
  histo_orig.Draw("same hist");
  histo_reweighted.Draw("same hist");

  canvas->BuildLegend();
  canvas->Update();
  canvas->Print((output_dir + "/validate_ff.png").c_str());

  delete canvas;

  delete data_tree;
  delete weight_tree;

  delete data_file;
  delete weight_file;
}
