// Author: Yipeng Sun
// License: GPLv2
// Description: Validation of FF reweighting from ISGW2 -> CLN
// Last Change: Mon Nov 02, 2020 at 03:42 PM +0100

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
TH1D q2_histo(BMeson b_type, FFType ff_type, Double_t m_lep,
              Long64_t normalization, const char* name, const char* title,
              Int_t nbinsx, Double_t xlow, Double_t xup) {
  auto histo   = TH1D(name, title, nbinsx, xlow, xup);
  auto ff_calc = BToDstaunu{};

  ff_calc.SetMasses(b_type);

  for (auto bin = 1; bin <= nbinsx; bin++) {
    auto q2      = histo.GetXaxis()->GetBinCenter(bin);
    auto q2_dist = ff_calc.Compute(q2, ff_type, m_lep);

    histo.SetBinContent(bin, q2_dist);
  }

  histo.ComputeIntegral();
  cout << title << " has a " << *(histo.GetIntegral()) << " integral." << endl;

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

  auto normalization = data_tree->GetEntries();

  // Reference CLN
  auto histo_ref_cln_B0ToDstTauNu =
      q2_histo(BMeson::Neutral, FFType::CLN, m_Tau, normalization, "CLN",
               "Reference CLN", 800, 2.5, 12);

  histo_ref_cln_B0ToDstTauNu.SetLineWidth(2);
  histo_ref_cln_B0ToDstTauNu.SetLineColor(kRed);

  // Reference ISGW2
  auto histo_ref_isgw2_B0ToDstTauNu =
      q2_histo(BMeson::Neutral, FFType::ISGW2, m_Tau, normalization, "ISGW2",
               "Reference ISGW2", 800, 2.5, 12);
  histo_ref_isgw2_B0ToDstTauNu.SetLineWidth(2);
  histo_ref_isgw2_B0ToDstTauNu.SetLineColor(kBlue);

  // Original ISGW2
  auto histo_orig =
      fill_histo(data_tree, "q2", "q2_orig", "q2 original", 80, 2.5, 12);
  histo_orig.SetLineColor(kMagenta);
  histo_orig.SetLineWidth(0);

  // Reweighted CLN
  auto histo_reweighted = fill_histo(data_tree, "q2", "w_ff", "q2_reweighted",
                                     "q2 reweighted", 80, 2.5, 12);

  // Original plot
  auto canvas = new TCanvas("canvas", "FF validation", 4000, 3000);
  histo_ref_cln_B0ToDstTauNu.Draw();
  histo_ref_isgw2_B0ToDstTauNu.Draw("same");
  histo_orig.Draw("same *H");

  // Legend
  auto legend = new TLegend(0.1, 0.8, 0.3, 0.9);
  legend->SetTextSize(0.03);
  legend->AddEntry(&histo_ref_cln_B0ToDstTauNu, "CLN", "f");
  legend->AddEntry(&histo_ref_isgw2_B0ToDstTauNu, "ISGW2", "f");
  legend->AddEntry(&histo_orig, "Original", "f");
  legend->Draw();
  // canvas->BuildLegend();

  canvas->Update();
  canvas->Print((output_dir + "/ff_cln_original.png").c_str());
  canvas->Clear();

  // Reweighted plot
  histo_reweighted.Draw("*H");

  canvas->Update();
  canvas->Print((output_dir + "/ff_cln_reweighted.png").c_str());

  delete canvas;
  delete legend;

  delete data_tree;
  delete weight_tree;

  delete data_file;
  delete weight_file;
}
