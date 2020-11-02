// Author: Yipeng Sun
// License: GPLv2
// Description: Validation of FF reweighting from ISGW2 -> CLN
// Last Change: Mon Nov 02, 2020 at 12:11 PM +0100

#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>

#include <ff_dstaunu.hpp>

using namespace std;

enum BMeson { Charged = 1, Neutral = 0 };

enum FFType { ISGW2 = 0, CLN = 1 };

// in GeV
const Double_t m_B0    = 5.2792;
const Double_t m_Bplus = 5.2792;
const Double_t m_Tau   = 1.7768;

// q2 distributions with a particular FF parameterization
TH1D q2_histo(BMeson b_type, FFType ff_type, Double_t m_lep, Int_t nbinsx,
              Double_t xlow, Double_t xup, const char* name,
              const char* title) {
  auto histo   = TH1D(name, title, nbinsx, xlow, xup);
  auto ff_calc = BToDstaunu{};

  ff_calc.SetMasses(b_type);

  for (auto bin = 1; bin <= nbinsx; bin++) {
    auto q2 = histo.GetXaxis()->GetBinCenter(bin);
    histo.SetBinContent(bin, ff_calc.Compute(q2, ff_type, m_lep * m_lep));
  }

  return histo;
}

int main(int, char** argv) {
  TFile* input_file = new TFile(argv[1], "read");
  string output_dir = argv[2];

  gROOT->SetBatch(kTRUE);

  auto canvas = new TCanvas("canvas", "FF validation", 4000, 3000);

  // Reference CLN vs reweighted
  auto histo_ref_cln_B0ToDstTauNu =
      q2_histo(BMeson::Neutral, FFType::CLN, m_Tau, 400, 3, 12, "histo_ref_cln",
               "Reference CLN");

  histo_ref_cln_B0ToDstTauNu.Draw();
  canvas->Update();

  canvas->Print((output_dir + "/ff_cln.png").c_str());

  delete input_file;
  delete canvas;
}
