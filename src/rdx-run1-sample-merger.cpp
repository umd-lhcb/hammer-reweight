// Author: Yipeng Sun
// License: GPLv2
// Description: FF reweighting for R(D(*)) run 1, step 1 ntuples.
// Based on:
//   https://github.com/ZishuoYang/my-hammer-reweighting/blob/master/Bc2JpsiMuNu.cc
// Last Change: Mon Oct 12, 2020 at 09:17 PM +0800

#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <string>

using namespace std;

void merge_data_weight(TFile* data_file, TFile* weight_file, TFile* output_file,
                       const char* data_tree_name   = "mc_dst_tau_aux",
                       const char* weight_tree_name = "mc_dst_tau_ff_w") {
  auto data_tree   = static_cast<TTree*>(data_file->Get(data_tree_name));
  auto weight_tree = static_cast<TTree*>(weight_file->Get(weight_tree_name));

  // Build index for weights so that the fit variables in data tree have correct
  // weights assigned to them
  weight_tree->BuildIndex("runNumber", "evnetNumber");
  data_tree->AddFriend(weight_tree);

  TTreeReader reader(data_tree);
  TTree       output(data_tree_name, data_tree_name);

  // Read input branches ///////////////////////////////////////////////////////
  // General
  TTreeReaderValue<ULong64_t> eventNumber(reader, "eventNumber");
  TTreeReaderValue<UInt_t>    runNumber(reader, "runNumber");

  // Fit variables
  TTreeReaderValue<Double_t> mm2(reader, "mm2");
  TTreeReaderValue<Double_t> q2(reader, "q2");
  TTreeReaderValue<Double_t> el(reader, "el");

  // MC weight
  TTreeReaderValue<Double_t> w_ff(reader, "w_ff");

  // Define output branches ////////////////////////////////////////////////////
  ULong64_t eventNumber_out;
  output.Branch("eventNumber", &eventNumber_out);
  UInt_t runNumber_out;
  output.Branch("runNumber", &runNumber_out);

  Double_t mm2_out;
  output.Branch("mm2", &mm2_out);
  Double_t q2_out;
  output.Branch("q2", &q2_out);
  Double_t el_out;
  output.Branch("el", &el_out);

  Double_t w_ff_out;
  output.Branch("w_ff", &w_ff_out);

  while (reader.Next()) {
    if (w_ff.Get() != nullptr) {
      eventNumber_out = *eventNumber;
      runNumber_out   = *runNumber;

      mm2_out = *mm2;
      q2_out  = *q2;
      el_out  = *el;

      w_ff_out = *w_ff;

      output.Fill();
    }
  }

  output_file->Write("", TObject::kOverwrite);
}

int main(int, char** argv) {
  TFile* data_file   = new TFile(argv[1], "read");
  TFile* weight_file = new TFile(argv[2], "read");
  TFile* output_file = new TFile(argv[3], "recreate");

  merge_data_weight(data_file, weight_file, output_file);

  delete data_file;
  delete weight_file;
  delete output_file;

  return 0;
}
