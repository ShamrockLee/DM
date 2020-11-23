/// dumpCompare.C by Yueh-Shun Li (ShamrockLee)
// Based on the dumpComparePDF.C by Prof. Shin-Shan Eiko Yu (syuvivida)
// The function chi2NbinsCompare is borrowed directly
// and the function dumpCompareMerged is a reimplementation
// of the function dumpComparePDF
// that passes the arguments to dumpCompare

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TSystem.h>

#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>

void chi2NbinsCompare(const TH1* h1, const TH1* h2, double& chi2, int& nbins,
                      int binLo = -1, int binHi = -1) {
  printf("Chi2 Calculation: \n");

  chi2 = 0;
  nbins = 0;

  int nBins = 0;
  double evtThresh = 1e-6;

  double h1Err, h1Value, h2Err, h2Value, lowEdge, highEdge, binChi2;

  double binWidth = h1->GetXaxis()->GetBinWidth(1);

  if (binLo < 0 || binHi < 0) {
    binLo = 1;
    binHi = h1->GetNbinsX();
  }

  for (int i = binLo; i <= binHi; i++) {
    h1Value = h1->GetBinContent(i);
    h1Err = h1->GetBinError(i);

    h2Value = h2->GetBinContent(i);
    h2Err = h2->GetBinError(i);

    if (h1Value < evtThresh && h2Value < evtThresh) continue;
    if (h1Err < 1e-6 && h2Err < 1e-6) continue;
    lowEdge = h1->GetXaxis()->GetBinLowEdge(i);
    highEdge = h1->GetXaxis()->GetBinUpEdge(i);

    binChi2 = (h1Value - h2Value) / sqrt(h1Err * h1Err + h2Err * h2Err);
    binChi2 *= binChi2;

    chi2 += binChi2;

    // printf( "%d) [%d, %d] [%f, %f] h1: %f h2: %f h1Err: %f h2Err: %f chi2: %f
    // total: %f \n ", 	    nBins, i, i, lowEdge, highEdge,h1Value, h2Value,
    // h1Err, h2Err, binChi2, chi2);

    nBins++;
  }
  int NDF = nBins;

  printf("Fit chi2/NDF = %f/%d, prob: %f\n", chi2, NDF,
         TMath::Prob(chi2, NDF) * 100);
  nbins = NDF;
}

void drawHistsCompared(TH1* h1, TH1* h2, const TString endfix1, const TString endfix2, Bool_t normalize = true) {
  Int_t nEntries1 = h1->GetEntries();
  Int_t nEntries2 = h2->GetEntries();
  TLegend* leg = new TLegend(0.3333, 0.7027, 0.8333, 0.9023);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);

  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);
  h1->SetLineWidth(3);
  h1->SetLineColor(4);

  h2->SetMarkerStyle(25);
  h2->SetMarkerSize(1);
  h2->SetLineWidth(3);
  h2->SetLineColor(2);

  gStyle->SetOptStat(0);
  h1->Sumw2();
  h2->Sumw2();
  if (normalize) {
    if (nEntries1) h1->Scale(1.0 / h1->Integral());
    if (nEntries2) h2->Scale(1.0 / h2->Integral());
  }

  Double_t chi2 = 0.;
  Int_t nbins = 0;
  if (nEntries1 && nEntries2) {
    Float_t max =
    h1->GetMaximum() > h2->GetMaximum() ? h1->GetMaximum() : h2->GetMaximum();
    h1->SetMaximum(1.1 * max);
    h2->SetMaximum(1.1 * max);

    chi2NbinsCompare(h1, h2, chi2, nbins, 1, h1->GetNbinsX());
  }

  h1->Draw("hist");
  h2->Draw("histsame");

  leg->Clear();
  leg->SetHeader("");
  if (nEntries1 && nEntries2) {
    leg->AddEntry((TObject*)0, Form("#chi^{2}/NDF=%.1f/%d", chi2, nbins), "");
    leg->AddEntry((TObject*)0,
                  Form("#chi^{2} Prob=%.2f", TMath::Prob(chi2, nbins)), "");
    leg->AddEntry((TObject*)0, "", "");
  }

  leg->AddEntry(h1, endfix1.Data(), "l");
  leg->AddEntry(h2, endfix2.Data(), "l");
  leg->Draw("same");
}

void refgetIdxHistOpenableRecent(TDirectory *tdir, std::vector<TString>& vName, std::vector<UInt_t>& vIdx) {
  UInt_t nKeys = tdir->GetNkeys();
  vName.clear();
  vIdx.clear();
  std::vector<Short_t> vCycle(nKeys);
  UInt_t iKey=0;
  for (const auto&& keyRaw: *tdir->GetListOfKeys()) {
    Short_t cycle = 0;
    if (keyRaw != nullptr) {
      TKey *key = (TKey *) keyRaw;
      vName.push_back(key->GetName());
      cycle = key->GetCycle();
      if (key->ReadObj()->IsA()->InheritsFrom( "TH1" )) {
        vIdx.push_back(iKey);
      }
    } else {
      vName.push_back("");
    }
    iKey++;
  }
  std::stable_sort(vIdx.begin(), vIdx.end(), [&vName, &vCycle](UInt_t idx1, UInt_t idx2)->Bool_t{
    const TString &name1 = vName[idx1];
    const TString &name2 = vName[idx2];
    return (name1 < name2) || (name1 == name2 && vCycle[idx1] < vCycle[idx2]);
  });
  for (UInt_t iIdx=vIdx.size()-1; iIdx>1; iIdx--) {
    if (vName[vIdx[iIdx]] == vName[vIdx[iIdx-1]]) {
      vIdx.erase(vIdx.begin()+(iIdx-1));
    }
  }
}

void dumpCompare(TDirectory *tdir1, TDirectory *tdir2, Bool_t normalize, TString endfix1, TString endfix2, std::function<void (TCanvas*, TString, UInt_t)> funPrint, std::function<void (TCanvas*, UInt_t)> funEndPrint, Int_t indexOrder=1) {
  const UInt_t nTDir = 2;
  TDirectory *arrTDirIn[nTDir] = {tdir1, tdir2};
  std::vector<TString> vNameMain;
  std::vector<UInt_t> vIdxMain;
  UInt_t indexChooseLeaves = (indexOrder == 0) ? 0 : indexOrder - 1;
  refgetIdxHistOpenableRecent(arrTDirIn[indexChooseLeaves], vNameMain, vIdxMain);
  std::vector<UInt_t> vpermIdx(vIdxMain.size());
  std::iota(vpermIdx.begin(), vpermIdx.end(), 0);
  if (indexOrder == 0) {
    std::sort(vpermIdx.begin(), vpermIdx.end(), [&vIdxMain](UInt_t i, UInt_t j)->Bool_t {return vIdxMain[i] < vIdxMain[j];});
  }
  TCanvas *c1 = new TCanvas;
  UInt_t nPage = 0;
  for (UInt_t perm:vpermIdx) {
    TString nameHist = vNameMain[perm];
    TH1 *arrHist[nTDir];
    Bool_t isLeafCommon = true;
    for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
      TObject *histRaw = arrTDirIn[iTDir]->Get(nameHist);
      if (histRaw == nullptr || !histRaw->IsA()->InheritsFrom("TH1")) {
        isLeafCommon = false;
        std::cout << "Leaf " << nameHist << " not found in " << arrTDirIn[iTDir]->GetName() << "(" << arrTDirIn[iTDir] << ")" << std::endl;
      } else {
        arrHist[iTDir] = (TH1 *) histRaw;
      }
    }
    if (!isLeafCommon) continue;
    drawHistsCompared(arrHist[0], arrHist[1], endfix1, endfix2, normalize);
    funPrint(c1, nameHist, nPage);
    nPage++;
  }
  funEndPrint(c1, nPage);
}

void dumpCompareMerged(TDirectory *tdir1, TDirectory *tdir2, Bool_t normalize, TString endfix1, TString endfix2, TString pathDir="fig", TString nameFile="", TString optionPrint="pdf", Int_t indexOrder=1) {
  // const UInt_t nTDir = 2;
  // TDirectory *arrTDirIn[nTDir] = {tdir1, tdir2};
  // std::vector<TString> vNameMain;
  // std::vector<UInt_t> vIdxMain;
  // UInt_t indexChooseLeaves = (indexOrder == 0) ? 0 : indexOrder - 1;
  // refgetIdxHistOpenableRecent(arrTDirIn[indexChooseLeaves], vNameMain, vIdxMain);
  // std::vector<UInt_t> vpermIdx(vIdxMain.size());
  // std::iota(vpermIdx.begin(), vpermIdx.end(), 0);
  // if (indexOrder == 0) {
  //   std::sort(vpermIdx.begin(), vpermIdx.end(), [&vIdxMain](UInt_t i, UInt_t j)->Bool_t {return vIdxMain[i] < vIdxMain[j];});
  // }
  if (nameFile == "") {
    nameFile = endfix1 + "_" + endfix2 + "." + optionPrint;
    if (nameFile.EndsWith("+")) {
      nameFile.Resize(nameFile.Length()-1);
    }
  }
  
  if (pathDir.Length() >= 2 && pathDir.EndsWith("/")) {
    pathDir.Resize(pathDir.Length()-1);
  }
  TString filenamePrint = (pathDir == "" ? "" : pathDir + "/") + nameFile;
  TCanvas *c1 = new TCanvas;
  std::function<void (TCanvas *, TString, UInt_t)> funPrint = [filenamePrint](TCanvas *c, TString nameHist, UInt_t iPage) {
    if (iPage == 0) {
      c->Print(filenamePrint + "(");
    } else {
      c->Print(filenamePrint);
    }
  };
  std::function<void (TCanvas *, UInt_t)> funEndPrint = [filenamePrint](TCanvas *c, UInt_t nPage) {
    c->Print(filenamePrint + ")");
  };
  dumpCompare(tdir1, tdir2, normalize, endfix1, endfix2, funPrint, funEndPrint, indexOrder);
  // UInt_t nPage = 0;
  // for (UInt_t perm:vpermIdx) {
  //   TString nameLeaf = vNameMain[perm];
  //   TH1 *arrHist[nTDir];
  //   Bool_t isLeafCommon = true;
  //   for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
  //     TObject *histRaw = arrTDirIn[iTDir]->Get(nameLeaf);
  //     if (histRaw == nullptr || !histRaw->IsA()->InheritsFrom("TH1")) {
  //       isLeafCommon = false;
  //       std::cout << "Leaf " << nameLeaf << " not found in " << arrTDirIn[iTDir]->GetName() << "(" << arrTDirIn[iTDir] << ")" << std::endl;
  //     } else {
  //       arrHist[iTDir] = (TH1 *) histRaw;
  //     }
  //   }
  //   if (!isLeafCommon) continue;
  //   drawHistsCompared(arrHist[0], arrHist[1], endfix1, endfix2, normalize);
  //   c1->Print(filenamePrint + ((nPage == 0) ? "(" : ""), optionPrint);
  // }
  // c1->Print(filenamePrint + ")", optionPrint);
}

void dumpCompareMerged(TString pathFileIn1, TString pathFileIn2, Bool_t normalize=true, TString endfix1="", TString endfix2="", TString nameDir="fig", TString nameFile="", TString optionPrint="pdf", Int_t indexOrder=1) {
  const UInt_t nTDir = 2;
  TString arrPathFileIn[nTDir] = {pathFileIn1, pathFileIn2};
  TString arrEndfix[nTDir] = {endfix1, endfix2};
  TFile *arrTFileIn[nTDir];
  for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
    if (arrEndfix[iTDir] == "") {
      arrEndfix[iTDir]=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"${test%%.root}\"", arrPathFileIn[iTDir].Data()));
      arrTFileIn[iTDir] = TFile::Open(arrPathFileIn[iTDir], "read");
    }
  }
  dumpCompareMerged(arrTFileIn[0], arrTFileIn[1], normalize, arrEndfix[0], arrEndfix[1], nameDir, nameFile, optionPrint, indexOrder);
  for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
    arrTFileIn[iTDir]->Close();
  }
}

void dumpCompareSeperated(TDirectory *tdir1, TDirectory *tdir2, Bool_t normalize, TString endfix1, TString endfix2, TString pathDir="", std::function<TString (TString, UInt_t)> funNameFile=nullptr, TString optionPrint="svg", Int_t indexOrder=1) {
  if (pathDir == "") {
    pathDir = endfix1 + "_" + endfix2;
  }
  if (funNameFile == nullptr) {
    TString suffixFromOption = optionPrint;
    if (suffixFromOption.EndsWith("+")) {
      suffixFromOption.Resize(suffixFromOption.Length()-1);
    }
    funNameFile = [suffixFromOption](TString nameHist, UInt_t iPage)->TString{return nameHist + "." + suffixFromOption;};
  }
  if (pathDir.Length() >= 2 && pathDir.EndsWith("/")) {
    pathDir.Resize(pathDir.Length()-1);
  }
  gSystem->mkdir(pathDir);
  std::function<void (TCanvas *, TString, UInt_t)> funPrint = [pathDir, endfix1, endfix2, funNameFile, optionPrint](TCanvas *c, TString nameHist, UInt_t iPage) {
    c->Print((pathDir=="" ? "" : pathDir + "/") + funNameFile(nameHist, iPage), optionPrint);
  };
  std::function<void (TCanvas *, UInt_t)> funEndPrint = [](TCanvas *c, UInt_t nPage){};
  dumpCompare(tdir1, tdir2, normalize, endfix1, endfix2, funPrint, funEndPrint, indexOrder);
}

void dumpCompareSeperated(TString pathFileIn1, TString pathFileIn2, Bool_t normalize=true, TString endfix1="", TString endfix2="", TString pathDir="", std::function<TString (TString, UInt_t)> funNameFile=nullptr, TString optionPrint="svg", Int_t indexOrder=1) {
  const UInt_t nTDir = 2;
  TString arrPathFileIn[nTDir] = {pathFileIn1, pathFileIn2};
  TString arrEndfix[nTDir] = {endfix1, endfix2};
  TFile *arrTFileIn[nTDir];
  for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
    if (arrEndfix[iTDir] == "") {
      arrEndfix[iTDir]=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"${test%%.root}\"", arrPathFileIn[iTDir].Data()));
      arrTFileIn[iTDir] = TFile::Open(arrPathFileIn[iTDir], "read");
    }
  }
  dumpCompareSeperated(arrTFileIn[0], arrTFileIn[1], normalize, arrEndfix[0], arrEndfix[1], pathDir, funNameFile, optionPrint, indexOrder);
  for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
    arrTFileIn[iTDir]->Close();
  }
}

void dumpCompareSeperated(TDirectory *tdir1, TDirectory *tdir2, Bool_t normalize, TString endfix1, TString endfix2, std::function<TString (TString, TString)> funPathDir, std::function<TString (TString, UInt_t)> funNameFile=nullptr, TString optionPrint="svg", Int_t indexOrder=1) {
  TString pathDir = funPathDir(endfix1, endfix2);
  dumpCompareSeperated(tdir1, tdir2, normalize, endfix1, endfix2, pathDir, funNameFile, optionPrint, indexOrder);
}

void dumpCompareSeperated(TString nameFileIn1, TString nameFileIn2, Bool_t normalize, TString endfix1, TString endfix2, std::function<TString (TString, TString)> funPathDir, std::function<TString (TString, UInt_t)> funNameFile=nullptr, TString optionPrint="svg", Int_t indexOrder=1) {
  TString pathDir = funPathDir(endfix1, endfix2);
  dumpCompareSeperated(nameFileIn1, nameFileIn2, normalize, endfix1, endfix2, pathDir, funNameFile, optionPrint, indexOrder);
}
