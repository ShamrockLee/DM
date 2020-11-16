/// dumpCompare.C by Yueh-Shun Li (ShamrockLee)
// Based on the dumpComparePDF.C by Prof. Shin-Shan Eiko Yu (syuvivida)
// The function chi2NbinsCompare is borrowed directly
// and the function dumpCompare is a reimplementation
// of the function dumpComparePDF

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <THashList.h>

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
  Int_t nBins1 = h1->GetNbinsX();
  Int_t nBins2 = h2->GetNbinsX();
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
    if (nBins1) h1->Scale(1.0 / h1->Integral());
    if (nBins2) h2->Scale(1.0 / h2->Integral());
  }

  Double_t chi2 = 0.;
  Int_t nbins = 0;
  if (nBins1 && nBins2) {
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
  if (nBins1 && nBins2) {
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

template<UInt_t nTDir>
std::vector<UInt_t[nTDir]> getVArrIdxCommon(std::vector<TString> arrVName[nTDir], std::vector<UInt_t> arrVIdx[nTDir]) {
  std::vector<UInt_t[nTDir]> vArrIdxCommon;
  vArrIdxCommon.clear();
  UInt_t arrIdxIdx[nTDir];
  UInt_t arrSizeIdx[nTDir];
  for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
    arrIdxIdx[iTDir] = 0;
    arrSizeIdx[iTDir] = arrVIdx[iTDir].size();
    if (arrSizeIdx[iTDir] == 0) {
      return vArrIdxCommon;
    }
  };
  auto getName = [&arrVName, &arrVIdx](UInt_t iTDir, UInt_t ii)->TString {
    return arrVName[iTDir][arrVIdx[iTDir][ii]];
  };
  Bool_t finish = false;
  while (true) {
    TString nameMax = getName(0, arrIdxIdx[0]);
    for (UInt_t iTDir=1; iTDir<nTDir; iTDir++) {
      if (nameMax < getName(iTDir, arrIdxIdx[iTDir])) {
        nameMax = getName(iTDir, arrIdxIdx[iTDir]);
      }
    }
    Bool_t findCommon = false;
    while (!findCommon) {
      UInt_t iTDir = 0;
      for (iTDir=0; iTDir<nTDir; iTDir++) {
        while (getName(iTDir, arrIdxIdx[iTDir]) < nameMax) {
          arrIdxIdx[iTDir]++;
          if (arrIdxIdx[iTDir] >= arrSizeIdx[iTDir]) {
            return vArrIdxCommon; // finish
          }
        }
        if (getName(iTDir, arrIdxIdx[iTDir]) > nameMax) {
          nameMax = getName(iTDir, arrIdxIdx[iTDir]);
          break;
        }
      }
      findCommon = (iTDir == nTDir);
    }
    UInt_t arrIdxCommon[nTDir];
    for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
      arrIdxCommon[iTDir] = arrVIdx[iTDir][arrIdxIdx[iTDir]];
    }
    vArrIdxCommon.push_back(arrIdxCommon);
    for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
      arrIdxIdx[iTDir]++;
      if (arrIdxIdx[iTDir] >= arrSizeIdx[iTDir]) {
        return vArrIdxCommon; // finish
      }
    }
    // for (UInt_t iTdir=0; iTdir<nTDir; iTdir++) {
    //   if (arrIdxIdx[iTdir] >= arrSizeIdx[iTdir]) {
    //     finish = true;
    //     break;
    //   }
    // }
    // if (finish) break;
  }
}

void dumpCompareMerged(TDirectory *tdir1, TDirectory *tdir2, TString endfix1, TString endfix2, TString nameDir="fig", TString nameFile="", TString optionPrint="pdf", Bool_t normalize=true, Int_t indexOrder=1) {
  const UInt_t nTDir = 2;
  TDirectory *arrTDirIn[nTDir] = {tdir1, tdir2};
  std::vector<TString> arrVName[nTDir];
  std::vector<UInt_t> arrVIdx[nTDir];
  THashList *arrTListKeys[nTDir];
  for (UInt_t i=0; i<nTDir; i++) {
    refgetIdxHistOpenableRecent(arrTDirIn[i], arrVName[i], arrVIdx[i]);
    arrTListKeys[i] = (THashList *) arrTDirIn[i]->GetListOfKeys();
  }
  std::vector<UInt_t[nTDir]> vArrIdxCommon = getVArrIdxCommon<nTDir>(arrVName, arrVIdx);
  if (indexOrder != 0) {
    std::sort(vArrIdxCommon.begin(), vArrIdxCommon.end(), [&arrVName, &arrVIdx, &indexOrder](UInt_t arrIdxCommon[nTDir])->UInt_t {
      return arrIdxCommon[indexOrder-1];
    });
  }
  if (nameFile = "") {
    nameFile = endfix1 + "_" + endfix2;
  }
  TString filenamePrint = (nameDir == "" ? "" : nameDir + "/") + nameFile;
  TCanvas *c1 = new TCanvas;
  UInt_t nPage = 0;
  for (auto arrIdxCommon: vArrIdxCommon) {
    TH1 *arrHist[nTDir];
    for (UInt_t iTDir=0; iTDir<nTDir; iTDir++) {
      arrHist[iTDir] = (TH1 *) arrTListKeys[iTDir]->At(arrIdxCommon[iTDir]);
    }
    drawHistsCompared(arrHist[0], arrHist[1], endfix1, endfix2, normalize);
    c1->Print(nPage ? filenamePrint : (filenamePrint + "("), optionPrint);
  }
  c1->Print(filenamePrint + ")", optionPrint);
}