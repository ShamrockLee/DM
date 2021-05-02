#include <Rtypes.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TTree.h>
#include <TString.h>

#include "dumpCompare.C"

void plot_20210311(const Bool_t debug=false) {
  const TString namesLepton[] = {"Electron", "Muon", "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon", "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};

  const TString pathCondorPacks = "/home/shamrock/Documents/Project_HEP/ROOT_assignment_20191210/ctau-proper/lxplus_HTcondor";
  const TString nameCondorPack = "preselect";

  const TString pathOutImages = "/home/shamrock/Documents/Project_HEP/ROOT_assignment_20191210/out_images";
  const TString nameFolderToday = "output_" + nameCondorPack + "_20210311";

  auto funOpenFile = [pathCondorPacks, nameCondorPack](const TString nameDatagroup, const TString nameClusterID, const TString nameTT)->TFile*{
    return TFile::Open(pathCondorPacks + "/" + nameCondorPack + "/" + "output_" + nameDatagroup + "_" + nameTT + "_" + nameClusterID + "_hist.root");
  };

  const TString nameDatagroupSignalHeavy = "signal_Mx2-150_Mv-500_Mx1-1_ctau-1";
  const TString nameClusterIDSignalHeavy = "20210311";

  const TString nameDatagroupDYJets = "DYJets";
  const TString nameClusterIDDYJets = "1325718";

  const TString suffixImage = "svg";

  auto funRatioFromBoolHist = [](const TH1 *const hist)->Double_t {
    const Double_t valFalse =  static_cast<Double_t>(hist->GetBinContent(1));
    const Double_t valTrue = static_cast<Double_t>(hist->GetBinContent(2));
    if (valFalse == 0 && valTrue == 0) {
      return -1;
    }
    return valTrue / (valFalse + valTrue);
  };

  TFile *const fileSignalHeavyGenElectron                 =  funOpenFile(nameDatagroupSignalHeavy, nameClusterIDSignalHeavy, "Gen" + namesLepton[0]);
  TFile *const fileSignalHeavyNumCorrectElectron          =  funOpenFile(nameDatagroupSignalHeavy, nameClusterIDSignalHeavy, "NumCorrect" + namesLepton[0]);
  TFile *const fileSignalHeavyZMassCuttedElectron         =  funOpenFile(nameDatagroupSignalHeavy, nameClusterIDSignalHeavy, "ZMassCutted" + namesLepton[0]);
  TFile *const arrFileSignalHeavyPreselectedElectron[2]   = {funOpenFile(nameDatagroupSignalHeavy, nameClusterIDSignalHeavy, "PreselectedMatchingJet" + namesLepton[0]),
                                                             funOpenFile(nameDatagroupSignalHeavy, nameClusterIDSignalHeavy, "PreselectedMatchingFatJet" + namesLepton[0])};
  TFile *const arrFileSignalHeavyAllMatchedElectron[2]    = {funOpenFile(nameDatagroupSignalHeavy, nameClusterIDSignalHeavy, "AllMatchedJet" + namesLepton[0]),
                                                             funOpenFile(nameDatagroupSignalHeavy, nameClusterIDSignalHeavy, "AllMatchedFatJet" + namesLepton[0])};

  // TFile *const fileDYJetsGenElectron                      =  funOpenFile(nameDatagroupDYJets,      nameClusterIDDYJets,      "Gen" + namesLepton[0]);
  TFile *const fileDYJetsNumCorrectElectron               =  funOpenFile(nameDatagroupDYJets,      nameClusterIDDYJets,      "NumCorrect" + namesLepton[0]);
  TFile *const fileDYJetsZMassCuttedElectron              =  funOpenFile(nameDatagroupDYJets,      nameClusterIDDYJets,      "ZMassCutted" + namesLepton[0]);
  TFile *const arrFileDYJetsPreselectedEletron[2]         = {funOpenFile(nameDatagroupDYJets,      nameClusterIDDYJets,      "PreselectedMatchingJet" + namesLepton[0]),
                                                             funOpenFile(nameDatagroupDYJets,      nameClusterIDDYJets,      "PreselectedMatchingFatJet" + namesLepton[0])};

  gSystem->mkdir(pathOutImages + "/" + nameFolderToday);

  if (true) {
    TString nameLeafMain = "nElectron";
    std::vector<Bool_t> vIsHist = {};
    std::vector<TH1 *> vHist = {};
    std::vector<TString> vNameLeg = {};
    for (Byte_t isPassedPtEta = 0; isPassedPtEta < 2; isPassedPtEta++) {
      const TString nameLeafToGet = TString::Format("nElectron%s", isPassedPtEta ? "PassedPtEta" : "");
      vHist.push_back(static_cast<TH1*>(fileSignalHeavyGenElectron->Get("h" + nameLeafToGet)));
      vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
      vNameLeg.push_back("SignalHeavy_" + nameLeafToGet);
    }
    if (gPad) gPad->SetLogy(false);
    drawHistsMultiple(vIsHist, vHist, vNameLeg, "hNElectron", "nElectron", vHist[0]->GetXaxis()->GetTitle());
    gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
                    "Preselections_" + nameLeafMain + "." + suffixImage);
    gPad->SetLogy(true);
    gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
                    "Preselections_" + nameLeafMain + "_logy." + suffixImage);
  }

  if (false) {
    std::vector<TString> vName = {};
    std::vector<UInt_t> vIdx = {};
    refgetIdxHistOpenableRecent(fileSignalHeavyGenElectron, vName, vIdx);
    for (const auto &nameHist : vName) {
      const TString nameLeaf = nameHist(1, nameHist.Length());
      if (nameLeaf.BeginsWith("Electron_")) {
        if (debug) Info("plotting", "Plotting %s", nameLeaf.Data());
        const TString tstrUnderscoreAndSuffix = nameLeaf(8, nameLeaf.Length());
        std::vector<Bool_t> vIsHist = {};
        std::vector<TH1 *> vHist = {};
        std::vector<TString> vNameLeg = {};
        {
          const TString nameDatagroupShorter = "SignalHeavy";
          {
            const TString nameLeafToGet = nameLeaf;
            vHist.push_back(static_cast<TH1 *>(
                fileSignalHeavyGenElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_Gen_%s", nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
            vHist.push_back(static_cast<TH1 *>(
                fileSignalHeavyNumCorrectElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_NumCorrect_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
            vHist.push_back(static_cast<TH1 *>(
                fileSignalHeavyZMassCuttedElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_ZMassCutted_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
          }
          {
            const TString nameLeafToGet =
                namesLepton[0] + "NumCorrect" + tstrUnderscoreAndSuffix;
            vHist.push_back(static_cast<TH1 *>(
                fileSignalHeavyNumCorrectElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_NumCorrect_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
            vHist.push_back(static_cast<TH1 *>(
                fileSignalHeavyZMassCuttedElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_ZMassCutted_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
          }
        }
        {
          const TString nameDatagroupShorter = nameDatagroupDYJets;
          {
            const TString nameLeafToGet = nameLeaf;
            vHist.push_back(static_cast<TH1 *>(
                fileDYJetsZMassCuttedElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_NumCorrect_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
            vHist.push_back(static_cast<TH1 *>(
                fileDYJetsZMassCuttedElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_ZMassCutted_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
          }
        }
        if (gPad) gPad->SetLogy(false);
        drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeaf, true, false);
        gPad->SetTitle(nameLeaf);
        gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
                    "Preselections_" + nameLeaf + "." + suffixImage);
        gPad->SetLogy(true);
        gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
                    "Preselections_" + nameLeaf + "_logy." + suffixImage);
      }
      if (nameLeaf.BeginsWith("ElectronPair_")) {
        if (debug) Info("plotting", "Plotting %s", nameLeaf.Data());
        const TString tstrUnderscoreAndSuffix = nameLeaf(12, nameLeaf.Length());
        std::vector<Bool_t> vIsHist = {};
        std::vector<TH1 *> vHist = {};
        std::vector<TString> vNameLeg = {};
        {
          const TString nameDatagroupShorter = nameDatagroupSignalHeavy;
          {
            const TString nameLeafToGet = nameLeaf;
            vHist.push_back(static_cast<TH1 *>(
                fileSignalHeavyNumCorrectElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_NumCorrect_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
            vHist.push_back(static_cast<TH1 *>(
                fileSignalHeavyZMassCuttedElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_ZMassCutted_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
          }
        }
        {
          const TString nameDatagroupShorter = nameDatagroupDYJets;
          {
            const TString nameLeafToGet = nameLeaf;
            vHist.push_back(static_cast<TH1 *>(
                fileDYJetsNumCorrectElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_NumCorrect_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
            vHist.push_back(static_cast<TH1 *>(
                fileDYJetsZMassCuttedElectron->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(Form("%s_ZMassCutted_%s",
                                    nameDatagroupShorter.Data(),
                                    nameLeafToGet.Data()));
          }
        }
        if (gPad) gPad->SetLogy(false);
        drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeaf, nameLeaf,
                          vHist[0]->GetXaxis()->GetTitle());
        gPad->SetTitle(nameLeaf);
        gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
                    "Preselections_" + nameLeaf + "." + suffixImage);
        gPad->SetLogy(true);
        gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
                    "Preselections_" + nameLeaf + "_logy." + suffixImage);
      }
    }
  }

  for (Byte_t iAK = 0; iAK<2; iAK++) {
    {
      if (debug) Info("plotting", "Getting names from file %s", arrFileSignalHeavyPreselectedElectron[iAK]->GetName());
      std::vector<TString> vName = {};
      std::vector<UInt_t> vIdx = {};
      refgetIdxHistOpenableRecent(arrFileSignalHeavyPreselectedElectron[iAK], vName, vIdx);
      for (const auto& nameHist: vName) {
        const TString nameLeaf = nameHist(1, nameHist.Length());
        if (nameLeaf.BeginsWith(iAK ? "FatJet_" : "Jet_")) {
          if (debug) Info("plotting", "Plotting %s", nameLeaf.Data());
          const TString tstrUnderscoreAndSuffix = nameLeaf(iAK ? 6 : 3, nameLeaf.Length());
          std::vector<Bool_t> vIsHist = {};
          std::vector<TH1*> vHist = {};
          std::vector<TString> vNameLeg = {};
          TString tstrXTitle = "";
          {
            const TString nameDatagroupShorter = "SignalHeavy";
            {
              const TString nameLeafToGet = TString::Format("%sMatched%s", iAK ? "FatJet" : "Jet", tstrUnderscoreAndSuffix.Data());
              vHist.push_back(static_cast<TH1*>(arrFileSignalHeavyAllMatchedElectron[iAK]->Get("h" + nameLeafToGet)));
              if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
              vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
              vNameLeg.push_back(Form("AllMatched%s_%s", iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
            }
            // for (Int_t nFirst=1; nFirst <= (iAK ? 6 : 8); nFirst++) {
            {
              Int_t nFirst = (iAK ? 2 : 6);
              const TString nameLeafToGet = TString::Format("%sFirst%d%s", iAK ? "FatJet" : "Jet", nFirst, tstrUnderscoreAndSuffix.Data());
              vHist.push_back(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[iAK]->Get("h" + nameLeafToGet)));
              if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
              vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
              vNameLeg.push_back(Form("Preselected%s_%s", iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
            }

            {
              const TString& nameLeafToGet = nameLeaf;
              vHist.push_back(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[iAK]->Get("h" + nameLeafToGet)));
              if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
              vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
              vNameLeg.push_back(Form("Preselected%s_%s", iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
              if (vIsHist.back()) {
                tstrXTitle = vHist.back()->GetXaxis()->GetTitle();
              }
            }
          }
          if (gPad) gPad->SetLogy(false);
          drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeaf, nameLeaf, tstrXTitle);
          gPad->SetTitle(nameLeaf);
          gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedAndAllMatched_" + nameLeaf + "." + suffixImage);
          gPad->SetLogy(true);
          gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedAndAllMatched_" + nameLeaf + "_logy." + suffixImage);
          {
            const TString nameDatagroupShorter = nameDatagroupDYJets;
            
            // for (Int_t nFirst=1; nFirst <= (iAK ? 6 : 8); nFirst++) {
            {
              Int_t nFirst = (iAK ? 2 : 6);
              const TString nameLeafToGet = TString::Format("%sFirst%d%s", iAK ? "FatJet" : "Jet", nFirst, tstrUnderscoreAndSuffix.Data());
              vHist.push_back(static_cast<TH1*>(arrFileDYJetsPreselectedEletron[iAK]->Get("h" + nameLeafToGet)));
              if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
              vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
              vNameLeg.push_back(Form("%s_Preselected%s_%s", nameDatagroupShorter.Data(), iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
            }

            {
              const TString nameLeafToGet = nameLeaf;
              vHist.push_back(static_cast<TH1*>(arrFileDYJetsPreselectedEletron[iAK]->Get("h" + nameLeafToGet)));
              if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
              vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
              vNameLeg.push_back(Form("%s_Preselected%s_%s", nameDatagroupShorter.Data(), iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
            }
          }
          if (gPad) gPad->SetLogy(false);
          drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeaf, nameLeaf, tstrXTitle);
          gPad->SetTitle(nameLeaf);
          gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedAndAllMatched_" + nameLeaf + "vsDYJets" + "." + suffixImage);
          gPad->SetLogy(true);
          gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedAndAllMatched_" + nameLeaf + "vsDYJets" + "_logy." + suffixImage);
        }
        if (nameLeaf.EqualTo("SoftActivityJetHT")) {
          for (auto cstrNumber: {"", "2", "5", "10"}) {
            TString nameLeafToGet = nameLeaf + cstrNumber;
            std::vector<Bool_t> vIsHist = {};
            std::vector<TH1 *> vHist = {};
            std::vector<TString> vNameLeg = {};
            vHist.push_back(static_cast<TH1 *>(arrFileSignalHeavyPreselectedElectron[iAK]->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(TString::Format("SignalHeavy_Electron_%s",
                                               iAK ? "FatJet" : "Jet"));
            vHist.push_back(static_cast<TH1 *>(
                arrFileDYJetsPreselectedEletron[iAK]->Get("h" + nameLeafToGet)));
            vIsHist.push_back(vHist.back() != nullptr &&
                              !vHist.back()->IsZombie());
            vNameLeg.push_back(TString::Format("DYJets_Electron_%s",
                                               iAK ? "FatJet" : "Jet"));
            if (gPad) gPad->SetLogy(false);
            drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeafToGet,
                              nameLeafToGet + " (Preselected)",
                              vHist.back()->GetXaxis()->GetTitle());
            gPad->SetTitle(nameLeaf);
            gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
                        "PreselectedJetElectron_" + nameLeaf + "vsDYJets" +
                        "." + suffixImage);
            gPad->SetLogy(false);
            gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
                        "PreselectedJetElectron_" + nameLeaf + "vsDYJets" +
                        "_logy." + suffixImage);
          }
        }
        if (nameLeaf.Contains("MET_")) {
          std::vector<Bool_t> vIsHist = {};
          std::vector<TH1*> vHist = {};
          std::vector<TString> vNameLeg = {};
          vHist.push_back(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[iAK]->Get(nameHist)));
          vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
          vNameLeg.push_back(TString::Format("SignalHeavy_Electron_%s", iAK ? "FatJet" : "Jet"));
          vHist.push_back(static_cast<TH1*>(arrFileDYJetsPreselectedEletron[iAK]->Get(nameHist)));
          vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
          vNameLeg.push_back(TString::Format("DYJets_Electron_%s", iAK ? "FatJet" : "Jet"));
          if (gPad) gPad->SetLogy(false);
          drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeaf, nameLeaf + " (Preselected)", TString(vHist[0]->GetXaxis()->GetTitle()) + " (Pt > #)");
          gPad->SetTitle(nameLeaf);
          gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedJetElectron_" + nameLeaf + "vsDYJets" + "." + suffixImage);
          gPad->SetLogy(false);
          gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedJetElectron_" + nameLeaf + "vsDYJets" + "_logy." + suffixImage);
        }
      }
    }
  }

  if (true) {
    std::cout
    << "| " << "Cuts" << " | " << nameDatagroupSignalHeavy + "(" + namesLeptonNota[0] + "-pair)"<< " |\n"
    << "|-" << "|-" << "|\n"
    << "| " << "At lease 2 electrons in each event" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyGenElectron->Get("hhave2MoreElectron"))) << " |\n"
    << "| " << "Exactly 2 electrons in each event" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyGenElectron->Get("hhave2Electron"))) << " |\n"
    << "| " << "At lease 2 electrons passing Pt-Eta cut" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyGenElectron->Get("hhave2MoreElectronPassedPtEta"))) << " |\n"
    << "| " << "Exactly 2 electrons passing Pt-Eta cut" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyGenElectron->Get("hhave2ElectronPassedPtEta"))) << " |\n"
    << "| " << "Exactly 2 loose e (Pt>20GeV, Abs(Eta)<4.5)" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyGenElectron->Get("his" + namesLepton[0] + "NumCorrect"))) << " |\n"
    << "| " << "Z mass cut on e-pair +-20GeV (continued)" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyNumCorrectElectron->Get("his" + namesLepton[0] + "PairPassingZMassCut"))) << " |\n"
    // << "| " << " | " << " |\n"
    // << "| " << "without prior cuts:" << " |\n"
    // << "| " << "All GEN-level d's Pt>30GeV && Abs(Eta)<3" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyGenElectron->Get("hareAllGenDMatchingPassed"))) << " |\n"
    // << "| " << "All GEN-level ddbar's Pt>200GeV && Abs(Eta)<4.5" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyGenElectron->Get("hareAllGenDPairMatchingPassed"))) << " |\n"
    // << "| " << "after 2-e and Z mass cuts:" << " |\n"
    // << "| " << "All GEN-level d's Pt>30GeV && Abs(Eta)<3" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyZMassCuttedElectron->Get("hareAllGenDMatchingPassed"))) << " |\n"
    // << "| " << "All GEN-level ddbar's Pt>200GeV && Abs(Eta)<4.5" << " | " << funRatioFromBoolHist(static_cast<TH1*>(fileSignalHeavyZMassCuttedElectron->Get("hareAllGenDPairMatchingPassed"))) << " |\n"
    // << "| " << " | " << " |\n"
    // << "| " << "after preselections:" << " |\n"
    // << "| " << "after the above preselections:" << " |\n"
    // << "| " << "At least two d-AK4 match" << " | " << funRatioFromBoolHist(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[0]->Get("hhave2MoreJetPassed"))) << " |\n"
    // << "| " << "At least two ddbar-AK8 match" << " | " << funRatioFromBoolHist(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[1]->Get("hhave2MoreFatJetPassed"))) << " |\n"
    // << "| " << "All d has an AK4 match" << " | " << funRatioFromBoolHist(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[0]->Get("hareAllGenDMatched"))) << " |\n"
    // << "| " << "All ddbar has an AK8 match" << " | " << funRatioFromBoolHist(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[1]->Get("hareAllGenDPairMatched"))) << " |\n"
    << std::flush;
  }
}


int main(int argc, char* argv[]) {
  plot_20210311();
}
