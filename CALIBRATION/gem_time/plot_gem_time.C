#include <TChain.h>
#include <TFile.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

void plot_gem_time() {
  const TString inputList = "filelist_optics1_17deg.txt";
  // Full dataset for actual analysis: "filelist_23367_23370.txt"
  const TString outputName = "gem_time.root";
  const char *spectrometers[] = {"H", "P"};
  constexpr int nSpectrometers = 2;
  constexpr int nModules = 2;
  constexpr int nSamples = 6;

  TChain chain("T");
  int nFiles = 0;
  if (inputList.EndsWith(".root")) {
    if (chain.Add(inputList.Data()) == 0) {
      std::cerr << "Error: unable to add input ROOT file " << inputList << std::endl;
      return;
    }
    ++nFiles;
  } else {
    std::ifstream fileList(inputList.Data());
    if (!fileList.is_open()) {
      std::cerr << "Error: unable to open file list " << inputList << std::endl;
      return;
    }

    std::string line;
    while (std::getline(fileList, line)) {
      const std::string::size_type first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos || line[first] == '#') {
        continue;
      }
      const std::string::size_type last = line.find_last_not_of(" \t\r");
      const std::string fileName = line.substr(first, last - first + 1);
      if (chain.Add(fileName.c_str()) == 0) {
        std::cerr << "Error: unable to add input file " << fileName << std::endl;
        return;
      }
      ++nFiles;
    }
  }
  if (nFiles == 0) {
    std::cerr << "Error: file list contains no ROOT files" << std::endl;
    return;
  }

  const Long64_t nEntries = chain.GetEntries();
  if (nEntries <= 0 || chain.LoadTree(0) < 0) {
    std::cerr << "Error: input chain contains no readable entries" << std::endl;
    return;
  }

  for (int iSpec = 0; iSpec < nSpectrometers; ++iSpec) {
    for (int iModule = 0; iModule < nModules; ++iModule) {
      const TString stripBranch =
          TString::Format("%s.gem.m%d.strip.istrip", spectrometers[iSpec], iModule);
      const TString isUBranch =
          TString::Format("%s.gem.m%d.strip.IsU", spectrometers[iSpec], iModule);
      const TString isVBranch =
          TString::Format("%s.gem.m%d.strip.IsV", spectrometers[iSpec], iModule);
      const TString sampleBranch =
          TString::Format("%s.gem.m%d.strip.ADCsamples", spectrometers[iSpec], iModule);
      const TString adcSumBranch =
          TString::Format("%s.gem.m%d.strip.ADCsum", spectrometers[iSpec], iModule);
      const TString tMeanBranch =
          TString::Format("%s.gem.m%d.strip.Tmean", spectrometers[iSpec], iModule);
      const TString maxSampleBranch =
          TString::Format("%s.gem.m%d.strip.isampmax", spectrometers[iSpec], iModule);
      const TString nStripBranch =
          TString::Format("Ndata.%s.gem.m%d.strip.istrip", spectrometers[iSpec], iModule);
      const TString nIsUBranch =
          TString::Format("Ndata.%s.gem.m%d.strip.IsU", spectrometers[iSpec], iModule);
      const TString nIsVBranch =
          TString::Format("Ndata.%s.gem.m%d.strip.IsV", spectrometers[iSpec], iModule);
      const TString nSampleBranch =
          TString::Format("Ndata.%s.gem.m%d.strip.ADCsamples", spectrometers[iSpec], iModule);
      const TString nAdcSumBranch =
          TString::Format("Ndata.%s.gem.m%d.strip.ADCsum", spectrometers[iSpec], iModule);
      const TString nTMeanBranch =
          TString::Format("Ndata.%s.gem.m%d.strip.Tmean", spectrometers[iSpec], iModule);
      const TString nMaxSampleBranch =
          TString::Format("Ndata.%s.gem.m%d.strip.isampmax", spectrometers[iSpec], iModule);
      if (!chain.GetBranch(stripBranch.Data()) || !chain.GetBranch(isUBranch.Data()) ||
          !chain.GetBranch(isVBranch.Data()) || !chain.GetBranch(sampleBranch.Data()) ||
          !chain.GetBranch(adcSumBranch.Data()) || !chain.GetBranch(tMeanBranch.Data()) ||
          !chain.GetBranch(maxSampleBranch.Data()) ||
          !chain.GetBranch(nStripBranch.Data()) || !chain.GetBranch(nIsUBranch.Data()) ||
          !chain.GetBranch(nIsVBranch.Data()) || !chain.GetBranch(nSampleBranch.Data()) ||
          !chain.GetBranch(nAdcSumBranch.Data()) || !chain.GetBranch(nTMeanBranch.Data()) ||
          !chain.GetBranch(nMaxSampleBranch.Data())) {
        std::cerr << "Error: required LADlib branch is missing for " << spectrometers[iSpec]
                  << " GEM module " << iModule << std::endl;
        return;
      }
    }
  }

  struct StripBranches {
    std::unique_ptr<TTreeReaderArray<Double_t>> strip;
    std::unique_ptr<TTreeReaderArray<Double_t>> isU;
    std::unique_ptr<TTreeReaderArray<Double_t>> isV;
    std::unique_ptr<TTreeReaderArray<Double_t>> samples;
    std::unique_ptr<TTreeReaderArray<Double_t>> adcSum;
    std::unique_ptr<TTreeReaderArray<Double_t>> tMean;
    std::unique_ptr<TTreeReaderArray<Double_t>> maxSample;
    std::unique_ptr<TTreeReaderValue<Int_t>> nStripEntries;
    std::unique_ptr<TTreeReaderValue<Int_t>> nIsUEntries;
    std::unique_ptr<TTreeReaderValue<Int_t>> nIsVEntries;
    std::unique_ptr<TTreeReaderValue<Int_t>> nSampleEntries;
    std::unique_ptr<TTreeReaderValue<Int_t>> nAdcSumEntries;
    std::unique_ptr<TTreeReaderValue<Int_t>> nTMeanEntries;
    std::unique_ptr<TTreeReaderValue<Int_t>> nMaxSampleEntries;
  };
  StripBranches branches[nSpectrometers][nModules];
  constexpr int nViews = 2;
  const char *viewNames[nViews] = {"U", "V"};

  TTreeReader reader(&chain);
  for (int iSpec = 0; iSpec < nSpectrometers; ++iSpec) {
    for (int iModule = 0; iModule < nModules; ++iModule) {
      branches[iSpec][iModule].strip.reset(new TTreeReaderArray<Double_t>(
          reader, TString::Format("%s.gem.m%d.strip.istrip", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].isU.reset(new TTreeReaderArray<Double_t>(
          reader, TString::Format("%s.gem.m%d.strip.IsU", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].isV.reset(new TTreeReaderArray<Double_t>(
          reader, TString::Format("%s.gem.m%d.strip.IsV", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].samples.reset(new TTreeReaderArray<Double_t>(
          reader, TString::Format("%s.gem.m%d.strip.ADCsamples", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].adcSum.reset(new TTreeReaderArray<Double_t>(
          reader, TString::Format("%s.gem.m%d.strip.ADCsum", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].tMean.reset(new TTreeReaderArray<Double_t>(
          reader, TString::Format("%s.gem.m%d.strip.Tmean", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].maxSample.reset(new TTreeReaderArray<Double_t>(
          reader, TString::Format("%s.gem.m%d.strip.isampmax", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].nStripEntries.reset(new TTreeReaderValue<Int_t>(
          reader, TString::Format("Ndata.%s.gem.m%d.strip.istrip", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].nIsUEntries.reset(new TTreeReaderValue<Int_t>(
          reader, TString::Format("Ndata.%s.gem.m%d.strip.IsU", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].nIsVEntries.reset(new TTreeReaderValue<Int_t>(
          reader, TString::Format("Ndata.%s.gem.m%d.strip.IsV", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].nSampleEntries.reset(new TTreeReaderValue<Int_t>(
          reader, TString::Format("Ndata.%s.gem.m%d.strip.ADCsamples", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].nAdcSumEntries.reset(new TTreeReaderValue<Int_t>(
          reader, TString::Format("Ndata.%s.gem.m%d.strip.ADCsum", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].nTMeanEntries.reset(new TTreeReaderValue<Int_t>(
          reader, TString::Format("Ndata.%s.gem.m%d.strip.Tmean", spectrometers[iSpec], iModule)));
      branches[iSpec][iModule].nMaxSampleEntries.reset(new TTreeReaderValue<Int_t>(
          reader, TString::Format("Ndata.%s.gem.m%d.strip.isampmax", spectrometers[iSpec], iModule)));
    }
  }

  std::unique_ptr<TH1D> hEventFraction[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hNegativeStripAdcSum[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hPositiveStripAdcSum[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hNegativeStripIndex[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hNegativeStripTMean[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hPositiveStripTMean[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hNegativeStripMaxSample[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hPositiveStripMaxSample[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hMostNegativeSample[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH1D> hNegativeSamplePattern[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH2D> hNegativeTMeanVsAdcSum[nSpectrometers][nModules][nViews];
  std::unique_ptr<TH2D> hPositiveTMeanVsAdcSum[nSpectrometers][nModules][nViews];
  Long64_t nFiredStrips[nSpectrometers][nModules][nViews] = {};
  Long64_t nNegativeSampleStrips[nSpectrometers][nModules][nViews] = {};
  Long64_t nEventsWithFiredStrips[nSpectrometers][nModules][nViews] = {};
  for (int iSpec = 0; iSpec < nSpectrometers; ++iSpec) {
    for (int iModule = 0; iModule < nModules; ++iModule) {
      for (int iView = 0; iView < nViews; ++iView) {
        const TString suffix =
            TString::Format("%s_m%d_%s", spectrometers[iSpec], iModule, viewNames[iView]);
        hEventFraction[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_event_fraction_negative_sample_%s", suffix.Data()),
                     TString::Format("%s per-event fraction of fired strips with a negative ADC sample;"
                                     "Negative-sample fired strips / all fired strips per event;Events",
                                     suffix.Data()),
                     101, -0.005, 1.005));
        hNegativeStripAdcSum[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_negative_sample_strip_adcsum_%s", suffix.Data()),
                     TString::Format("%s CM and ped subtracted ADCsum for fired strips with a negative ADC sample;"
                                     "ADCsum;Fired strips",
                                     suffix.Data()),
                     2000, 0, 20000));
        hPositiveStripAdcSum[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_all_positive_strip_adcsum_%s", suffix.Data()),
                     TString::Format("%s CM and ped subtracted ADCsum for fired strips with all six ADC samples positive;"
                                     "ADCsum;Fired strips",
                                     suffix.Data()),
                     2000, 0, 20000));
        const int nViewStrips = iView == 0 ? 12 * 128 : 24 * 128;
        hNegativeStripIndex[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_negative_sample_strip_index_%s", suffix.Data()),
                     TString::Format("%s strip index frequency for fired strips with any negative ADC sample;"
                                     "Strip index;Negative-sample strip occurrences",
                                     suffix.Data()),
                     nViewStrips, -0.5, nViewStrips - 0.5));
        hNegativeStripTMean[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_negative_sample_strip_tmean_%s", suffix.Data()),
                     TString::Format("%s Tmean for fired strips with any negative ADC sample;"
                                     "Tmean;Fired strips",
                                     suffix.Data()),
                     800, -400, 400));
        hPositiveStripTMean[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_all_positive_strip_tmean_%s", suffix.Data()),
                     TString::Format("%s Tmean for fired strips with all six ADC samples positive;"
                                     "Tmean;Fired strips",
                                     suffix.Data()),
                     800, -400, 400));
        hNegativeStripMaxSample[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_negative_sample_strip_maxsamp_%s", suffix.Data()),
                     TString::Format("%s max-sample index for fired strips with any negative ADC sample;"
                                     "Max-sample index;Fired strips",
                                     suffix.Data()),
                     nSamples, -0.5, nSamples - 0.5));
        hPositiveStripMaxSample[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_all_positive_strip_maxsamp_%s", suffix.Data()),
                     TString::Format("%s max-sample index for fired strips with all six ADC samples positive;"
                                     "Max-sample index;Fired strips",
                                     suffix.Data()),
                     nSamples, -0.5, nSamples - 0.5));
        hMostNegativeSample[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_most_negative_sample_adc_%s", suffix.Data()),
                     TString::Format("%s most negative corrected ADC sample per fired strip;"
                                     "Most negative ADC sample;Fired strips",
                                     suffix.Data()),
                     1000, -500, 0));
        hNegativeSamplePattern[iSpec][iModule][iView].reset(
            new TH1D(TString::Format("h_negative_sample_pattern_%s", suffix.Data()),
                     TString::Format("%s negative-sample pattern per fired strip;"
                                     "Negative-sample pattern;Fired strips",
                                     suffix.Data()),
                     5, 0.5, 5.5));
        hNegativeSamplePattern[iSpec][iModule][iView]->GetXaxis()->SetBinLabel(1, "Beginning only");
        hNegativeSamplePattern[iSpec][iModule][iView]->GetXaxis()->SetBinLabel(2, "End only");
        hNegativeSamplePattern[iSpec][iModule][iView]->GetXaxis()->SetBinLabel(
            3, "Both ends; positive middle");
        hNegativeSamplePattern[iSpec][iModule][iView]->GetXaxis()->SetBinLabel(4, "All six negative");
        hNegativeSamplePattern[iSpec][iModule][iView]->GetXaxis()->SetBinLabel(5, "Other/mixed");
        hNegativeTMeanVsAdcSum[iSpec][iModule][iView].reset(
            new TH2D(TString::Format("h_negative_sample_tmean_vs_adcsum_%s", suffix.Data()),
                     TString::Format("%s Tmean vs ADCsum for fired strips with any negative ADC sample;"
                                     "ADCsum;Tmean",
                                     suffix.Data()),
                     1000, 0, 20000, 800, -400, 400));
        hPositiveTMeanVsAdcSum[iSpec][iModule][iView].reset(
            new TH2D(TString::Format("h_all_positive_tmean_vs_adcsum_%s", suffix.Data()),
                     TString::Format("%s Tmean vs ADCsum for fired strips with all six ADC samples positive;"
                                     "ADCsum;Tmean",
                                     suffix.Data()),
                     1000, 0, 20000, 800, -400, 400));
        hEventFraction[iSpec][iModule][iView]->SetDirectory(nullptr);
        hNegativeStripAdcSum[iSpec][iModule][iView]->SetDirectory(nullptr);
        hPositiveStripAdcSum[iSpec][iModule][iView]->SetDirectory(nullptr);
        hNegativeStripIndex[iSpec][iModule][iView]->SetDirectory(nullptr);
        hNegativeStripTMean[iSpec][iModule][iView]->SetDirectory(nullptr);
        hPositiveStripTMean[iSpec][iModule][iView]->SetDirectory(nullptr);
        hNegativeStripMaxSample[iSpec][iModule][iView]->SetDirectory(nullptr);
        hPositiveStripMaxSample[iSpec][iModule][iView]->SetDirectory(nullptr);
        hMostNegativeSample[iSpec][iModule][iView]->SetDirectory(nullptr);
        hNegativeSamplePattern[iSpec][iModule][iView]->SetDirectory(nullptr);
        hNegativeTMeanVsAdcSum[iSpec][iModule][iView]->SetDirectory(nullptr);
        hPositiveTMeanVsAdcSum[iSpec][iModule][iView]->SetDirectory(nullptr);
      }
    }
  }

  Long64_t nProcessed = 0;
  while (reader.Next()) {
    Long64_t eventFiredStrips[nSpectrometers][nModules][nViews] = {};
    Long64_t eventNegativeSampleStrips[nSpectrometers][nModules][nViews] = {};

    for (int iSpec = 0; iSpec < nSpectrometers; ++iSpec) {
      for (int iModule = 0; iModule < nModules; ++iModule) {
        const StripBranches &module = branches[iSpec][iModule];
        const Int_t nStripEntries = **module.nStripEntries;
        const Int_t nIsUEntries = **module.nIsUEntries;
        const Int_t nIsVEntries = **module.nIsVEntries;
        const Int_t nSampleEntries = **module.nSampleEntries;
        const Int_t nAdcSumEntries = **module.nAdcSumEntries;
        const Int_t nTMeanEntries = **module.nTMeanEntries;
        const Int_t nMaxSampleEntries = **module.nMaxSampleEntries;
        const std::size_t nStrips = module.strip->GetSize();
        if (nStripEntries < 0 || nIsUEntries < 0 || nIsVEntries < 0 ||
            nSampleEntries < 0 || nAdcSumEntries < 0 || nTMeanEntries < 0 ||
            nMaxSampleEntries < 0 ||
            static_cast<std::size_t>(nStripEntries) != nStrips ||
            static_cast<std::size_t>(nIsUEntries) != nStrips ||
            static_cast<std::size_t>(nIsVEntries) != nStrips ||
            static_cast<std::size_t>(nSampleEntries) != module.samples->GetSize() ||
            static_cast<std::size_t>(nAdcSumEntries) != nStrips ||
            static_cast<std::size_t>(nTMeanEntries) != nStrips ||
            static_cast<std::size_t>(nMaxSampleEntries) != nStrips ||
            static_cast<std::size_t>(nSampleEntries) != nStrips * nSamples) {
          std::cerr << "Error: LADlib Ndata counts and strip-array sizes do not match at entry "
                    << nProcessed << " for " << spectrometers[iSpec] << " GEM module " << iModule
                    << " (Ndata istrip=" << nStripEntries << ", IsU=" << nIsUEntries
                    << ", IsV=" << nIsVEntries << ", ADCsamples=" << nSampleEntries
                    << ", ADCsum=" << nAdcSumEntries << ", Tmean=" << nTMeanEntries
                    << ", isampmax=" << nMaxSampleEntries << ")"
                    << std::endl;
          return;
        }

        for (std::size_t iStrip = 0; iStrip < nStrips; ++iStrip) {
          const Double_t isU = (*module.isU)[iStrip];
          const Double_t isV = (*module.isV)[iStrip];
          if ((isU != 0.0 && isU != 1.0) || (isV != 0.0 && isV != 1.0) ||
              isU + isV != 1.0) {
            std::cerr << "Error: strip is not classified as exactly one view at entry "
                      << nProcessed << " for " << spectrometers[iSpec] << " GEM module "
                      << iModule << ", strip-array index " << iStrip << std::endl;
            return;
          }
          const int iView = isU ? 0 : 1;
          ++eventFiredStrips[iSpec][iModule][iView];
          bool hasNegativeSample = false;
          bool allSamplesPositive = true;
          bool allSamplesNegative = true;
          bool hasZeroSample = false;
          Double_t mostNegativeSample = 0;
          Double_t samples[nSamples];
          for (int iSample = 0; iSample < nSamples; ++iSample) {
            const Double_t sample = (*module.samples)[iStrip * nSamples + iSample];
            samples[iSample] = sample;
            if (sample < 0) {
              hasNegativeSample = true;
              if (sample < mostNegativeSample) {
                mostNegativeSample = sample;
              }
            } else {
              allSamplesNegative = false;
            }
            if (sample <= 0) {
              allSamplesPositive = false;
            }
            if (sample == 0) {
              hasZeroSample = true;
            }
          }
          if (hasNegativeSample) {
            ++eventNegativeSampleStrips[iSpec][iModule][iView];
            const Double_t adcSum = (*module.adcSum)[iStrip];
            const Double_t tMean = (*module.tMean)[iStrip];
            hNegativeStripAdcSum[iSpec][iModule][iView]->Fill(adcSum);
            hNegativeStripIndex[iSpec][iModule][iView]->Fill((*module.strip)[iStrip]);
            hNegativeStripTMean[iSpec][iModule][iView]->Fill(tMean);
            hNegativeStripMaxSample[iSpec][iModule][iView]->Fill((*module.maxSample)[iStrip]);
            hMostNegativeSample[iSpec][iModule][iView]->Fill(mostNegativeSample);
            hNegativeTMeanVsAdcSum[iSpec][iModule][iView]->Fill(adcSum, tMean);
            int patternBin = 5;
            if (allSamplesNegative) {
              patternBin = 4;
            } else if (!hasZeroSample) {
              bool beginningOnly = samples[0] < 0;
              bool seenPositive = false;
              for (int iSample = 0; iSample < nSamples; ++iSample) {
                if (samples[iSample] > 0) {
                  seenPositive = true;
                } else if (seenPositive) {
                  beginningOnly = false;
                }
              }
              beginningOnly = beginningOnly && seenPositive;

              bool endOnly = samples[nSamples - 1] < 0;
              seenPositive = false;
              for (int iSample = nSamples - 1; iSample >= 0; --iSample) {
                if (samples[iSample] > 0) {
                  seenPositive = true;
                } else if (seenPositive) {
                  endOnly = false;
                }
              }
              endOnly = endOnly && seenPositive;

              bool bothEnds = samples[0] < 0 && samples[nSamples - 1] < 0;
              bool seenPositiveMiddle = false;
              bool negativeAfterMiddle = false;
              for (int iSample = 1; iSample < nSamples - 1; ++iSample) {
                if (samples[iSample] > 0) {
                  if (negativeAfterMiddle) {
                    bothEnds = false;
                  }
                  seenPositiveMiddle = true;
                } else if (seenPositiveMiddle) {
                  negativeAfterMiddle = true;
                }
              }
              bothEnds = bothEnds && seenPositiveMiddle;

              if (beginningOnly) {
                patternBin = 1;
              } else if (endOnly) {
                patternBin = 2;
              } else if (bothEnds) {
                patternBin = 3;
              }
            }
            hNegativeSamplePattern[iSpec][iModule][iView]->Fill(patternBin);
          } else if (allSamplesPositive) {
            const Double_t adcSum = (*module.adcSum)[iStrip];
            const Double_t tMean = (*module.tMean)[iStrip];
            hPositiveStripAdcSum[iSpec][iModule][iView]->Fill(adcSum);
            hPositiveStripTMean[iSpec][iModule][iView]->Fill(tMean);
            hPositiveStripMaxSample[iSpec][iModule][iView]->Fill((*module.maxSample)[iStrip]);
            hPositiveTMeanVsAdcSum[iSpec][iModule][iView]->Fill(adcSum, tMean);
          }
        }
      }
    }

    for (int iSpec = 0; iSpec < nSpectrometers; ++iSpec) {
      for (int iModule = 0; iModule < nModules; ++iModule) {
        for (int iView = 0; iView < nViews; ++iView) {
          if (eventFiredStrips[iSpec][iModule][iView] > 0) {
            hEventFraction[iSpec][iModule][iView]->Fill(
                static_cast<Double_t>(eventNegativeSampleStrips[iSpec][iModule][iView]) /
                eventFiredStrips[iSpec][iModule][iView]);
            ++nEventsWithFiredStrips[iSpec][iModule][iView];
            nFiredStrips[iSpec][iModule][iView] += eventFiredStrips[iSpec][iModule][iView];
            nNegativeSampleStrips[iSpec][iModule][iView] +=
                eventNegativeSampleStrips[iSpec][iModule][iView];
          }
        }
      }
    }
    ++nProcessed;
  }

  if (nProcessed != nEntries) {
    std::cerr << "Error: read " << nProcessed << " of " << nEntries << " input entries" << std::endl;
    return;
  }

  TFile outputFile(outputName, "RECREATE");
  if (outputFile.IsZombie()) {
    std::cerr << "Error: cannot create output ROOT file " << outputName << std::endl;
    return;
  }
  for (int iSpec = 0; iSpec < nSpectrometers; ++iSpec) {
    for (int iModule = 0; iModule < nModules; ++iModule) {
      for (int iView = 0; iView < nViews; ++iView) {
        hEventFraction[iSpec][iModule][iView]->Write();
        hNegativeStripAdcSum[iSpec][iModule][iView]->Write();
        hPositiveStripAdcSum[iSpec][iModule][iView]->Write();
        hNegativeStripIndex[iSpec][iModule][iView]->Write();
        hNegativeStripTMean[iSpec][iModule][iView]->Write();
        hPositiveStripTMean[iSpec][iModule][iView]->Write();
        hNegativeStripMaxSample[iSpec][iModule][iView]->Write();
        hPositiveStripMaxSample[iSpec][iModule][iView]->Write();
        hMostNegativeSample[iSpec][iModule][iView]->Write();
        hNegativeSamplePattern[iSpec][iModule][iView]->Write();
        hNegativeTMeanVsAdcSum[iSpec][iModule][iView]->Write();
        hPositiveTMeanVsAdcSum[iSpec][iModule][iView]->Write();
      }
    }
  }
  outputFile.Close();

  std::cout << "Processed " << nProcessed << " events." << std::endl;
  for (int iSpec = 0; iSpec < nSpectrometers; ++iSpec) {
    for (int iModule = 0; iModule < nModules; ++iModule) {
      for (int iView = 0; iView < nViews; ++iView) {
        std::cout << spectrometers[iSpec] << " GEM m" << iModule << " " << viewNames[iView] << ": "
                  << nNegativeSampleStrips[iSpec][iModule][iView] << " / "
                  << nFiredStrips[iSpec][iModule][iView] << " fired strips had a negative sample";
        if (nFiredStrips[iSpec][iModule][iView] > 0) {
          std::cout << " ("
                    << static_cast<Double_t>(nNegativeSampleStrips[iSpec][iModule][iView]) /
                           nFiredStrips[iSpec][iModule][iView]
                    << ")";
        }
        std::cout << "; " << nEventsWithFiredStrips[iSpec][iModule][iView]
                  << " events contributed to the fraction histogram." << std::endl;
      }
    }
  }
  std::cout << "Wrote per-spectrometer/module/view fraction and ADCsum histograms to " << outputName
            << "." << std::endl;
}
