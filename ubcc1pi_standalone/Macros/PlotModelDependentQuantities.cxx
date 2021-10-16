/**
 *  @file  ubcc1pi_standalone/Macros/PlotModelDependentQuantities.cxx
 *
 *  @brief The implementation file of the PlotModelDependentQuantities macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotModelDependentQuantities(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    // Setup the plots
    PlottingHelper::MultiPlot neutrinoEnergyPlot("Neutrino energy / GeV", "Events", 10u, 0.f, 2.f);
    PlottingHelper::MultiPlot invariantMassPlot("Hadronic invariant mass / GeV", "Events", 10u, 0.f, 2.f);

    TH2F *hRes_nuEnergyEstRecoEstTrue = new TH2F("hRes_nuEnergyEstRecoEstTrue", "", 10, 0.f, 2.f, 10, 0.f, 2.f);
    TH2F *hRes_nuEnergyEstRecoTrue = new TH2F("hRes_nuEnergyEstRecoTrue", "", 10, 0.f, 2.f, 10, 0.f, 2.f);
    TH2F *hRes_nuEnergyEstTrueTrue = new TH2F("hRes_nuEnergyEstTrueTrue", "", 10, 0.f, 2.f, 10, 0.f, 2.f);

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetDefaultSelection();

    // Loop over the events
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // Run the event selection and store which cuts are passed
            const auto &[isSelectedGolden, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto isSelectedGeneric = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            // Only use events that at least pass the golden selection (as we want the pion momentum)
            if (!isSelectedGolden)
                continue;

            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            // Get the reco analysis data
            const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, isSelectedGolden);

            // Estimate the neutrino energy using the muon-pion kinematics from reco and truth
            const float muonMass = 0.1056583755f;
            const float pionMass = 0.13957039f;
            const float protonMass = 0.93827208816f;
            const float neutronMass = 0.93956542052f;
            const float nucleonMass = 0.5f * (protonMass + neutronMass); // Assume we don't know which nucleon is involved, so average

            const float muonEnergyReco = std::pow(recoData.muonMomentum*recoData.muonMomentum + muonMass*muonMass, 0.5f);
            const float pionEnergyReco = std::pow(recoData.pionMomentum*recoData.pionMomentum + pionMass*pionMass, 0.5f);

            // This formula is from the MiniBooNE CC1pi paper (this assumes a CC-Res interaction on a stationary target with no FSI)
            const float nuEnergyEstimateReco = ( muonMass * muonMass
                                               + pionMass * pionMass
                                               - 2 * nucleonMass * (muonEnergyReco + pionEnergyReco)
                                               + 2 * muonEnergyReco * pionEnergyReco
                                               - 2 * recoData.muonMomentum * recoData.pionMomentum * std::cos(recoData.muonPionAngle)
                                               ) / 2 * (
                                                 muonEnergyReco + pionEnergyReco
                                               - recoData.muonMomentum * recoData.muonCosTheta
                                               - recoData.pionMomentum * recoData.pionCosTheta
                                               - nucleonMass
                                               );

            // Estimate the hadronic invariant mass
            const float invariantMassReco = std::pow(
                std::pow(nucleonMass + nuEnergyEstimateReco - muonEnergyReco, 2)
                - (nuEnergyEstimateReco*nuEnergyEstimateReco + recoData.muonMomentum*recoData.muonMomentum - 2*nuEnergyEstimateReco*recoData.muonMomentum*recoData.muonCosTheta)
                , 0.5f);

            neutrinoEnergyPlot.Fill(nuEnergyEstimateReco, plotStyle, weight);
            invariantMassPlot.Fill(invariantMassReco, plotStyle, weight);

            // Now only deal with true CC1Pi events
            if (sampleType != AnalysisHelper::Overlay || !AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
                continue;

            // Estimate the neutrino energy using truth info
            const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
            const float muonEnergyTruth = std::pow(truthData.muonMomentum*truthData.muonMomentum + muonMass*muonMass, 0.5f);
            const float pionEnergyTruth = std::pow(truthData.pionMomentum*truthData.pionMomentum + pionMass*pionMass, 0.5f);

            const float nuEnergyEstimateTruth = ( muonMass * muonMass
                                               + pionMass * pionMass
                                               - 2 * nucleonMass * (muonEnergyTruth + pionEnergyTruth)
                                               + 2 * muonEnergyTruth * pionEnergyTruth
                                               - 2 * truthData.muonMomentum * truthData.pionMomentum * std::cos(truthData.muonPionAngle)
                                               ) / 2 * (
                                                 muonEnergyTruth + pionEnergyTruth
                                               - truthData.muonMomentum * truthData.muonCosTheta
                                               - truthData.pionMomentum * truthData.pionCosTheta
                                               - nucleonMass
                                               );

            // Estimate the hadronic invariant mass using truth info
            const float invariantMassTruth = std::pow(
                std::pow(nucleonMass + nuEnergyEstimateTruth - muonEnergyTruth, 2)
                - (nuEnergyEstimateTruth*nuEnergyEstimateTruth + truthData.muonMomentum*truthData.muonMomentum - 2*nuEnergyEstimateTruth*truthData.muonMomentum*truthData.muonCosTheta)
                , 0.5f);

            // Now get the true neutrino energy from the simulation
            const auto nuEnergyTruth = pEvent->truth.nuEnergy();

            // Fill the 2D plots
            hRes_nuEnergyEstRecoEstTrue->Fill(nuEnergyEstimateReco, nuEnergyEstimateTruth, weight);
            hRes_nuEnergyEstRecoTrue->Fill(nuEnergyEstimateReco, nuEnergyTruth, weight);
            hRes_nuEnergyEstTrueTrue->Fill(nuEnergyEstimateTruth, nuEnergyTruth, weight);
        }
    }

    neutrinoEnergyPlot.SaveAsStacked("modelDepentent_neutrinoEnergy");
    invariantMassPlot.SaveAsStacked("modelDepentent_hadronicInvariantMass");

    auto pCanvas = PlottingHelper::GetCanvas();
    hRes_nuEnergyEstRecoEstTrue->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "modelDepentent_neutrinoEnergy_estReco-vs-estTrue");

    hRes_nuEnergyEstRecoTrue->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "modelDepentent_neutrinoEnergy_estReco-vs-true");

    hRes_nuEnergyEstTrueTrue->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "modelDepentent_neutrinoEnergy_estTrue-vs-true");
}

} // namespace ubcc1pi_macros
