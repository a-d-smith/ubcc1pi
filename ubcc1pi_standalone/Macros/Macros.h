/**
 *  @file  ubcc1pi_standalone/Macros/Macros.h
 *
 *  @brief The header file for the macros
 */

#ifndef UBCC1PI_STANDALONE_MACROS_MACROS
#define UBCC1PI_STANDALONE_MACROS_MACROS

#include "ubcc1pi_standalone/Objects/Config.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

    /**
     *  @brief  Macro to print the configuration
     *
     *  @param  config the input configuration
     */
    void PrintConfig(const Config &config = Config());

    /**
     *  @brief  Count the POT for the overlays and dirt samples and print the results
     *
     *  @param  config the input configuration
     */
    void CountPOT(const Config &config = Config());

    /**
     *  @brief  Get the run-subrun list that's needed to run Zarco's POT / trigger / spill counting tool
     *
     *  @param  config the input configuration
     */
    void GetRunSubrunList(const Config &config = Config());

    /**
     *  @brief  Look at the truth information related to signal events
     *
     *  @param  config the input configuration
     */
    void TruthStudy(const Config &config = Config());

    /**
     *  @brief  Count the number of true CC1Pi events and count those that have a numubar
     *
     *  @param  config the input configuration
     */
    void CountNuMuBars(const Config &config = Config());

    /**
     *  @brief  Look at the truth information relating to secondary interactions
     *
     *  @param  config the input configuration
     */
    void SecondaryInteractionsStudy(const Config &config = Config());

    /**
     *  @brief  Get the range->momentum curve fit parameters
     *
     *  @param  config the input configuration
     */
    void FitRangeCurves(const Config &config = Config());

    /**
     *  @brief  Determine the momentum thresholds to apply
     *
     *  @param  config the input configuration
     */
    void MomentumThresholdsStudy(const Config &config = Config());

    /**
     *  @brief  Look at the accuracy of the muon PID from the CC inclusive as a function of the muon kinematics
     *
     *  @param  config the input configuration
     */
    void CCInclusiveMuonPIDStudy(const Config &config = Config());

    /**
     *  @brief  Demonstrate the impact of track direction on the PID, and how it can be improved
     *
     *  @param  config the input configuration
     */
    void MultiPlanePIDDemo(const Config &config = Config());

    /**
     *  @brief  Plot the variables that form the input to the particle ID BDTs after the CC inclusive selection
     *
     *  @param  config the input configuration
     */
    void PlotInputVariables(const Config &config = Config());

    /**
     *  @brief  Plot the variables that form the input to the particle ID BDTs after the CC inclusive selection broken down for pions
     *
     *  @param  config the input configuration
     */
    void PlotPionInputVariables(const Config &config = Config());

    /**
     *  @brief  Get some event IDs for pions with various values for wigglines
     *
     *  @param  config the input configuration
     */
    void FindWigglyPions(const Config &config = Config());

    /**
     *  @brief  Get the correlation plots between the input variables to the BDT
     *
     *  @param  config the input configuration
     */
    void GetCorrelationPlots(const Config &config = Config());

    /**
     *  @brief  Train the BDTs
     *
     *  @param  config the input configuration
     */
    void TrainBDTs(const Config &config = Config());

    /**
     *  @brief  Run the N-1 BDT study
     *
     *  @param  config the input configuration
     */
    void NMinusOneBDTStudy(const Config &config = Config());

    /**
     *  @brief  Make data/MC plots for one iteration of the N-1 BDT study
     *
     *  @param  config the input configuration
     */
    void NMinusOneBDTDataMCStudy(const Config &config = Config());

    /**
    *  @brief  Run multiple iterations of the N-1 BDT study following a hard-coded ordering
    *
    *  @param  config the input configuration
    */
    void NMinusOneBDTStudyFull(const Config &config = Config());

    /**
     *  @brief  Study the impact of the muon PID
     *
     *  @param  config the input configuration
     */
    void MuonPIDStudy(const Config &config = Config());

    /**
     *  @brief  Make the event selection table
     *
     *  @param  config the input configuration
     */
    void MakeEventSelectionTable(const Config &config = Config());

    /**
     *  @brief  Make the PID table
     *
     *  @param  config the input configuration
     */
    void MakeSelectedPIDTable(const Config &config = Config());

    /**
     *  @brief  Make the selection efficiency plots
     *
     *  @param  config the input configuration
     */
    void MakeEventSelectionEfficiencyPlots(const Config &config = Config());

    /**
     *  @brief  Plot the variables that are used by the selection at the point of usage
     *
     *  @param  config the input configuration
     */
    void PlotEventSelectionCuts(const Config &config = Config());

    /**
     *  @brief  Plot the reconstructed variables for the muon candidate
     *
     *  @param  config the input configuration
     */
    void PlotMuonRecoVariables(const Config &config = Config());

    /**
     *  @brief  Plot the recontructed variables
     *
     *  @param  config the input configuration
     */
    void PlotReconstructedVariables(const Config &config = Config());

    /**
    *  @brief  Dump the truth and reconstructed information about selected events to a ROOT file
    *
    *  @param  config the input configuration
    */
    void DumpSelectedEventInfo(const Config &config = Config());

    /**
    *  @brief  Plot the reconstructed neutrino energy and the hadronic invariant mass
    *
    *  @param  config the input configuration
    */
    void PlotModelDependentQuantities(const Config &config = Config());

    /**
     *  @brief  Make the plots to motivate the binning choice
     *
     *  @param  config the input configuration
     */
    void MakeBinningPlots(const Config &config = Config());

    /**
     *  @brief  Plot the flux distribution
     *
     *  @param  config the input configuration
     */
    void PlotFlux(const Config &config = Config());

    /**
     *  @brief  Plot the variations of the flux for each systematic parameter
     *
     *  @param  config the input configuration
     */
    void PlotFluxVariations(const Config &config = Config());

    /**
     *  @brief  Extract the cross-sections
     *
     *  @param  config the input configuration
     */
    void ExtractXSecs(const Config &config = Config());

    /**
     *  @brief  Print a summary of the uncertainties for the total cross-section
     *
     *  @param  config the input configuration
     */
    void PrintUncertaintiesSummary(const Config &config = Config());

    /**
    *  @brief  Produce the plots for the previously extracted cross-section data
    *
    *  @param  config the input configuration
    */
    void MakeXSecPlots(const Config &config = Config());

    /**
     *  @brief  Top-level macro that runs the full analysis chain
     *
     *  @param  config the input config
     */
    void RunFullAnalysis(const Config &config = Config());


} // namespace ubcc1pi_macros

#endif

