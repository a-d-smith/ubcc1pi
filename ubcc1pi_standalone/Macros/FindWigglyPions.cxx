/**
 *  @file  ubcc1pi_standalone/Macros/FindWigglyPions.cxx
 *
 *  @brief The implementation file of the FindWigglyPions macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void FindWigglyPions(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();

    // Setup a container for the pions we find
    //                       run, subRun, event, isGoldenPion, wiggliness
    std::vector< std::tuple< int, int,    int,   bool,         float      > > pionData;

    // Loop over the events
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Ge the event weight
        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);

        // Only consider true CC1pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        // Only consider events passing the CC1pi generic selection
        const auto &[isSelectedTotal, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
        const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);
        if (!passedGenericSelection)
            continue;

        // Get the reconstructed pion
        const auto pionIndex = AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211);
        const auto pion = pEvent->reco.particles.at(pionIndex);

        // Check if it's really a true pion
        unsigned int truthParticleIndex = std::numeric_limits<unsigned int>::max();
        bool hasTruthMatch = false;
        try
        {
            truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(pion, pEvent->truth.particles);
            hasTruthMatch = true;
        }
        catch (const std::exception &)
        {
        }

        if (!hasTruthMatch)
            continue;

        const auto truthParticle = pEvent->truth.particles.at(truthParticleIndex);
        const auto truePdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();
        if (truePdgCode != 211)
            continue;

        // Now we have a true CC1Pi event that's passed the selection for which the reconstructed pion matches to a true pion
        // Store the info
        const auto isGoldenPion = AnalysisHelper::IsGolden(truthParticle);
        pionData.emplace_back(pEvent->metadata.run(), pEvent->metadata.subRun(), pEvent->metadata.event(), isGoldenPion, pion.wiggliness());
    }

    // Sort the pion data
    std::sort(pionData.begin(), pionData.end(), [](const auto &a, const auto &b){

        // Unpack the tuples
        const auto &[runA, subRunA, eventA, isGoldenPionA, wigglinessA] = a;
        const auto &[runB, subRunB, eventB, isGoldenPionB, wigglinessB] = b;

        // Put golden pions first
        if (isGoldenPionA != isGoldenPionB)
            return isGoldenPionA;

        // Then sort by wiggliness
        if (std::abs(wigglinessA - wigglinessB) > std::numeric_limits<float>::epsilon())
            return wigglinessA < wigglinessB;

        // If two entries have identical wiggliness, then sort by run, subRun and event
        if (runA != runB)
            return runA < runB;

        if (subRunA != subRunB)
            return subRunA < subRunB;

        if (eventA != eventB)
            return eventA < eventB;

        // At this point the entries are identical, so just stick with whatever order they were initially filled
        return true;
    });

    // Print the results in a table
    FormattingHelper::Table table({"Golden", "Wiggliness", "", "Run", "SubRun", "Event"});
    for (const auto &[run, subRun, event, isGoldenPion, wiggliness] : pionData)
    {
        table.AddEmptyRow();
        table.SetEntry("Golden", isGoldenPion);
        table.SetEntry("Wiggliness", wiggliness);
        table.SetEntry("Run", run);
        table.SetEntry("SubRun", subRun);
        table.SetEntry("Event", event);
    }

    table.WriteToFile("findWigglyPions.md");
}

} // namespace ubcc1pi_macros
