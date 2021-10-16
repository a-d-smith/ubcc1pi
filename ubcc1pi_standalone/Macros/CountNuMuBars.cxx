/**
 *  @file  ubcc1pi_standalone/Macros/CountNuMuBars.cxx
 *
 *  @brief The implementation file of the CountNuMuBars macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void CountNuMuBars(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    float nCC1PiEvents = 0.f;
    float nCC1PiNuMuBarEvents = 0.0f;

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

        // Work out if this event has a numubar
        const auto hasNuMuBar = AnalysisHelper::CountParticlesWithPdgCode(pEvent->truth.particles, -13, false);

        nCC1PiEvents += weight;
        if (hasNuMuBar) nCC1PiNuMuBarEvents += weight;
    }

    std::cout << "When finding CC1Pi events, useAbsPdg = " << config.global.useAbsPdg << std::endl;
    std::cout << "CC1Pi events: " << nCC1PiEvents << std::endl;
    std::cout << "  ... with numubar: " << nCC1PiNuMuBarEvents << std::endl;
}

} // namespace ubcc1pi_macros
