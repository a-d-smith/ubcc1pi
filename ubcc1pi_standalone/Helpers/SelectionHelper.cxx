/**
 *  @file  ubcc1pi_standalone/Helpers/SelectionHelper.cxx
 *
 *  @brief The implementation file for the selection helper class
 */

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

namespace ubcc1pi
{

SelectionHelper::EventSelection::Cut::Cut() :
    m_name(""),
    m_hasValue(false),
    m_value(-std::numeric_limits<float>::max())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::Cut::Cut(const std::string &name) :
    m_name(name),
    m_hasValue(false),
    m_value(-std::numeric_limits<float>::max())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::Cut::Cut(const std::string &name, const float value) :
    m_name(name),
    m_hasValue(true),
    m_value(value)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string SelectionHelper::EventSelection::Cut::GetName() const
{
    if (m_name.empty())
        throw std::logic_error("Cut::GetName - Name has not been set");

    return m_name;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::Cut::HasValue() const
{
    return m_hasValue;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::Cut::GetValue() const
{
    if (!m_hasValue)
        throw std::logic_error("Cut::GetValue - Can't get value of cut: \"" + m_name + "\" - it has no value");

    return m_value;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::Cut::SetValue(const float value)
{
    if (!m_hasValue)
        throw std::logic_error("Cut::SetValue - Can't set value of cut: \"" + m_name + "\" - it has no value");

    m_value = value;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::CutTracker::CutTracker(const std::vector<Cut> &cuts, const unsigned int nRecoParticles) :
    m_cuts(cuts),
    m_assignedPdgs(nRecoParticles, 0) // ATTN here we use a PDG code of zero to mean "unassigned"
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::CutTracker::GetCutValue(const std::string &name) const
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("CutTracker::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->GetValue();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> SelectionHelper::EventSelection::CutTracker::GetCutsPassed() const
{
    return m_cutsPassed;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> SelectionHelper::EventSelection::CutTracker::GetAssignedPDGCodes() const
{
    return m_assignedPdgs;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::CutTracker::MarkCutAsPassed(const std::string &name)
{
    // Count the number of cuts that have already passed
    const auto nCutsPassed = m_cutsPassed.size();

    // Check we haven't passed all cuts already
    if (nCutsPassed == m_cuts.size())
        throw std::logic_error("CutTracker::MarkCutAsPassed - All cuts are already marked as passed");

    // Check that the cut we are marking as "passed" matches up with the next cut in the vector
    if (name != m_cuts.at(nCutsPassed).GetName())
        throw std::logic_error("CutTracker::MarkCutAsPassed - The cut \"" + name + "\" is not the next cut in the sequence");

    m_cutsPassed.push_back(name);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::CutTracker::AssignPDGCode(const unsigned int recoParticleIndex, const int pdg)
{
    if (recoParticleIndex >= m_assignedPdgs.size())
        throw std::out_of_range("CutTracker::AssignPDGCode - The input recoParticleIndex is out of range");

    // Set the PDG code
    m_assignedPdgs.at(recoParticleIndex) = pdg;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::EventSelection(const std::vector<Cut> &cuts, BDTMap &bdtMap, const SelectionLogic &logic) :
    m_cuts(cuts),
    m_bdtMap(bdtMap),
    m_logic(logic)
{
    // Make sure none of the input cut names are not repeated
    for (const auto &cut : m_cuts)
    {
        const auto cutName = cut.GetName();

        if (std::count_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) { return x.GetName() == cutName; }) != 1)
            throw std::invalid_argument("EventSelection::EventSelection - Repeated cut name: \"" + cutName + "\"");
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> SelectionHelper::EventSelection::GetCuts() const
{
    std::vector<std::string> cutNames;
    for (const auto &cut : m_cuts)
        cutNames.push_back(cut.GetName());

    return cutNames;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::CutHasValue(const std::string &name) const
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("EventSelection::CutHasValue - Unknown cut: \"" + name + "\"");

    return iter->HasValue();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::GetCutValue(const std::string &name) const
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("EventSelection::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->GetValue();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::SetCutValue(const std::string &name, const float value)
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("EventSelection::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->SetValue(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::SelectionResult SelectionHelper::EventSelection::Execute(const std::shared_ptr<Event> &pEvent)
{
    // Setup the cut tracker
    CutTracker cutTracker(m_cuts, pEvent->reco.particles.size());

    // Pass the event, BDTs and cut tracker through to the selection logic and run the selection
    const auto passed = m_logic(pEvent, m_bdtMap, cutTracker);

    // Return the selection result
    return {
        passed,
        cutTracker.GetCutsPassed(),
        cutTracker.GetAssignedPDGCodes()
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCCInclusiveSelection()
{
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        {"passesCCInclusive"}
    };

    // We don't need any BDTs so just use an empty map
    EventSelection::BDTMap bdtMap;

    // Define the actual selection logic
    const auto logic = [](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        // We must pass the CC inclusive selection
        if (!pEvent->reco.passesCCInclusive())
            return false;

        // Mark the cut "passesCCInclusive" as passed
        cutTracker.MarkCutAsPassed("passesCCInclusive");

        // Identify the muon candidate
        const auto &recoParticles = pEvent->reco.particles;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (particle.isCCInclusiveMuonCandidate())
            {
                // Assign the muon candidate a PDG code of 13
                cutTracker.AssignPDGCode(i, 13);
            }
        }

        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetDefaultSelection()
{
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        {"passesCCInclusive"},
        {"min2Tracks"},
        {"max1Uncontained"},
        {"2NonProtons", -0.06f},
        {"pionNotInGap"},
        {"muonNotInGap"},
        {"openingAngle", 2.65f},
        {"topologicalScore", 0.67f},
        {"startNearVertex", 9.5f},
        {"likelyGoldenPion", -0.03f}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)},
        {"goldenPion", std::make_shared<BDTHelper::BDT>("goldenPion", BDTHelper::GoldenPionBDTFeatureNames)}
    };

    // Define the actual selection logic
    const auto logic = [](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        // ----------------------------------------------------------------------------------
        // passesCCInclusive
        // ----------------------------------------------------------------------------------

        // Insist the event passes the CC inclusive preselection
        if (!pEvent->reco.passesCCInclusive())
            return false;

        // Mark the cut "passesCCInclusive" as passed
        cutTracker.MarkCutAsPassed("passesCCInclusive");

        // ----------------------------------------------------------------------------------
        // min2Tracks
        // ----------------------------------------------------------------------------------

        // Count the particles with a track fit and check if they are contained
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;

        const auto &recoParticles = pEvent->reco.particles;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            nTrackParticles++;

            if (!AnalysisHelper::IsContained(particle))
                nUncontainedParticles++;
        }

        // Insist that we have at least two tracks
        if (nTrackParticles < 2)
            return false;

        // Mark the cut "min2Tracks" as passed
        cutTracker.MarkCutAsPassed("min2Tracks");

        // ----------------------------------------------------------------------------------
        // max1Uncontained
        // ----------------------------------------------------------------------------------

        // Insist that at most one particle is uncontained
        if (nUncontainedParticles > 1)
            return false;

        // Mark the cut "max1Uncontained" as passed
        cutTracker.MarkCutAsPassed("max1Uncontained");

        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // ----------------------------------------------------------------------------------
        // 2NonProtons
        // ----------------------------------------------------------------------------------

        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThreshold = cutTracker.GetCutValue("2NonProtons");

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");

        // Keep track of the number of protons and pions we have identifies
        unsigned int nProtons = 0;
        std::vector<unsigned int> pionIndices;

        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
                continue;

            // Assume particles without a track fit are just small protons
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("DefaultSelection - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            if (!hasFeatures)
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThreshold)
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // If we've got to this point then we haven't identified the particle as a muon or a proton.
            // Instead let's identify the particle as a pion
            pionIndices.push_back(i);
            cutTracker.AssignPDGCode(i, 211);
        }

        // Sanity check that we have identified every particle
        const auto nMuons = 1u;
        const auto nPions = pionIndices.size();
        if (nProtons + nPions + nMuons != recoParticles.size())
            throw std::logic_error("DefaultSelection - Identified the wrong number of particles. This shouldn't happen!");

        // Insist that we exacly one pion (i.e we have have 2 non-protons)
        if (nPions != 1)
            return false;

        // Mark the cut "2NonProtons" as passed
        cutTracker.MarkCutAsPassed("2NonProtons");

        // ----------------------------------------------------------------------------------
        // pionNotInGap
        // ----------------------------------------------------------------------------------

        // Get the pion reco particle
        const auto pionIndex = pionIndices.front();
        const auto &pion = recoParticles.at(pionIndex);

        // Sanity check that our muon and pion are not the same particle
        if (muonIndex == pionIndex)
            throw std::logic_error("DefaultSelection - The muon and the pion candidates are the same particle. This shouldn't happen!");

        // Insist that the pion has at least one hit in each view (i.e. not in a gap)
        if (pion.nHitsU() == 0 || pion.nHitsV() == 0 || pion.nHitsW() == 0)
            return false;

        // Mark the cut "pionNotInGap" as passed
        cutTracker.MarkCutAsPassed("pionNotInGap");

        // ----------------------------------------------------------------------------------
        // muonNotInGap
        // ----------------------------------------------------------------------------------

        // Get the muon reco particle
        const auto &muon = recoParticles.at(muonIndex);

        // Insist that the muon has at least one hit in each view (i.e. not in a gap)
        if (muon.nHitsU() == 0 || muon.nHitsV() == 0 || muon.nHitsW() == 0)
            return false;

        // Mark the cut "muonNotInGap" as passed
        cutTracker.MarkCutAsPassed("muonNotInGap");

        // ----------------------------------------------------------------------------------
        // openingAngle
        // ----------------------------------------------------------------------------------

        // Get the opening angle cut value
        const auto maxOpeningAngle = cutTracker.GetCutValue("openingAngle");

        // Get the opening angle between the muon and pion
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(pionDir);

        // Insist that the opening angle isn't too wide
        if (openingAngle >= maxOpeningAngle)
            return false;

        // Mark the cut "openingAngle" as passed
        cutTracker.MarkCutAsPassed("openingAngle");

        // ----------------------------------------------------------------------------------
        // topologicalScore
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScore = cutTracker.GetCutValue("topologicalScore");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScore)
            return false;

        // Mark the cut "topologicalScore" as passed
        cutTracker.MarkCutAsPassed("topologicalScore");

        // ----------------------------------------------------------------------------------
        // startNearVertex
        // ----------------------------------------------------------------------------------

        // Get the start near vertex cut value
        const auto maxVertexDist = cutTracker.GetCutValue("startNearVertex");
        const auto maxVertexDist2 = maxVertexDist*maxVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
            const auto vertexDist2 = (start - recoVertex).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxVertexDist2)
                return false;
        }

        // Mark the cut "startNearVertex" as passed
        cutTracker.MarkCutAsPassed("startNearVertex");

        // ----------------------------------------------------------------------------------
        // likelyGoldenPion
        // ----------------------------------------------------------------------------------

        // Get the likely golden pion cut value
        const auto goldenPionBDTThreshold = cutTracker.GetCutValue("likelyGoldenPion");

        // Get the golden pion BDT
        auto &pGoldenPionBDT = bdtMap.at("goldenPion");

        // Get the features of the pion
        std::vector<float> features;
        if (!BDTHelper::GetBDTFeatures(pion, BDTHelper::GoldenPionBDTFeatureNames, features))
            throw std::logic_error("DefaultSelection - Can't get golden pion BDT features of pion candidate");

        // Insist that the BDT response is greater than the cut value to identify the pion as a golden pion
        const auto bdtResponse = pGoldenPionBDT->GetResponse(features);
        if (bdtResponse <= goldenPionBDTThreshold)
            return false;

        // Mark the cut "likelyGoldenPion" as passed
        cutTracker.MarkCutAsPassed("likelyGoldenPion");

        // We passed all cuts!
        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::IsCutPassed(const std::vector<string> &cutsPassed, const std::string &cut)
{
    return std::find(cutsPassed.begin(), cutsPassed.end(), cut) != cutsPassed.end();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int SelectionHelper::GetMuonCandidateIndex(const std::vector<Event::Reco::Particle> &particles, const std::vector<std::string> &featureNames, BDTHelper::BDT &muonBDT)
{
    bool foundCCInclusiveMuon = false;

    unsigned int ccInclusiveMuonIndex = std::numeric_limits<unsigned int>::max();
    unsigned int muonIndex = std::numeric_limits<unsigned int>::max();

    unsigned int nUncontainedParticles = 0;

    for (unsigned int index = 0; index < particles.size(); ++index)
    {
        const auto &particle = particles.at(index);

        // Check if this particle is the CC inclusive muon candidate
        if (particle.isCCInclusiveMuonCandidate())
        {
            if (foundCCInclusiveMuon)
                throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found multiple CC inclusive muon candidates");

            foundCCInclusiveMuon = true;
            ccInclusiveMuonIndex = index;
        }

        if (!AnalysisHelper::HasTrackFit(particle))
            continue;

        // For now make the last escaping particle the muon candidate, if there are multiple we will fall back on the CC inclusive candidate
        if (!AnalysisHelper::IsContained(particle))
        {
            muonIndex = index;
            nUncontainedParticles++;
        }
    }

    if (!foundCCInclusiveMuon)
        throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found no CC inclusive muon candidate");

    if (nUncontainedParticles > 1)
        return ccInclusiveMuonIndex;

    if (nUncontainedParticles == 1)
        return muonIndex;

    // If we are here all particles are contained, choose the muon using the BDT
    float maxMuonBDTResponse = -std::numeric_limits<float>::max();
    bool foundMuon = false;

    for (unsigned int index = 0; index < particles.size(); ++index)
    {
        const auto &particle = particles.at(index);

        if (!AnalysisHelper::HasTrackFit(particle))
            continue;

        if (!AnalysisHelper::IsContained(particle))
            throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found escaping particle when not expecting to!");

        std::vector<float> features;
        const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, featureNames, features);

        if (!hasFeatures)
            continue;

        const auto muonBDTResponse = muonBDT.GetResponse(features);

        if (muonBDTResponse < maxMuonBDTResponse)
            continue;

        maxMuonBDTResponse = muonBDTResponse;
        muonIndex = index;
        foundMuon = true;
    }

    // If no muon can be found, then default to the CC inclusive candidate
    if (!foundMuon)
        return ccInclusiveMuonIndex;

    return muonIndex;
}

} // namespace ubcc1pi
