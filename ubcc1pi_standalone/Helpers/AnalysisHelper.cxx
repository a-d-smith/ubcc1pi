/**
 *  @file  ubcc1pi_standalone/Helpers/AnalysisHelper.cxx
 *
 *  @brief The implementation of the analysis helper class
 */

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include "ubcc1pi_standalone/Helpers/GeometryHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <stdexcept>
#include <algorithm>
#include <numeric>

#include <TGraph.h>
#include <TFitResult.h>

namespace ubcc1pi
{

std::string AnalysisHelper::GetSampleTypeName(const SampleType &sampleType)
{
    switch (sampleType)
    {
        case DataBNB:
            return "Data BNB";
        case Overlay:
            return "Overlay";
        case DataEXT:
            return "Data EXT";
        case Dirt:
            return "Dirt";
        case DetectorVariation:
            return "Detector variation";
        default: break;
    }

    throw std::invalid_argument("AnalysisHelper::GetSampleTypeName - unknown sample type");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::EventCounter::EventCounter()
{
    // The "all" is treated differently when printing out
    m_tags.push_back("all");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetWeight(const std::string &tag, const SampleType &sampleType, const std::string &classification) const
{
    float weight = 0.f;
    if (this->GetWeight(tag, sampleType, classification, weight))
        return weight;

    return 0.f;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::EventCounter::GetWeight(const std::string &tag, const SampleType &sampleType, const std::string &classification, float &weight) const
{
    const auto tagIter = m_eventWeightMap.find(tag);
    if (tagIter == m_eventWeightMap.end())
        return false;

    const auto tagMap = tagIter->second;

    const auto sampleIter = tagMap.find(sampleType);
    if (sampleIter == tagMap.end())
        return false;

    const auto sampleMap = sampleIter->second;

    const auto classificationIter = sampleMap.find(classification);
    if (classificationIter == sampleMap.end())
        return false;

    weight = classificationIter->second;
    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetTotalMCWeight(const std::string &tag) const
{
    // Get the total weight for this tag
    float totalWeight = 0.f;
    for (const auto &type : AnalysisHelper::AllSampleTypes)
    {
        if (type == DataBNB)
            continue;

        for (const auto &classification : m_classifications)
            totalWeight += this->GetWeight(tag, type, classification);
    }

    return totalWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetBNBDataWeight(const std::string &tag) const
{
    float weight = 0.f;

    for (const auto &classification : m_classifications)
        weight += this->GetWeight(tag, DataBNB, classification);

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::EventCounter::IsSignalClassification(const std::string &classification, const std::string &searchQuery) const
{
    if (classification.empty())
        return false;

    // Signal classifications always begin with 'S'
    if (classification.at(0) != 'S')
        return false;

    if (searchQuery.empty())
        return true;

    // Insist that there is also the additional query
    return (classification.find(searchQuery) != std::string::npos);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> AnalysisHelper::EventCounter::GetSignalClassifications(const std::string &searchQuery) const
{
    std::vector<std::string> signalClassifications;

    for (const auto &classification : m_classifications)
    {
        if (!this->IsSignalClassification(classification, searchQuery))
            continue;

        signalClassifications.push_back(classification);
    }

    return signalClassifications;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetSignalWeight(const std::string &tag, const std::string &searchQuery) const
{
    float weight = 0.f;

    for (const auto &classification : m_classifications)
    {
        if (!this->IsSignalClassification(classification, searchQuery))
            continue;

        weight += this->GetWeight(tag, Overlay, classification);
    }

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetBackgroundWeight(const std::string &tag, const std::string &searchQuery) const
{
    float weight = 0.f;

    for (const auto &classification : m_classifications)
    {
        if (this->IsSignalClassification(classification, searchQuery))
            continue;

        for (const auto &type : AnalysisHelper::AllSampleTypes)
        {
            if (type == DataBNB)
                continue;

            weight += this->GetWeight(tag, type, classification);
        }
    }

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetSignalEfficiency(const std::string &tag, const std::string &searchQuery) const
{
    const auto allWeight = this->GetSignalWeight("all", searchQuery);

    if (allWeight <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    const auto weight = this->GetSignalWeight(tag, searchQuery);

    return weight / allWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetSignalPurity(const std::string &tag, const std::string &searchQuery) const
{
    if (tag == "all")
        throw std::invalid_argument("AnalysisHelper::EventCounter::GetSignalPurity - Can't get purity for tag \"all\", as it won't be available for off-beam data");

    const auto totalWeight = this->GetTotalMCWeight(tag);

    if (totalWeight <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    const auto weight = this->GetSignalWeight(tag, searchQuery);

    return weight / totalWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetEfficiency(const std::string &tag, const SampleType &sampleType, const std::string &classification) const
{
    if (sampleType != Overlay)
        throw std::invalid_argument("AnalysisHelper::EventCounter::GetEfficiency - Can't get efficiency of data!");

    const auto allWeight = this->GetWeight("all", sampleType, classification);

    if (allWeight <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    const auto weight = this->GetWeight(tag, sampleType, classification);

    return weight / allWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetPurity(const std::string &tag, const SampleType &sampleType, const std::string &classification) const
{
    if (sampleType == DataBNB)
        throw std::invalid_argument("AnalysisHelper::EventCounter::GetPurity - Can't get purity of BNB data!");

    if (tag == "all")
        throw std::invalid_argument("AnalysisHelper::EventCounter::GetPurity - Can't get purity for tag \"all\", as it won't be available for off-beam data");

    const auto totalWeight = this->GetTotalMCWeight(tag);

    if (totalWeight <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::epsilon();

    const auto weight = this->GetWeight(tag, sampleType, classification);

    return weight / totalWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::CountEvent(const std::string &tag, const SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const float weight)
{
    // Keep track of this tag if we haven't seen it before
    if (std::find(m_tags.begin(), m_tags.end(), tag) == m_tags.end())
        m_tags.push_back(tag);

    // Classify the event
    const auto useAbsPdg = true; // TODO make this configurable
    const auto countProtonsInclusively = true; // TODO make this configurable
    const auto classification = AnalysisHelper::GetClassificationString(pEvent, useAbsPdg, countProtonsInclusively);

    // Keep track of this classification if we haven't seen it before
    if (std::find(m_classifications.begin(), m_classifications.end(), classification) == m_classifications.end())
        m_classifications.push_back(classification);

    std::sort(m_classifications.begin(), m_classifications.end());

    // Get the mapping from the classification string to the weight, making it if it doesn't exist for this tag and sampleType
    auto &stringToWeightMap = m_eventWeightMap[tag][sampleType];
    auto weightIter = stringToWeightMap.find(classification);
    if (weightIter == stringToWeightMap.end())
    {
        stringToWeightMap.emplace(classification, weight);
    }
    else
    {
        weightIter->second += weight;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> AnalysisHelper::EventCounter::GetTags() const
{
    return m_tags;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::PrintBreakdownSummary(const std::string &outputFileName) const
{
    FormattingHelper::Table table({"Tag", "", "Signal", "Background", "Efficiency", "Purity", "E*P", "Golden fraction", "", "BNB Data", "Data/MC ratio"});
    for (const auto &tag : m_tags)
    {
        table.AddEmptyRow();
        table.SetEntry("Tag", tag);

        table.SetEntry("Signal", this->GetSignalWeight(tag));
        table.SetEntry("Background", this->GetBackgroundWeight(tag));

        const auto efficiency = this->GetSignalEfficiency(tag);
        table.SetEntry("Efficiency", efficiency);

        const auto isAll = (tag == "all");
        if (tag == "all")
        {
            table.SetEntry("Purity", "?");
            table.SetEntry("E*P", "?");
        }
        else
        {
            const auto purity = this->GetSignalPurity(tag);

            table.SetEntry("Purity", purity);
            table.SetEntry("E*P", efficiency * purity);
        }

        table.SetEntry("Golden fraction", this->GetSignalWeight(tag, "S G") / this->GetSignalWeight(tag));

        const auto bnbDataWeight = this->GetBNBDataWeight(tag);
        const auto totalMCWeight = this->GetTotalMCWeight(tag);

        table.SetEntry("BNB Data", bnbDataWeight);

        if (totalMCWeight <= std::numeric_limits<float>::epsilon())
        {
            table.SetEntry("Data/MC ratio", "?");
        }
        else
        {
            table.SetEntry("Data/MC ratio", bnbDataWeight / totalMCWeight);
        }
    }

    table.WriteToFile(outputFileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::PrintBreakdownDetails(const std::string &outputFileName, const unsigned int nEntries) const
{
    FormattingHelper::Table table({"Tag", "", "Sample", "Classification", "", "Weight", "Efficiency", "Purity"});

    for (const auto &tag : m_tags)
    {
        // Print a separator row
        if (tag != m_tags.front())
            table.AddEmptyRow();

        // The total MC weight we have and haven't printed
        float totalMCWeight = 0.f;
        float totalMCWeightPrinted = 0.f;

        for (const auto &sample : AnalysisHelper::AllSampleTypes)
        {
            const auto sampleName = AnalysisHelper::GetSampleTypeName(sample);
            std::vector<std::pair<std::string, float> > classificationWeightVector;

            for (const auto &classification : m_classifications)
            {
                float weight = 0.f;
                if (!this->GetWeight(tag, sample, classification, weight))
                    continue;

                classificationWeightVector.emplace_back(classification, weight);
            }

            std::sort(classificationWeightVector.begin(), classificationWeightVector.end(), [&](const auto &a, const auto &b){
                const bool isASignal = this->IsSignalClassification(a.first);
                const bool isBSignal = this->IsSignalClassification(b.first);

                // Put signal before background
                if (isASignal != isBSignal)
                    return isASignal;

                // Order in numerically, largest numbers first
                return a.second > b.second;
            });

            // Print the entries
            unsigned int entriesPrinted = 0;
            for (const auto &entry : classificationWeightVector)
            {
                const auto &classification = entry.first;
                const auto weight = entry.second;

                if (sample != AnalysisHelper::DataBNB)
                    totalMCWeight += weight;

                if (entriesPrinted == nEntries)
                    continue;

                if (sample != AnalysisHelper::DataBNB)
                    totalMCWeightPrinted += weight;

                table.AddEmptyRow();
                table.SetEntry("Tag", tag);
                table.SetEntry("Sample", sampleName);
                table.SetEntry("Classification", classification);
                table.SetEntry("Weight", weight);

                if (sample == Overlay)
                {
                    table.SetEntry("Efficiency", this->GetEfficiency(tag, sample, classification));
                }

                if (sample != DataBNB && tag != "all")
                {
                    table.SetEntry("Purity", this->GetPurity(tag, sample, classification));
                }

                entriesPrinted++;
            }
        }

        // Print the other catagory
        const auto mcWeightNotPrinted = totalMCWeight - totalMCWeightPrinted;
        table.AddEmptyRow();
        table.SetEntry("Tag", tag);
        table.SetEntry("Classification", "Other");
        table.SetEntry("Weight", mcWeightNotPrinted);

        if (totalMCWeight > std::numeric_limits<float>::epsilon())
        {
            const auto otherPurity = mcWeightNotPrinted / totalMCWeight;
            table.SetEntry("Purity", otherPurity);
        }

    }

    table.WriteToFile(outputFileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetNominalEventWeight(const std::shared_ptr<Event> &pEvent)
{
    float weight = 1.f;

    const auto &truth = pEvent->truth;

    // ATTN unfortunately when the event weights are applied a small fraction of events are given an infinite weight, I assume due to a
    // failure in the event weighting code - for now we just skip infinite weights
    if (truth.splineEventWeight.IsSet() && std::abs(truth.splineEventWeight()) < std::numeric_limits<float>::max())
        weight *= truth.splineEventWeight();

    if (truth.genieTuneEventWeight.IsSet() && std::abs(truth.genieTuneEventWeight()) < std::numeric_limits<float>::max())
        weight *= truth.genieTuneEventWeight();

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsTrueCCInclusive(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg)
{
    if (!pEvent->metadata.hasTruthInfo())
        throw std::invalid_argument("AnalysisHelper::IsTrueCCInclusive - Input event doesn't have truth information!");

    const auto truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return false;

    // Count the visible particles
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
    const auto nMu = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 13, useAbsPdg);

    // Insist the CCInclusive topology
    return (nMu == 1);

}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsTrueCC1Pi(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg)
{
    if (!pEvent->metadata.hasTruthInfo())
        throw std::invalid_argument("AnalysisHelper::IsTrueCC1Pi - Input event doesn't have truth information!");

    const auto truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return false;

    // Count the visible particles
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
    const auto nMu = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 13, useAbsPdg);
    const auto nProton = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 2212, useAbsPdg);
    const auto nPion = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 211, useAbsPdg);
    const auto nOther = visibleParticles.size() - (nMu + nProton + nPion);

    // Insist the CC1Pi topology
    return (nMu == 1 && nPion == 1 && nOther == 0);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::PassesVisibilityThreshold(const Event::Truth::Particle &particle)
{
    // First only consider specific particle types that can lead to hits
    const auto absPDG = std::abs(particle.pdgCode());

    return (absPDG == 11   || // Electron
            absPDG == 13   || // Muon
            absPDG == 2212 || // Proton
            absPDG == 22   || // Photon
            absPDG == 111  || // Pi0 (visible via decay products)
            absPDG == 211  || // Charged pion
            absPDG == 321);   // Charged kaon
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<Event::Truth::Particle> AnalysisHelper::SelectVisibleParticles(const std::vector<Event::Truth::Particle> &particles)
{
    std::vector<Event::Truth::Particle> visibleParticles;

    for (const auto &particle : particles)
    {
        if (AnalysisHelper::PassesVisibilityThreshold(particle))
            visibleParticles.push_back(particle);
    }

    return visibleParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::CountParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode, const bool useAbsPdg)
{
    unsigned int count = 0;

    for (const auto &particle : particles)
    {
        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();
        if (pdg == pdgCode)
            count++;
    }

    return count;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::CountGoldenParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode, const bool useAbsPdg)
{
    std::vector<Event::Truth::Particle> goldenParticles;
    for (const auto &particle : particles)
    {
        if (AnalysisHelper::IsGolden(particle))
            goldenParticles.push_back(particle);
    }

    return AnalysisHelper::CountParticlesWithPdgCode(goldenParticles, pdgCode, useAbsPdg);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::GetPdgCodeCountMap(const std::vector<Event::Truth::Particle> &particles, const bool useAbsPdg, std::vector<int> &foundPdgs, std::unordered_map<int, unsigned int> &pdgCodeCountMap)
{
    if (!foundPdgs.empty())
        throw std::invalid_argument("AnalysisHelper::GetPdgCodeCountMap - input foundPdgs vector isn't empty");

    if (!pdgCodeCountMap.empty())
        throw std::invalid_argument("AnalysisHelper::GetPdgCodeCountMap - input pdgCodeCountMap isn't empty");

    // Get the mapping
    for (const auto &particle : particles)
    {
        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();
        auto iter = pdgCodeCountMap.find(pdg);

        if (iter == pdgCodeCountMap.end())
        {
            // Add a new entry
            pdgCodeCountMap.emplace(pdg, 1);
            foundPdgs.push_back(pdg);
        }
        else
        {
            // Increment the counter
            iter->second++;
        }
    }

    // Sort the vector for reproducibility. Here we sort numerically in increasing order, but place particles and antiparticles next to each other
    std::sort(foundPdgs.begin(), foundPdgs.end(), [](const int &a, const int &b) {

        const auto aAbs = std::abs(a);
        const auto bAbs = std::abs(b);

        return (aAbs == bAbs) ? a > b : aAbs < bAbs;
    });
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string AnalysisHelper::GetTopologyString(const std::vector<Event::Truth::Particle> &particles, const bool useAbsPdg, const bool countProtonsInclusively)
{
    // Count the particles by PDG code
    std::vector<int> foundPdgs;
    std::unordered_map<int, unsigned int> pdgCodeCountMap;
    AnalysisHelper::GetPdgCodeCountMap(particles, useAbsPdg, foundPdgs, pdgCodeCountMap);

    unsigned int nOther = 0;
    unsigned int nProtons = 0;
    std::string topology = "";

    for (const auto &pdg : foundPdgs)
    {
        const auto count = pdgCodeCountMap.at(pdg);
        const auto countStr = std::to_string(count);

        switch (pdg)
        {
            case 11:
                topology += countStr + " e" + (useAbsPdg ? "" : "-") + "  ";
                break;
            case -11:
                topology += countStr + " e" + (useAbsPdg ? "" : "+") + "  ";
                break;
            case 13:
                topology += countStr + " Mu" + (useAbsPdg ? "" : "-") + "  ";
                break;
            case -13:
                topology += countStr + " Mu" + (useAbsPdg ? "" : "+") + "  ";
                break;
            case 22:
                topology += countStr + " Gamma  ";
                break;
            case 211:
                topology += countStr + " Pi" + (useAbsPdg ? "" : "+") + "  ";
                break;
            case -211:
                topology += countStr + " Pi" + (useAbsPdg ? "" : "-") + "  ";
                break;
            case 111:
                topology += countStr + " Pi0  ";
                break;
            case 321:
                topology += countStr + " K" + (useAbsPdg ? "" : "+") + "  ";
                break;
            case -321:
                topology += countStr + " K" + (useAbsPdg ? "" : "-") + "  ";
                break;
            case 2112:
                topology += countStr + " n  ";
                break;
            case 2212:
                nProtons = count;
                if (!countProtonsInclusively)
                    topology += countStr + " p  ";
                break;
            default:
                nOther++;
                break;
        }
    }

    if (countProtonsInclusively)
        topology += "X p  ";

    if (nOther != 0)
        topology += std::to_string(nOther) + " other  ";


    return topology;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string AnalysisHelper::GetClassificationString(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg, const bool countProtonsInclusively)
{
    // Check if data
    if (!pEvent->metadata.hasTruthInfo())
        return "D,  ";

    const auto truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return "NF, ";

    std::string classification = "";
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(truth.particles);

    // Signal or background
    const auto isTrueCC1Pi = AnalysisHelper::IsTrueCC1Pi(pEvent, useAbsPdg);
    classification += isTrueCC1Pi ? "S" : "B,  ";

    if (isTrueCC1Pi)
    {
        // Check if we have a golden pion
        const auto hasGoldenPion = (AnalysisHelper::CountGoldenParticlesWithPdgCode(visibleParticles, 211, useAbsPdg) != 0);

        if (hasGoldenPion)
            classification += " G,";
        else
            classification += ",  ";
    }

    // Add the topology classification
    classification += "  " + AnalysisHelper::GetTopologyString(visibleParticles, useAbsPdg, countProtonsInclusively);

    return classification;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsPointWithinMargins(const TVector3 &point, const float lowXMargin, const float highXMargin, const float lowYMargin, const float highYMargin, const float lowZMargin, const float highZMargin)
{
    return ((point.x() > GeometryHelper::lowX + lowXMargin)   &&
            (point.x() < GeometryHelper::highX - highXMargin) &&
            (point.y() > GeometryHelper::lowY + lowYMargin)   &&
            (point.y() < GeometryHelper::highY - highYMargin) &&
            (point.z() > GeometryHelper::lowZ + lowZMargin)   &&
            (point.z() < GeometryHelper::highZ - highZMargin) );
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetFiducialVolume()
{
    // Get the width of the FV in each direction - here we take off the fiducial margins (see AnalysisHelper::IsFiducial)
    const auto widthX = (GeometryHelper::highX - GeometryHelper::lowX) - 10.f - 10.f;
    const auto widthY = (GeometryHelper::highY - GeometryHelper::lowY) - 10.f - 10.f;
    const auto widthZ = (GeometryHelper::highZ - GeometryHelper::lowZ) - 10.f - 50.f;

    return widthX * widthY * widthZ;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsFiducial(const TVector3 &point)
{
    return AnalysisHelper::IsPointWithinMargins(point, 10.f, 10.f, 10.f, 10.f, 10.f, 50.f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const TVector3 &point)
{
    return AnalysisHelper::IsPointWithinMargins(point, 5.f, 5.f, 5.f, 5.f, 5.f, 5.f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const Event::Truth::Particle &particle)
{
    const auto start = TVector3(particle.startX(), particle.startY(), particle.startZ());
    const auto end = TVector3(particle.endX(), particle.endY(), particle.endZ());

    return (AnalysisHelper::IsContained(start) && AnalysisHelper::IsContained(end));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::HasTrackFit(const Event::Reco::Particle &particle)
{
    return (particle.startX.IsSet() && particle.startY.IsSet() && particle.startZ.IsSet() && particle.endX.IsSet() && particle.endY.IsSet() && particle.endZ.IsSet());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const Event::Reco::Particle &particle)
{
    if (!AnalysisHelper::HasTrackFit(particle))
        throw std::invalid_argument("AnalysisHelper::IsContained - input reco particle doesn't have a fitted track start-end points");

    const auto start = TVector3(particle.startX(), particle.startY(), particle.startZ());
    const auto end = TVector3(particle.endX(), particle.endY(), particle.endZ());

    return (AnalysisHelper::IsContained(start) && AnalysisHelper::IsContained(end));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsInActiveVolume(const TVector3 &point)
{
    return AnalysisHelper::IsPointWithinMargins(point, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetBestMatchedTruthParticleIndex(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    if (!recoParticle.hasMatchedMCParticle.IsSet())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle is from an event without truth info");

    if (!recoParticle.hasMatchedMCParticle())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle has no matched truth particle");

    if (truthParticles.size() != recoParticle.truthMatchPurities().size())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - wrong number of input truth particles");

    if (truthParticles.empty())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - no truth particles supplied but reco particle has a match!");


    unsigned int bestMatchedParticleId = std::numeric_limits<unsigned int>::max();
    float bestMatchScore = -std::numeric_limits<float>::max();
    bool foundMatch = false;

    for (unsigned int i = 0; i < truthParticles.size(); ++i)
    {
        if (applyVisibilityThreshold && !AnalysisHelper::PassesVisibilityThreshold(truthParticles.at(i)))
            continue;

        const auto score = recoParticle.truthMatchPurities().at(i);
        if (score < bestMatchScore)
            continue;

        bestMatchScore = score;
        bestMatchedParticleId = i;
        foundMatch = true;
    }

    if (!foundMatch)
        throw std::logic_error("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle has no matched truth particle");

    return bestMatchedParticleId;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

Event::Truth::Particle AnalysisHelper::GetBestMatchedTruthParticle(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    return truthParticles.at(AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles, applyVisibilityThreshold));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsGolden(const Event::Truth::Particle &particle)
{
    return (particle.nElasticScatters() == 0 &&
            particle.nInelasticScatters() == 0 &&
            particle.isStopping() &&
            AnalysisHelper::IsContained(particle));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TF1> AnalysisHelper::GetRangeToMomentumFunction()
{
    // Work out the maximum straight-line range a particle can have, so we can set sensible limits on the fit
    const auto minRange = 0.f;
    const auto maxRange = std::pow( std::pow(GeometryHelper::highX - GeometryHelper::lowX, 2) +
                                    std::pow(GeometryHelper::highY - GeometryHelper::lowY, 2) +
                                    std::pow(GeometryHelper::highZ - GeometryHelper::lowZ, 2) , 0.5f);

    std::shared_ptr<TF1> pFunc(new TF1(("fitFunc_" + std::to_string(m_fitFunctionIndex++)).c_str(), "[0] + [1]*x - [2]*pow(x, -[3])", minRange, maxRange));
    //std::shared_ptr<TF1> pFunc(new TF1(("fitFunc_" + std::to_string(m_fitFunctionIndex++)).c_str(), "[0] + [1]*x", minRange, maxRange));
    return pFunc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TF1> AnalysisHelper::GetRangeToMomentumFunctionMuon()
{
    RangeToMomentumFitParameters params;
    params.a = 0.19581;
    params.b = 0.002107;
    params.c = 0.15896;
    params.d = 0.19468;

    auto pFunc = AnalysisHelper::GetRangeToMomentumFunction();
    AnalysisHelper::SetRangeToMomentumFunctionParameters(params, pFunc);

    return pFunc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TF1> AnalysisHelper::GetRangeToMomentumFunctionPion()
{
    RangeToMomentumFitParameters params;
    params.a = 0.25798;
    params.b = 0.0024088;
    params.c = 0.18828;
    params.d = 0.11687;

    auto pFunc = AnalysisHelper::GetRangeToMomentumFunction();
    AnalysisHelper::SetRangeToMomentumFunctionParameters(params, pFunc);

    return pFunc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::SetRangeToMomentumFunctionParameters(const RangeToMomentumFitParameters &params, std::shared_ptr<TF1> &pFunc)
{
    pFunc->SetParameter(0, params.a);
    pFunc->SetParameter(1, params.b);
    pFunc->SetParameter(2, params.c);
    pFunc->SetParameter(3, params.d);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::GetRangeToMomentumFunctionParameters(const std::shared_ptr<TF1> &pFunc, RangeToMomentumFitParameters &params)
{
    params.a = pFunc->GetParameter(0);
    params.b = pFunc->GetParameter(1);
    params.c = pFunc->GetParameter(2);
    params.d = pFunc->GetParameter(3);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::GetRangeToMomentumFitParameters(const std::vector<float> &ranges, const std::vector<float> &momenta, RangeToMomentumFitParameters &params)
{
    if (ranges.size() != momenta.size())
        throw std::invalid_argument("AnalysisHelper::GetRangeToMomentumFitParameters - input ranges and momenta are of different size");

    // In general for a given range, there can be a significant (and asymmetric) spread of momenta around the most probable value.
    // The idea here is to first bin the ranges, and get the most probable value at each range, then fit those most probable values

    // The number of data points to put in each bin
    const unsigned int nPointsPerBin = 100u;

    // The maximum allowed bin width
    const float maxBinWidth = 10.f;

    // The momentum resolution scale to use when finding the MPV
    const float momentumScale = 0.01f;

    // First turn the two input vectors into a vector of pairs so they can be sorted
    std::vector< std::pair<float, float> > dataPoints;
    for (unsigned int i = 0; i < ranges.size(); ++i)
    {
        dataPoints.emplace_back(ranges.at(i), momenta.at(i));
    }

    // Now sort the points by range
    std::sort(dataPoints.begin(), dataPoints.end(), [](const auto &a, const auto &b){ return a.first < b.first; });

    // Now organise the data points into bins, each bin corresponds to a set of nearby ranges
    std::vector< std::vector< std::pair<float, float> > > binnedPoints;
    for (unsigned int i = 0; i < dataPoints.size(); ++i)
    {
        // If we have enough points in the current bin, then add another
        if (i % nPointsPerBin == 0)
            binnedPoints.emplace_back();

        // Add this data point to the last bin
        binnedPoints.back().push_back(dataPoints.at(i));
    }

    // Next find the most probable value in each bin
    std::vector<float> fitRanges, fitMomenta;
    for (const auto &bin : binnedPoints)
    {
        // Ignore the last bin if it doesn't have enough points
        if (bin.size() < nPointsPerBin)
            continue;

        // Extract the ranges and momenta in this bin
        std::vector<float> rangesInBin, momentaInBin;
        for (const auto &point : bin)
        {
            rangesInBin.push_back(point.first);
            momentaInBin.push_back(point.second);
        }

        // Get the spread of ranges in the bin
        const auto minRange = *std::min_element(rangesInBin.begin(), rangesInBin.end());
        const auto maxRange = *std::max_element(rangesInBin.begin(), rangesInBin.end());
        const auto binWidth = maxRange - minRange;

        // Don't use bins that have a large spread of values
        if (binWidth > maxBinWidth)
            continue;

        // Get the mean of the ranges in the bin
        const auto meanRange = std::accumulate(rangesInBin.begin(), rangesInBin.end(), 0.f) / static_cast<float>(rangesInBin.size());

        // Get the minimum momenta in the bin
        const auto minMomentum = *std::min_element(momentaInBin.begin(), momentaInBin.end());

        // Organise the momenta into their own bins
        std::map<unsigned int, std::vector<float> > momentumBins;
        for (const auto &momentum : momentaInBin)
        {
            const unsigned int binIndex = std::floor((momentum - minMomentum) / momentumScale);
            momentumBins[binIndex].push_back(momentum);
        }

        // Sort the momentum bins by the number of entries
        std::vector< std::pair<unsigned int, std::vector<float> > > momentumBinsVector;
        for (const auto &momentumBin : momentumBins)
        {
            momentumBinsVector.emplace_back(momentumBin.first, momentumBin.second);
        }
        std::sort(momentumBinsVector.begin(), momentumBinsVector.end(), [](const auto &a, const auto &b){ return a.second.size() > b.second.size(); });

        // Find the mean momenta in the most probable bin
        const auto mostProbableMomenta = momentumBinsVector.front().second;
        const auto meanMomentum = std::accumulate(mostProbableMomenta.begin(), mostProbableMomenta.end(), 0.f) / static_cast<float>(mostProbableMomenta.size());

        // Save this as a data point to fit
        fitRanges.push_back(meanRange);
        fitMomenta.push_back(meanMomentum);
    }

    // Determine the range over which to fit
    const auto minRange = *std::min_element(fitRanges.begin(), fitRanges.end());
    const auto maxRange = *std::max_element(fitRanges.begin(), fitRanges.end());


    // Make a TGraph to fit
    const auto nPoints = fitRanges.size();
    std::shared_ptr<TGraph> pGraph(new TGraph(nPoints, fitRanges.data(), fitMomenta.data()));

    // Get the function we want to fit
    auto pFunc = AnalysisHelper::GetRangeToMomentumFunction();

    // Set some initial values to help the fit converge quickly. For reference, I just eyeballed these numbers from a plot of muons
    params.a = 2e-1f;
    params.b = 2e-3f;
    params.c = 1e-1f;
    params.d = 2e-1f;
    AnalysisHelper::SetRangeToMomentumFunctionParameters(params, pFunc);

    // Run the fit
    //  - M = do it gooder (about as much information as the ROOT documentation gives you!)
    //  - N = don't draw the result
    //  - Q = quiet mode
    //  - S = return the fit result
    const auto fitResult = pGraph->Fit(pFunc->GetName(), "MNQS", "", minRange, maxRange);
    if (!fitResult->IsValid())
        std::cout << "AnalysisHelper::GetRangeToMomentumFitParameters - WARNING. Fit was deemed problematic, check the result!" << std::endl;

    // Store the fitted parameters
    AnalysisHelper::GetRangeToMomentumFunctionParameters(pFunc, params);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetPionMomentumFromRange(const float &range)
{
    auto pFunc = AnalysisHelper::GetRangeToMomentumFunctionPion();
    return pFunc->Eval(range);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetMuonMomentumFromRange(const float &range)
{
    auto pFunc = AnalysisHelper::GetRangeToMomentumFunctionMuon();
    return pFunc->Eval(range);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetMuonMomentumFromMCS(const Event::Reco::Particle &muon)
{
    if (muon.mcsMomentumForwardMuon.IsSet())
        return muon.mcsMomentumForwardMuon();

    if (muon.mcsMomentumBackwardMuon.IsSet())
        return muon.mcsMomentumBackwardMuon();

    // ATTN could throw here instead - would need to add event selection cut to remove ents in which muon momentum isn't available
    return -std::numeric_limits<float>::max();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetMuonMomentum(const Event::Reco::Particle &muon)
{
    if (AnalysisHelper::IsContained(muon))
        return AnalysisHelper::GetMuonMomentumFromRange(muon.range());

    return AnalysisHelper::GetMuonMomentumFromMCS(muon);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::AnalysisData AnalysisHelper::GetRecoAnalysisData(const Event::Reco &reco, const std::vector<int> &assignedPdgCodes, const bool passesGoldenPionSelection)
{
    AnalysisData data;

    // Sanity check
    const auto recoParticles = reco.particles;
    if (assignedPdgCodes.size() != recoParticles.size())
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisData - The assigned PDG codes is the wrong size");

    auto muonIndex = std::numeric_limits<unsigned int>::max();
    auto pionIndex = std::numeric_limits<unsigned int>::max();
    unsigned int nMuons = 0u;
    unsigned int nPions = 0u;
    unsigned int nProtons = 0u;
    unsigned int nOther = 0u;

    for (unsigned int index = 0; index < assignedPdgCodes.size(); ++index)
    {
        const auto recoPdg = assignedPdgCodes.at(index);

        switch (recoPdg)
        {
            case 13:
                muonIndex = index;
                nMuons++;
                break;
            case 211:
                pionIndex = index;
                nPions++;
                break;
            case 2212:
                nProtons++;
                break;
            default:
                nOther++;
                break;
        }
    }

    // Make sure the reco PDGs make sense
    if (nMuons != 1)
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisData - Reconstructed " + std::to_string(nMuons) + " muons");

    if (nPions != 1)
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisData - Reconstructed " + std::to_string(nMuons) + " pions");

    // Get the muon and pion reconstructed particles
    const auto muon = recoParticles.at(muonIndex);
    const auto pion = recoParticles.at(pionIndex);

    // Find the reconstructed variables
    const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
    data.muonCosTheta = muonDir.Z();
    data.muonPhi = std::atan2(muonDir.Y(), muonDir.X());
    data.muonMomentum = AnalysisHelper::GetMuonMomentum(muon);

    const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
    data.pionCosTheta = pionDir.Z();
    data.pionPhi = std::atan2(pionDir.Y(), pionDir.X());
    data.pionMomentum = AnalysisHelper::GetPionMomentumFromRange(pion.range());

    data.muonPionAngle = std::acos(muonDir.Dot(pionDir));
    data.nProtons = nProtons;
    data.hasGoldenPion = passesGoldenPionSelection;

    return data;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetParticleIndexWithPdg(const std::vector<int> &assignedPdgCodes, const int pdgCode)
{
    std::vector<unsigned int> indices;

    for (unsigned int i = 0; i < assignedPdgCodes.size(); ++i)
    {
        if (assignedPdgCodes.at(i) != pdgCode)
            continue;

        indices.push_back(i);
    }

    if (indices.size() != 1)
        throw std::logic_error("AnalysisHelper::GetParticleIndexWithPdg - Found " + std::to_string(indices.size()) + " particles width assigned PDG code: " + std::to_string(pdgCode) + ", expected 1");

    return indices.front();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::AnalysisData AnalysisHelper::GetTruthAnalysisData(const Event::Truth &truth, const bool useAbsPdg, const float protonMomentumThreshold)
{
    AnalysisData data;

    bool foundPion = false;
    bool foundMuon = false;

    data.nProtons = 0u;
    TVector3 pionDir, muonDir;

    for (const auto &particle : AnalysisHelper::SelectVisibleParticles(truth.particles))
    {
        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();

        const auto dir = TVector3(particle.momentumX(), particle.momentumY(), particle.momentumZ()).Unit();
        const auto cosTheta = dir.Z();
        const auto phi = std::atan2(dir.Y(), dir.X());

        if (pdg == 211)
        {
            if (foundPion)
                throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Found multiple pions! Are you sure this is a signal event?");

            foundPion = true;

            data.pionMomentum = particle.momentum();
            data.pionCosTheta = cosTheta;
            data.pionPhi = phi;
            data.hasGoldenPion = AnalysisHelper::IsGolden(particle);

            pionDir = dir;

            continue;
        }

        if (pdg == 13)
        {
            if (foundMuon)
                throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Found multiple muons! Are you sure this is a signal event?");

            foundMuon = true;

            data.muonMomentum = particle.momentum();
            data.muonCosTheta = cosTheta;
            data.muonPhi = phi;

            muonDir = dir;

            continue;
        }

        if (pdg == 2212)
        {
            const auto passesMomentumThreshold = (particle.momentum() > protonMomentumThreshold);
            if (!passesMomentumThreshold)
                continue;

            data.nProtons++;
            continue;
        }

        throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Found particle with unexpected PDG code = " + std::to_string(particle.pdgCode()));
    }

    if (!foundPion || !foundMuon)
        throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Input event doesn't contain both a muon and a pion!");


    data.muonPionAngle = std::acos(muonDir.Dot(pionDir));

    return data;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::AnalysisData AnalysisHelper::GetDummyAnalysisData()
{
    AnalysisData data;

    data.muonMomentum = -std::numeric_limits<float>::max();
    data.muonCosTheta = -std::numeric_limits<float>::max();
    data.muonPhi = -std::numeric_limits<float>::max();
    data.pionMomentum = -std::numeric_limits<float>::max();
    data.pionCosTheta = -std::numeric_limits<float>::max();
    data.pionPhi = -std::numeric_limits<float>::max();
    data.muonPionAngle = -std::numeric_limits<float>::max();
    data.nProtons = std::numeric_limits<unsigned int>::max();
    data.hasGoldenPion = false;

    return data;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::GetLogLikelihoodRatio(const Member<float> &numerator, const Member<float> &denominator, float &ratio)
{
    ratio = -std::numeric_limits<float>::max();

    if (!numerator.IsSet() || !denominator.IsSet())
        return false;

    if (numerator() < 0.f)
        return false;

    if (denominator() <= std::numeric_limits<float>::epsilon())
        return false;

    ratio = std::log(numerator() / denominator());
    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::GetSoftmax(const Member<float> &signal, const Member<float> &background, float &softmax)
{
    softmax = -std::numeric_limits<float>::max();

    if (!signal.IsSet() || !background.IsSet())
        return false;

    const auto denominator = std::exp(signal()) + std::exp(background());
    if (std::abs(denominator) <= std::numeric_limits<float>::epsilon())
        return false;

    softmax = std::exp(signal()) / denominator;
    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetCountUncertainty(const float &count)
{
    if (count < 0)
        throw std::logic_error("AnalysisHelper::GetCountUncertainty - count is < 0");

    return std::pow(count + 1, 0.5f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetEfficiencyUncertainty(const float &numerator, const float &denominator)
{
    if (denominator <= std::numeric_limits<float>::epsilon())
        throw std::logic_error("AnalysisHelper::GetEfficiencyUncertainty - denominator is <= 0");

    if (numerator > denominator)
        throw std::logic_error("AnalysisHelper::GetEfficiencyUncertainty - numerator > denomintaor");

    const auto n1 = numerator + 1;
    const auto n2 = numerator + 2;
    const auto d2 = denominator + 2;
    const auto d3 = denominator + 3;

    return std::pow((n1*n2)/(d2*d3) - (n1*n1)/(d2*d2), 0.5f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::PrintLoadingBar(const unsigned int numerator, const unsigned int denominator)
{
    if (denominator == 0)
        return;

    // Get the current stage of the loading bar
    const auto width = 85u;
    const auto loadedWidth = static_cast<unsigned int>(std::floor(width * (static_cast<float>(numerator) / static_cast<float>(denominator))));
    const auto unloadedWidth = static_cast<unsigned int>(width - loadedWidth);

    // Get the fraction of the way through the current stage
    const auto stageLength = static_cast<float>(denominator) / static_cast<float>(width);
    const unsigned int tenth = std::floor(10 * ((static_cast<float>(numerator) - stageLength * std::floor(static_cast<float>(numerator) / stageLength)) / stageLength));

    // Now do the same but for the last value of the numerator
    const auto lastLoadedWidth = static_cast<unsigned int>(std::floor(width * (static_cast<float>(numerator - 1) / static_cast<float>(denominator))));
    const unsigned int lastTenth = std::floor(10 * ((static_cast<float>(numerator - 1) - stageLength * std::floor(static_cast<float>(numerator - 1) / stageLength)) / stageLength));

    // Work out if it's worth printing
    const bool shouldPrint = (numerator == 0) || (lastLoadedWidth != loadedWidth) || (lastTenth != tenth);

    if (!shouldPrint)
        return;

    // Jump back 2 lines
    if (numerator != 0)
        std::cout << "\033[F\033[F";

    // Print the fraction and loading bar
    std::cout << numerator << " / " << denominator << std::endl;
    std::cout << "|" << std::string(loadedWidth, '=') << tenth << std::string(unloadedWidth, ' ') << "|" << std::endl;
}

} // namespace ubcc1pi
