/**
 *  @file  ubcc1pi/Producers/HitSlicer_module.cc
 *
 *  @brief The implementation of the hit slicer producer
 */

#include "ubcc1pi/Producers/HitSlicer.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"

#include "art/Persistency/Common/PtrMaker.h"

namespace ubcc1pi
{

HitSlicer::HitSlicer(const art::EDProducer::Table<Config> &config) :
    m_config(config)
{
    produces< std::vector<recob::Hit> >();

    if (m_config().HasMC())
        produces< art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> >();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void HitSlicer::produce(art::Event &event)
{
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << " Run    : " << event.run() << std::endl;
    std::cout << " SubRun : " << event.subRun() << std::endl;
    std::cout << " Event  : " << event.event() << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    // Setup the output collections
    std::unique_ptr< std::vector<recob::Hit> > outputHits( new std::vector<recob::Hit> );
    std::unique_ptr< art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> > outputMCPsToHits( new art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> );

    // Get the input reco collections
    const auto slices = CollectionHelper::GetCollection<recob::Slice>(event, m_config().SliceLabel());
    const auto slicesToPFPs = CollectionHelper::GetAssociation<recob::Slice, recob::PFParticle>(event, m_config().SliceToPFParticlesLabel());
    const auto slicesToHits = CollectionHelper::GetAssociation<recob::Slice, recob::Hit>(event, m_config().SliceToHitsLabel());

    // Get the input truth collections (if have MC)
    AssociationData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> mcpToHitsData;
    if (m_config().HasMC())
    {
        mcpToHitsData = CollectionHelper::GetAssociationWithData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>(event, m_config().MCParticleLabel(), m_config().BacktrackerLabel());
    }
    const auto hitsToMCPData = CollectionHelper::GetReversedAssociation(mcpToHitsData);

    // Find the reconstructed neutrino slice (if it exists)
    const art::PtrMaker<recob::Hit> makePtrHit(event);
    bool foundNeutrino = false;
    for (const auto &slice : slices)
    {
        const auto pfps = CollectionHelper::GetManyAssociated(slice, slicesToPFPs);
        const auto neutrinos = RecoHelper::GetNeutrinos(pfps);

        if (neutrinos.empty())
            continue;

        if (neutrinos.size() != 1)
            throw cet::exception("HitSlicer::produce") << " - Found multiple neutrino PFPs in a slice." << std::endl;

        if (foundNeutrino)
            throw cet::exception("HitSlicer::produce") << " - Found multiple neutrino PFPs in the event." << std::endl;

        foundNeutrino = true;

        // Add the hits to the output
        const auto hits = CollectionHelper::GetManyAssociated(slice, slicesToHits);
        for (const auto &hit : hits)
        {
            if (m_config().HasMC())
            {
                const auto mcpDataVector = CollectionHelper::GetManyAssociatedWithData(hit, hitsToMCPData);

                for (const auto &mcpData : mcpDataVector)
                    outputMCPsToHits->addSingle(mcpData.first, makePtrHit(outputHits->size()), mcpData.second);

            }

            outputHits->push_back(*hit);
        }
    }

    // Save the output
    event.put(std::move(outputHits));

    if (m_config().HasMC())
        event.put(std::move(outputMCPsToHits));
}

} // namespace ubcc1pi
