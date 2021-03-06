/**
 *  @file  ubcc1pi/Helpers/BacktrackHelper.h
 *
 *  @brief The header file for the back tracking helper class
 */

#ifndef UBCC1PI_HELPERS_BACKTRACK_HELPER
#define UBCC1PI_HELPERS_BACKTRACK_HELPER

#include "ubcc1pi/Helpers/CollectionHelper.h"

#include <unordered_map>
#include <tuple>

namespace ubcc1pi
{

/**
 *  @brief  The backtrack helper class
 *
 *          Helper functions for matching PFParticles to MCParticles
 */
class BacktrackHelper
{
    public:
        /**
         *  @brief  Class holding backtracking data
         */
        class BacktrackerData
        {
            public:
                typedef std::unordered_map<art::Ptr<simb::MCParticle>, float> MCParticleToFloatMap;   ///< Mapping from MCParticles to floats
                typedef std::unordered_map<art::Ptr<recob::PFParticle>, float> PFParticleToFloatMap;  ///< Mapping from PFParticles to floats
                typedef AssociationData<recob::PFParticle, simb::MCParticle, float> MatchMap;         ///< Mapping from PFParticles to MCParticles along with the shared weight

                /**
                 *  @brief  Constructor
                 *
                 *  @param  pfParticles the input pfParticles
                 *  @param  mcParticles the input mcParticles
                 *  @param  hitsToPfps the input mapping from hits to PFParticles
                 *  @param  hitsToMcps the input mapping from hits to MCParticles
                 */
                BacktrackerData(const PFParticleVector &pfParticles, const MCParticleVector &mcParticles, const HitsToPFParticles &hitsToPfps, const HitsToMCParticleWeights &hitsToMcps);

                /**
                 *  @brief  Get the MCParticles
                 *
                 *  @return the MCParticles
                 */
                MCParticleVector GetMCParticles() const;

                /**
                 *  @brief  Get the PFParticles
                 *
                 *  @return the PFParticles
                 */
                PFParticleVector GetPFParticles() const;

                /**
                 *  @brief  Get the number of hits associated with a given PFParticle
                 *
                 *  @param  pfParticle the pfParticle
                 *
                 *  @return the number of hits
                 */
                unsigned int GetNHits(const art::Ptr<recob::PFParticle> &pfParticle) const;

                /**
                 *  @brief  Get the weight associated with a given PFParticle (= the number of hits)
                 *
                 *  @param  pfParticle the pfParticle
                 *
                 *  @return the weight
                 */
                float GetWeight(const art::Ptr<recob::PFParticle> &pfParticle) const;

                /**
                 *  @brief  Get the weight associated with a given MCParticle (= the weighted number of hits)
                 *
                 *  @param  mcParticle the mcParticle
                 *
                 *  @return the weight
                 */
                float GetWeight(const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get all hits that are associated the MCParticle with any weight
                 *
                 *  @param  mcParticle the mcParticle
                 *
                 *  @return the hits
                 */
                HitVector GetHits(const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get all hits that are associated the PFParticle
                 *
                 *  @param  pfParticle the pfParticle
                 *
                 *  @return the hits
                 */
                HitVector GetHits(const art::Ptr<recob::PFParticle> &pfParticle) const;

                /**
                 *  @brief  Get the weight of a given PFParticle-MCParticle pair
                 *
                 *  @param  pfParticle the pfParticle
                 *  @param  mcParticle the mcParticle
                 *
                 *  @return the match weight
                 */
                float GetMatchWeight(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the purity of a given PFParticle-MCParticle pair = shared weight / PFParticle weight
                 *          Purity is the fraction of the PFParticle that represents the MCParticle
                 *
                 *  @param  pfParticle the pfParticle
                 *  @param  mcParticle the mcParticle
                 *
                 *  @return the match purity
                 */
                float GetMatchPurity(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the completeness a given PFParticle-MCParticle pair = shared weight / MCParticle weight
                 *          Completeness is the fraction of the MCParticle that is represented by the PFParticle
                 *
                 *  @param  pfParticle the pfParticle
                 *  @param  mcParticle the mcParticle
                 *
                 *  @return the match completeness
                 */
                float GetMatchCompleteness(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the MCParticle with the highest weighted match to the given PFParticle
                 *
                 *  @param  pfParticle the pfParticle
                 *
                 *  @return the MCParticle
                 */
                art::Ptr<simb::MCParticle> GetBestMatchedMCParticle(const art::Ptr<recob::PFParticle> &pfParticle) const;

                /**
                 *  @brief  Get the PFParticles which have the given MCParticle as their strongest match
                 *
                 *  @param  mcParticle the mcParticle
                 *
                 *  @return the PFParticles
                 */
                PFParticleVector GetBestMatchedPFParticles(const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the weight associated with the input Hit and MCParticle
                 *          weight = fraction of the energy of the hit that was provided by the supplied MCParticle
                 *
                 *  @param  hit the input hit
                 *  @param  mcParticle the input MCParticle
                 *
                 *  @return the weight
                 */
                float GetWeight(const art::Ptr<recob::Hit> &hit, const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the summed weight associated with the input collection of Hits and MCParticle
                 *
                 *  @param  hits the input hits
                 *  @param  mcParticle the input MCParticle
                 *
                 *  @return the weight
                 */
                float GetWeight(const HitVector &hits, const art::Ptr<simb::MCParticle> &mcParticle) const;

            private:

                /**
                 *  @brief  Check that PFParticles in the input map from hits -> PFParticles all exist in the input PFPartice vector
                 *
                 *  @param  pfParticles the input PFParticles
                 *  @param  hitsToPfps the input mapping from hits to PFPs
                 *
                 *  @throws if an inconsistency is found
                 */
                void CheckConsistency(const PFParticleVector &pfParticles, const HitsToPFParticles &hitsToPfps) const;

                /**
                 *  @brief  Check that MCParticles in the input map from hits -> MCParticles all exist in the input MCPartice vector
                 *
                 *  @param  mcParticles the input MCParticles
                 *  @param  hitsToMcps the input mapping from hits to MCPs
                 *
                 *  @throws if an inconsistency is found
                 */
                void CheckConsistency(const MCParticleVector &mcParticles, const HitsToMCParticleWeights &hitsToMcps) const;

                /**
                 *  @brief  Collect all of the hits stored in the input maps ensuring we don't double count
                 *
                 *  @param  hitsToPfps the input mapping from hits to PFParticles
                 *  @param  hitsToMcps the input mapping from hits to MCParticles
                 *
                 *  @return an ordered vector of hits that exist in the input maps
                 */
                HitVector CollectHits(const HitsToPFParticles &hitsToPfps, const HitsToMCParticleWeights &hitsToMcps) const;

                /**
                 *  @brief  Get the PFParticle associated with a given hit
                 *
                 *  @param  hit the input hit
                 *  @param  hitsToPfps the input mapping from hits to PFParticles
                 *  @param  outputParticle the output particle associated to the hit
                 *
                 *  @return if a PFParticle could be found (this will be false if the hit is unclustered)
                 */
                bool CollectPFParticle(const art::Ptr<recob::Hit> &hit, const HitsToPFParticles &hitsToPfps, art::Ptr<recob::PFParticle> &outputParticle) const;

                /**
                 *  @brief  Collect the MCParticles that are associated with a given hit along with their weights (= fraction of the hits charge contributed by the MCParticle)
                 *
                 *  @param  hit the input hit
                 *  @param  hitsToMcps the input mapping from hits to MCParticle
                 *  @param  outputParticleWeights the output collection of associated MCParticle-weight pairs
                 *
                 *  @return if any associated MCParticles could be found (this will be false if the hit is from an external background, overlays, noise)
                 */
                bool CollectMCParticleWeights(const art::Ptr<recob::Hit> &hit, const HitsToMCParticleWeights &hitsToMcps, CollectionData<simb::MCParticle, float> &outputParticleWeights) const;

                MCParticleVector        m_mcParticles;         ///< The MCParticles
                PFParticleVector        m_pfParticles;         ///< The PFParticles
                MCParticleToFloatMap    m_mcParticleWeightMap; ///< The mapping from MCParticle to the total weight (= weighted number of hits)
                PFParticleToFloatMap    m_pfParticleWeightMap; ///< The mapping from PFParticle to the total weight (= number of hits)
                MatchMap                m_matchMap;            ///< The matching between PFParticles and MCParticles along with the shared weight
                HitsToMCParticleWeights m_hitsToMcps;          ///< The mapping from hits to MCParticle along with the associated weights
                HitsToPFParticles       m_hitsToPfps;          ///< The mapping from hits to PFParticles
        };

        /**
         *  @brief  Class holding the backtracking details of the slices
         */
        class SliceMetadata
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  slices the input slices
                 *  @param  sliceToIsSelectedAsNu the input map from a slice to whether it's been selected as the neutrino - There should only ever be 1 or zero slices selected
                 *  @param  slicesToHits the input mapping from hits to slices
                 *  @param  hitsToIsNuInduced the input mapping from a hit to whether it's neutrino induced or not
                 */
                SliceMetadata(const SliceVector &slices, const SlicesToBool &sliceToIsSelectedAsNu, const SlicesToHits &slicesToHits, const HitsToBool &hitsToIsNuInduced);

                /**
                 *  @brief  Get the slices
                 *
                 *  @return the slices
                 */
                SliceVector GetSlices() const;

                /**
                 *  @brief  Get the purity of a given slice
                 *          purity = fraction of the hits in the slice that are neutrino induced
                 *
                 *  @param  slice the input slice
                 *
                 *  @return the purity
                 */
                float GetPurity(const art::Ptr<recob::Slice> &slice) const;

                /**
                 *  @brief  Get the completeness of a given slice
                 *          completeness = fraction of the neutrino induced hits in the event that are in the slice
                 *
                 *  @param  slice the input slice
                 *
                 *  @return the completeness
                 */
                float GetCompleteness(const art::Ptr<recob::Slice> &slice) const;

                /**
                 *  @brief  Get the number of hits in the slice
                 *
                 *  @param  slice the input slice
                 *
                 *  @return the number of hits
                 */
                unsigned int GetNumberOfHits(const art::Ptr<recob::Slice> &slice) const;

                /**
                 *  @brief  Get the hits in the given slice
                 *
                 *  @param  slice the input slice
                 *
                 *  @return the hits
                 */
                HitVector GetHits(const art::Ptr<recob::Slice> &slice) const;

                /**
                 *  @brief  Get the number of neutrino induced hits in the slice
                 *
                 *  @param  slice the input slice
                 *
                 *  @return the number of neutrino induced hits
                 */
                unsigned int GetNumberOfNuInducedHits(const art::Ptr<recob::Slice> &slice) const;

                /**
                 *  @brief  Get the total number of neutrino induced hits
                 *
                 *  @return the total number of neutrino induced hits
                 */
                unsigned int GetTotalNumberOfNuInducedHits() const;

                /**
                 *  @brief  Get the slice with the largest completeness
                 *
                 *  @return the slice
                 */
                art::Ptr<recob::Slice> GetMostCompleteSlice() const;

                /**
                 *  @brief  Get the slices that have been selected a neutrinos
                 *
                 *  @return the slices
                 */
                SliceVector GetSelectedNeutrinoSlices() const;

                /**
                 *  @brief  Determine if the most complete slice has been selected as a neutrino
                 *
                 *  @return bool, true if most complete slicec is selected
                 */
                bool IsMostCompleteSliceSelected() const;

            private:

                /**
                 *  @brief  Get all hits in the event through the slices
                 *
                 *  @return the hits
                 */
                HitVector GetAllHits() const;

                /**
                 *  @brief  Count the number of hits in the input vector that are neutrino induced
                 *
                 *  @param  hits the input hits
                 *
                 *  @return the number of neutrino induced hits
                 */
                unsigned int CountNuInducedHits(const HitVector &hits) const;

                SliceVector  m_slices;                  ///< The slices
                SlicesToBool m_sliceToIsSelectedAsNu;   ///< The mapping from slices to a boolean noting if the slice was selected as a neutrino
                SlicesToHits m_slicesToHits;            ///< The mapping from slices to hits in the slice
                HitsToBool   m_hitsToIsNuInduced;       ///< The mapping from hits to a boolean noting if the hit has any of it's energy contributed by a true neutrino induced particle
        };

    /**
     *  @brief  Get the mapping from the input neutrino final state PFParticle to their hits, optionally folding in downstream PFParticles
     *
     *  @param  event the art event
     *  @param  pfParticleLabel the PFParticle producer label
     *  @param  finalStates the input vector of final state PFParticles
     *  @param  useDaughters whether to fold in downstream particles
     *
     *  @return the mapping
     */
    static HitsToPFParticles GetHitToPFParticleMap(const art::Event &event, const art::InputTag &pfParticleLabel, const PFParticleVector &finalStates, const bool &useDaughters);

    /**
     *  @brief  Get the mapping from neutrino final state MCParticles to hits, folding in all downstream particles
     *
     *  @param  event the art event
     *  @param  mcParticleLabel the MCParticle producer label
     *  @param  backtrackerLabel the MCParticle to hit producer label
     *  @param  finalStates the input vector of final state MCParticles
     *
     *  @return the hit to MCParticle weight map
     */
    static HitsToMCParticleWeights GetHitToMCParticleWeightMap(const art::Event &event, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const MCParticleVector &finalStates);

    /**
     *  @brief  Convert the input map from having backtracker hit matching data, to just a float for the MCPartile -> Hit weight
     *
     *  @param  mcParticleToHits the input mapping from MCParticles to hits
     *
     *  @return the output mapping in the required format
     */
    static MCParticlesToHitWeights GetMCParticleToHitWeightsMap(const AssociationData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> &mcParticleToHits);

    /**
     *  @brief  Get the mapping from all hits, to whether the hit is induced (at least in part) by a neutrino induced MCParticle
     *
     *  @param  event the art event
     *  @param  hitLabel the Hit producer label
     *  @param  mcParticleLabel the MCParticle producer label
     *  @param  backtrackerLabel the MCParticle to hit producer label
     *  @param  nuMCParticles the input vector of neutrino induced MCParticles
     *
     *  @return the hits to is neutrino induced map
     */
    static HitsToBool GetHitsToIsNuInducedMap(const art::Event &event, const art::InputTag &hitLabel, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const MCParticleVector &nuMCParticles);

    /**
     *  @brief  Count the number of hits associated with a given MCParticle in a given view
     *
     *  @param  mcParticle the MCParticle in question
     *  @param  mcParticleToHits the mapping from MCParticles to hits
     *  @param  view the view to use
     *
     *  @return the number of hits
     */
    static int CountHitsInView(const art::Ptr<simb::MCParticle> &mcParticle, const AssociationData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> &mcParticleToHits, const geo::View_t &view);

    /**
     *  @brief  Count the number of "good" hits associated with a given MCParticle in a given view for which the MCParticle contributed at least 1/2 of the charge of the hit
     *
     *  @param  mcParticle the MCParticle in question
     *  @param  mcParticleToHits the mapping from MCParticles to hits
     *  @param  view the view to use
     *
     *  @return the number of good hits
     */
    static int CountGoodHitsInView(const art::Ptr<simb::MCParticle> &mcParticle, const AssociationData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> &mcParticleToHits, const geo::View_t &view);

    /**
     *  @brief  Count the total hit weight associated with a given MCParticle in a given view
     *
     *  @param  mcParticle the MCParticle in question
     *  @param  mcParticleToHits the mapping from MCParticles to hits
     *  @param  view the view to use
     *
     *  @return the hit weight
     */
    static float GetHitWeightInView(const art::Ptr<simb::MCParticle> &mcParticle, const AssociationData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> &mcParticleToHits, const geo::View_t &view);

    /**
     *  @brief  Count the number of hits associated with a given MCParticle in a given view
     *
     *  @param  mcParticle the MCParticle in question
     *  @param  mcParticleToHits the mapping from MCParticles to hits
     *  @param  view the view to use
     *
     *  @return the number of hits
     */
    static int CountHitsInView(const art::Ptr<simb::MCParticle> &mcParticle, const MCParticlesToHitWeights &mcParticleToHits, const geo::View_t &view);

    /**
     *  @brief  Count the number of "good" hits associated with a given MCParticle in a given view for which the MCParticle contributed at least 1/2 of the charge of the hit
     *
     *  @param  mcParticle the MCParticle in question
     *  @param  mcParticleToHits the mapping from MCParticles to hits
     *  @param  view the view to use
     *
     *  @return the number of good hits
     */
    static int CountGoodHitsInView(const art::Ptr<simb::MCParticle> &mcParticle, const MCParticlesToHitWeights &mcParticleToHits, const geo::View_t &view);

    /**
     *  @brief  Count the total hit weight associated with a given MCParticle in a given view
     *
     *  @param  mcParticle the MCParticle in question
     *  @param  mcParticleToHits the mapping from MCParticles to hits
     *  @param  view the view to use
     *
     *  @return the hit weight
     */
    static float GetHitWeightInView(const art::Ptr<simb::MCParticle> &mcParticle, const MCParticlesToHitWeights &mcParticleToHits, const geo::View_t &view);
};

} // namespace ubcc1pi

#endif
