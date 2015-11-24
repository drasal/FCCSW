#include "JetClustering.h"

// #include "DataObjects/LorentzVector.h"


// #include "TLorentzVector.h"

// #include <vector>

// using namespace fastjet;

#include "datamodel/ParticleCollection.h"
#include "datamodel/BareParticle.h"
#include "datamodel/JetCollection.h"
#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenJetCollection.h"
#include "datamodel/JetParticleAssociationCollection.h"
#include "datamodel/GenJetParticleAssociationCollection.h"

typedef JetClustering<fccedm::ParticleCollection, fccedm::JetCollection, fccedm::JetParticleAssociationCollection> ParticleJetClustering;
typedef JetClustering<fccedm::MCParticleCollection, fccedm::GenJetCollection, fccedm::GenJetParticleAssociationCollection> MCParticleGenJetClustering;

DECLARE_COMPONENT(ParticleJetClustering)
DECLARE_COMPONENT_WITH_ID(ParticleJetClustering,"ParticleJetClustering")
DECLARE_COMPONENT(MCParticleGenJetClustering)
DECLARE_COMPONENT_WITH_ID(MCParticleGenJetClustering,"MCParticleGenJetClustering")

