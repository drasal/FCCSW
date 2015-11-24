#ifndef _DELPHESSIMULATION_H_
#define _DELPHESSIMULATION_H_

#include "TObjArray.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"
//#include "datamodel/GenJetCollection.h"
#include "datamodel/ParticleCollection.h"
#include "datamodel/ParticleMCAssociationCollection.h"
#include "datamodel/JetCollection.h"
#include "datamodel/JetParticleAssociationCollection.h"
//#include "datamodel/MissingEnergyCollection.h"

#include "TFile.h"
#include "FWCore/DataHandle.h"
#include "HepMC/GenEvent.h"

// for Delphes

// ROOT
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TTree.h"
#include "TApplication.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "Math/GenVector/LorentzVector.h"

#include "classes/DelphesModule.h"
//#include "SortableObject.h"
#include "classes/DelphesClasses.h"
#include "modules/Delphes.h"
#include "classes/DelphesFactory.h"
#include "DelphesExtHepMCReader.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTask.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "modules/ParticlePropagator.h"

// STL
#include <iostream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <signal.h>

/** @class DelphesSimulation
 *
 *  Gaudi algorithm: reads-in HepMC event (from file or memory), calls Dephes and performs Delphes simulation,
 *                   output to Delphes-based root file or FCC-EDM-based root file (recommended!).
 *
 *  @author: Z. Drasal (CERN), M. de Gruttola (CERN)
 *
 */
class DelphesSimulation: public GaudiAlgorithm {
  friend class AlgFactory<DelphesSimulation> ;

public:

  // Constructor.
  DelphesSimulation(const std::string& name, ISvcLocator* svcLoc);

  // Initialize.
  virtual StatusCode initialize();

  // Execute. This function actually does no simulation,
  // and simply converts the stable MCParticles in the input collection
  // into Particles that are written to the output collection.
  virtual StatusCode execute();

  // Finalize.
  virtual StatusCode finalize();

private:

  // Delphes detector card to be read in
  std::string              m_DelphesCard;

  // Delphes
  Delphes               *m_Delphes ;
  DelphesFactory        *m_DelphesFactory ;
  DelphesExtHepMCReader *m_HepMCReader ;

  // Handle for the HepMC to be read in from data store
  DataHandle<HepMC::GenEvent> m_hepmcHandle;

  // Or read a HepMC file
  FILE*        m_inHepMCFile ;
  std::string  m_inHepMCFileName;
  Long64_t     m_inHepMCFileLength;

  // Handle for the generated or reconstructed objects to be written out
  DataHandle<fccedm::MCParticleCollection> m_handleGenParticles;    // Generated particles handle
  DataHandle<fccedm::GenVertexCollection>  m_handleGenProdVertices; // Generated particle production vertices handle
  DataHandle<fccedm::GenVertexCollection>  m_handleGenDecVertices;  // Generated particle decay vertices handle
  DataHandle<fccedm::ParticleCollection>   m_handleRecMuons;        // Reconstructed muons
  DataHandle<fccedm::ParticleCollection>   m_handleRecElectrons;    // Reconstructed electrons
  DataHandle<fccedm::ParticleCollection>   m_handleRecPhotons;      // Reconstructed photons
  DataHandle<fccedm::JetCollection>        m_handleRecJets;         // Reconstructed jets

  DataHandle<fccedm::ParticleMCAssociationCollection>  m_handleRecMuonsToMC;     // Relation between muons & MC particle
  DataHandle<fccedm::ParticleMCAssociationCollection>  m_handleRecElectronsToMC; // Relation between electrons & MC particle
  DataHandle<fccedm::ParticleMCAssociationCollection>  m_handleRecPhotonsToMC;   // Relation between photons & MC particle
  DataHandle<fccedm::JetParticleAssociationCollection> m_handleRecJetsToPart;    // Relation between jets & particle
  //DataHandle<ParticleCollection> m_recMETs;               // MET
  //DataHandle<ParticleCollection> m_recHTSs;               // Scalar HT

  // Convert internal Delphes objects to FCC EDM objects - UPDATE
  void ConvertMCParticles(const TObjArray*  Input, fccedm::MCParticleCollection* colMCParticles, fccedm::GenVertexCollection* colProdVertices, fccedm::GenVertexCollection* colDecVertices);
  template <class T> void ConvertParticles(const TObjArray*  Input, fccedm::ParticleCollection* colParticles, fccedm::MCParticleCollection* colMCParticles, fccedm::ParticleMCAssociationCollection* ascColParticlesToMC);
  void ConvertPhotons(    const TObjArray*  Input, fccedm::ParticleCollection* colParticles, fccedm::MCParticleCollection* colMCParticles, fccedm::ParticleMCAssociationCollection* ascColParticlesToMC);
  void ConvertJets(       const TObjArray*  Input, fccedm::JetCollection* colJets);
  //void  ConvertMET(     TObjArray *  Input, ParticleCollection * coll );
  //void  ConvertHT(      TObjArray *  Input, ParticleCollection * coll );

  int m_eventCounter;

  // ROOT output
  TFile            *m_outRootFile;
  std::string       m_outRootFileName;
  ExRootTreeWriter *m_treeWriter;
  ExRootTreeBranch *m_branchEvent;
  ExRootConfReader *m_confReader;

  // Internal Delphes object arrays
  TObjArray *m_stablePartOutArray,
            *m_allPartOutArray ,
            *m_partonOutArray ,
            *m_muonOutArray ,
            *m_electronOutArray ,
            *m_photonOutArray ,
            *m_jetOutArray ,
            *m_metOutArray ,
            *m_htOutArray ;
};

#endif
