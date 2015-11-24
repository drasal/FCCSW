#ifndef _DUMMYSIMULATION_H_
#define _DUMMYSIMULATION_H_

#include "GaudiAlg/GaudiAlgorithm.h"
#include "FWCore/DataHandle.h"
#include "datamodel/ParticleCollection.h"
#include "datamodel/MCParticleCollection.h"


class DummySimulation: public GaudiAlgorithm {
  friend class AlgFactory<DummySimulation> ;

public:
  /// Constructor.
  DummySimulation(const std::string& name, ISvcLocator* svcLoc);
  /// Initialize.
  virtual StatusCode initialize();
  /// Execute. This function actually does no simulation,
  /// and simply converts the stable MCParticles in the input collection
  /// into Particles that are written to the output collection. 
  virtual StatusCode execute();
  /// Finalize.
  virtual StatusCode finalize();
private:

  /// Handle for the MCParticleCollection to be read
  DataHandle<fccedm::MCParticleCollection> m_genphandle;
  /// Handle for the "reconstructed" to be written
  DataHandle<fccedm::ParticleCollection> m_recphandle;

};

#endif
