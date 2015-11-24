#ifndef _GENPARTICLEFILTER_H_
#define _GENPARTICLEFILTER_H_

#include "GaudiAlg/GaudiAlgorithm.h"
#include "FWCore/DataHandle.h"
#include "datamodel/MCParticleCollection.h"

class GenParticleFilter: public GaudiAlgorithm {
  friend class AlgFactory<GenParticleFilter> ;

public:
  /// Constructor.
  GenParticleFilter(const std::string& name, ISvcLocator* svcLoc);
  /// Initialize.
  virtual StatusCode initialize();
  /// Execute.
  virtual StatusCode execute();
  /// Finalize.
  virtual StatusCode finalize();
private:
  /// Handle for the ParticleCollection to be read
  DataHandle<fccedm::MCParticleCollection> m_igenphandle;
  /// Handle for the genparticles to be written
  DataHandle<fccedm::MCParticleCollection> m_ogenphandle;
};

#endif
