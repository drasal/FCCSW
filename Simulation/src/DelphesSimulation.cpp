#include "DelphesSimulation.h"

#include <limits>
#include <ParticleStatus.h>

using namespace std;

DECLARE_COMPONENT(DelphesSimulation)

DelphesSimulation::DelphesSimulation(const std::string& name, ISvcLocator* svcLoc):
GaudiAlgorithm(name, svcLoc) ,
  m_DelphesCard(),
  m_Delphes(nullptr),
  m_DelphesFactory(nullptr),
  m_HepMCReader(nullptr),
  m_inHepMCFile(nullptr),
  m_inHepMCFileName(""),
  m_inHepMCFileLength(0),
  m_outRootFile(nullptr),
  m_outRootFileName(""),
  m_treeWriter(nullptr),
  m_branchEvent(nullptr),
  m_confReader(nullptr)
 {
   //declareProperty("filename", m_filename="" , "Name of the HepMC file to read");
   declareProperty("DelphesCard"      , m_DelphesCard    , "Name of Delphes tcl config file with detector and simulation parameters");
   declareProperty("HepMCInputFile"   , m_inHepMCFileName, "Name of HepMC input file, if defined file read in / if not data read in directly from the transient data store");
   declareProperty("ROOTOutputFile"   , m_outRootFileName, "Name of Root output file, if defined file write out / if not data written to the transient data store");

   declareInput("hepmc", m_hepmcHandle);
   
   declareOutput("genParticles"      , m_handleGenParticles);
   declareOutput("genProdVertices"   , m_handleGenProdVertices);
   declareOutput("genDecVertices"    , m_handleGenDecVertices);
   declareOutput("recMuons"          , m_handleRecMuons);
   declareOutput("recElectrons"      , m_handleRecElectrons);
   declareOutput("recPhotons"        , m_handleRecPhotons);
   declareOutput("recJets"           , m_handleRecJets);

   declareOutput("recMuonsToMC"      , m_handleRecMuonsToMC);
   declareOutput("recElectronsToMC"  , m_handleRecElectronsToMC);
   declareOutput("recPhotonsToMC"    , m_handleRecPhotonsToMC);
   declareOutput("recJetsToPart"     , m_handleRecJetsToPart);
   /*declareOutput("recMETs"           , m_recMETHandle);
   declareOutput("recHTs"            , m_recHTSHandle);*/

   m_stablePartOutArray = nullptr;
   m_allPartOutArray    = nullptr;
   m_partonOutArray     = nullptr;
   m_muonOutArray       = nullptr;
   m_electronOutArray   = nullptr;
   m_photonOutArray     = nullptr;
   m_jetOutArray        = nullptr;
   m_metOutArray        = nullptr;
   m_htOutArray         = nullptr;

   m_eventCounter = 0;
}

StatusCode DelphesSimulation::initialize() {

  // Open HepMC file if defined
  if (m_inHepMCFileName!="") {

    info()  << "Reading in HepMC file: " << m_inHepMCFileName << endmsg;
    m_inHepMCFile = fopen(m_inHepMCFileName.c_str(), "r");

    if (m_inHepMCFile==nullptr) {

      error() << "Can't open " << m_inHepMCFileName << endmsg;
      return Error ("ERROR, can't open defined HepMC input file.");
    }
  
    fseek(m_inHepMCFile, 0L, SEEK_END);
    m_inHepMCFileLength = ftello(m_inHepMCFile);
    fseek(m_inHepMCFile, 0L, SEEK_SET);
    info() << "Length of HepMC input file: " << m_inHepMCFileLength << endmsg;
    if (m_inHepMCFileLength<=0) {
  
      fclose(m_inHepMCFile);
      return Error ("ERROR, zero length HepMC input file.");
    }
  }
  
  // If required, export output directly to root file
  if (m_outRootFileName!="") {

    info()  << "Opening ROOT output file: " << m_outRootFileName << endmsg;
    m_outRootFile = new TFile(m_outRootFileName.c_str(), "RECREATE");
    if (m_outRootFile->IsZombie()) {

      error() << "Can't open " << m_outRootFileName << endmsg;
      return Error ("ERROR, can't open defined ROOT output file.");
    }
  }

  // Read Delphes configuration card
  m_confReader = new ExRootConfReader;
  m_confReader->ReadFile(m_DelphesCard.c_str());
   
  // Instance of Delphes
  m_Delphes = new Delphes("Delphes");
  m_Delphes->SetConfReader(m_confReader);

  // Get standard Delphes factory
  m_DelphesFactory = m_Delphes->GetFactory();

  // Delphes needs data structure to be defined (ROOT tree)
  m_treeWriter  = new ExRootTreeWriter( m_outRootFile , "DelphesSim");
  m_branchEvent = m_treeWriter->NewBranch("Event", HepMCEvent::Class());
  m_Delphes->SetTreeWriter(m_treeWriter);

  // Define event readers
  //
  //  HepMC reader --> reads either from a file or directly from data store
  m_HepMCReader    = new DelphesExtHepMCReader;
  if (m_inHepMCFile) m_HepMCReader->SetInputFile(m_inHepMCFile);
  
  // Create following arrays of Delphes objects --> starting objects
  m_allPartOutArray    = m_Delphes->ExportArray("allParticles");
  m_stablePartOutArray = m_Delphes->ExportArray("stableParticles");
  m_partonOutArray     = m_Delphes->ExportArray("partons");

  // Init Delphes - read in configuration & define modules to be executed
  m_Delphes->InitTask();

  // Print Delphes modules to be used
  ExRootConfParam param = m_confReader->GetParam("::ExecutionPath");
  Long_t          size  = param.GetSize();
  info()  << "Delphes simulation will use the following modules: " << endmsg;
  for( Long_t k = 0; k < size; ++k) {

    TString name = param[k].GetString();
    info()  << "-- Module: " <<  name << endmsg;
  }
  
  // Initialize all variables
  m_eventCounter = 0;
  if (m_outRootFile!=nullptr) m_treeWriter->Clear();
  m_Delphes->Clear();
  m_HepMCReader->Clear();
 
  return StatusCode::SUCCESS;
}


StatusCode DelphesSimulation::execute() {

  //
  // Read event & initialize event variables
  TStopwatch readStopWatch;
  readStopWatch.Start();

  bool isEventReady = false;

  if (m_inHepMCFile) {

    // Test end-of-file
    if ( ftello(m_inHepMCFile) == m_inHepMCFileLength) {

      info() << "End of file reached at lenght " << m_inHepMCFileLength << endmsg;
      return StatusCode::SUCCESS;
    }

    // Read event - read line-by-line until event complete
    isEventReady = m_HepMCReader->ReadEventFromFile(m_DelphesFactory, m_allPartOutArray, m_stablePartOutArray, m_partonOutArray);
  }
  else {

    // Read event
    const HepMC::GenEvent *hepMCEvent = m_hepmcHandle.get();
    isEventReady = m_HepMCReader->ReadEventFromStore(hepMCEvent, m_DelphesFactory, m_allPartOutArray, m_stablePartOutArray, m_partonOutArray);

    // Print HepMC event info
    /*for(auto ipart=hepMCEvent->particles_begin(); ipart!=hepMCEvent->particles_end(); ++ipart) {

      int motherID        = -1;
      int motherIDRange   = 0;
      int daughterID      = -1;
      int daughterIDRange = 0;
      if ((*ipart)->production_vertex()!=nullptr) {

        motherID      = (*((*ipart)->production_vertex()->particles_in_const_begin()))->barcode();//(*((*ipart)->production_vertex()->particles_begin()))->barcode();
        motherIDRange = (*ipart)->production_vertex()->particles_in_size() -1;
      }
      if ((*ipart)->end_vertex()!=nullptr) {

        daughterID      = (*((*ipart)->end_vertex()->particles_out_const_begin()))->barcode();//(*((*ipart)->production_vertex()->particles_begin()))->barcode();
        daughterIDRange = (*ipart)->end_vertex()->particles_out_size() -1;
      }

      std::cout << "Delphes HepMC: "
                << " Id: "       << setw(3)  << (*ipart)->barcode()
                << " Pdg: "      << setw(5)  << (*ipart)->pdg_id()
                << " Mothers: "  << setw(3)  << motherID << " -> " << setw(3) << motherID+motherIDRange
                << " Daughters: "<< setw(3)  << daughterID << " -> " << setw(3) << daughterID+daughterIDRange
                << " Stat: "     << setw(2)  << (*ipart)->status()
                << " Px: "       << setw(10) << (*ipart)->momentum().px()
                << " Py: "       << setw(10) << (*ipart)->momentum().py()
                << " Pz: "       << setw(10) << (*ipart)->momentum().pz()
                << " E: "        << setw(10) << (*ipart)->momentum().e()
                << " M: "        << setw(10) << (*ipart)->momentum().m();
      if ((*ipart)->production_vertex()!=nullptr) {
      std::cout << " Vx: "       << setw(10) << (*ipart)->production_vertex()->position().x()
                << " Vy: "       << setw(10) << (*ipart)->production_vertex()->position().y()
                << " Vz: "       << setw(10) << (*ipart)->production_vertex()->position().z()
                << " T: "        << setw(10) << (*ipart)->production_vertex()->position().t();
      }
      std::cout << std::endl;
    }*/
  }

  if (!isEventReady) return StatusCode::FAILURE;

  // Print Delphes event info
  for (auto i=0; i<m_allPartOutArray->GetEntries(); i++) {

    Candidate *candidate = static_cast<Candidate *>(m_allPartOutArray->At(i));

    std::cout << "Delphes Object: "
              << " Id: "       << setw(3)  << i+1
              << " Pdg: "      << setw(5)  << candidate->PID
              << " Mothers: "  << setw(3)  << candidate->M1+1 << " -> " << setw(3) << candidate->M2+1
              << " Daughters: "<< setw(3)  << candidate->D1+1 << " -> " << setw(3) << candidate->D2+1
              << " Stat: "     << setw(2)  << candidate->Status
              << " Px: "       << setw(10) << candidate->Momentum.Px()
              << " Py: "       << setw(10) << candidate->Momentum.Py()
              << " Pz: "       << setw(10) << candidate->Momentum.Pz()
              << " E: "        << setw(10) << candidate->Momentum.E()
              << " M: "        << setw(10) << candidate->Mass
              << " Vx: "       << setw(10) << candidate->Position.X()
              << " Vy: "       << setw(10) << candidate->Position.Y()
              << " Vz: "       << setw(10) << candidate->Position.Z()
              << " T: "        << setw(10) << candidate->Position.T()
              << std::endl;
  }

  m_eventCounter++;
  readStopWatch.Stop();

  //
  // Process event
  TStopwatch procStopWatch;

  // Delphes process
  procStopWatch.Start();
  m_Delphes->ProcessTask();
  procStopWatch.Stop();

  // Generate Delphes branch: Event
  m_HepMCReader->MakeEventBranch(m_branchEvent, &readStopWatch, &procStopWatch);
  if (m_outRootFile!=nullptr) m_treeWriter->Fill();

  // FCC EDM (event-data model) based output
  fccedm::MCParticleCollection* genParticles     = new fccedm::MCParticleCollection();
  fccedm::GenVertexCollection*  genProdVertices  = new fccedm::GenVertexCollection();
  fccedm::GenVertexCollection*  genDecVertices   = new fccedm::GenVertexCollection();
  fccedm::ParticleCollection*   recMuons         = new fccedm::ParticleCollection();
  fccedm::ParticleCollection*   recElectrons     = new fccedm::ParticleCollection();
  fccedm::ParticleCollection*   recPhotons       = new fccedm::ParticleCollection();
  fccedm::JetCollection*        recJets          = new fccedm::JetCollection();

  fccedm::ParticleMCAssociationCollection*  recMuonsToMC     = new fccedm::ParticleMCAssociationCollection();
  fccedm::ParticleMCAssociationCollection*  recElectronsToMC = new fccedm::ParticleMCAssociationCollection();
  fccedm::ParticleMCAssociationCollection*  recPhotonsToMC   = new fccedm::ParticleMCAssociationCollection();
  fccedm::JetParticleAssociationCollection* recJetsToPart    = nullptr; //new fccedm::JetParticleAssociationCollection();
  //ParticleCollection* recMETs            = nullptr;
  //ParticleCollection* recHTs             = nullptr;

  // Fill FCC collections
  m_muonOutArray     = m_Delphes->ImportArray("MuonMomentumSmearing/muons");
  m_electronOutArray = m_Delphes->ImportArray("ElectronEnergySmearing/electrons");
  m_photonOutArray   = m_Delphes->ImportArray("Ecal/eflowPhotons");
  m_jetOutArray      = m_Delphes->ImportArray("JetEnergyScale/jets");

  DelphesSimulation::ConvertMCParticles(          m_allPartOutArray , genParticles, genProdVertices, genDecVertices);
  DelphesSimulation::ConvertParticles<Muon>(      m_muonOutArray    , recMuons    , genParticles, recMuonsToMC    );
  DelphesSimulation::ConvertParticles<Electron>(  m_electronOutArray, recElectrons, genParticles, recElectronsToMC);
  DelphesSimulation::ConvertPhotons(              m_photonOutArray  , recPhotons  , genParticles, recPhotonsToMC  );
  DelphesSimulation::ConvertJets(                 m_jetOutArray     , recJets);

  // Save FCC-EDM collections to FCCSw data store
  m_handleGenParticles.put(genParticles);
  m_handleGenProdVertices.put(genProdVertices);
  m_handleGenDecVertices.put(genDecVertices);

  m_handleRecMuons.put(recMuons);
  //m_handleRecMuonsToMC.put(recMuonsToMC);
  m_handleRecElectrons.put(recElectrons);
  //m_handleRecElectronsToMC.put(recElectronsToMC);
  m_handleRecPhotons.put(recPhotons);
  //m_handleRecPhotonsToMC.put(recPhotonsToMC);
  m_handleRecJets.put(recJets);
  //m_handleRecJetsToPart.put(recJetsToPart);

  // Initialize for next event reading (Will also zero Delphes arrays)
  if (m_outRootFile!=nullptr) m_treeWriter->Clear();
  m_Delphes->Clear();
  m_HepMCReader->Clear();


  return StatusCode::SUCCESS;
}

StatusCode DelphesSimulation::finalize() {

  // Finish Delphes task
  m_Delphes->FinishTask();

  // Close HepMC input file if defined
  if (m_inHepMCFile!=nullptr) {

    fclose(m_inHepMCFile);
  }

  // Write output to Root file
  if (m_outRootFile!=nullptr) {

    m_treeWriter->Write();
    m_outRootFile->Close();
    if (m_outRootFile){delete m_outRootFile; m_outRootFile = nullptr;}
  }
  
  info() << "Exiting Delphes..." << endmsg;
  
  // Clear memory
  if (m_HepMCReader) {delete m_HepMCReader; m_HepMCReader = nullptr; } // Releases also the memory allocated by inHepMCFile
  if (m_Delphes)     {delete m_Delphes;     m_Delphes     = nullptr; } // Releases also the memory allocated by treeWriter
  if (m_confReader)  {delete m_confReader;  m_confReader  = nullptr; }
  
  return GaudiAlgorithm::finalize();
}

void DelphesSimulation::ConvertMCParticles(const TObjArray* Input , fccedm::MCParticleCollection* colMCParticles, fccedm::GenVertexCollection* colProdVertices, fccedm::GenVertexCollection* colDecVertices)
{
  Candidate* cand = nullptr;
  for(int j=0; j<Input->GetEntries(); j++) {

    cand = static_cast<Candidate *>(Input->At(j));

    fccedm::MCParticleHandle particleHandle   = colMCParticles->create();
    fccedm::GenVertexHandle  prodVertexHandle = colProdVertices->create();
    fccedm::GenVertexHandle  decVertexHandle  = colDecVertices->create();

    auto particle    = particleHandle.mod().Core;
    auto prodVertex  = prodVertexHandle.mod();
    auto decayVertex = decVertexHandle.mod();

    particleHandle.mod().StartVertex = prodVertexHandle;
    particleHandle.mod().EndVertex   = decVertexHandle;

    particle.Type     = cand->PID;
    particle.Status   = cand->Status;
  	particle.P4.Px    = cand->Momentum.X();
    particle.P4.Py    = cand->Momentum.Y();
    particle.P4.Pz    = cand->Momentum.Z();
    particle.P4.Mass  = cand->Mass;
    particle.Charge   = cand->Charge;
    particle.Vertex.X = cand->Position.X();
    particle.Vertex.Y = cand->Position.Y();
    particle.Vertex.Z = cand->Position.Z();

    prodVertex.Position.X = cand->Position.X();
    prodVertex.Position.Y = cand->Position.Y();
    prodVertex.Position.Z = cand->Position.Z();
    prodVertex.Ctau       = cand->Position.T();

    // Colliding particle
    if (cand->M1==-1) {

      particle.Bits       = ParticleStatus::Beam;
      Candidate* daughter = nullptr;
      if (cand->D1!=-1) {

        daughter = static_cast<Candidate *>(Input->At(cand->D1));
        decayVertex.Position.X = daughter->Position.X();
        decayVertex.Position.Y = daughter->Position.Y();
        decayVertex.Position.Z = daughter->Position.Z();
        decayVertex.Ctau       = daughter->Position.T();
      }
      else {

        particle.Bits = ParticleStatus::Stable;
        decayVertex.Position.X = -std::numeric_limits<float>::max();
        decayVertex.Position.Y = -std::numeric_limits<float>::max();
        decayVertex.Position.Z = -std::numeric_limits<float>::max();
        decayVertex.Ctau       = -std::numeric_limits<float>::max();
      }
    }
    // Stable particle
    else if (cand->D1==-1) {

      particle.Bits = ParticleStatus::Stable;
      decayVertex.Position.X = -std::numeric_limits<float>::max();
      decayVertex.Position.Y = -std::numeric_limits<float>::max();
      decayVertex.Position.Z = -std::numeric_limits<float>::max();
      decayVertex.Ctau       = -std::numeric_limits<float>::max();
    }
    // Intermediate state particle
    else {

      particle.Bits = ParticleStatus::Decayed;
      Candidate* daughter = static_cast<Candidate *>(Input->At(cand->D1));
      decayVertex.Position.X = daughter->Position.X();
      decayVertex.Position.Y = daughter->Position.Y();
      decayVertex.Position.Z = daughter->Position.Z();
      decayVertex.Ctau       = daughter->Position.T();
    }

    /*std::cout << "Delphes FCCEDM: "
              << " Id: "       << setw(3)  << j
              << " Pdg: "      << setw(5)  << particle.Type
              << " Stat: "     << setw(2)  << particle.Status
              << " Px: "       << setw(10) << particle.P4.Px
              << " Py: "       << setw(10) << particle.P4.Py
              << " Pz: "       << setw(10) << particle.P4.Pz
              << " E: "        << setw(10) << sqrt(particle.P4.Px*particle.P4.Px + particle.P4.Py*particle.P4.Py + particle.P4.Pz*particle.P4.Pz + particle.P4.Mass*particle.P4.Mass)
              << " M: "        << setw(10) << particle.P4.Mass;
    std::cout << " Vx: "       << setw(10) << prodVertex.Position.X
              << " Vy: "       << setw(10) << prodVertex.Position.Y
              << " Vz: "       << setw(10) << prodVertex.Position.Z
              << " T: "        << setw(10) << prodVertex.Ctau;
    std::cout << std::endl;*/
  }
}   

template <class T> void DelphesSimulation::ConvertParticles(const TObjArray*  Input, fccedm::ParticleCollection* colParticles, fccedm::MCParticleCollection* colMCParticles, fccedm::ParticleMCAssociationCollection* ascColParticlesToMC)
{
  T* part = nullptr;
  for(int j=0; j<Input->GetEntries(); j++) {

    part = (T *)(Input->At(j));

    fccedm::ParticleHandle               particleHandle = colParticles->create();
    fccedm::ParticleMCAssociationHandle  relationHandle = ascColParticlesToMC->create();

    auto particle = particleHandle.mod().Core;
    auto relation = relationHandle.mod();

    particle.Charge   = part->Charge;
    particle.Status   = -1;
    particle.Type     = -1;
    particle.P4.Mass  = part->P4().M();
    particle.P4.Px    = part->P4().Px();
    particle.P4.Py    = part->P4().Py();
    particle.P4.Pz    = part->P4().Pz();
    particle.Vertex.X = -1;
    particle.Vertex.Y = -1;
    particle.Vertex.Z = -1;

    // Reference to MC - Delphes holds references to all objects related to the <T> object, only one relates to MC particle
    TObjArray* refParticles  = ((Candidate *)(Input->At(j)))->GetCandidates();
    int        nRefParticles = 0;
    int        indexMCPart   = -1;
    int        nMCParticles  = colMCParticles->getHandles().size();
    for(int k=0; k<refParticles->GetEntries(); k++) {

      int index = refParticles->At(k)->GetUniqueID()-1; // To start with zero value
      if (index<nMCParticles) {

        indexMCPart = index;
        nRefParticles++;
      }
    }

    if (nRefParticles!=1) {

      particle.Bits = ParticleStatus::Unmatched;
      std::cout << "WARNING: Can't build relation to MC particle!" << std::endl;
    }
    else {

      //std::cout << ">> Index reference: " << indexMCPart+1 << std::endl;
      relation.Rec = particleHandle;
      //relation.Sim = colMCParticles->get(indexMCPart);
    }
  }
}

void DelphesSimulation::ConvertPhotons(const TObjArray*  Input, fccedm::ParticleCollection* colParticles, fccedm::MCParticleCollection* colMCParticles, fccedm::ParticleMCAssociationCollection* ascColParticlesToMC)
{
  Photon* part = nullptr;
  for(int j=0; j<Input->GetEntries(); j++) {

    part = static_cast<Photon *>(Input->At(j));

    fccedm::ParticleHandle               particleHandle = colParticles->create();
    fccedm::ParticleMCAssociationHandle  relationHandle = ascColParticlesToMC->create();

    auto particle = particleHandle.mod().Core;
    auto relation = relationHandle.mod();

    particle.Charge   =  0;
    particle.Status   = -1;
    particle.Type     = -1;
    particle.P4.Mass  = part->P4().M();
    particle.P4.Px    = part->P4().Px();
    particle.P4.Py    = part->P4().Py();
    particle.P4.Pz    = part->P4().Pz();
    particle.Vertex.X = -1;
    particle.Vertex.Y = -1;
    particle.Vertex.Z = -1;

    // Reference to MC - Delphes holds references to all objects related to the Photon object, only one relates to MC particle
    // (relation can be a cascade of photons, only the mother of the cascade relates to the MC particle)
    TObjArray* refParticles   = ((Candidate *)(Input->At(j)))->GetCandidates();
    int        nRefParticles  = 0;
    int        indexMCPart    = -1;
    int        nMCParticles   = colMCParticles->getHandles().size();
    int        nRefInCascade  = 0;
    for(int k=0; k<refParticles->GetEntries(); k++) {

      int index = refParticles->At(k)->GetUniqueID()-1; // To start with zero value
      if (index<nMCParticles) {

        indexMCPart = index;
        nRefParticles++;
      }
      else {
        TObjArray* ref2Particles = ((Candidate *)(refParticles->At(k)))->GetCandidates();
        for(int l=0; l<ref2Particles->GetEntries(); l++) {

          int index = ref2Particles->At(l)->GetUniqueID()-1; //
          if (index<nMCParticles) {

            indexMCPart = index;
            nRefInCascade++;
          }
        }
      }
    }

    if ((nRefParticles+nRefInCascade)!=1) {

      particle.Bits = ParticleStatus::Unmatched;
      std::cout << "WARNING: Can't build relation to MC particle!" << std::endl;
    }
    else {

      //std::cout << ">> Index reference: " << indexMCPart+1 << std::endl;
      relation.Rec = particleHandle;
      if (nRefInCascade!=0) particle.Bits = ParticleStatus::MatchInCascade;
      //relation.Sim = colMCParticles->get(indexMCPart);
    }
  }
}

void DelphesSimulation::ConvertJets(const TObjArray*  Input, fccedm::JetCollection* colJets)
{
  Jet* jet;
  for(int j = 0; j < Input->GetEntries(); ++j) {
      
    jet = static_cast<Jet *>(Input->At(j));

    fccedm::JetHandle jetHandle = colJets->create();

    fccedm::BareJet& core = jetHandle.mod().Core;
    core.Area     = 0;
    core.P4.Px    = (double) jet->P4().X();
    core.P4.Py    = (double) jet->P4().Y();
    core.P4.Pz    = (double) jet->P4().Z();
    core.P4.Mass  = (double) jet->Mass ;
  }
}   

/*void DelphesSimulation::ConvertMET(   TObjArray *  Input , ParticleCollection *  coll  ){

  Candidate * cand;
  for(int j = 0; j < Input->GetEntries(); ++j) {
      
    cand = static_cast<Candidate *>(Input->At(j));
    ParticleHandle outptc = coll->create();
    BareParticle& core = outptc.mod().Core;
    core.P4.Px         = (double) cand->Momentum.X();
    core.P4.Py         = (double) cand->Momentum.Y();
  }
}   

void DelphesSimulation::ConvertHT(   TObjArray *  Input , ParticleCollection *  coll  ){

  Candidate * cand;
  for(int j = 0; j < Input->GetEntries(); ++j) {

    cand = static_cast<Candidate *>(Input->At(j));
    ParticleHandle outptc = coll->create();
    BareParticle& core = outptc.mod().Core;
    core.P4.Px         = (double) cand->Momentum.X();
    core.P4.Py         = (double) cand->Momentum.Y();
  }
}*/
