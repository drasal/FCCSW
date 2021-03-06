// local
#include "SimG4ConstantMagneticFieldTool.h"

// FCCSW
#include "SimG4Common/ConstantField.h"

// Geant 4
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"

#include "G4PropagatorInField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4NystromRK4.hh"
#include "G4ClassicalRK4.hh"

#include "G4SystemOfUnits.hh"

// Declaration of the Tool
DECLARE_COMPONENT(SimG4ConstantMagneticFieldTool)

SimG4ConstantMagneticFieldTool::SimG4ConstantMagneticFieldTool(const std::string& type, const std::string& name,
                                                         const IInterface* parent)
    : GaudiTool(type, name, parent),
      m_field(nullptr),
      m_fieldOn(false),
      m_minEps(0),
      m_maxEps(0),
      m_deltaChord(0),
      m_deltaOneStep(0),
      m_maxStep(1. * m),
      m_integratorStepper("NystromRK4"),
      m_fieldComponentX(0),
      m_fieldComponentY(0),
      m_fieldComponentZ(-6. * tesla),
      m_fieldRadMax(6 * m),
      m_fieldZMax(6. * m) {
  declareInterface<ISimG4MagneticFieldTool>(this);

  declareProperty("FieldOn", m_fieldOn, "Switch to turn field off");
  declareProperty("DeltaChord", m_deltaChord, "Missing distance for the chord finder");
  declareProperty("DeltaOneStep", m_deltaOneStep, "Delta(one-step)");
  declareProperty("MinimumEpsilon", m_minEps, "Minimum epsilon (see G4 documentation)");
  declareProperty("MaximumEpsilon", m_maxEps, "Maximum epsilon (see G4 documentation)");
  declareProperty("MaximumStep", m_maxStep, "Maximum step length in field (see G4 documentation)");

  declareProperty("IntegratorStepper", m_integratorStepper = "NystromRK4", "Integrator stepper name");
  declareProperty("FieldComponentX", m_fieldComponentX, "Field X component");
  declareProperty("FieldComponentY", m_fieldComponentY, "Field Y component");
  declareProperty("FieldComponentZ", m_fieldComponentZ, "Field Z component");
  declareProperty("FieldRMax", m_fieldRadMax, "Field max radius");
  declareProperty("FieldZMax", m_fieldZMax, "field max Z");
}

SimG4ConstantMagneticFieldTool::~SimG4ConstantMagneticFieldTool() {
  if (nullptr != m_field) delete m_field;
}

StatusCode SimG4ConstantMagneticFieldTool::initialize() {
  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) return sc;

  if (m_fieldOn) {
    G4TransportationManager* transpManager = G4TransportationManager::GetTransportationManager();
    G4FieldManager* fieldManager = transpManager->GetFieldManager();
    G4PropagatorInField* propagator = transpManager->GetPropagatorInField();

    // The field manager keeps an observing pointer to the field, ownership stays with this tool. (Cleaned up in dtor)
    m_field =
        new sim::ConstantField(m_fieldComponentX, m_fieldComponentY, m_fieldComponentZ, m_fieldRadMax, m_fieldZMax);
    fieldManager->SetDetectorField(m_field);

    fieldManager->CreateChordFinder(m_field);
    G4ChordFinder* chordFinder = fieldManager->GetChordFinder();
    chordFinder->GetIntegrationDriver()->RenewStepperAndAdjust(stepper(m_integratorStepper, m_field));

    propagator->SetLargestAcceptableStep(m_maxStep);

    if (m_deltaChord > 0) fieldManager->GetChordFinder()->SetDeltaChord(m_deltaChord);
    if (m_deltaOneStep > 0) fieldManager->SetDeltaOneStep(m_deltaOneStep);
    if (m_minEps > 0) fieldManager->SetMinimumEpsilonStep(m_minEps);
    if (m_maxEps > 0) fieldManager->SetMaximumEpsilonStep(m_maxEps);
  }
  return sc;
}

StatusCode SimG4ConstantMagneticFieldTool::finalize() {
  StatusCode sc = GaudiTool::finalize();
  return sc;
}

const G4MagneticField* SimG4ConstantMagneticFieldTool::field() const { return m_field; }

G4MagIntegratorStepper* SimG4ConstantMagneticFieldTool::stepper(const std::string& name, G4MagneticField* field) const {
  G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(field);
  if (name == "HelixImplicitEuler")
    return new G4HelixImplicitEuler(fEquation);
  else if (name == "HelixSimpleRunge")
    return new G4HelixSimpleRunge(fEquation);
  else if (name == "HelixExplicitEuler")
    return new G4HelixExplicitEuler(fEquation);
  else if (name == "NystromRK4")
    return new G4NystromRK4(fEquation);
  else if (name == "ClassicalRK4")
    return new G4ClassicalRK4(fEquation);
  else {
    error() << "Stepper " << name << " not available! returning NystromRK4!" << endmsg;
    return new G4NystromRK4(fEquation);
  }
}
