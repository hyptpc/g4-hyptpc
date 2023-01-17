#include "SACEventAction.hh"
#include "SACAnalysisManager.hh"

SACEventAction::SACEventAction(SACAnalysisManager *analysisManager)
  : G4UserEventAction(), anaMan(analysisManager)
{
}

SACEventAction::~SACEventAction()
{
}

void SACEventAction::BeginOfEventAction(const G4Event *anEvent)
{
  if (anaMan) anaMan->BeginOfEvent(anEvent);
}

void SACEventAction::EndOfEventAction(const G4Event *anEvent)
{
  if (anaMan) anaMan->EndOfEvent(anEvent);
}
