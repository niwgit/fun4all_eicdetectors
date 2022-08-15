#include "G4EicDircOpBoundaryProcess.h"

#include <Geant4/G4Step.hh>
#include <Geant4/G4TouchableHistory.hh>
#include <Geant4/G4Track.hh>

#include <set>

G4VParticleChange* G4EicDircOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  G4StepPoint* pPreStepPoint = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4VParticleChange* pParticleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep);
  
  double endofbar = -(1)*(0.5*(4235 + 4*0.05)) - 437.5; 
    
  if(pPostStepPoint->GetPosition().z() > pPreStepPoint->GetPosition().z())
    {
      pParticleChange->ProposeTrackStatus(fStopAndKill);
      if(pPreStepPoint->GetPosition().z() > endofbar) pParticleChange->ProposeTrackStatus(fStopAndKill);
    }

 if ((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName().contains("lLens3")) && pPostStepPoint->GetPosition().z() > pPreStepPoint->GetPosition().z())
  {
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }

  // kill photons outside bar and prizm

  if (GetStatus() == FresnelRefraction && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World")
  {
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }

  if (((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName().contains("lLens1")) || (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName().contains("lLens2"))) && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World")
  {
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }

  // // black edge of the lens3
  // if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
  //     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc")
  //    || (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3"
  // 	 &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wLens3")){
  //   pParticleChange->ProposeTrackStatus(fStopAndKill);
  // }

  if (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName().contains("lLens1") && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName().contains("lLens1"))
  {
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }
  if (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName().contains("lLens2") && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName().contains("lLens2"))
  {
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }

  return pParticleChange;
}
