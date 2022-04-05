#include "G4EicDircSteppingAction.h"
#include "G4EicDircDetector.h"
#include "PrtHit.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4StepStatusDecode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <TSystem.h>
#include <TVector3.h>

#include <Geant4/G4NavigationHistory.hh>
#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>  // for G4TouchableHandle
#include <Geant4/G4Track.hh>            // for G4Track
#include <Geant4/G4TrackStatus.hh>      // for fStopAndKill
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation
#include <Geant4/Randomize.hh>

#include <cmath>  // for isfinite
#include <iostream>
#include <string>  // for operator<<, string

class PHCompositeNode;
//____________________________________________________________________________..
G4EicDircSteppingAction::G4EicDircSteppingAction(
    G4EicDircDetector *detector, const PHParameters *parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_BlackHoleFlag(m_Params->get_int_param("blackhole"))
{
}

G4EicDircSteppingAction::~G4EicDircSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
  delete m_PrtHit;
}

//____________________________________________________________________________..
bool G4EicDircSteppingAction::UserSteppingAction(const G4Step *aStep,
                                                 bool was_used)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume *volume = touch->GetVolume();
  G4VPhysicalVolume *volume_post = touchpost->GetVolume();
  // IsInDetector(volume) returns
  //  == 0 outside of detector
  //   > 0 for hits in active volume
  //  < 0 for hits in passive material
  int whichactive_int = m_Detector->IsInDetector(volume);
  int whichactive_int_post = m_Detector->IsInDetector(volume_post);
  bool whichactive = (whichactive_int > 0 && whichactive_int < 12);

  if (!whichactive)
  {
    return false;
  }

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion =
      (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) /
      GeV;
  const G4Track *aTrack = aStep->GetTrack();

  if (aTrack->GetCurrentStepNumber() > 50000 || aTrack->GetTrackLength() > 30000)
  {
    //std::cout<<"WRN: too many steps or track length > 30 m  N=" << aTrack->GetCurrentStepNumber()<<" Len = "<< aTrack->GetTrackLength()/1000. <<std::endl;
    return false;
  }

  if (whichactive_int == 10 && volume_post->GetName() == "World") return false;  // 10-Mcp
  if (whichactive_int == 1 && volume_post->GetName() == "World") return false;   // 1-PMT plane

  /*if ((whichactive_int == 6 || whichactive_int == 7) && (volume_post->GetName()=="World")) return false; // 6-lens1, 7-lens2

  if(whichactive_int == 6 && whichactive_int_post == 6) return false;
  if(whichactive_int == 7 && whichactive_int_post == 7) return false;

  if (whichactive_int == 8 && (aStep->GetPostStepPoint()->GetPosition().z() > aStep->GetPreStepPoint()->GetPosition().z()))
    {
      // 8 - lens layer 3
      return false;
    }
  */

  // if this block stops everything, just put all kinetic energy into edep
  if (m_BlackHoleFlag)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track *killtrack = const_cast<G4Track *>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  // make sure we are in a volume
  if (m_ActiveFlag)
  {
    bool geantino = false;
      // the check for the pdg code speeds things up, I do not want to make
      // an expensive string compare for every track when we know
      // geantino or chargedgeantino has pid=0
      if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
          aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != std::string::npos)
      {
        geantino = true;
      }

      const G4DynamicParticle *dynParticle = aTrack->GetDynamicParticle();                                                    
      G4ParticleDefinition *particle = dynParticle->GetDefinition();                                                          
      G4String ParticleName = particle->GetParticleName();                                                                    
         
      G4StepPoint *prePoint = aStep->GetPreStepPoint();
      G4StepPoint *postPoint = aStep->GetPostStepPoint();
      //       std::cout << "track id " << aTrack->GetTrackID() << std::endl;
      //       std::cout << "time prepoint: " << prePoint->GetGlobalTime() << std::endl;
      //       std::cout << "time postpoint: " << postPoint->GetGlobalTime() << std::endl;
      int prepointstatus = prePoint->GetStepStatus();
      if (prepointstatus == fGeomBoundary ||
          prepointstatus == fUndefined ||
          (prepointstatus == fPostStepDoItProc && m_SavePostStepStatus == fGeomBoundary))
      {
        if (prepointstatus == fPostStepDoItProc && m_SavePostStepStatus == fGeomBoundary)
        {
          std::cout << GetName() << ": New Hit for  " << std::endl;
          std::cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
                    << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
                    << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
                    << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
          std::cout << "last track: " << m_SaveTrackId
                    << ", current trackid: " << aTrack->GetTrackID() << std::endl;
          std::cout << "phys pre vol: " << volume->GetName()
                    << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
          std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
                    << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
        }
        if (!m_Hit)
        {
          m_Hit = new PHG4Hitv1();
        }
    
	
	if((whichactive_int == 2 || whichactive_int == 3) && (prepointstatus == fGeomBoundary) && (aTrack->GetParentID() == 0) && (ParticleName != "opticalphoton")) //pre step in a long piece of bar(2) or short piece of bar(3) ------
	  {
	//here we set the entrance values in cm
        m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
        m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
        m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

        m_Hit->set_px(0, prePoint->GetMomentum().x() / GeV);
        m_Hit->set_py(0, prePoint->GetMomentum().y() / GeV);
        m_Hit->set_pz(0, prePoint->GetMomentum().z() / GeV);
        // time in ns
        m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
        //set the track ID
        m_Hit->set_trkid(aTrack->GetTrackID());
	  }
        m_SaveTrackId = aTrack->GetTrackID();

        //set the initial energy deposit
        m_Hit->set_edep(0);
        if (!geantino && !m_BlackHoleFlag && m_ActiveFlag)
        {
          m_Hit->set_eion(0);
        }
        if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
        {
          if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
          {
            m_Hit->set_trkid(pp->GetUserTrackId());
            m_Hit->set_shower_id(pp->GetShower()->get_id());
          }
        }
	  
      }

      // some sanity checks for inconsistencies
      // check if this hit was created, if not print out last post step status
      if (!m_Hit || !isfinite(m_Hit->get_x(0)))
      {
        std::cout << GetName() << ": hit was not created" << std::endl;
        std::cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
                  << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
                  << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
                  << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
        std::cout << "last track: " << m_SaveTrackId
                  << ", current trackid: " << aTrack->GetTrackID() << std::endl;
        std::cout << "phys pre vol: " << volume->GetName()
                  << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
        std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
                  << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
        exit(1);
      }
      m_SavePostStepStatus = postPoint->GetStepStatus();
      // check if track id matches the initial one when the hit was created
      if (aTrack->GetTrackID() != m_SaveTrackId)
      {
        std::cout << "hits do not belong to the same track" << std::endl;
        std::cout << "saved track: " << m_SaveTrackId
                  << ", current trackid: " << aTrack->GetTrackID()
                  << std::endl;
        exit(1);
      }
      m_SavePreStepStatus = prePoint->GetStepStatus();
      m_SavePostStepStatus = postPoint->GetStepStatus();
      m_SaveVolPre = volume;
      m_SaveVolPost = touchpost->GetVolume();

      // here we just update the exit values, it will be overwritten
      // for every step until we leave the volume or the particle
      // ceases to exist
      
      if ((whichactive_int_post == 2 || whichactive_int_post == 3) && (postPoint->GetStepStatus() == fGeomBoundary) && (aTrack->GetParentID() == 0) && (ParticleName != "opticalphoton")) //post step in a long piece of bar(2) or short piece of bar(3) ------ 
	{
      
	  m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
	  m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
	  m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

	  m_Hit->set_px(1, postPoint->GetMomentum().x() / GeV);
	  m_Hit->set_py(1, postPoint->GetMomentum().y() / GeV);
	  m_Hit->set_pz(1, postPoint->GetMomentum().z() / GeV);

	  m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
	  //sum up the energy to get total deposited
	  m_Hit->set_edep(m_Hit->get_edep() + edep);
	

      // update ionization energy only for active volumes, not for black holes or geantinos
      // if the hit is created without eion, get_eion() returns NAN
      // if you see NANs check the creation of the hit
	  if (m_ActiveFlag && !m_BlackHoleFlag && !geantino)
	    {
	      m_Hit->set_eion(m_Hit->get_eion() + eion);
	    }
	  if (geantino)
	    {
	      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
	    }
	  if (edep > 0)
	    {
	      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
		{
		  if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
		    {
		      pp->SetKeep(1);  // we want to keep the track
		    }
		}
	    }
      // if any of these conditions is true this is the last step in
      // this volume and we need to save the hit
      // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
      // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
      // (happens when your detector goes outside world volume)
      // postPoint->GetStepStatus() == fAtRestDoItProc: track stops (typically
      // aTrack->GetTrackStatus() == fStopAndKill is also set)
      // aTrack->GetTrackStatus() == fStopAndKill: track ends
      if (postPoint->GetStepStatus() == fGeomBoundary ||
          postPoint->GetStepStatus() == fWorldBoundary ||
          postPoint->GetStepStatus() == fAtRestDoItProc ||
          aTrack->GetTrackStatus() == fStopAndKill)
      {
        // save only hits with energy deposit (or -1 for geantino)
        if (m_Hit->get_edep())
        {
          m_AbsorberHitContainer->AddHit(0, m_Hit);
          // ownership has been transferred to container, set to null
          // so we will create a new hit for the next track
          m_Hit = nullptr;
        }
        else
        {
          // if this hit has no energy deposit, just reset it for reuse
          // this means we have to delete it in the dtor. If this was
          // the last hit we processed the memory is still allocated
          m_Hit->Reset();
        }
      }
	} // end save track information	
        
    
    //-----------------------------

      int detector_id = 0;  // we use here only one detector in this simple example
      
      G4String prePointVolName = prePoint->GetPhysicalVolume()->GetName();
      G4String postPointVolName = postPoint->GetPhysicalVolume()->GetName();
      
      if (!m_PrtHit)
        {
          m_PrtHit = new PrtHit();
        }

      
      if (postPoint->GetStepStatus() == fGeomBoundary ||
          postPoint->GetStepStatus() == fWorldBoundary ||
          postPoint->GetStepStatus() == fAtRestDoItProc ||
          aTrack->GetTrackStatus() == fStopAndKill)
      {
	if (whichactive_int_post == 11)  // post step in Pixel ---------------
          {
            G4ThreeVector globalpos = aStep->GetPostStepPoint()->GetPosition();
            G4ThreeVector localpos = touchpost->GetHistory()->GetTopTransform().TransformPoint(globalpos);
            G4ThreeVector translation = touchpost->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0, 0, 0));
            G4ThreeVector inPrismpos = touchpost->GetHistory()->GetTransform(1).TransformPoint(globalpos);
            G4ThreeVector g4mom = aTrack->GetVertexMomentumDirection();  //GetMomentum();
            G4ThreeVector g4pos = aTrack->GetVertexPosition();

            G4ThreeVector localvec = touchpost->GetHistory()->GetTopTransform().TransformAxis(g4mom);

            TVector3 globalPos(inPrismpos.x(), inPrismpos.y(), inPrismpos.z());
            TVector3 localPos(localpos.x(), localpos.y(), localpos.z());
            TVector3 digiPos(translation.x(), translation.y(), translation.z());
            TVector3 momentum(g4mom.x(), g4mom.y(), g4mom.z());
            TVector3 position(g4pos.x(), g4pos.y(), g4pos.z());

            double time = aStep->GetPreStepPoint()->GetLocalTime();

            int mcp = touchpost->GetReplicaNumber(1);
            int pix = touchpost->GetReplicaNumber(0);

            Double_t wavelength = 1.2398 / (aTrack->GetMomentum().mag() * 1E6) * 1000;

            // transport efficiency ----

            double pi(4 * atan(1));
            double roughness(0.5);  // nm
            double angleX = localvec.angle(G4ThreeVector(1, 0, 0));
            double angleY = localvec.angle(G4ThreeVector(0, 1, 0));
            if (angleX > 0.5 * pi) angleX = pi - angleX;
            if (angleY > 0.5 * pi) angleY = pi - angleY;
            double length = aTrack->GetTrackLength() - 400;  // 400 - average path in EV
            double lengthx = fabs(length * localvec.x());    // along the bar
            double lengthy = fabs(length * localvec.y());

            //cout << "track length = " << aTrack->GetTrackLength() << endl;

            int nBouncesX = (int) (lengthx) / 17.25;
            int nBouncesY = (int) (lengthy) / 35.0;

            double ll = wavelength * wavelength;
            double n_quartz = sqrt(1. + (0.696 * ll / (ll - pow(0.068, 2))) +
                                   (0.407 * ll / (ll - pow(0.116, 2))) + 0.897 * ll / (ll - pow(9.896, 2)));
            double bounce_probX = 1 - pow(4 * pi * cos(angleX) * roughness * n_quartz / wavelength, 2);
            double bounce_probY = 1 - pow(4 * pi * cos(angleY) * roughness * n_quartz / wavelength, 2);

            double totalProb = pow(bounce_probX, nBouncesX) * pow(bounce_probY, nBouncesY);

            if (G4UniformRand() > totalProb) return false;

            // time since track created
            m_PrtHit->SetLeadTime(time);
            m_PrtHit->SetTotTime(wavelength);  //set photon wavelength
            m_PrtHit->SetMcpId(mcp);
            m_PrtHit->SetPixelId(pix);
            m_PrtHit->SetGlobalPos(globalPos);
            m_PrtHit->SetLocalPos(localPos);
            m_PrtHit->SetDigiPos(digiPos);
            m_PrtHit->SetPosition(position);
            m_PrtHit->SetMomentum(momentum);

            
            m_SaveHitContainer->AddHit(detector_id, m_PrtHit);

            // ownership has been transferred to container, set to null
            // so we will create a new hit for the next track
            //m_Hit = nullptr;
            m_PrtHit = nullptr;
          }
      }

        
      // return true to indicate the hit was used
      return true;
  }

   return false;
}

//____________________________________________________________________________..
void G4EicDircSteppingAction::SetInterfacePointers(PHCompositeNode *topNode)
{
  if (!m_HitNodeName.empty())
  {
    m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  }
  if (!m_AbsorberNodeName.empty())
  {
    m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_AbsorberNodeName);
    if (!m_AbsorberHitContainer)
    {
      if (Verbosity() > 0)
      {
        std::cout << "G4EicDircSteppingAction::SetTopNode - unable to find " << m_AbsorberNodeName << std::endl;
      }
    }
  }
  if (!m_SupportNodeName.empty())
  {
    m_SupportHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_SupportNodeName);
    if (!m_SupportHitContainer)
    {
      if (Verbosity() > 0)
      {
        std::cout << "G4EicDircSteppingAction::SetTopNode - unable to find " << m_SupportNodeName << std::endl;
      }
    }
  }
  if (!m_HitContainer)
  {
    std::cout << "G4EicDircSteppingAction::SetTopNode - unable to find " << m_HitNodeName << std::endl;
  }
}
