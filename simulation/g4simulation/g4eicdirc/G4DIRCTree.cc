#include "G4DIRCTree.h"

#include "PrtHit.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>  // for PHG4HitContainer, PHG4Hit...
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phparameter/PHParameters.h>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <cmath>     // for atan2, sqrt
#include <cstring>   // for strcmp
#include <iostream>  // for ostringstream, operator<<
#include <sstream>
#include <utility>  // for pair

G4DIRCTree::G4DIRCTree(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , _filename(filename)
  , g4tree(nullptr)
  , outfile(nullptr)
{
}

int G4DIRCTree::Init(PHCompositeNode *)
{
  outfile = new TFile(_filename.c_str(), "RECREATE");
  g4tree = new TTree("mG4EvtTree", "g4tree");
  g4tree->SetAutoSave(1000000);

  std::cout << "Initialize Geant 4 Tree ... << " << std::endl;

  /// Event level
  g4tree->Branch("momentum", &mG4EvtTree.momentum, "p/D");
  g4tree->Branch("theta", &mG4EvtTree.theta, "theta/D");
  g4tree->Branch("phi", &mG4EvtTree.phi, "phi/D");
  g4tree->Branch("px", &mG4EvtTree.px, "px/D");
  g4tree->Branch("py", &mG4EvtTree.py, "py/D");
  g4tree->Branch("pz", &mG4EvtTree.pz, "pz/D");
  g4tree->Branch("pid", &mG4EvtTree.pid, "pid/I");

  /// Hit level
  g4tree->Branch("nhits", &mG4EvtTree.nhits, "nhits/I");
  g4tree->Branch("detid", mG4EvtTree.detid, "detid[nhits]/I");
  g4tree->Branch("hitid", mG4EvtTree.hitid, "hitid[nhits]/I");
  g4tree->Branch("trkid", mG4EvtTree.trkid, "trkid[nhits]/I");
  g4tree->Branch("x0", mG4EvtTree.x0, "x0[nhits]/D");
  g4tree->Branch("y0", mG4EvtTree.y0, "y0[nhits]/D");
  g4tree->Branch("z0", mG4EvtTree.z0, "z0[nhits]/D");
  g4tree->Branch("x1", mG4EvtTree.x1, "x1[nhits]/D");
  g4tree->Branch("y1", mG4EvtTree.y1, "y1[nhits]/D");
  g4tree->Branch("z1", mG4EvtTree.z1, "z1[nhits]/D");
  g4tree->Branch("edep", mG4EvtTree.edep, "edep[nhits]/D");

  g4tree->Branch("mcp_id", mG4EvtTree.mcp_id, "mcp_id[nhits]/I");
  g4tree->Branch("pixel_id", mG4EvtTree.pixel_id, "pixel_id[nhits]/I");
  g4tree->Branch("lead_time", mG4EvtTree.lead_time,"lead_time[nhits]/D");
  g4tree->Branch("wavelength", mG4EvtTree.wavelength,"wavelength[nhits]/D");
  g4tree->Branch("hit_pathId", mG4EvtTree.hit_pathId, "hit_pathId[nhits]/L");
  g4tree->Branch("nrefl", mG4EvtTree.nrefl, "nrefl[nhits]/I");
  g4tree->Branch("parent_pid", mG4EvtTree.parent_pid, "parent_pid[nhits]/I");
  g4tree->Branch("parent_momentum", mG4EvtTree.parent_momentum, "parent_momentum[nhits]/D");

  g4tree->Branch("hit_globalPos", mG4EvtTree.hit_globalPos, "hit_globalPos[nhits][3]/D");
  g4tree->Branch("hit_localPos", mG4EvtTree.hit_localPos, "hit_localPos[nhits][3]/D");
  g4tree->Branch("hit_digiPos", mG4EvtTree.hit_digiPos, "hit_digiPos[nhits][3]/D");
  g4tree->Branch("hit_mom", mG4EvtTree.hit_mom, "hit_mom[nhits][3]/D");
  g4tree->Branch("hit_pos", mG4EvtTree.hit_pos, "hit_pos[nhits][3]/D");

  return 0;
}

int G4DIRCTree::process_event(PHCompositeNode *topNode)
{
  // get the primary particle which did this to us....
  PHG4TruthInfoContainer *truthInfoList = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  const PHG4TruthInfoContainer::Range primRange = truthInfoList->GetPrimaryParticleRange();
  double px = primRange.first->second->get_px();
  double py = primRange.first->second->get_py();
  double pz = primRange.first->second->get_pz();
  double p = sqrt(px * px + py * py + pz * pz);
  double pt = sqrt(px * px + py * py);
  double phi = atan2(py, px);
  double theta = atan2(pt, pz);
  double pid = primRange.first->second->get_pid();
  TVector3 dir_vec(px, py, pz);

  mG4EvtTree.momentum = p;
  mG4EvtTree.theta = theta;
  mG4EvtTree.phi = phi;
  mG4EvtTree.px = px;
  mG4EvtTree.py = py;
  mG4EvtTree.pz = pz;
  mG4EvtTree.pid = pid;

  int nhits = 0;

  std::ostringstream nodename;
  std::set<std::string>::const_iterator iter;

  for (iter = _node_postfix.begin(); iter != _node_postfix.end(); ++iter)
  {
    int detid = (_detid.find(*iter))->second;
    nodename.str("");
    nodename << "G4HIT_" << *iter;
    PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());

    if (!strcmp("G4HIT_hpDIRC", nodename.str().c_str()))  // DIRC
    {
      process_hit(hits, "G4HIT_hpDIRC", detid, nhits, dir_vec);
    }
  }

  mG4EvtTree.nhits = nhits;
  //std::cout << "number of hits = " << nhits << std::endl;
  if (g4tree) g4tree->Fill();

  evt_num++;

  return 0;
}

int G4DIRCTree::End(PHCompositeNode *topNode)
{
  outfile->cd();
  g4tree->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  return 0;
}

void G4DIRCTree::AddNode(const std::string &name, const int detid)
{
  _node_postfix.insert(name);
  _detid[name] = detid;
  return;
}

int G4DIRCTree::process_hit(PHG4HitContainer *hits, const std::string &dName, int detid, int &nhits, TVector3 dir_vec)
{
  if (hits)
  {
    PHG4HitContainer::ConstRange hit_range_0 = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter_0 = hit_range_0.first; hit_iter_0 != hit_range_0.second; hit_iter_0++)
    {
      PrtHit *dirc_hit = dynamic_cast<PrtHit *>(hit_iter_0->second);

      mG4EvtTree.detid[nhits] = detid;
      mG4EvtTree.hitid[nhits] = dirc_hit->get_hit_id();
      mG4EvtTree.trkid[nhits] = dirc_hit->get_trkid();

      mG4EvtTree.x0[nhits] = dirc_hit->get_x(0);
      mG4EvtTree.y0[nhits] = dirc_hit->get_y(0);
      mG4EvtTree.z0[nhits] = dirc_hit->get_z(0);
      mG4EvtTree.x1[nhits] = dirc_hit->get_x(1);
      mG4EvtTree.y1[nhits] = dirc_hit->get_y(1);
      mG4EvtTree.z1[nhits] = dirc_hit->get_z(1);
      mG4EvtTree.edep[nhits] = dirc_hit->get_edep();

      mG4EvtTree.mcp_id[nhits] = dirc_hit->GetMcpId();
      mG4EvtTree.pixel_id[nhits] = dirc_hit->GetPixelId();
      mG4EvtTree.lead_time[nhits] = dirc_hit->GetLeadTime();
      mG4EvtTree.wavelength[nhits] = dirc_hit->GetTotTime();
      mG4EvtTree.hit_pathId[nhits] = dirc_hit->GetPathInPrizm();
      mG4EvtTree.nrefl[nhits] = dirc_hit->GetNreflectionsInPrizm();
      mG4EvtTree.parent_pid[nhits] = dirc_hit->GetParentParticleId();
      mG4EvtTree.parent_momentum[nhits] = dirc_hit->GetParentParticleMomentum().Mag();

      for (int i = 0; i < 3; i++)
      {
        mG4EvtTree.hit_globalPos[nhits][i] = dirc_hit->GetGlobalPos()(i);
	mG4EvtTree.hit_localPos[nhits][i] = dirc_hit->GetLocalPos()(i);
	mG4EvtTree.hit_digiPos[nhits][i] = dirc_hit->GetDigiPos()(i);
	mG4EvtTree.hit_mom[nhits][i] = dirc_hit->GetMomentum()(i);
	mG4EvtTree.hit_pos[nhits][i] = dirc_hit->GetPosition()(i);	    
	  
      }

      nhits++;
    }

  }

  return 0;
}
