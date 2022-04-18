#include "G4DIRCTree.h"

#include "PrtHit.h"
//#include "PrtLutNode.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>  // for PHG4HitContainer, PHG4Hit...
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>

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
  , nblocks(0)
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

  //g4tree->Branch("track_mom_bar", mG4EvtTree.track_mom_bar, "track_mom_bar[3]/D");
  //g4tree->Branch("track_hit_pos_bar", mG4EvtTree.track_hit_pos_bar, "track_hit_pos_bar[3]/D");

  /// Hit level
  g4tree->Branch("nhits", &mG4EvtTree.nhits, "nhits/I");
  g4tree->Branch("nbarhits", &mG4EvtTree.nbarhits, "nbarhits/I");
  g4tree->Branch("detid", mG4EvtTree.detid, "detid[nhits]/I");
  g4tree->Branch("hitid", mG4EvtTree.hitid, "hitid[nhits]/I");
  g4tree->Branch("trkid", mG4EvtTree.trkid, "trkid[nhits]/I");
  /*g4tree->Branch("x0", mG4EvtTree.x0, "x0[nhits]/D");
  g4tree->Branch("y0", mG4EvtTree.y0, "y0[nhits]/D");
  g4tree->Branch("z0", mG4EvtTree.z0, "z0[nhits]/D");
  g4tree->Branch("x1", mG4EvtTree.x1, "x1[nhits]/D");
  g4tree->Branch("y1", mG4EvtTree.y1, "y1[nhits]/D");
  g4tree->Branch("z1", mG4EvtTree.z1, "z1[nhits]/D");*/
  g4tree->Branch("x0", mG4EvtTree.x0, "x0[nbarhits]/D");                                                                       
  g4tree->Branch("y0", mG4EvtTree.y0, "y0[nbarhits]/D");                                                                        
  g4tree->Branch("z0", mG4EvtTree.z0, "z0[nbarhits]/D");                                                                        
  g4tree->Branch("x1", mG4EvtTree.x1, "x1[nbarhits]/D");                                                                        
  g4tree->Branch("y1", mG4EvtTree.y1, "y1[nbarhits]/D");                                                                        
  g4tree->Branch("z1", mG4EvtTree.z1, "z1[nbarhits]/D");

  g4tree->Branch("edep", mG4EvtTree.edep, "edep[nhits]/D");

  /*g4tree->Branch("track_angle_at_bar", mG4EvtTree.track_angle_at_bar, "track_angle_at_bar[nhits]/D");
  g4tree->Branch("track_phi", mG4EvtTree.track_phi, "track_phi[nhits]/D");
  g4tree->Branch("bar_hit_time", mG4EvtTree.bar_hit_time, "bar_hit_time[nhits]/D");
  g4tree->Branch("track_mom", mG4EvtTree.track_mom, "track_mom[nhits][3]/D");
  g4tree->Branch("track_pos", mG4EvtTree.track_pos, "track_pos[nhits][3]/D");*/

  g4tree->Branch("mcp_id", mG4EvtTree.mcp_id, "mcp_id[nhits]/I");
  g4tree->Branch("pixel_id", mG4EvtTree.pixel_id, "pixel_id[nhits]/I");
  g4tree->Branch("lead_time", mG4EvtTree.lead_time,"lead_time[nhits]/D");
  g4tree->Branch("wavelength", mG4EvtTree.wavelength,"wavelength[nhits]/D");
  //g4tree->Branch("hit_pathId", mG4EvtTree.hit_pathId, "hit_pathId[nhits]/L");
  //g4tree->Branch("nrefl", mG4EvtTree.nrefl, "nrefl[nhits]/I");

  g4tree->Branch("hit_globalPos", mG4EvtTree.hit_globalPos, "hit_globalPos[nhits][3]/D");
  g4tree->Branch("hit_localPos", mG4EvtTree.hit_localPos, "hit_localPos[nhits][3]/D");
  g4tree->Branch("hit_digiPos", mG4EvtTree.hit_digiPos, "hit_digiPos[nhits][3]/D");
  g4tree->Branch("hit_mom", mG4EvtTree.hit_mom, "hit_mom[nhits][3]/D");
  g4tree->Branch("hit_pos", mG4EvtTree.hit_pos, "hit_pos[nhits][3]/D");
  //g4tree->Branch("track_mom_bar", mG4EvtTree.track_mom_bar, "track_mom_bar[nhits][3]/D");
  //g4tree->Branch("track_hit_pos_bar", mG4EvtTree.track_hit_pos_bar, "track_hit_pos_bar[nhits][3]/D");

  //-------- LUT --------
  
  /*fLut = new TClonesArray("PrtLutNode");
  fLutTree = new TTree("prtlut","Look-up table for the geometrical reconstruction.");
  fLutTree->Branch("LUT",&fLut,256000,2); 
  Int_t Nnodes = 100000;
    
  TClonesArray &fLuta = *fLut; 
  for (Long64_t n=0; n<Nnodes; n++) {
    new((fLuta)[n]) PrtLutNode(n);
    } */   

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
  int nbarhits = 0;

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

    if (!strcmp("G4HIT_ABSORBER_hpDIRC", nodename.str().c_str()))  // BAR
      {
	process_track_bar_hit(hits, "G4HIT_ABSORBER_hpDIRC", detid, nbarhits);
      }
			      
  }

  mG4EvtTree.nhits = nhits;
  mG4EvtTree.nbarhits = nbarhits;

  if (g4tree) g4tree->Fill();
  //if (fLutTree) fLutTree->Fill();

  evt_num++;

  return 0;
}

int G4DIRCTree::End(PHCompositeNode *topNode)
{
  outfile->cd();
  g4tree->Write();
  //fLutTree->Write();
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


int G4DIRCTree::process_track_bar_hit(PHG4HitContainer *hits, const std::string &dName, int detid, int &nbarhits)
{
  if (hits)
    {
      PHG4HitContainer::ConstRange bar_hit_range = hits->getHits();
      for (PHG4HitContainer::ConstIterator bar_hit_iter = bar_hit_range.first; bar_hit_iter != bar_hit_range.second; bar_hit_iter++)
	{
	  PHG4Hit* bar_Hit = dynamic_cast<PHG4Hit *>(bar_hit_iter->second);

	  mG4EvtTree.x0[nbarhits] = bar_Hit->get_x(0);
	  mG4EvtTree.y0[nbarhits] = bar_Hit->get_y(0);
	  mG4EvtTree.z0[nbarhits] = bar_Hit->get_z(0);
	  mG4EvtTree.x1[nbarhits] = bar_Hit->get_x(1);
	  mG4EvtTree.y1[nbarhits] = bar_Hit->get_y(1);
	  mG4EvtTree.z1[nbarhits] = bar_Hit->get_z(1);

	  nbarhits++;
	}
    }
  return 0;
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

      /*mG4EvtTree.x0[nhits] = dirc_hit->get_x(0);
      mG4EvtTree.y0[nhits] = dirc_hit->get_y(0);
      mG4EvtTree.z0[nhits] = dirc_hit->get_z(0);
      mG4EvtTree.x1[nhits] = dirc_hit->get_x(1);
      mG4EvtTree.y1[nhits] = dirc_hit->get_y(1);
      mG4EvtTree.z1[nhits] = dirc_hit->get_z(1);*/
      mG4EvtTree.edep[nhits] = dirc_hit->get_edep();

      //mG4EvtTree.track_angle_at_bar[nhits] = dirc_hit->GetAngleTrack();
      //mG4EvtTree.track_phi[nhits] = dirc_hit->GetMomentum().Phi();
      //mG4EvtTree.bar_hit_time[nhits] = dirc_hit->GetBarHitTime();
      
      mG4EvtTree.mcp_id[nhits] = dirc_hit->GetMcpId();
      mG4EvtTree.pixel_id[nhits] = dirc_hit->GetPixelId();
      mG4EvtTree.lead_time[nhits] = dirc_hit->GetLeadTime();
      mG4EvtTree.wavelength[nhits] = dirc_hit->GetTotTime();
      //mG4EvtTree.hit_pathId[nhits] = dirc_hit->GetPathInPrizm();
      //mG4EvtTree.nrefl[nhits] = dirc_hit->GetNreflectionsInPrizm();

      for (int i = 0; i < 3; i++)
      {
        //mG4EvtTree.track_mom[nhits][i] = dirc_hit->GetMomentum()(i);
        //mG4EvtTree.track_pos[nhits][i] = dirc_hit->GetPosition()(i);
        mG4EvtTree.hit_globalPos[nhits][i] = dirc_hit->GetGlobalPos()(i);
	mG4EvtTree.hit_localPos[nhits][i] = dirc_hit->GetLocalPos()(i);
	mG4EvtTree.hit_digiPos[nhits][i] = dirc_hit->GetDigiPos()(i);
	mG4EvtTree.hit_mom[nhits][i] = dirc_hit->GetMomentum()(i);
	mG4EvtTree.hit_pos[nhits][i] = dirc_hit->GetPosition()(i);	    
	  
      }

      /*int id = 300*dirc_hit->GetMcpId() + dirc_hit->GetPixelId();
      ((PrtLutNode*)(fLut->At(id)))->
	AddEntry(id, dir_vec, dirc_hit->GetPathInPrizm(),
		 dirc_hit->GetNreflectionsInPrizm(),
		 dirc_hit->GetLeadTime(), dirc_hit->GetGlobalPos(), dirc_hit->GetDigiPos());
      */
      nhits++;
    }

  }

  return 0;
}
