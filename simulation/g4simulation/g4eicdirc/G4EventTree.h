#ifndef G4EVENTTREE_H
#define G4EVENTTREE_H

#define MAXHIT 500000

#include <Rtypes.h>

typedef struct
{
  // Event Level
  int nhits;
  
  // Hit level
  int detid[MAXHIT];
  int hitid[MAXHIT];
  int trkid[MAXHIT];
  Double_t x0[MAXHIT];
  Double_t y0[MAXHIT];
  Double_t z0[MAXHIT];
  Double_t x1[MAXHIT];
  Double_t y1[MAXHIT];
  Double_t z1[MAXHIT];
  Double_t edep[MAXHIT];
  Int_t bar_id[MAXHIT];
  Int_t mcp_id[MAXHIT];
  Int_t pixel_id[MAXHIT];
  Double_t lead_time[MAXHIT];
  Double_t wavelength[MAXHIT];
  Double_t hit_globalPos[MAXHIT][3];
  Double_t hit_localPos[MAXHIT][3];
  Double_t hit_digiPos[MAXHIT][3];
  Double_t hit_mom[MAXHIT][3];
  Double_t hit_pos[MAXHIT][3];
  Double_t track_momentum_at_bar[MAXHIT][3];
  Long64_t hit_pathId[MAXHIT];
  Int_t nrefl[MAXHIT];
  Int_t parent_pid[MAXHIT];
  Int_t parent_track_id[MAXHIT];
  Double_t parent_momentum[MAXHIT];


} G4EventTree;

#endif
