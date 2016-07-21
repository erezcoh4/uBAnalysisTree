#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"

int ana2csv_data()
{
  Int_t   trackid; 
  Int_t   flip;
  Int_t   primary;     
  Int_t   nhits;
  Float_t length;
  Float_t startx;
  Float_t starty;
  Float_t startz;
  Float_t endx;
  Float_t endy;
  Float_t endz;
  Float_t ttheta;
  Float_t tphi;
  Float_t distlenratio;
  Float_t startdqdx;
  Float_t enddqdx;
  Float_t dqdxdiff;
  Float_t dqdxratio;
  Float_t totaldqdx;
  Float_t averagedqdx;
  Float_t cosmicscore;
  Float_t coscontscore;
  Float_t pidpida;
  Float_t pidchi;
  Float_t cftime;
  Float_t cftimewidth;
  Float_t cfzcenter;
  Float_t cfzwidth;
  Float_t cfycenter;
  Float_t cfywidth;
  Float_t cftotalpe; 
  Float_t cfdistance;

  Short_t startwire;
  Short_t endwire;
  Float_t starttime;
  Float_t endtime;

  Float_t gstartx;
  Float_t gstarty;
  Float_t gstartz;
  Float_t gendx;
  Float_t gendy;
  Float_t gendz;
  Float_t gtheta;
  Float_t gphi;
  Float_t gstartdqdx;
  Float_t genddqdx;

  static Int_t   run;
  static Int_t   subrun;
  static Int_t   event;
  static Int_t   ntracks_pandoraNu;
  static Int_t   trkId_pandoraNu[50];
  static Int_t   trkg4id_pandoraNu[50];
  static Float_t trklen_pandoraNu[50];
  static Float_t trkstartx_pandoraNu[50];
  static Float_t trkstarty_pandoraNu[50];
  static Float_t trkstartz_pandoraNu[50];
  static Float_t trkendx_pandoraNu[50];
  static Float_t trkendy_pandoraNu[50];
  static Float_t trkendz_pandoraNu[50];
  static Float_t trktheta_pandoraNu[50];
  static Float_t trkphi_pandoraNu[50];
  static Short_t ntrkhits_pandoraNu[50][3];
  static Float_t trkdqdx_pandoraNu[50][3][5000];
  static Float_t trkresrg_pandoraNu[50][3][5000];
  static Float_t trkxyz_pandoraNu[50][3][5000][3];
  static Short_t trkncosmictags_tagger_pandoraNu[50]; 
  static Float_t trkcosmicscore_tagger_pandoraNu[50][10];
  static Short_t trkcosmictype_tagger_pandoraNu[50][10];
  static Short_t trkncosmictags_containmenttagger_pandoraNu[50];
  static Float_t trkcosmicscore_containmenttagger_pandoraNu[50][10];
  static Short_t trkcosmictype_containmenttagger_pandoraNu[50][10];
  static Float_t trkpidchi_pandoraNu[50][3];
  static Float_t trkpidpida_pandoraNu[50][3];
  static Short_t trkpidbestplane_pandoraNu[50];

  static Int_t   no_hits;
  static Short_t hit_plane[5000];
  static Short_t hit_wire[5000];
  static Float_t hit_peakT[5000];
  static Short_t hit_trkid[5000];
  static Short_t hit_trkKey[5000];

  static Int_t   no_flashes;
  static Float_t flash_time[5000];
  static Float_t flash_timewidth[5000];
  static Float_t flash_pe[5000];
  static Float_t flash_ycenter[5000];
  static Float_t flash_ywidth[5000];
  static Float_t flash_zcenter[5000];
  static Float_t flash_zwidth[5000];

  // open csv
  std::ofstream ofile("features_ana.csv");
  ofile << "run,subrun,event,trackid,flip,nhits,length," 
        << "startx,starty,startz,endx,endy,endz," 
        << "theta,phi,distlenratio,startdqdx,enddqdx," 
        << "dqdxdiff,dqdxratio,totaldqdx,averagedqdx," 
        << "cosmicscore,coscontscore,pidpida,pidchi,"  
        << "cftime,cftimewidth,cfzcenter,cfzwidth,cfycenter,cfywidth,cftotalpe,cfdistance,"
        << "startwire,endwire,starttick,endtick"
        << std::endl;


  // open file
  std::ifstream infile("filesana_BNBEXT.list");

  std::string filename;
  while(infile >> filename)
  {
    //make filename into char array
    char* a = new char[filename.size() + 1];
    a[filename.size()] = 0;
    memcpy(a,filename.c_str(),filename.size());
    std::cout << a << std::endl;

    TFile* f = new TFile(a,"read");
    // open tree
    TTree* tree = (TTree*)f->Get("analysistree/anatree");

    // Get Branches
    TBranch* b_run                                         = tree->GetBranch("run");
    TBranch* b_subrun                                      = tree->GetBranch("subrun");
    TBranch* b_event                                       = tree->GetBranch("event");
    TBranch* b_ntracks_pandoraNu                           = tree->GetBranch("ntracks_pandoraNu");
    TBranch* b_trkId_pandoraNu                             = tree->GetBranch("trkId_pandoraNu");
    TBranch* b_trkg4id_pandoraNu                           = tree->GetBranch("trkg4id_pandoraNu");
    TBranch* b_trklen_pandoraNu                            = tree->GetBranch("trklen_pandoraNu");
    TBranch* b_trkstartx_pandoraNu                         = tree->GetBranch("trkstartx_pandoraNu");
    TBranch* b_trkstarty_pandoraNu                         = tree->GetBranch("trkstarty_pandoraNu");
    TBranch* b_trkstartz_pandoraNu                         = tree->GetBranch("trkstartz_pandoraNu");
    TBranch* b_trkendx_pandoraNu                           = tree->GetBranch("trkendx_pandoraNu");
    TBranch* b_trkendy_pandoraNu                           = tree->GetBranch("trkendy_pandoraNu");
    TBranch* b_trkendz_pandoraNu                           = tree->GetBranch("trkendz_pandoraNu");
    TBranch* b_trktheta_pandoraNu                          = tree->GetBranch("trktheta_pandoraNu");
    TBranch* b_trkphi_pandoraNu                            = tree->GetBranch("trkphi_pandoraNu");
    TBranch* b_ntrkhits_pandoraNu                          = tree->GetBranch("ntrkhits_pandoraNu");
    TBranch* b_trkdqdx_pandoraNu                           = tree->GetBranch("trkdqdx_pandoraNu");
    TBranch* b_trkresrg_pandoraNu                          = tree->GetBranch("trkresrg_pandoraNu");
    TBranch* b_trkxyz_pandoraNu                            = tree->GetBranch("trkxyz_pandoraNu");
    TBranch* b_trkncosmictags_tagger_pandoraNu             = tree->GetBranch("trkncosmictags_tagger_pandoraNu");
    TBranch* b_trkcosmicscore_tagger_pandoraNu             = tree->GetBranch("trkcosmicscore_tagger_pandoraNu");
    TBranch* b_trkcosmictype_tagger_pandoraNu              = tree->GetBranch("trkcosmictype_tagger_pandoraNu");
    TBranch* b_trkncosmicscore_containmenttagger_pandoraNu = tree->GetBranch("trkcosmicscore_containmenttagger_pandoraNu");
    TBranch* b_trkcosmicscore_containmenttagger_pandoraNu  = tree->GetBranch("trkcosmicscore_containmenttagger_pandoraNu");
    TBranch* b_trkcosmictype_containmenttagger_pandoraNu   = tree->GetBranch("trkcosmictype_containmenttagger_pandoraNu");
    TBranch* b_trkpidchi_pandoraNu                         = tree->GetBranch("trkpidchi_pandoraNu");
    TBranch* b_trkpidpida_pandoraNu                        = tree->GetBranch("trkpidpida_pandoraNu");
    TBranch* b_trkpidbestplane_pandoraNu                   = tree->GetBranch("trkpidbestplane_pandoraNu");

    TBranch* b_no_hits                                     = tree->GetBranch("no_hits");
    TBranch* b_hit_plane                                   = tree->GetBranch("hit_plane");
    TBranch* b_hit_wire                                    = tree->GetBranch("hit_wire");
    TBranch* b_hit_peakT                                   = tree->GetBranch("hit_peakT");
    TBranch* b_hit_trkid                                   = tree->GetBranch("hit_trkid");
    TBranch* b_hit_trkKey                                  = tree->GetBranch("hit_trkKey");

    TBranch* b_no_flashes                                  = tree->GetBranch("no_flashes");
    TBranch* b_flash_time                                  = tree->GetBranch("flash_time");
    TBranch* b_flash_timewidth                             = tree->GetBranch("flash_timewidth");
    TBranch* b_flash_pe                                    = tree->GetBranch("flash_pe");
    TBranch* b_flash_ycenter                               = tree->GetBranch("flash_ycenter");
    TBranch* b_flash_ywidth                                = tree->GetBranch("flash_ywidth");
    TBranch* b_flash_zcenter                               = tree->GetBranch("flash_zcenter");
    TBranch* b_flash_zwidth                                = tree->GetBranch("flash_zwidth");
  
    // Set addresses
    b_run->SetAddress(&run);
    b_subrun->SetAddress(&subrun);
    b_event->SetAddress(&event);
    b_ntracks_pandoraNu->SetAddress(&ntracks_pandoraNu);
    b_trkId_pandoraNu->SetAddress(&trkId_pandoraNu);
    b_trkg4id_pandoraNu->SetAddress(&trkg4id_pandoraNu);
    b_trklen_pandoraNu->SetAddress(&trklen_pandoraNu);
    b_trkstartx_pandoraNu->SetAddress(&trkstartx_pandoraNu);
    b_trkstarty_pandoraNu->SetAddress(&trkstarty_pandoraNu);
    b_trkstartz_pandoraNu->SetAddress(&trkstartz_pandoraNu);
    b_trkendx_pandoraNu->SetAddress(&trkendx_pandoraNu);
    b_trkendy_pandoraNu->SetAddress(&trkendy_pandoraNu);
    b_trkendz_pandoraNu->SetAddress(&trkendz_pandoraNu);
    b_trktheta_pandoraNu->SetAddress(&trktheta_pandoraNu);
    b_trkphi_pandoraNu->SetAddress(&trkphi_pandoraNu);
    b_ntrkhits_pandoraNu->SetAddress(&ntrkhits_pandoraNu);
    b_trkdqdx_pandoraNu->SetAddress(&trkdqdx_pandoraNu);
    b_trkresrg_pandoraNu->SetAddress(&trkresrg_pandoraNu);
    b_trkxyz_pandoraNu->SetAddress(&trkxyz_pandoraNu);
    b_trkncosmictags_tagger_pandoraNu->SetAddress(&trkncosmictags_tagger_pandoraNu); 
    b_trkcosmicscore_tagger_pandoraNu->SetAddress(&trkcosmicscore_tagger_pandoraNu);
    b_trkcosmictype_tagger_pandoraNu->SetAddress(&trkcosmictype_tagger_pandoraNu);
    b_trkncosmicscore_containmenttagger_pandoraNu->SetAddress(&trkcosmicscore_containmenttagger_pandoraNu);
    b_trkcosmicscore_containmenttagger_pandoraNu->SetAddress(&trkcosmicscore_containmenttagger_pandoraNu);
    b_trkcosmictype_containmenttagger_pandoraNu->SetAddress(&trkcosmictype_containmenttagger_pandoraNu);
    b_trkpidchi_pandoraNu->SetAddress(&trkpidchi_pandoraNu);
    b_trkpidpida_pandoraNu->SetAddress(&trkpidpida_pandoraNu);
    b_trkpidbestplane_pandoraNu->SetAddress(&trkpidbestplane_pandoraNu);

    b_no_hits->SetAddress(&no_hits);
    b_hit_plane->SetAddress(&hit_plane);
    b_hit_wire->SetAddress(&hit_wire);
    b_hit_peakT->SetAddress(&hit_peakT);
    b_hit_trkid->SetAddress(&hit_trkid);
    b_hit_trkKey->SetAddress(&hit_trkKey);

    b_no_flashes->SetAddress(&no_flashes);
    b_flash_time->SetAddress(&flash_time);
    b_flash_timewidth->SetAddress(&flash_timewidth);
    b_flash_pe->SetAddress(&flash_pe);
    b_flash_ycenter->SetAddress(&flash_ycenter);
    b_flash_ywidth->SetAddress(&flash_ywidth);
    b_flash_zcenter->SetAddress(&flash_zcenter);
    b_flash_zwidth->SetAddress(&flash_zwidth);
    
    Int_t nentries  = (Int_t)tree->GetEntries();
    std::cout << "Looping over " << nentries << " events." << std::endl;
    for(Int_t i = 0; i<nentries;i++)
    {
      b_run->GetEntry(i);
      b_subrun->GetEntry(i);
      b_event->GetEntry(i);
      b_ntracks_pandoraNu->GetEntry(i);
      if(ntracks_pandoraNu == 0) continue;
    
      b_trkId_pandoraNu->GetEntry(i);  
      b_trkg4id_pandoraNu->GetEntry(i);
      b_trklen_pandoraNu->GetEntry(i);
      b_trkstartx_pandoraNu->GetEntry(i);
      b_trkstarty_pandoraNu->GetEntry(i);
      b_trkstartz_pandoraNu->GetEntry(i);
      b_trkendx_pandoraNu->GetEntry(i);
      b_trkendy_pandoraNu->GetEntry(i);
      b_trkendz_pandoraNu->GetEntry(i);
      b_trktheta_pandoraNu->GetEntry(i);
      b_trkphi_pandoraNu->GetEntry(i);
      b_ntrkhits_pandoraNu->GetEntry(i);
      b_trkdqdx_pandoraNu->GetEntry(i);
      b_trkresrg_pandoraNu->GetEntry(i);
      b_trkxyz_pandoraNu->GetEntry(i);
      b_trkncosmictags_tagger_pandoraNu->GetEntry(i); 
      b_trkcosmicscore_tagger_pandoraNu->GetEntry(i);
      b_trkcosmictype_tagger_pandoraNu->GetEntry(i);
      b_trkncosmicscore_containmenttagger_pandoraNu->GetEntry(i);
      b_trkcosmicscore_containmenttagger_pandoraNu->GetEntry(i);
      b_trkcosmictype_containmenttagger_pandoraNu->GetEntry(i);
      b_trkpidchi_pandoraNu->GetEntry(i);
      b_trkpidpida_pandoraNu->GetEntry(i);
      b_trkpidbestplane_pandoraNu->GetEntry(i);

      b_no_hits->GetEntry(i);
      b_hit_plane->GetEntry(i);
      b_hit_wire->GetEntry(i);
      b_hit_peakT->GetEntry(i);
      b_hit_trkid->GetEntry(i);
      b_hit_trkKey->GetEntry(i);
  
      b_no_flashes->GetEntry(i);
      b_flash_time->GetEntry(i);
      b_flash_timewidth->GetEntry(i);
      b_flash_pe->GetEntry(i);
      b_flash_ycenter->GetEntry(i);
      b_flash_ywidth->GetEntry(i);
      b_flash_zcenter->GetEntry(i);
      b_flash_zwidth->GetEntry(i);

      // Get list of in-time flashes for event
      std::vector<int> goodflashidx;
      if(no_flashes > 0)
      {
        for(Int_t iflash=0; iflash < no_flashes; iflash++)
        {
          if(flash_time[iflash] > 0 && flash_time[iflash] < 10 && flash_pe[iflash] > 6.5)
          {
            goodflashidx.push_back(iflash);
          }
        }
      }

      // loop over all reconstructed tracks
      for(Int_t j=0; j < ntracks_pandoraNu; j++)
      {
        // check if contained
        if(trkstartx_pandoraNu[j] < 3.45    | trkstartx_pandoraNu[j] > 255.8 ) continue;
        if(trkendx_pandoraNu[j]   < 3.45    | trkendx_pandoraNu[j]   > 255.8 ) continue;
        if(trkstarty_pandoraNu[j] < -110.53 | trkstarty_pandoraNu[j] > 112.47) continue;
        if(trkendy_pandoraNu[j]   < -110.53 | trkendy_pandoraNu[j]   > 112.47) continue;
        if(trkstartz_pandoraNu[j] < 5.1     | trkstartz_pandoraNu[j] > 1031.9) continue;
        if(trkendz_pandoraNu[j]   < 5.1     | trkendz_pandoraNu[j]   > 1031.9) continue;

        trackid    = trkg4id_pandoraNu[j];
        length     = trklen_pandoraNu[j];
        startx     = trkstartx_pandoraNu[j];
        starty     = trkstarty_pandoraNu[j];
        startz     = trkstartz_pandoraNu[j];
        endx       = trkendx_pandoraNu[j];
        endy       = trkendy_pandoraNu[j];
        endz       = trkendz_pandoraNu[j];
        ttheta     = trktheta_pandoraNu[j];
        tphi       = trkphi_pandoraNu[j];
       

        // get flash info
        // compare reconstructed track to list of flashes in beam to find closest
        float tzcenter = (startz + endz)/2.;
        if(goodflashidx.size() > 0)
        {
          int   minzi    = goodflashidx.at(0);
          float minzdiff = TMath::Abs(flash_zcenter[minzi] - tzcenter);
          for(size_t k=0; k < goodflashidx.size(); k++)
          {
            int   fidx     = goodflashidx.at(k);
            float fzcenter = flash_zcenter[fidx];
            if(TMath::Abs(fzcenter - tzcenter) < minzdiff)
            {
              minzi    = fidx;
              minzdiff = TMath::Abs(fzcenter = tzcenter);
            }
          }
          cftime      = flash_time[minzi];
          cftimewidth = flash_timewidth[minzi];
          cfzcenter   = flash_zcenter[minzi];
          cfzwidth    = flash_zwidth[minzi];
          cfycenter   = flash_ycenter[minzi];
          cfywidth    = flash_ywidth[minzi];
          cftotalpe   = flash_pe[minzi];
          cfdistance  = tzcenter - cfzcenter;
        }
        else
        {
          cftime      = -9999; 
          cftimewidth = -9999; 
          cfzcenter   = -9999; 
          cfzwidth    = -9999; 
          cfycenter   = -9999; 
          cfywidth    = -9999; 
          cftotalpe   = -9999; 
          cfdistance  = -9999; 
        }

        // get cosmic scores
        cosmicscore  = trkcosmicscore_tagger_pandoraNu[j][0];
        coscontscore = trkcosmicscore_containmenttagger_pandoraNu[j][0];

        // get pid info
        pidpida = trkpidpida_pandoraNu[j][trkpidbestplane_pandoraNu[j]];
        pidchi  = trkpidchi_pandoraNu[j][trkpidbestplane_pandoraNu[j]];
        
        // get dqdx info
        // loop over range from end of track to find start and end
        int   rmin[3];
        int   rmax[3];
        float totaldqdx = 0.;
        float startdqdx = 0.;
        float enddqdx   = 0.;
        int   nhits     = 0;
        for(Int_t fr=0; fr<3;fr++)
        {
          if(ntrkhits_pandoraNu[j][fr] >= 0)
          {
            rmin[fr]   = trkresrg_pandoraNu[j][fr][0];
            rmax[fr]   = trkresrg_pandoraNu[j][fr][0];
            totaldqdx += trkdqdx_pandoraNu[j][fr][0];
            int minidx = 0;
            int maxidx = 0;
            for(Int_t ridx=0; ridx < ntrkhits_pandoraNu[j][fr]; ridx++)
            {
              if(trkresrg_pandoraNu[j][fr][ridx] < rmin[fr] && trkdqdx_pandoraNu[j][fr][ridx] != 0)
              {
                rmin[fr] = trkresrg_pandoraNu[j][fr][ridx];
                minidx   = ridx;
              }
              if(trkresrg_pandoraNu[j][fr][ridx] > rmax[fr])
              {
                rmax[fr] = trkresrg_pandoraNu[j][fr][ridx];
                maxidx   = ridx;
              }
              totaldqdx += trkdqdx_pandoraNu[j][fr][ridx];
            }
            if(rmax[fr] != rmin[fr])
            {
              nhits       += ntrkhits_pandoraNu[j][fr];
            }
            if(maxidx >= 3)
            {
              startdqdx   += (trkdqdx_pandoraNu[j][fr][maxidx] + trkdqdx_pandoraNu[j][fr][maxidx-1]
                              + trkdqdx_pandoraNu[j][fr][maxidx-2]);
              enddqdx     += (trkdqdx_pandoraNu[j][fr][minidx] + trkdqdx_pandoraNu[j][fr][minidx+1]
                              + trkdqdx_pandoraNu[j][fr][minidx+2]);
            }
            else
            {
              startdqdx   += trkdqdx_pandoraNu[j][fr][maxidx];
              enddqdx     += trkdqdx_pandoraNu[j][fr][minidx];
            }
          }
        }
        // do our own pida
        //
        
        if(nhits > 0) averagedqdx = totaldqdx/nhits;
        else averagedqdx = 0;

        // fix start and end points by charge
        // If start charge is greater, flip track
        if(startdqdx > enddqdx)
        {
          gstartx    = endx;
          gstarty    = endy;
          gstartz    = endz;
          gendx      = startx;
          gendy      = starty;
          gendz      = startz;
          gtheta     = 3.1416 - ttheta;
          gphi       = (tphi > 0 ? tphi-3.1416 : tphi+3.1416);
          gstartdqdx = enddqdx;
          genddqdx   = startdqdx;
          flip       = 1;
        }
        else
        {
          gstartx    = startx;
          gstarty    = starty;
          gstartz    = startz;
          gendx      = endx;
          gendy      = endy;
          gendz      = endz;
          gtheta     = ttheta;
          gphi       = tphi;
          gstartdqdx = startdqdx;
          genddqdx   = enddqdx;
          flip       = 0;
        }

        // get wire and time tick estimates
        // to help find tracks in event displays
        startwire    = gstartz*3455./1036.;
        endwire      = gendz*3455./1036.;
        starttime    = (gstartx + 45.)*6400./357.;
        endtime      = (gendx + 45.)*6400./357.;

        // make some caloritric variables
        dqdxdiff     = genddqdx - gstartdqdx;
        dqdxratio    = genddqdx/gstartdqdx;

        // get straightness (distance from start to end of track over reco'd length)
        float dist   = TMath::Sqrt(TMath::Power(gendx-gstartx,2) + TMath::Power(gendy-gstarty,2) + TMath::Power(gendz-gstartz,2));
        distlenratio = dist/length;

        // dump all variable to a line of csv
        ofile << run << "," << subrun << "," << event << "," << trackid << "," << flip << "," << nhits << "," << length << "," 
              << gstartx << "," << gstarty << "," << gstartz << "," << gendx << "," << gendy << "," << gendz << "," 
              << gtheta << "," << gphi << "," << distlenratio << "," << gstartdqdx << "," << genddqdx << "," 
              << dqdxdiff << "," << dqdxratio << "," << totaldqdx << "," << averagedqdx << "," 
              << cosmicscore << "," << coscontscore << "," << pidpida << "," << pidchi << ","  
              << cftime << "," << cftimewidth << "," << cfzcenter << "," << cfzwidth << "," << cfycenter << "," << cfywidth << "," << cftotalpe << "," << cfdistance << ","
              << startwire << "," << endwire << "," << starttime << "," << endtime
              << std::endl;
      }
    }
    f->Close();
  }
  ofile.close();
  return 0;
}
