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

int ana2csv_example()
{
  // In analysistree, ntracks is just an int, but all of the variables about the tracks are arrays
  // You have to make these arrays in the largest size that they could be or root will error
  // I'm assuming here that there won't be more than 100 reconstructed pandoraNu tracks per event
  static Int_t   ntracks_pandoraNu;
  static Float_t trklen_pandoraNu[100];
  static Float_t trkstarty_pandoraNu[100];
  static Float_t trkendy_pandoraNu[100];

  // open csv
  /*
  This part creates opens a csv file and writes the headers to the first line
  */
  std::ofstream ofile("ylength.csv");
  ofile << "ystart,yend,ydiff,length" 
        << std::endl;

  // open file
  /*
  These next few lines reads files from a list.
  If you just want to use one files, you would cut out everything between here and "TFile* f = ..."
  */
  std::ifstream infile("filesana.list");

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
    /*
    This part creates branch objects
    */
    TBranch* b_ntracks_pandoraNu                           = tree->GetBranch("ntracks_pandoraNu");
    TBranch* b_trklen_pandoraNu                            = tree->GetBranch("trklen_pandoraNu");
    TBranch* b_trkstarty_pandoraNu                         = tree->GetBranch("trkstarty_pandoraNu");
    TBranch* b_trkendy_pandoraNu                           = tree->GetBranch("trkendy_pandoraNu");

    // Set addresses
    /*
    This connects the branches to the variables you use to access them
    */
    b_ntracks_pandoraNu->SetAddress(&ntracks_pandoraNu);
    b_trklen_pandoraNu->SetAddress(&trklen_pandoraNu);
    b_trkstarty_pandoraNu->SetAddress(&trkstarty_pandoraNu);
    b_trkendy_pandoraNu->SetAddress(&trkendy_pandoraNu);
    
    Int_t nentries  = (Int_t)tree->GetEntries();
    std::cout << "Looping over " << nentries << " events." << std::endl;
    for(Int_t i = 0; i<nentries;i++)
    {
      /*
      Now we're looping over each event.
      This part gets a new entry for each variable.
      Now we can actually do stuff with the anatree variables.
      */
      b_ntracks_pandoraNu->GetEntry(i);
      b_trklen_pandoraNu->GetEntry(i);
      b_trkstarty_pandoraNu->GetEntry(i);
      b_trkendy_pandoraNu->GetEntry(i);

      /*
      Each event has multiple reconstructed tracks, so here we will loop over each track.
      */ 
      for(Int_t j = 0; j<ntracks_pandoraNu; j++)
      {
        // The length in each dimension isn't stored in analysistree for the track objects
        // so we can just take the difference in y start and y end
        // I'm also outputting the 3D length 
        float trkydiff = TMath::Abs(trkstarty_pandoraNu[j] - trkendy_pandoraNu[j]);

        // This writes each the variables that you want to save to a line in the csv file
        // You have to be careful to keep the same order as the headers that we wrote above
        ofile << trkstarty_pandoraNu[j] << "," << trkendy_pandoraNu[j] << "," 
              << trkydiff << "," << trklen_pandoraNu[j]
              << std::endl;
      }
    }
    f->Close();
  }
  ofile.close();
  return 0;
}
