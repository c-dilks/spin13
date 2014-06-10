void AddLumi(const char * filename="lum_perrun_FMSJP2.txt")
{
  Int_t runnum; // run number
  Double_t rstart; // seconds since 12/31/2011 for run start
  Double_t rend; // seconds since 12/31/2011 for run end
  Int_t fill; // fill no.
  Double_t lumi; // luminosity
  Double_t prescale; // trigger prescale
  Double_t live_time; // live time for this trigger from TCU
  char base_trig_name[16]; // name of trigger for lumi estimation (base trigger)
  Double_t base_trig_live_time; // live time for base trigger
  Double_t FOM1; // figure of merit as (P_B*P_B+P_Y*P_Y)/2*lum
  Double_t FOM2; // figure of merit as P_B*P_B*P_Y*P_Y*lum

  TTree * tr = new TTree();
  tr->ReadFile(filename,
   "runnum/I:rstart/D:rend/D:fill/I:lumi/D:prescale/D:live_time/D:base_trig_name/C:base_trig_live_time/D:FOM1/D:FOM2/D");

  TTree * rr = new TTree();
  rr->ReadFile("../runlist","index/I:runnum/I:fill/I:rellum/F");
  Int_t runnum_good;
  rr->SetBranchAddress("runnum",&runnum_good);

  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("rstart",&rstart);
  tr->SetBranchAddress("rend",&rend);
  tr->SetBranchAddress("fill",&fill);
  tr->SetBranchAddress("lumi",&lumi);
  tr->SetBranchAddress("prescale",&prescale);
  tr->SetBranchAddress("live_time",&live_time);
  tr->SetBranchAddress("base_trig_name",base_trig_name);
  tr->SetBranchAddress("base_trig_live_time",&base_trig_live_time);
  tr->SetBranchAddress("FOM1",&FOM1);
  tr->SetBranchAddress("FOM2",&FOM2);

  Double_t int_lumi=0;
  Bool_t add_lumi;
  /*
  for(Int_t i=0; i<rr->GetEntries(); i++)
  {
    rr->GetEntry(i);
    add_lumi=0;
    for(Int_t j=0; j<tr->GetEntries(); j++)
    {
      tr->GetEntry(j);
      if(runnum==runnum_good)
      {
        add_lumi=1;
        int_lumi+=lumi;
      };
    };
    if(!add_lumi) printf("%d in pi0 run list but not in lumi run list\n",runnum_good);
  };
  */
  for(Int_t j=0; j<tr->GetEntries(); j++)
  {
    tr->GetEntry(j);
    int_lumi+=lumi;
  };

  printf("integrated luminosity: %f\n",int_lumi);
};



