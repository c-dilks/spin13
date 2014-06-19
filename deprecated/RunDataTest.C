// tests calling methods of RunData for seg faults

void RunDataTest(Bool_t callDefault=1)
{
  gROOT->Reset();
  gSystem->Load("src/RunData13.so");
  RunData13 * r;
  if(callDefault) r = new RunData13();
  else r = new RunData13("/home/dilks/h5/root12fms/spin13");
  Int_t runnum = 14158017;
  Int_t bx=1;
  printf("runnum=%d\n",runnum);
  printf("fill=%d\n",r->GetFill(runnum));
  printf("rellum=%f\n",r->Rellum(runnum,3,"zdc"));
  printf("rellumerr=%f\n",r->RellumErr(runnum,3,"zdc"));
  printf("BluePol=%f\n",r->BluePol(runnum));
  printf("YellPol=%f\n",r->YellPol(runnum));
  printf("BlueSpin=%d\n",r->BlueSpin(runnum,bx));
  printf("YellSpin=%d\n",r->YellSpin(runnum,bx));


  
};
