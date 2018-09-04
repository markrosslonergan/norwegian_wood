void split_file(){

TFile *f = new TFile("DUNE_ntuple.root","read");

f->GetListOfKeys()->Print(); 

for(const auto&& obj: *(f->GetListOfKeys())){
   TTree *t = (TTree*)f->Get(obj->GetName());
   
   TFile* fout = new TFile(Form("%s%s%s","split/",obj->GetName(),".root"),"recreate");
   TTree *tnew = t->CloneTree();
   fout->cd();


   fout->cd();
   tnew->Write();
   fout->Close();

   f->cd();
}


f->Close();
}
