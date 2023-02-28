#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include <iostream>
#include <fstream>


void GetEnergy(Double_t M,Double_t IZ,Double_t BRO,Double_t &E);


std::vector<AtHitCluster> fHitClusterArray; //< Clusterized hits container

std::vector<AtHitCluster> *GetHitClusterArray() { 
 

         return &fHitClusterArray;


 }



void Bepp_ana(Int_t nEvents = 1000)
{


   TH2F *angle_vs_energy = new TH2F("angle_vs_energy", "angle_vs_energy", 720, 0, 179, 1000, 0, 100.0);
   TH2F *Z_vs_Y = new TH2F("Z_vs_Y", "energy_vs_clusterangle", 720, 0, -3, 1000, 0, 2.0);
   TH2F *energy_vs_Zorb = new TH2F("energy_vs_Zorb", "energy_vs_Zorb", 720, 0, -5, 1000, 0, 100.0);
   TCanvas *c1 = new TCanvas();
   c1->Divide(2, 2);
   c1->Draw();



   FairRunAna *run = new FairRunAna();
   std::vector<TString> files = {"data/output_digi"};

   for (auto iFile = 0; iFile < files.size(); ++iFile) {

       TString mcFileNameHead = files[iFile];
       TString mcFileNameTail = ".root";
       TString mcFileName = mcFileNameHead + mcFileNameTail;
  // std:
    //  cout << " Analysis of simulation file  " << mcFileName << endl;

       TFile *file = new TFile(mcFileName.Data(), "READ");
       TTree *tree = (TTree *)file->Get("cbmsim");

       tree = (TTree *)file->Get("cbmsim");
       // TBranch *branch = tree->GetBranch("AtTpcPoint");
       TTreeReader Reader1("cbmsim", file);
       TTreeReaderValue<TClonesArray> eventArray(Reader1, "AtPatternEvent");
     // TTreeReaderValue<TClonesArray> SimPointArray(Reader1, "AtMCPoint");

       Double_t deltaPhi = 0.0;
       Double_t POCAOrbZ = 1E6;
       Double_t firstOrbZ = 0.0;
       Double_t phiOrbZ = 0.0;
       Double_t lengthOrbZ = 0.0;
       Double_t length = 0;
       Double_t phiClus =0;
       Double_t zIniCal = 0;


       ROOT::Math::XYZPoint iniPos;
       ROOT::Math::XYZPoint secPos;

       for (Int_t i = 0; i < nEvents; i++) {

           //std::cout << " Event Number : " << i << "\n";

           Reader1.Next();

           AtPatternEvent *patternEvent = (AtPatternEvent *)eventArray->At(0);

           if (patternEvent) {
              std::vector<AtTrack> &patternTrackCand = patternEvent->GetTrackCand();
              //std::cout << " Number of pattern tracks " << patternTrackCand.size() << "\n";
              for (auto track : patternTrackCand) {

	         // std::cout << " === Track " << track.GetTrackID() << " with "
                        // << track.GetHitClusterArray()->size() << " clusters "  << "\n";


		  if ( track.GetHitClusterArray()->size() < 5) {
                     //std::cout << " Track is noise or has less than 5 clusters! "  << "\n";
                     continue;
                  }
		  auto hitClusterArray = track.GetHitClusterArray();
                  AtHitCluster iniCluster;
                  AtHitCluster SecCluster;

                  Double_t theta = track.GetGeoTheta();
                  Double_t rad   = track.GetGeoRadius();
                  //Double_t phi = track.SetGeoPhi(-phiClus);
                  //std::vector<AtHitCluster> *hitClusterArray = track.GetHitClusterArray();
                  Double_t B_f = 3.0;
                  double bro = B_f * rad / TMath::Sin(theta) / 1000.0;
                  double ener = 0;
                  Double_t Am = 1.0;

                 // GetEnergy(Am, 1.0, bro, ener);

                  //angle_vs_energy->Fill(180.0 - theta * TMath::RadToDeg(), ener * Am);

                  if (i == 17 || i == 141 || i== 263 || i==347 || i== 395 || i== 807){
		     std::cout<< "Processing event " << i  << "with " << track.GetHitClusterArray()->size() << " clusters" << endl;

                     for (auto iclus = 1; iclus < hitClusterArray->size(); ++iclus) {
                         SecCluster = hitClusterArray->at(iclus-1);
                         iniCluster = hitClusterArray->at(iclus);

		         auto cluster_inipos = iniCluster.GetPosition();
                         auto cluster_secpos = SecCluster.GetPosition();

                         auto dir = cluster_secpos - cluster_inipos;
                         length += dir.mag2();

                         Double_t phiInc = TMath::ATan2(cluster_secpos.Y()-cluster_inipos.Y(),cluster_inipos.X()-cluster_secpos.X());//calculate angle between the two clusters.
 
                         Double_t distance = TMath::Sqrt(cluster_secpos.X() * cluster_secpos.X() + cluster_secpos.Y() * cluster_secpos.Y());
                         phiInc = (phiInc > 0) ? phiInc : 2.0 * TMath::Pi() + phiInc;
                         deltaPhi += phiInc;

                      //std::cout<<phiInc*TMath::RadToDeg()<<"    "<<deltaPhi*TMath::RadToDeg()<<"\n";
                      //std::cout<<" Distance to z axis : "<<distance << " " << POCAOrbZ << endl;;
                      //std::cout<<" Length : "<<length<<"\n";
                      //std::cout << cluster_secpos.X()<< " "<<" Z projection : "<<cluster_secpos.Z()<<"\n";
                      // std::cout<<" momLoss "<<momLoss<<"\n";
                        // energy_vs_Zorb->Fill(distance,ener*Am);

                      // NB: Find poca after at least half turn.
	                 if (distance < POCAOrbZ && deltaPhi > TMath::Pi()){
			    POCAOrbZ = distance;
		            firstOrbZ=cluster_secpos.Z();
			    phiOrbZ = deltaPhi;
                            lengthOrbZ = length;

                         }
 
		  }
                  double bro = B_f * rad / TMath::Sin(theta) / 1000.0;
                  double ener = 0;
                  Double_t Am = 1.0;
                  Double_t Velocity = 1.0;
                  Double_t T_cycle = 21.9;   //in nano seconds. 
                  GetEnergy(Am, 1.0, bro,ener);
                  std::cout<< firstOrbZ<<endl;
                  //Double_t Energy_lab = ener*Am - (1/2)*Am*pow(Velocity,2) + ((Am *Velocity)/T_cycle)*firstOrbZ;
                  energy_vs_Zorb->Fill(firstOrbZ,ener*Am);


                  }

	      }
           }
       }
   }


   Double_t *ThetaCMS = new Double_t[20000];
   Double_t *ThetaLabRec = new Double_t[20000];
   Double_t *EnerLabRec = new Double_t[20000];
   Double_t *ThetaLabSca = new Double_t[20000];
   Double_t *EnerLabSca = new Double_t[20000];
   Double_t *MomLabRec = new Double_t[20000];

   
/*   TString fileKine = "O16_aa_el_kine.txt";
   std::ifstream *kineStr = new std::ifstream(fileKine.Data());
   Int_t numKin = 0;

   if (!kineStr->fail()) {
      while (!kineStr->eof()) {
         *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >> ThetaLabSca[numKin] >>
            EnerLabSca[numKin];
         // numKin++;

         // MomLabRec[numKin] =( pow(EnerLabRec[numKin] + M_Ener,2) - TMath::Power(M_Ener, 2))/1000.0;
         // std::cout<<" Momentum : " <<MomLabRec[numKin]<<"\n";
         // Double_t E = TMath::Sqrt(TMath::Power(p, 2) + TMath::Power(M_Ener, 2)) - M_Ener;
         numKin++;
      }
   } else if (kineStr->fail())
      std::cout << " Warning : No Kinematics file found for this reaction!" << std::endl;

    TGraph *Kine_AngRec_EnerRec = new TGraph(numKin, ThetaLabRec, EnerLabRec);
*/   
   //TCanvas *c1 = new TCanvas();
   c1->cd(1);
   angle_vs_energy->Draw();
   //c1->cd(2);
   //energy_vs_clusterangle->Draw();
   //energy_vs_clusterangle->SetMarkerStyle(20);
   c1->cd(2);
   energy_vs_Zorb->Draw();
   energy_vs_Zorb->SetMarkerStyle(20);

}


void  GetEnergy(Double_t M,Double_t IZ,Double_t BRO,Double_t &E){

  //Energy per nucleon

 Double_t Velocity = 1.0;
 Double_t T_cycle = 21.9; 

  Float_t  AM=931.5;
  Float_t X=BRO/0.1439*IZ/M;
  X=pow(X,2);
  X=2.*AM*X;
  X=X+pow(AM,2);
  E=TMath::Sqrt(X)-AM;

 // Double_t Energy_lab = E -(1/2)*Am*pow(Velocity,2) +((Am *Velocity)/T_cycle)*firstOrbZ;


  //std::cout<<E<< endl;
  }


