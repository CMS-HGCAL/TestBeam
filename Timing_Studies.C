void Timing_Studies(){

     TCanvas *c1 = new TCanvas("c1", "c1", 800,600);   
     TFile *f1 = TFile::Open("test_DigiAndRechitPlotter_TB_4130_Ped_Profile.root");
     TFile *f2 = TFile::Open("test_DigiAndRechitPlotter_TB_4120_Ped_Profile.root");
//     TFile *f3 = TFile::Open("test_DigiAndRechitPlotter_TB_4098_Ped_Profile.root");
//     TFile *f3 = TFile::Open("test_DigiAndRechitPlotter_TB_2051_Ped_Profile.root");
     TFile *f3 = TFile::Open("test_DigiAndRechitPlotter_TB_4132_Ped_Profile.root");

     TFile *f4 = TFile::Open("test_DigiAndRechitPlotter_TB_4100_Ped_Profile.root");
     TFile *f5 = TFile::Open("test_DigiAndRechitPlotter_TB_4112_Ped_Profile.root");
     TFile *f6 = TFile::Open("test_DigiAndRechitPlotter_TB_4114_Ped_Profile.root");
     TFile *f7 = TFile::Open("test_DigiAndRechitPlotter_TB_4118_Ped_Profile.root");


     double time[7]={0., 16., 32., 48., 64., 80., 96.};
     double adc_mean[7]={0.}; 
     double adc_error[7]={0.}; 


     f1->cd("hgcaltbdigisplotter_ped_profile");
     TH1F* h1= (TH1F*) f1->FindObjectAny("Sum_Cluster_ADC");     
     adc_mean[0] = h1->GetMean();
     adc_error[0] = h1->GetRMS()/sqrt(512);
/*
     h1->Rebin(4);
     h1->SetLineColor(2);
     h1->Draw();
*/
     f2->cd("hgcaltbdigisplotter_ped_profile");
     TH1F* h2= (TH1F*) f2->FindObjectAny("Sum_Cluster_ADC"); 
     adc_mean[1] = h2->GetMean();
     adc_error[1] = h2->GetRMS()/sqrt(512);
     f3->cd("hgcaltbdigisplotter_ped_profile");
     TH1F* h3= (TH1F*) f3->FindObjectAny("Sum_Cluster_ADC"); 
     adc_mean[2] = h3->GetMean();
     adc_error[2] = h3->GetRMS()/sqrt(512);
     f1->cd("hgcaltbdigisplotter_ped_profile");
     TH1F* h4= (TH1F*) f4->FindObjectAny("Sum_Cluster_ADC"); 
     adc_mean[3] = h4->GetMean();
     adc_error[3] = h4->GetRMS()/sqrt(512);
     f5->cd("hgcaltbdigisplotter_ped_profile");
     TH1F* h5= (TH1F*) f5->FindObjectAny("Sum_Cluster_ADC"); 
     adc_mean[4] = h5->GetMean();
     adc_error[4] = h5->GetRMS()/sqrt(512);
/*
     h5->Rebin(4);
     h5->Draw("sames");
*/
     f6->cd("hgcaltbdigisplotter_ped_profile");
     TH1F* h6= (TH1F*) f6->FindObjectAny("Sum_Cluster_ADC"); 
     adc_mean[5] = h6->GetMean();
     adc_error[5] = h6->GetRMS()/sqrt(512);
     f1->cd("hgcaltbdigisplotter_ped_profile");
     TH1F* h7= (TH1F*) f7->FindObjectAny("Sum_Cluster_ADC"); 
     adc_mean[6] = h7->GetMean();
     adc_error[6] = h7->GetRMS()/sqrt(512);

     TGraphErrors* gr = new TGraphErrors(7,time,adc_mean,0,adc_error);
                  gr->Draw("*ACE");
                  gr->GetXaxis()->SetTitle("Delay(ns)");
                  gr->GetYaxis()->SetTitle("<Cluster(adc counts)>");

     c1->SaveAs("Timing_Max_Cluster.png");
}
