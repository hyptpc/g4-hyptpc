/*parameters
par[0]:dark current probability
par[1]:pedestal peak sigma
par[2]:pedestal peak mean
par[3]:coefficient of the exponentially decay of the dark current

par[4]:lambda, mean value of NPE
par[5]:single photon peak(first peak except for pedestal) sigma
par[6]:number of events
par[7]:average charge at the PM output when one photoelectron (difference of ADC channel between a peak and next peak)
par[8]:effecive spectrum shift due to background
*/


double Gain_Fitting(double *x, double *par){
  // Background        pedestal peak + dark current
  double model = (1-par[0])*TMath::Gaus(x[0], par[2], par[1], kTRUE)+par[0]*TMath::Floor(x[0]-par[2])*par[3]*TMath::Exp((-par[3])*(x[0]-par[2]));
  // We sum convolution of Poisson and Gaussian
  double D=3; //number of peaks to consider (except for pedestal)
  for(int n=1;n<D+1;n++){
    model += ((TMath::Power(par[4],n))/(TMath::Factorial(n)))*(TMath::Gaus(x[0], par[2]/*+par[8]*/+n*par[7], (TMath::Sqrt(n))*par[5], kTRUE));
  }
  return par[6]*(TMath::Exp(-par[4]))*model;
}

double Gain_Fitting_1(double *x, double *par){
  // Background        pedestal peak + dark current
  double model = (1-par[0])*TMath::Gaus(x[0], par[2], par[1], kTRUE)+par[0]*TMath::Floor(x[0]-par[2])*par[3]*TMath::Exp((-par[3])*(x[0]-par[2]));
  // We sum convolution of Poisson and Gaussian
  double D=3; //number of peaks to consider (except for pedestal)
  for(int n=1;n<D+1;n++){
    model += ((TMath::Power(par[4],n))/(TMath::Factorial(n)))*(TMath::Gaus(x[0], par[2]/*+par[8]*/+n*par[7], (TMath::Sqrt(n))*par[5], kTRUE));
  }
  return par[6]*(TMath::Exp(-par[4]))*model;
}



void gain_check(){

  //Pedestal
  Int_t N=4;
  TFile *file_pe[N];
  TTree *data_pe[N];
  Double_t ADC_pe[N][16];
  Double_t total_pe[N];
  Double_t to_min=1000000;
  TH1D* hist_pe[N][16];
  TF1* fit_pe[N][16];
  
  for(int i=0;i<N;i++){
    file_pe[i] = new TFile(Form("LED_modified/BAC%d_pedestal.root",i+1),"read");
    data_pe[i] = (TTree*)file_pe[i]->Get("tree");
    data_pe[i]->SetBranchAddress("ADC",ADC_pe[i]);
    total_pe[i] = data_pe[i]->GetEntries();
    if(total_pe[i]<to_min)to_min = total_pe[i];
    for(int j=0;j<16;j++){
      hist_pe[i][j] = new TH1D(Form("hist_pe_BAC%d_Ch%d",i+1,j),Form("hist_pe_BAC%d_Ch%d",i+1,j),160,90,250);
      fit_pe[i][j] = new TF1(Form("fit_pe_BAC%d_Ch%d",i+1,j),"gaus(0)",90,250);
    }
  }



  for(int n=0;n<to_min;n++){
    for(int i=0;i<N;i++){
      data_pe[i]->GetEntry(n);
      for(int j=0;j<16;j++){
	hist_pe[i][j]->Fill(ADC_pe[i][j]);
      }
    }
  }

  Double_t parameter_pe[N][16][3];
  TCanvas *c_pe[N];

  for(int i=0;i<N;i++){
    c_pe[i] = new TCanvas(Form("c_pe_BAC%d",i+1),Form("c_pe_BAC%d",i+1),800,650);
    c_pe[i]->Divide(4,4);
    for(int j=0;j<16;j++){
      c_pe[i]->cd(j+1);
      hist_pe[i][j]->Fit(fit_pe[i][j],"Q","",90,250);
      fit_pe[i][j]->GetParameters(parameter_pe[i][j]);
    }
  }

  //Gain for each channel
  TFile *file[N];
  TTree *data[N];
  Double_t ADC[N][16];
  Double_t total[N];
  TH1D* hist[N][16];
  TF1* fit[N][16];
  to_min = 100000;

  for(int i=0;i<N;i++){
    file[i] = new TFile(Form("LED_modified/BAC%d_LED.root",i+1),"read");
    data[i] = (TTree*)file[i]->Get("tree");
    data[i]->SetBranchAddress("ADC",ADC[i]);
    total[i] = data[i]->GetEntries();
    if(total[i]<to_min)to_min = total[i];
    for(int j=0;j<16;j++){
      hist[i][j] = new TH1D(Form("hist_BAC%d_Ch%d",i+1,j),Form("hist_BAC%d_Ch%d",i+1,j),310,90,400);
      fit[i][j] = new TF1(Form("fit_BAC%d_Ch%d",i+1,j),Gain_Fitting,90,400,8);

      //best condition

      fit[i][j]->SetParLimits(0,0.1,0.8);
      fit[i][j]->SetParLimits(1,1,parameter_pe[i][j][2]+2);
      fit[i][j]->SetParLimits(2,parameter_pe[i][j][1]-5,parameter_pe[i][j][1]+3);
      fit[i][j]->SetParLimits(3,0,1);

      fit[i][j]->SetParLimits(5,0.1,5);
      fit[i][j]->SetParLimits(7,13,18.5);


      
    }
  }


  Double_t mean[N][16];

  for(int n=0;n<to_min;n++){
    for(int i=0;i<N;i++){
      data[i]->GetEntry(n);
      for(int j=0;j<16;j++){
	hist[i][j]->Fill(ADC[i][j]);
      }
    }
  }


  TCanvas *c_gain[N];

  for(int i=0;i<N;i++){
    c_gain[i] = new TCanvas(Form("c_gain_BAC%d",i+1),Form("c_gain_BAC%d",i+1),800,650);
    c_gain[i]->Divide(4,4);
    for(int j=0;j<16;j++){
      c_gain[i]->cd(j+1);
      mean[i][j] = hist[i][j]->GetMean();
      //std::cout<<mean[i][j]<<std::endl;
      fit[i][j]->SetParLimits(4,(mean[i][j]-parameter_pe[i][j][1])/25-2,(mean[i][j]-parameter_pe[i][j][1])/25+2);
      fit[i][j]->SetParLimits(6,to_min-1500,to_min+1500);
      hist[i][j]->Fit(fit[i][j],"","",parameter_pe[i][j][1]-20,parameter_pe[i][j][1]+300);
      //fit_pe[i][j]->GetParameters(parameter_pe[i][j]);
    }
  }

  Double_t channel[N][16];
  Double_t error_ch[N][16];
  Double_t parameter_gain[N][16];
  Double_t parameter_fps[N][16];
  Double_t error_gain[N][16];
  TGraphErrors* gain[N];
  //Double_t plus = 0;
  Double_t average[N];
  for(int i=0;i<N;i++){
    for(int j=0;j<16;j++){
      parameter_gain[i][j] = fit[i][j]->GetParameter(7);
      parameter_fps[i][j] = fit[i][j]->GetParameter(5);
      average[i]+=parameter_gain[i][j];
      error_gain[i][j] = fit[i][j]->GetParError(7);
      channel[i][j] = 0.0+j;
      error_ch[i][j] = 0.;
    }
    gain[i] = new TGraphErrors(16,channel[i],parameter_gain[i],error_ch[i],error_gain[i]);
    //gain[i]->SetMarkerStyle(21+i);
    gain[i]->SetMarkerStyle(24);
    gain[i]->SetMarkerColor(1+i);
    gain[i]->SetLineColor(1+i);
    gain[i]->SetMarkerSize(0.7);
    average[i] = average[i]/16;
    std::cout<<"BAC"<<i+1<<" gain : "<<average[i]<<std::endl;
  }

  TMultiGraph *mg = new TMultiGraph();
  TLegend *le1 = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<N;i++){
    mg->Add(gain[i]);
    le1->AddEntry(gain[i],Form("BAC%d",i+1));
    
  }

  TCanvas *c_result = new TCanvas("c_result","c_result",1000,800);
  c_result->cd();
  mg->SetTitle("Gain for each MPPC cell;channel;Gain [ADC Ch.]");
  mg->Draw("AP");
  le1->Draw();

  
  

  
  

  //ELPH LED
  TFile *file_elph_pe[N];
  TTree *data_elph_pe[N];
  Double_t ADC_elph_pe[N][3];
  Double_t total_elph_pe[N];
  to_min=1000000;
  TH1D* hist_elph_pe[N][3];
  TF1* fit_elph_pe[N][3];
  
  for(int i=0;i<N;i++){
    file_elph_pe[i] = new TFile(Form("LED_modified/BAC%d_ELPH_pedestal.root",i+1),"read");
    data_elph_pe[i] = (TTree*)file_elph_pe[i]->Get("tree");
    data_elph_pe[i]->SetBranchAddress("ADC",ADC_elph_pe[i]);
    total_elph_pe[i] = data_elph_pe[i]->GetEntries();
    if(total_elph_pe[i]<to_min)to_min = total_elph_pe[i];
    for(int j=0;j<3;j++){
      hist_elph_pe[i][j] = new TH1D(Form("hist_elph_pe_elph_BAC%d_Ch%d",i+1,j),Form("hist_elph_pe_elph_BAC%d_Ch%d",i+1,j+14),200,50,250);
      fit_elph_pe[i][j] = new TF1(Form("fit_elph_pe_elph_BAC%d_Ch%d",i+1,j+14),"gaus(0)",50,250);
    }
  }



  for(int n=0;n<to_min;n++){
    for(int i=0;i<N;i++){
      data_elph_pe[i]->GetEntry(n);
      for(int j=0;j<3;j++){
	hist_elph_pe[i][j]->Fill(ADC_elph_pe[i][j]);
      }
    }
  }

  Double_t parameter_elph_pe[N][3][3];
  TCanvas *c_elph_pe[N];

  for(int i=0;i<N;i++){
    c_elph_pe[i] = new TCanvas(Form("c_elph_pe_BAC%d",i+1),Form("c_elph_pe_BAC%d",i+1),800,650);
    c_elph_pe[i]->Divide(3);
    for(int j=0;j<3;j++){
      c_elph_pe[i]->cd(j+1);
      hist_elph_pe[i][j]->Fit(fit_elph_pe[i][j],"Q","",50,250);
      fit_elph_pe[i][j]->GetParameters(parameter_elph_pe[i][j]);
    }
  }

  //Gain for each channel
  TFile *file_elph[N];
  TTree *data_elph[N];
  Double_t ADC_elph[N][3];
  Double_t total_elph[N];
  TH1D* hist_elph[N][3];
  TF1* fit_elph[N][3];
  to_min = 100000;

  for(int i=0;i<N;i++){
    file_elph[i] = new TFile(Form("LED_modified/BAC%d_ELPH_LED.root",i+1),"read");
    data_elph[i] = (TTree*)file_elph[i]->Get("tree");
    data_elph[i]->SetBranchAddress("ADC",ADC_elph[i]);
    total_elph[i] = data_elph[i]->GetEntries();
    if(total_elph[i]<to_min)to_min = total_elph[i];
    for(int j=0;j<3;j++){
      hist_elph[i][j] = new TH1D(Form("hist_BAC%d_Ch%d",i+1,j),Form("hist_BAC%d_Ch%d",i+1,j+14),200,50,250);
      fit_elph[i][j] = new TF1(Form("fit_BAC%d_Ch%d",i+1,j+14),Gain_Fitting_1,50,250,8);

      //best condition


      fit_elph[i][j]->SetParLimits(0,0.0001,1);
      fit_elph[i][j]->SetParLimits(1,1,parameter_elph_pe[i][j][2]+3);
      fit_elph[i][j]->SetParLimits(2,parameter_elph_pe[i][j][1]-2,parameter_elph_pe[i][j][1]+2);
      fit_elph[i][j]->SetParLimits(3,0,2);
	
      fit_elph[i][j]->SetParLimits(5,1,parameter_fps[i][j+13]+1);
      fit_elph[i][j]->SetParLimits(7,parameter_gain[i][j+13]*0.8-3,parameter_gain[i][j+13]);
      



      
    }
  }


  Double_t mean_elph[N][3];

  for(int n=0;n<to_min;n++){
    for(int i=0;i<N;i++){
      data_elph[i]->GetEntry(n);
      for(int j=0;j<3;j++){
	hist_elph[i][j]->Fill(ADC_elph[i][j]);
      }
    }
  }


  TCanvas *c_elph_gain[N];

  for(int i=0;i<N;i++){
    c_elph_gain[i] = new TCanvas(Form("c_elph_gain_BAC%d",i+1),Form("c_elph_gain_BAC%d",i+1),800,650);
    c_elph_gain[i]->Divide(3);
    for(int j=0;j<3;j++){
      c_elph_gain[i]->cd(j+1);
      mean_elph[i][j] = hist_elph[i][j]->GetMean();
      //std::cout<<mean_elph[i][j]<<std::endl;
      fit_elph[i][j]->SetParLimits(4,(mean_elph[i][j]-parameter_elph_pe[i][j][1])/15-3,(mean_elph[i][j]-parameter_elph_pe[i][j][1])/15+3);
      fit_elph[i][j]->SetParLimits(6,to_min-2000,to_min+2000);
      hist_elph[i][j]->Fit(fit_elph[i][j],"","",parameter_elph_pe[i][j][1]-20,parameter_elph_pe[i][j][1]+300);
      //fit_elph_pe[i][j]->GetParameters(parameter_elph_pe[i][j]);
    }
  }

  Double_t channel_elph[N][3];
  Double_t error_ch_elph[N][3];
  Double_t parameter_gain_elph[N][3];
  Double_t error_gain_elph[N][3];
  TGraphErrors* gain_elph[N];
  TGraph* compare[N];
  Double_t ratio[N][3];
  for(int i=0;i<N;i++){
    for(int j=0;j<3;j++){
      parameter_gain_elph[i][j] = fit_elph[i][j]->GetParameter(7);
      std::cout<<parameter_gain_elph[i][j]<<std::endl;
      ratio[i][j] = parameter_gain_elph[i][j]/parameter_gain[i][13+j];

      error_gain_elph[i][j] = fit[i][j]->GetParError(7);
      channel_elph[i][j] = 14.0+j;
      error_ch_elph[i][j] = 0.;
    }
    gain_elph[i] = new TGraphErrors(3,channel_elph[i],parameter_gain_elph[i],error_ch_elph[i],error_gain_elph[i]);
    compare[i] = new TGraph(3,channel_elph[i],ratio[i]);
    //gain_elph[i]->SetMarkerStyle(21+i);
    gain_elph[i]->SetMarkerStyle(24);
    gain_elph[i]->SetMarkerColor(1+i);
    gain_elph[i]->SetLineColor(1+i);
    gain_elph[i]->SetMarkerSize(0.7);
    compare[i]->SetMarkerStyle(24);
    compare[i]->SetMarkerColor(1+i);
    compare[i]->SetLineColor(1+i);
    compare[i]->SetMarkerSize(0.7);
  }

  Double_t ave_ratio = 0;
  Double_t gain_ind[N];
  for(int i=0;i<N;i++){
    gain_ind[i] = 0;
    for(int j=0;j<3;j++){
      ave_ratio +=ratio[i][j];
      gain_ind[i]+=ratio[i][j];
    }
    std::cout<<"Gain average of BAC"<<i+1<<" is "<<gain_ind[i]/3<<std::endl;
  }
  std::cout<<"Gain Average is "<<ave_ratio/(N*3)<<std::endl;


  TMultiGraph *mg_elph = new TMultiGraph();
  TMultiGraph *mg_com = new TMultiGraph();
  TLegend *le2 = new TLegend(0.8,0.5,0.48,0.6);
  TLegend *le3 = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<N;i++){
    mg_elph->Add(gain_elph[i]);
    mg_com->Add(compare[i]);
    le2->AddEntry(gain_elph[i],Form("BAC%d",i+1));
    le3->AddEntry(compare[i],Form("BAC%d",i+1));
    
  }

  TCanvas *c_result_elph = new TCanvas("c_result_elph","c_result_elph",1000,800);
  c_result_elph->Divide(2);
  c_result_elph->cd(1);
  mg_elph->SetTitle("Gain for each MPPC cell;channel;Gain [ADC Ch.]");
  mg_elph->Draw("AP");
  le2->Draw();

  c_result_elph->cd(2);
  mg_com->SetTitle("Gain comparison for each MPPC cell;channel;Gain [ADC Ch.]");
  mg_com->Draw("AP");
  le3->Draw();



 
  
      
      
    

  
  
  
  
}
