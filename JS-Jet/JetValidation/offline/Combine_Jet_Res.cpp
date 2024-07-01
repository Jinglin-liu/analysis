// Convert parameter into bin number
int findBin(double para, vector<double> bins){ // Parameters are cent or pt
	int bin = 0;
	int nbins = bins.size();
	for(int i = 0; i < nbins; i++){
		if(para > bins[i])bin++;
	}
	if(bin >= nbins)return 0;
	return bin;
}

unsigned long getNEvents(string fileName, string treeName){
	auto _file = TFile::Open(fileName.c_str());
	auto _tree = (TTree*)_file->Get(treeName.c_str());
	return _tree->GetEntriesFast();
}

void Combine_Jet_Res(){
    /////////////////////////////////
    // Defining Reco Pt Bins
	vector<double> ptBins = {-1, 0, 5, 15, 25, 35, 45, 70};
	int nsize_pt = 8;
	TH1D* hRecoPt = new TH1D("hRecoPt", "Reco Pt", ptBins.size()-1, ptBins.data());
    
    // Resolution histograms over pt cuts
    vector<TH1D*> h1ResPT;
	for(unsigned int iPt = 1; iPt < nsize_pt; iPt++){
        //cout<<"Defining res histograms for pt bin:"<<iPt<<", "<<ptBins[iPt-1]<<" to "<<ptBins[iPt]<<endl;
        h1ResPT.push_back(new TH1D(Form("h1ResPT_%d", iPt), Form("Jet Resolution, %.f < Reco pT < %.f ", ptBins[iPt-1], ptBins[iPt]), 30,-5,10));
        h1ResPT[iPt-1]->GetXaxis()->SetTitle("res (p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}");
	}
    //////////////////////////////////

    // Defining Centrality Bins
	vector<double> centBins = {-1, 0, 10, 30, 50, 80};
	int nsize_cent = 6;
	TH1D* hCent = new TH1D("hCent", "Centrality", centBins.size()-1, centBins.data());
	//cout<<"Centrality hist defined"<<endl;
    
    // Reco pt vs. gen pt histograms, resolution histograms over centrality cuts
    vector<TH2D*> h2Pt;
    vector<TH1D*> h1ResCent;
    vector<TH1D*> h1ResPToCent;
	for(unsigned int iCent = 1; iCent < nsize_cent; iCent++){
	//cout<<"Defining pt histograms for centrality bin:"<<iCent<<", "<<centBins[iCent-1]<<" to "<<centBins[iCent]<<endl;
        h2Pt.push_back(new TH2D(Form("h2Pt_%d", iCent), Form("RecoPt vs TruthPt, %.f < cent < %.f ", centBins[iCent-1], centBins[iCent]), 80, 0, 80, 80, 0, 80));
        h2Pt[iCent-1]->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
        h2Pt[iCent-1]->GetYaxis()->SetTitle("p_{T}^{reco} [GeV]");

        //cout<<"Defining res histograms for centrality bin:"<<iCent<<", "<<centBins[iCent-1]<<" to "<<centBins[iCent]<<endl;
        h1ResCent.push_back(new TH1D(Form("h1ResCent_%d", iCent), Form("Jet Resolution, %.f < cent < %.f ", centBins[iCent-1], centBins[iCent]), 30,-5,10));
        h1ResCent[iCent-1]->GetXaxis()->SetTitle("res (p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}");

        // Resolution histograms over pt cuts for each centrality bin
	    for(unsigned int iPt = 1; iPt < nsize_pt; iPt++){
            //cout<<"Defining res histograms for pt bin:"<<iPt<<", "<<ptBins[iPt-1]<<" to "<<ptBins[iPt]<<endl;
            h1ResPToCent.push_back(new TH1D(Form("h1ResPT_%d_Cent_%d_%d", iPt, iCent, iPt+7*(iCent-1)), Form("Jet Resolution, %.f < Reco pT < %.f, %.f < cent < %.f ", ptBins[iPt-1], ptBins[iPt], centBins[iCent-1], centBins[iCent]), 30,-5,10));
            h1ResPToCent[iPt+7*(iCent-1)-1]->GetXaxis()->SetTitle("res (p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}");
	    }
	}
    //TH2D *h2Eta = new TH2D("h2Eta", "truthEta vs eta ", 20, -1, 1, 20, -1, 1);
    //TH2D *h2Phi = new TH2D("h2Phi", "truthPhi vs Phi ", 60, -3.2, 3.2, 60, -3.2, 3.2);
    TH2D *h2_Pt = new TH2D("h2_Pt", "RecoPt vs TruthPt, All Centrality", 80, 0, 80, 80, 0, 80);
    h2_Pt->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
    h2_Pt->GetYaxis()->SetTitle("p_{T}^{reco} [GeV]");

    // Resolution histograms
    TH1D *h_res = new TH1D("h_res","Jet Resolution",30,-5,10);
    h_res->GetXaxis()->SetTitle("res (p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}");

        // res bins
    const int res_N = 30;
    Double_t res_range[res_N+1];
    for(int i = 0; i < res_N+1; i++){
        res_range[i] = i*0.5 -5;
    }
        // Resolution histograms with Centrality bins
    TH2D *h_res_cent = new TH2D("h_res_cent","Jet Resolution over Centrality bins",res_N, res_range, centBins.size()-1, centBins.data());
    h_res_cent->GetXaxis()->SetTitle("res (p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}");
    h_res_cent->GetYaxis()->SetTitle("Centrality");

        // Reconstructed pt bins
    const int reco_pt_N = 10;
    Double_t reco_pt_bins[reco_pt_N+1];
    for(int i = 0; i < reco_pt_N+1; i++){
        reco_pt_bins[i] = i*10;
    }
        // Resolution histograms with matchPt bins
    TH2D *h_res_pt = new TH2D("h_res_pt","Jet Resolution over RecoPt bins", res_N, res_range, reco_pt_N, reco_pt_bins);
    h_res_pt->GetXaxis()->SetTitle("res (p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}");
    h_res_pt->GetYaxis()->SetTitle("p_{T}^{reco} [GeV]");


    string file1 = "../macro/simdst_10M/Jet30/output/sim_TTree_Jet30_1M.root"; //type 11 directory
    string file2 = "../macro/simdst_10M/Jet10/output/sim_TTree_Jet10_1M.root"; //type 12 directory
    string tree = "T";
    vector<string> files = {file1, file2};
    vector<double> weights = {0.002718, 3.21}; //cross section: 1 type 12 10GeV, 3.210microb; 0 type 11 30GeV, 2.718nb
	for(int iFile = 0; iFile < files.size(); iFile++){
	 	auto nEvents = getNEvents(files[iFile], tree);
		weights[iFile] /= double(nEvents);
        cout<<"Number of Events: "<<nEvents<<endl;
	}
    vector<string> columnNames = {"cent", "pt", "eta", "phi", "truthPt", "truthEta", "truthPhi"}; //Parameters in fillHist

    ROOT::RDataFrame df(tree, files, columnNames); //strings together TTrees from the files into one RDataFrame

    auto makeLabels = [&](unsigned int slot, const ROOT::RDF::RSampleInfo &id) { //This function will later be used to define the file number from which the event comes
        //Static variables are initialized only once and keep their values between function calls
        static unsigned int iFile = 0; //Static variable to keep track of the file number
        static string currentTree = id.AsString(); //Static variable to keep track of the current tree
        if(currentTree != id.AsString()){
            iFile++;
            currentTree = id.AsString();
        }
        cout<<"Tree#: "<<iFile<<"; currentTree: "<<currentTree<<endl;
        return iFile;
    };

    //This function will be called every event to find the jet matches and fill the histograms
    auto fillHist = [&](const int cent, const ROOT::RVecF &pt, const ROOT::RVecF &eta, const ROOT::RVecF &phi, const ROOT::RVecF &truthPt, const ROOT::RVecF &truthEta, const ROOT::RVecF &truthPhi, const unsigned int ifile) {
        static unsigned int nEventsProcessed = 0;
	    hCent->Fill(cent, weights[ifile]);
	    int centBin = findBin(cent, centBins);
	    //cout<<cent<<" "<<centBin<<endl;
	    if(centBin == 0)return;
	
	// Discard overlapping events
        auto leadingTruthPt = ROOT::VecOps::Max(truthPt);
        if(ifile==0 && leadingTruthPt < 30) return;
        if(ifile==1 && (leadingTruthPt < 10 || leadingTruthPt > 30)) return;

        // Jet matching and filling histograms
        for(int iTruth = 0; iTruth < truthPt.size(); iTruth++){
            if(truthPt[iTruth] < 10) continue;
            int iRecoMatch = -1;
            double dRmin = 100;
            for(int iReco = 0; iReco < pt.size(); iReco++){
                //if(truthPt[j] < 10) continue;
                double dEta = eta[iReco] - truthEta[iTruth];
                double dPhi = abs(phi[iReco] - truthPhi[iTruth]);
                if(dPhi > TMath::Pi()) dPhi = 2.*TMath::Pi() - dPhi;
                double dR = sqrt(dEta*dEta + dPhi*dPhi);
                if(dR > 0.4) continue;
                if(dR < dRmin){ // Look for min dR, which is the match Reco/Truth
                    dRmin = dR;
                    iRecoMatch = iReco;
                }
            }
            if(iRecoMatch == -1) continue; // Discard the truth/gen if no match

            hRecoPt->Fill(pt[iRecoMatch], weights[ifile]);
            int ptBin = findBin(pt[iRecoMatch], ptBins);
	        if(ptBin == 0)return;

            h2_Pt->Fill(truthPt[iTruth], pt[iRecoMatch],  weights[ifile]);
            h2Pt[centBin-1]->Fill(truthPt[iTruth], pt[iRecoMatch],  weights[ifile]);
            //h2Eta->Fill(truthEta[iTruth], eta[iRecoMatch],  weights[ifile]);
            //h2Phi->Fill(truthPhi[iTruth], phi[iRecoMatch],  weights[ifile]);
        
            // Fill resolution histograms
	        double res = 0.0;
	        res = (pt[iRecoMatch] - truthPt[iTruth]) / (truthPt[iTruth]);

	        h_res->Fill(res, weights[ifile]);
            h_res_cent->Fill(res, cent, weights[ifile]);
            h_res_pt->Fill(res, pt[iRecoMatch], weights[ifile]);
            h1ResCent[centBin-1]->Fill(res, weights[ifile]);

            h1ResPT[ptBin-1]->Fill(res, weights[ifile]);

            h1ResPToCent[ptBin+7.*(centBin-1.)-1]->Fill(res, weights[ifile]);
	
	    }

        if(++nEventsProcessed%10000==0) cout<<"Processed "<<nEventsProcessed<<" events"<<endl;

    };
    columnNames.push_back("ifile");
    df.DefinePerSample("ifile", makeLabels).Foreach(fillHist, columnNames); // df.Lazy.Lazy.....Instant Action (start a event loop)

    gStyle->SetOptStat(0);
    //TCanvas *c = new TCanvas("c", "c", 800, 700);
    //gPad->SetLogz();
    //h2Pt->Draw("colz");
    //c->SaveAs("jinglin_hist_reweight.png");

    // Fit resolution histograms
    cout<<"Fitting Jet Resolution with all Centrality"<<endl;
    //cout<<"Mean:"<<h_res->GetMean()<<", RMS"<<h_res->GetRMS()<<endl;
    h_res->Fit("gaus");
    for(auto& hist : h1ResCent){
        cout<<"Fitting "<<hist->GetTitle()<<endl;
		hist->Fit("gaus");
	}
    for(auto& hist : h1ResPT){
        cout<<"Fitting "<<hist->GetTitle()<<endl;
		hist->Fit("gaus");
	}
    for(auto& hist : h1ResPToCent){
        cout<<"Fitting "<<hist->GetTitle()<<endl;
		hist->Fit("gaus");
	}

    // Get mean and sigma from fitting curve
    TF1 *fitFunction = h_res->GetFunction("gaus");
    cout<<"h_res Fitted para: mean: "<<fitFunction->GetParameter(1)<<", RMS: "<<fitFunction->GetParameter(2)<<endl;
    cout<<"h_res Fitted para: mean Error: "<<fitFunction->GetParError(1)<<", RMS Error: "<<fitFunction->GetParError(2)<<endl;

    // For each cent bin, create Tgraph of mean/sigma, vs pt
    // h1ResPToCent # loop, #/7 = cent bin, #%7 = pt bin, array of mean, sigma, array of pt, then make a TGraph
    int grCentBin;
    int grPtBin;
    double mean[nsize_cent-1][nsize_pt-2];
    double sigma[nsize_cent-1][nsize_pt-2];
    double Errmean[nsize_cent-1][nsize_pt-2];
    double Errsigma[nsize_cent-1][nsize_pt-2];
    for(unsigned int i = 0; i < (nsize_cent-1)*(nsize_pt-1); i++){ //i=0~34
        fitFunction = h1ResPToCent[i]->GetFunction("gaus");
        if (!fitFunction) continue; // skipped all PtBin=0 case
        grCentBin = i/7; // 0~6/7=0; 7~13/7=1... 0-4
        grPtBin = i%7; // 0~6%7=0~6; 7~13%7=0~6... 1-6 
        mean[grCentBin][grPtBin-1] = fitFunction->GetParameter(1);
            cout<<"Cent Bin "<<grCentBin+1<<" mean: "<<mean[grCentBin][grPtBin-1]<<endl;
        Errmean[grCentBin][grPtBin-1] = fitFunction->GetParError(1);
            cout<<"Cent Bin "<<grCentBin+1<<" mean error: "<<Errmean[grCentBin][grPtBin-1]<<endl;

        sigma[grCentBin][grPtBin-1] = fitFunction->GetParameter(2);
            cout<<"Cent Bin "<<grCentBin+1<<" sigma: "<<sigma[grCentBin][grPtBin-1]<<endl;
        Errsigma[grCentBin][grPtBin-1] = fitFunction->GetParError(2);
            cout<<"Cent Bin "<<grCentBin+1<<" sigma error: "<<Errsigma[grCentBin][grPtBin-1]<<endl;
    }

    const int n = nsize_pt-2;
    double grPt[nsize_cent-1][nsize_pt-2];
    double grPtRg[nsize_cent-1][nsize_pt-2];
    for(unsigned int i=0; i < nsize_cent-1; i++){
        for(unsigned int j=0; j < nsize_pt-2; j++){
            grPt[i][j]=(ptBins[j+1]+ptBins[j+2])/2.;
            //cout<<"Pt bins are: "<<grPt[i][j]<<endl;
            grPtRg[i][j]=(ptBins[j+2]-ptBins[j+1])/2.;
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 2400, 1200);
    c1->Divide(4, 2);
    // Graph with error bars
    vector<TGraphErrors*> meangr;
    vector<TGraphErrors*> sigmagr;
	for(unsigned int iCent = 0; iCent < nsize_cent-1; iCent++){
        meangr.push_back(new TGraphErrors(n,grPt[iCent],mean[iCent],grPtRg[iCent],Errmean[iCent]));
        meangr[iCent]->SetTitle(Form("Res Mean vs Pt, %.f < cent < %.f; p_{T}^{reco} [GeV]; Mean", centBins[iCent], centBins[iCent+1]));
        sigmagr.push_back(new TGraphErrors(n,grPt[iCent],sigma[iCent],grPtRg[iCent],Errsigma[iCent]));
        sigmagr[iCent]->SetTitle(Form("Res Sigma vs Pt, %.f < cent < %.f; p_{T}^{reco} [GeV]; Sigma", centBins[iCent], centBins[iCent+1]));

        if(iCent == 0)continue;

        c1->cd(iCent);
        meangr[iCent]->SetMarkerStyle(8);
        meangr[iCent]->Draw("AP");
        //c1->SaveAs(Form("MeanVsPtCent%.f_to_%.f.png", centBins[iCent], centBins[iCent+1]));

        c1->cd(iCent+4);
        sigmagr[iCent]->SetMarkerStyle(8);
        sigmagr[iCent]->Draw("AP");
        //c1->SaveAs(Form("SigmaVsPtCent%.f_to_%.f.png", centBins[iCent], centBins[iCent+1]));
    }
    c1->SaveAs("MeanSigmaVsPtCent.png");
    
/*
    TCanvas *c1 = new TCanvas("c1", "c1", 1300, 700);
	for(unsigned int iCent = 1; iCent < nsize_cent; iCent++){
        meangr[iCent-1]->Draw("AB");
        sigmagr[iCent-1]->Draw("AB");
	}
*/

/*
    TCanvas *c1 = new TCanvas("c1", "c1", 1300, 700);
    c1->Divide(2, 1);
    c1->cd(1);
    gPad->SetLogz();
    h2Eta->Draw("colz");
    c1->cd(2);
    gPad->SetLogz();
    h2Phi->Draw("colz");
*/

	TFile* outf = TFile::Open("output.root", "RECREATE");
	hCent->Write();
    hRecoPt->Write();
    h2_Pt->Write();
	for(auto& hist : h2Pt){
		hist->Write();
	}
	h_res->Write();
    h_res_cent->Write();
    h_res_pt->Write();
    for(auto& hist : h1ResCent){
		hist->Write();
	}
    for(auto& hist : h1ResPT){
		hist->Write();
	}
    for(auto& hist : h1ResPToCent){
		hist->Write();
	}
    for(auto& gr : meangr){
		gr->Write();
	}
    for(auto& gr : sigmagr){
		gr->Write();
	}
    c1->Write();


}

