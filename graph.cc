#include <fstream>
#include <iostream>
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include <vector>
#include "TPad.h"
int graph () {


	double xb,col, error;
	std::vector< double > vxb;
	std::vector< double > vcol;
	std::vector< double > verror;
    gStyle->SetOptLogx();
//	TFile* out = new TFile("Asym.root", "RECREATE");
       	TCanvas *c1 = new TCanvas("canv", "Asymmetry", 800, 800);
/*  	TH2D *hist = new TH2D("h", "Asymmetry",5, 0.0, 0.4,5, -0.1, 0.3);
	std::fstream f;
	f.open("logpiminus");
	while (!f.eof() ) {
		f >> xb;
		f >> col;
		f >> error;
		vxb.push_back(xb);
		vcol.push_back(col);
		verror.push_back(error);

		
		cout << xb << " " << col << " " << error << endl;
		//gr->AddPoint(xb, col);	
		hist->Fill(xb, col);
	}
	c->Divide(2,2, 0.05, 0.05, 0);*/

	TGraphErrors* gr1 = new TGraphErrors("positron_collins_Kplus_27.5_1", "%lg %lg %lg");
	gr1->SetTitle("K^{+} Collins asymmetry; x;A_{Col}");
	gr1->GetYaxis()->SetTitleOffset(0.5);
	gr1->SetMarkerStyle(21);
	gr1->SetMarkerColor(4);
	gr1->Draw("AP");
    gr1->SetMinimum(-0.2);
    gr1->SetMaximum(0.2);
	c1->cd();
	c1->Print("acol_Kplus.png");
       	TCanvas *c2 = new TCanvas("canv", "Asymmetry", 800, 800);
	TGraphErrors* gr2 = new TGraphErrors("positron_collins_Kminus_27.5_1", "%lg %lg %lg");
	gr2->SetTitle("K^{-} Collins asymmetry; x;Acol");
	gr2->SetMarkerStyle(21);
	gr2->GetYaxis()->SetTitleOffset(0.5);
	gr2->SetMarkerColor(4);
	gr2->Draw("AP");
    gr2->SetMinimum(-0.2);
    gr2->SetMaximum(0.2);
	c2->cd();
	c2->Print("acol_Kminus.png");
       	TCanvas *c3 = new TCanvas("canv", "Asymmetry", 800, 800);
	TGraphErrors* gr3 = new TGraphErrors("positron_collins_piplus_27.5_1", "%lg %lg %lg");
	gr3->SetTitle("#pi^{+} Collins asymmetry; x;Acol");
	gr3->SetMarkerStyle(21);
	gr3->GetYaxis()->SetTitleOffset(0.5);
	gr3->SetMarkerColor(3);

    gr3->SetMinimum(-0.2);
    gr3->SetMaximum(0.2);
	gr3->Draw("AP");
	c3->cd();
	c3->Print("acol_piplus.png");
       	TCanvas *c4 = new TCanvas("canv", "Asymmetry", 800, 800);
	TGraphErrors* gr4 = new TGraphErrors("positron_collins_piminus_27.5_1", "%lg %lg %lg");
	gr4->SetTitle("#pi^{-} Collins asymmetry; x;Acol");
	gr4->GetYaxis()->SetTitleOffset(0.5);
	gr4->SetMarkerStyle(21);
	gr4->SetMarkerColor(3);

    gr4->SetMinimum(-0.2);
    gr4->SetMaximum(0.2);
	gr4->Draw("AP");
	c4->cd();
	c4->Print("acol_piminus.png");

    
    TCanvas *c5 = new TCanvas("canv", "Asymmetry", 800, 800);
	TGraphErrors* grpi0 = new TGraphErrors("positron_collins_pi0_27.5_1", "%lg %lg %lg");
	grpi0->SetTitle("#pi^{0} Collins asymmetry; x;Acol");
	grpi0->SetMarkerStyle(21);
	grpi0->GetYaxis()->SetTitleOffset(0.5);
	grpi0->SetMarkerColor(4);
	grpi0->Draw("AP");
    grpi0->SetMinimum(-0.2);
    grpi0->SetMaximum(0.2);
	c5->cd();
	c5->Print("acol_pi0.png");
    
    TCanvas *c6 = new TCanvas("canv", "Asymmetry", 800, 800);
	TGraphErrors* grrho0 = new TGraphErrors("positron_collins_rho0_27.5_1", "%lg %lg %lg");
	grrho0->SetTitle("#rho^{0} Collins asymmetry; x;Acol");
	grrho0->SetMarkerStyle(21);
	grrho0->GetYaxis()->SetTitleOffset(0.5);
	grrho0->SetMarkerColor(4);
	grrho0->Draw("AP");
    grrho0->SetMinimum(-0.2);
    grrho0->SetMaximum(0.2);
	c6->cd();
	c6->Print("acol_rho0.png");

    TCanvas *c7 = new TCanvas("canv", "Asymmetry", 800, 800);
	TGraphErrors* grphi = new TGraphErrors("positron_collins_phi_27.5_1", "%lg %lg %lg");
	grphi->SetTitle("#phi^{0} Collins asymmetry; x;Acol");
	grphi->SetMarkerStyle(21);
	grphi->GetYaxis()->SetTitleOffset(0.5);
	grphi->SetMarkerColor(4);
	grphi->Draw("AP");
    grphi->SetMinimum(-0.2);
    grphi->SetMaximum(0.2);
	c7->cd();
	c7->Print("acol_phi.png");


    TCanvas *c8 = new TCanvas("canv", "Asymmetry", 800, 800);
	TGraphErrors* grrhop = new TGraphErrors("positron_collins_rho+_27.5_1", "%lg %lg %lg");
	grrhop->SetTitle("#rho^{+} Collins asymmetry; x;Acol");
	grrhop->SetMarkerStyle(21);
	grrhop->GetYaxis()->SetTitleOffset(0.5);
	grrhop->SetMarkerColor(4);
	grrhop->Draw("AP");
    grrhop->SetMinimum(-0.1);
    grrhop->SetMaximum(0.2);
	c8->cd();
	c8->Print("acol_rho+.png");


    TCanvas *c9 = new TCanvas("canv", "Asymmetry", 800, 800);
	TGraphErrors* grrhom = new TGraphErrors("positron_collins_rho-_27.5_1", "%lg %lg %lg");
	grrhom->SetTitle("#rho^{-} Collins asymmetry; x;Acol");
	grrhom->SetMarkerStyle(21);
	grrhom->GetYaxis()->SetTitleOffset(0.5);
	grrhom->SetMarkerColor(4);
	grrhom->Draw("AP");
    grrhom->SetMinimum(-0.1);
    grrhom->SetMaximum(0.2);
	c9->cd();
	c9->Print("acol_rho-.png");
    //	c->Update();
    //	c->Print("asymmetry.png");
//	hist->Write();
//	out->Close();
//	f.close();
	return 0;
}
