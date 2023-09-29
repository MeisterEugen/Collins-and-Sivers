#include <fstream>
#include <iostream>
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include <vector>
#include "TPad.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
int graph () {

	// Compass data

	double x[9] = {0.006452, 0.01054, 0.01628, 0.02533, 0.0398, 0.06279, 0.1008, 0.1608, 0.2847};
	double z[8] = {0.2237, 0.2737, 0.3237, 0.3737, 0.4448, 0.5646, 0.7164, 0.8787};
	double pT[9] = {0.1547, 0.2519, 0.3496, 0.4482, 0.5475, 0.6681, 0.817, 1.042, 1.549};

	double xpiplus[9] = {-0.0647, -0.0039, 0.0078, -0.0007, 0.0131, 0.0085, 0.0018, 0.0214, -0.0085};
	double expiplus[9] = {0.0252, 0.0124, 0.0093, 0.0079, 0.0088, 0.0111, 0.0146, 0.0218, 0.0376};
	double zpiplus[8] = {0.0013, -0.0055, 0.0011, 0.0134, 0.0089, 0.0106, -0.0048, -0.0107};
	double ezpiplus[8] = {0.0078, 0.009, 0.0106, 0.0122, 0.0109, 0.0124, 0.0179, 0.0243};
	double pTpiplus[9] = {0.0214, 0.0, -0.005, -0.0028, -0.0044, -0.0044, 0.0245, 0.0256, -0.0456};
	double epTpiplus[9] = {0.0114, 0.0095, 0.0092, 0.0099, 0.0114, 0.0114, 0.0157, 0.0173, 0.0401};

	double xpimin[9] = {0.0261, 0.013, 0.0192, -0.0015, -0.012, 0.0154, 0.0021, 0.0059, 0.0076};
	double expimin[9] = {0.026, 0.0128, 0.0098, 0.0084, 0.0095, 0.0121, 0.0164, 0.0251, 0.0453};
	double zpimin[8] = {0.014, 0.014, 0.0001, -0.0061, -0.0107, 0.0086, 0.0121, -0.0214};
	double ezpimin[8] = {0.0083, 0.0096, 0.0113, 0.0132, 0.012, 0.0138, 0.0199, 0.0254};
	double pTpimin[9] = {-0.0056, -0.0007, -0.011, 0.0142, 0.0034, 0.0136, 0.0368, 0.0143, 0.1183};
	double epTpimin[9] = {0.0123, 0.0102, 0.0099, 0.0107, 0.0122, 0.0171, 0.0185, 0.0433};

	double xKplus[9] = {};
	double zKplus[8] = {};
	double pTKplus[9] = {};

	

    gStyle->SetOptLogx();
    TCanvas *c1 = new TCanvas("canv", "Asymmetry", 800, 800);

	// c->Divide(2,2, 0.05, 0.05, 0);


	TMultiGraph* mgpiplus = new TMultiGraph();
	mgpiplus->SetTitle("#pi^{-} Collins asymmetry; x;A_{Col}");
	mgpiplus->SetMinimum(-0.2);
	mgpiplus->SetMaximum(0.2);

	TGraphErrors* grpiplus = new TGraphErrors("Collins_211", "%lg %lg %lg");
	grpiplus->GetYaxis()->SetTitleOffset(0.5);
	grpiplus->SetMarkerStyle(21);
	grpiplus->SetMarkerColor(4);


	TGraphErrors* grhpiplus = new TGraphErrors();
	grhpiplus->SetMinimum(-0.2);
	grhpiplus->SetMaximum(0.2);
    
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
