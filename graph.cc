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



void CheckInputFile(std::string name)
{
	std::ifstream file(name);
	if (!file.is_open())
	{
		std::cout << "File " << name << " does not exist" << std::endl;
		exit(1);
	}
}

template <typename ...T>
void Print(T... args)
{
	((std::cout << args << " "), ...);
	std::cout << endl;
}

class mgraph
{	
private:
	double xmin, xmax, ymin, ymax;
	std::vector<TGraphErrors> vecgraph;
	TLegend legend = TLegend(0.7, 0.7, 0.9, 0.9);
	std::unique_ptr<TCanvas> canv;
public:

	mgraph(double minx, double miny, double maxx, double maxy, std::string canv_name)
	{
		xmin = minx;
		xmax = maxx;
		ymax = maxy;
		ymin = miny;
		canv = std::unique_ptr<TCanvas>(new TCanvas(canv_name.c_str(), "a", 1680, 1200));
		//Print(xmin, xmax, ymin, ymax);

		legend.SetLineColorAlpha(0., 0.);
		legend.SetFillColorAlpha(0., 0.);
	};

	void fast_divide(int ncolumns, int nrows)
	{
		canv->Divide( ncolumns, 1, 0., 0.);
		for (int i = 1; i <= ncolumns; i++)
		{
			canv->cd(i);
			gPad->Divide(1,nrows, 0., 0);
		}
	}

	void SetGraphStyle(TGraphErrors *gr, Color_t color, Style_t style)
	{
		gr->SetMarkerColor(color);
		gr->SetMarkerStyle(style);
		gr->SetMarkerSize(1.4);
		gr->SetLineColor(color);
	}

	void add_graph(double *x, double *y, double *x_err, double *y_err, std::string legend_entry, Color_t color, Style_t style){
		vecgraph.push_back(TGraphErrors(9,x,y, x_err, y_err));  
		SetGraphStyle(&vecgraph.back(), color, style);
		//Print(vecgraph.back().GetN());
	
		legend.AddEntry(vecgraph.back().Clone(), legend_entry.c_str(), "P");
	}

	void add_graph(double *x, double *y, double *x_err, double *y_err, Color_t color, Style_t style){
		vecgraph.push_back(TGraphErrors(9,x,y, x_err, y_err));  
		SetGraphStyle(&vecgraph.back(), color, style);
		
	
		// legend.AddEntry(vecgraph.back().Clone(), legend_entry.c_str(), "P");
	}

	void add_graph(std::string file_name, std::string legend_entry, Color_t color, Style_t style){
		CheckInputFile(file_name);

		vecgraph.push_back(TGraphErrors(file_name.c_str(), "%lg %lg %lg"));
		SetGraphStyle(&vecgraph.back(), color, style);
		legend.AddEntry(vecgraph.back().Clone(), legend_entry.c_str(), "P");

	}

	void add_graph(std::string file_name, Color_t color, Style_t style){
		CheckInputFile(file_name);

		vecgraph.push_back(TGraphErrors(file_name.c_str(), "%lg %lg %lg"));
		SetGraphStyle(&vecgraph.back(), color, style);
		//legend.AddEntry(vecgraph.back().Clone(), legend_entry.c_str(), "P");

	}
	
	

	void draw(std::string title_text, int icanv, int igpad, std::string x_label, std::string y_label){
		canv->cd(icanv);
		gPad->cd(igpad);
		if (icanv == 1){
			gPad->SetLogx();
		}
	
		TH1F* hist = gPad->DrawFrame(xmin, ymin, xmax, ymax);
		hist->GetXaxis()->SetTitle(x_label.c_str());
		hist->GetYaxis()->SetTitle(y_label.c_str());
		hist->Draw("SAME AXIS X+ Y+");
		for (TGraphErrors gr: vecgraph){
			//gr.Clone()->DrawClone("P");
			gr.Clone()->Draw("P");
		}

		hist->GetXaxis()->SetTitleSize(0.07);
		hist->GetYaxis()->SetTitleSize(0.07);
		hist->GetXaxis()->SetTitleOffset(0.8);
		//hist->GetYaxis()->SetTitleOffset(1);

		hist->GetXaxis()->SetLabelSize(0.05);
		hist->GetYaxis()->SetLabelSize(0.05);

		TLatex text = TLatex(0.01, 0.15, title_text.c_str());
		text.SetTextSize(0.1);

		legend.DrawClone();

		text.DrawClone();
		vecgraph.clear();
		legend.Clear();
	}

	void FullDraw(std::string fig_name)
	{
		canv->SaveAs(fig_name.c_str());
	}
};

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

	double xKplus[9] = {0.0194,  -0.0031, -0.0027, 0.0198, -0.0166, -0.0287, -0.0245, -0.1305, -0.0564};
	double exKplus[9] = {0.0497, 0.025, 0.0211, 0.0207, 0.0244, 0.0284, 0.0357, 0.0514, 0.0846};
	double zKplus[8] = {-0.0116, -0.0142, 0.0039, 0.0002, -0.0134, 0.0013, -0.009, -0.0546};
	double ezKplus[8] = {0.0266, 0.026, 0.0266, 0.0282, 0.0219, 0.0232, 0.0336, 0.0519};
	double pTKplus[9] = {-0.1139, 0.031, 0.0026, -0.0232, -0.0034, 0.0262, -0.0322, 0.0219, -0.0416};
	double epTKplus[9] = {0.0306, 0.0262, 0.0251, 0.0259, 0.0276, 0.0255, 0.0317, 0.0314, 0.0622};

	double xKmin[9] = {0.0557, -0.0229, 0.013, 0.0265, 0.0181, -0.0406, -0.079, -0.1907, 0.0149};
	double exKmin[9] = {0.0541, 0.0286, 0.0257, 0.026, 0.0317, 0.0394, 0.0518, 0.0822, 0.1516};
	double zKmin[8] = {-0.0066, -0.0331, 0.0325, 0.0267, -0.0462, 0.0197, 0.0751, 0.0809};
	double ezKmin[8] = {0.0305, 0.0305, 0.0315, 0.0354, 0.0283, 0.0312, 0.0515, 0.0897};
	double pTKmin[9] = {-0.0153, -0.0486, 0.0408, -0.0225, -0.0109, 0.0098, 0.0023, 0.0765, -0.0747};
	double epTKmin[9] = {0.039, 0.0328, 0.0324, 0.0328, 0.035, 0.0322, 0.0402, 0.0392, 0.0756};

	double zero[9];
	for (double &arr_val: zero){
		arr_val = 0.0;
	}


    //gStyle->SetOptLogx();
	gROOT->SetBatch(1);
	mgraph multi_canv = mgraph(0.005, -0.32, 1, 0.31, "allah");

	multi_canv.fast_divide(3, 4);
	multi_canv.add_graph(x, xpiplus, zero, expiplus, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("Collins_211", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("#pi^{+}", 1, 1, "x", "A_{Col}");
	multi_canv.add_graph(x, xpimin, zero, expimin, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("Collins_-211", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("#pi^{-}", 1, 2, "x", "A_{Col}");
	multi_canv.add_graph(x, xKplus, zero, exKplus, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("Collins_311", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("K^{+}", 1, 3, "x", "A_{Col}");
	
	multi_canv.add_graph(x, xKmin, zero, exKmin, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("Collins_-311", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("K^{-}", 1, 4, "x", "A_{Col}");

	multi_canv.add_graph(z, zpiplus, zero, ezpiplus, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("zCollins_211", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("#pi^{+}", 2, 1, "z", "A_{Col}");

	multi_canv.add_graph(z, zpimin, zero, ezpimin, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("zCollins_-211", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("#pi^{-}", 2, 2, "z", "A_{Col}");

	multi_canv.add_graph(z, zKplus, zero, ezKplus, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("zCollins_311", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("K^{+}", 2, 3, "z", "A_{Col}");

	multi_canv.add_graph(z, zKmin, zero, ezKmin, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("zCollins_-311", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("K^{-}", 2, 4, "z", "A_{Col}");

	multi_canv.add_graph(pT, pTpiplus, zero, epTpiplus, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("ptCollins_211", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("#pi^{+}", 3, 1, "p_{T}", "A_{Col}");

	multi_canv.add_graph(pT, pTpimin, zero, epTpimin, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("ptCollins_-211", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("#pi^{-}", 3, 2, "p_{T}", "A_{Col}");

	multi_canv.add_graph(pT, pTKplus, zero, epTKplus, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("ptCollins_311", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("K^{+}", 3, 3, "p_{T}", "A_{Col}");

	multi_canv.add_graph(pT, pTKmin, zero, epTKmin, "COMPASS", kRed-3, 55);
	multi_canv.add_graph("ptCollins_-311", "PYTHIA8", kAzure-3, 59);
	multi_canv.draw("K^{-}", 3, 4, "p_{T}", "A_{Col}");


	multi_canv.FullDraw("COMPASS_Collins.png");
	
	mgraph sivers = mgraph(0.005, -0.21, 1, 0.21, "Zoroaster");
	sivers.fast_divide(3, 4);

	sivers.add_graph("Expsivers_211.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("Sivers_211", "PYTHIA8", kAzure-3, 59);
	sivers.draw("#pi^{+}", 1, 1, "x", "A_{Siv}");

	sivers.add_graph("Expsivers_-211.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("Sivers_-211", "PYTHIA8", kAzure-3, 59);
	sivers.draw("#pi^{-}", 1, 2, "x", "A_{Siv}");

	sivers.add_graph("Expsivers_311.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("Sivers_311", "PYTHIA8", kAzure-3, 59);
	sivers.draw("K^{+}", 1, 3, "x", "A_{Siv}");

	sivers.add_graph("Expsivers_-311.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("Sivers_-311", "PYTHIA8", kAzure-3, 59);
	sivers.draw("K^{-}", 1, 4, "x", "A_{Siv}");

	sivers.add_graph("zExpsivers_211.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("zSivers_211", "PYTHIA8", kAzure-3, 59);
	sivers.draw("#pi^{+}", 2, 1, "z", "A_{Siv}");

	sivers.add_graph("zExpsivers_-211.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("zSivers_-211", "PYTHIA8", kAzure-3, 59);
	sivers.draw("#pi^{-}", 2, 2, "z", "A_{Siv}");

	sivers.add_graph("zExpsivers_311.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("zSivers_311", "PYTHIA8", kAzure-3, 59);
	sivers.draw("K^{+}", 2, 3, "z", "A_{Siv}");	

	sivers.add_graph("zExpsivers_-311.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("zSivers_-311", "PYTHIA8", kAzure-3, 59);
	sivers.draw("K^{-}", 2, 4, "z", "A_{Siv}");	

	sivers.add_graph("ptExpsivers_211.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("ptSivers_211", "PYTHIA8", kAzure-3, 59);
	sivers.draw("#pi^{+}", 3, 1, "p_{T}", "A_{Siv}");

	sivers.add_graph("ptExpsivers_-211.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("ptSivers_-211", "PYTHIA8", kAzure-3, 59);
	sivers.draw("#pi^{-}", 3, 2, "p_{T}", "A_{Siv}");

	sivers.add_graph("ptExpsivers_311.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("ptSivers_311", "PYTHIA8", kAzure-3, 59);
	sivers.draw("K^{+}", 3, 3, "p_{T}", "A_{Siv}");	

	sivers.add_graph("ptExpsivers_-311.log", "COMPASS", kRed-3, 55);
	sivers.add_graph("ptSivers_-311", "PYTHIA8", kAzure-3, 59);
	sivers.draw("K^{-}", 3, 4, "p_{T}", "A_{Siv}");	

	sivers.FullDraw("COMPASS_Sivers.png");

	// HERMES 

	mgraph hermes_col = mgraph(0.005, -0.21, 1, 0.21, "christ");
	hermes_col.fast_divide(3,4);

	hermes_col.add_graph("HERMES/exp_xCollins_211.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/xCollins_211", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("#pi^{+}", 1, 1, "x", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_xCollins_-211.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/xCollins_-211", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("#pi^{-}", 1, 2, "x", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_xCollins_311.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/xCollins_311", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("K^{+}", 1, 3, "x", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_xCollins_-311.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/xCollins_-311", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("K^{-}", 1, 4, "x", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_zCollins_211.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/zCollins_211", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("#pi^{+}", 2, 1, "z", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_zCollins_-211.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/zCollins_-211", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("#pi^{-}", 2, 2, "z", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_zCollins_311.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/zCollins_311", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("K^{+}", 2, 3, "z", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_zCollins_-311.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/zCollins_-311", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("K^{-}", 2, 4, "z", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_ptCollins_211.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/ptCollins_211", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("#pi^{+}", 3, 1, "p_{T}", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_ptCollins_-211.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/ptCollins_211", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("#pi^{-}", 3, 2, "p_{T}", "A_{Col}");	

	hermes_col.add_graph("HERMES/exp_ptCollins_311.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/ptCollins_311", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("K^{+}", 3, 3, "p_{T}", "A_{Col}");

	hermes_col.add_graph("HERMES/exp_ptCollins_-311.log", "HERMES", kRed-3, 55);
	hermes_col.add_graph("HERMES/ptCollins_-311", "PYTHIA8", kAzure-3, 59);
	hermes_col.draw("K^{-}", 3, 4, "p_{T}", "A_{Col}");

	hermes_col.FullDraw("HERMES_Collins.png");


	mgraph hermes_siv = mgraph(0.005, -0.21, 1, 0.21, "Thor");
	hermes_siv.fast_divide(3,4);
	hermes_siv.add_graph("HERMES/exp_xSivers_211.log", "HERMES", kRed-3, 55);
	hermes_siv.add_graph("HERMES/xSivers_211", "PYTHIA8", kAzure-3, 59);
	hermes_siv.draw("#pi^{+}", 1, 1, "x", "A_{Siv}");

	hermes_siv.add_graph("HERMES/exp_xSivers_-211.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/xSivers_211", kAzure-3, 59);
	hermes_siv.draw("#pi^{-}", 1, 2, "", "");

	hermes_siv.add_graph("HERMES/exp_xSivers_311.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/xSivers_311", kAzure-3, 59);
	hermes_siv.draw("K^{+}", 1, 3, "", "");

	hermes_siv.add_graph("HERMES/exp_xSivers_-311.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/xSivers_-311", kAzure-3, 59);
	hermes_siv.draw("K^{-}", 1, 4, "x", "");

	hermes_siv.add_graph("HERMES/exp_zSivers_211.log", "HERMES", kRed-3, 55);
	hermes_siv.add_graph("HERMES/zSivers_211", "PYTHIA8", kAzure-3, 59);
	hermes_siv.draw("#pi^{+}", 2, 1, "z", "A_{Siv}");

	hermes_siv.add_graph("HERMES/exp_zSivers_-211.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/zSivers_211", kAzure-3, 59);
	hermes_siv.draw("#pi^{-}", 2, 2, "", "");

	hermes_siv.add_graph("HERMES/exp_zSivers_311.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/zSivers_311", kAzure-3, 59);
	hermes_siv.draw("K^{+}", 2, 3, "", "");	

	hermes_siv.add_graph("HERMES/exp_zSivers_-311.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/zSivers_-311", kAzure-3, 59);
	hermes_siv.draw("K^{-}", 2, 4, "z", "");	

	hermes_siv.add_graph("HERMES/exp_ptSivers_211.log", "HERMES", kRed-3, 55);
	hermes_siv.add_graph("HERMES/ptSivers_211", "PYTHIA8", kAzure-3, 59);
	hermes_siv.draw("#pi^{+}", 3, 1, "p_{T}", "A_{Siv}");

	hermes_siv.add_graph("HERMES/exp_ptSivers_-211.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/ptSivers_-211", kAzure-3, 59);
	hermes_siv.draw("#pi^{+}", 3, 2, "", "");

	hermes_siv.add_graph("HERMES/exp_ptSivers_311.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/ptSivers_311", kAzure-3, 59);
	hermes_siv.draw("K^{+}", 3, 3, "", "");		

	hermes_siv.add_graph("HERMES/exp_ptSivers_-311.log", kRed-3, 55);
	hermes_siv.add_graph("HERMES/ptSivers_-311", kAzure-3, 59);
	hermes_siv.draw("K^{-}", 3, 4, "p_{T}", "");	



	hermes_siv.FullDraw("HERMES_Sivers.png");

 
	return 0;
}
