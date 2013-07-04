
{
	gROOT->ProcessLine(".L lhcbStyle.C");
	lhcbStyle();

	double* XValues = new double[6];
	XValues[0] = 990+(1008-990)*0.5;
	XValues[1] = 1008+(1016-1008)*0.5;
	XValues[2] = 1016+(1020-1016)*0.5;
	XValues[3] = 1020+(1024-1020)*0.5;
	XValues[4] = 1024+(1032-1024)*0.5;
	XValues[5] = 1032+(1050-1032)*0.5;
	double* XErrors = new double[6];
	XErrors[0] = (1008-990)*0.5;
	XErrors[1] = (1016-1008)*0.5;
	XErrors[2] = (1020-1016)*0.5;
	XErrors[3] = (1024-1020)*0.5;
	XErrors[4] = (1032-1024)*0.5;
	XErrors[5] = (1050-1032)*0.5;

	double offset = 0.5;
	double* XValues2 = new double[6];
	XValues2[0] = 990+(1008-990)*0.5+offset;
	XValues2[1] = 1008+(1016-1008)*0.5+offset;
	XValues2[2] = 1016+(1020-1016)*0.5+offset;
	XValues2[3] = 1020+(1024-1020)*0.5+offset;
	XValues2[4] = 1024+(1032-1024)*0.5+offset;
	XValues2[5] = 1032+(1050-1032)*0.5+offset;
	double* XErrors2_l = new double[6];
	double* XErrors2_h = new double[6];
	XErrors2_l[0] = (1008-990)*0.5+offset;  XErrors2_h[0] = (1008-990)*0.5-offset;
	XErrors2_l[1] = (1016-1008)*0.5+offset; XErrors2_h[1] = (1016-1008)*0.5-offset;
	XErrors2_l[2] = (1020-1016)*0.5+offset; XErrors2_h[2] = (1020-1016)*0.5-offset;
	XErrors2_l[3] = (1024-1020)*0.5+offset; XErrors2_h[3] = (1024-1020)*0.5-offset;
	XErrors2_l[4] = (1032-1024)*0.5+offset; XErrors2_h[4] = (1032-1024)*0.5-offset;
	XErrors2_l[5] = (1050-1032)*0.5+offset; XErrors2_h[5] = (1050-1032)*0.5-offset;

	double* YValues_Phi_zero      = new double[6];
	double* YErrors_low_Phi_zero  = new double[6];
	double* YErrors_high_Phi_zero = new double[6];
	YValues_Phi_zero[0] = 1.3244;   YErrors_low_Phi_zero[0] = sqrt(0.50782*0.50782+0.007*0.007); YErrors_high_Phi_zero[0] = sqrt(0.76539*0.76539+0.007*0.007);
	YValues_Phi_zero[1] = 0.77121;  YErrors_low_Phi_zero[1] = sqrt(0.23524*0.23524+0.028*0.028); YErrors_high_Phi_zero[1] = sqrt(0.38682*0.38682+0.028*0.028);
	YValues_Phi_zero[2] = 0.4911;   YErrors_low_Phi_zero[2] = sqrt(0.46911*0.46911+0.049*0.049); YErrors_high_Phi_zero[2] = sqrt(1.3086*1.3086+0.049*0.049);
	YValues_Phi_zero[3] = -0.51623; YErrors_low_Phi_zero[3] = sqrt(0.36407*0.36407+0.025*0.025); YErrors_high_Phi_zero[3] = sqrt(0.2173*0.2173+0.025*0.025);
	YValues_Phi_zero[4] = -0.45231; YErrors_low_Phi_zero[4] = sqrt(0.25606*0.25606+0.021*0.021); YErrors_high_Phi_zero[4] = sqrt(0.18422*0.18422+0.021*0.021);
	YValues_Phi_zero[5] = -0.65291; YErrors_low_Phi_zero[5] = sqrt(0.22554*0.22554+0.020*0.020); YErrors_high_Phi_zero[5] = sqrt(0.19914*0.19914+0.020*0.020);

	double* YValues_Phi_pi      = new double[6];
	double* YErrors_Phi_pi_low  = new double[6];
	double* YErrors_Phi_pi_high = new double[6];
	YValues_Phi_pi[0] = 1.818;  YErrors_Phi_pi_low[0] = sqrt(0.76618*0.76618+0.007*0.007); YErrors_Phi_pi_high[0] = sqrt(0.50701*0.50701+0.007*0.007);
	YValues_Phi_pi[1] = 2.3706; YErrors_Phi_pi_low[1] = sqrt(0.38721*0.38721+0.028*0.028); YErrors_Phi_pi_high[1] = sqrt(0.235*0.235+0.028*0.028);
	YValues_Phi_pi[2] = 2.6506; YErrors_Phi_pi_low[2] = sqrt(1.3087*1.3087+0.049*0.049);   YErrors_Phi_pi_high[2] = sqrt(0.29762*0.29762+0.049*0.049);
	YValues_Phi_pi[3] = 3.658;  YErrors_Phi_pi_low[3] = sqrt(0.21753*0.21753+0.025*0.025); YErrors_Phi_pi_high[3] = sqrt(0.36374*0.36374+0.025*0.025);
	YValues_Phi_pi[4] = 3.594;  YErrors_Phi_pi_low[4] = sqrt(0.18432*0.18432+0.021*0.021); YErrors_Phi_pi_high[4] = sqrt(0.25597*0.25597+0.021*0.021);
	YValues_Phi_pi[5] = 3.7946; YErrors_Phi_pi_low[5] = sqrt(0.18395*0.18395+0.020*0.020); YErrors_Phi_pi_high[5] = sqrt(0.2254*0.2254+0.020*0.020);

	TCanvas* c1 = new TCanvas( "AssymPlot", "AssymPlot" );
	c1->SetTitle("");
	TGraphAsymmErrors* soln_Phi_zero = new TGraphAsymmErrors( 6, XValues2, YValues_Phi_zero, XErrors2_l, XErrors2_h, YErrors_low_Phi_zero, YErrors_high_Phi_zero );
	soln_Phi_zero->SetMarkerColor( kBlue );
	soln_Phi_zero->SetLineColor( kBlue );
	soln_Phi_zero->SetLineWidth( 2 );
	TGraphAsymmErrors* soln_Phi_pi   = new TGraphAsymmErrors( 6, XValues, YValues_Phi_pi, XErrors, XErrors, YErrors_Phi_pi_low, YErrors_Phi_pi_high );
	soln_Phi_pi->SetMarkerColor( kRed );
	soln_Phi_pi->SetMarkerStyle( 23 );
	soln_Phi_pi->SetLineColor( kRed );
	soln_Phi_pi->SetLineWidth( 2 );

	soln_Phi_zero->Draw("AP9");
	soln_Phi_pi->Draw("P9SAME");

	c1->Update();
	c1->SetTitle("");

	soln_Phi_zero->GetYaxis()->SetRangeUser(-1.5,4.5);

	soln_Phi_zero->GetXaxis()->SetTitle( "m(K^{+}K^{-})   [MeV/c^{2}]" );
	soln_Phi_zero->GetYaxis()->SetTitle( "#delta_{S} - #delta_{#perp}   [rad]" );

	// add LHCb label
	TPaveText* lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
			0.87 - gStyle->GetPadTopMargin(),
			gStyle->GetPadLeftMargin() + 0.20,
			0.95 - gStyle->GetPadTopMargin(),
			"BRNDC");
	lhcbName->AddText("LHCb");
	lhcbName->SetFillColor(0);
	lhcbName->SetTextAlign(12);
	lhcbName->SetBorderSize(0);

	lhcbName->Draw("SAME");

	c1->Update();

	c1->Print("AssymPlot.pdf");

}


