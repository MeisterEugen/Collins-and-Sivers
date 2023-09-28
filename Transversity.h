	// This file stores the parametrizations of
	// the transversity PDF.


	// Transversity PDF for valence u quarks in the proton.
	inline double uTransvPDF(double xB, double Q2){

		double au = 3.2;
		double bu = 1.28;
		double xhu = au * pow( xB, bu) * pow( 1. - xB, 4);
		
		return xhu;
	
	}
	
	
	// Transversity PDF for valence d quarks in the proton.
	inline double dTransvPDF(double xB, double Q2) {
				
		double ad = 4.6;
		double bd = 1.44;
		double xhd = - ad * pow( xB, bd) * pow( 1. - xB, 4);
		
		return xhd;
	
	
	}

