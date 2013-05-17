//============================================================================
// Name        : VKM.cpp
// Author      : Constantin Maslov
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "set_const.h"
#include <iostream>
#include "EoS.h"
#include <math.h>
#include "constants.h"
#include <fstream>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>

using namespace std;
double func_enermin(double x, void * params) {
	set_const C = *(set_const *) params;
	return E(x, C) / (D * x);
}
double enermin(set_const C) {
	cout << "HI, I AM THE ENERGY MINIMIZER" << endl;
	gsl_function F;
	F.function = &func_enermin;
	F.params = &C;

	const gsl_min_fminimizer_type *T;
	T = gsl_min_fminimizer_brent;
	gsl_min_fminimizer *s;
	s = gsl_min_fminimizer_alloc(T);
	cout << "WOW" << endl;
	gsl_min_fminimizer_set(s, &F, 0.5, 0.1, 2.0);
	double m, a, b, m_expected = 0.5;
	int status;
	int iter = 0, max_iter = 100;
	do {
		iter++;
		status = gsl_min_fminimizer_iterate(s);
		m = gsl_min_fminimizer_x_minimum(s);
		a = gsl_min_fminimizer_x_lower(s);
		b = gsl_min_fminimizer_x_upper(s);

		status = gsl_min_test_interval(a, b, 0.001, 0.0);

		if (status == GSL_SUCCESS)
			printf("Converged:\n");

		printf("%5d [%.7f, %.7f] "
			"%.7f %+.7f %.7f\n", iter, a, b, m, m - m_expected, b - a);
	} while (status == GSL_CONTINUE && iter < max_iter);
	cout << "DONE" << endl;
	gsl_min_fminimizer_free(s);

	cout << pow(m * 0.321 * 6 * pow(pi, 2.0) / g, 1.0 / 3.0) << endl;
	cout << func_enermin(m, NULL) - m_n << endl;
	return m;
}


int plot_E(set_const C){
	ofstream ofs("plot_E.dat");
	for (int i = 1; i< 1000; i++){
		double k = 0.01*i;
		ofs << k << "        " << E(k, C) << endl;
	}

	ofs.close();
	ofstream os("plot_E.plt");
	os << "plot 'plot_E.dat' w l" << endl;
	os.close();
	int res = system("gnuplot 'plot_E.plt'");

	return 0;
}
int plot_Ener(set_const C){
	ofstream ofs("plot_Ener.dat");
	for (int i = 100; i< 1100; i++){
		double k = 0.001*i;

		//OF N/N0

		ofs << k/(0.16/0.321) << "        " << t_E(0.5*k,0.5*k, C)/(D*k) - m_n << endl;
//		ofs << k/(0.16/0.321) << "        " << E(k, C)/(D*k) - m_n << endl;

		// OF K_F_0

		//ofs << pow(0.321*k*6.0*pow(pi,2)/g, 1.0/3.0) << "        " << E(k, C)/(D*k) - m_n << endl;
	}

	ofs.close();
	ofstream os("plot_Ener.plt");
	os << "plot 'plot_Ener.dat' w l" << endl << "pause -1" << endl;
	os.close();
	int res = system("gnuplot 'plot_Ener.plt'");

	return 0;
}

double eps_second(double x, set_const C){
	double h = 1e-5;
	//return (E(x+h) - 2.0*E(x) + E(x-h))/(pow(h,2));
	return (func_enermin(x+h,&C) - 2.0*func_enermin(x,&C) + func_enermin(x-h,&C))/(pow(h,2));
}


double compr_modulus(set_const C){
	double rho_eq = 0.498;
	return 9.0*pow(rho_eq,2.0) * eps_second(rho_eq, C);
}



int main(void) {
	puts("Hello World!!!");
  set_const A(sqrt(329.7), sqrt(249.40), sqrt(68.09),0.0,0.0);// WALECKA
//	set_const A(sqrt(189.94), sqrt(90.768), sqrt(100.18), 6.3714e-3, 1.6288e-2);//MOD WALECKA
//    set_const A(sqrt(169.36), sqrt(59.055), sqrt(104.56), 0.0, 0.0); // ZM
//    cout << A.b << "       " << A.c << endl;
//    cout << A.U(1.0) << endl;
	cout << A.C_o << endl;
	cout << A.phi_n(0.5) << endl;
	cout << p_f(0.25) << endl;
	cout << "NP   " <<  np_eq(0.5, A) << endl;
	cout << f_eq(0.25, 0.25, A) << endl;
	cout << A.phi_n(f_eq(0.25,0.25, A)) << endl;
	cout << A.diff_phi_n(10.0) << endl;
	cout << t_E(0.25,0.25,A) << endl;
	cout << t_E(0.25,0.25,A)/(D*0.5) - m_n << endl;
	cout << "REAL E DENSITY = " << E(0.5, A) << endl;
	cout << "REAL E = " << E(0.5, A)/(D*0.5) - m_n << endl;
	cout << 1 - f_eq(0.5 -  np_eq(0.5, A), np_eq(0.5,A), A)<<endl;
	cout << "COMRPESS. MODULUS = " << compr_modulus(A) << endl;
	plot_Ener(A);
	cout << np_eq(6.0, A) << endl;
	return EXIT_SUCCESS;
}
