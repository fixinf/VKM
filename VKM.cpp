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
#include "fitting.h"
#include <gsl/gsl_vector.h>

using namespace std;
double func_enermin(double x, void * params) {
	set_const C = *(set_const *) params;
	return t_E(x*0.5, x*0.5, C) / (D * x);
}
double enermin(set_const C) {
	bool debug = 0;
	if (debug){
	cout << "HI, I AM THE ENERGY MINIMIZER" << endl;
	}
	gsl_function F;
	F.function = &func_enermin;
	F.params = &C;

	const gsl_min_fminimizer_type *T;
	T = gsl_min_fminimizer_brent;
	gsl_min_fminimizer *s;
	s = gsl_min_fminimizer_alloc(T);
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
			if (debug)
			printf("Converged:\n");

		if (debug)
		printf("%5d [%.7f, %.7f] "
			"%.7f %+.7f %.7f\n", iter, a, b, m, m - m_expected, b - a);
	} while (status == GSL_CONTINUE && iter < max_iter);
	cout << "DONE" << endl;
	gsl_min_fminimizer_free(s);

	cout << m << endl;
//	cout << func_enermin(m, NULL) - m_n << endl;
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
	double rho_eq = enermin(C);
	return 9.0*pow(rho_eq,2.0) * eps_second(rho_eq, C);
}



int main(void) {
	puts("Hello World!!!");
//  set_const A(sqrt(329.7), sqrt(249.40), sqrt(68.09),0.0,0.0,0.0);// WALECKA
//	set_const A(sqrt(189.94), sqrt(90.768), sqrt(100.18), 6.3714e-3, 1.6288e-2,0.0);//MOD WALECKA
//    set_const A(sqrt(169.36), sqrt(59.055), sqrt(104.56), 0.0, 0.0,0.0); // ZM
	set_const A(sqrt(179.56), sqrt(87.600), sqrt(100.64), 7.7346e-3, 3.4462e-4,0.65);//MW(nu)
//	set_const A(sqrt(184.36), sqrt(87.600), sqrt(100.64), 5.5387e-3, 2.2976e-2,0.0);//MW(u)

//    cout << A.b << "       " << A.c << endl;
//    cout << A.U(1.0) << endl;




//	cout << A.C_o << endl;
//	cout << A.phi_n(0.5) << endl;
//	cout << p_f(0.25) << endl;
//	cout << "NP   " <<  np_eq(0.5, A) << endl;
//	cout << f_eq(0.25, 0.25, A) << endl;
//	cout << A.phi_n(f_eq(0.25,0.25, A)) << endl;
//	cout << A.diff_phi_n(10.0) << endl;
//	cout << t_E(0.25,0.25,A) << endl;
//	cout << t_E(0.25,0.25,A)/(D*0.5) - m_n << endl;
//
//	cout << 1 - f_eq(0.5 -  np_eq(0.5, A), np_eq(0.5,A), A)<<endl;
//	cout << "COMRPESS. MODULUS = " << compr_modulus(A) << endl;
////
//	cout << enermin(A) << endl;
//	cout << np_eq(6.0, A) << endl;
//	plot_Ener(A);
//	cout << "ALL DONE" << endl;
//	return EXIT_SUCCESS;

	printf("%f %f %f %f %f %f \n", A.C_s, A.C_o, A.C_r, A.b, A.c, A.z);

	gsl_vector * X = gsl_vector_alloc(5);
	gsl_vector * F = gsl_vector_alloc(5);
	gsl_vector_set(X, 0, A.C_s);
	cout << "SET 0 " << endl;
	gsl_vector_set(X, 1, A.C_o);
	cout << "SET 1 " << endl;
	gsl_vector_set(X, 2, A.C_r);
	cout << "SET 2 " << endl;
	gsl_vector_set(X, 3, A.b);
	cout << "SET 3 " << endl;
	gsl_vector_set(X, 4, A.c);
	cout << "SET 4 " << endl;

	cout << "SET 5 " << endl;
	cout << "PARAMS SET" << endl;
	fit_eq_params params = {-16.0, 0.5, 275, 28.0, 0.8, A.z};
	cout << "PASSED TO FUNC_FIT" << endl;
	func_fit(X, &params, F);
	cout << gsl_vector_get(F, 0) << endl;
	cout << gsl_vector_get(F, 1) << endl;
	cout << gsl_vector_get(F, 2) << endl;
	cout << gsl_vector_get(F, 3) << endl;
	cout << gsl_vector_get(F, 4) << endl;

	cout << "TRYING TO SOLVE" << endl;
	set_const C(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	fit(params, &C);

	printf("%f %f %f %f %f %f \n", C.C_s, C.C_o, C.C_r, C.b, C.c, C.z);
	printf("%f %f %f %f %f %f \n", C.C_s*C.C_s, C.C_o*C.C_o, C.C_r*C.C_r, C.b, C.c, C.z);

	cout << "COMRPESS. MODULUS = " << compr_modulus(C) << endl;

	double n0 = enermin(C);
	cout << n0 << endl;
	cout << np_eq(6.0, C) << endl;
	cout<<func_enermin(n0, &C) - m_n<< endl ;


	plot_Ener(C);



	cout << "DONE" << endl;

	return 0;
}
