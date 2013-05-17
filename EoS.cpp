/*
 * EoS.cpp
 *	Provides EoS described in arXiv:nucl-th/0410063 v1 Pt.2
 *
 *
 *
 *  Created on: 16.05.2013
 *      Author: const.maslov@gmail.com
 */
#include <math.h>
#include "constants.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "set_const.h"
#include <iostream>

#include <gsl/gsl_integration.h>

using namespace std;
double p_f(double n) {
	return pow(3.0 * pi * pi * D * n, 1.0 / 3.0);
}

struct func_f_eq_params {
	double nn, np;
	set_const C;
};

double intf_func_f_eq(double p, void * params) {
	double m_eff = *(double *) params;
	return p * p / sqrt(m_eff * m_eff + p * p);
}

double func_f_eq(double f, void * params) {
	bool debug = false;
	struct func_f_eq_params *p = (struct func_f_eq_params *) params;
	double nn = p->nn;
	double np = p->np;
	set_const C = p->C;
	if (debug) {
		cout << "func_f : f = " << f << endl;
	}
	double res = pow(m_n / C.C_s, 2.0) * f * C.eta_s(f);
	if (debug) {
		cout << "func_f : 1 = " << pow(m_n / C.C_s, 2.0) * f * C.eta_s(f)
				<< endl;
	}
	double dx = 1.e-7;
	res += (C.U(f + dx) - C.U(f - dx)) / (2.0*dx * m_n * m_n);

	//Integration

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	double result, error;
	gsl_function F;
	F.function = &intf_func_f_eq;
	double m_eff = m_n * C.phi_n(f);
	F.params = &m_eff;

	//Over the protons' Fermi-sea
	gsl_integration_qags(&F, 0.0, p_f(np), 0, 1e-10, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	if (debug) {
		cout << "func_f : 2 : result = " << result << endl;
	}
	res += result * C.phi_n(f) * C.diff_phi_n(f) / (pi * pi);
	if (debug) {
		cout << "func_f : 2 = " << result * C.phi_n(f) * C.diff_phi_n(f) / (pi
				* pi) << endl;
	}
	result = 0.0;
	error = 0.0;

	//Over the neutrons' Fermi-sea
	w = gsl_integration_workspace_alloc(1000);
	gsl_integration_qags(&F, 0.0, p_f(nn), 0, 1e-10, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);
	res += result * C.phi_n(f) * C.diff_phi_n(f) / (pi * pi);
	if (debug) {
		cout << "func_f : 3: result = " << result << endl;
		cout << "func_f : 2 = " << result * C.phi_n(f) * C.diff_phi_n(f) / (pi
				* pi) << endl;
	}

	//Additions from \eta_i
	dx = 1e-7;
	double d_eta_o = (C.eta_o(f+dx) - C.eta_o(f-dx))/(2.0*dx);
	double d_eta_r = (C.eta_r(f+dx) - C.eta_r(f-dx))/(2.0*dx);

	res -= 0.5*pow(C.C_o * (nn + np) / m_n * C.eta_o(f),2.0) * d_eta_o;
	res -= pow(C.C_r * (nn - np) / m_n * C.eta_r(f),2.0) * d_eta_r / 8.0;

	return res;
}

double f_eq(double nn, double np, set_const C) {
	int status;
	int iter = 0, max_iter = 300;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	double x = 0.0, x_expect = 0.5;
	double xmin = 0.0, xmax = 2.0;

	gsl_function F;
	F.function = &func_f_eq;
	func_f_eq_params params = { nn, np, C };
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, xmin, xmax);

	do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		x = gsl_root_fsolver_root(s);
		xmin = gsl_root_fsolver_x_lower(s);
		xmax = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(xmin, xmax, 0, 1e-10);
		//		if (status == GSL_SUCCESS)
		//		printf("Hi, i'm done");
	} while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free(s);
	if (status != 0) {
		cout << "error, status = " << status << endl;
	}
	return x;
}

//Eq. (55.1) for \mu_e
double mu_e(double n, double np, double f, set_const C) {
	double res = pow(C.C_r, 2.0) * D *(n - 2.0 * np) / (2.0 * m_n *m_n* C.eta_r(f));
	res -= sqrt(pow(m_n * C.phi_n(f), 2.0) + pow(p_f(np), 2.0));
	res += sqrt(pow(m_n * C.phi_n(f), 2.0) + pow(p_f(n - np), 2.0));
	return res;
}

double theta(double x) {
	if (x > 0.0) {
		return 1.0;
	} else {
		return 0.0;
	}
}
struct func_np_params {
	double n;
	set_const C;
};

//Eq(55.2) for n_p
double func_np(double np, void * params) {
	struct func_np_params *p = (struct func_np_params *) params;
	double n = p->n;
	set_const C = p->C;
	double f = f_eq(n-np, np, C);
	double mue = mu_e(n, np, f,C);
	//cout << "MUE = " << mue << " N = "<< n << " NP = " << np << " F = " << f << endl;
	double result = D*np;
	if (mue*mue - m_e*m_e >= 0){
		result -= pow(mue*mue - m_e*m_e, 3.0/2.0)/(3.0*pi*pi);
	}
	if (mue*mue - m_mu*m_mu >= 0){
		result -= pow(mue*mue - m_mu*m_mu, 3.0/2.0)/(3.0*pi*pi);
	}
	//cout << "RESULT = " << result << endl;
	return result;
}

double np_eq(double n, set_const C) {
	int status;
	int iter = 0, max_iter = 300;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double x = 0.0, x_expect = 0.1;
	double xmin = 0.000001, xmax = 0.5 * n;
	gsl_function F;
	F.function = &func_np;
	func_np_params par = { n, C };
	//std::cout<<rho_s<<std::endl;
	F.params = &par;
	//gsl_set_error_handler_off();

	T = gsl_root_fsolver_bisection;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, xmin, xmax);

	do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		x = gsl_root_fsolver_root(s);
		xmin = gsl_root_fsolver_x_lower(s);
		xmax = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(xmin, xmax, 0, 1e-10);
		//		if (status == GSL_SUCCESS)
		//		printf("Hi, i'm done");
	} while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free(s);
	if (status != 0) {
		cout << "error, status = " << status << endl;
	}
	return x;
}

double func_e(double p, void * params) {
	double m_eff = *(double *) params;
	return p * p * sqrt(p * p + m_eff * m_eff);
}

double t_E(double nn, double np, set_const C) {
	double f = f_eq(nn, np, C);
	double res = 0.5 * pow(m_n * m_n * f / C.C_s, 2.0);
	res += C.U(f);
	//Integration
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

	double result, error;
	gsl_function F;
	F.function = &func_e;
	double m_eff = m_n * C.phi_n(f);
	F.params = &m_eff;

	//Over the protons' Fermi-sea
	gsl_integration_qags(&F, 0.0, p_f(np), 0, 1e-10, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	res += result / (pi * pi);

	result = 0;
	error = 0;
	//Neutrons'

	w = gsl_integration_workspace_alloc(1000);
	gsl_integration_qags(&F, 0.0, p_f(nn), 0, 1e-10, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	res += result / (pi * pi);

	//Omega

	res += 0.5 * pow(C.C_o * D * (nn + np) / m_n, 2.0) / (C.eta_o(f));

	//Rho

	res += pow(C.C_r * D * (nn - np) / m_n, 2.0) / (8.0 * C.eta_r(f));

	//TODO electrons & muons
	w = gsl_integration_workspace_alloc(1000);
	double me = m_e;
	F.params = &me;
	double pf_e = 0;
	if (pow(mu_e(nn+np, np, f, C),2.0) - me*me >= 0){
		pf_e = sqrt(pow(mu_e(nn+np, np, f, C),2.0) - me*me );
	}
//	cout << "MU_E = " <<  mu_e(nn+np, np, f, C) << endl;
	gsl_integration_qags(&F, 0.0, pf_e, 0, 1e-10, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	res += result/(pi*pi);
//	cout << "RESULT   " << result << endl;
	w = gsl_integration_workspace_alloc(1000);
	double mmu = m_mu;
	F.params = &mmu;
	double pf_mu =0;
	if (pow(mu_e(nn+np, np, f, C),2.0) - mmu*mmu >= 0){
		pf_mu = sqrt(pow(mu_e(nn+np, np, f, C),2.0) - mmu*mmu );
	}
	gsl_integration_qags(&F, 0.0, pf_mu, 0, 1e-10, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	res += result / (pi*pi);



	return res;
}

double t_P(double nn, double np, set_const C) {
	return 0.0;
}

double E(double n, set_const C) {
	double np = np_eq(n, C);
	return t_E(n-np, np, C);
}

double P(double n, set_const C) {
	return 0.0;
}
