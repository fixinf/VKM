/*
 * set_const.h
 *
 *  Created on: 16.05.2013
 *      Author: fixinf
 */

#ifndef SET_CONST_H_
#define SET_CONST_H_

class set_const {
public:
	set_const(double, double, double, double, double);
	double diff_phi_n(double);
	double U(double);
	double phi_n(double);
	double eta_s(double);
	double eta_o(double);
	double eta_r(double);
	double C_s;
	double C_o;
	double C_r;
	double b;
	double c;
};

#endif /* SET_CONST_H_ */
