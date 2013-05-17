/*
 * set_const.cpp
 *
 *  Created on: 16.05.2013
 *      Author: fixinf
 */

#include "set_const.h"
#include <math.h>
#include "constants.h"
#include <iostream>


double set_const::phi_n(double f){
	return 1.0 - f;
//	return 1.0/(1.0 + f);
}

double set_const::diff_phi_n(double f){
	double df = 0.0001;
	return (this->phi_n(f+df) - this->phi_n(f-df))/(2.0*df);
}

double set_const::eta_r(double f){
	return 1.0;
}


double set_const::eta_o(double f){
	return 1.0;
}


double set_const::eta_s(double f){
	return 1.0;
}

double set_const::U(double f){

	//std::cout << "Cs = " << C_s << "B = " << this->b << " C = " << this->c << std::endl;
	return pow(m_n,4.0)*(b * pow(f,3.0)/3.0 + c*pow(f,4.0)/4.0);
}

set_const::set_const(double C_s, double C_o, double C_r, double b, double c) {
	// TODO Auto-generated constructor stub
	set_const::C_o = C_o;
	set_const::C_r = C_r;
	set_const::C_s = C_s;
	set_const::b = b;
	set_const::c = c;
	//std::cout << "BBBBBB" << this->b <<"CCCC" << this->c << std::endl;
	// TODO использование переданной функции как метода или же нафиг
}
