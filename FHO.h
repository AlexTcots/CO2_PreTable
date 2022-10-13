//
// Created by alext on 11/20/2021.
//

#ifndef FHO_PRETABLE_FHO_H
#define FHO_PRETABLE_FHO_H
#include "vib_const.h"
#include "species.h"
#include <cmath>
#include "random.h"
#include <iostream>
#include <omp.h>
#include <vector>




double Vibr_eng( int vib_lev);// Kustova's book

double   s_func(int i,int f);//Adamovich

double Ns(int i ,int f);//Adamovich

double ave_vib_qua(double e1,double e2,double s);//Adamovich

double theta_p(double w);


double u_from_E(double E);



std::pair<double,double> Epsilon_Generator(RanPark& rd);

double gam(double theta1,double theta2,double phi1,double phi2,std::pair<double,double> e12, double y);

double Q_func(double theta_prime,double theta1,double phi1,double theta,double w,double u,double gamma);



//modify
double Compute_pro_if_deex(int n, int i, int f,const std::vector<double>& Ebins,RanPark& rd, int MS,const double THETA_mole);
double Compute_pro_if_exci(int n, int i, int f,const std::vector<double>& Ebins,RanPark& rd, int MS,const double THETA_mole);

#endif //FHO_PRETABLE_FHO_H
