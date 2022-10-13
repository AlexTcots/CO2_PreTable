//
// Created by alext on 11/10/2021.
//

#ifndef FHO_PRETABLE_SHOCKWAVE_H
#define FHO_PRETABLE_SHOCKWAVE_H
#pragma once
#include <cmath>
class Shockwave {
private:
   double Rho_R;
   double Temp_R;
   double Vs_L;
   double Vs_R;
public:
    Shockwave();
    Shockwave(double rho_l,double Temp_l,double avm,double Max,double gamma);
    double get_T()const;
    double get_Rho()const;
    double get_V()const;
    double get_V_L()const;

};


#endif //FHO_PRETABLE_SHOCKWAVE_H
