//
// Created by alext on 11/10/2021.
//

#include "shockwave.h"

Shockwave::Shockwave(double rho_l, double Temp_l, double avm, double Max, double gamma) {
        /*--G is the specific heat ratio
        WRITE (*,*) ' SPECIFIC HEAT RATIO = ',G
        GM1=G-1.
        GP1=G+1.
        GR=GM1/GP1
        USS=SQRT(G*BOLTZ*FTMP/AVM)
            *--USS is the speed of sound
        SMN=FVEL/USS
            *--SMN is the shock MACH NUMBER
        WRITE (*,*) ' UPSTREAM SPEED OF SOUND = ',USS
        WRITE (*,*) ' SHOCK MACH NUMBER = ',SMN
                                            *
        UDN=FND*AVM
        SM2=SMN*SMN
        A=GP1*SM2/(2.+GM1*SM2)
        WRITE (*,*) ' DENSITY RATIO = ',A
        DDN=UDN*A
        WRITE (*,*) ' DOWNSTREAM DENSITY = ',DDN
                                             *      DVEL=226.138
        DVEL=FVEL/A
        WRITE (*,*) ' DOWNSTREAM SPEED = ',DVEL
        B=2.*G*SM2/GP1-GR
        DTMP=FTMP*B/A
        WRITE (*,*) ' DOWNSTREAM TEMPERATURE = ',DTMP
         */
        double GR = (gamma-1)/(gamma+1);
        double uss = sqrt(gamma * 1.38064852e-23 * Temp_l / avm);
        Vs_L = Max*uss;
        double sm2 = Max*Max;
        double a = (gamma+1)*sm2/(2+(gamma-1)*sm2);
        Rho_R = rho_l * a;
        Vs_R = (uss*Max) / a;
        double B = 2 * gamma * sm2 / (gamma+1) - GR;
        Temp_R = Temp_l * B / a;



}

double Shockwave::get_T() const {
    return Temp_R;
}

double Shockwave::get_Rho() const {
    return Rho_R;
}

double Shockwave::get_V() const {
    return Vs_R;
}

double Shockwave::get_V_L() const {
    return Vs_L;
}
