//
// Created by alext on 11/20/2021.
//

#ifndef FHO_PRETABLE_SPECIES_H
#define FHO_PRETABLE_SPECIES_H

 const double Pi = 3.14159265359;
 const double BOLTZ = 1.38064852e-23;
 const double THETA_N2 = 3395;
 const double m_N = 2.3258671e-26;
 const double m_N2 = 4.65e-26;
 const double m_r_o_N2 = 1.1629336e-26;
 const double m_r_N2  = 2.3258671e-26;
 const double Xi_N2 = 1.0;
// from Adamovich's paper;
 const double d_n2 = 4.467e-10;
 const double Morse_N2 = 3.7e10;
 const double CAP = 5.0;
 const double ETA = 0.25;
const double ETA_HS = 0.25;
// From Kustova's paper Calculation of vibrational relaxation times
//in carbon dioxide using forced harmonic
//oscillator model
 const double d_co2 = 5.54e-10;
 const double Morse_CO2 = 4.3e10;
 const double m_r_CO2 = 3.655e-26;
 const double m_r_o_CO2 =0.725e-26;
 const double Xi_CO2 = 2.5;
 const double THETA_CO2_symetric = 1918.6;
 const double THETA_CO2_bending = 959;
 const double THETA_CO2_asymtric = 3382.6;





 // define
 const double THETA_mole = THETA_CO2_symetric;

const double m_r_mole = m_r_CO2;
const double Morse_mole = Morse_CO2;
const double Xi_mole = Xi_CO2;


#endif //FHO_PRETABLE_SPECIES_H
