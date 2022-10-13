//
// Created by alext on 11/20/2021.
//

#ifndef FHO_PRETABLE_VIB_CONST_H
#define FHO_PRETABLE_VIB_CONST_H

//const Planck
const double h = 6.6162e-34;
// speed of light
const double c = 2.99792458e8;//in meter
//omega_e_c for Molecule in m-1
const double omega_N2 = 235857;// in m-1
//omegaX for Molecule in m-1
const double omegaX_N2 = 1432;
//h_bar;
const double h_bar = 1.0545718e-34;

const double omega_CO2_symtric = 148000;
const double omega_CO2_bending = 52600;
const double omega_co2_asymtric = 256500;


const double omega_mole = omega_CO2_symtric;
const double omegaX_mole = 0;
#endif //FHO_PRETABLE_VIB_CONST_H
