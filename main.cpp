#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <omp.h>
#include "shockwave.h"
#include "species.h"
#include "random.h"
#include "FHO.h"


double Enb(double e_min, double e_max, int Nbins, int nbin) {
    double n = static_cast<double>(nbin);
    double Nb = static_cast<double>(Nbins - 1);
    double expent = n / Nb;

    return e_min * pow(e_max / e_min, expent);
}


int main() {
    //the parameter for N2 gas
    RanPark rd(1991.1112);
    const int SEED = 891040435;
    const int MS = 1000'000;
    //std::fstream outfile;
    //outfile.open("table.txt", std::fstream::out | std::fstream::trunc);
    //outfile<<500<<'-'<<50<<'-'<<51<<'\n';

    std::fstream outfile2;
    outfile2.open("Ebins.txt", std::fstream::out | std::fstream::trunc);
    //outfile2<<500<<'\n';

     std::fstream outfile3;
    outfile3.open("Ebins_omega.txt", std::fstream::out | std::fstream::trunc);

    //std::fstream outfile4;
    //outfile4.open("table_test.txt", std::fstream::out | std::fstream::trunc);
    std::fstream outfile5;
    outfile5.open("table_dsmc_co2_sym.txt", std::fstream::out | std::fstream::trunc);

    //test_ugamma(rd);
    //test_epsilon();
    //test_rd(rd);
    double temp = 300;
    double rho_l = 6.66 / (296.8 * temp);//rho at 6.66 pa at 300K;
    double avm = 7.31e-26;// mass of molecule;
    // computing the shock-wave parameters
    Shockwave n2m5(rho_l, temp, avm, 5, 1.4);
    //T_max is the temperature behind shock wave;
    double T_max = n2m5.get_T();
    double Vs_L = n2m5.get_V_L();
    double Vs_R = n2m5.get_V();
    // Energy bins Adamovich'paper;

    //double E_min = 2.5 * BOLTZ * temp;
    //double E_max = 2.5 * BOLTZ * 10000;
    double E_min = 1.5 * h * c * 1e5;
    double E_max = 1.5 * h * c * 1e8;

    //double E_min = 1e-21;
    //double E_max = 1e-18;


    int Nbins = 500;
    // Generate energy bins;
    std::vector<double> Ebins;
    std::vector<double> E_w;
    Ebins.reserve(Nbins);
    E_w.reserve(Nbins);
    for (int i = 0; i < Nbins; ++i) {
        double e = Enb(E_min, E_max, Nbins, i);
        double w = 0.01 * e / (1.5 * h * c);// in 1\cm
        Ebins.push_back(std::move(e));
        E_w.push_back(std::move(w));

    }

    //Energy bins done;
    for (const auto &E: Ebins) {
        outfile2 << std::setprecision(10) << E << std::endl;

    }

    for (const auto &Ew: E_w) {
        outfile3 << std::setprecision(10) << Ew << std::endl;
    }
    int max_level = 50;
    int count =0;
    for (int n = 0; n < Nbins; ++n) {

        // consider max vib level is 50; i varify (0,...50)  51 levels

        //Calculating the look-up table


       for (int i = 0; i <= max_level; ++i) {

           double psum = 0;
           for (int f = 0; f <= max_level; ++f) {
               ++count;
               int s = s_func(i, f);
               if (s <= 10 && s != 0) {
                   //chose excitation or deexcitation
                   if (i > f) { //deexcitation
                       double p = Compute_pro_if_deex(n, i, f, Ebins, rd, MS,THETA_mole);
                       psum += p;
                       outfile5 << n << ' ' << i << ' ' << f << ' '
                                << std::setprecision(10)
                                << p << ' '
                                << std::endl;
                   } else if (i < f) {// excitation
                       double e1 = Vibr_eng(i);
                       double e2 = Vibr_eng(f);
                       double e_t_r = Ebins.at(n);
                       // Trans+rot should greater than evib1 -evib2;
                       if (e_t_r > std::fabs(e1 - e2)) {
                           double p = Compute_pro_if_exci(n, i, f, Ebins, rd, MS,THETA_mole);
                           psum += p;
                           outfile5 << n << ' ' << i << ' ' << f << ' '
                                    << std::setprecision(10)
                                    << p << ' '
                                    << std::endl;
                       } else {
                           double p = 0;
                           outfile5 << n << ' ' << i << ' ' << f << ' '
                                   << std::setprecision(10)
                                   << p << ' '
                                   << std::endl;
                       }
                   }


               } else {
                   outfile5 << n << ' ' << i << ' ' << f << ' '
                            << std::setprecision(10)
                            << 0.0 << ' '
                            << std::endl;
               }

              std::cout<<n<<' '<<i<<' '<<f<<' '<<"complete"<<' '<<count<<std::endl;
           }
           outfile5 << n << ' ' << i << ' ' << max_level + 1 << ' ' << std::setprecision(10) << psum << std::endl;
           std::cout<<n<<' '<<i<<' '<<max_level + 1<<' '<<"Psum complete"<<std::endl;


       }

        //testing
      /*
       int i = 0;
       int f = 1;
       double p = Compute_pro_if_exci(n,i,f,Ebins,rd,MS);

       outfile5 << n << ' ' << i << ' ' << f << ' '
                << std::setprecision(10)
                << E_w.at(n) << ' '
                << p << ' '
                << '\n';*/


    }


    return 0;

}