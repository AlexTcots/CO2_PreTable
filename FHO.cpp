//
// Created by alext on 11/20/2021.
//

#include "FHO.h"


double Vibr_eng(int i) {
    //return h * c * (omega_mole * (i + 0.5) - omegaX_mole * pow(i + 0.5, 2));
    return h * c * (omega_mole * (i + 0.5));
}

double s_func(int i, int f) {

    return static_cast<double>(abs(i - f));
}

double Ns(int i, int f) {
    int max = std::max(i, f);
    int min = std::min(i, f);
    double s = s_func(i, f);

    return pow(tgamma(max + 1) / tgamma(min + 1), 1. / s);

}

double ave_vib_qua(double e1, double e2, double s) {

    double deltaE_if = fabs(e1 - e2);

    return deltaE_if / (s * h_bar);


}


std::pair<double, double> Epsilon_Generator(RanPark &rd) {

    double eta = ETA_HS;
    double xm = 1 / (2 - eta);
    double Fm = xm * pow(1 - xm, 1 - eta);
    double k1;
    double cond;
    bool flag = true;
    while (flag) {
        k1 = rd.uniform();
        cond = (k1 * pow(1 - k1, 1 - eta) / Fm);
        if (rd.uniform() < cond) {
            flag = false;
        }
    }
    double epsilon1 = k1 * rd.uniform();
    double epsilon2 = k1 - epsilon1;
    return std::make_pair(epsilon1, epsilon2);
}

double u_from_E(double E) {
    return sqrt(2 * E / m_r_mole);
}

double theta_p(double w) {
    return (4 *  pow(Pi * w, 2) * m_r_mole)
           / (pow(Morse_mole, 2) * BOLTZ);
}

double gam(double theta1, double theta2, double phi1, double phi2, std::pair<double, double> e12, double y) {


    double tmp = -0.5 * sin(2 * theta1) * cos(phi1)
                 * sqrt(e12.first) - 0.5 * sin(2 * theta2) * cos(phi2)
                                     * sqrt(e12.second) + sqrt((1 - e12.first - e12.second) * (1 - y));

    return std::max(0.0, tmp);
}

double Q_func(double theta_prime, double theta1, double phi1, double theta,double w, double u, double gamma) {
    double Q;
    Q = (theta_prime * Xi_mole * pow(cos(theta1), 2) * pow(cos(phi1), 2))
        / (4 * theta * pow(sinh(Pi * w / (Morse_mole * u * gamma)), 2));
    return Q;
}

//*------------------------------------------------





double Compute_pro_if_deex(int n, int i, int f, const std::vector<double> &Ebins, RanPark &rd, int MS,const double THETA_mole) {
    double E_t_r = Ebins.at(n);
    double p;
    double u1, y1, s1, ns1, w1, gamma1, theta_prime1,theta, q1, p1;
    p = 0;
    theta = THETA_mole;
    u1 = y1 = s1 = ns1 = w1 = gamma1 = theta_prime1 = q1 = p1 = 0;
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
#pragma omp parallel default(none) private (u1, y1, s1, ns1, w1, gamma1, theta_prime1, q1, p1) shared(E_t_r, p, rd, MS, i, f,theta)
#pragma omp  for reduction (+:p)
    for (int j = 0; j < MS; ++j) {

        auto e12 = Epsilon_Generator(rd);
        double u1 = u_from_E(E_t_r);

        double e1 = Vibr_eng(i);
        double e2 = Vibr_eng(f);
        double y1 = rd.uniform();

        double theta1, theta2, phi1, phi2;
        theta1 = rd.uniform() * 2 * Pi;
        theta2 = rd.uniform() * 2 * Pi;
        phi1 = rd.uniform() * 2 * Pi;
        phi2 = rd.uniform() * 2 * Pi;

        double s1 = s_func(i, f);

        double ns1 = Ns(i, f);

        double w1 = ave_vib_qua(e1, e2, s1);

        double gamma1 = gam(theta1, theta2, phi1, phi2, e12, y1);

        double theta_prime1 = theta_p(w1);

        double q1 = Q_func(theta_prime1, theta1, phi1, theta,w1, u1, gamma1);

        double p1 = (pow(ns1, s1) / pow(tgamma(s1 + 1), 2)) * pow(q1, s1)
                    * exp((-2 * ns1 * q1) / (s1 + 1) -
                          (pow(ns1, 2) * pow(q1, 2)) / (pow(s1 + 1, 2) * (s1 + 2)));
        p += p1;

    }
    p /= MS;
    return p;

}

double Compute_pro_if_exci(int n, int i, int f, const std::vector<double> &Ebins, RanPark &rd, int MS,const double THETA_mole) {
    double E_t_r = Ebins.at(n);
    double p;
    double u1, y1, s1, ns1, w1, gamma1, theta_prime1,theta, q1, p1, E_t_f, E_prime, g_i, g_f, em1, em2, factor;
    std::pair<double, double> new_e12;
    int count = 0;
    p = 0;
    theta = THETA_mole;
    u1 = y1 = s1 = ns1 = w1 = gamma1 = theta_prime1 = q1 = p1 = E_t_f = E_prime = g_i = g_f = em1 = em2 = factor = 0;
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
#pragma omp parallel default(none) private (u1, y1, s1, ns1, w1, gamma1, theta_prime1, q1, p1, E_t_f, E_prime, g_i,g_f, em1, em2, factor, new_e12) shared(E_t_r, p, rd, MS, i, f,theta, count)
#pragma omp  for reduction (+:p, count)

    for (int j = 0; j < MS; ++j) {

        auto e12 = Epsilon_Generator(rd);


        double e1 = Vibr_eng(i);
        double e2 = Vibr_eng(f);
        double y1 = rd.uniform();

        double theta1, theta2, phi1, phi2;
        theta1 = rd.uniform() * 2 * Pi;
        theta2 = rd.uniform() * 2 * Pi;
        phi1 = rd.uniform() * 2 * Pi;
        phi2 = rd.uniform() * 2 * Pi;

        //from deex to excitation process

        double E_t_f = E_t_r * (1 - e12.first - e12.second) + e1 - e2;
        if (E_t_f > 0) {

            ++count;


            double E_prime = E_t_f + (e12.first + e12.second) * E_t_r;
            double g_i = sqrt(2 * (1 - e12.first - e12.second) * E_t_r / m_r_N2);
            double g_f = sqrt(2 * E_t_f / m_r_N2);
            double em1 = e12.first * (E_t_r / E_prime);
            double em2 = e12.second * (E_t_r / E_prime);
            new_e12 = std::make_pair(em1, em2);
            double factor = pow(g_f / g_i, 2 - 2 * ETA);
            double u1 = u_from_E(E_prime);
            double s1 = s_func(i, f);

            double ns1 = Ns(i, f);

            double w1 = ave_vib_qua(e1, e2, s1);

            double gamma1 = gam(theta1, theta2, phi1, phi2, new_e12, y1);

            double theta_prime1 = theta_p(w1);

            double q1 = Q_func(theta_prime1, theta1, phi1, theta, w1,u1, gamma1);

            double p1 = (pow(ns1, s1) / pow(tgamma(s1 + 1), 2)) * pow(q1, s1)
                        * exp((-2 * ns1 * q1) / (s1 + 1) -
                              (pow(ns1, 2) * pow(q1, 2)) / (pow(s1 + 1, 2) * (s1 + 2)));
            p1 *= factor;
            p += p1;
        } else {

            p += 0;
        }


    }
    if(count!=0) {
        p /= count;
        std::cout<<count<<' ';
    }else{
        p = 0;
    }
    return p;
}
