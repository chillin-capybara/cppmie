/**
 * @file cppmie.h
 * @author codesaurus97
 * @brief Implementation of Mie scattering calculation using Hong Du's algorithm.
 * @version 1.0.1
 * @date 2021-10-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef CPPMIE_CPPMIE_H
#define CPPMIE_CPPMIE_H

#include <complex>
#include <vector>
#include <iostream>

/** Default value of N_STAR for the downward recurrence of the rate function
 * r_n(x). */
#define CPPMIE_NSTAR_DEFAULT (60000U)

/**
 * @brief The namespace `cppmie::helpers` contains the helper functions and templates for the
 * calculation of Mie scattering efficiencies.
 */
namespace cppmie::helpers {

/* --- Concepts for template argument requirements ------------------------------------------------------------------ */

template<typename T>
concept FloatingType = std::is_floating_point<T>::value;

template<class T>
struct is_complex_t : public std::false_type {
};

template<class T> requires FloatingType<T>
struct is_complex_t<std::complex<T>> : public std::true_type {
};

}

namespace cppmie {

using std::complex, std::vector;

template<typename T>
    requires helpers::FloatingType<T>
void MieScattering(const T& x, const complex<T>& m,
                   T& qext, T& qsca, T& qback, int n_star = CPPMIE_NSTAR_DEFAULT)
{
    const complex<T> mx = m * x;
    const int n = static_cast<int>(std::ceil(x + 4.0 * pow(x, 1.0 / 3.0) + 2.0 + 10.0));

    vector<complex<T>> rate(n_star + 1);
    rate[n_star] = static_cast<T>(2 * n_star + 1) / mx;

    for (int i = n_star - 1; i >= 0; i--) {
        rate[i] = static_cast<T>(2 * i + 1) / mx - 1. / rate[i + 1];
    }

    // Approximate the Ricatti-Bessel functions
    T psi_0, psi_1, psi_2;
    T chi_0, chi_1, chi_2;
    complex<T> zeta_1, zeta_2;
    complex<T> an, bn;

    psi_0 = std::cos(x);
    psi_1 = std::sin(x);

    chi_0 = std::sin(x);
    chi_1 = - std::cos(x);

    zeta_1 = {psi_1, chi_1};

    qsca = 0.;
    qext = 0.;
    complex<T> qback_stack = {0.0, 0.0};

    for (int i = 1; i <= n; i++) {
        chi_2 = static_cast<T>(2 * i - 1) * chi_1 / x - chi_0;
        psi_2 = static_cast<T>(2 * i - 1) * psi_1 / x - psi_0;
        zeta_2 = {psi_2, chi_2};

        complex<T> rf = (rate[i] / m + static_cast<T>(i) * (1.0 - 1.0 / (m * m)) / x);

        an = (rf * psi_2 - psi_1) / (rf * zeta_2 - zeta_1);
        bn = (rate[i] * m * psi_2 - psi_1)
            / (rate[i] * m * zeta_2 - zeta_1);

        qext += static_cast<T>(2 * i + 1) * (an.real() + bn.real());
        qsca += static_cast<T>(2 * i + 1) * (abs(an) * abs(an) + abs(bn) * abs(bn));
        qback_stack = static_cast<T>(2 * i + 1) * (an - bn) * std::pow(-1.0, i);

        psi_0 = psi_1;
        psi_1 = psi_2;
        chi_0 = chi_1;
        chi_1 = chi_2;
        zeta_1 = zeta_2;
    }

    qext *= 2.0 / (x * x);
    qsca *= 2.0 / (x * x);
    qback = abs(qback_stack) * abs(qback_stack) / (x * x);
}

template<typename T>
    requires helpers::FloatingType<T>
void MieScatteringBaseLine(const T& x, const complex<T>& m,
                   T& qext, T& qsca, T& qback, int n_star = CPPMIE_NSTAR_DEFAULT)
{
    const complex<T> mx = m * x;

    const int n = static_cast<int>(std::ceil(x + 4.0 * pow(x, 1.0 / 3.0) + 2.0 + 10.0));

    vector<complex<T>> rate(n_star + 1);
    rate[n_star] = static_cast<T>(2 * n_star + 1) / mx;

    for (int i = n_star - 1; i >= 0; i--) {
        rate[i] = static_cast<T>(2 * i + 1) / mx - 1. / rate[i + 1];
    }

    // Approximate the Ricatti-Bessel functions
    vector<T> psi(n_star + 2);
    vector<T> chi(n_star + 2);
    vector<complex<T>> zeta(n_star + 2);

    psi[0] = std::cos(x);
    psi[1] = std::sin(x);

    chi[0] = std::sin(x);
    chi[1] = - std::cos(x);

    zeta[0] = {psi[0], chi[0]};
    zeta[1] = {psi[1], chi[1]};

    for (int i = 1; i <= n+1; i++) {
        chi[i + 1] = static_cast<T>(2 * i - 1) * chi[i] / x - chi[i - 1];
        psi[i + 1] = static_cast<T>(2 * i - 1) * psi[i] / x - psi[i - 1];
        zeta[i + 1] = {psi[i + 1], chi[i + 1]};
    }

    vector<complex<T>> a(n + 2);
    vector<complex<T>> b(n + 2);
    for (int i = 0; i <= n+1; i++) {
        complex<T> rf = (rate[i] / m + static_cast<T>(i) * (1.0 - 1.0 / (m * m)) / x);

        a[i] = (rf * psi[i+1] - psi[i]) / (rf * zeta[i+1] - zeta[i]);
        b[i] = (rate[i] * m * psi[i+1] - psi[i])
            / (rate[i] * m * zeta[i+1] - zeta[i]);
    }

    qsca = 0.;
    qext = 0.;
    qback = 0.;
    complex<T> qback_stack {0.0, 0.0};

    for (int i = 1; i <= n; i++) {
        qext += static_cast<T>(2 * i + 1) * (a[i].real() + b[i].real());
        qsca += static_cast<T>(2 * i + 1) * (pow(abs(a[i]), 2.) + pow(abs(b[i]), 2.));
        qback_stack += static_cast<T>(2 * i + 1) * (a[i] - b[i]) * std::pow(-1.0, i);
    }

    qext *= 2.0 / (x * x);
    qsca *= 2.0 / (x * x);
    qback = abs(qback_stack) * abs(qback_stack) / (x * x);
}


}

#endif //CPPMIE_CPPMIE_H
