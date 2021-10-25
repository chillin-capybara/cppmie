//
// Created by Marcell Pigniczki on 24.10.21.
//

#ifndef CPPMIE_CPPMIE_H
#define CPPMIE_CPPMIE_H

#include <complex>
#include <vector>
#include <iostream>

// TODO: make some template requirements

#define CPPMIE_NSTAR_DEFAULT (60000U)

namespace cppmie::helpers {

template<typename T>
static inline std::vector<T> calc_r(const T& mx, size_t n_star)
{
	T mx_inv = T(1) / mx;  // Pre-Calculate 1/mx to speed up the computation

	std::vector<T> r(n_star + 1);
	r[n_star] = static_cast<T>(2 * n_star + 1) * mx_inv;

	/**
	 * Calculate the complex ratio \[ r_n(x) = \Psi_{n}(x) / \Psi_{n-1}(x) \] and the its downward recurrence
	 * from \[ N_* \] until 0.
	 * \f[
		 r_n(mx) = \frac{2n+1}{mx} - \frac{1}{r_{n+1}(mx)}.
	 * \f]
	 */
	for (size_t i = n_star - 1; i >= 1; i--) {
		r[i - 1] = static_cast<T>(2 * i + 1) * mx_inv - 1.0 / (r[i]);
	}

	return std::move(r);
}

template<typename TIntercept, typename TRefractive>
static inline void mie_core(const TIntercept& x, const TRefractive& m, size_t n_star)
{
	TRefractive mx     = m * x;
	TRefractive mx_inv = TIntercept(1) / mx;

	// Determine the number of elements in for a_n and b_n
	auto n = static_cast<size_t>(x + 4.0 * std::pow(x, 1.0 / 3.0) + 2.0 + 10.0);

	// Calculate the rate of the Psi(x) function using recursion
	auto r = calc_r(mx, n_star);

	std::vector<TIntercept> psi(n + 1);
	std::vector<TIntercept> chi(n + 1);

	/**
	 * \f[\Psi_{-1}(x) = sin(x) f\]
	 * \f[\Psi_{0}(x) = \Psi_{-1} / x - cos(x) f\]
	 */
	psi[0] = std::sin(x);
	psi[1] = psi[0] / x - std::cos(x);

	/**
	 * \f[\Chi_{-1}(x) = cos(x) f\]
	 * \f[\Chi_{0}(x) = \Chi_{-1} / x + sin(x) f\]
	 */
	chi[0] = std::cos(x);
	chi[1] = chi[0] / x + std::sin(x);

	for (size_t i = 1; i < n; i++) {
		psi[i + 1] = static_cast<TIntercept>(2 * i + 1) * psi[i] / x - psi[i - 1];
		chi[i + 1] = static_cast<TIntercept>(2 * i + 1) * chi[i] / x - chi[i - 1];
	}

	std::vector<std::complex<TIntercept>> a(n);
	std::vector<std::complex<TIntercept>> b(n);

	for (size_t i = 1; i < n + 1; i++) {  // TODO: Optimize the memory usage by calculating it simultaneously with the previous loop
		std::complex<TIntercept> factor = r[i - 1] / m + static_cast<TIntercept>(i) * (1.0 - 1.0 / (m * m)) / x;
		std::complex<TIntercept> zeta_im1{psi[i - 1], chi[i - 1]};  // zeta[i-1]
		std::complex<TIntercept> zeta_i{psi[i], chi[i]};            // zeta[i]

		a[i - 1] = (factor * psi[i] - psi[i - 1]) / (factor * zeta_i - zeta_im1);
		b[i - 1] = (r[i - 1] * m * psi[i] - psi[i - 1]) / (r[i - 1] * m * zeta_i - zeta_im1);
	}

	TIntercept               qsca  = TIntercept(0);
	TIntercept               qext  = TIntercept(0);
	TIntercept               qback = TIntercept(0);
	std::complex<TIntercept> qback_pre{0.0, 0.0};

	for (size_t i = 1; i < n + 1; i++) {
		qext += static_cast<TIntercept>(2 * i + 1) * (a[i - 1].real() + b[i - 1].real());
		qsca += static_cast<TIntercept>(2 * i + 1)
				* (std::abs(a[i - 1]) * std::abs(a[i - 1]) + std::abs(b[i - 1]) * std::abs(b[i - 1]));
		qback_pre += static_cast<TIntercept>(2 * i + 1) * (a[i - 1] - b[i - 1]) * std::pow(-1.0, i - 1);
	}

	qext *= 2.0 / (x * x);
	qsca *= 2.0 / (x * x);
	qback         = std::abs(qback_pre) * std::abs(qback_pre) / (x * x);
}

/**
 * @brief Microoptimized algorithm to calculate the mie scattering values.
 * @details This implementation of the algorithm aims to minimize the memory usage and maximize the calculation speed
 * by keeping most of the values in CPU registers.
 * @tparam TIntercept
 * @tparam TRefractive
 * @param x
 * @param m
 * @param n_star
 */
template<typename TIntercept, typename TRefractive>
static inline void mie_core_microopt(const TIntercept& x, const TRefractive& m, size_t n_star)
{
	TRefractive mx     = m * x;
	TRefractive mx_inv = TIntercept(1) / mx;

	// Determine the number of elements in for a_n and b_n
	auto n = static_cast<size_t>(x + 4.0 * std::pow(x, 1.0 / 3.0) + 2.0 + 10.0);

	// Calculate the rate of the Psi(x) function using recursion
	auto r = calc_r(mx, n_star);

	/**
	 * \f[\Psi_{-1}(x) = sin(x) f\]
	 * \f[\Psi_{0}(x) = \Psi_{-1} / x - cos(x) f\]
	 */
	TIntercept psi_0 = std::sin(x);
	TIntercept psi_1 = psi_0 / x - std::cos(x);
	TIntercept psi_2;

	/**
	 * \f[\Chi_{-1}(x) = cos(x) f\]
	 * \f[\Chi_{0}(x) = \Chi_{-1} / x + sin(x) f\]
	 */
	TIntercept chi_0 = std::cos(x);
	TIntercept chi_1 = chi_0 / x + std::sin(x);
	TIntercept chi_2;

	std::complex<TIntercept> a;
	std::complex<TIntercept> b;

	TIntercept               qsca  = TIntercept(0);
	TIntercept               qext  = TIntercept(0);
	TIntercept               qback = TIntercept(0);
	std::complex<TIntercept> qback_pre{0.0, 0.0};

	std::complex<TIntercept> factor;
	std::complex<TIntercept> zeta_0;
	std::complex<TIntercept> zeta_1;

	for (size_t i = 1; i < n + 1; i++) {
		if (i < n) {
			psi_2 = static_cast<TIntercept>(2 * i + 1) * psi_1 / x - psi_0;
			chi_2 = static_cast<TIntercept>(2 * i + 1) * chi_1 / x - chi_0;
		}

		factor = r[i - 1] / m + static_cast<TIntercept>(i) * (1.0 - 1.0 / (m * m)) / x;
		zeta_0 = {psi_0, chi_0};  // zeta[i-1]
		zeta_1 = {psi_1, chi_1};  // zeta[i]

		a = (factor * psi_1 - psi_0) / (factor * zeta_1 - zeta_0);
		b = (r[i - 1] * m * psi_1 - psi_0) / (r[i - 1] * m * zeta_1 - zeta_0);

		qext += static_cast<TIntercept>(2 * i + 1) * (a.real() + b.real());
		qsca += static_cast<TIntercept>(2 * i + 1)
				* (std::abs(a) * std::abs(a) + std::abs(b) * std::abs(b));
		qback_pre += static_cast<TIntercept>(2 * i + 1) * (a - b) * std::pow(-1.0, i - 1);

		// Shift the register values into the correct index
		psi_0 = psi_1;
		psi_1 = psi_2;
		chi_0 = chi_1;
		chi_1 = chi_2;
	}

	qext *= 2.0 / (x * x);
	qsca *= 2.0 / (x * x);
	qback         = std::abs(qback_pre) * std::abs(qback_pre) / (x * x);
}

}

namespace cppmie {

template<typename T>
void mie(const T& x, const T& m, size_t n_star = CPPMIE_NSTAR_DEFAULT)
{
	helpers::mie_core_microopt(x, m, n_star);
}

template<typename T>
void mie(const T& x, const std::complex<T>& m, size_t n_star = CPPMIE_NSTAR_DEFAULT)
{
	if (m.imag() == T(0)) {  // If complex type was given, but it's still a real number
		helpers::mie_core_microopt(x, m.real(), n_star);
	} else {
		helpers::mie_core_microopt(x, m, n_star);
	}
}

}


#endif //CPPMIE_CPPMIE_H
