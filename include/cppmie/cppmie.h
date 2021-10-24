//
// Created by Marcell Pigniczki on 24.10.21.
//

#ifndef CPPMIE_CPPMIE_H
#define CPPMIE_CPPMIE_H

#include <complex>
#include <vector>

void mie_scattering(double x, std::complex<double> m, int n_star = 60000)
{
	//
	std::vector<std::complex<double>> r(n_star + 1);
	std::complex<double> mx = m * x;

	// Initialize the downward recurrence
	/**
	 * Calculate the complex ratio \[ r_n(x) = \Psi_{n}(x) / \Psi_{n-1}(x) \] and the its downward recurrence
	 * \f[
	 	r_n(mx) = \frac{2n+1}{mx} - \frac{1}{r_{n+1}(mx)}.
	 * \f]
	 */
	r[n_star] = static_cast<double>(2 *  n_star + 1) / mx;
	for (int i = n_star - 1; i >= 1; i--) {
		r[i - 1] = static_cast<double>(2 * i + 1) / mx - 1.0 / (r[i]);
	}

	auto n = static_cast<size_t>(x+4.0*pow(x,1.0/3.0)+2.0+10.0);
	std::vector<double> psi(n + 1);
	std::vector<double> chi(n + 1);
	std::vector<std::complex<double>> a(n);
	std::vector<std::complex<double>> b(n);

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
		psi[i+1] = static_cast<double>(2 * i + 1) * psi[i] / x - psi[i-1];
		chi[i+1] = static_cast<double>(2 * i + 1) * chi[i] / x - chi[i-1];
	}

	for (size_t i = 1; i < n + 1; i++) {
		std::complex<double> factor = r[i - 1] / m + static_cast<double>(i) * (1.0 - 1.0 / (m * m)) / x;
		std::complex<double> zeta_im1{psi[i-1], chi[i-1]};  // zeta[i-1]
		std::complex<double> zeta_i{psi[i], chi[i]};         // zeta[i]

		a[i-1] = (factor * psi[i] - psi[i - 1]) / (factor * zeta_i - zeta_im1);
		b[i-1] = (r[i-1] * m * psi[i] - psi[i-1]) / (r[i-1] * m * zeta_i - zeta_im1);

		if (a[i-1] != a[i-1] || b[i-1] != b[i-1]) {
			break;
		}

	}

	double qsca = 0.0;
	double qext = 0.0;
	double qback = 0.0;
	std::complex<double> qback_pre{0.0, 0.0};

	for (size_t i = 1; i < n + 1; i++) {
		qext += static_cast<double>(2 * i + 1) * (a[i-1].real() + b[i-1].real());
		qsca += static_cast<double>(2 * i + 1) * (std::abs(a[i-1]) * std::abs(a[i-1]) + std::abs(b[i-1]) * std::abs(b[i-1]));
		qback_pre += static_cast<double>(2 * i + 1) * (a[i-1] - b[i-1]) * std::pow(-1.0, i-1);
	}

	qext *= 2.0 / (x * x);
	qsca *= 2.0 / (x * x);
	qback = std::abs(qback_pre) * std::abs(qback_pre) / (x * x);
}

#endif //CPPMIE_CPPMIE_H
