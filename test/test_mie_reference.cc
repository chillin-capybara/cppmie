#include <cppmie/cppmie.h>
#include <gtest/gtest.h>

#include <complex>
using std::complex;

class MieScatteringTestSet : public ::testing::TestWithParam<std::tuple<complex<double>, double, double, double>> {
 public:
  double dummy;
};

constexpr double ApplyPrecision(const double& value, int precision)
{
    const auto exponent = static_cast<int>(std::floor(std::log10(value)))+1;
    const double fixed_value = value * std::pow(10.,  -exponent + precision);
    const auto round_value = static_cast<int64_t>(std::round(fixed_value));
    return static_cast<double>(round_value) * std::pow(10., - (-exponent+precision));

}

TEST_P(MieScatteringTestSet, QbackQext)
{
    const auto m = std::get<0>(GetParam());
    const auto x = std::get<1>(GetParam());
    const auto qext_gold = std::get<2>(GetParam());
    const auto qsca_gold = std::get<3>(GetParam());

    double qext, qsca, qback;
    cppmie::MieScattering(x, m, qext, qsca, qback);

    // test for 6 value precisions (NOT DECIMALS)
    ASSERT_FLOAT_EQ(ApplyPrecision(qext,6), qext_gold);
    ASSERT_FLOAT_EQ(ApplyPrecision(qsca, 6), qsca_gold);
}

INSTANTIATE_TEST_SUITE_P(
    MieScatteringValidation,
    MieScatteringTestSet,
    ::testing::Values(
        std::make_tuple(complex<double>{0.75, 0.}, 0.099, 7.41786e-6, 7.41786e-6),
        std::make_tuple(complex<double>{0.75, 0.}, 0.101, 8.03354e-6, 8.03354e-6),
        std::make_tuple(complex<double>{0.75, 0.}, 10., 2.23226, 2.23226),
        std::make_tuple(complex<double>{0.75, 0.}, 1000., 1.99791, 1.99791),
        std::make_tuple(complex<double>{1.33, 1.e-5}, 100., 2.10132, 2.09659),
        std::make_tuple(complex<double>{1.33, 1.e-5}, 10000., 2.00409, 1.72386),
        std::make_tuple(complex<double>{1.5, 1.}, 0.055, 0.101491, 1.13169e-5),
        std::make_tuple(complex<double>{1.5, 1.}, 0.056, 0.103347, 1.21631e-5),
        std::make_tuple(complex<double>{1.5, 1.}, 100., 2.09750, 1.28370),
        std::make_tuple(complex<double>{1.5, 1.}, 10000., 2.00437, 1.23657),
        std::make_tuple(complex<double>{10., 10.}, 1., 2.53299, 2.04941),
        std::make_tuple(complex<double>{10., 10.}, 100., 2.07112, 1.83679),
        std::make_tuple(complex<double>{10., 10.}, 10000., 2.00591, 1.79539)
    ));
