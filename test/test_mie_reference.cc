#include <cppmie/cppmie.h>
#include <gtest/gtest.h>

#include <complex>
using std::complex;

class MieScatteringTestSet : public ::testing::TestWithParam<std::tuple<complex<double>, double, double, double>> {
 public:
  double dummy;
};

class MieScatteringWiscombe : public ::testing::TestWithParam<std::tuple<complex<double>, double, double>> {
 public:
  double dummy;
};

static double ApplyPrecision(const double& value, int precision)
{
    const auto exponent = static_cast<int>(std::floor(std::log10(value)))+1;
    const double fixed_value = value * std::pow(10.,  -exponent + precision);
    const auto round_value = static_cast<int64_t>(std::round(fixed_value));
    return static_cast<double>(round_value) * std::pow(10., - (-exponent+precision));

}

TEST_P(MieScatteringTestSet, QbackQext_opt)
{
    const auto m = std::get<0>(GetParam());
    const auto x = std::get<1>(GetParam());
    const auto qext_gold = std::get<2>(GetParam());
    const auto qsca_gold = std::get<3>(GetParam());

    double qext;
    double qsca;
    double qback;
    cppmie::MieScattering(x, m, qext, qsca, qback);

    // test for 6 value precisions (NOT DECIMALS)
    ASSERT_FLOAT_EQ(ApplyPrecision(qext,6), qext_gold);
    ASSERT_FLOAT_EQ(ApplyPrecision(qsca, 6), qsca_gold);
}

TEST_P(MieScatteringTestSet, QbackQext_baseline)
{
    const auto m = std::get<0>(GetParam());
    const auto x = std::get<1>(GetParam());
    const auto qext_gold = std::get<2>(GetParam());
    const auto qsca_gold = std::get<3>(GetParam());

    double qext;
    double qsca;
    double qback;
    cppmie::MieScatteringBaseLine(x, m, qext, qsca, qback);

    // test for 6 value precisions (NOT DECIMALS)
    ASSERT_DOUBLE_EQ(ApplyPrecision(qext,6), qext_gold);
    ASSERT_DOUBLE_EQ(ApplyPrecision(qsca, 6), qsca_gold);
}

TEST_P(MieScatteringWiscombe, QbackWiscombe)  // Wiscombe test for backscattering
{
    const auto m = std::get<0>(GetParam());
    const auto x = std::get<1>(GetParam());
    const auto qback_wiscombe = std::get<2>(GetParam());

    double qext;
    double qsca;
    double qback;
    cppmie::MieScatteringBaseLine(x, m, qext, qsca, qback);

    ASSERT_DOUBLE_EQ(qback, qback_wiscombe);
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


INSTANTIATE_TEST_SUITE_P(
    MieScatteringBaseline,
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

// Wiscombe Tests for Mie Scattering Qback
// Source for Wiscombe Test Values: https://miepython.readthedocs.io/en/latest/09_backscattering.html
INSTANTIATE_TEST_SUITE_P(
    WiscombePerfectlyConducting,
    MieScatteringWiscombe,
    ::testing::Values(
        std::make_tuple(complex<double>{1.5500, 0.0}, 5.213, 2.925340e+00),
        std::make_tuple(complex<double>{0.0000, -1000.00}, 0.099, 8.630064e-04),
        std::make_tuple(complex<double>(0.0000 -1000.00), 0.101, 9.347732e-04),
        std::make_tuple(complex<double>(0.0000 -1000.00), 100., 9.990256e-01),
        std::make_tuple(complex<double>(0.0000 -1000.00), 10000., 9.999997e-01)));

INSTANTIATE_TEST_SUITE_P(
    SmallerRefractiveThanEnvironment,
    MieScatteringWiscombe,
    ::testing::Values(
        std::make_tuple(complex<double>{0.75, 0.0}, 0.099, 1.108554e-05),
        std::make_tuple(complex<double>{0.75, 0.0}, 0.101, 1.200382e-05),
        std::make_tuple(complex<double>(0.75, 0.0), 10.000, 4.658462e-02),
        std::make_tuple(complex<double>(0.75, 0.0), 1000.000, 9.391600e-01)));

INSTANTIATE_TEST_SUITE_P(
    NonAbsorbing,
    MieScatteringWiscombe,
    ::testing::Values(
        std::make_tuple(complex<double>{1.5, 0.0}, 10., 1.695084e+00),
        std::make_tuple(complex<double>{1.5, 0.0}, 100., 1.736102e+00),
        std::make_tuple(complex<double>(1.5, 0.0), 1000., 1.030182e+01)));

INSTANTIATE_TEST_SUITE_P(
    WaterDroplets,
    MieScatteringWiscombe,
    ::testing::Values(
        std::make_tuple(complex<double>{1.3300, -0.00001}, 10., 8.462494e-02),
        std::make_tuple(complex<double>{1.3300, -0.00001}, 100., 2.146327e+00),
        std::make_tuple(complex<double>(1.3300, -0.00001), 1000., 3.757215e-02)));

INSTANTIATE_TEST_SUITE_P(
    ModeratelyAbsorbingSpheres,
    MieScatteringWiscombe,
    ::testing::Values(
        std::make_tuple(complex<double>{1.5, -1.}, 0.055, 1.695493e-05),
        std::make_tuple(complex<double>{1.5, -1.}, 0.056, 1.822197e-05),
        std::make_tuple(complex<double>{1.5, -1.}, 1., 5.730036e-01),
        std::make_tuple(complex<double>{1.5, -1.}, 100., 1.724214e-01),
        std::make_tuple(complex<double>(1.5, -1.), 1000., 1.724138e-01)));

INSTANTIATE_TEST_SUITE_P(
    BigRefraction,
    MieScatteringWiscombe,
    ::testing::Values(
        std::make_tuple(complex<double>{10., -10.}, 1., 3.308998e+00),
        std::make_tuple(complex<double>{10., -10.}, 100., 8.201267e-01),
        std::make_tuple(complex<double>(10., -10.), 1000., 8.190052e-01)));