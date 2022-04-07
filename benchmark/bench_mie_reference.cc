#include <cppmie/cppmie.h>
#include <complex>
#include <vector>
#include <tuple>

#include <benchmark/benchmark.h>

using std::vector, std::complex, std::tuple;
const static vector<tuple<complex<double>, double, double, double>> parameters{
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
};

static void BM_MieReference(benchmark::State& state) {
    const int id = state.range(0);
    const auto& [m,x,gold_qext,gold_qsca] = parameters[id];
    double qext, qback, qsca;
    for (auto _ : state) {
        cppmie::MieScattering(x, m, qext, qsca, qback);
        benchmark::DoNotOptimize(qext);
        benchmark::DoNotOptimize(qback);
        benchmark::DoNotOptimize(qsca);
        benchmark::ClobberMemory();
    }
}

static void BM_MieReferenceBaseline(benchmark::State& state) {
    const int id = state.range(0);
    const auto& [m,x,gold_qext,gold_qsca] = parameters[id];
    double qext, qback, qsca;
    for (auto _ : state) {
        cppmie::MieScattering(x, m, qext, qsca, qback);
        benchmark::DoNotOptimize(qext);
        benchmark::DoNotOptimize(qback);
        benchmark::DoNotOptimize(qsca);
        benchmark::ClobberMemory();
    }
}

BENCHMARK(BM_MieReference)->DenseRange(0,parameters.size(), 1);
BENCHMARK(BM_MieReferenceBaseline)->DenseRange(0,parameters.size(), 1);
BENCHMARK_MAIN();
