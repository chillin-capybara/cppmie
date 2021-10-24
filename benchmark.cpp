//
// Created by Marcell Pigniczki on 25.10.21.
//

#include <benchmark/benchmark.h>
#include "cppmie/cppmie.h"

static void bm_mie_real_refractive(benchmark::State& state) {
	double numIterations = 0;

	// Perform setup here
	for (auto _ : state) {
		// This code gets timed
		cppmie::mie(4.96, 1.5);
		numIterations += 1.0;
	}

	state.counters["exec_rate"] = benchmark::Counter(numIterations, benchmark::Counter::kIsRate);
}

static void bm_mie_complex_refractive(benchmark::State& state) {
	double numIterations = 0;

	// Perform setup here
	for (auto _ : state) {
		// This code gets timed
		cppmie::mie(4.96, {1.5, 0.1});
		numIterations += 1.0;
	}

	state.counters["exec_rate"] = benchmark::Counter(numIterations, benchmark::Counter::kIsRate);
}

static void bm_mie_complex_refractive_opt(benchmark::State& state) {
	double numIterations = 0;

	// Perform setup here
	for (auto _ : state) {
		// This code gets timed
		cppmie::mie(4.96, {1.5, 0.0});
		numIterations += 1.0;
	}

	state.counters["exec_rate"] = benchmark::Counter(numIterations, benchmark::Counter::kIsRate);
}

// Register the function as a benchmark
BENCHMARK(bm_mie_real_refractive);
BENCHMARK(bm_mie_complex_refractive);
BENCHMARK(bm_mie_complex_refractive_opt);

// Run the benchmark
BENCHMARK_MAIN();