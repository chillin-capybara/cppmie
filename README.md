[![CMake](https://github.com/codesaurus97/cppmie/actions/workflows/cmake.yml/badge.svg)](https://github.com/codesaurus97/cppmie/actions/workflows/cmake.yml)

# miecpp

A **header-only** library in **C++20** for calculating Mie scattering efficiencies using the algorithm presented by Hong Du [[1]](#1).

*I am open for suggestions and will accept pull requests that contribute to this library in a useful way. Feel free to create issues and come forward with feature requests.*

## 💽 Installation

### 🐱 As GIT Submodule

Clone the repository into your project folder or add it as a submodule.

```bash
git submodule add https://github.com/codesaurus97/cppmie.git
cd cppmie
git checkout main
```

If you are using CMake for your project, you can simply add the library as a subdirectory and link the header-only library to your project.

```cmake
add_subdirectory(cppmie)
target_link_libraries(myapp cppmie)
```

### 💪 As a global library for CMake

The library can also be installed via CMake as  a global package. To do so, run the following code.

```bash
git clone https://github.com/codesaurus97/cppmie.git
cd cppmie && mkdir build && cd build
cmake ..
cmake --build . --config Release --target install -- -j $(nproc)
```

After the installation is finished the library can be accessed from `CMakeLists.txt` with the `find_package(` directive. To add **cppmie** to your project, you should link against the library `cppmie::cppmie`.

```cmake
find_package(cppmie CONFIG REQUIRED)
target_link_libraries(myapp cppmie::cppmie)
```

## 👨‍💻 Usage

To access the library's functions include the header-file `cppmie.h` from the library as follows.

```c++
#include "cppmie/cppmie.h"
```

After the include the namespace `cppmie::` which containe all the relevant function of the library will be available.

> The namespace `cppmie::helpers` only contains templated helper functions that are used for internal calculations only.

### 😎 Mie scattering efficiencies as data structure

To retrieve the calculated scattering efficiencies as a data structure `cppmie::MieResult` call the function `cppmie::mie()` with two mandatory paramters (intercept parameter + refractive index) and one optional parameter N_STAR.

```c++
const double pi = 3.14159265359;
double diameter = 1.0;             // Diameter in microns
double wavelength = 0.6328;        // Wavelength in microns
double x = diameter * pi / wavelength;
std::complex<double> m{1.5, 0.0};  // Complex refractive index

MieResult result = cppmie::mie(x, m);         // Without N_star
MieResult result = cppmie::mie(x, m, 60000);  // With N_star

std::cout << "Qext  = " << result.Qext  << std::endl;  // Print the results to the console
std::cout << "Qsca  = " << result.Qsca  << std::endl;
std::cout << "Qback = " << result.Qback << std::endl;
```

### ⚙️ Mie scattering efficiencies written back to reference (advanced)

If you are designing your software for maximum efficiency, then return an data structure each time will not give the best performance. The function `cppmie::mie()` offers an alternative call, such that all the values will be written back to a given reference.

```c++
const double pi = 3.14159265359;
double diameter = 1.0;             // Diameter in microns
double wavelength = 0.6328;        // Wavelength in microns
double x = diameter * pi / wavelength;
std::complex<double> m{1.5, 0.0};  // Complex refractive index

double Qext, Qsca, Qback;

cppmie::mie(x, m, Qext, Qsca, Qback);         // Without N_star
cppmie::mie(x, m, Qext, Qsca, Qback, 60000);  // With N_star
```

### 📊 Benchmark

The benchmarks were conducted on a MacBook Pro (13-inch, M1, 2020) with **Apple M1** CPU and 16GB RAM using [google benchmark](https://github.com/google/benchmark).

- Each function uses the write-back method with references.
- The benchmarks with the suffix _opt use the micro-optimized functions.
- All benchmarks are single-threaded.

```
----------------------------------------------------------------------------------------
Benchmark                              Time             CPU   Iterations UserCounters...
----------------------------------------------------------------------------------------
bm_mie_real_refractive            254523 ns       254447 ns         2698 exec_rate=3.9k/s
bm_mie_complex_refractive         255559 ns       255251 ns         2728 exec_rate=3.9k/s
bm_mie_real_refractive_opt        253977 ns       253854 ns         2749 exec_rate=3.9k/s
bm_mie_complex_refractive_opt     253288 ns       253243 ns         2754 exec_rate=3.9k/s
```

## 🚧 Roadmap

✅ --  Implement Mie scattering calculation for `float` and `double` types.

🔲 -- Implement unit tests with pre-calculated results.

🦆 -- Consider optimizing / inlining the rate function r_n(x) for better performance.

## 📔 References

<a id="1">[1]</a> H. Du, "Mie-scattering calculation," Appl. Opt.  43, 1951-1956 (2004).

