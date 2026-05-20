# Building the Project

> **Quick start:** A pre-compiled `langevin_mex.mexw64` for Windows x64 is included in the project root. If you are on Windows x64 and just want to run the examples, skip directly to [Building the MATLAB MEX File](#building-the-matlab-mex-file) — you only need to run `build_mex.m` inside MATLAB.

## Requirements

| Tool | Minimum Version |
|------|----------------|
| CMake | 3.18 |
| C++ compiler | MSVC 2019 / GCC 9 / Clang 10 (must support C++17) |
| Internet access | Required at configure time (FetchContent downloads Eigen3 and NLopt) |

Eigen3 and NLopt are downloaded automatically — no manual installation needed.

---

## Windows (MSVC)

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

The output is `build/Release/langevin_core.lib` (static library).

**AVX2 note:** The build enables `/arch:AVX2`.  If your CPU does not support AVX2, remove that flag from `CMakeLists.txt` before building.

---

## Linux / macOS (GCC or Clang)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Output: `build/liblangevin_core.a`.

---

## CMake Targets

| Target | Type | Description |
|--------|------|-------------|
| `gsl_mini` | Static library | Custom GSL implementation (splines + quasi-random) |
| `langevin_core` | Static library | Main library; links `gsl_mini` and NLopt |

---

## Compiler Flags

| Compiler | Flags |
|----------|-------|
| MSVC | `/O2 /W3 /MP /EHsc /arch:AVX2 /DEIGEN_NO_DEBUG` |
| GCC/Clang | `-O3 -march=native -Wall -Wextra -DEIGEN_NO_DEBUG` |

`EIGEN_NO_DEBUG` disables Eigen bounds-checking assertions for maximum performance.

---

## Building the MATLAB MEX File

The MEX file is not part of the CMake build.  Compile it from within MATLAB after building the core library:

```matlab
% Example — adjust paths for your system
mex -O ...
    CXXFLAGS='$CXXFLAGS -std=c++17' ...
    -I/path/to/eigen3/include ...
    langevin_mex.cpp langevin_core.cpp gsl_impl.cpp ...
    build/Release/langevin_core.lib ...     % Windows
    -lnlopt ...
    -o langevin_mex
```

See the project `README.md` for a ready-to-run build script.

---

## Verifying the Build

A quick smoke-test from C++:

```cpp
#include "langevin_core.h"
#include <cstdio>
#include <cmath>

int main() {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = 0.1 * i;

    LC_Result *r = euler_reconstruct(
        data, 100, 0.1,
        NULL, 3, NULL, 2,
        NULL, NULL, NAN, NAN,
        LC_CUBIC, LC_ADDITIVE,
        LC_SLSQP, 1, 1000, 1e-6
    );

    if (r) {
        printf("OK — NLL = %f\n", r->negloglik);
        lc_result_free(r);
    } else {
        printf("FAIL\n");
        return 1;
    }
    return 0;
}
```

Link against `langevin_core` (and its transitive dependencies `gsl_mini`, `nlopt`).
