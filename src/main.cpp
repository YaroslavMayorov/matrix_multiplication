#include "../include/Algorithms.h"
#include "../include/Matrix.h"
#include "../include/Operations.h"
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
using namespace std::chrono;

template <typename T>
double calculate_memory_kb(const string &algorithm, int n) {
  double base_size = sizeof(T);
  double n_squared = (double)n * n;

  if (algorithm == "Naive" || algorithm == "BLAS") {
    // A + B + C = 3N²
    return (3.0 * n_squared * base_size) / 1024.0;
  } else if (algorithm == "Winograd") {
    // A + B + C + row_factors(N) + col_factors(N) = 3N² + 2N
    return ((3.0 * n_squared + 2.0 * n) * base_size) / 1024.0;
  } else if (algorithm == "Strassen") {
    // Strassen recursively allocates:
    // - 19 submatrices of size (N/2)² each at each level
    // - 3 matrices of size N²
    return (7.75 * n_squared * base_size) / 1024.0;
  }

  // Default fallback
  return (3.0 * n_squared * base_size) / 1024.0;
}

template <typename T>
double calculate_error(const Matrix<T> &Ref, const Matrix<T> &Result) {
  double error = 0.0;
  for (int i = 0; i < Ref.rows * Ref.cols; ++i) {
    if constexpr (std::is_same_v<T, int>) {
      error += std::abs(Ref.data[i] - Result.data[i]);
    } else {
      error += std::abs(Ref.data[i] - Result.data[i]);
    }
  }
  return error;
}

template <typename T>
void run_benchmark_for_type(const string &type_name_prefix,
                            ofstream &csv_file) {
  vector<int> sizes = {4, 8, 16, 32, 64, 128, 256, 512};
  vector<string> types = {"Random", "Symmetric", "Identity"};
  const int TARGET_DURATION_MS = 100;

  for (const string &type : types) {
    for (int n : sizes) {
      Matrix<T> A(n, n), B(n, n);
      if (type == "Random") {
        A = Matrix<T>::random(n, n);
        B = Matrix<T>::random(n, n);
      } else if (type == "Symmetric") {
        A = Matrix<T>::symmetric(n);
        B = Matrix<T>::symmetric(n);
      } else if (type == "Identity") {
        A = Matrix<T>::identity(n);
        B = Matrix<T>::identity(n);
      }

      // Calculate Reference for Error Check (Naive)
      // Save stats because naive_multiply increments them
      long long pre_adds = stats.additions;
      long long pre_mults = stats.multiplications;
      Matrix<T> Reference = naive_multiply(A, B);
      stats.additions = pre_adds;
      stats.multiplications = pre_mults;

      auto benchmark_algo = [&](const string &name, auto func) {
        stats.reset();

        // 1. Run for Result and Error check
        Matrix<T> Result = func(A, B);

        double error = 0.0;
        if (name == "Naive" || name == "BLAS") {
          error = 0.0;
        } else {
          error = calculate_error(Reference, Result);
        }

        // 2. Run for Timing
        auto start = high_resolution_clock::now();
        func(A, B);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        long long single_run_duration = duration.count();
        int reps = 1;

        if (single_run_duration < TARGET_DURATION_MS) {
          if (single_run_duration == 0)
            reps = 100;
          else
            reps = TARGET_DURATION_MS / single_run_duration + 1;

          if (reps > 1000)
            reps = 1000;
          if (reps < 5)
            reps = 5;

          stats.reset();
          func(A, B);
          long long stored_adds = stats.additions;
          long long stored_mults = stats.multiplications;

          start = high_resolution_clock::now();
          for (int i = 0; i < reps; ++i) {
            func(A, B);
          }
          stop = high_resolution_clock::now();
          duration = duration_cast<milliseconds>(stop - start);

          stats.additions = stored_adds;
          stats.multiplications = stored_mults;
        } else {
          stats.reset();
          func(A, B);
        }

        double avg_time = (double)duration.count() / reps;
        double memory_kb = calculate_memory_kb<T>(name, n);

        cout << type_name_prefix << " " << type << "," << n << "," << name
             << "," << stats.additions << "," << stats.multiplications << ","
             << avg_time << "," << reps << "," << memory_kb << "," << error
             << endl;

        csv_file << type_name_prefix << " " << type << "," << n << "," << name
                 << "," << stats.additions << "," << stats.multiplications
                 << "," << avg_time << "," << reps << "," << memory_kb << ","
                 << error << endl;
      };

      benchmark_algo("Naive", naive_multiply<T>);
      benchmark_algo("Strassen", strassen_multiply<T>);
      benchmark_algo("Winograd", winograd_multiply<T>);

#ifdef HAS_ACCELERATE
      benchmark_algo("BLAS", blas_multiply<T>);
#endif
    }
  }
}

void run_benchmark() {
  ofstream csv_file("benchmark_results.csv");
  if (!csv_file.is_open()) {
    cerr << "Failed to open benchmark_results.csv for writing" << endl;
    return;
  }

  cout << "Starting Benchmark..." << endl;
  string header = "DataType Type,Size,Algorithm,Additions,Multiplications,Time("
                  "ms),Reps,Memory(KB),Error";
  cout << header << endl;
  csv_file << header << endl;

  run_benchmark_for_type<int>("Int", csv_file);

  cout << "--- Switching to Double ---" << endl;
  run_benchmark_for_type<double>("Double", csv_file);

  cout << "--- Switching to Complex ---" << endl;
  run_benchmark_for_type<std::complex<double>>("Complex", csv_file);

  csv_file.close();
  cout << "Benchmark results written to benchmark_results.csv" << endl;
}

int main() {
  try {
    run_benchmark();
  } catch (const exception &e) {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }
  return 0;
}
