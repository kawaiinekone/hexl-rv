#include "hexl_rv/ntt.hpp"
#include "hexl_rv/eltwise_mult_mod.hpp"
#include "hexl_rv/number_theory.hpp"

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <cassert>

// O(N^2) Schoolbook Negacyclic Polynomial Multiplication in Z_q[X]/(X^N + 1)
void PolyMultSchoolbookNegacyclic(uint64_t* result,
                                  const uint64_t* op1,
                                  const uint64_t* op2,
                                  uint64_t n,
                                  uint64_t modulus) {
    std::vector<unsigned __int128> acc(2 * n, 0);

    // Standard acyclic polynomial multiplication
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            acc[i + j] += static_cast<unsigned __int128>(op1[i]) * op2[j];
        }
    }

    // Negacyclic folding: X^N = -1 (mod X^N + 1)
    for (uint64_t i = 0; i < n; ++i) {
        unsigned __int128 low = acc[i];
        unsigned __int128 high = acc[i + n];
        uint64_t low_mod = static_cast<uint64_t>(low % modulus);
        uint64_t high_mod = static_cast<uint64_t>(high % modulus);

        result[i] = (low_mod >= high_mod) 
            ? (low_mod - high_mod) 
            : (low_mod + modulus - high_mod);
    }
}

// O(N log N) Fast Convolution using NTT
void PolyMultNTT(uint64_t* result,
                 const uint64_t* op1,
                 const uint64_t* op2,
                 const hexl_rv::NTT& ntt) {
    uint64_t n = ntt.GetDegree();
    uint64_t modulus = ntt.GetModulus();

    std::vector<uint64_t> buf1(op1, op1 + n);
    std::vector<uint64_t> buf2(op2, op2 + n);

    // 1. Forward NTT on both operands
    ntt.ComputeForward(buf1.data());
    ntt.ComputeForward(buf2.data());

    // 2. Element-wise Dyadic Modular Multiplication (uses RV64 assembly kernel)
    hexl_rv::EltwiseMultMod(result, buf1.data(), buf2.data(), n, modulus, 1);

    // 3. Inverse NTT to transform back to polynomial representation
    ntt.ComputeInverse(result);
}

void BenchmarkDegree(uint64_t degree, uint64_t modulus, uint64_t root, size_t iterations) {
    std::cout << "\n==================================================================" << std::endl;
    std::cout << " Evaluating Degree N = " << degree 
              << ", Modulus q = " << modulus << " (" << hexl_rv::MSB(modulus) + 1 << "-bit)" << std::endl;
    std::cout << "==================================================================" << std::endl;

    hexl_rv::NTT ntt(degree, modulus, root);
    std::mt19937_64 rng(42);

    std::vector<uint64_t> op1(degree);
    std::vector<uint64_t> op2(degree);
    std::vector<uint64_t> res_schoolbook(degree, 0);
    std::vector<uint64_t> res_ntt(degree, 0);

    for (size_t i = 0; i < degree; ++i) {
        op1[i] = rng() % modulus;
        op2[i] = rng() % modulus;
    }

    // Correctness Verification: Schoolbook vs NTT
    PolyMultSchoolbookNegacyclic(res_schoolbook.data(), op1.data(), op2.data(), degree, modulus);
    PolyMultNTT(res_ntt.data(), op1.data(), op2.data(), ntt);

    for (size_t i = 0; i < degree; ++i) {
        if (res_schoolbook[i] != res_ntt[i]) {
            std::cerr << "[FAIL] Output mismatch at index " << i << ": Schoolbook = "
                      << res_schoolbook[i] << ", NTT = " << res_ntt[i] << std::endl;
            exit(1);
        }
    }
    std::cout << "[PASS] Mathematical equivalence verified: INTT(NTT(A) * NTT(B)) == A * B (mod X^N + 1, q)" << std::endl;

    // Time Schoolbook
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t it = 0; it < iterations; ++it) {
        PolyMultSchoolbookNegacyclic(res_schoolbook.data(), op1.data(), op2.data(), degree, modulus);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double time_schoolbook_ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / iterations;

    // Time NTT
    auto t2 = std::chrono::high_resolution_clock::now();
    for (size_t it = 0; it < iterations; ++it) {
        PolyMultNTT(res_ntt.data(), op1.data(), op2.data(), ntt);
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    double time_ntt_ms = std::chrono::duration<double, std::milli>(t3 - t2).count() / iterations;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << " • Schoolbook O(N^2) Latency : " << time_schoolbook_ms << " ms / mult" << std::endl;
    std::cout << " • Accelerated NTT Latency   : " << time_ntt_ms << " ms / mult" << std::endl;
    std::cout << " • Speedup Factor           : " << (time_schoolbook_ms / time_ntt_ms) << "x faster" << std::endl;
}

int main() {
    std::cout << ">>> HEXL-RV RING POLYNOMIAL MULTIPLICATION BENCHMARK <<<" << std::endl;

    // 50-bit prime supporting 2N-th roots of unity (modulus = 1 mod 2N)
    const uint64_t modulus = 1125899906842625ULL;

    // N = 1024
    uint64_t root_1024 = hexl_rv::PowerMod(3, (modulus - 1) / (2 * 1024), modulus);
    BenchmarkDegree(1024, modulus, root_1024, 5);

    // N = 2048
    uint64_t root_2048 = hexl_rv::PowerMod(3, (modulus - 1) / (2 * 2048), modulus);
    BenchmarkDegree(2048, modulus, root_2048, 2);

    return 0;
}
