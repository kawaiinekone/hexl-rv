#include "hexl_rv/ntt.hpp"
#include "hexl_rv/number_theory.hpp"
#include <iostream>
#include <vector>
#include <random>

int main() {
    std::cout << "[TEST] Running NTT & INTT roundtrips..." << std::endl;
    const uint64_t degree = 8;
    const uint64_t modulus = 1153;
    uint64_t root = hexl_rv::PowerMod(5, (modulus - 1) / (2 * degree), modulus);

    hexl_rv::NTT ntt(degree, modulus, root);
    std::vector<uint64_t> poly = {1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<uint64_t> orig = poly;

    ntt.ComputeForward(poly.data());
    ntt.ComputeInverse(poly.data());

    for (size_t i = 0; i < degree; ++i) {
        if (poly[i] != orig[i]) return 1;
    }

    const uint64_t large_degree = 1024;
    const uint64_t large_modulus = 1125899906842625ULL;
    uint64_t large_root = hexl_rv::PowerMod(3, (large_modulus - 1) / (2 * large_degree), large_modulus);

    hexl_rv::NTT large_ntt(large_degree, large_modulus, large_root);
    std::vector<uint64_t> large_orig(large_degree);
    std::mt19937_64 rng(777);
    for (size_t i = 0; i < large_degree; ++i) {
        large_orig[i] = rng() % large_modulus;
    }

    std::vector<uint64_t> large_poly = large_orig;
    large_ntt.ComputeForward(large_poly.data());
    large_ntt.ComputeInverse(large_poly.data());

    for (size_t i = 0; i < large_degree; ++i) {
        if (large_poly[i] != large_orig[i]) return 1;
    }

    std::cout << "[PASS] NTT / INTT roundtrips passed." << std::endl;
    return 0;
}
