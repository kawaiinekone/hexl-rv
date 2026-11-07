#include "hexl_rv/eltwise_mult_mod.hpp"
#include <iostream>
#include <vector>
#include <random>

int main() {
    std::cout << "[TEST] Running EltwiseMultMod (V-V & V-S)..." << std::endl;
    const uint64_t n = 8;
    const uint64_t q = 1153;
    std::vector<uint64_t> op1 = {0, 1, 10, 500, 1152, 200, 1000, 7};
    std::vector<uint64_t> op2 = {10, 1, 5, 2, 2, 5, 3, 1152};
    std::vector<uint64_t> res(n, 0);

    hexl_rv::EltwiseMultMod(res.data(), op1.data(), op2.data(), n, q, 1);
    for (size_t i = 0; i < n; ++i) {
        uint64_t exp = (static_cast<unsigned __int128>(op1[i]) * op2[i]) % q;
        if (res[i] != exp) return 1;
    }

    const size_t test_size = 1024;
    uint64_t modulus = (1ULL << 45) + 113;
    std::mt19937_64 rng(42);
    std::vector<uint64_t> a(test_size), b(test_size), out(test_size);

    for (size_t i = 0; i < test_size; ++i) {
        a[i] = rng() % modulus;
        b[i] = rng() % modulus;
    }

    hexl_rv::EltwiseMultMod(out.data(), a.data(), b.data(), test_size, modulus, 1);
    for (size_t i = 0; i < test_size; ++i) {
        uint64_t exp = (static_cast<unsigned __int128>(a[i]) * b[i]) % modulus;
        if (out[i] != exp) return 1;
    }

    uint64_t scalar = 0xDEADBEEFCAFEULL % modulus;
    hexl_rv::EltwiseMultMod(out.data(), a.data(), scalar, test_size, modulus, 1);
    for (size_t i = 0; i < test_size; ++i) {
        uint64_t exp = (static_cast<unsigned __int128>(a[i]) * scalar) % modulus;
        if (out[i] != exp) return 1;
    }

    std::cout << "[PASS] EltwiseMultMod passed." << std::endl;
    return 0;
}
