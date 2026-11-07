#include "hexl_rv/eltwise_add_sub_mod.hpp"
#include <iostream>
#include <vector>
#include <random>

int main() {
    std::cout << "[TEST] Running EltwiseAddMod & EltwiseSubMod..." << std::endl;
    const size_t n = 2048;
    const uint64_t q = (1ULL << 40) - 87;
    std::mt19937_64 rng(999);

    std::vector<uint64_t> op1(n), op2(n), res_add(n), res_sub(n);
    for (size_t i = 0; i < n; ++i) {
        op1[i] = rng() % q;
        op2[i] = rng() % q;
    }

    hexl_rv::EltwiseAddMod(res_add.data(), op1.data(), op2.data(), n, q);
    hexl_rv::EltwiseSubMod(res_sub.data(), op1.data(), op2.data(), n, q);

    for (size_t i = 0; i < n; ++i) {
        uint64_t exp_add = (op1[i] + op2[i]) % q;
        uint64_t exp_sub = (op1[i] >= op2[i]) ? (op1[i] - op2[i]) : (op1[i] + q - op2[i]);
        if (res_add[i] != exp_add || res_sub[i] != exp_sub) return 1;
    }

    std::cout << "[PASS] EltwiseAddMod & EltwiseSubMod passed." << std::endl;
    return 0;
}
