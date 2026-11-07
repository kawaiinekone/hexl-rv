#include "hexl_rv/eltwise_add_sub_mod.hpp"

namespace hexl_rv {

void EltwiseAddMod(uint64_t* result,
                   const uint64_t* operand1,
                   const uint64_t* operand2,
                   uint64_t n,
                   uint64_t modulus) {
    if (!result || !operand1 || !operand2 || n == 0 || modulus <= 1) return;

    for (uint64_t i = 0; i < n; ++i) {
        uint64_t a = operand1[i];
        uint64_t b = operand2[i];
        if (a >= modulus) a -= modulus;
        if (b >= modulus) b -= modulus;

        uint64_t sum = a + b;
        result[i] = (sum >= modulus) ? (sum - modulus) : sum;
    }
}

void EltwiseAddMod(uint64_t* result,
                   const uint64_t* operand1,
                   uint64_t scalar,
                   uint64_t n,
                   uint64_t modulus) {
    if (!result || !operand1 || n == 0 || modulus <= 1) return;

    uint64_t b = (scalar >= modulus) ? (scalar % modulus) : scalar;
    for (uint64_t i = 0; i < n; ++i) {
        uint64_t a = operand1[i];
        if (a >= modulus) a -= modulus;

        uint64_t sum = a + b;
        result[i] = (sum >= modulus) ? (sum - modulus) : sum;
    }
}

void EltwiseSubMod(uint64_t* result,
                   const uint64_t* operand1,
                   const uint64_t* operand2,
                   uint64_t n,
                   uint64_t modulus) {
    if (!result || !operand1 || !operand2 || n == 0 || modulus <= 1) return;

    for (uint64_t i = 0; i < n; ++i) {
        uint64_t a = operand1[i];
        uint64_t b = operand2[i];
        if (a >= modulus) a -= modulus;
        if (b >= modulus) b -= modulus;

        result[i] = (a >= b) ? (a - b) : (a + modulus - b);
    }
}

void EltwiseSubMod(uint64_t* result,
                   const uint64_t* operand1,
                   uint64_t scalar,
                   uint64_t n,
                   uint64_t modulus) {
    if (!result || !operand1 || n == 0 || modulus <= 1) return;

    uint64_t b = (scalar >= modulus) ? (scalar % modulus) : scalar;
    for (uint64_t i = 0; i < n; ++i) {
        uint64_t a = operand1[i];
        if (a >= modulus) a -= modulus;

        result[i] = (a >= b) ? (a - b) : (a + modulus - b);
    }
}

} // namespace hexl_rv
