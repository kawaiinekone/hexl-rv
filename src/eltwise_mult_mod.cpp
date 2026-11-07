#include "hexl_rv/eltwise_mult_mod.hpp"
#include "hexl_rv/number_theory.hpp"

namespace hexl_rv {

namespace {

template <int InputModFactor>
void EltwiseMultModKernel(uint64_t* result,
                          const uint64_t* operand1,
                          const uint64_t* operand2,
                          uint64_t n,
                          uint64_t modulus) {
    const uint64_t twice_modulus = 2 * modulus;
    BarrettFactor barrett(modulus);
    const uint64_t mu = barrett.barrett_factor;
    const uint64_t shift = barrett.prod_right_shift;

    for (uint64_t i = 0; i < n; ++i) {
        uint64_t x = operand1[i];
        uint64_t y = operand2[i];

        if constexpr (InputModFactor == 2) {
            if (x >= modulus) x -= modulus;
            if (y >= modulus) y -= modulus;
        } else if constexpr (InputModFactor == 4) {
            if (x >= twice_modulus) x -= twice_modulus;
            if (x >= modulus) x -= modulus;
            if (y >= twice_modulus) y -= twice_modulus;
            if (y >= modulus) y -= modulus;
        }

        result[i] = MultiplyModBarrett(x, y, modulus, mu, shift);
    }
}

template <int InputModFactor>
void EltwiseMultModScalarKernel(uint64_t* result,
                                const uint64_t* operand1,
                                uint64_t scalar,
                                uint64_t n,
                                uint64_t modulus) {
    const uint64_t twice_modulus = 2 * modulus;
    BarrettFactor barrett(modulus);
    const uint64_t mu = barrett.barrett_factor;
    const uint64_t shift = barrett.prod_right_shift;

    uint64_t y = scalar;
    if constexpr (InputModFactor == 2) {
        if (y >= modulus) y -= modulus;
    } else if constexpr (InputModFactor == 4) {
        if (y >= twice_modulus) y -= twice_modulus;
        if (y >= modulus) y -= modulus;
    }

    for (uint64_t i = 0; i < n; ++i) {
        uint64_t x = operand1[i];

        if constexpr (InputModFactor == 2) {
            if (x >= modulus) x -= modulus;
        } else if constexpr (InputModFactor == 4) {
            if (x >= twice_modulus) x -= twice_modulus;
            if (x >= modulus) x -= modulus;
        }

        result[i] = MultiplyModBarrett(x, y, modulus, mu, shift);
    }
}

} // namespace

void EltwiseMultMod(uint64_t* result,
                    const uint64_t* operand1,
                    const uint64_t* operand2,
                    uint64_t n,
                    uint64_t modulus,
                    uint64_t input_mod_factor) {
    if (!result || !operand1 || !operand2 || n == 0 || modulus <= 1) return;

    switch (input_mod_factor) {
        case 1:  EltwiseMultModKernel<1>(result, operand1, operand2, n, modulus); break;
        case 2:  EltwiseMultModKernel<2>(result, operand1, operand2, n, modulus); break;
        case 4:  EltwiseMultModKernel<4>(result, operand1, operand2, n, modulus); break;
        default: EltwiseMultModKernel<1>(result, operand1, operand2, n, modulus); break;
    }
}

void EltwiseMultMod(uint64_t* result,
                    const uint64_t* operand1,
                    uint64_t scalar,
                    uint64_t n,
                    uint64_t modulus,
                    uint64_t input_mod_factor) {
    if (!result || !operand1 || n == 0 || modulus <= 1) return;

    switch (input_mod_factor) {
        case 1:  EltwiseMultModScalarKernel<1>(result, operand1, scalar, n, modulus); break;
        case 2:  EltwiseMultModScalarKernel<2>(result, operand1, scalar, n, modulus); break;
        case 4:  EltwiseMultModScalarKernel<4>(result, operand1, scalar, n, modulus); break;
        default: EltwiseMultModScalarKernel<1>(result, operand1, scalar, n, modulus); break;
    }
}

} // namespace hexl_rv
