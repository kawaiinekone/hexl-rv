#include "hexl_rv/ntt.hpp"
#include "hexl_rv/number_theory.hpp"
#include <stdexcept>

namespace hexl_rv {

NTT::NTT(uint64_t degree, uint64_t modulus, uint64_t root_of_unity)
    : m_degree(degree), m_modulus(modulus) {
    if (m_degree == 0 || (m_degree & (m_degree - 1)) != 0) {
        throw std::invalid_argument("NTT degree must be a power of 2");
    }
    if (m_modulus <= 1) {
        throw std::invalid_argument("NTT modulus must be > 1");
    }

    m_degree_bits = MSB(m_degree);
    m_inv_degree = InverseMod(m_degree, m_modulus);

    PrecomputeTwiddleFactors(root_of_unity);
}

void NTT::PrecomputeTwiddleFactors(uint64_t root_of_unity) {
    m_root_of_unity_powers.resize(m_degree);
    m_inv_root_of_unity_powers.resize(m_degree);

    uint64_t root = root_of_unity;
    uint64_t inv_root = InverseMod(root, m_modulus);

    for (uint64_t i = 0; i < m_degree; ++i) {
        uint64_t rev_i = ReverseBits(i, m_degree_bits);
        m_root_of_unity_powers[i] = PowerMod(root, rev_i, m_modulus);
        m_inv_root_of_unity_powers[i] = PowerMod(inv_root, rev_i, m_modulus);
    }
}

void NTT::ComputeForward(uint64_t* operand, uint64_t input_mod_factor) const {
    if (!operand) return;

    if (input_mod_factor > 1) {
        for (uint64_t i = 0; i < m_degree; ++i) {
            operand[i] = ReduceMod(operand[i], m_modulus);
        }
    }

    uint64_t t = m_degree;
    for (uint64_t m = 1; m < m_degree; m <<= 1) {
        t >>= 1;
        for (uint64_t i = 0; i < m; ++i) {
            uint64_t j1 = 2 * i * t;
            uint64_t j2 = j1 + t - 1;
            uint64_t W = m_root_of_unity_powers[m + i];

            for (uint64_t j = j1; j <= j2; ++j) {
                uint64_t u = operand[j];
                uint64_t v = static_cast<uint64_t>((static_cast<unsigned __int128>(operand[j + t]) * W) % m_modulus);

                uint64_t add = u + v;
                operand[j] = (add >= m_modulus) ? (add - m_modulus) : add;

                operand[j + t] = (u >= v) ? (u - v) : (u + m_modulus - v);
            }
        }
    }
}

void NTT::ComputeInverse(uint64_t* operand) const {
    if (!operand) return;

    uint64_t t = 1;
    for (uint64_t m = m_degree >> 1; m >= 1; m >>= 1) {
        for (uint64_t i = 0; i < m; ++i) {
            uint64_t j1 = 2 * i * t;
            uint64_t j2 = j1 + t - 1;
            uint64_t W = m_inv_root_of_unity_powers[m + i];

            for (uint64_t j = j1; j <= j2; ++j) {
                uint64_t u = operand[j];
                uint64_t v = operand[j + t];

                uint64_t add = u + v;
                operand[j] = (add >= m_modulus) ? (add - m_modulus) : add;

                uint64_t sub = (u >= v) ? (u - v) : (u + m_modulus - v);
                operand[j + t] = static_cast<uint64_t>((static_cast<unsigned __int128>(sub) * W) % m_modulus);
            }
        }
        t <<= 1;
    }

    for (uint64_t i = 0; i < m_degree; ++i) {
        operand[i] = static_cast<uint64_t>((static_cast<unsigned __int128>(operand[i]) * m_inv_degree) % m_modulus);
    }
}

} // namespace hexl_rv
