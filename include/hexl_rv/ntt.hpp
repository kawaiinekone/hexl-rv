#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>

namespace hexl_rv {

class NTT {
public:
    NTT(uint64_t degree, uint64_t modulus, uint64_t root_of_unity = 0);

    void ComputeForward(uint64_t* operand, uint64_t input_mod_factor = 1) const;
    void ComputeInverse(uint64_t* operand) const;

    uint64_t GetDegree() const { return m_degree; }
    uint64_t GetModulus() const { return m_modulus; }

private:
    uint64_t m_degree;
    uint64_t m_modulus;
    uint64_t m_degree_bits;
    uint64_t m_inv_degree;
    std::vector<uint64_t> m_root_of_unity_powers;
    std::vector<uint64_t> m_inv_root_of_unity_powers;

    void PrecomputeTwiddleFactors(uint64_t root_of_unity);
};

} // namespace hexl_rv
