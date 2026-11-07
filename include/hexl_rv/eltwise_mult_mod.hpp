#pragma once

#include <cstdint>
#include <cstddef>

namespace hexl_rv {

void EltwiseMultMod(uint64_t* result,
                    const uint64_t* operand1,
                    const uint64_t* operand2,
                    uint64_t n,
                    uint64_t modulus,
                    uint64_t input_mod_factor = 1);

void EltwiseMultMod(uint64_t* result,
                    const uint64_t* operand1,
                    uint64_t scalar,
                    uint64_t n,
                    uint64_t modulus,
                    uint64_t input_mod_factor = 1);

} // namespace hexl_rv
