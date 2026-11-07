#pragma once

#include <cstdint>
#include <cstddef>

namespace hexl_rv {

void EltwiseAddMod(uint64_t* result,
                   const uint64_t* operand1,
                   const uint64_t* operand2,
                   uint64_t n,
                   uint64_t modulus);

void EltwiseAddMod(uint64_t* result,
                   const uint64_t* operand1,
                   uint64_t scalar,
                   uint64_t n,
                   uint64_t modulus);

void EltwiseSubMod(uint64_t* result,
                   const uint64_t* operand1,
                   const uint64_t* operand2,
                   uint64_t n,
                   uint64_t modulus);

void EltwiseSubMod(uint64_t* result,
                   const uint64_t* operand1,
                   uint64_t scalar,
                   uint64_t n,
                   uint64_t modulus);

} // namespace hexl_rv
