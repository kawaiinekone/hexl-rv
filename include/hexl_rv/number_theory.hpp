#pragma once

#include <cstdint>
#include <cstddef>
#include <limits>

namespace hexl_rv {

inline void MultiplyUInt64(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo) {
#if defined(__SIZEOF_INT128__)
    unsigned __int128 prod = static_cast<unsigned __int128>(a) * b;
    *lo = static_cast<uint64_t>(prod);
    *hi = static_cast<uint64_t>(prod >> 64);
#else
    uint64_t a_lo = static_cast<uint32_t>(a);
    uint64_t a_hi = a >> 32;
    uint64_t b_lo = static_cast<uint32_t>(b);
    uint64_t b_hi = b >> 32;

    uint64_t p0 = a_lo * b_lo;
    uint64_t p1 = a_lo * b_hi;
    uint64_t p2 = a_hi * b_lo;
    uint64_t p3 = a_hi * b_hi;

    uint64_t cy = ((p0 >> 32) + static_cast<uint32_t>(p1) + static_cast<uint32_t>(p2)) >> 32;
    *lo = p0 + (p1 << 32) + (p2 << 32);
    *hi = p3 + (p1 >> 32) + (p2 >> 32) + cy;
#endif
}

inline uint64_t MultiplyUInt64Hi(uint64_t a, uint64_t b) {
#if defined(__SIZEOF_INT128__)
    unsigned __int128 prod = static_cast<unsigned __int128>(a) * b;
    return static_cast<uint64_t>(prod >> 64);
#else
    uint64_t hi, lo;
    MultiplyUInt64(a, b, &hi, &lo);
    return hi;
#endif
}

inline uint64_t MSB(uint64_t x) {
    if (x == 0) return 0;
    return 63 - __builtin_clzll(x);
}

struct BarrettFactor {
    uint64_t modulus;
    uint64_t barrett_factor;
    uint64_t prod_right_shift;

    explicit BarrettFactor(uint64_t mod = 0) 
        : modulus(mod), barrett_factor(0), prod_right_shift(0) {
        if (modulus > 1) {
            uint64_t ceil_log_mod = MSB(modulus) + 1;
            prod_right_shift = (ceil_log_mod >= 2) ? (ceil_log_mod - 2) : 0;
            unsigned __int128 numerator = static_cast<unsigned __int128>(1) << (ceil_log_mod + 62);
            barrett_factor = static_cast<uint64_t>(numerator / modulus);
        }
    }
};

inline uint64_t MultiplyModBarrett(uint64_t a, uint64_t b, uint64_t modulus, 
                                   uint64_t barrett_factor, uint64_t prod_right_shift) {
    uint64_t prod_hi, prod_lo;
    MultiplyUInt64(a, b, &prod_hi, &prod_lo);

    uint64_t c1 = (prod_right_shift == 0) 
        ? prod_lo 
        : (prod_lo >> prod_right_shift) | (prod_hi << (64 - prod_right_shift));

    uint64_t c2_hi = MultiplyUInt64Hi(c1, barrett_factor);
    uint64_t q_hat = c2_hi;

    uint64_t rem = prod_lo - q_hat * modulus;
    return (rem >= modulus) ? rem - modulus : rem;
}

inline uint64_t ReduceMod(uint64_t x, uint64_t modulus) {
    return x >= modulus ? x % modulus : x;
}

inline uint64_t PowerMod(uint64_t base, uint64_t exp, uint64_t modulus) {
    uint64_t res = 1;
    base = base % modulus;
    while (exp > 0) {
        if (exp & 1) {
            res = static_cast<uint64_t>((static_cast<unsigned __int128>(res) * base) % modulus);
        }
        base = static_cast<uint64_t>((static_cast<unsigned __int128>(base) * base) % modulus);
        exp >>= 1;
    }
    return res;
}

inline uint64_t InverseMod(uint64_t a, uint64_t modulus) {
    int64_t t = 0, newt = 1;
    int64_t r = static_cast<int64_t>(modulus), newr = static_cast<int64_t>(a % modulus);

    while (newr != 0) {
        int64_t q = r / newr;
        int64_t temp_t = t - q * newt;
        t = newt;
        newt = temp_t;

        int64_t temp_r = r - q * newr;
        r = newr;
        newr = temp_r;
    }

    if (r > 1) return 0;
    if (t < 0) t += modulus;
    return static_cast<uint64_t>(t);
}

inline uint64_t ReverseBits(uint64_t x, uint64_t bit_width) {
    uint64_t rev = 0;
    for (uint64_t i = 0; i < bit_width; ++i) {
        rev = (rev << 1) | (x & 1);
        x >>= 1;
    }
    return rev;
}

} // namespace hexl_rv
