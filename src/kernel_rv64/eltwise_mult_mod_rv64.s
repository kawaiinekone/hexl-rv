.text
.balign 4
.global EltwiseMultModKernel_rv64
.type EltwiseMultModKernel_rv64, @function

# void EltwiseMultModKernel_rv64(
#     uint64_t* result,        // a0
#     const uint64_t* op1,     // a1
#     const uint64_t* op2,     // a2
#     uint64_t n,              // a3
#     uint64_t modulus,        // a4
#     uint64_t barrett_factor, // a5
#     uint64_t shift           // a6
# );

EltwiseMultModKernel_rv64:
    beqz    a3, .Ldone

    li      t1, 64
    sub     t1, t1, a6         # t1 = 64 - shift

.Lloop:
    ld      t2, 0(a1)          # t2 = op1[i]
    ld      t3, 0(a2)          # t3 = op2[i]
    addi    a1, a1, 8
    addi    a2, a2, 8

    # 128-bit product
    mul     t4, t2, t3         # prod_lo
    mulhu   t5, t2, t3         # prod_hi

    # Combine shifted components
    srl     t6, t4, a6         # prod_lo >> shift
    sll     t5, t5, t1         # prod_hi << (64 - shift)
    or      t6, t6, t5         # c1

    # Quotient estimation
    mulhu   t6, t6, a5         # q_hat = mulhu(c1, barrett_factor)
    mul     t6, t6, a4         # q_hat * modulus
    sub     t4, t4, t6         # rem = prod_lo - q_hat * modulus

    # Branchless reduction (rem - (rem >= modulus ? modulus : 0))
    sltu    t0, t4, a4         # t0 = (rem < modulus) ? 1 : 0
    addi    t0, t0, -1         # mask = (rem < modulus) ? 0 : ~0ULL
    and     t0, t0, a4         # sub_val = mask & modulus
    sub     t4, t4, t0         # rem = rem - sub_val

    sd      t4, 0(a0)          # store result[i]
    addi    a0, a0, 8

    addi    a3, a3, -1
    bnez    a3, .Lloop

.Ldone:
    ret
