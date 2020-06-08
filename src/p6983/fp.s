.intel_syntax noprefix

.set pbytes,32
.set plimbs,4
.global p
.global _p
p: _p:
    .quad 0x50ca4291ffffffff, 0xd1b004f94a5952c9, 0xb25c76a437728f3b, 0xa3091565b678a990

.inv_min_p_mod_r: /* -p^-1 mod 2^64 */
    .quad 0x50ca429200000001

.global fp_0
.global _fp_0
fp_0: _fp_0:
    .zero pbytes

.global fp_1
.global _fp_1
fp_1: _fp_1: /* 2^256 mod p */
    .quad 0xaf35bd6e00000001, 0x2e4ffb06b5a6ad36, 0x4da3895bc88d70c4, 0x5cf6ea9a4987566f

.r_squared_mod_p: /* (2^256)^2 mod p */
    .quad 0xce1eebd1f8f7c152, 0x3d0e272c4d770926, 0x868e6f58f9294cfb, 0x578b875020f84a87

.p_minus_2:
    .quad 0x50ca4291fffffffd, 0xd1b004f94a5952c9, 0xb25c76a437728f3b, 0xa3091565b678a990

.p_minus_1_halves:
    .quad 0xa8652148ffffffff, 0xe8d8027ca52ca964, 0x592e3b521bb9479d, 0x51848ab2db3c54c8

.p_plus_1_quarter:
	.quad 0x543290a480000000, 0xf46c013e529654b2, 0x2c971da90ddca3ce, 0x28c245596d9e2a64

.data
.global fp_mul_count
.global _fp_mul_count
fp_mul_count: _fp_mul_count:
    .quad 0

.text
.p2align 4,,15

.global fp_copy
.global _fp_copy
fp_copy: _fp_copy:
    cld
    mov rcx, plimbs
    rep movsq
    ret

.global fp_set
.global _fp_set
fp_set: _fp_set:
    push rdi
    call uintbig_set
    pop rdi
    mov rsi, rdi
    jmp fp_enc

.global fp_cswap
.global _fp_cswap
fp_cswap: _fp_cswap:
    movzx rax, dl
    neg rax
    .set k, 0
    .rept plimbs
        mov rcx, [rdi + 8*k]
        mov rdx, [rsi + 8*k]

        mov r8, rcx
        xor r8, rdx
        and r8, rax

        xor rcx, r8
        xor rdx, r8

        mov [rdi + 8*k], rcx
        mov [rsi + 8*k], rdx

        .set k, k+1
    .endr
    ret

.reduce_once:
    push rbp
    mov rbp, rdi

    mov rdi, [rbp +  0]
    sub rdi, [rip + p +  0]
    mov rsi, [rbp +  8]
    sbb rsi, [rip + p +  8]
    mov rdx, [rbp + 16]
    sbb rdx, [rip + p + 16]
    mov rcx, [rbp + 24]
    sbb rcx, [rip + p + 24]
    sbb rax, 0 /* handle carry from caller */

    setnc al
    movzx rax, al
    neg rax

.macro cswap2, r, m
    xor \r, \m
    and \r, rax
    xor \m, \r
.endm

    cswap2 rdi, [rbp +  0]
    cswap2 rsi, [rbp +  8]
    cswap2 rdx, [rbp + 16]
    cswap2 rcx, [rbp + 24]

    pop rbp
    ret

.global fp_add2
.global _fp_add2
fp_add2: _fp_add2:
    mov rdx, rdi
.global fp_add3
.global _fp_add3
fp_add3: _fp_add3:
    push rdi
    call uintbig_add3 /* may leave a carry in rax */
    pop rdi
    jmp .reduce_once

.global fp_sub2
.global _fp_sub2
fp_sub2: _fp_sub2:
  mov rdx, rdi
  xchg rsi, rdx
.global fp_sub3
.global _fp_sub3
fp_sub3: _fp_sub3:
    push rdi
    call uintbig_sub3
    pop rdi
    xor rsi, rsi
    xor rdx, rdx
    xor rcx, rcx
    test rax, rax
    cmovnz rax, [rip + p +  0]
    cmovnz rsi, [rip + p +  8]
    cmovnz rdx, [rip + p + 16]
    cmovnz rcx, [rip + p + 24]
    add [rdi +  0], rax
    adc [rdi +  8], rsi
    adc [rdi + 16], rdx
    adc [rdi + 24], rcx
    ret


/* Montgomery arithmetic */

.global fp_enc
.global _fp_enc
fp_enc: _fp_enc:
    lea rdx, [rip + .r_squared_mod_p]
    jmp fp_mul3

.global fp_dec
.global _fp_dec
fp_dec: _fp_dec:
    lea rdx, [rip + uintbig_1]
    jmp fp_mul3

.global fp_mul2
.global _fp_mul2
fp_mul2: _fp_mul2:
  mov rdx, rdi
.global fp_mul3
.global _fp_mul3
fp_mul3: _fp_mul3:
    push rbp
    push rbx
    push r12

    push rdi

    /* inc qword ptr fp_mul_count */

    mov rdi, rsi
    mov rsi, rdx

    xor r8,  r8
    xor r9,  r9
    xor r10, r10
    xor r11, r11
    xor r12, r12
    xor rbp, rbp

    /* flags are already cleared */

.macro MULSTEP, k, r0, r1, r2, r3, r4, r5

    mov rdx, [rsi +  0]
    mulx rcx, rdx, [rdi + 8*\k]
    add rdx, \r0
    mulx rcx, rdx, [rip + .inv_min_p_mod_r]

    xor rax, rax /* clear flags */

    mulx rbx, rax, [rip + p +  0]
    adox \r0, rax

    mulx rcx, rax, [rip + p +  8]
    adcx \r1, rbx
    adox \r1, rax

    mulx rbx, rax, [rip + p + 16]
    adcx \r2, rcx
    adox \r2, rax

    mulx rcx, rax, [rip + p + 24]
    adcx \r3, rbx
    adox \r3, rax

    mov rax, 0
    adcx \r4, rcx
    adox \r4, rax

    adcx \r5, rax
    adox \r5, rax

    mov rdx, [rdi + 8*\k]

    xor rax, rax /* clear flags */

    mulx rbx, rax, [rsi +  0]
    adox \r0, rax

    mulx rcx, rax, [rsi +  8]
    adcx \r1, rbx
    adox \r1, rax

    mulx rbx, rax, [rsi + 16]
    adcx \r2, rcx
    adox \r2, rax

    mulx rcx, rax, [rsi + 24]
    adcx \r3, rbx
    adox \r3, rax

    mov rax, 0
    adcx \r4, rcx
    adox \r4, rax

    adcx \r5, rax
    adox \r5, rax

.endm

    MULSTEP 0, r8,  r9,  r10, r11, r12, rbp
    MULSTEP 1, r9,  r10, r11, r12, rbp, r8
    MULSTEP 2, r10, r11, r12, rbp, r8,  r9
    MULSTEP 3, r11, r12, rbp, r8,  r9,  r10

    pop rdi

    mov [rdi +  0], r12
    mov [rdi +  8], rbp
    mov [rdi + 16], r8
    mov [rdi + 24], r9
    mov rax, r10

    pop r12
    pop rbx
    pop rbp
    jmp .reduce_once

.global fp_sq1
.global _fp_sq1
fp_sq1: _fp_sq1:
    mov rsi, rdi
.global fp_sq2
.global _fp_sq2
fp_sq2: _fp_sq2:
    /* TODO implement optimized Montgomery squaring */
    mov rdx, rsi
    jmp fp_mul3

/* (obviously) not constant time in the exponent! */
.fp_pow:
    push rbx
    mov rbx, rsi
    push r12
    push r13
    push rdi
    sub rsp, pbytes

    mov rsi, rdi
    mov rdi, rsp
    call fp_copy

    mov rdi, [rsp + pbytes]
    lea rsi, [rip + fp_1]
    call fp_copy

.macro POWSTEP, k
        mov r13, [rbx + 8*\k]
        xor r12, r12

        2:
        test r13, 1
        jz 1f

        mov rdi, [rsp + pbytes]
        mov rsi, rsp
        call fp_mul2

        1:
        mov rdi, rsp
        call fp_sq1

        shr r13

        inc r12
        test r12, 64
        jz 2b
.endm

    POWSTEP 0
    POWSTEP 1
    POWSTEP 2
    POWSTEP 3

    add rsp, pbytes+8
    pop r13
    pop r12
    pop rbx
    ret

/* TODO use a better addition chain? */
.global fp_inv
.global _fp_inv
fp_inv: _fp_inv:
    lea rsi, [rip + .p_minus_2]
    jmp .fp_pow

.global fp_sqrt
.global _fp_sqrt
fp_sqrt: _fp_sqrt:
    lea rsi, [rip + .p_plus_1_quarter]
    jmp .fp_pow

/* TODO use a better addition chain? */
.global fp_issquare
.global _fp_issquare
fp_issquare: _fp_issquare:
    push rdi
    lea rsi, [rip + .p_minus_1_halves]
    call .fp_pow
    pop rdi

    xor rax, rax
    .set k, 0
    .rept plimbs
        mov rsi, [rdi + 8*k]
        xor rsi, [rip + fp_1 + 8*k]
        or rax, rsi
        .set k, k+1
    .endr
    test rax, rax
    setz al
    movzx rax, al
    ret

/* not constant time (but this shouldn't leak anything of importance) */
.global fp_random
.global _fp_random
fp_random: _fp_random:

    push rdi
    mov rsi, pbytes
    call _randombytes
    pop rdi

    .set k, plimbs-1
    .rept plimbs
        mov rax, [rip + p + 8*k]
        cmp [rdi + 8*k], rax
        ja fp_random
        jb 0f
        .set k, k-1
    .endr
    jmp fp_random
    0:
    ret
