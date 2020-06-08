.intel_syntax noprefix

.global uintbig_1
.global _uintbig_1
uintbig_1: _uintbig_1:
    .quad 1, 0, 0, 0


.text

.global uintbig_set
.global _uintbig_set
uintbig_set: _uintbig_set:
    cld
    mov rax, rsi
    stosq
    xor rax, rax
    mov rcx, 3
    rep stosq
    ret


.global uintbig_bit
.global _uintbig_bit
uintbig_bit: _uintbig_bit:
    mov rcx, rsi
    and rcx, 0x3f
    shr rsi, 6
    mov rax, [rdi + 8*rsi]
    shr rax, cl
    and rax, 1
    ret


.global uintbig_add3
.global _uintbig_add3
uintbig_add3: _uintbig_add3:
    mov rax, [rsi +  0]
    add rax, [rdx +  0]
    mov [rdi +  0], rax
    .set k, 1
    .rept 3
        mov rax, [rsi + 8*k]
        adc rax, [rdx + 8*k]
        mov [rdi + 8*k], rax
        .set k, k+1
    .endr
    setc al
    movzx rax, al
    ret

.global uintbig_sub3
.global _uintbig_sub3
uintbig_sub3: _uintbig_sub3:
    mov rax, [rsi +  0]
    sub rax, [rdx +  0]
    mov [rdi +  0], rax
    .set k, 1
    .rept 3
        mov rax, [rsi + 8*k]
        sbb rax, [rdx + 8*k]
        mov [rdi + 8*k], rax
        .set k, k+1
    .endr
    setc al
    movzx rax, al
    ret


.global uintbig_mul3_64
.global _uintbig_mul3_64
uintbig_mul3_64: _uintbig_mul3_64:

    mulx r10, rax, [rsi +  0]
    mov [rdi +  0], rax

    mulx r11, rax, [rsi +  8]
    add  rax, r10
    mov [rdi +  8], rax

    mulx r10, rax, [rsi + 16]
    adcx rax, r11
    mov [rdi + 16], rax

    mulx r11, rax, [rsi + 24]
    adcx rax, r10
    mov [rdi + 24], rax

    ret

.global uintbig_div3_64
.global _uintbig_div3_64
uintbig_div3_64: _uintbig_div3_64:
    mov r10, rdx
    mov rdx, 0

    mov rax, [rsi + 24]
    div r10
    mov [rdi + 24], rax

    mov rax, [rsi + 16]
    div r10
    mov [rdi + 16], rax

    mov rax, [rsi +  8]
    div r10
    mov [rdi +  8], rax

    mov rax, [rsi +  0]
    div r10
    mov [rdi +  0], rax

    mov rax, rdx

    ret
