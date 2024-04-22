// Macros also used in htslib, very useful.
#if defined __GNUC__
#define GCC_AT_LEAST(major, minor) \
    (__GNUC__ > (major) || (__GNUC__ == (major) && __GNUC_MINOR__ >= (minor)))
#else 
# define GCC_AT_LEAST(major, minor) 0
#endif

#if defined(__clang__) && defined(__has_attribute)
#define CLANG_COMPILER_HAS(attribute) __has_attribute(attribute)
#else
#define CLANG_COMPILER_HAS(attribute) 0
#endif

#define COMPILER_HAS_TARGET (GCC_AT_LEAST(4, 8) || CLANG_COMPILER_HAS(__target__))
#define COMPILER_HAS_OPTIMIZE (GCC_AT_LEAST(4,4) || CLANG_COMPILER_HAS(optimize))

#if defined(__x86_64__) || defined(_M_X64)
#define BUILD_IS_X86_64 1
#include "immintrin.h"
#else
#define BUILD_IS_X86_64 0
#endif

#include <stdint.h>
#include <string.h>
#include <stddef.h>

static void 
decode_bam_sequence_default(uint8_t *dest, const uint8_t *encoded_sequence, size_t length)  {
    static const char code2base[512] =
        "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";
    static const uint8_t *nuc_lookup = (uint8_t *)"=ACMGRSVTWYHKDBN";
    size_t length_2 = length / 2; 
    for (size_t i=0; i < length_2; i++) {
        memcpy(dest + i*2, code2base + ((size_t)encoded_sequence[i] * 2), 2);
    }
    if (length & 1) {
        uint8_t encoded = encoded_sequence[length_2] >> 4;
        dest[(length - 1)] = nuc_lookup[encoded];
    }
}

#if COMPILER_HAS_TARGET && BUILD_IS_X86_64
__attribute__((__target__("ssse3")))
static void 
decode_bam_sequence_ssse3(uint8_t *dest, const uint8_t *encoded_sequence, size_t length) 
{

    static const uint8_t *nuc_lookup = (uint8_t *)"=ACMGRSVTWYHKDBN";
    const uint8_t *dest_end_ptr = dest + length;
    uint8_t *dest_cursor = dest;
    const uint8_t *encoded_cursor = encoded_sequence;
    const uint8_t *dest_vec_end_ptr = dest_end_ptr - (2 * sizeof(__m128i));
    __m128i first_upper_shuffle = _mm_setr_epi8(
        0, 0xff, 1, 0xff, 2, 0xff, 3, 0xff, 4, 0xff, 5, 0xff, 6, 0xff, 7, 0xff);
    __m128i first_lower_shuffle = _mm_setr_epi8(
        0xff, 0, 0xff, 1, 0xff, 2, 0xff, 3, 0xff, 4, 0xff, 5, 0xff, 6, 0xff, 7);
    __m128i second_upper_shuffle = _mm_setr_epi8(
        8, 0xff, 9, 0xff, 10, 0xff, 11, 0xff, 12, 0xff, 13, 0xff, 14, 0xff, 15, 0xff);
    __m128i second_lower_shuffle = _mm_setr_epi8(
        0xff, 8, 0xff, 9, 0xff, 10, 0xff, 11, 0xff, 12, 0xff, 13, 0xff, 14, 0xff, 15);
    __m128i nuc_lookup_vec = _mm_lddqu_si128((__m128i *)nuc_lookup);
    /* Work on 16 encoded characters at the time resulting in 32 decoded characters
       Examples are given for 8 encoded characters A until H to keep it readable.
        Encoded stored as |AB|CD|EF|GH|
        Shuffle into |AB|00|CD|00|EF|00|GH|00| and
                     |00|AB|00|CD|00|EF|00|GH|
        Shift upper to the right resulting into
                     |0A|B0|0C|D0|0E|F0|0G|H0| and
                     |00|AB|00|CD|00|EF|00|GH|
        Merge with or resulting into (X stands for garbage)
                     |0A|XB|0C|XD|0E|XF|0G|XH|
        Bitwise and with 0b1111 leads to:
                     |0A|0B|0C|0D|0E|0F|0G|0H|
        We can use the resulting 4-bit integers as indexes for the shuffle of
        the nucleotide lookup. */
    while (dest_cursor < dest_vec_end_ptr) {
        __m128i encoded = _mm_lddqu_si128((__m128i *)encoded_cursor);

        __m128i first_upper = _mm_shuffle_epi8(encoded, first_upper_shuffle);
        __m128i first_lower = _mm_shuffle_epi8(encoded, first_lower_shuffle);
        __m128i shifted_first_upper = _mm_srli_epi64(first_upper, 4);
        __m128i first_merged = _mm_or_si128(shifted_first_upper, first_lower);
        __m128i first_indexes = _mm_and_si128(first_merged, _mm_set1_epi8(0b1111));
        __m128i first_nucleotides = _mm_shuffle_epi8(nuc_lookup_vec, first_indexes);
        _mm_storeu_si128((__m128i *)dest_cursor, first_nucleotides);

        __m128i second_upper = _mm_shuffle_epi8(encoded, second_upper_shuffle);
        __m128i second_lower = _mm_shuffle_epi8(encoded, second_lower_shuffle);
        __m128i shifted_second_upper = _mm_srli_epi64(second_upper, 4);
        __m128i second_merged = _mm_or_si128(shifted_second_upper, second_lower);
        __m128i second_indexes = _mm_and_si128(second_merged, _mm_set1_epi8(0b1111));
        __m128i second_nucleotides = _mm_shuffle_epi8(nuc_lookup_vec, second_indexes);
        _mm_storeu_si128((__m128i *)(dest_cursor + 16), second_nucleotides);

        encoded_cursor += sizeof(__m128i);
        dest_cursor += 2 * sizeof(__m128i);
    }
    decode_bam_sequence_default(dest_cursor, encoded_cursor, dest_end_ptr - dest_cursor);
}

static void (*decode_bam_sequence)(
    uint8_t *dest, const uint8_t *encoded_sequence, size_t length);

/* Simple dispatcher function, updates the function pointer after testing the
   CPU capabilities. After this, the dispatcher function is not needed anymore. */
static void decode_bam_sequence_dispatch(
        uint8_t *dest, const uint8_t *encoded_sequence, size_t length) {
    if (__builtin_cpu_supports("ssse3")) {
        decode_bam_sequence = decode_bam_sequence_ssse3;
    }
    else {
        decode_bam_sequence = decode_bam_sequence_default;
    }
    decode_bam_sequence(dest, encoded_sequence, length);
}

static void (*decode_bam_sequence)(
    uint8_t *dest, const uint8_t *encoded_sequence, size_t length
) = decode_bam_sequence_dispatch;

#else
static inline void decode_bam_sequence(uint8_t *dest, const uint8_t *encoded_sequence, size_t length) 
{
    decode_bam_sequence_default(dest, encoded_sequence, length);
}
#endif 

// Code is simple enough to be auto vectorized.
#if COMPILER_HAS_OPTIMIZE
__attribute__((optimize("O3")))
#endif
static void 
decode_bam_qualities(
    uint8_t *restrict dest,
    const uint8_t *restrict encoded_qualities,
    size_t length)
{
    for (size_t i=0; i<length; i++) {
        dest[i] = encoded_qualities[i] + 33;
    }
}
