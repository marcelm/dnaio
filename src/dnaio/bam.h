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
#define COMPILER_HAS_CONSTRUCTOR (__GNUC__ || CLANG_COMPILER_HAS(constructor))
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

static void (*decode_bam_sequence)(
    uint8_t *dest, const uint8_t *encoded_sequence, size_t length
) = decode_bam_sequence_default;

#if COMPILER_HAS_TARGET && COMPILER_HAS_CONSTRUCTOR && BUILD_IS_X86_64
__attribute__((__target__("ssse3")))
static void 
decode_bam_sequence_ssse3(uint8_t *dest, const uint8_t *encoded_sequence, size_t length) 
{

    static const uint8_t *nuc_lookup = (uint8_t *)"=ACMGRSVTWYHKDBN";
    const uint8_t *dest_end_ptr = dest + length;
    uint8_t *dest_cursor = dest;
    const uint8_t *encoded_cursor = encoded_sequence;
    const uint8_t *dest_vec_end_ptr = dest_end_ptr - (2 * sizeof(__m128i) - 1);
    __m128i nuc_lookup_vec = _mm_lddqu_si128((__m128i *)nuc_lookup);
    /* Nucleotides are encoded 4-bits per nucleotide and stored in 8-bit bytes
       as follows: |AB|CD|EF|GH|. The 4-bit codes (going from 0-15) can be used
       together with the pshufb instruction as a lookup table. The most efficient
       the upper codes (|A|C|E|G|) and one with the lower codes (|B|D|F|H|).
       The lookup can then be performed and the resulting vectors can be
       interleaved again using the unpack instructions. */
    while (dest_cursor < dest_vec_end_ptr) {
        __m128i encoded = _mm_lddqu_si128((__m128i *)encoded_cursor);
        __m128i encoded_upper = _mm_srli_epi64(encoded, 4);
        encoded_upper = _mm_and_si128(encoded_upper, _mm_set1_epi8(15));
        __m128i encoded_lower = _mm_and_si128(encoded, _mm_set1_epi8(15));
        __m128i nucs_upper = _mm_shuffle_epi8(nuc_lookup_vec, encoded_upper);
        __m128i nucs_lower = _mm_shuffle_epi8(nuc_lookup_vec, encoded_lower);
        __m128i first_nucleotides = _mm_unpacklo_epi8(nucs_upper, nucs_lower);
        __m128i second_nucleotides = _mm_unpackhi_epi8(nucs_upper, nucs_lower);
        _mm_storeu_si128((__m128i *)dest_cursor, first_nucleotides);
        _mm_storeu_si128((__m128i *)(dest_cursor + 16), second_nucleotides);
        encoded_cursor += sizeof(__m128i);
        dest_cursor += 2 * sizeof(__m128i);
    }
    decode_bam_sequence_default(dest_cursor, encoded_cursor, dest_end_ptr - dest_cursor);
}

/* Constructor functions run at dynamic link time. This checks the CPU capabilities 
   and updates the function pointer accordingly. */
__attribute__((constructor))
static void decode_bam_sequence_dispatch(void) {
    if (__builtin_cpu_supports("ssse3")) {
        decode_bam_sequence = decode_bam_sequence_ssse3;
    }
    else {
        decode_bam_sequence = decode_bam_sequence_default;
    }
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
