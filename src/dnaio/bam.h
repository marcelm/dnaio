#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>

#ifdef __SSE2__
#include "emmintrin.h"
#endif

#ifdef __SSSE3__
#include "tmmintrin.h"
#endif

static void
decode_bam_sequence(uint8_t *dest, const uint8_t *encoded_sequence, size_t length)
{
    /* Reuse a trick from sam_internal.h in htslib. Have a table to lookup
       two characters simultaneously.*/
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
    const uint8_t *dest_end_ptr = dest + length;
    uint8_t *dest_cursor = dest;
    const uint8_t *encoded_cursor = encoded_sequence;
    #ifdef __SSSE3__
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
    #endif
    /* Do two at the time until it gets to the last even address. */
    const uint8_t *dest_end_ptr_twoatatime = dest + (length & (~1ULL));
    while (dest_cursor < dest_end_ptr_twoatatime) {
        /* According to htslib, size_t cast helps the optimizer.
           Code confirmed to indeed run faster. */
        memcpy(dest_cursor, code2base + ((size_t)*encoded_cursor * 2), 2);
        dest_cursor += 2;
        encoded_cursor += 1;
    }
    assert((dest_end_ptr - dest_cursor) < 2);
    if (dest_cursor != dest_end_ptr) {
        /* There is a single encoded nuc left */
        uint8_t encoded_nucs = *encoded_cursor;
        uint8_t upper_nuc_index = encoded_nucs >> 4;
        dest_cursor[0] = nuc_lookup[upper_nuc_index];
    }
}

static void
decode_bam_qualities(uint8_t *dest, const uint8_t *encoded_qualities, size_t length)
{
    const uint8_t *end_ptr = encoded_qualities + length;
    const uint8_t *cursor = encoded_qualities;
    uint8_t *dest_cursor = dest;
    #ifdef __SSE2__
    const uint8_t *vec_end_ptr = end_ptr - sizeof(__m128i);
    while (cursor < vec_end_ptr) {
        __m128i quals = _mm_loadu_si128((__m128i *)cursor);
        __m128i phreds = _mm_add_epi8(quals, _mm_set1_epi8(33));
        _mm_storeu_si128((__m128i *)dest_cursor, phreds);
        cursor += sizeof(__m128i);
        dest_cursor += sizeof(__m128i);
    }
    #endif
    while (cursor < end_ptr) {
        *dest_cursor = *cursor + 33;
        cursor += 1;
        dest_cursor += 1;
    }
}