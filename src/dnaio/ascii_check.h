#include <stddef.h>
#include <stdint.h>
#ifdef __SSE2__
#include "emmintrin.h"
#endif

#define ASCII_MASK_8BYTE 0x8080808080808080ULL
#define ASCII_MASK_1BYTE 0x80

/**
 * @brief Check if a string of given length only contains ASCII characters.
 *
 * @param string A char pointer to the start of the string.
 * @param length The length of the string. This funtion does not check for
 *               terminating NULL bytes.
 * @returns 1 if the string is ASCII-only, 0 otherwise.
 */
static int
string_is_ascii(const char * string, size_t length) {
    // By performing bitwise OR on all characters in 8-byte chunks (16-byte
    // with SSE2) we can
    // determine ASCII status in a non-branching (except the loops) fashion.
    uint64_t all_chars = 0;
    const char *cursor = string;
    const char *string_end_ptr = string + length;
    const char *string_8b_end_ptr = string_end_ptr - sizeof(uint64_t);
    int non_ascii_in_vec = 0;
    #ifdef __SSE2__
    const char *string_16b_end_ptr = string_end_ptr - sizeof(__m128i);
    __m128i vec_all_chars = _mm_setzero_si128();
    while (cursor < string_16b_end_ptr) {
        __m128i loaded_chars = _mm_loadu_si128((__m128i *)cursor);
        vec_all_chars = _mm_or_si128(loaded_chars, vec_all_chars);
        cursor += sizeof(__m128i);
    }
    non_ascii_in_vec = _mm_movemask_epi8(vec_all_chars);
    #endif

    while (cursor < string_8b_end_ptr) {
        all_chars |= *(uint64_t *)cursor;
        cursor += sizeof(uint64_t);
    }
    while (cursor < string_end_ptr) {
        all_chars |= *cursor;
        cursor += 1;
    }
    return !(non_ascii_in_vec + (all_chars & ASCII_MASK_8BYTE));
}
