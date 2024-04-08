#include <stddef.h>

#define ASCII_MASK_8BYTE 0x8080808080808080ULL
#define ASCII_MASK_1BYTE 0x80

static inline int string_is_ascii_fallback(const char *string, size_t length)
{
    /* Combining all characters with OR allows for only one bit check at the end */
    size_t all_chars = 0;
    for (size_t i=0; i<length; i++) {
        all_chars |= string[i];
    }
    return !(all_chars & ASCII_MASK_1BYTE);
}

/**
 * @brief Check if a string of given length only contains ASCII characters.
 *
 * @param string A char pointer to the start of the string.
 * @param length The length of the string. This funtion does not check for
 *               terminating NULL bytes.
 * @returns 1 if the string is ASCII-only, 0 otherwise.
 */
static int
string_is_ascii(const char *string, size_t length)
{
    if (length < sizeof(size_t)) {
        return string_is_ascii_fallback(string, length);
    }
    size_t number_of_chunks = length / sizeof(size_t);
    size_t *chunks = (size_t *)string;
    size_t number_of_unrolls = number_of_chunks / 4;
    size_t remaining_chunks = number_of_chunks - (number_of_unrolls * 4);
    size_t *chunk_ptr = chunks;
    size_t all_chars0 = 0;
    size_t all_chars1 = 0;
    size_t all_chars2 = 0;
    size_t all_chars3 = 0;
    for (size_t i=0; i < number_of_unrolls; i++) {
        /* Performing indepedent OR calculations allows the compiler to use
           vectors. It also allows out of order execution. */
        all_chars0 |= chunk_ptr[0];
        all_chars1 |= chunk_ptr[1];
        all_chars2 |= chunk_ptr[2];
        all_chars3 |= chunk_ptr[3];
        chunk_ptr += 4;
    }
    size_t all_chars = all_chars0 | all_chars1 | all_chars2 | all_chars3;
    for (size_t i=0; i<remaining_chunks; i++) {
        all_chars |= chunk_ptr[i];
    }
    /* Load the last few bytes left in a single integer for fast operations.
       There is some overlap here with the work done before, but for a simple
       ascii check this does not matter. */
    size_t last_chunk = *(size_t *)(string + length - sizeof(size_t));
    all_chars |= last_chunk;
    return !(all_chars & ASCII_MASK_8BYTE);
}
