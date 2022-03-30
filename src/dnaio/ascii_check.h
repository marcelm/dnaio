#define ASCII_MASK_8BYTE 0x8080808080808080ULL
#define ASCII_MASK_1BYTE 0x80

#include <stddef.h>
#include <stdint.h>

static int
string_is_ascii(char * string, size_t length) {
    size_t n = length;
    uint64_t all_chars = 0;
    char * char_ptr = string;
    // The first loop aligns the memory address. Char_ptr is cast to a size_t
    // to return the memory address. Uint64_t is 8 bytes long, and the processor
    // handles this better when its address is a multiplier of 8. This loops
    // handles the first few bytes that are not on such a multiplier boundary.
    while ((size_t)char_ptr % sizeof(uint64_t) && n != 0) {
        all_chars |= *char_ptr;
        char_ptr += 1;
        n -= 1;
    }
    uint64_t *longword_ptr = (uint64_t *)char_ptr;
    while (n >= sizeof(uint64_t)) {
        all_chars |= *longword_ptr;
        longword_ptr += 1;
        n -= sizeof(uint64_t);
    }
    char_ptr = (char *)longword_ptr;
    while (n != 0) {
        all_chars |= *char_ptr;
        char_ptr += 1;
        n -= 1;
    }
    return !(all_chars & ASCII_MASK_8BYTE);
}
