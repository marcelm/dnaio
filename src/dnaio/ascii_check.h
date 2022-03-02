#define ASCII_MASK_8BYTE 0x8080808080808080ULL
#define ASCII_MASK_1BYTE 0x80

#include <stddef.h>
#include <stdint.h>

static int
string_is_ascii(char * string, size_t length) {
    size_t n = length;
    char * char_ptr = string;
    // The first loop aligns the memory address. Char_ptr is cast to a size_t
    // to return the memory address. Uint64_t is 8 bytes long, and the processor
    // handles this better when its address is a multiplier of 8. This loops
    // handles the first few bytes that are not on such a multiplier boundary.
    while ((size_t)char_ptr % sizeof(uint64_t) && n != 0) {
        if (*char_ptr & ASCII_MASK_1BYTE) {
            return 0;
        }
        char_ptr += 1;
        n -= 1;
    }
    uint64_t *longword_ptr = (uint64_t *)char_ptr;
    while (n >= sizeof(uint64_t)) {
        if (*longword_ptr & ASCII_MASK_8BYTE){
            return 0;
        }
        longword_ptr += 1;
        n -= sizeof(uint64_t);
    }
    char_ptr = (char *)longword_ptr;
    while (n != 0) {
        if (*char_ptr & ASCII_MASK_1BYTE) {
            return 0;
        }
        char_ptr += 1;
        n -= 1;
    }
    return 1;
}
