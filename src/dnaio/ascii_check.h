#define ASCII_MASK_8BYTE 0x8080808080808080ULL
#define ASCII_MASK_1BYTE 0x80

#include <stddef.h>
#include <stdint.h>

static inline int 
string_is_ascii(char * string, size_t length) {
    size_t n = length;
    uint64_t *longword_ptr = (uint64_t *)string;
    while (n >= sizeof(uint64_t)) {
        if (*longword_ptr & ASCII_MASK_8BYTE){
            return 0;
        }
        longword_ptr += 1;
        n -= sizeof(uint64_t);
    }
    char * char_ptr = (char *)longword_ptr;
    while (n != 0) {
        if (*char_ptr & ASCII_MASK_1BYTE) {
            return 0;
        }
        char_ptr += 1;
        n -= 1;
    }
    return 1;
}
