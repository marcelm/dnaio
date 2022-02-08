#define ASCII_MASK_8BYTE 0x8080808080808080ULL
#define ASCII_MASK_1BYTE 0x80

#include <stddef.h>
#include <stdint.h>

static inline int 
string_is_ascii(char * string, ssize_t length) {
    uint64_t *longword_ptr = (uint64_t *)string;
    while (((char *)longword_ptr - string) < length) {
        if (*longword_ptr & ASCII_MASK_8BYTE){
            return 0;
        }
        longword_ptr += 1;
    }
    char * char_ptr = (char *)longword_ptr;
    while (char_ptr - string < length) {
        if (*char_ptr & ASCII_MASK_1BYTE) {
            return 0;
        }
    }
    return 1;
}
