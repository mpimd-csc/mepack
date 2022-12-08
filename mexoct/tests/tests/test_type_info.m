function test_type_info()

assert(strcmp(type_info(int8(1)) , 'int8_t') == 1, 'int8_t');
assert(strcmp(type_info(uint8(1)) , 'uint8_t') == 1, 'uint8_t');
assert(strcmp(type_info(int16(1)) , 'int16_t') == 1, 'int16_t');
assert(strcmp(type_info(uint16(1)) , 'uint16_t') == 1, 'uint16_t');
assert(strcmp(type_info(int32(1)) , 'int32_t') == 1, 'int32_t');
assert(strcmp(type_info(uint32(1)) , 'uint32_t') == 1, 'uint32_t');
assert(strcmp(type_info(int64(1)) , 'int64_t') == 1, 'int64_t');
assert(strcmp(type_info(uint64(1)) , 'uint64_t') == 1, 'uint64_t');

assert(strcmp(type_info(1.1) , 'double') == 1, 'double');
assert(strcmp(type_info(single(1)) , 'float') == 1, 'float');

assert(strcmp(type_info(1.1+sqrt(-1)) , 'double complex') == 1, 'double complex');
assert(strcmp(type_info(single(1+sqrt(-1))) , 'float complex') == 1, 'float complex');

end
