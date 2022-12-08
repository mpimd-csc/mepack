function test_copy_input()

x = int8(2); assert(copy_input(x) == x, 'int8_t');
x = uint8(2); assert(copy_input(x) == x, 'uint8_t');
x = int16(2); assert(copy_input(x) == x, 'int16_t');
x = uint16(2); assert(copy_input(x) == x, 'uint16_t');
x = int32(2); assert(copy_input(x) == x, 'int32_t');
x = uint32(2); assert(copy_input(x) == x, 'uint32_t');
x = int64(2); assert(copy_input(x) == x, 'int64_t');
x = uint64(2); assert(copy_input(x) == x, 'uint64_t');

x = 2.0; assert(copy_input(x) == x, 'double');
x = single(2.0) ; assert(copy_input(x) == x, 'single');

x = 2.0+sqrt(-1); assert(copy_input(x) == x, 'double complex');
x = single(2.0+sqrt(-1)) ; assert(copy_input(x) == x, 'single complex');

A = randi(20,20,20);
assert(norm(double(copy_input(int8(A)))-A,1) == 0, 'int8_t');
assert(norm(double(copy_input(uint8(A)))-A,1) == 0, 'uint8_t');
assert(norm(double(copy_input(int16(A)))-A,1) == 0, 'int16_t');
assert(norm(double(copy_input(uint16(A)))-A,1) == 0, 'uint16_t');
assert(norm(double(copy_input(int32(A)))-A,1) == 0, 'int32_t');
assert(norm(double(copy_input(uint32(A)))-A,1) == 0, 'uint32_t');
assert(norm(double(copy_input(int64(A)))-A,1) == 0, 'int64_t');
assert(norm(double(copy_input(uint64(A)))-A,1) == 0, 'uint64_t');

assert(norm(copy_input(A)-A,1) == 0, 'double');
assert(norm(double(copy_input(single(A))-A),1) == 0,'single');

A = randi(20,20,20) + sqrt(-1) * randi(20,20,20);
assert(norm(copy_input(A)-A,1) == 0, 'double');
assert(norm(double(copy_input(single(A))-A),1) == 0,'single');


end
