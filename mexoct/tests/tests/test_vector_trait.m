function test_vector_trait()
rv = [ 1 2 3 4];
cv = [ 1; 2; 3; 4];
a  = [ 1 2; 3 4];

assert(strcmp(vector_trait(rv), 'vector') == 1, 'vector');
assert(strcmp(vector_trait(cv), 'vector') == 1, 'vector');
fail = 0;
try
	assert(strcmp(vector_trait(a), 'vector') ~= 1, 'unknown');
catch
	fail = 1;
end
if ( ~ fail )
	error('fail');
end


end
