function test_vector_trait()
rv = [ 1 2 3 4];
cv = [ 1; 2; 3; 4];
a  = [ 1 2; 3 4];

assert(strcmp(vector_trait_row(rv), 'rowvector') == 1, 'rowvector');
fail = 0;
try
	assert(strcmp(vector_trait_row(cv), 'rowvector') == 1, 'colvector');
catch
	fail = 1;
end
if ( ~ fail )
	error('fail');
end


fail = 0;
try
	assert(strcmp(vector_trait_row(a), 'rowvector') ~= 1, 'unknown');
catch
	fail = 1;
end
if ( ~ fail )
	error('fail');
end


end
