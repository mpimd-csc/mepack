function test_vector_trait_column()
rv = [ 1 2 3 4];
cv = [ 1; 2; 3; 4];
a  = [ 1 2; 3 4];

assert(strcmp(vector_trait_column(cv), 'colvector') == 1, 'colvector');
fail = 0;
try
	assert(strcmp(vector_trait_column(rv), 'colvector') == 1, 'colvector');
catch
	fail = 1;
end
if ( ~ fail )
	error('fail');
end


fail = 0;
try
	assert(strcmp(vector_trait_column(a), 'colvector') ~= 1, 'unknown');
catch
	fail = 1;
end
if ( ~ fail )
	error('fail');
end


end
