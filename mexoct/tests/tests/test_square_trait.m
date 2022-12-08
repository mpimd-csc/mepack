function test_square_trait()
rv = [ 1 2 3 4];
cv = [ 1; 2; 3; 4];
a  = [ 1 2; 3 4];

assert(strcmp(square_trait(a), 'square') == 1, 'square');
fail = 0;
try
	assert(strcmp(square_trait(rv), 'square') ~= 1, 'unknown');
catch
	fail = 1;
end
if ( ~ fail )
	error('fail1');
end

fail = 0;
try
	assert(strcmp(square_trait(cv), 'square') ~= 1, 'unknown');
catch
	fail = 1;
end
if ( ~ fail )
	error('fail2');
end


end
