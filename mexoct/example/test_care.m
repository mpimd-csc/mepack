
A = rand(8);
B = rand(8);
C = rand(8);
E = eye(8);

% valid calling variants 
X = example_care(A, B, C);
X = example_care(A, B, C, E);
X = example_care(sparse(A), sparse(B), sparse(C), 'foo');
X = example_care(sparse(A), sparse(B), sparse(C), 'bar');

% errors
fprintf('\n\ntest: traits invalid \n');
try X = example_care(A, sparse(B), C); catch end
fprintf('\n\ntest: wrong type \n');
try X = example_care(sparse(A), sparse(B), sparse(C), 1); catch end
fprintf('\n\ntest: to few return values \n');
try example_care(A, B, C); catch end