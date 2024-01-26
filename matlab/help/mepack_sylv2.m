%% 
% X = mepack_sylv2(A, B, Y)
% X = mepack_sylv2(A, B, Y, OpA, OpB)
% X = mepack_sylv2(A, B, Y, OptSet)
% X = mepack_sylv2(A, B, Y, OpA, OpB, OptSet)
%  Input:
%     A      -- first coefficient matrix
%     B      -- second coefficient matrix
%     Y      -- right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
%  Return:
%     X      -- Solution of the Sylvester equation
% 
% 
% [X, S, T, Q, Z] = mepack_sylv2(A, B, Y)
% [X, S, T, Q, Z] = mepack_sylv2(A, B, Y, OpA, OpB)
% [X, S, T, Q, Z] = mepack_sylv2(A, B, Y, OptSet)
% [X, S, T, Q, Z] = mepack_sylv2(A, B, Y, OpA, OpB, OptSet)
%   Input:
%     A      -- first coefficient matrix
%     B      -- second coefficient matrix
%     Y      -- right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     X      -- Solution of the Sylvester  equation
%     S, Q   -- Schur decomposition of A, A = Q * S * Q'
%     T, Z   -- Schur decomposition of B, B = Z * T * B'
% 
% X = mepack_sylv2(S, T, Q, Z, Y)
% X = mepack_sylv2(S, T, Q, Z, Y, OpA, OpB)
% X = mepack_sylv2(S, T, Q, Z, Y, OptSet)
% X = mepack_sylv2(S, T, Q, Z, Y, OpA, OpB, OptSet)
%   Input:
%     S      -- (quasi) upper triangular coefficient matrix
%     T      -- (quasi) upper triangular coefficient matrix
%     Q      -- unitary transformation matrix, A = Q*S*Q'
%     Z      -- unitary transformation matrix, B = Z*T*Z'
%     Y      -- right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     X      -- Solution of the Sylvester equation
% 
% The mepack_sylv2 function solves the Sylvester equation
% 
%    opA(A) * X * opB(B) + X = Y,
% 
% or
% 
%    opA(A) * X * opB(B) - X = Y,
% 
% where A and B are general matrices. The operation opA (or opB) either returns the matrix or
% if OpA == 'T' (or opB == 'T') its transpose. In this way, the transposed equation
% can be solved without transposing the coefficient matrix before. The function
% can also be called with A and B decomposed into their Schur decompositions,
% A = Q*S*Q' and B = Q*T*Q', where S is the Schur form of A and T is the Schur form
% of B. If the function is called with five return values,
% the Schur forms of A and B and their unitary transformation matrices are returned.
% 
% The function uses the level-3 solver D/SLA_GESYLV2 of MEPACK.
% 
% The involved Schur decomposition is aware of matrices in Hessenberg form.
% In this case the initial Hessenberg reduction of the Schur decomposition is skipped.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_sylv2
% 
%    - mb      -- row block size of the level-3 algorithm (mb >=2)
%    - bignb   -- block size of two-stage solver (bignb >=2)
%    - nb      -- column block size of the level-3 algorithm (nb >=2)
%    - isolver -- internal level-2 solver in the algorithm
%        = 0   -- Default/Automatic selection
%        = 1   -- Local Copies, with alignments (if supported by the compiler)
%        = 2   -- Local Copies, without alignments
%        = 3   -- Reordering and Fortran Vectorization
%        = 4   -- Classic level-2 solver
%        = 5   -- Naive level-2 solver
%        = 6   -- Recursive Blocking
%    - sign    -- swap the sign between both terms
%        = 1   -- solve the equation with '+'
%        = -1  -- solve the equation with '-'
%    - fsolver -- determine the internal level 3 solver
%        = 0   -- automatic or default selection
%        = 1   -- standard level 3 solver
%        = 2   -- standard level 2 solver
%        = 3   -- DAG OpenMP 4 solver
%        = 4   -- Two-Stage level 3 solver
%        = 5   -- Recursive Blocking Solver
%    - openmp  -- enable/disable OpenMP support. By default
%                 the support is enabled GNU Octave and
%                 disabled in MATLAB
% 
% 
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 