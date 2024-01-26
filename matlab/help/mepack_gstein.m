%% 
% X = mepack_Stein(A, B, Y)
% X = mepack_Stein(A, B, Y, Op)
% X = mepack_gstein(A, B, Y, OptSet)
% X = mepack_gstein(A, B, Y, Op, OptSet)
%  Input:
%     A      -- first coefficient matrix
%     B      -- second coefficient matrix
%     Y      -- right hand side
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
%  Return:
%     X      -- Solution of the generalized Stein equation
% 
% 
% [X, S, T, Q, Z] = mepack_gstein(A, B, Y)
% [X, S, T, Q, Z] = mepack_gstein(A, B, Y, Op)
% [X, S, T, Q, Z] = mepack_gstein(A, B, Y, OptSet)
% [X, S, T, Q, Z] = mepack_gstein(A, B, Y, Op, OptSet)
%   Input:
%     A      -- first coefficient matrix
%     B      -- second coefficient matrix
%     Y      -- right hand side
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     X      -- Solution of the generalized Stein equation
%     S, T,
%     Q, Z   -- generalized Schur decomposition of (A,B),
%               A = Q*S*Z', B = Q*T*Z'
% 
% X = mepack_gstein(S, T, Q, Z, Y)
% X = mepack_gstein(S, T, Q, Z, Y, Op)
% X = mepack_gstein(S, T, Q, Z, Y, OptSet)
% X = mepack_gstein(S, T, Q, Z, Y, Op, OptSet)
%   Input:
%     S      -- (quasi) upper triangular coefficient matrix
%     T      -- upper triangular coefficient matrix
%     Q      -- unitary transformation matrix, A = Q*S*Z'
%     Z      -- unitary transformation matrix, B = Q*T*Z'
%     Y      -- right hand side
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     X      -- Solution of the generalized Stein equation
% 
% The mepack_gstein function solves the generalized Lyapunov equation
% 
%    op(A) * X * op(A)' - op(B) * X * op(A)' = Y,
% 
% where A and B are general matrices. The operation op either returns the matrix or
% if Op == 'T' its transpose. In this way, the transposed equation
% can be solved without transposing the coefficient matrix before. The function
% can also be called with (A,B) decomposed in to A = Q*S*Z' and B = Q*T*Z',
% where (S,T) is the generalized Schur form of (A, B).
% If the function is called with five return values,
% the generalized Schur form of (A,B) and its unitary transformation matrices are returned.
% 
% The matrices A, B, S, T, Q, Z  and Y must be of the same size and Y must be symmetric.
% 
% The function uses the level-3 solver D/SLA_GGSTEIN of MEPACK.
% 
% The involved generalized Schur decomposition is aware of matrix pairs in
% Hessenberg-Triangular form. In this case the initial Hessenberg reduction
% of the generalized Schur decomposition is skipped.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_gstein:
% 
%    - mb      -- block size of the level-3 algorithm (mb >=2)
%    - bignb   -- block size of two-stage solver (bignb >=2)
%    - isolver -- internal level-2 solver in the algorithm
%        = 0   -- Default/Automatic selection
%        = 1   -- Local Copies, with alignments (if supported by the compiler)
%        = 2   -- Local Copies, without alignments
%        = 3   -- Reordering and Fortran Vectorization
%        = 4   -- Classic level-2 solver
%        = 5   -- Naive level-2 solver
%        = 6   -- Recursive Blocking
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