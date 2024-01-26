%% 
% X = mepack_lyap(A, Y)
% X = mepack_lyap(A, Y, Op)
% X = mepack_lyap(A, Y, OptSet)
% X = mepack_lyap(A, Y, Op, OptSet)
%  Input:
%     A      -- coefficient matrix
%     Y      -- right hand side
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
%  Return:
%     X      -- Solution of the Lyapunov equation
% 
% 
% [X, S, Q] = mepack_lyap(A, Y)
% [X, S, Q] = mepack_lyap(A, Y, Op)
% [X, S, Q] = mepack_lyap(A, Y, OptSet)
% [X, S, Q] = mepack_lyap(A, Y, Op, OptSet)
%   Input:
%     A      -- coefficient matrix
%     Y      -- right hand side
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     X      -- Solution of the Lyapunov equation
%     S, Q   -- Schur decomposition of A, A = Q*S*Q'
% 
% X = mepack_lyap(S, Q, Y)
% X = mepack_lyap(S, Q, Y, Op)
% X = mepack_lyap(S, Q, Y, OptSet)
% X = mepack_lyap(S, Q, Y, Op, OptSet)
%   Input:
%     S      -- (quasi) upper triangular coefficient matrix, Schur form of A
%     Q      -- unitary transformation matrix, A = Q*S*Q'
%     Y      -- right hand side
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     X      -- Solution of the Lyapunov equation
% 
% The mepack_lyap function solves the Lyapunov equation
% 
%    op(A) * X + X * op(A)' = Y,
% 
% where A is a general matrix. The operation op either returns the matrix or
% if Op == 'T' its transpose. In this way, the transposed equation
% can be solved without transposing the coefficient matrix before. The function
% can also be called with A decomposed in to A = Q*S*Q', where S is the
% Schur form of A. If the function is called with three return values,
% the Schur form of A and its unitary transformation matrix is returned.
% 
% The matrices A and Y must be of the same size and Y must be symmetric.
% 
% The function uses the level-3 solver D/SLA_GELYAP of MEPACK.
% 
% The involved Schur decomposition is aware of matrices in Hessenberg form.
% In this case the initial Hessenberg reduction of the Schur decomposition is skipped.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_lyap:
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
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 