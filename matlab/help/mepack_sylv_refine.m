%% 
% X = mepack_sylv_refine(A, B, S, T, Q, U, Y)
% X = mepack_sylv_refine(A, B, S, T, Q, U, Y, X0)
% X = mepack_sylv_refine(A, B, S, T, Q, U, Y, OpA, OpB)
% X = mepack_sylv_refine(A, B, S, T, Q, U, Y, X0, OpA, OpB)
% X = mepack_sylv_refine(A, B, S, T, Q, U, Y, OptSet)
% X = mepack_sylv_refine(A, B, S, T, Q, U, Y, X0, OptSet)
% X = mepack_sylv_refine(A, B, S, T, Q, U, Y, Op, OptSet)
% X = mepack_sylv_refine(A, B, S, T, Q, U, Y, X0, OpA, OpB, OptSet)
% 
% Input:
%     A      -- coefficient matrix
%     B      -- coefficient matrix
%     S      -- Schur factor of the coefficient matrix A
%     T      -- Schur factor of the coefficient matrix B
%     Q      -- Unitary transformation of the Schur decomposition of A
%     U      -- Unitary transformation  of the Schur decomposition of B
%     Y      -- right hand side
%     X0     -- Initial Guess
%     OpA    -- Transpose operation on the coefficient matrix A
%     OpB    -- Transpose operation on the coefficient matrix B
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% Return:
%     X      -- Solution of the Lyapunov equation
% 
% [X, convlog] = mepack_sylv_refine(A, B, S, T, Q, U, Y)
% [X, convlog] = mepack_sylv_refine(A, B, S, T, Q, U, Y, X0)
% [X, convlog] = mepack_sylv_refine(A, B, S, T, Q, U, Y, OpA, OpB)
% [X, convlog] = mepack_sylv_refine(A, B, S, T, Q, U, Y, X0, OpA, OpB)
% [X, convlog] = mepack_sylv_refine(A, B, S, T, Q, U, Y, OptSet)
% [X, convlog] = mepack_sylv_refine(A, B, S, T, Q, U, Y, X0, OptSet)
% [X, convlog] = mepack_sylv_refine(A, B, S, T, Q, U, Y, Op, OptSet)
% [X, convlog] = mepack_sylv_refine(A, B, S, T, Q, U, Y, X0, OpA, OpB, OptSet)
% 
% Input:
%     A      -- coefficient matrix
%     B      -- coefficient matrix
%     S      -- Schur factor of the coefficient matrix A
%     T      -- Schur factor of the coefficient matrix B
%     Q      -- Unitary transformation of the Schur decomposition of A
%     U      -- Unitary transformation  of the Schur decomposition of B
%     Y      -- right hand side
%     X0     -- Initial Guess
%     OpA    -- Transpose operation on the coefficient matrix A
%     OpB    -- Transpose operation on the coefficient matrix B
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% Return:
%     X        -- Solution of the Lyapunov equation
%     convlog  -- convergence history
% 
% 
% The mepack_sylv_refine function solves the Sylvester equation
% 
%    opA(A) * X +/- X * opB(B) = Y,
% 
% where A and B are general matrices employing the iterative refinement
% strategy. The operation opA (or opB) either returns the matrix or
% if Op == 'T' its transpose. In this way, the transposed equation
% can be solved without transposing the coefficient matrix before.
% The function requires the matrices  A  and B and their Schur decompositions
% 
%    [ Q, S ] = schur (A, 'real')
%    [ U, T ] = schur (B, 'real')
% 
% in order to work. The Schur decomposition can also be obtain from
% the output of a previous run to mepack_lyap.
% The iterative refinement stops once
% 
%   || opA(A) * X_k +/- X_k * opB(B) - Y ||_F
%   -----------------------------------------  < sqrt(m*n) * (||A||_F+||B||_F) * u * tau
%              || Y ||_F
% 
% is fulfilled. Thereby m is the order of the matrix A, u is the
% machine precision in double or single precision, and tau is a
% security factor. Tau can be changed in the OptSet.
% 
% The function uses the level-3 solver D/SLA_GESYLV_REFINE of MEPACK.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_sylv_refine:
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
% 
% 
%    - maxit   -- maximum number of iterations in the iterative
%                 refinement, default: 15
%    - tau     -- security in the iterative refinement, default: 1.0
%    - openmp  -- enable/disable OpenMP support. By default
%                 the support is enabled GNU Octave and
%                 disabled in MATLAB
% 
% 
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 