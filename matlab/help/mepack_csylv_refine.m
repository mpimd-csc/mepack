%% 
% [R, L] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F)
% [R, L] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0)
% [R, L] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, OpA, OpB)
% [R, L] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OpA, OpB)
% [R, L] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, OptSet)
% [R, L] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OptSet)
% [R, L] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, Op, OptSet)
% [R, L] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OpA, OpB, OptSet)
% [R, L, convlog] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F)
% [R, L, convlog] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0)
% [R, L, convlog] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, OpA, OpB)
% [R, L, convlog] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OpA, OpB)
% [R, L, convlog] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, OptSet)
% [R, L, convlog] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OptSet)
% [R, L, convlog] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, Op, OptSet)
% [R, L, convlog] = mepack_csylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OpA, OpB, OptSet)
% 
% 
% Input:
%     A      -- coefficient matrix
%     B      -- coefficient matrix
%     C      -- coefficient matrix
%     B      -- coefficient matrix
%     AS     -- Schur factor of the matrix pair (A, C)
%     BS     -- Schur factor of the matrix pair (B, D)
%     CS     -- Schur factor of the matrix pair (A, C)
%     DS     -- Schur factor of the matrix pair (B, D)
%     Q      -- Unitary transformation of the Schur decomposition of (A,C)
%     Z      -- Unitary transformation of the Schur decomposition of (A,C)
%     U      -- Unitary transformation  of the Schur decomposition of (B,D)
%     V      -- Unitary transformation  of the Schur decomposition of (B,D)
%     E      -- first right hand side
%     F      -- second right hand side
%     R0     -- Initial Guess for the first solution
%     L0     -- Initial Guess for the second solution
%     OpA    -- Transpose operation on the coefficient matrix A
%     OpB    -- Transpose operation on the coefficient matrix B
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% Return:
%     R       -- first solution of the coupled Sylvester equation
%     L       -- second solution of the coupled Sylvester equation
%     convlog -- convergence history
% 
% The mepack_csylv_refine function solves the coupled generalized Sylvester equation
% 
%     opA(A) * R  +/- L  * opB(B) =  E
%     opA(C) * R  +/- L  * opB(D) =  F
% 
% where A, B, C, and D are general matrices employing the iterative refinement
% strategy. The operation opA (or opB) either returns the matrix or
% if Op == 'T' its transpose. In this way, the transposed equation
% can be solved without transposing the coefficient matrix before.
% The function requires the matrix pairs (A,C) and (B,D) and their Schur decompositions
% 
%     if (exist('OCTAVE_VERSION', 'builtin'))
%         [AS, CS, Q, Z] = qz(A, C);
%         [BS, DS, U, V] = qz(B, D);
%     else
%         [AS, CS, Q, Z] = qz(A, C, 'real');
%         [BS, DS, U, V] = qz(B, D, 'real');
%     end
%     Q = Q';
%     U = U';
% 
% in order to work. The Schur decomposition can also be obtain from
% the output of a previous run to mepack_csylv.
% The iterative refinement stops once
% 
%   || opA(A) * R_k * +/- L_k * opB(D) - E ||_F
%   ------------------------------------------- < tol
%                  || E ||_F
% 
% and
% 
%   || opA(A) * R_k * +/- L_k * opB(D) - F||_F
%   ------------------------------------------ < tol
%                  || F ||_F
% 
% with
% 
%   tol = sqrt ( m * n) * (||A||_F + ||B||_F + ||C||_F + ||D||_F) * u * tau
% 
% are fulfilled. Thereby m is the order of the matrix A, n the order
% of the matrix B, u is the machine precision in double or single precision,
% and tau is a security factor. Tau can be changed in the OptSet.
% 
% The function uses the level-3 solver D/SLA_GGCSYLV_REFINE of MEPACK.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_csylv_refine:
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
%    - sign    -- sign in the first of the two equations
%        = 1   -- solve the equation with '+'
%        = -1  -- solve the equation with '-'
%    - sign2   -- sign in the second of the two equations
%        = 1   -- solve the equation with '+'
%        = -1  -- solve the equation with '-'
%    - fsolver -- determine the internal level 3 solver
%        = 0   -- automatic or default selection
%        = 1   -- standard level 3 solver
%        = 2   -- standard level 2 solver
%        = 3   -- DAG OpenMP 4 solver
%        = 4   -- Two-Stage level 3 solver
%        = 5   -- Recursive Blocking Solver
%    - maxit   -- maximum number of iterations in the iterative
%                 refinement, default: 15
%    - tau     -- security in the iterative refinement, default: 1.0
%    - openmp  -- enable/disable OpenMP support. By default
%                 the support is enabled GNU Octave and
%                 disabled in MATLAB
% 
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 