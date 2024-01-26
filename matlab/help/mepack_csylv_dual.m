%% 
% [R, L] = mepack_csylv_dual(A, B, C, D, E, F)
% [R, L] = mepack_csylv_dual(A, B, C, D, E, F, OpA, OpB)
% [R, L] = mepack_csylv_dual(A, B, C, D, E, F, OptSet)
% [R, L] = mepack_csylv_dual(A, B, C, D, E, F, OpA, OpB, OptSet)
%  Input:
%     A      -- first coefficient matrix
%     B      -- second coefficient matrix
%     C      -- third coefficient matrix
%     D      -- fourth coefficient matrix
%     E      -- first right hand side
%     F      -- second right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
%  Return:
%     R      -- first solution of the Sylvester equation
%     L      -- second solution of the Sylvester equation
% 
% [R, L, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E, F)
% [R, L, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E, F, OpA, OpB)
% [R, L, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E, F, OptSet)
% [R, L, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E, F, OpA, OpB, OptSet)
%   Input:
%     A      -- first coefficient matrix
%     B      -- second coefficient matrix
%     C      -- third coefficient matrix
%     D      -- fourth coefficient matrix
%     E      -- first right hand side
%     F      -- second right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     R      -- first solution of the Sylvester equation
%     L      -- second solution of the Sylvester equation
%     AS, CS -- generalized Schur decomposition of (A, C)
%     BS, DS -- generalized Schur decomposition of (B, D)
%     QA, ZA -- transformation matrices for (A,C) = (QA*AS*ZA',QA*CS*ZA')
%     QB, ZB -- transformation matrices for (B,D) = (QB*BS*ZB',QB*DS*ZB')
% 
% 
% [R, L] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E, F)
% [R, L] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E, F, OpA, OpB)
% [R, L] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E, F, OptSet)
% [R, L] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E, F, OpA, OpB, OptSet)
%   Input:
%     AS, CS -- Schur form of (A,C)
%     BS, DS -- Schur form of (B,D)
%     QA, ZA -- transformation matrices for (A,C) = (QA*AS*ZA',QA*CS*ZA')
%     QB, ZB -- transformation matrices for (B,D) = (QB*BS*ZB',QB*DS*ZB')
%     E      -- first right hand side
%     F      -- second right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
%   Return:
%     R      -- first solution of the Sylvester equation
%     L      -- second solution of the Sylvester equation
% 
% The mepack_csylv_dual function solves the coupled Sylvester equation
% 
%    opA(A)**T * R + opA(C)**T * L = E,
%    +/- R * opB(B)**T +/- L * opB(D)**T = F,
% 
% where A, B, C, D are general matrices. The operation opA (or opB) either returns the matrix or
% if OpA == 'T' (or opB == 'T') its transpose. In this way, the transposed equation
% can be solved without transposing the coefficient matrix before. The function
% can also be called with (A, C) and (B, D)  decomposed into their generalized
% Schur decompositions, (A,C) = (QA*AS*ZA',QA*CS*ZA') and (QB*BS*ZB',QB*DS*ZB')
% If the function is called with nine return values,
% the generalized Schur forms of (A, C) and (B,D) including their
% unitary transformation matrices are returned.
% 
% The function uses the level-3 solver D/SLA_GGCSYLV_DUAL of MEPACK.
% 
% The involved generalized Schur decomposition is aware of matrix pairs in
% Hessenberg-Triangular form. In this case the initial Hessenberg reduction
% of the generalized Schur decomposition is skipped.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_csylv_dual
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
%    - sign    -- swap the sign in the first equation
%        = 1   -- solve the equation with '+'
%        = -1  -- solve the equation with '-'
%    - sign2   -- swap the sign in the second equation
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
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 