%% 
% X = mepack_gsylv(A, B, C, D, Y)
% X = mepack_gsylv(A, B, C, D, Y, OpA, OpB)
% X = mepack_gsylv(A, B, C, D, Y, OptSet)
% X = mepack_gsylv(A, B, C, D, Y, OpA, OpB, OptSet)
%  Input:
%     A      -- first coefficient matrix
%     B      -- second coefficient matrix
%     C      -- third coefficient matrix
%     D      -- fourth coefficient matrix
%     Y      -- right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
%  Return:
%     X      -- Solution of the Sylvester equation
% 
% 
% [X, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_gsylv(A, B, C, D, Y)
% [X, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_gsylv(A, B, C, D, Y, OpA, OpB)
% [X, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_gsylv(A, B, C, D, Y, OptSet)
% [X, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_gsylv(A, B, C, D, Y, OpA, OpB, OptSet)
%   Input:
%     A      -- first coefficient matrix
%     B      -- second coefficient matrix
%     C      -- third coefficient matrix
%     D      -- fourth coefficient matrix
%     Y      -- right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     X      -- Solution of the Sylvester  equation
%     AS, CS -- generalized Schur decomposition of (A, C)
%     BS, DS -- generalized Schur decomposition of (B, D)
%     QA, ZA -- transformation matrices for (A,C) = (QA*AS*ZA',QA*CS*ZA')
%     QB, ZB -- transformation matrices for (B,D) = (QB*BS*ZB',QB*DS*ZB')
% 
% 
% X = mepack_gsylv(AS, BS, CS, DS, QA, ZA, QB, ZB, Y)
% X = mepack_gsylv(AS, BS, CS, DS, QA, ZA, QB, ZB, Y, OpA, OpB)
% X = mepack_gsylv(AS, BS, CS, DS, QA, ZA, QB, ZB, Y, OptSet)
% X = mepack_gsylv(AS, BS, CS, DS, QA, ZA, QB, ZB, Y, OpA, OpB, OptSet)
%   Input:
%     AS, CS -- Schur form of (A,C)
%     BS, DS -- Schur form of (B,D)
%     QA, ZA -- transformation matrices for (A,C) = (QA*AS*ZA',QA*CS*ZA')
%     QB, ZB -- transformation matrices for (B,D) = (QB*BS*ZB',QB*DS*ZB')
%     Y      -- right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
%   Return:
%     X      -- Solution of the Sylvester equation
% 
% The mepack_gsylv function solves the Sylvester equation
% 
%    opA(A) * X * opB(B) + opA(C) * X * opB(D) = Y,
% 
% or
% 
%    opA(A) * X * opB(B) - opA(C) * X * opB(D) = Y,
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
% The function uses the level-3 solver D/SLA_GGSYLV of MEPACK.
% 
% The involved generalized Schur decomposition is aware of matrix pairs in
% Hessenberg-Triangular form. In this case the initial Hessenberg reduction
% of the generalized Schur decomposition is skipped.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_gsylv
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
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 