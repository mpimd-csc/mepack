%% 
% X = mepack_tgsylv(A, B, C, D, Y)
% X = mepack_tgsylv(A, B, C, D, Y, OpA, OpB)
% X = mepack_tgsylv(A, B, C, D, Y, OptSet)
% X = mepack_tgsylv(A, B, C, D, Y, OpA, OpB, OptSet)
%     A      -- double precision (quasi) upper triangular coefficient matrix
%     B      -- double precision (quasi) upper triangular coefficient matrix
%     C      -- double precision upper triangular coefficient matrix
%     D      -- double precision upper triangular coefficient matrix,
%               can be quasi upper triangular if B is upper triangular
%     Y      -- double precision right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
% Xs = mepack_tgsylv(As, Bs, Cs, Ds, Ys)
% Xs = mepack_tgsylv(As, Bs, Cs, Ds, Ys, Op)
% Xs = mepack_tgsylv(As, Bs, Cs, Ds, Ys, OptSet)
% Xs = mepack_tgsylv(As, Bs, Cs, Ds, Ys, Op, OptSet)
%     As     -- single precision (quasi) upper triangular coefficient matrix
%     Bs     -- single precision (quasi) upper triangular coefficient matrix
%     Cs     -- single precision upper triangular coefficient matrix
%     Ds     -- single precision upper triangular coefficient matrix,
%               can be quasi upper triangular if Bs is upper triangular
%     Ys     -- single precision right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
% The mepack_trsylv function solves the Sylvester equation
% 
%    opA(A) * X * opB(B) + opA(C) * X * opB(D) = Y,
% 
% or
% 
%    opA(A) * X * op(B) - opA(C) * X * opB(D)  = Y,
% 
% where (A, C) and (B,D) (or (D,B)) are (quasi) upper triangular matrix pairs,
% both coming from a generalized Schur decomposition. The operation op
% either returns the matrix or if Op == 'T' its transpose. In this way, the
% transposed equation can be solved without transposing the coefficient
% matrix before.
% 
% The function uses the level-3 solver D/SLA_TGSYLV_* of MEPACK.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_tgsylv:
% 
%    - mb      -- row block size of the level-3 algorithm (mb >=2)
%    - nb      -- column block size of the level-3 algorithm (nb >=2)
%    - isolver -- internal level-2 solver in the algorithm
%        = 0   -- Default/Automatic selection
%        = 1   -- Local Copies, with alignments (if supported by the compiler)
%        = 2   -- Local Copies, without alignments
%        = 3   -- Reordering and Fortran Vectorization
%        = 4   -- Classic level-2 solver
%        = 5   -- Naive level-2 solver
%        = 6   -- Recursive Blocking
%    - openmp  -- OpenMP selection
%        = 1   -- OpenMP - DAG solver enabled
%        = 0   -- Standard level-3 solver
%    - sign    -- swap the sign between both terms
%        = 1   -- solve the equation with '+'
%        = -1  -- solve the equation with '-'
% 
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 