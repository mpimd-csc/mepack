%% 
% X = mepack_trsylv2(A, B, Y)
% X = mepack_trsylv2(A, B, Y, OpA, OpB)
% X = mepack_trsylv2(A, B, Y, OptSet)
% X = mepack_trsylv2(A, B, Y, OpA, OpB, OptSet)
%     A      -- double precision (quasi) upper triangular coefficient matrix
%     B      -- double precision (quasi) upper triangular coefficient matrix
%     Y      -- double precision right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
% Xs = mepack_trsylv2(As, Bs, Ys)
% Xs = mepack_trsylv2(As, Bs, Ys, Op)
% Xs = mepack_trsylv2(As, Bs, Ys, OptSet)
% Xs = mepack_trsylv2(As, Bs, Ys, Op, OptSet)
%     As     -- single precision (quasi) upper triangular coefficient matrix
%     BS     -- single precision (quasi) upper triangular coefficient matrix
%     Ys     -- single precision right hand side
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
% The mepack_trsylv function solves the Sylvester equation
% 
%    opA(A) * X * opB(B) + X = Y,
% 
% or
% 
%    opA(A) * X * op(B) - X  = Y,
% 
% where A and B are (quasi) upper triangular matrices, both coming from
% a generalized Schur decomposition. The operation op
% either returns the matrix or if Op == 'T' its transpose. In this way, the
% transposed equation can be solved without transposing the coefficient
% matrix before.
% 
% The function uses the level-3 solver D/SLA_TRSYLV2_* of MEPACK.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_trsylv2:
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
%    - sign    -- swap the sign between both terms
%        = 1   -- solve the equation with '+'
%        = -1  -- solve the equation with '-'
%    - openmp  -- OpenMP selection
%        = 1   -- OpenMP - DAG solver enabled and OpenMP accelerated operations
%        = 0   -- Standard level-3 solver without OpenMP acceleration
% 
% 
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 