%% 
% X = mepack_trlyap(A, Y)
% X = mepack_trlyap(A, Y, Op)
% X = mepack_trlyap(A, Y, OptSet)
% X = mepack_trlyap(A, Y, Op, OptSet)
%     A      -- double precision (quasi) upper triangular coefficient matrix
%     Y      -- double precision right hand side
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
% Xs = mepack_trlyap(As, Ys)
% Xs = mepack_trlyap(As, Ys, Op)
% Xs = mepack_trlyap(As, Ys, OptSet)
% Xs = mepack_trlyap(As, Ys, Op, OptSet)
%     As     -- single precision (quasi) upper triangular coefficient matrix
%     Ys     -- single precision right hand side
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
% The mepack_trlyap function solves the Lyapunov equation
% 
%    op(A) * X + X * op(A)' = Y,
% 
% where A is a (quasi) upper triangular matrix coming from a Schur
% decomposition. The operation op either returns the matrix or
% if Op == 'T' its transpose. In this way, the transposed equation
% can be solved without transposing the coefficient matrix before.
% 
% The matrices A and Y must be of the same size and Y must be symmetric.
% 
% The function uses the level-3 solver D/SLA_TRLYAP_LEVEL3 of MEPACK.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_trlyap:
% 
%    - mb      -- block size of the level-3 algorithm (mb >=2)
%    - isolver -- internal level-2 solver in the algorithm
%        = 0   -- Default/Automatic selection
%        = 1   -- Local Copies, with alignments (if supported by the compiler)
%        = 2   -- Local Copies, without alignments
%        = 3   -- Reordering and Fortran Vectorization
%        = 4   -- Classic level-2 solver
%        = 5   -- Naive level-2 solver
%        = 6   -- Recursive Blocking
%    - openmp  -- OpenMP selection
%        = 1   -- OpenMP - DAG solver enabled and OpenMP accelerated operations
%        = 0   -- Standard level-3 solver without OpenMP acceleration
% 
% 
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 