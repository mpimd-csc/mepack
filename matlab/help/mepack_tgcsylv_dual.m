%% 
% [R, L] = mepack_tgcsylv_dual(A, B, C, D, E, F)
% [R, L] = mepack_tgcsylv_dual(A, B, C, D, E, F, OpA, OpB)
% [R, L] = mepack_tgcsylv_dual(A, B, C, D, E, F, OptSet)
% [R, L] = mepack_tgcsylv_dual(A, B, C, D, E, F, OpA, OpB, OptSet)
%     A      -- double precision (quasi) upper triangular coefficient matrix
%     B      -- double precision (quasi) upper triangular coefficient matrix
%     C      -- double precision upper triangular coefficient matrix
%     D      -- double precision upper triangular coefficient matrix,
%               can be quasi upper triangular if B is upper triangular
%     E, F   -- double precision right hand sides
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
% [Rs, Ls] = mepack_tgcsylv_dual(As, Bs, Cs, Ds, Es, Fs)
% [Rs, Ls] = mepack_tgcsylv_dual(As, Bs, Cs, Ds, Es, Fs, Op)
% [Rs, Ls] = mepack_tgcsylv_dual(As, Bs, Cs, Ds, Es, Fs, OptSet)
% [Rs, Ls] = mepack_tgcsylv_dual(As, Bs, Cs, Ds, Es, Fs, Op, OptSet)
%     As     -- single precision (quasi) upper triangular coefficient matrix
%     Bs     -- single precision (quasi) upper triangular coefficient matrix
%     Cs     -- single precision upper triangular coefficient matrix
%     Ds     -- single precision upper triangular coefficient matrix,
%               can be quasi upper triangular if Bs is upper triangular
%     Es, Fs -- single precision right hand sides
%     OpA    -- Transpose operation on the left coefficient matrix
%     OpB    -- Transpose operation on the right coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
% 
% The mepack_tgcsylv_dual function solves the dual coupled generalized Sylvester equation
% 
%     opA(A)' * R            + opA(C)' * R     = E
%           +/- R * opB(B)'  +/- L * opB(D)'   = F
% 
% where (A, C) and (B,D) (or (D,B)) are (quasi) upper triangular matrix pairs,
% both coming from a generalized Schur decomposition. The operation op
% either returns the matrix or if Op == 'T' its transpose. In this way, the
% transposed equation can be solved without transposing the coefficient
% matrix before.
% 
% The function uses the level-3 solver D/SLA_TGCSYLV_DUAL_* of MEPACK.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_tgcsylv_dual:
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
%    - sign    -- swap the sign between both terms in the first half of the second equation
%        = 1   -- solve the equation with '+' in the first half of the second equation
%        = -1  -- solve the equation with '-' in the first half of the second equation
%    - sign2   -- swap the sign between both terms in the second half of the second equation
%        = 1   -- solve the equation with '+' in the second half of the second equation
%        = -1  -- solve the equation with '-' in the second half of the second equation
%    - openmp  -- OpenMP selection
%        = 1   -- OpenMP - DAG solver enabled and OpenMP accelerated operations
%        = 0   -- Standard level-3 solver without OpenMP acceleration
% 
% This file file is part of MEPACK, the Matrix Equation Package
% by Martin Koehler.
% 