%% 
% X = mepack_lyap_refine(A, S, Q, Y)
% X = mepack_lyap_refine(A, S, Q, Y, X0)
% X = mepack_lyap_refine(A, S, Q, Y, Op)
% X = mepack_lyap_refine(A, S, Q, Y, X0, Op)
% X = mepack_lyap_refine(A, S, Q, Y, OptSet)
% X = mepack_lyap_refine(A, S, Q, Y, X0, OptSet)
% X = mepack_lyap_refine(A, S, Q, Y, Op, OptSet)
% X = mepack_lyap_refine(A, S, Q, Y, X0, Op, OptSet)
% 
% Input:
%     A      -- coefficient matrix
%     S      -- Schur factor of the coefficient matrix
%     Q      -- Unitary transformation of the Schur decomposition
%     Y      -- right hand side
%     X0     -- Initial Guess
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
%  Return:
%     X      -- Solution of the Lyapunov equation
% 
% [X, convlog] = mepack_lyap_refine(A, S, Q, Y)
% [X, convlog] = mepack_lyap_refine(A, S, Q, Y, X0)
% [X, convlog] = mepack_lyap_refine(A, S, Q, Y, Op)
% [X, convlog] = mepack_lyap_refine(A, S, Q, Y, X0, Op)
% [X, convlog] = mepack_lyap_refine(A, S, Q, Y, OptSet)
% [X, convlog] = mepack_lyap_refine(A, S, Q, Y, X0, OptSet)
% [X, convlog] = mepack_lyap_refine(A, S, Q, Y, Op, OptSet)
% [X, convlog] = mepack_lyap_refine(A, S, Q, Y, X0, Op, OptSet)
% 
% Input:
%     A      -- coefficient matrix
%     S      -- Schur factor of the coefficient matrix
%     Q      -- Unitary transformation of the Schur decomposition
%     Y      -- right hand side
%     X0     -- Initial Guess
%     Op     -- Transpose operation on the coefficient matrix
%     OptSet -- Options structure to override the default settings
%               of MEPACK.
%  Return:
%     X        -- Solution of the Lyapunov equation
%     convlog  -- convergence history
% 
% 
% The mepack_lyap_refine function solves the Lyapunov equation
% 
%    op(A) * X + X * op(A)' = Y,
% 
% where A is a general matrix employing the iterative refinement
% strategy. The operation op either returns the matrix or
% if Op == 'T' its transpose. In this way, the transposed equation
% can be solved without transposing the coefficient matrix before.
% The function requires the matrix A and its Schur decomposition
% 
%    [ Q, S ] = schur (A, 'real')
% 
% in order to work. The Schur decomposition can also be obtain from
% the output of a previous run to mepack_lyap.
% The iterative refinement stops once
% 
%   || op(A) * X_k + X_k * op(A)' - Y ||_F
%   --------------------------------------  < 2 * m * ||A||_F * u * tau
%              || Y ||_F
% 
% is fulfilled. Thereby m is the order of the matrix A, u is the
% machine precision in double or single precision, and tau is a
% security factor. Tau can be changed in the OptSet.
% 
% The matrices A and Y must be of the same size and Y must be symmetric.
% 
% The function uses the level-3 solver D/SLA_GELYAP_REFINE of MEPACK.
% 
% If the OptSet argument is given, the default settings of MEPACK can
% be overwritten this structure can have the following members for
% mepack_lyap_refine:
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