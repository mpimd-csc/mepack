% Copyright (C) Martin Koehler, 2017-2023
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.
%

function test_lyap_refine()
     mepack_test_init_random();
    is_octave =  exist('OCTAVE_VERSION', 'builtin');
    if ( is_octave )
        openmp = 1;
        isolver = 3;
    else
        isolver = 1;
        openmp = 0;
    end
    optset = struct();
    optset.mb = 16;
    optset.nb = 32;
    optset.isolver = isolver;
    optset.openmp = openmp;

    n = 257;
    tol = sqrt(eps) * n;
    ierr = 0;

    A = rand( n , n );
    X = ones ( n , n );
    [Q, S] = schur(A, 'real');

    Y1 = A*X + X*A';
    Y2 = A'*X + X*A;

    X0 = X + rand(n);

    X1 = mepack_lyap_refine(A, S, Q, Y1);
    X2 = mepack_lyap_refine(A, S, Q, Y1, X0);
    X3 = mepack_lyap_refine(A, S, Q, Y2, 'T');
    X4 = mepack_lyap_refine(A, S, Q, Y2, X0, 'T');
    X5 = mepack_lyap_refine(A, S, Q, Y1, optset);
    X6 = mepack_lyap_refine(A, S, Q, Y1, X0, optset);
    X7 = mepack_lyap_refine(A, S, Q, Y2, 'T', optset);
    X8 = mepack_lyap_refine(A, S, Q, Y2, X0, 'T', optset);

    assert(norm(A * X1 + X1 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1 ] A * X1 + X1 * A**T = Y1' );
    assert(norm(A * X2 + X2 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 2 ] A * X2 + X2 * A**T = Y1' );
    assert(norm(A' * X3 + X3 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X3 + X3 * A = Y2' );
    assert(norm(A' * X4 + X4 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 4 ] A**T * X4 + X4 * A = Y2' );

    assert(norm(A * X5 + X5 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 5 ] A * X5 + X5 * A**T = Y1' );
    assert(norm(A * X6 + X6 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 6 ] A * X6 + X6 * A**T = Y1' );
    assert(norm(A' * X7 + X7 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 7 ] A**T * X7 + X7 * A = Y2' );
    assert(norm(A' * X8 + X8 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 8 ] A**T * X8 + X8 * A = Y2' );


    [ X1, convlog1 ] = mepack_lyap_refine(A, S, Q, Y1);
    [ X2, convlog2 ] = mepack_lyap_refine(A, S, Q, Y1, X0);
    [ X3, convlog3 ] = mepack_lyap_refine(A, S, Q, Y2, 'T');
    [ X4, convlog4 ] = mepack_lyap_refine(A, S, Q, Y2, X0, 'T');
    [ X5, convlog5 ] = mepack_lyap_refine(A, S, Q, Y1, optset);
    [ X6, convlog6 ] = mepack_lyap_refine(A, S, Q, Y1, X0, optset);
    [ X7, convlog7 ] = mepack_lyap_refine(A, S, Q, Y2, 'T', optset);
    [ X8, convlog8 ] = mepack_lyap_refine(A, S, Q, Y2, X0, 'T', optset);

    assert(norm(A * X1 + X1 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1c ] A * X1 + X1 * A**T = Y1' );
    assert(norm(A * X2 + X2 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 2c ] A * X2 + X2 * A**T = Y1' );
    assert(norm(A' * X3 + X3 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3c ] A**T * X3 + X3 * A = Y2' );
    assert(norm(A' * X4 + X4 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 4c ] A**T * X4 + X4 * A = Y2' );

    assert(norm(A * X5 + X5 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 5c ] A * X5 + X5 * A**T = Y1' );
    assert(norm(A * X6 + X6 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 6c ] A * X6 + X6 * A**T = Y1' );
    assert(norm(A' * X7 + X7 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 7c ] A**T * X7 + X7 * A = Y2' );
    assert(norm(A' * X8 + X8 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 8c ] A**T * X8 + X8 * A = Y2' );


end
