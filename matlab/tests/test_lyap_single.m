% Copyright (C) Martin Koehler, 2017-2022
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

function test_lyap_single()
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
    tol = single(sqrt(eps) * n);
    ierr = 0;

    A = single(rand( n , n ));
    X = single(ones ( n , n ));

    Y1 = A*X + X*A';
    Y2 = A'*X + X*A;

    X1 = mepack_lyap(A, Y1);
    X2 = mepack_lyap(A, Y2, 'T');
    X3 = mepack_lyap(A, Y1, optset);
    X4 = mepack_lyap(A, Y2, 'T', optset);

    assert(norm(A * X1 + X1 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1 ] A * X1 + X1 * A**T = Y1' );
    assert(norm(A' * X2 + X2 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 1 ] A**T * X2 + X2 * A = Y2' );
    assert(norm(A * X3 + X3 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1 ] A * X3 + X3 * A**T = Y1' );
    assert(norm(A' * X4 + X4 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 1 ] A**T * X4 + X4 * A = Y2' );

    [ X1, S1, Q1 ] = mepack_lyap(A, Y1);
    [ X2, S2, Q2 ] = mepack_lyap(A, Y2, 'T');
    [ X3, S3, Q3 ] = mepack_lyap(A, Y1, optset);
    [ X4, S4, Q4 ] = mepack_lyap(A, Y2, 'T', optset);

    assert(norm(A * X1 + X1 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 2 ] A * X1 + X1 * A**T = Y1' );
    assert(norm(A' * X2 + X2 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 2 ] A**T * X2 + X2 * A = Y2' );
    assert(norm(A * X3 + X3 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 2 ] A * X3 + X3 * A**T = Y1' );
    assert(norm(A' * X4 + X4 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 2 ] A**T * X4 + X4 * A = Y2' );

    X1 = mepack_lyap(S1, Q1 , Y1);
    X2 = mepack_lyap(S2, Q2 , Y2, 'T');
    X3 = mepack_lyap(S3, Q3 , Y1, optset);
    X4 = mepack_lyap(S4, Q4 , Y2, 'T', optset);

    assert(norm(A * X1 + X1 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 3 ] A * X1 + X1 * A**T = Y1' );
    assert(norm(A' * X2 + X2 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X2 + X2 * A = Y2' );
    assert(norm(A * X3 + X3 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 3 ] A * X3 + X3 * A**T = Y1' );
    assert(norm(A' * X4 + X4 * A - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X4 + X4 * A = Y2' );


end
