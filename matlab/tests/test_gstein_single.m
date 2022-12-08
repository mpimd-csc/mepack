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

function test_gstein_single()
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
    tol = sqrt(eps("single")) * n;
    ierr = 0;

    A = single(rand( n , n ));
    B = single(rand( n , n ));
    X = single(ones ( n , n ));

    Y1 = A*X*A' - B*X*B';
    Y2 = A'*X*A - B'*X*B;

    X1 = mepack_gstein(A, B, Y1);
    X2 = mepack_gstein(A, B, Y2, 'T');
    X3 = mepack_gstein(A, B, Y1, optset);
    X4 = mepack_gstein(A, B, Y2, 'T', optset);

    assert(norm(A  * X1 * A' - B  * X1 * B' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1 ] A * X1 * A**T - B * X1 * B**T = Y1' );
    assert(norm(A' * X2 * A  - B' * X2 * B  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 1 ] A**T * X2 * A - B**T * X2 * B = Y2' );
    assert(norm(A  * X3 * A' - B  * X3 * B' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1 ] A * X3 * A**T - B*X3 * B**T = Y1' );
    assert(norm(A' * X4 * A  - B' * X4 * B  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 1 ] A**T * X4 * A - B**T * X4 * B = Y2' );

    [ X1, S1, T1, Q1, Z1 ] = mepack_gstein(A, B, Y1);
    [ X2, S2, T2, Q2, Z2 ] = mepack_gstein(A, B, Y2, 'T');
    [ X3, S3, T3, Q3, Z3 ] = mepack_gstein(A, B, Y1, optset);
    [ X4, S4, T4, Q4, Z4 ] = mepack_gstein(A, B, Y2, 'T', optset);

    assert(norm(A  * X1 * A' - B  * X1 * B' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 2 ] A * X1 * A**T - B * X1 * B**T = Y1' );
    assert(norm(A' * X2 * A  - B' * X2 * B  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 2 ] A**T * X2 * A - B**T * X2 * B = Y2' );
    assert(norm(A  * X3 * A' - B  * X3 * B' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 2 ] A * X3 * A**T - B*X3 * B**T = Y1' );
    assert(norm(A' * X4 * A  - B' * X4 * B  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 2 ] A**T * X4 * A - B**T * X4 * B = Y2' );

    X1 = mepack_gstein(S1, T1, Q1, Z1, Y1);
    X2 = mepack_gstein(S2, T2, Q2, Z2, Y2, 'T');
    X3 = mepack_gstein(S3, T3, Q3, Z3, Y1, optset);
    X4 = mepack_gstein(S4, T4, Q4, Z4, Y2, 'T', optset);

    assert(norm(A  * X1 * A' - B  * X1 * B' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 3 ] A * X1 * A**T - B * X1 * B**T = Y1' );
    assert(norm(A' * X2 * A  - B' * X2 * B  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X2 * A - B**T * X2 * B = Y2' );
    assert(norm(A  * X3 * A' - B  * X3 * B' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 3 ] A * X3 * A**T - B*X3 * B**T = Y1' );
    assert(norm(A' * X4 * A  - B' * X4 * B  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X4 * A - B**T * X4 * B = Y2' );



end
