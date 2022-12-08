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

function test_sylv_single()
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
    optset.sign = -1;


    n = 257;
    m = 199;

    tol = sqrt(eps('single')) * max(m,n);
    ierr = 0;

    A = single(rand( m , m ));
    B = single(rand( n,  n ));

    X = single(ones ( m , n ));

    Y1 = A*X + X*B;
    Y2 = A*X + X*B';
    Y3 = A'*X + X*B;
    Y4 = A'*X + X*B';

    Y5 = A*X - X*B;
    Y6 = A*X - X*B';
    Y7 = A'*X - X*B;
    Y8 = A'*X - X*B';


    X1 = mepack_sylv(A, B, Y1);
    X2 = mepack_sylv(A, B, Y2, 'N', 'T');
    X3 = mepack_sylv(A, B, Y3, 'T', 'N');
    X4 = mepack_sylv(A, B, Y4, 'T', 'T');

    assert(norm(A * X1 + X1 * B - Y1, 'fro')/norm(Y1, 'fro') < tol,   '[ 1 ] A * X1 + X1 * B = Y1' );
    assert(norm(A * X2 + X2 * B' - Y2, 'fro')/norm(Y2, 'fro') < tol,  '[ 1 ] A * X2 + X2 * B**T = Y2' );
    assert(norm(A' * X3 + X3 * B - Y3, 'fro')/norm(Y3, 'fro') < tol,  '[ 1 ] A**T * X3 + X3 * B = Y3' );
    assert(norm(A' * X4 + X4 * B' - Y4, 'fro')/norm(Y4, 'fro') < tol, '[ 1 ] A**T * X4 + X4 * B**T = Y4' );

    [X1, S1, T1, Q1, Z1 ] = mepack_sylv(A, B, Y5, optset);
    [X2, S2, T2, Q2, Z2 ] = mepack_sylv(A, B, Y6, 'N', 'T', optset);
    [X3, S3, T3, Q3, Z3 ] = mepack_sylv(A, B, Y7, 'T', 'N', optset);
    [X4, S4, T4, Q4, Z4 ] = mepack_sylv(A, B, Y8, 'T', 'T', optset);

    assert(norm(A * X1  - X1 * B  - Y5, 'fro')/norm(Y5, 'fro') < tol,   '[ 2 ] A * X1 - X1 * B = Y5 opt' );
    assert(norm(A * X2  - X2 * B' - Y6, 'fro')/norm(Y6, 'fro') < tol,  '[ 2 ] A * X2 - X2 * B**T = Y6 opt' );
    assert(norm(A' * X3 - X3 * B  - Y7, 'fro')/norm(Y7, 'fro') < tol,   '[ 2 ] A**T * X3 - X3 * B = Y7 opt' );
    assert(norm(A' * X4 - X4 * B' - Y8, 'fro')/norm(Y8, 'fro') < tol,  '[ 2 ] A**T * X4 - X4 * B**T = Y8 opt' );

    X1 = mepack_sylv(S1, T1, Q1, Z1, Y5, optset);
    X2 = mepack_sylv(S2, T2, Q2, Z2, Y6, 'N', 'T', optset);
    X3 = mepack_sylv(S3, T3, Q3, Z3, Y7, 'T', 'N', optset);
    X4 = mepack_sylv(S4, T4, Q4, Z4, Y8, 'T', 'T', optset);

    assert(norm(A * X1  - X1 * B  - Y5, 'fro')/norm(Y5, 'fro') < tol,   '[ 3 ] A * X1 - X1 * B = Y5 opt' );
    assert(norm(A * X2  - X2 * B' - Y6, 'fro')/norm(Y6, 'fro') < tol,  '[ 3 ] A * X2 - X2 * B**T = Y6 opt' );
    assert(norm(A' * X3 - X3 * B  - Y7, 'fro')/norm(Y7, 'fro') < tol,   '[ 3 ] A**T * X3 - X3 * B = Y7 opt' );
    assert(norm(A' * X4 - X4 * B' - Y8, 'fro')/norm(Y8, 'fro') < tol,  '[ 3 ] A**T * X4 - X4 * B**T = Y8 opt' );



end
