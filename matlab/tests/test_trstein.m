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

function test_trstein()
     mepack_test_init_random();

    is_octave =  exist('OCTAVE_VERSION', 'builtin');
    if ( is_octave )
        openmp = 1;
    else
        openmp = 0;
    end
    optset = struct();
    optset.mb = 16;
    optset.nb = 32;
    optset.isolver = 3;
    optset.openmp = openmp;

    n = 257;
    tol = sqrt(eps) * n;
    tols = sqrt(eps('single')) * n;
    ierr = 0;

    A = rand( n , n );
    A = schur(A);
    X = ones ( n , n );

    Y1 = A*X + X*A';
    Y2 = A'*X + X*A;

    As  = single(A);
    XS  = single(X);
    Y1s = single(Y1);
    Y2s = single(Y2);

    X1 = mepack_trstein(A, Y1);
    X2 = mepack_trstein(A, Y2, 'T');
    X3 = mepack_trstein(A, Y1, optset);
    X4 = mepack_trstein(A, Y2, 'T', optset);


    X1s = mepack_trstein(As, Y1s);
    X2s = mepack_trstein(As, Y2s, 'T');
    X3s = mepack_trstein(As, Y1s, optset );
    X4s = mepack_trstein(As, Y2s, 'T', optset);


    assert(norm(A * X1 * A' - X1 - Y1, 'fro')/norm(Y1, 'fro') < tol, 'A * X1 * A**T - X1 = Y1' );
    assert(norm(A' * X2 * A - X2 - Y2, 'fro')/norm(Y2, 'fro') < tol, 'A**T * X2 * A - X2 = Y2' );
    assert(norm(A * X3 * A' - X3 - Y1, 'fro')/norm(Y1, 'fro') < tol, 'A * X3 * A**T - X3 = Y1' );
    assert(norm(A' * X4 * A - X4 - Y2, 'fro')/norm(Y2, 'fro') < tol, 'A**T * X4 * A - X4 = Y2' );

    assert(norm(As * X1s * As' - X1s - Y1s, 'fro')/norm(Y1s, 'fro') < tols, 'As * X1s * As**T - X1s = Y1s' );
    assert(norm(As' * X2s * As - X2s - Y2s, 'fro')/norm(Y2s, 'fro') < tols, 'As**T * X2s * As - X2s = Y2s' );
    assert(norm(As * X3s * As' - X3s - Y1s, 'fro')/norm(Y1s, 'fro') < tols, 'As * X3s * As**T - X3s = Y1s' );
    assert(norm(As' * X4s * As - X4s - Y2s, 'fro')/norm(Y2s, 'fro') < tols, 'As**T * X4s * As - X4s = Y2s' );

