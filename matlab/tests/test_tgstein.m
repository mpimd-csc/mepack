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

function test_tgstein()
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
    B = rand( n,  n );

    is_octave =  exist('OCTAVE_VERSION', 'builtin');
    if (is_octave)
        [A, B] = qz(A, B);
    else
        [A, B] = qz(A, B, 'real');
    end

    X = ones ( n , n );

    Y1 = A*X*A' - B*X*B';
    Y2 = A'*X*A - B'*X*B;

    As  = single(A);
    Bs  = single(B);
    XS  = single(X);
    Y1s = single(Y1);
    Y2s = single(Y2);

    X1 = mepack_tgstein(A, B, Y1);
    X2 = mepack_tgstein(A, B, Y2, 'T');
    X3 = mepack_tgstein(A, B, Y1, optset);
    X4 = mepack_tgstein(A, B, Y2, 'T', optset);


    X1s = mepack_tgstein(As, Bs, Y1s);
    X2s = mepack_tgstein(As, Bs, Y2s, 'T');
    X3s = mepack_tgstein(As, Bs, Y1s, optset );
    X4s = mepack_tgstein(As, Bs, Y2s, 'T', optset);


    assert(norm(A * X1 * A'  - B*X1 * B' - Y1, 'fro')/norm(Y1, 'fro') < tol, 'A * X1 * A**T - B * X1 * B**T = Y1' );
    assert(norm(A' * X2 * A  - B'*X2 * B - Y2, 'fro')/norm(Y2, 'fro') < tol, 'A**T * X2 * A - B**T X2 * B   = Y2' );
    assert(norm(A * X3 * A'  - B*X3 * B' - Y1, 'fro')/norm(Y1, 'fro') < tol, 'A * X3 * A**T - B* X3 * B**T = Y1' );
    assert(norm(A' * X4 * A  - B'*X4 * B - Y2, 'fro')/norm(Y2, 'fro') < tol, 'A**T * X4 * A - B**T * X4 * B = Y2' );

    assert(norm(As * X1s * As' - Bs * X1s * Bs' - Y1s, 'fro')/norm(Y1s, 'fro') < tols, 'As * X1s * As**T - Bs * X1s * Bs**T = Y1s' );
    assert(norm(As' * X2s * As - Bs'* X2s * Bs  - Y2s, 'fro')/norm(Y2s, 'fro') < tols, 'As**T * X2s * As - Bs**T * X2s * Bs = Y2s' );
    assert(norm(As * X3s * As' - Bs * X3s * Bs' - Y1s, 'fro')/norm(Y1s, 'fro') < tols, 'As * X3s * As**T - Bs * X3s * Bs**T = Y1s' );
    assert(norm(As' * X4s * As - Bs' *X4s * Bs  - Y2s, 'fro')/norm(Y2s, 'fro') < tols, 'As**T * X4s * As - Bs ** T * X4s * Bs = Y2s' );




end
