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

function test_gsylv_single()
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
    C = single(rand( m,  m ));
    B = single(rand( n,  n ));
    D = single(rand( n,  n ));

    X = single(ones ( m , n ));

    Y1 = A*X*B   + C*X*D;
    Y2 = A*X*B'  + C*X*D';
    Y3 = A'*X*B  + C'*X*D;
    Y4 = A'*X*B' + C'*X*D';

    Y5 = A*X*B   - C*X*D;
    Y6 = A*X*B'  - C*X*D';
    Y7 = A'*X*B  - C'*X*D;
    Y8 = A'*X*B' - C'*X*D';

    X1 = mepack_gsylv(A, B, C, D, Y1);
    X2 = mepack_gsylv(A, B, C, D, Y2, 'N', 'T');
    X3 = mepack_gsylv(A, B, C, D, Y3, 'T', 'N');
    X4 = mepack_gsylv(A, B, C, D, Y4, 'T', 'T');

    assert(norm(A * X1 * B   + C * X1 * D   - Y1, 'fro')/norm(Y1, 'fro') < tol,   '[ 1 ] A * X1 *B + C * X1 * D = Y1' );
    assert(norm(A * X2 * B'  + C * X2 * D'  - Y2, 'fro')/norm(Y1, 'fro') < tol,   '[ 1 ] A * X2 *B**T + C * X2 * D**T = Y2' );
    assert(norm(A' * X3 * B  + C' * X3 * D  - Y3, 'fro')/norm(Y1, 'fro') < tol,   '[ 1 ] A**T * X3 *B + C**T * X3 * D = Y3' );
    assert(norm(A' * X4 * B' + C' * X4 * D' - Y4, 'fro')/norm(Y1, 'fro') < tol,   '[ 1 ] A**T * X4 *B**T + C**T * X4 * D**T = Y4' );


    [X1, AS, BS, CS, DS, QA, ZA, QB, ZB ] = mepack_gsylv(A, B, C, D, Y5, optset);
    [X2, AS, BS, CS, DS, QA, ZA, QB, ZB ] = mepack_gsylv(A, B, C, D, Y6, 'N', 'T', optset);
    [X3, AS, BS, CS, DS, QA, ZA, QB, ZB ] = mepack_gsylv(A, B, C, D, Y7, 'T', 'N', optset);
    [X4, AS, BS, CS, DS, QA, ZA, QB, ZB ] = mepack_gsylv(A, B, C, D, Y8, 'T', 'T', optset);

    assert(norm(A * X1 * B   - C * X1 * D   - Y5, 'fro')/norm(Y1, 'fro') < tol,   '[ 2 ] A * X1 *B + C * X1 * D = Y5' );
    assert(norm(A * X2 * B'  - C * X2 * D'  - Y6, 'fro')/norm(Y1, 'fro') < tol,   '[ 2 ] A * X2 *B**T + C * X2 * D**T = Y6' );
    assert(norm(A' * X3 * B  - C' * X3 * D  - Y7, 'fro')/norm(Y1, 'fro') < tol,   '[ 2 ] A**T * X3 *B + C**T * X3 * D = Y7' );
    assert(norm(A' * X4 * B' - C' * X4 * D' - Y8, 'fro')/norm(Y1, 'fro') < tol,   '[ 2 ] A**T * X4 *B**T + C**T * X4 * D**T = Y8' );

    X1 = mepack_gsylv(AS, BS, CS, DS, QA, ZA, QB, ZB, Y5, optset);
    X2 = mepack_gsylv(AS, BS, CS, DS, QA, ZA, QB, ZB, Y6, 'N', 'T', optset);
    X3 = mepack_gsylv(AS, BS, CS, DS, QA, ZA, QB, ZB, Y7, 'T', 'N', optset);
    X4 = mepack_gsylv(AS, BS, CS, DS, QA, ZA, QB, ZB, Y8, 'T', 'T', optset);

    assert(norm(A * X1 * B   - C * X1 * D   - Y5, 'fro')/norm(Y1, 'fro') < tol,   '[ 3 ] A * X1 *B + C * X1 * D = Y5' );
    assert(norm(A * X2 * B'  - C * X2 * D'  - Y6, 'fro')/norm(Y1, 'fro') < tol,   '[ 3 ] A * X2 *B**T + C * X2 * D**T = Y6' );
    assert(norm(A' * X3 * B  - C' * X3 * D  - Y7, 'fro')/norm(Y1, 'fro') < tol,   '[ 3 ] A**T * X3 *B + C**T * X3 * D = Y7' );
    assert(norm(A' * X4 * B' - C' * X4 * D' - Y8, 'fro')/norm(Y1, 'fro') < tol,   '[ 3 ] A**T * X4 *B**T + C**T * X4 * D**T = Y8' );

end
