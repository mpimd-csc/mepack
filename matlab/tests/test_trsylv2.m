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

function test_trsylv2()
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

    optset2 = struct();
    optset2.sign = -1;

    optset3 = optset;
    optset3.sign = -1;

    n = 257;
    m = 198;

    tol = sqrt(eps) * max(m,n);
    tols = sqrt(eps('single')) * max(m,n);
    ierr = 0;

    A = rand( m , m );
    B = rand( n,  n );

    A = schur(A);
    B = schur(B);

    X = ones ( m , n );

    for i = 1:2
        Y1 = A*X*B + X;
        Y2 = A*X*B' + X;
        Y3 = A'*X*B + X;
        Y4 = A'*X*B' + X;

        Y5 = A*X*B - X;
        Y6 = A*X*B' - X;
        Y7 = A'*X*B - X;
        Y8 = A'*X*B'- X;


        if i == 2
            A = single(A);
            B = single(B);
            X = single(X);
            Y1 = single(Y1);
            Y2 = single(Y2);
            Y3 = single(Y3);
            Y4 = single(Y4);

            Y5 = single(Y5);
            Y6 = single(Y6);
            Y7 = single(Y7);
            Y8 = single(Y8);


            tol = tols;
        end

        X1 = mepack_trsylv2(A, B, Y1);
        X2 = mepack_trsylv2(A, B, Y2, 'N', 'T');
        X3 = mepack_trsylv2(A, B, Y3, 'T', 'N');
        X4 = mepack_trsylv2(A, B, Y4, 'T', 'T');

        nrm = norm(A * X1 * B  + X1 - Y1, 'fro')/norm(Y1, 'fro');
        assert(norm(A * X1 * B  + X1 - Y1, 'fro')/norm(Y1, 'fro') < tol, sprintf('A * X1 * B + X1 = Y1 %g',nrm) );
        assert(norm(A * X2 * B' + X2 - Y2, 'fro')/norm(Y2, 'fro') < tol, 'A * X2 * B**T+ X2 = Y2' );
        assert(norm(A' * X3* B  + X3 - Y3, 'fro')/norm(Y3, 'fro') < tol, 'A**T * X3 * B + X3 = Y3' );
        assert(norm(A' * X4* B' + X4 - Y4, 'fro')/norm(Y4, 'fro') < tol, 'A**T * X4 * B**T + X4 = Y4' );

        X1 = mepack_trsylv2(A, B, Y1, optset);
        X2 = mepack_trsylv2(A, B, Y2, 'N', 'T', optset);
        X3 = mepack_trsylv2(A, B, Y3, 'T', 'N', optset);
        X4 = mepack_trsylv2(A, B, Y4, 'T', 'T', optset);

        assert(norm(A * X1 * B  + X1 - Y1, 'fro')/norm(Y1, 'fro') < tol, 'A * X1 * B+ X1 = Y1 opt' );
        assert(norm(A * X2 * B' + X2 - Y2, 'fro')/norm(Y2, 'fro') < tol, 'A * X2 * B**T+ X2 = Y2 opt' );
        assert(norm(A' * X3* B  + X3 - Y3, 'fro')/norm(Y3, 'fro') < tol, 'A**T * X3 * B + X3 = Y3 opt' );
        assert(norm(A' * X4* B' + X4 - Y4, 'fro')/norm(Y4, 'fro') < tol, 'A**T * X4 * B**T + X4 = Y4 opt' );


        X1 = mepack_trsylv2(A, B, Y5, optset2);
        X2 = mepack_trsylv2(A, B, Y6, 'N', 'T', optset2);
        X3 = mepack_trsylv2(A, B, Y7, 'T', 'N', optset2);
        X4 = mepack_trsylv2(A, B, Y8, 'T', 'T', optset2);

        assert(norm(A * X1 * B  - X1 - Y5, 'fro')/norm(Y1, 'fro') < tol, 'A * X1 * B - X1 = Y5 ' );
        assert(norm(A * X2 * B' - X2 - Y6, 'fro')/norm(Y2, 'fro') < tol, 'A * X2 * B**T - X2 = Y6' );
        assert(norm(A' * X3* B  - X3 - Y7, 'fro')/norm(Y3, 'fro') < tol, 'A**T * X3 * B - X3 = Y7' );
        assert(norm(A' * X4* B' - X4 - Y8, 'fro')/norm(Y4, 'fro') < tol, 'A**T * X4 * B**T - X4 = Y8' );

        X1 = mepack_trsylv2(A, B, Y5, optset3);
        X2 = mepack_trsylv2(A, B, Y6, 'N', 'T', optset3);
        X3 = mepack_trsylv2(A, B, Y7, 'T', 'N', optset3);
        X4 = mepack_trsylv2(A, B, Y8, 'T', 'T', optset3);

        assert(norm(A * X1 * B  - X1 - Y5, 'fro')/norm(Y1, 'fro') < tol, 'A * X1 * B - X1 = Y5 opt' );
        assert(norm(A * X2 * B' - X2 - Y6, 'fro')/norm(Y2, 'fro') < tol, 'A * X2 * B**T - X2 = Y6 opt' );
        assert(norm(A' * X3* B  - X3 - Y7, 'fro')/norm(Y3, 'fro') < tol, 'A**T * X3 * B - X3 = Y7 opt' );
        assert(norm(A' * X4* B' - X4 - Y8, 'fro')/norm(Y4, 'fro') < tol, 'A**T * X4 * B**T - X4 = Y8 opt' );


    end
end
