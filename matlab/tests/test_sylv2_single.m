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

function test_sylv2_single()
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

    hess_comp = { 'NN', 'NH', 'HN', 'HH'};
    Ao = single(rand( m , m ));
    Bo = single(rand( n , n ));

    X = single(ones ( m , n ));

    for h = 1:4
        current_hess = hess_comp{h};

        fprintf(1,"Tesing Hessenberg Setup %s\n", current_hess);

        if ( current_hess(1) == 'N')
            A = Ao;
        else
            A = hess(Ao);
        end
        if ( current_hess(2) == 'N')
            B = Bo;
        else
            B = hess(Bo);
        end


        Y1 = A*X*B + X;
        Y2 = A*X*B' + X;
        Y3 = A'*X*B + X;
        Y4 = A'*X*B' + X;

        Y5 = A*X*B - X;
        Y6 = A*X*B' - X;
        Y7 = A'*X*B - X;
        Y8 = A'*X*B'- X;

        X1 = mepack_sylv2(A, B, Y1);
        X2 = mepack_sylv2(A, B, Y2, 'N', 'T');
        X3 = mepack_sylv2(A, B, Y3, 'T', 'N');
        X4 = mepack_sylv2(A, B, Y4, 'T', 'T');

        assert(norm(A * X1 * B  + X1 - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1 ] A * X1 * B    + X1 = Y1' );
        assert(norm(A * X2 * B' + X2 - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 1 ] A * X2 * B**T + X2 = Y2' );
        assert(norm(A' * X3 * B + X3 - Y3, 'fro')/norm(Y3, 'fro') < tol, '[ 1 ] A**T * X3 * B + X3 = Y3' );
        assert(norm(A' * X4 * B'+ X4 - Y4, 'fro')/norm(Y4, 'fro') < tol, '[ 1 ] A**T * X4 * B**T +  X4 = Y4' );

        [X1, S1, T1, Q1, Z1 ] = mepack_sylv2(A, B, Y5, optset);
        [X2, S2, T2, Q2, Z2 ] = mepack_sylv2(A, B, Y6, 'N', 'T', optset);
        [X3, S3, T3, Q3, Z3 ] = mepack_sylv2(A, B, Y7, 'T', 'N', optset);
        [X4, S4, T4, Q4, Z4 ] = mepack_sylv2(A, B, Y8, 'T', 'T', optset);

        assert(norm(A * X1  * B  - X1 - Y5, 'fro')/norm(Y5, 'fro') < tol,  '[ 2 ] A * X1 * B    - X1 = Y5 opt' );
        assert(norm(A * X2  * B' - X2 - Y6, 'fro')/norm(Y6, 'fro') < tol,  '[ 2 ] A * X2 * B**T - X2 = Y6 opt' );
        assert(norm(A' * X3 * B  - X3 - Y7, 'fro')/norm(Y7, 'fro') < tol,  '[ 2 ] A**T * X3 * B - X3 = Y7 opt' );
        assert(norm(A' * X4 * B' - X4 - Y8, 'fro')/norm(Y8, 'fro') < tol,  '[ 2 ] A**T * X4 * B**T - X4 = Y8 opt' );

        X1 = mepack_sylv2(S1, T1, Q1, Z1, Y5, optset);
        X2 = mepack_sylv2(S2, T2, Q2, Z2, Y6, 'N', 'T', optset);
        X3 = mepack_sylv2(S3, T3, Q3, Z3, Y7, 'T', 'N', optset);
        X4 = mepack_sylv2(S4, T4, Q4, Z4, Y8, 'T', 'T', optset);

        assert(norm(A * X1  * B  - X1 - Y5, 'fro')/norm(Y5, 'fro') < tol,  '[ 3 ] A * X1 * B    - X1 = Y5 opt' );
        assert(norm(A * X2  * B' - X2 - Y6, 'fro')/norm(Y6, 'fro') < tol,  '[ 3 ] A * X2 * B**T - X2 = Y6 opt' );
        assert(norm(A' * X3 * B  - X3 - Y7, 'fro')/norm(Y7, 'fro') < tol,  '[ 3 ] A**T * X3 * B - X3 = Y7 opt' );
        assert(norm(A' * X4 * B' - X4 - Y8, 'fro')/norm(Y8, 'fro') < tol,  '[ 3 ] A**T * X4 * B**T - X4 = Y8 opt' );

    end
end
