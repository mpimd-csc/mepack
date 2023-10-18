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

function test_glyap()
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

    hess_comp = { 'N', 'H'};

    Ao = rand( n , n );
    Bo = rand( n , n );
    X = ones ( n , n );

    for h = 1:2
        current_hess = hess_comp{h};

        fprintf(1,"Tesing Hessenberg Setup %s\n", current_hess);

        if ( is_octave )
            hess_fun = @qzhess;
        else
            hess_fun = @hess;
        end

        if ( current_hess(1) == 'N')
            A = Ao;
            B = Bo;
        else
            [A, B, ~, ~ ] = hess_fun(Ao, Bo);
        end

        Y1 = A*X*B' + B*X*A';
        Y2 = A'*X*B + B'*X*A;

        X1 = mepack_glyap(A, B, Y1);
        X2 = mepack_glyap(A, B, Y2, 'T');
        X3 = mepack_glyap(A, B, Y1, optset);
        X4 = mepack_glyap(A, B, Y2, 'T', optset);

        assert(norm(A  * X1 * B' + B  * X1 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1 ] A * X1 * B**T + B * X1 * A**T = Y1' );
        assert(norm(A' * X2 * B  + B' * X2 * A  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 1 ] A**T * X2 * B + B**T * X2 * A = Y2' );
        assert(norm(A  * X3 * B' + B  * X3 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 1 ] A * X3 * B**T + B*X3 * A**T = Y1' );
        assert(norm(A' * X4 * B  + B' * X4 * A  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 1 ] A**T * X4 * B+ B**T * X4 * A = Y2' );

        [ X1, S1, T1, Q1, Z1 ] = mepack_glyap(A, B, Y1);
        [ X2, S2, T2, Q2, Z2 ] = mepack_glyap(A, B, Y2, 'T');
        [ X3, S3, T3, Q3, Z3 ] = mepack_glyap(A, B, Y1, optset);
        [ X4, S4, T4, Q4, Z4 ] = mepack_glyap(A, B, Y2, 'T', optset);

        assert(norm(A  * X1 * B' + B  * X1 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 3 ] A * X1 * B**T + B * X1 * A**T = Y1' );
        assert(norm(A' * X2 * B  + B' * X2 * A  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X2 * B + B**T * X2 * A = Y2' );
        assert(norm(A  * X3 * B' + B  * X3 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 3 ] A * X3 * B**T + B*X3 * A**T = Y1' );
        assert(norm(A' * X4 * B  + B' * X4 * A  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X4 * B+ B**T * X4 * A = Y2' );

        X1 = mepack_glyap(S1, T1, Q1, Z1, Y1);
        X2 = mepack_glyap(S2, T2, Q2, Z2, Y2, 'T');
        X3 = mepack_glyap(S3, T3, Q3, Z3, Y1, optset);
        X4 = mepack_glyap(S4, T4, Q4, Z4, Y2, 'T', optset);

        assert(norm(A  * X1 * B' + B  * X1 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 3 ] A * X1 * B**T + B * X1 * A**T = Y1' );
        assert(norm(A' * X2 * B  + B' * X2 * A  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X2 * B + B**T * X2 * A = Y2' );
        assert(norm(A  * X3 * B' + B  * X3 * A' - Y1, 'fro')/norm(Y1, 'fro') < tol, '[ 3 ] A * X3 * B**T + B*X3 * A**T = Y1' );
        assert(norm(A' * X4 * B  + B' * X4 * A  - Y2, 'fro')/norm(Y2, 'fro') < tol, '[ 3 ] A**T * X4 * B+ B**T * X4 * A = Y2' );

    end

end
