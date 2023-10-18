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
function test_cgsylv_dual()
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
    m = 198;

    tol = sqrt(eps) * max(m,n);
    ierr = 0;

    hess_comp = { 'NN', 'NH', 'HN', 'HH'};

    Ao = rand( m , m );
    Co = rand( m , m );
    Bo = rand( n , n );
    Do = rand( n , n );

    R = ones ( m , n );
    L = ones ( m , n );

    for h = 1:4
        current_hess = hess_comp{h};

        fprintf(1,"Tesing Hessenberg Setup %s\n", current_hess);

        if ( is_octave )
            hess_fun = @qzhess;
        else
            hess_fun = @hess;
        end

        if ( current_hess(1) == 'N')
            A = Ao;
            C = Co;
        else
            [A, C, ~, ~ ] = hess_fun(Ao, Co);
        end
        if ( current_hess(2) == 'N')
            B = Bo;
            D = Do;
        else
            [B, D, ~, ~ ] = hess_fun(Bo, Do);
        end

        for i = 1:4
            if i == 1
                optset.sign = 1;
                optset.sign2 = 1;
            elseif i == 2
                optset.sign = 1;
                optset.sign2 = -1;
            elseif i == 3
                optset.sign = -1;
                optset.sign2 = 1;
            elseif i == 4
                optset.sign = -1;
                optset.sign2 = -1;
            end

            E1 = A'*R + C'*L;
            F1 = optset.sign * R * B' + optset.sign2 * L * D';
            E2 = A*R + C*L;
            F2 = optset.sign * R * B' + optset.sign2 * L * D';
            E3 = A'*R + C'*L;
            F3 = optset.sign * R * B + optset.sign2 * L * D;
            E4 = A*R + C*L;
            F4 = optset.sign * R * B + optset.sign2 * L * D;

            if i == 1
                [R1, L1] = mepack_csylv_dual(A, B, C, D, E1, F1);
                [R2, L2] = mepack_csylv_dual(A, B, C, D, E2, F2, 'N', 'T');
                [R3, L3] = mepack_csylv_dual(A, B, C, D, E3, F3, 'T', 'N');
                [R4, L4] = mepack_csylv_dual(A, B, C, D, E4, F4, 'T', 'T');
            else
                [R1, L1] = mepack_csylv_dual(A, B, C, D, E1, F1, optset);
                [R2, L2] = mepack_csylv_dual(A, B, C, D, E2, F2, 'N', 'T', optset);
                [R3, L3] = mepack_csylv_dual(A, B, C, D, E3, F3, 'T', 'N', optset);
                [R4, L4] = mepack_csylv_dual(A, B, C, D, E4, F4, 'T', 'T', optset);
            end
            assert(max( norm ( A' * R1  +  C' * L1 - E1, 'fro'), norm( optset.sign * R * B' + optset.sign2 * L * D' - F1, 'fro')) ...
                /max(norm(E1, 'fro'), norm(F1,'fro')) < tol, sprintf('[%d-1] E1 - F1 failed', i) );
            assert(max( norm ( A * R1  +  C * L1 - E2, 'fro'), norm( optset.sign * R * B' + optset.sign2 * L * D' - F2, 'fro')) ...
                /max(norm(E2, 'fro'), norm(F2,'fro')) < tol, sprintf('[%d-1] E2 - F2 failed', i) );
            assert(max( norm ( A' * R1  +  C' * L1 - E3, 'fro'), norm( optset.sign * R * B + optset.sign2 * L * D - F3, 'fro')) ...
                /max(norm(E3, 'fro'), norm(F3,'fro')) < tol, sprintf('[%d-1] E3 - F3 failed', i) );
            assert(max( norm ( A * R1  +  C * L1 - E4, 'fro'), norm( optset.sign * R * B + optset.sign2 * L * D - F4, 'fro')) ...
                /max(norm(E4, 'fro'), norm(F4,'fro')) < tol, sprintf('[%d-1] E4 - F4 failed', i) );



            [R1, L1, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E1, F1, optset);
            [R2, L2, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E2, F2, 'N', 'T', optset);
            [R3, L3, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E3, F3, 'T', 'N', optset);
            [R4, L4, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E4, F4, 'T', 'T', optset);

            assert(max( norm ( A' * R1  +  C' * L1 - E1, 'fro'), norm( optset.sign * R * B' + optset.sign2 * L * D' - F1, 'fro')) ...
                /max(norm(E1, 'fro'), norm(F1,'fro')) < tol, sprintf('[%d-1] E1 - F1 failed', i) );
            assert(max( norm ( A * R1  +  C * L1 - E2, 'fro'), norm( optset.sign * R * B' + optset.sign2 * L * D' - F2, 'fro')) ...
                /max(norm(E2, 'fro'), norm(F2,'fro')) < tol, sprintf('[%d-1] E2 - F2 failed', i) );
            assert(max( norm ( A' * R1  +  C' * L1 - E3, 'fro'), norm( optset.sign * R * B + optset.sign2 * L * D - F3, 'fro')) ...
                /max(norm(E3, 'fro'), norm(F3,'fro')) < tol, sprintf('[%d-1] E3 - F3 failed', i) );
            assert(max( norm ( A * R1  +  C * L1 - E4, 'fro'), norm( optset.sign * R * B + optset.sign2 * L * D - F4, 'fro')) ...
                /max(norm(E4, 'fro'), norm(F4,'fro')) < tol, sprintf('[%d-1] E4 - F4 failed', i) );


            [R1, L1] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E1, F1, optset);
            [R2, L2] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E2, F2, 'N', 'T', optset);
            [R3, L3] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E3, F3, 'T', 'N', optset);
            [R4, L4] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E4, F4, 'T', 'T', optset);

            assert(max( norm ( A' * R1  +  C' * L1 - E1, 'fro'), norm( optset.sign * R * B' + optset.sign2 * L * D' - F1, 'fro')) ...
                /max(norm(E1, 'fro'), norm(F1,'fro')) < tol, sprintf('[%d-1] E1 - F1 failed', i) );
            assert(max( norm ( A * R1  +  C * L1 - E2, 'fro'), norm( optset.sign * R * B' + optset.sign2 * L * D' - F2, 'fro')) ...
                /max(norm(E2, 'fro'), norm(F2,'fro')) < tol, sprintf('[%d-1] E2 - F2 failed', i) );
            assert(max( norm ( A' * R1  +  C' * L1 - E3, 'fro'), norm( optset.sign * R * B + optset.sign2 * L * D - F3, 'fro')) ...
                /max(norm(E3, 'fro'), norm(F3,'fro')) < tol, sprintf('[%d-1] E3 - F3 failed', i) );
            assert(max( norm ( A * R1  +  C * L1 - E4, 'fro'), norm( optset.sign * R * B + optset.sign2 * L * D - F4, 'fro')) ...
                /max(norm(E4, 'fro'), norm(F4,'fro')) < tol, sprintf('[%d-1] E4 - F4 failed', i) );

        end
    end
end
