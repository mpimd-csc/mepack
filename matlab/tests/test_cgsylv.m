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
function test_csylv()
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

    A = rand( m , m );
    C = rand( m , m );
    B = rand( n , n );
    D = rand( n , n );

    R = ones ( m , n );
    L = ones ( m , n );

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

        E1 = A*R + optset.sign * L*B;
        F1 = C*R + optset.sign2 * L*D;
        E2 = A*R + optset.sign * L*B';
        F2 = C*R + optset.sign2 * L*D';
        E3 = A'*R + optset.sign * L*B;
        F3 = C'*R + optset.sign2 * L*D;
        E4 = A'*R + optset.sign * L*B';
        F4 = C'*R + optset.sign2 * L*D';

        if i == 1
            [R1, L1] = mepack_csylv(A, B, C, D, E1, F1);
            [R2, L2] = mepack_csylv(A, B, C, D, E2, F2, 'N', 'T');
            [R3, L3] = mepack_csylv(A, B, C, D, E3, F3, 'T', 'N');
            [R4, L4] = mepack_csylv(A, B, C, D, E4, F4, 'T', 'T');
        else
            [R1, L1] = mepack_csylv(A, B, C, D, E1, F1, optset);
            [R2, L2] = mepack_csylv(A, B, C, D, E2, F2, 'N', 'T', optset);
            [R3, L3] = mepack_csylv(A, B, C, D, E3, F3, 'T', 'N', optset);
            [R4, L4] = mepack_csylv(A, B, C, D, E4, F4, 'T', 'T', optset);
        end
        assert(max( norm ( A * R1  +  optset.sign * L1 * B  - E1, 'fro'), norm( C * R1  +  optset.sign2 * L1 * D  - F1, 'fro')) ...
                /max(norm(E1, 'fro'), norm(F1,'fro')) < tol, sprintf('[%d-1] E1 - F1 failed', i) );
        assert(max( norm ( A * R2  +  optset.sign * L2 * B' - E2, 'fro'), norm( C * R2  +  optset.sign2 * L2 * D' - F2, 'fro')) ...
                /max(norm(E2, 'fro'), norm(F2,'fro')) < tol, sprintf('[%d-1] E2 - F2 failed', i) );
        assert(max( norm ( A' * R3 +  optset.sign * L3 * B  - E3, 'fro'), norm( C' * R3 +  optset.sign2 * L3 * D  - F3, 'fro')) ...
                /max(norm(E3, 'fro'), norm(F3,'fro')) < tol, sprintf('[%d-1] E3 - F3 failed', i) );
        assert(max( norm ( A' * R4 +  optset.sign * L4 * B' - E4, 'fro'), norm( C' * R4 +  optset.sign2 * L4 * D' - F4, 'fro')) ...
                /max(norm(E4, 'fro'), norm(F4,'fro')) < tol, sprintf('[%d-1] E4 - F4 failed', i) );


        [R1, L1, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv(A, B, C, D, E1, F1, optset);
        [R2, L2, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv(A, B, C, D, E2, F2, 'N', 'T', optset);
        [R3, L3, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv(A, B, C, D, E3, F3, 'T', 'N', optset);
        [R4, L4, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv(A, B, C, D, E4, F4, 'T', 'T', optset);

        assert(max( norm ( A * R1  +  optset.sign * L1 * B  - E1, 'fro'), norm( C * R1  +  optset.sign2 * L1 * D  - F1, 'fro')) ...
                /max(norm(E1, 'fro'), norm(F1,'fro')) < tol, sprintf('[%d-2] E1 - F1 failed', i) );
        assert(max( norm ( A * R2  +  optset.sign * L2 * B' - E2, 'fro'), norm( C * R2  +  optset.sign2 * L2 * D' - F2, 'fro')) ...
                /max(norm(E2, 'fro'), norm(F2,'fro')) < tol, sprintf('[%d-2] E2 - F2 failed', i) );
        assert(max( norm ( A' * R3 +  optset.sign * L3 * B  - E3, 'fro'), norm( C' * R3 +  optset.sign2 * L3 * D  - F3, 'fro')) ...
                /max(norm(E3, 'fro'), norm(F3,'fro')) < tol, sprintf('[%d-2] E3 - F3 failed', i) );
        assert(max( norm ( A' * R4 +  optset.sign * L4 * B' - E4, 'fro'), norm( C' * R4 +  optset.sign2 * L4 * D' - F4, 'fro')) ...
                /max(norm(E4, 'fro'), norm(F4,'fro')) < tol, sprintf('[%d-2] E4 - F4 failed', i) );

        [R1, L1] = mepack_csylv(AS, BS, CS, DS, QA, ZA, QB, ZB, E1, F1, optset);
        [R2, L2] = mepack_csylv(AS, BS, CS, DS, QA, ZA, QB, ZB, E2, F2, 'N', 'T', optset);
        [R3, L3] = mepack_csylv(AS, BS, CS, DS, QA, ZA, QB, ZB, E3, F3, 'T', 'N', optset);
        [R4, L4] = mepack_csylv(AS, BS, CS, DS, QA, ZA, QB, ZB, E4, F4, 'T', 'T', optset);

        assert(max( norm ( A * R1  +  optset.sign * L1 * B  - E1, 'fro'), norm( C * R1  +  optset.sign2 * L1 * D  - F1, 'fro')) ...
                /max(norm(E1, 'fro'), norm(F1,'fro')) < tol, sprintf('[%d-3] E1 - F1 failed', i) );
        assert(max( norm ( A * R2  +  optset.sign * L2 * B' - E2, 'fro'), norm( C * R2  +  optset.sign2 * L2 * D' - F2, 'fro')) ...
                /max(norm(E2, 'fro'), norm(F2,'fro')) < tol, sprintf('[%d-3] E2 - F2 failed', i) );
        assert(max( norm ( A' * R3 +  optset.sign * L3 * B  - E3, 'fro'), norm( C' * R3 +  optset.sign2 * L3 * D  - F3, 'fro')) ...
                /max(norm(E3, 'fro'), norm(F3,'fro')) < tol, sprintf('[%d-2] E3 - F3 failed', i) );
        assert(max( norm ( A' * R4 +  optset.sign * L4 * B' - E4, 'fro'), norm( C' * R4 +  optset.sign2 * L4 * D' - F4, 'fro')) ...
                /max(norm(E4, 'fro'), norm(F4,'fro')) < tol, sprintf('[%d-2] E4 - F4 failed', i) );
    end
end
