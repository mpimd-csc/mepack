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
function test_tgcsylv()
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
    m = 198;

    tol = sqrt(eps) * max(m,n);
    tols = sqrt(eps('single')) * max(m,n);
    ierr = 0;

    A = rand( m , m );
    C = rand( m , m );
    B = rand( n , n );
    D = rand( n , n );

    is_octave =  exist('OCTAVE_VERSION', 'builtin');

    if (is_octave)
        [A, C] = qz(A, C);
        [B, D] = qz(B, D);
    else
        [A, C] = qz(A, C, 'real');
        [B, D] = qz(B, D, 'real');
    end



    R = ones ( m , n );
    L = ones ( m , n );

    for i = 1:4
        E1 = A*R + L*B;
        F1 = C*R + L*D;
        E2 = A*R + L*B';
        F2 = C*R + L*D';
        E3 = A'*R + L*B;
        F3 = C'*R + L*D;
        E4 = A'*R + L*B';
        F4 = C'*R + L*D';

        E5 = A*R  - L*B;
        F5 = C*R  + L*D;
        E6 = A*R  - L*B';
        F6 = C*R  + L*D';
        E7 = A'*R - L*B;
        F7 = C'*R + L*D;
        E8 = A'*R - L*B';
        F8 = C'*R + L*D';

        E9 = A*R   + L*B;
        F9 = C*R   - L*D;
        E10 = A*R  + L*B';
        F10 = C*R  - L*D';
        E11 = A'*R + L*B;
        F11 = C'*R - L*D;
        E12 = A'*R + L*B';
        F12 = C'*R - L*D';

        E13 = A*R  - L*B;
        F13 = C*R  - L*D;
        E14 = A*R  - L*B';
        F14 = C*R  - L*D';
        E15 = A'*R - L*B;
        F15 = C'*R - L*D;
        E16 = A'*R - L*B';
        F16 = C'*R - L*D';



        if i >= 3
            A = single(A);
            B = single(B);
            C = single(C);
            D = single(D);
            R = single(R);
            L = single(L);

            E1 = single ( E1 );
            E2 = single ( E2 );
            E3 = single ( E3 );
            E4 = single ( E4 );
            E5 = single ( E5 );
            E6 = single ( E6 );
            E7 = single ( E7 );
            E8 = single ( E8 );
            E9 = single ( E9 );
            E10 = single ( E10 );
            E11 = single ( E11 );
            E12 = single ( E12 );
            E13 = single ( E13 );
            E14 = single ( E14 );
            E15 = single ( E15 );
            E16 = single ( E16 );

            F1 = single ( F1 );
            F2 = single ( F2 );
            F3 = single ( F3 );
            F4 = single ( F4 );
            F5 = single ( F5 );
            F6 = single ( F6 );
            F7 = single ( F7 );
            F8 = single ( F8 );
            F9 = single ( F9 );
            F10 = single ( F10 );
            F11 = single ( F11 );
            F12 = single ( F12 );
            F13 = single ( F13 );
            F14 = single ( F14 );
            F15 = single ( F15 );
            F16 = single ( F16 );

            tol = tols;
        end

        if i == 2 || i == 4
            optset.openmp = 1;
        else
            optset.openmp = 0;
        end

        [R1, L1] = mepack_tgcsylv(A, B, C, D, E1, F1);
        [R2, L2] = mepack_tgcsylv(A, B, C, D, E2, F2, 'N', 'T');
        [R3, L3] = mepack_tgcsylv(A, B, C, D, E3, F3, 'T', 'N');
        [R4, L4] = mepack_tgcsylv(A, B, C, D, E4, F4, 'T', 'T');

        assert(max( norm ( A * R1  +  L1 * B  - E1, 'fro'), norm( C * R1  +  L1 * D  - F1, 'fro'))/max(norm(E1, 'fro'), norm(F1,'fro')) < tol, ' E1 - F1 failed' );
        assert(max( norm ( A * R2  +  L2 * B' - E2, 'fro'), norm( C * R2  +  L2 * D' - F2, 'fro'))/max(norm(E2, 'fro'), norm(F2,'fro')) < tol, ' E2 - F2 failed' );
        assert(max( norm ( A' * R3 +  L3 * B  - E3, 'fro'), norm( C' * R3 +  L3 * D  - F3, 'fro'))/max(norm(E3, 'fro'), norm(F3,'fro')) < tol, ' E3 - F3 failed' );
        assert(max( norm ( A' * R4 +  L4 * B' - E4, 'fro'), norm( C' * R4 +  L4 * D' - F4, 'fro'))/max(norm(E4, 'fro'), norm(F4,'fro')) < tol, ' E4 - F4 failed' );

        ox = optset;
        ox.sign = -1;
        ox.sign2 = 1;

        [R1, L1] = mepack_tgcsylv(A, B, C, D, E5, F5, ox);
        [R2, L2] = mepack_tgcsylv(A, B, C, D, E6, F6, 'N', 'T', ox);
        [R3, L3] = mepack_tgcsylv(A, B, C, D, E7, F7, 'T', 'N', ox);
        [R4, L4] = mepack_tgcsylv(A, B, C, D, E8, F8, 'T', 'T', ox);

        assert(max( norm ( A * R1   - L1 * B  - E5, 'fro'), norm( C * R1  + L1 * D  - F5, 'fro'))/max(norm(E5, 'fro'), norm(F5,'fro')) < tol, sprintf(' E5 - F5 failed %d', i) );
        assert(max( norm ( A * R2   - L2 * B' - E6, 'fro'), norm( C * R2  + L2 * D' - F6, 'fro'))/max(norm(E6, 'fro'), norm(F6,'fro')) < tol, sprintf(' E6 - F6 failed %d', i) );
        assert(max( norm ( A' * R3  - L3 * B  - E7, 'fro'), norm( C' * R3 + L3 * D  - F7, 'fro'))/max(norm(E7, 'fro'), norm(F7,'fro')) < tol, sprintf(' E7 - F7 failed %d', i) );
        assert(max( norm ( A' * R4  - L4 * B' - E8, 'fro'), norm( C' * R4 + L4 * D' - F8, 'fro'))/max(norm(E8, 'fro'), norm(F8,'fro')) < tol, sprintf(' E8 - F8 failed %d', i) );

        ox = optset;
        ox.sign = 1;
        ox.sign2 = -1;
        [R1, L1] = mepack_tgcsylv(A, B, C, D, E9, F9, ox);
        [R2, L2] = mepack_tgcsylv(A, B, C, D, E10, F10, 'N', 'T', ox);
        [R3, L3] = mepack_tgcsylv(A, B, C, D, E11, F11, 'T', 'N', ox);
        [R4, L4] = mepack_tgcsylv(A, B, C, D, E12, F12, 'T', 'T', ox);

        assert(max( norm ( A * R1   + L1 * B  - E9, 'fro'),  norm( C * R1  - L1 * D  - F9, 'fro')) /max(norm(E9, 'fro'),  norm(F9,'fro')) < tol,  ' E9  - F9  failed' );
        assert(max( norm ( A * R2   + L2 * B' - E10, 'fro'), norm( C * R2  - L2 * D' - F10, 'fro'))/max(norm(E10, 'fro'), norm(F10,'fro')) < tol, ' E10 - F10 failed' );
        assert(max( norm ( A' * R3  + L3 * B  - E11, 'fro'), norm( C' * R3 - L3 * D  - F11, 'fro'))/max(norm(E11, 'fro'), norm(F11,'fro')) < tol, ' E11 - F11 failed' );
        assert(max( norm ( A' * R4  + L4 * B' - E12, 'fro'), norm( C' * R4 - L4 * D' - F12, 'fro'))/max(norm(E12, 'fro'), norm(F12,'fro')) < tol, ' E12 - F12 failed' );

        ox = optset;
        ox.sign = -1;
        ox.sign2 = -1;
        [R1, L1] = mepack_tgcsylv(A, B, C, D, E13, F13, ox);
        [R2, L2] = mepack_tgcsylv(A, B, C, D, E14, F14, 'N', 'T', ox);
        [R3, L3] = mepack_tgcsylv(A, B, C, D, E15, F15, 'T', 'N', ox);
        [R4, L4] = mepack_tgcsylv(A, B, C, D, E16, F16, 'T', 'T', ox);

        assert(max( norm ( A * R1  -  L1 * B  - E13, 'fro'), norm( C * R1  - L1 * D  - F13, 'fro'))/max(norm(E13, 'fro'), norm(F13,'fro')) < tol, ' E13 - F13  failed' );
        assert(max( norm ( A * R2  -  L2 * B' - E14, 'fro'), norm( C * R2  - L2 * D' - F14, 'fro'))/max(norm(E14, 'fro'), norm(F14,'fro')) < tol, ' E14 - F14 failed' );
        assert(max( norm ( A' * R3 -  L3 * B  - E15, 'fro'), norm( C' * R3 - L3 * D  - F15, 'fro'))/max(norm(E15, 'fro'), norm(F15,'fro')) < tol, ' E15 - F15 failed' );
        assert(max( norm ( A' * R4 -  L4 * B' - E16, 'fro'), norm( C' * R4 - L4 * D' - F16, 'fro'))/max(norm(E16, 'fro'), norm(F16,'fro')) < tol, ' E16 - F16 failed' );

    end
end
