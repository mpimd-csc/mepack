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
function test_tgcsylv_dual()
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
        E1 = A'*R + C'*L;
        F1 = R*B' + L*D';
        E2 = A'*R + C'*L;
        F2 = R*B  + L*D;
        E3 = A*R  + C*L;
        F3 = R*B' + L*D';
        E4 = A*R  + C*L;
        F4 = R*B  + L*D;

        E5 = A'*R + C'*L;
        F5 = -R*B' + L*D';
        E6 = A'*R + C'*L;
        F6 = -R*B  + L*D;
        E7 = A*R  + C*L;
        F7 = -R*B' + L*D';
        E8 = A*R  + C*L;
        F8 = -R*B  + L*D;

        E9  = A'*R + C'*L;
        F9  = R*B' - L*D';
        E10 = A'*R + C'*L;
        F10 = R*B  - L*D;
        E11 = A*R  + C*L;
        F11 = R*B' - L*D';
        E12 = A*R  + C*L;
        F12 = R*B  - L*D;

        E13 = A'*R + C'*L;
        F13 = - R*B' - L*D';
        E14 = A'*R + C'*L;
        F14 = - R*B  - L*D;
        E15 = A*R  + C*L;
        F15 = - R*B' - L*D';
        E16 = A*R  + C*L;
        F16 = - R*B  - L*D;



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

        [R1, L1] = mepack_tgcsylv_dual(A, B, C, D, E1, F1);
        [R2, L2] = mepack_tgcsylv_dual(A, B, C, D, E2, F2, 'N', 'T');
        [R3, L3] = mepack_tgcsylv_dual(A, B, C, D, E3, F3, 'T', 'N');
        [R4, L4] = mepack_tgcsylv_dual(A, B, C, D, E4, F4, 'T', 'T');

        assert(max( norm ( A' * R1  +  C' * L1  - E1, 'fro'), norm( R1 * B' +  L1 * D'  - F1, 'fro'))/max(norm(E1, 'fro'), norm(F1,'fro')) < tol, ' E1 - F1 failed' );
        assert(max( norm ( A' * R2  +  C' * L2  - E2, 'fro'), norm( R2 * B  +  L2 * D   - F2, 'fro'))/max(norm(E2, 'fro'), norm(F2,'fro')) < tol, ' E2 - F2 failed' );
        assert(max( norm ( A  * R3  +  C  * L3  - E3, 'fro'), norm( R3 * B' +  L3 * D'  - F3, 'fro'))/max(norm(E3, 'fro'), norm(F3,'fro')) < tol, ' E3 - F3 failed' );
        assert(max( norm ( A  * R4  +  C  * L4  - E4, 'fro'), norm( R4 * B  +  L4 * D   - F4, 'fro'))/max(norm(E4, 'fro'), norm(F4,'fro')) < tol, ' E4 - F4 failed' );

        ox = optset;
        ox.sign = -1;
        ox.sign2 = 1;
        [R1, L1] = mepack_tgcsylv_dual(A, B, C, D, E5, F5, ox);
        [R2, L2] = mepack_tgcsylv_dual(A, B, C, D, E6, F6, 'N', 'T', ox);
        [R3, L3] = mepack_tgcsylv_dual(A, B, C, D, E7, F7, 'T', 'N', ox);
        [R4, L4] = mepack_tgcsylv_dual(A, B, C, D, E8, F8, 'T', 'T', ox);

        assert(max( norm ( A' * R1  +  C' * L1  - E5, 'fro'), norm( - R1 * B' +  L1 * D'  - F5, 'fro'))/max(norm(E5, 'fro'), norm(F5,'fro')) < tol, ' E5 - F5 failed' );
        assert(max( norm ( A' * R2  +  C' * L2  - E6, 'fro'), norm( - R2 * B  +  L2 * D   - F6, 'fro'))/max(norm(E6, 'fro'), norm(F6,'fro')) < tol, ' E6 - F6 failed' );
        assert(max( norm ( A  * R3  +  C  * L3  - E7, 'fro'), norm( - R3 * B' +  L3 * D'  - F7, 'fro'))/max(norm(E7, 'fro'), norm(F7,'fro')) < tol, ' E7 - F7 failed' );
        assert(max( norm ( A  * R4  +  C  * L4  - E8, 'fro'), norm( - R4 * B  +  L4 * D   - F8, 'fro'))/max(norm(E8, 'fro'), norm(F8,'fro')) < tol, ' E8 - F8 failed' );


        ox = optset;
        ox.sign = 1;
        ox.sign2 = -1;
        [R1, L1] = mepack_tgcsylv_dual(A, B, C, D, E9, F9, ox);
        [R2, L2] = mepack_tgcsylv_dual(A, B, C, D, E10, F10, 'N', 'T', ox);
        [R3, L3] = mepack_tgcsylv_dual(A, B, C, D, E11, F11, 'T', 'N', ox);
        [R4, L4] = mepack_tgcsylv_dual(A, B, C, D, E12, F12, 'T', 'T', ox);

        assert(max( norm ( A' * R1  +  C' * L1  - E9, 'fro'),  norm(  R1 * B' -  L1 * D'  - F9, 'fro'))/ max(norm(E9, 'fro'),  norm(F9,'fro')) < tol , ' E9  - F9 failed' );
        assert(max( norm ( A' * R2  +  C' * L2  - E10, 'fro'), norm(  R2 * B  -  L2 * D   - F10, 'fro'))/max(norm(E10, 'fro'), norm(F10,'fro')) < tol, ' E10 - F10 failed' );
        assert(max( norm ( A  * R3  +  C  * L3  - E11, 'fro'), norm(  R3 * B' -  L3 * D'  - F11, 'fro'))/max(norm(E11, 'fro'), norm(F11,'fro')) < tol, ' E11 - F11 failed' );
        assert(max( norm ( A  * R4  +  C  * L4  - E12, 'fro'), norm(  R4 * B  -  L4 * D   - F12, 'fro'))/max(norm(E12, 'fro'), norm(F12,'fro')) < tol, ' E12 - F12 failed' );



        ox = optset;
        ox.sign = -1;
        ox.sign2 = -1;

        [R1, L1] = mepack_tgcsylv_dual(A, B, C, D, E13, F13, ox);
        [R2, L2] = mepack_tgcsylv_dual(A, B, C, D, E14, F14, 'N', 'T', ox);
        [R3, L3] = mepack_tgcsylv_dual(A, B, C, D, E15, F15, 'T', 'N', ox);
        [R4, L4] = mepack_tgcsylv_dual(A, B, C, D, E16, F16, 'T', 'T', ox);

        assert(max( norm ( A' * R1  +  C' * L1  - E13, 'fro'), norm( - R1 * B' -  L1 * D'  - F13, 'fro'))/max(norm(E13, 'fro'), norm(F13,'fro')) < tol ,' E13  -F13 failed' );
        assert(max( norm ( A' * R2  +  C' * L2  - E14, 'fro'), norm( - R2 * B  -  L2 * D   - F14, 'fro'))/max(norm(E14, 'fro'), norm(F14,'fro')) < tol, ' E14 - F14 failed' );
        assert(max( norm ( A  * R3  +  C  * L3  - E15, 'fro'), norm( - R3 * B' -  L3 * D'  - F15, 'fro'))/max(norm(E15, 'fro'), norm(F15,'fro')) < tol, ' E15 - F15 failed' );
        assert(max( norm ( A  * R4  +  C  * L4  - E16, 'fro'), norm( - R4 * B  -  L4 * D   - F16, 'fro'))/max(norm(E16, 'fro'), norm(F16,'fro')) < tol, ' E16 - F16 failed' );


    end
end
