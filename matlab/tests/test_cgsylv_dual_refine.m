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

function test_csylv_dual_refine()
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
    m = 199;

    tol = sqrt(eps()) * max(m,n);
    ierr = 0;

    A = rand( m , m );
    C = rand( m,  m );
    B = rand( n,  n );
    D = rand( n,  n );

    if (exist('OCTAVE_VERSION', 'builtin'))
        [AS, CS, Q, Z] = qz(A, C);
        [BS, DS, U, V] = qz(B, D);
    else
        [AS, CS, Q, Z] = qz(A, C, 'real');
        [BS, DS, U, V] = qz(B, D, 'real');
    end
    Q = Q';
    U = U';

    R = ones ( m , n );
    L = 2*ones ( m , n );
    R0 = R + rand( m, n);
    L0 = L + rand( m, n);

    problems = cell(16,1);
    problems{1} = {'N', 'N', 1, 1};
    problems{2} = {'N', 'T', 1, 1};
    problems{3} = {'T', 'N', 1, 1};
    problems{4} = {'N', 'N', 1, 1};
    problems{5} = {'N', 'N', 1, -1};
    problems{6} = {'N', 'T', 1, -1};
    problems{7} = {'T', 'N', 1, -1};
    problems{8} = {'N', 'N', 1, -1};
    problems{9} =  {'N', 'N', -1, 1};
    problems{10} = {'N', 'T', -1, 1};
    problems{11} = {'T', 'N', -1, 1};
    problems{12} = {'N', 'N', -1, 1};
    problems{13} = {'N', 'N', -1, -1};
    problems{14} = {'N', 'T', -1, -1};
    problems{15} = {'T', 'N', -1, -1};
    problems{16} = {'N', 'N', -1, -1};

    for p = 1:max(size(problems))

        prob = problems{p};

        E = op(A, prob{1})' * R + op(C, prob{1})' * L;
        F = prob{3} * R * op(B, prob{2})' + prob{4} * L * op(D, prob{2});

        if (prob{1} == 'N' && prob{2} == 'N')
            tsp = '';
        else
            tsp = sprintf(', ''%s'', ''%s''', prob{1}, prob{2});
        end
        optset.sign = prob{3};
        optset.sign2 = prob{4};

        for k = 1:2
            if k == 1
                convlog = '';
            else
                convlog = ', convlog';
            end

            if prob{3} == 1 && prob{4} == 1
                t1 = sprintf('[R1, L1 %s] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F %s);', convlog, tsp);
                t2 = sprintf('[R2, L2 %s] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0 %s);', convlog, tsp);
                t3 = sprintf('[R3, L3 %s] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F %s , optset );', convlog, tsp);
                t4 = sprintf('[R4, L4 %s] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0 %s , optset );', convlog, tsp);

                eval(t1);
                eval(t2);
                eval(t3);
                eval(t4);

                assert(max( norm ( op(A, prob{1})' * R1           +  op(C, prob{1})' * L1           - E, 'fro')/norm(E, 'fro'), ...
                            norm ( optset.sign * R1 * op(B, prob{2})' + optset.sign2 * L1 * op(D, prob{2})' - F, 'fro')/norm(F, 'fro')) ...
                            < tol,   '[ %d-%d-1 ] %s', p, k, t1 );
                assert(max( norm ( op(A, prob{1})' * R2           +   op(C, prob{1})' * L2           - E, 'fro')/norm(E, 'fro'), ...
                            norm ( optset.sign * R2 * op(B, prob{2})' + optset.sign2 * L2 * op(D, prob{2})' - F, 'fro')/norm(F, 'fro')) ...
                            < tol,   '[ %d-%d-2 ] %s', p, k, t1 );
                assert(max( norm ( op(A, prob{1})' * R3           +  op(C, prob{1})' * L3           - E, 'fro')/norm(E, 'fro'), ...
                            norm ( optset.sign * R3 * op(B, prob{2})' + optset.sign2 * L3 * op(D, prob{2})' - F, 'fro')/norm(F, 'fro')) ...
                            < tol,   '[ %d-%d-3 ] %s', p, k, t1 );
                assert(max( norm ( op(A, prob{1})' * R4           +  op(C, prob{1})' * L4           - E, 'fro')/norm(E, 'fro'), ...
                            norm ( optset.sign * R4 * op(B, prob{2})' + optset.sign2 * L4 * op(D, prob{2})' - F, 'fro')/norm(F, 'fro')) ...
                            < tol,   '[ %d-%d-4 ] %s', p, k, t1 );


            else
                t1 = sprintf('[R1, L1 %s] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F %s , optset );', convlog, tsp);
                t2 = sprintf('[R2, L2 %s] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0 %s , optset );', convlog, tsp);

                eval(t1);
                eval(t2);

                assert(max( norm ( op(A, prob{1})' * R1           +  op(C, prob{1})' * L1           - E, 'fro')/norm(E, 'fro'), ...
                            norm ( optset.sign * R1 * op(B, prob{2})' + optset.sign2 * L1 * op(D, prob{2})' - F, 'fro')/norm(F, 'fro')) ...
                            < tol,   '[ %d-%d-1 ] %s', p, k, t1 );
                assert(max( norm ( op(A, prob{1})' * R2           +  op(C, prob{1})' * L2           - E, 'fro')/norm(E, 'fro'), ...
                            norm ( optset.sign * R2 * op(B, prob{2})' + optset.sign2 * L2 * op(D, prob{2})' - F, 'fro')/norm(F, 'fro')) ...
                            < tol,   '[ %d-%d-2 ] %s', p, k, t1 );

            end
        end
    end
end

function X = op(A, op)
     mepack_test_init_random();
    if (op == 'N')
        X = A;
    elseif ( op == 'T' )
        X = A';
    end
end
