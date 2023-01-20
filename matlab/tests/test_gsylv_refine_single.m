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

function test_gsylv_refine_single()
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

    tol = sqrt(eps('single')) * max(m,n);
    ierr = 0;

    A = single(rand( m , m ));
    C = single(rand( m,  m ));
    B = single(rand( n,  n ));
    D = single(rand( n,  n ));

    if (exist('OCTAVE_VERSION', 'builtin'))
        [AS, CS, Q, Z] = qz(A, C);
        [BS, DS, U, V] = qz(B, D);
    else
        [AS, CS, Q, Z] = qz(A, C, 'real');
        [BS, DS, U, V] = qz(B, D, 'real');
    end
    Q = single(Q)';
    U = single(U)';
    Z = single(Z);
    V = single(V);
    AS = single(AS);
    BS = single(BS);
    CS = single(CS);
    DS = single(DS);

    X = single(ones ( m , n ));
    X0 = X + single(rand( m, n));

    problems = cell(8,1);
    problems{1} = {'N', 'N', 1};
    problems{2} = {'N', 'T', 1};
    problems{3} = {'T', 'N', 1};
    problems{4} = {'N', 'N', 1};
    problems{5} = {'N', 'N', -1};
    problems{6} = {'N', 'T', -1};
    problems{7} = {'T', 'N', -1};
    problems{8} = {'N', 'N', -1};

    for p = 1:max(size(problems))

        prob = problems{p};

        Y = op(A, prob{1}) * X * op(B, prob{2} )+ prob{3} * op(C, prob{1}) *X * op(D, prob{2});

        if (prob{1} == 'N' && prob{2} == 'N')
            tsp = '';
        else
            tsp = sprintf(', ''%s'', ''%s''', prob{1}, prob{2});
        end
        optset.sign = prob{3};

        for k = 1:2
            if k == 1
                convlog = '';
            else
                convlog = ', convlog';
            end

            if prob{3} == 1
                t1 = sprintf('[X1 %s] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y %s);', convlog, tsp);
                t2 = sprintf('[X2 %s] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0 %s);', convlog, tsp);
                t3 = sprintf('[X3 %s] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y %s , optset );', convlog, tsp);
                t4 = sprintf('[X4 %s] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0  %s , optset );', convlog, tsp);

                eval(t1);
                eval(t2);
                eval(t3);
                eval(t4);

                assert(norm(op(A, prob{1}) * X1 * op(B, prob{2}) + prob{3} * op(C, prob{1}) * X1 * op(D, prob{2}) - Y, 'fro')/norm(Y, 'fro') < tol,   '[ %d-%d-1 ] %s', p, k, t1 );
                assert(norm(op(A, prob{1}) * X2 * op(B, prob{2}) + prob{3} * op(C, prob{1}) * X2 * op(D, prob{2}) - Y, 'fro')/norm(Y, 'fro') < tol,   '[ %d-%d-2 ] %s', p, k, t2 );
                assert(norm(op(A, prob{1}) * X3 * op(B, prob{2}) + prob{3} * op(C, prob{1}) * X3 * op(D, prob{2}) - Y, 'fro')/norm(Y, 'fro') < tol,   '[ %d-%d-3 ] %s', p, k, t3 );
                assert(norm(op(A, prob{1}) * X4 * op(B, prob{2}) + prob{3} * op(C, prob{1}) * X4 * op(D, prob{2}) - Y, 'fro')/norm(Y, 'fro') < tol,   '[ %d-%d-4 ] %s', p, k, t4 );

            else

                t1 = sprintf('[X1 %s] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y %s , optset );', convlog, tsp);
                t2 = sprintf('[X2 %s] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0  %s , optset );', convlog, tsp);

                eval(t1);
                eval(t2);

                assert(norm(op(A, prob{1}) * X1 * op(B, prob{2}) + prob{3} * op(C, prob{1}) * X1 * op(D, prob{2}) - Y, 'fro')/norm(Y, 'fro') < tol,   '[ %d-%d-1 ] %s', p, k, t1 );
                assert(norm(op(A, prob{1}) * X2 * op(B, prob{2}) + prob{3} * op(C, prob{1}) * X2 * op(D, prob{2}) - Y, 'fro')/norm(Y, 'fro') < tol,   '[ %d-%d-2 ] %s', p, k, t2 );

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
