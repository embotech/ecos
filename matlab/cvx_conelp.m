function shim = cvx_conelp( shim )

% CVX_SOLVER_SHIM	Conelp interface for CVX.
%   This procedure returns a 'shim': a structure containing the necessary
%   information CVX needs to use this solver in its modeling framework.

if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    [ ver, isoctave, fs, ps ] = cvx_version; %#ok
    temp = strfind( shim.spath, fs );
    shim.name    = 'conelp';
    shim.dualize = true;
    shim.path    = [ shim.spath(1:temp(end-1)), 'conelp', ps ];
end
if isempty( shim.error ),
    shim.check = @check;
    shim.solve = @solve;
else
    shim.check = [];
    shim.solve = [];
end
    
function found_bad = check( nonls ) %#ok
found_bad = false;

function [ x, status, tol, iters, y, z ] = solve( At, b, c, nonls, quiet, prec, settings )

n = length( c );
m = length( b );
K = struct( 'f', 0, 'l', 0, 'q', [], 's', [], 'scomplex', [], 'ycomplex', [] );
reord = struct( 'n', 0, 'r', [], 'c', [], 'v', [] );
reord = struct( 'f', reord, 'l', reord, 'a', reord, 'q', reord, 's', reord, 'h', reord );
reord.f.n = n;
for k = 1 : length( nonls ),
    temp = nonls( k ).indices;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    nnv = nn * nv;
    tt = nonls( k ).type;
    reord.f.n = reord.f.n - nnv;
    if strncmp( tt, 'i_', 2 ),
        error( 'Conelp does not support integer variables.' );
    elseif nn == 1 || isequal( tt, 'nonnegative' ),
        reord.l.r = [ reord.l.r ; temp(:) ];
        reord.l.c = [ reord.l.c ; reord.l.n + ( 1 : nnv )' ];
        reord.l.v = [ reord.l.v ; ones( nnv, 1 ) ];
        reord.l.n = reord.l.n + nnv;
    elseif isequal( tt, 'lorentz' ),
        if nn == 2,
            rr = [ temp ; temp ];
            cc = reshape( floor( 1 : 0.5 : 2 * nv + 0.5 ), 4, nv );
            vv = [1;1;-1;1]; vv = vv(:,ones(1,nv));
            reord.a.r = [ reord.a.r ; rr(:) ];
            reord.a.c = [ reord.a.c ; cc(:) + reord.a.n ];
            reord.a.v = [ reord.a.v ; vv(:) ];
            reord.a.n = reord.a.n + nnv;
        else
            temp = temp( [ end, 1 : end - 1 ], : );
            reord.q.r = [ reord.q.r ; temp(:) ];
            reord.q.c = [ reord.q.c ; reord.q.n + ( 1 : nnv )' ];
            reord.q.v = [ reord.q.v ; ones(nnv,1) ];
            reord.q.n = reord.q.n + nnv;
            K.q = [ K.q, nn * ones( 1, nv ) ];
        end
    elseif isequal( tt, 'semidefinite' ),
        error( 'Conelp does not support semidefinite programming.' );
    elseif isequal( tt, 'hermitian-semidefinite' ) && 1,
        error( 'Conelp does not support semidefinite programming.' );
    elseif isequal( tt, 'hermitian-semidefinite' ),
        error( 'Conelp does not support semidefinite programming.' );
    else
        error( 'Unsupported nonlinearity: %s', tt );
    end
end
if reord.f.n > 0,
    reord.f.r = ( 1 : n )';
    reord.f.r( [ reord.l.r ; reord.a.r ; reord.q.r ; reord.s.r ] ) = [];
    reord.f.c = ( 1 : reord.f.n )';
    reord.f.v = ones(reord.f.n,1);
end
n_d = max( m - n - reord.f.n + 1, isempty( At ) );
if n_d,
    reord.l.n = reord.l.n + n_d;
end
K.f = reord.f.n;
K.l = reord.l.n + reord.a.n;
n_out = reord.f.n;
reord.l.c = reord.l.c + n_out; n_out = n_out + reord.l.n;
reord.a.c = reord.a.c + n_out; n_out = n_out + reord.a.n;
reord.q.c = reord.q.c + n_out; n_out = n_out + reord.q.n;
reord.s.c = reord.s.c + n_out; n_out = n_out + reord.s.n;
reord = sparse( ...
    [ reord.f.r ; reord.l.r ; reord.a.r ; reord.q.r ; reord.s.r ], ...
    [ reord.f.c ; reord.l.c ; reord.a.c ; reord.q.c ; reord.s.c ], ...
    [ reord.f.v ; reord.l.v ; reord.a.v ; reord.q.v ; reord.s.v ], ...
    n, n_out );

At = reord' * At;
c  = reord' * c;
pars.free = K.f > 1 && nnz( K.q );
pars.eps     = prec(1);
pars.bigeps  = prec(3);
if quiet,
    pars.fid = 0;
end
add_row = isempty( At );
if add_row,
    K.f = K.f + 1;
    At = sparse( 1, 1, 1, n_out + 1, 1 );
    b = 1;
    c = [ 0 ; c ];
end
m = K.l + sum(K.q);
if( K.f > 0 )
    G = sparse([zeros(m,K.f), -speye(m)]);
else
    G = -speye(m);
end
h = zeros(m,1);
if( ~isreal(At) || ~isreal(c) || ~isreal(b) )
    error('Paris does not handle complex data');
end
A = At';
dims = K;
save P c G h dims A b
disp('problem data dumped to P.mat');
[ xx, yy, info ] = cvx_run_solver( @conelp, c, G, h, dims, A, b, 'xx', 'yy', 'info', settings, 5 );
if add_row,
    xx = xx(2:end);
    yy = zeros(0,1);
    At = zeros(n_out,0);
    % b  = zeros(0,1);
    c  = c(2:end);
end
if ~isfield( info, 'r0' ) && info.pinf,
    info.r0 = 0;
    info.iter = 0;
    info.numerr = 0;
end
tol = info.r0;
iters = info.iter;
xx = full( xx );
yy = full( yy );
status = '';
if info.pinf ~= 0,
    status = 'Infeasible';
    x = NaN * ones( n, 1 );
    y = yy;
    z = - real( reord * ( At * yy ) );
    if add_row, y = zeros( 0, 1 ); end
elseif info.dinf ~= 0
    status = 'Unbounded';
    y = NaN * ones( m, 1 );
    z = NaN * ones( n, 1 );
    x = real( reord * xx );
else
    x = real( reord * xx );
    y = yy;
    z = real( reord * ( c - At * yy ) );
    if add_row, y = zeros( 0, 1 ); end
end
if info.numerr == 2,
    status = 'Failed';
    if any( K.q == 2 ),
        warning( 'CVX:Conelp', cvx_error_format( 'This solver failure may possibly be due to a known bug in the Conelp solver. Try switching to SDPT3 by inserting "cvx_solver sdpt3" into your model.', ...
            [66,75], false, '' ) );
    end
else
    if isempty( status ),
        status = 'Solved';
    end
    if info.numerr == 1 && info.r0 > prec(2),
        status = [ 'Inaccurate/', status ];
    end
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
