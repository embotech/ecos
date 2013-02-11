function a = conelp_linesearch(s,z,tau,kap,ds,dz,dtau,dkap,dims,lambda)
% Backtracking linesearch for conelp.

%% parameters
amin = 1e-8;
amax = 10;
% nbisect = 10;

%% LPcone
% au = amax;
% al = amin;
afact = 0.9;
a = amax;
if( dims.l > 0 )
    sl = s(1:dims.l); dsl = ds(1:dims.l);
    zl = z(1:dims.l); dzl = dz(1:dims.l);
    snew = sl + a*dsl;
    znew = zl + a*dzl;
    taunew = tau + a*dtau;
    kapnew = kap + a*dkap;
    while( any( snew <= 0 ) || any( znew <= 0 ) || taunew <= 0 || kapnew <= 0 )
        a = a*afact;
        assert( a > amin, 'No further progress possible along search direction, exiting.' );
        snew = sl + a*dsl;
        znew = zl + a*dzl;
        taunew = tau + a*dtau;
        kapnew = kap + a*dkap;
    end
    
    %     for i = 1:dims.l
    %         if( dsl(i) < 0 ), a = min([a, -sl(i)/dsl(i)-1e-6]); end
    %         if( dzl(i) < 0 ), a = min([a, -zl(i)/dzl(i)-1e-6]); end
    %     end
    
    %     % bisection on s and l
    %     for i = 1:nbisect
    %         ac = al + (au - al)/2;
    %         if( any( sl + ac*dsl < 0 ) || any( zl + ac*dzl < 0 ) || (tau+ac*dtau < 0) || (kap+ac*dkap <0) )
    %             au = ac;
    %         else
    %             al = ac;
    %         end
    %     end
    %
    %     a = min([a, al]);
    
    
    %     rhok = ds(1:dims.l) ./ lambda(1:dims.l);
    %     sigmak = dz(1:dims.l) ./ lambda(1:dims.l);
    %     a = min([a, 1/max([0, -min(rhok), -min(sigmak), -dtau/tau, -dkap/kap])]);
    %alpha = min([alpha, 1/max([0, -min(rhok), -min(sigmak)])]);
    
    assert( a >= amin, 'step length too short (LP cone)');
    assert( a <= amax, 'step length too long  (LP-cone)');
    assert( all(sl + a*dsl >= 0), 's leaving LP cone');
    assert( all(zl + a*dzl >= 0), 'z leaving LP cone');
end


%% second-order cone
if( ~isempty(dims.q) )
    for k = 1:length(dims.q)
        coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
        
        zk = z(coneidx);   dzk = dz(coneidx);
        sk = s(coneidx);   dsk = ds(coneidx);
        
        %zknew = zk + a*dzk;    sknew = sk + a*dsk;
        
        a = conelp_linesearch_soc(a,sk,dsk);
        a = conelp_linesearch_soc(a,zk,dzk);
        
        %         zres = zknew(1)^2 - zknew(2:end)'*zknew(2:end);
        %         sres = sknew(1)^2 - sknew(2:end)'*sknew(2:end);
        %         while( any(zres <= 0) || any(sres <= 0) )
        %             a = a*afact;
        %             assert( a > amin, 'No further progress possible along search direction, exiting.' );
        %             zknew = zk + a*dzk; sknew = sk + a*dsk;
        %             zres = zknew(1)^2 - zknew(2:end)'*zknew(2:end);
        %             sres = sknew(1)^2 - sknew(2:end)'*sknew(2:end);
        %         end
    end
end