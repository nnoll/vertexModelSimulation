N = 12;
L = sqrt(3);
% alpha = .07;
alpha = 0;
Tscale = .2;
deltaT = .3;

if (~exist('rV'))
    [ rV, d0, c2v, bndryVs, bndryEs, T, T_bc, rB_bc, Tri, q ] = generate.lattice( N, L, alpha );
    rB = .5*abs(d0)*rV;
    [~,ind] = min(sum(rB.^2,2));
%     T(ind) = 1 + deltaT;
    T = Tscale * T;
    T(bndryEs) = .5*T(bndryEs);
    m = T;
    T(ind) = T(ind) + deltaT;
    T_bc = Tscale * T_bc;
    Fbc = bsxfun(@times,T_bc,rB_bc);

end

%%
nu = .1;
omega = .09;
Lc = .01;
mode = 'shear';
F = 0;
P = Tscale/sqrt(3);
Fbc = zeros(size(Fbc));
deltaVec = 0;
clear pR meanTheta
for ii = 1:length(deltaVec)
    
    delta = deltaVec(ii);
    Tpull = delta + 5000;

    Fext = @( t, nB, vp, vm, dir ) sim.returnLeadingExt( t, F, delta, nB, vp, vm, dir );
    vMsim = sim.vertexModel( rV, bndryVs, q, d0, c2v, Tri, T, m, nu, omega, Fbc, Fext, P, Lc, mode );

    vMsim = vMsim.evolve(Tpull);
%     mod(ii) = vMsim.returnStrain();

%     Fext = @( t, nB, vp, vm, dir ) sim.returnEndingExt( t, Tpull, F, delta, nB, vp, vm, dir );
% 
%     vMsim = vMsim.updateFext(Fext);
%     vMsim = vMsim.evolve(Tpull);
    
%     [theta] = vMsim.fitIsogonal();
%     meanTheta(ii) = mean(theta-min(theta));
    
end

% clearvars -except vMsim delta F Tpull ii deltaVec mod meanTheta omega nu

% vMsim.compareInitialtoFinal();
