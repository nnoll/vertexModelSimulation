N = 10;
L = sqrt(3);
alpha = .13;
Tscale = .5;

if (~exist('rV'))
    [ rV, d0, c2v, bndryVs, T, T_bc, rB_bc, Tri, q ] = generate.rectLattice( N, L, alpha );
    T = Tscale * T;
    m = T;
    T_bc = Tscale * T_bc;
    Fbc = bsxfun(@times,T_bc,rB_bc);

end

%%
nu = .1;
omega = 10;
Lc = .01;
mode = 'shear';
F = 0;
P = Tscale/sqrt(3);
% Fbc = zeros(size(Fbc));
deltaVec = 0;
clear pR meanTheta
for ii = 1:length(deltaVec)
    
    delta = deltaVec(ii);
    Tpull = delta + 20000;

    Fext = @( t, nB, vp, vm, dir ) sim.returnLeadingExt( t, F, delta, nB, vp, vm, dir );
    vMsim = sim.vertexModel( rV, bndryVs, q, d0, c2v, Tri, T, m, nu, omega, Fbc, Fext, P, Lc, mode );

    vMsim = vMsim.evolve(Tpull);
    
end
