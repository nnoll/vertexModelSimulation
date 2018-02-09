N = 10;
L = sqrt(3);
alpha = 0;
Tscale = .5;

if (~exist('rV'))
    [ rV, d0, c2v, bndryVs, T, T_bc, rB_bc, Tri, q ] = generate.rectLattice( N, L, alpha );
    T = Tscale * T;
    m = T;
    T_bc = Tscale * T_bc;
    Fbc = bsxfun(@times,T_bc,rB_bc);

end

%%
Ltissue = max(rV(:,2)) - min(rV(:,2));
rV = rV/Ltissue;
T = T/Ltissue;
m = m/Ltissue;
Fbc = Fbc/Ltissue;

nu = .01;
alpha = .0001;
Lc = 0;
mode = 'pull';
F = .0002/Ltissue;
P = 0;
omegaVec = logspace(-5,0,100);
% omegaVec = .000001;

clear pR meanTheta
for ii = 1:length(omegaVec)
    ii
    omega = omegaVec(ii);
    f0 = omega/(2*pi);
    deltaT = 1/(100*f0);
    Tpull = max(round(10/f0),1);
    Fext = @( t, nB, vp, vm, dir ) sim.returnSinusoidExt( t, F, omega, nB, vp, vm, dir );
    
    vMsim = sim.vertexModel( rV, bndryVs, q, d0, c2v, Tri, T, m, nu, alpha, Fbc, Fext, P, Lc, mode );

    vMsim = vMsim.evolve(Tpull,deltaT);
    [phase,mag] = vMsim.computeFFT_rB(f0,deltaT);
    [phaseV(:,ii),magV(:,ii),y] = vMsim.produceAveragedFcn(phase,mag);

end
