N = 9;
L = 2;
alpha = .05;
Tscale = 1;
betaVec = linspace(0,.5,300);
gamma = .1;

[ rV, d0, c2v, bndryVs, bndryCs, T_bc, rB_bc ] = generate.hexLattice( N, L, alpha );
T_bc = Tscale * T_bc;
Fbc = bsxfun(@times,0*T_bc,rB_bc);

for ii = 1:length(betaVec)
    beta = betaVec(ii);
    a0 = abs((3*sqrt(3)/2 * L^2)*(1 + beta*randn(size(c2v,1),1)));
    Lc = 0;
    l0 = ones(size(d0,1),1);
    l0 = l0 + .25*randn(size(l0));
    % l0 = .9*l0 + exprnd(.1*ones(size(l0)));
    Tpull = 500;
    Pb = gamma*(mean(a0) - (3*sqrt(3)/2 * L^2)) - 1;

    vMsim = simP.vertexModel( rV, bndryVs, bndryCs, d0, c2v, gamma, a0, Pb, Fbc, l0, Lc );
    vMsim = vMsim.evolve(Tpull);
    [ kappa ] = vMsim.returnEndCurvature();

    Struct = vMsim.returnStruct();
    [T,P] = MI.invertMechABIC(Struct,0,1e-5);
    [Ta,Pa] = MI.returnActualMech(Struct,0);

    % mu = logspace(-7,2,10);
    % clear corT corP
    % for ii = 1:length(mu)
    %     [Tm,Pm] = MI.invertMechABIC(Struct,0,mu(ii));
    %     corT(ii) = corr(Tm,Ta);
    %     corP(ii) = corr(Pm,Pa);
    % end
    corT(ii) = corr(T,Ta);
    [~,p(ii)] = kstest((Ta-mean(Ta))/std(Ta));
end