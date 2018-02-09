N = 21;
L = 2;
alpha = .25;
Tscale = 1;

[ rV, d0, c2v, bndryVs, bndryCs, T_bc, rB_bc ] = generate.hexLattice( N, L, alpha );
T_bc = Tscale * T_bc;
Fbc = bsxfun(@times,0*T_bc,rB_bc);

%%
% alphaVec = linspace(1,5,8);
alphaVec = 1;
% gammaVec = linspace(0,40,8);
gammaVec = 40;
% betaVec = linspace(0,.35,20);
betaVec = .35;

nIter = 1;
for n = 1:nIter
    n
    for ii = 1:length(gammaVec)
        for jj = 1:length(betaVec)
            for kk = 1:length(alphaVec)
                beta = betaVec(jj);
                a0 = alphaVec(kk)*abs((3*sqrt(3)/2 * L^2)*(1 + beta*randn(size(c2v,1),1)));
                gamma = gammaVec(ii);
                Lc = 0;
                Tpull = 500;
                Pb = gamma*(mean(a0) - (3*sqrt(3)/2 * L^2)) - 5;

                vMsim = simP.vertexModel( rV, bndryVs, bndryCs, d0, c2v, gamma, a0, Pb, Fbc, Lc );
                vMsim = vMsim.evolve(Tpull);
                [ kappa ] = vMsim.returnEndCurvature();
                k(ii,jj,n) = median(kappa);
                
                Struct = vMsim.returnStruct();
%                 [ chi, chiEst, p ] = checkPressureCompt( Struct );
                [ chi, chiC, chiCV ] = isogonal.buildControlDistribution( Struct, 0 );
%                 c(ii,jj,n) = corr(chi,chiEst);
            end
        end
    end
end

% ind = abs(chiCV{1}) < .25;
chiT = {chiCV{1},chi{1}};
plot.niceHist(chiT,30)
    
