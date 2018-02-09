N = 9;
L = 2;
alpha = .25;
Tscale = 1;

[ rV, d0, c2v, bndryVs, bndryCs, T_bc, rB_bc ] = generate.hexLattice( N, L, alpha );
T_bc = Tscale * T_bc;
Fbc = bsxfun(@times,0*T_bc,rB_bc);

%%
% alphaVec = linspace(1,5,8);
% alphaVec = linspace(0,1,4);
alphaVec = 1;
% gammaVec = linspace(0,3,4);
gammaVec = .01;
% betaVec = linspace(.05,.4,4);
betaVec = .05;

nIter = 5;
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
                k(ii,jj,kk,n) = median(kappa);
                
                Struct = vMsim.returnStruct();
                [ chi, chiC, chiCV ] = isogonal.buildControlDistribution( Struct, 0 );
                sigma(ii,jj,kk,n) = std(chi{1});
                sigmaC(ii,jj,kk,n) = std(chiC{1});
                sigmaCV(ii,jj,kk,n) = std(chiCV{1});
            end
        end
    end
end

k = squeeze(k);
imagesc(gammaVec,betaVec,median(k,3))
    
