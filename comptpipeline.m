N = 7;
L = 2;
alpha = 0;
Tscale = 1;

[ rV, d0, c2v, bndryVs, bndryCs, T_bc, rB_bc ] = generate.hexLattice( N, L, alpha );
% [ rV, d0, c2v, bndryVs, bndryEs, bndryCs, T, T_bc, rB_bc ] = generate.hextLattice( N, L, alpha );
T_bc = 0 * Tscale * T_bc;
Fbc = bsxfun(@times,T_bc,rB_bc);

%%
% gammaVec = linspace(.05,5,5);
gammaVec = 3;
beta = .3;

nIter = 1;
for n = 1:nIter
    n
    for ii = 1:length(gammaVec)
            
        a0 = abs((3*sqrt(3)/2 * L^2)*(1 + beta*randn(size(c2v,1),1)));
        gamma = gammaVec(ii);
        Lc = 0;
        Tpull = 100;
        Pb = gamma*(mean(a0) - (3*sqrt(3)/2 * L^2)) - .15;
        l0 = abs(ones(size(d0,1),1) + .5*randn(size(d0,1),1));
        kappaL = abs(1 + 0*randn(size(l0)));
        vMsim = simP.vertexModel( rV, bndryVs, bndryCs, d0, c2v, gamma, a0, Pb, Fbc, l0, kappaL, Lc );
        vMsim = vMsim.evolve(Tpull);
%         [ kappa ] = vMsim.returnEndCurvature();
%         k(ii,n) = median(kappa);
% 
        Struct = vMsim.returnStruct();
%         R = median([Struct.Bdat.length]);
%         for jj = 1:length(alphaVec)
%             Struct = vMsim.returnStruct(alphaVec(jj)/R);
% %             [ chi, chiC, chiCV ] = isogonal.buildControlDistribution( Struct, 0 );
% %             c(ii,jj,n) = std(chi{1}(~isnan(chi{1})));
% %             cV(ii,jj,n) = std(chiCV{1}(~isnan(chiCV{1})));
% %             cC(ii,jj,n) = std(chiC{1}(~isnan(chiC{1})));
%         end
            
    end
end
    
