N = 10;
L = 10;
alpha = 0;

[ rV, d0, c2v, bndryVs, T, T_bc, rB_bc, Tri, q ] = generate.lattice( N, L, alpha );

rB = d0*rV;
R = sqrt(sum(rB.^2,2));
rV = rV/mean(R);
q = q/mean(R);

n = 1;
T0 = T;
T_bc0 = T_bc;
omega = .0001;
nu = .01;
Lc = 0;
mode = 'pull';
F = .1;
    
% for Tscale = linspace(.001,1,10)
Tscale = .5;
T = Tscale * T0;
T_bc = Tscale * T_bc0;
Fbc = bsxfun(@times,T_bc,rB_bc);

% nuV = linspace(.0001,.001,3);
% delta = [1000,5000,10000];
delta = sort([0,10000,250000,50000,75000,1000000]);
for d = delta
% Fvec = linspace(.00001,.15,10);
% for F = Fvec 
    Tpull = 100000 + d;
    Fext = @( t, nB, vp, vm ) sim.returnLeadingExt( t, F, d, nB, vp, vm );
    
    vMsim = sim.vertexModel( rV, bndryVs, q, d0, c2v, Tri, T, T, nu, omega, Fbc, Fext, Lc, mode );
    vMsim = vMsim.evolve(Tpull);
    
%     [ rDiffTraj(:,n) ] = vMsim.plotVertexDisplacement();
%     
%     plot(nu*(1:size(rDiffTraj,1)),rDiffTraj(:,n))
%     hold all
    [ mod(n,:) ] = vMsim.returnStrain( );
    n = n + 1;
end
    
clearvars -except vMsim delta F Fvec Tpull ii deltaVec mod meanTheta omega nu nuV rDiffTraj Tscale

