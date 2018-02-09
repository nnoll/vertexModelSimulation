N = 13;
L = sqrt(3);
alpha = 0;
Tscale = .5;

[ rV, d0, c2v, bndryVs, bndryEs, bndryCs, T, T_bc, rB_bc, Tri, q ] = generate.lattice( N, L, alpha );

rB = .5*abs(d0)*rV;
[~,ind] = min(sum(rB.^2,2));
T = Tscale * T;
m = T;

T_bc = Tscale * T_bc;

Fbc = bsxfun(@times,T_bc,rB_bc);

tauL = .5;
alphaBar = 0;
Pb = 0;
Lc = 0;
A = 1.5*Tscale;
Tpulse = 30;
sigma = .2;
yMax = 5;
Kappa = .1;
A0 = sqrt(3)*3/2;

vMsim = myoPulse.vertexModel( rV, T, m, bndryVs, bndryCs, d0, c2v, Pb, tauL, alphaBar, A, Tpulse, sigma, Kappa, A0, yMax, Fbc, Lc );

Tpull = 10;
vMsim = vMsim.evolve(Tpull);

