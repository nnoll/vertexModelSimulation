function [ nrv, nadj, nc2v ] = commit_T1s( rv, adj, c2v, niter )
%COMMIT_T1S Takes in a generated hexagonal lattice and performs T1s to
%perturb it off regular topology.

nrv = rv;
nadj = zeros(size(adj));
nc2v = c2v;

%% Seperate external verts for now.
nnverts = sum(adj,1);
ext_verts = find(nnverts==1);
bulk_verts = find(nnverts>1);
ext_cverts = zeros(size(ext_verts)); %Vertices we connect bc's to.
ii = 1;
for ev = ext_verts
    ext_cverts(ii) = find(adj(:,ev),1);
    ii = ii + 1;
end
clear ii

%% Clear out ext. verts from our data structures
rv(:,ext_verts) = [];
adj(:,ext_verts) = [];
adj(ext_verts,:) = [];
c2v(:,ext_verts) = [];

shift_cverts = zeros(size(ext_cverts));
for ev=ext_verts
    shift_cverts(ext_cverts > ev) = shift_cverts(ext_cverts > ev) - 1;
end
shext_cverts = ext_cverts + shift_cverts;

%% Create bond LUT and adjacency matrix.
nbonds = sum(adj(:))/2;
N = size(adj,1);
LUT = zeros(nbonds,N);
ii = 1;
for n=1:N
    nv = find(adj(n,:));
    nv = nv(find(nv>n));
    for nn=nv
        LUT(ii,n) = 1; LUT(ii,nn) = -1;
        ii = ii + 1;
    end
end
rb = LUT*rv'; 
D = sqrt(sum(rb.^2,2));
[adj_b] = generate.create_adjB(LUT,N,nbonds);
flipped_bonds = [];

for nn=1:niter
    contin = 1;
    while (contin)
        T1_bond = randi(nbonds,1);
%         T1_bond = 154;
        bverts = find(LUT(T1_bond,:));
        contin = any(ismember(bverts,shext_cverts)) || ismember(T1_bond,flipped_bonds);
    end
%     T1_bond
    flipped_bonds = horzcat(flipped_bonds,T1_bond);
    %% Change topology.
    [ rot_sign, adj, LUT, adj_b ] = generate.change_topology( adj, LUT, T1_bond, bverts, rv, D, adj_b );

    %% Compute new vertex positions after applying a T1
    rc_t1 = .5*abs(LUT(T1_bond,:))*rv'; 
    rc_t1 = rc_t1';
    delta_b = 1.2*[0, -1; 1, 0]*rb(T1_bond,:)'./(ones(2,1)*D(T1_bond)');
    delta_b = delta_b*rot_sign;
    rv_t1 = zeros(2,2);
    rv_t1(1,:) = rc_t1 + .5*ones(size(delta_b,1),1)*LUT(T1_bond,bverts(1)).*delta_b;
    rv_t1(2,:) = rc_t1 + .5*ones(size(delta_b,1),1)*LUT(T1_bond,bverts(2)).*delta_b;

    %% Update graph geometry.
    rv(:,bverts(1)) = rv_t1(1,:); 
    rv(:,bverts(2)) = rv_t1(2,:);
    
    %% Update cell data structure.
    [ c2v ] = generate.update_c2v( c2v, adj, bverts );
end

%% Update output variables
for ii=1:size(rv,2)
    nrv(:,bulk_verts(ii)) = rv(:,ii);
end

for ii=1:size(adj,1)
    for jj=1:size(adj,2)
        nadj(bulk_verts(ii),bulk_verts(jj)) = adj(ii,jj);        
    end
end

for ii=1:length(ext_verts)
    nadj(ext_verts(ii),ext_cverts(ii)) = 1;
    nadj(ext_cverts(ii),ext_verts(ii)) = 1;
end

for ii=1:size(c2v,2)
    nc2v(:,bulk_verts(ii)) = c2v(:,ii);
end

end

