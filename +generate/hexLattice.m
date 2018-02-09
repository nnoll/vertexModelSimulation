function [ rV, d0, c2v, bndryVs, bndryCs, T_bc, rB_bc ] = hexLattice( nLayers, L, alpha )
%HEXLATTICE 

    N = 1 + 3*nLayers*(nLayers+1);
    [ rv, adj ] = generate.seedLattice( N, L );

    [ ~, rv2, ~ ] = generate.find_boundary_edges( adj, rv );
    [ ~, c2v ] = generate.get_cells( rv2 );
    [ ~, rV, adj ] = generate.find_boundary_edges( adj, rv );
    [ rV, adj, c2v ] = generate.commit_T1s( rV, adj, c2v, 0 );
    
    rV = rV + alpha*randn(size(rV));
    
    Z = sum(adj,1);
    extVs = find(Z==1);
    
    [ adj ] = generate.clean_adj( adj );
    nbonds = sum(adj(:))/2;
    d0 = zeros(nbonds,size(adj,1));
    nb = 1;
    for v = 1:(size(adj,1)-1)
        for vv = v+1:size(adj,1)
            if (adj(v,vv) == 1)
                d0(nb,v) = 1;
                d0(nb,vv) = -1;
                nb = nb + 1;
            end
        end
    end
    
    rB_bc = zeros(size(extVs,1),2);
    bndryVs = zeros(size(extVs));
    n = 1;
    rV = rV';
    for v = extVs
       bndryVs(n) = find((adj(v,:)==1));
       delta = rV(v,:) - rV(bndryVs(n),:);
       delta = delta/sqrt(sum(delta.^2,2));
       rB_bc(n,:) = delta;
       n = n + 1;
    end
    
    T_bc = ones(size(rB_bc,1),1);
    
    rV(extVs,:) = [];
    d0(:,extVs) = [];
    
    badBonds = sum(abs(d0),2) < 2;
    d0(badBonds,:) = [];
    
    bndryCs = find(sum(c2v(:,extVs),2) >= 1)';
%     c2v(badCells,:) = [];
    c2v(:,extVs) = [];
    
    shiftBndryVs = zeros(size(bndryVs));
    for ii = 1:length(extVs)
       shiftBndryVs(bndryVs > extVs(ii)) = shiftBndryVs(bndryVs > extVs(ii))-1;
    end
    
    bndryVs = bndryVs + shiftBndryVs;
    
end

