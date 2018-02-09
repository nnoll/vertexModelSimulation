function [ rV, d0, c2v, bndryVs, bndryEs, bndryCs, T, T_ext, rB_ext, Tri, q ] = lattice( N, L, alpha )
    % LATTICE 
    
    if (nargin == 2)
        alpha = 0;
    end
    
    dim = N*1*[1,1] + .01;
    [ q0 ] = generate.triLattice( 1, dim );
    q0(abs(q0(:,1)-min(q0(:,1))) < 1e-3,:) = [];
    q0(abs(max(q0(:,1))-q0(:,1)) < 1e-3,:) = [];
    
    q0 = q0 + alpha*randn(size(q0));
        
    [rV,C] = voronoin(q0);
    q = L*q0;
    rV = L*rV;
        
    c2v = zeros(length(C),size(rV,1));
    for c = 1:length(C)
       c2v(c,C{c}) = 1; 
    end
    
    badCells = (c2v(:,1) == 1);
    c2v(badCells,:) = [];
    c2v(:,1) = [];
    rV(1,:) = [];
    
    rC = bsxfun(@rdivide,c2v*rV,sum(c2v,2));
    D = pdist2(rC,q);
    [ind,~] = generate.munkres(D);
    q = q(ind,:);

    bulkVerts = find(sum(c2v,1)==3);
    extVerts = find(sum(c2v,1)==2);
    
    rV_ext = rV(extVerts,:);
    
    Tri = zeros(length(bulkVerts),3);
    for v = 1:length(bulkVerts)
       nCells = find(sum(c2v(:,bulkVerts(v)),2) >= 1);
       deltaR = bsxfun(@minus,q(nCells,:),rV(bulkVerts(v),:));
       theta = atan2(deltaR(:,2),deltaR(:,1));
       [~,ind] = sort(theta);
       Tri(v,:) = nCells(ind);
    end
    
    % Build adjacency matrix for verts in cell array.
    adj = zeros(length(bulkVerts));
    for v = 1:length(bulkVerts)
        for vv = v+1:length(bulkVerts)
            if (sum(ismember(Tri(v,:),Tri(vv,:))) == 2)
                adj(v,vv) = 1;
                adj(vv,v) = 1;
            end
        end
    end
    
    % Find bndry verts.
    numZ = sum(adj,2);
    b0 = find(numZ == 2);
    Db = pdist2(rV(bulkVerts(b0),:),rV_ext);
    [bC,~] = generate.munkres(Db);
    rB_ext = rV_ext(bC,:) - rV(bulkVerts(b0),:);
    
    T_ext = zeros(length(b0),1);
    for b = 1:length(b0)
        bCells = find(prod(c2v(:,[bulkVerts(b0(b)),extVerts(bC(b))]),2)==1);
        T_ext(b) = norm(q(bCells(1),:)-q(bCells(2),:));
    end
    
    rB_ext = bsxfun(@rdivide,rB_ext,sqrt(sum(rB_ext.^2,2)));
    bndryVs = b0;
    
    % Store the tensions as given by edge lengths in the triangulation.
    numBulkEdges = sum(adj(:))/2;
    edgeVerts = zeros(numBulkEdges,2);
    T = zeros(numBulkEdges,1);
    
    n = 1;
    for v = 1:size(adj,1)
        nVerts = find(adj(v,:));
        nVerts = nVerts(nVerts>v);
        for nv = nVerts
           edgeVerts(n,:) = [v,nv]; 
           v1 = bulkVerts(v);
           v2 = bulkVerts(nv);
           bCells = find(prod(c2v(:,[v1,v2]),2) == 1);
           T(n) = norm(q(bCells(1),:) - q(bCells(2),:));
           n = n + 1;
        end
    end

    d0 = zeros(numBulkEdges,length(bulkVerts));
    d0((1:numBulkEdges)' + numBulkEdges*(edgeVerts(:,1)-1)) = 1;
    d0((1:numBulkEdges)' + numBulkEdges*(edgeVerts(:,2)-1)) = -1;
    
    c2v = c2v(:,bulkVerts);
    numV = sum(c2v,2);
    c2v(numV==0,:) = [];
    q(numV==0,:) = [];
    rV = rV(bulkVerts,:);
    bndryEs = (sum(abs(d0(:,bndryVs)),2) >= 1);
    extVs = (sum(abs(d0(bndryEs,:)),1) >= 1);
    bndryEs = find(sum(abs(d0(:,extVs)),2) == 2);
    bndryCs = find(sum(c2v(:,extVs),2) >= 1)';
    
    Tri = zeros(size(c2v,2),3);
    for v = 1:size(c2v,2)
       nCells = find(c2v(:,v));
       deltaR = bsxfun(@minus,q(nCells,:),rV(v,:));
       theta = atan2(deltaR(:,2),deltaR(:,1));
       [~,ind] = sort(theta);
       Tri(v,:) = nCells(ind);
    end
    
    Tmed = median(T);
    T = T/Tmed;
    T_ext = T_ext/Tmed;
    

end

