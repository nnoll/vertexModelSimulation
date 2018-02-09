function [ rV, d0, c2v, bndryVs, bndryEs, bndryCs, T, T_ext, rB_ext, Tri, q ] = rectLattice( N, L, alpha )
    % LATTICE 
    
    if (nargin == 2)
        alpha = 0;
    end
    
    dim = N*1*[1,1] + .01;
    [ q0 ] = generate.triLattice( 1, dim );
    
    % Modify q so that we form a `rectangle' in tension space.
    minX = min(q0(:,1)); maxX = max(q0(:,1));
    minY = min(q0(:,2)); maxY = max(q0(:,2));
    
    indL = (abs(q0(:,1) - minX) <= 1e-3);
    q0(indL,1) = minX + 1/2;
    
    indR = (abs(maxX - q0(:,1)) <= 1e-3);
    q0(indR,1) = maxX - 1/2;

    minX = minX + 1/2;
    maxX = maxX - 1/2;

    delta = 1/4;
    bndry = [minX-delta,minY-delta;minX-delta,maxY+delta;maxX+delta,maxY+delta;maxX+delta,minY-delta];
    
    ind = (abs(max(q0(:,1)) - q0(:,1)) <= 1e-3) | (abs(q0(:,1)-min(q0(:,1))) <= 1e-3) ...
        | (abs(max(q0(:,2)) - q0(:,2)) <= 1e-3) | (abs(q0(:,2)-min(q0(:,2))) <= 1e-3);
       
    qB = q0(ind,:);
    
    q0 = q0 + alpha*randn(size(q0));
%     q0(ind,:) = qB;
    
    [rV,C,q] = generate.voronoiLimit(q0(:,1),q0(:,2),'bs_ext',bndry,'figure','off');
    
    q = L*q;
    rV = L*rV;
    bndry = L*bndry;
        
    c2v = zeros(length(C),size(rV,1));
    for c = 1:length(C)
       c2v(c,C{c}) = 1; 
    end
    
    % Obtain list of boundary verts.
    bndryVs = (abs(rV(:,1) - min(bndry(:,1))) <= 1e-3) | (abs(rV(:,1) - max(bndry(:,1))) <= 1e-3) | ...
              (abs(rV(:,2) - min(bndry(:,2))) <= 1e-3) | (abs(rV(:,2) - max(bndry(:,2))) <= 1e-3);
    rV_ext = rV(bndryVs,:);
    
    bulkVerts = 1:size(rV,1);
    bulkVerts(bndryVs) = [];
    bndryVs = find(bndryVs);
    
    
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
    
    T_ext1 = zeros(length(b0),1);
    for b = 1:length(b0)
        bCells = find(prod(c2v(:,[b0(b),bndryVs(bC(b))]),2)==1);
        T_ext1(b) = norm(q(bCells(1),:)-q(bCells(2),:));
    end
    
    extVert = find(numZ == 1);
    T_ext2 = zeros(length(extVert),1);
    n = 1;
    for v = extVert'
        nVert = find(adj(v,:));
        rB_ext = [rB_ext;rV(v,:) - rV(nVert,:)];
        bCells = find(prod(c2v(:,[nVert,v]),2)==1);
        T_ext2(n) = norm(q(bCells(1),:)-q(bCells(2),:));

        b0 = [b0;nVert];
        n = n + 1;
    end
    
    rB_ext = bsxfun(@rdivide,rB_ext,sqrt(sum(rB_ext.^2,2)));
    adj(extVert,:) = [];
    adj(:,extVert) = [];
    
    Q = [q,zeros(size(q,1),1)];
    
    % Find new indices of indL and indR
    DL = pdist2(L*q0(indL,:),q);
    DR = pdist2(L*q0(indR,:),q);
    indL = generate.munkres(DL);
    indR = generate.munkres(DR);
% 
%     if (alpha < .1)
%         Q(indL,3) = -L/2;
%         Q(indR,3) = -L/2;
%     end
    
    [ rV ] = generate.atnVerts( Tri, Q );
    rV(extVert,:) = [];
    bulkVerts(extVert) = [];
    
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
    
    T_ext = [T_ext1;T_ext2]; 
    bndryVs = b0;
    shiftBndryVs = zeros(size(b0));
    for ii = 1:length(extVert)
        shiftBndryVs(bndryVs >= extVert(ii)) = shiftBndryVs(bndryVs >= extVert(ii)) - 1;
    end
    bndryVs = bndryVs + shiftBndryVs;
    
    c2v = c2v(:,bulkVerts);
    numV = sum(c2v,2);
    c2v(numV==0,:) = [];
    q(numV==0,:) = [];
    
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

