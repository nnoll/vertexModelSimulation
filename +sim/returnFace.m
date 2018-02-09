function [ face ] = returnFace( c2v, rV )
    % RETURN FACE 
    
    z = sum(c2v,2);
    zMax = max(z);
    
    face = nan*ones(size(c2v,1),zMax);
    for c = 1:size(c2v,1)
        cVerts = find(c2v(c,:));
        Rc = mean(rV(cVerts,:));
        deltaR = bsxfun(@minus,rV(cVerts,:),Rc);
        theta = atan2(deltaR(:,2),deltaR(:,1));
        theta = mod(theta,2*pi);
        [~,ind] = sort(theta);
        face(c,1:length(cVerts)) = cVerts(ind);
    end


end

