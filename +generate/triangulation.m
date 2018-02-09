function [ Tri ] = triangulation( rV, c2v, bulkVerts )
    % TRIANGULATION 

    q = zeros(size(c2v,1),2);
    for c = 1:size(c2v,1)
       q(c,:) = mean(rV(c2v(c,:)==1,:),1); 
    end
    
    Tri = zeros(length(bulkVerts),3);
    for v = 1:length(bulkVerts)
       nCells = find(sum(c2v(:,bulkVerts(v)),2) >= 1);
       deltaR = bsxfun(@minus,q(nCells,:),rV(bulkVerts(v),:));
       theta = atan2(deltaR(:,2),deltaR(:,1));
       [~,ind] = sort(theta);
       Tri(v,:) = nCells(ind);
    end


end

