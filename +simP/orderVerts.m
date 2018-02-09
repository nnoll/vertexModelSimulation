function [ order ] = orderVerts( r )
    % ORDER VERTS 
    
    R = mean(r,1);
    r = bsxfun(@minus,r,R);
    
    theta = atan2(r(:,2),r(:,1));
    theta(theta<0) = theta(theta<0) + 2*pi;
    [~,order] = sort(theta);

end

