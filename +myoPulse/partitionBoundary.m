function [ l, r, u, d ] = partitionBoundary( rV, bV, L )
    %PARTITION BOUNDARY 

    rBV = rV(bV,:);
    l = find(rBV(:,1) - min(rBV(:,1)) < L);
    [~,ind] = sort(rBV(l,2));
    l = l(ind);
    r = find(max(rBV(:,1))-rBV(:,1) < L);
    [~,ind] = sort(rBV(r,2));
    r = r(ind);
    u = find(max(rBV(:,2))-rBV(:,2) < L);
    [~,ind] = sort(rBV(u,1));
    u = u(ind);
    d = find(rBV(:,2) - min(rBV(:,2)) < L);
    [~,ind] = sort(rBV(d,1));
    d = d(ind);
    
end

