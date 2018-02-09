function [ bulk1 ] = returnIntBonds( d0 )
%RETURNINTBONDS 

    Z = sum(abs(d0),1);
    bndry0 = Z < 3;
    bndry1 = (sum(abs(d0(:,bndry0)),2) >= 1);
    badVerts = (sum(abs(d0(bndry1,:)),1) >= 1);
    bulk1 = find(sum(abs(d0(:,badVerts)),2) < 1);
    
end

