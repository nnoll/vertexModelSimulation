function [ r ] = triLattice( a, dim )
    % GENERATE TRIANGULAR LATTICE 

    b1 = a*[1;0];
    b2 = a*[1/2;sqrt(3)/2];

    N = round(dim(2)/a);
    r = zeros(2,(2*N+1)^2);
    n = 1;
    for ii = -N:N
        for jj = -N:N
            r(:,n) = ii*b1 + jj*b2; 
            n = n + 1;
        end
    end

    r = r + repmat((dim+1)',[1,size(r,2)])/2;
    r = bsxfun(@minus,r,mean(r,2));
    
    indX = ( (r(1,:) >= -dim(1)/2) .* (r(1,:) <= dim(1)/2) .* (r(2,:) >= -dim(2)/2) .* (r(2,:) <= dim(2)/2) );
    r = r(:,indX==1);
    r = r';
    
end

