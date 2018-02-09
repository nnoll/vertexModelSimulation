function [value,terminal,direction] = T1event(t, X, Ts, d0, Lc, bulkEdges)
    %T1 EVENT
    
    terminal = 1;
    direction = -1;
    
    %Compute value
    N = size(d0,2);
    rv = reshape(X(2*size(d0,1)+1:end),N,2); %Vertex positions sit on bottom.
    rb = d0*rv; 
    rB = .5*abs(d0)*rv;
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    
    if (sum((D(bulkEdges) < Lc) & (sqrt(sum(rB(bulkEdges,:).^2,2)) <= 6)) < 1 || t <= Ts + .001)
        value = 1;
    else
        value = 0;
    end

end

