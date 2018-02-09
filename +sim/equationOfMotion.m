function [ Xdot ] = equationOfMotion(t, X, nu, omega, d0, d0T, bndryVs, Fbc, Fext)
%EQUATIONOFMOTION 

    nBonds = size(d0,1);
    nVerts = size(d0,2);
    
    T = X(1:nBonds,:);
    m = X(nBonds+1:2*nBonds,:); %Myosins sandwiched in the middle.
    rv = reshape(X(2*nBonds+1:end,:),nVerts,2,size(X,2)); %Vertex positions sit on bottom.
    rb = mtimesx(d0,rv); 
    
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    rb = bsxfun(@rdivide,rb,D); %Norm lengths
    
    %Compute instantaneous tension vectors
    Tvec = permute(bsxfun(@times,T,permute(rb,[1,3,2])),[1,3,2]);
    
    %Get vertex movement here.
    rv_dot = -mtimesx(d0T,Tvec);

    %Apply boundary conditions.
    rv_dot(bndryVs,:,:) = bsxfun(@plus,rv_dot(bndryVs,:,:),Fext(t));
    rv_dot(bndryVs,:,:) = bsxfun(@plus,rv_dot(bndryVs,:,:),Fbc);
    
    %Get bond length change here.
    deltaR_dot = mtimesx(d0,rv_dot);
    rdot = squeeze(dot(deltaR_dot,rb,2));
    
    %Flatten vertex movement.
    rv_dot = reshape(rv_dot,2*size(rv_dot,1),size(rv_dot,3));
    
    l = squeeze(D) - T;
    T_dot = rdot - nu*l.*(T./m - 1); 
    
    %Get myosin dot here.
    m_dot = omega*(T - m);
    
    %Vertically concatenant all time derivatives.
    Xdot = vertcat(T_dot, m_dot, rv_dot);
    
end

