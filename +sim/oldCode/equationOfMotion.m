function [ Xdot ] = equationOfMotion(t, X, nu, omega, d0, bndryVs, l, r, d, u, Lx, Ly, Fbc, Fext, extVs, P)
%EQUATIONOFMOTION 

    nBonds = size(d0,1);
    nVerts = size(d0,2);

    T = X(1:nBonds);
    m = X(nBonds+1:2*nBonds); %Myosins sandwiched in the middle.
    rv = reshape(X(2*nBonds+1:end),nVerts,2); %Vertex positions sit on bottom.
    
    rb = d0*rv; 
    
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    rb = bsxfun(@rdivide,rb,D); %Norm lengths
    
    %Compute instantaneous tension vectors
    Tvec = bsxfun(@times,T,rb);
    
    %Get vertex movement here.
    rv_dot = -d0'*(Tvec);
    
    %Apply boundary conditions.
    rv_dot(bndryVs,:) = rv_dot(bndryVs,:) + Fbc;
%     rv_dot(bndryVs,:) = rv_dot(bndryVs,:) + Fext(t);
%     rv_dot(extVs,:) = rv_dot(extVs,:) - .5 * P * (rv(circshift(extVs,[0,-1]),:)-rv(circshift(extVs,[0,1]),:)) * [0,-1;1,0];
    
%     kappa = 0;
%     rv_dot(bndryVs(l),1) = rv_dot(bndryVs(l),1) - kappa*(rv(bndryVs(l),1)+Lx).*heaviside(-(rv(bndryVs(l),1)+Lx));
%     rv_dot(bndryVs(r),1) = rv_dot(bndryVs(r),1) - kappa*(rv(bndryVs(r),1)-Lx).*heaviside(rv(bndryVs(r),1)-Lx);
%     rv_dot(bndryVs(d),2) = rv_dot(bndryVs(d),2) - kappa*(rv(bndryVs(d),2)+Ly).*heaviside(-(rv(bndryVs(d),2)+Ly));
%     rv_dot(bndryVs(u),2) = rv_dot(bndryVs(u),2) - kappa*(rv(bndryVs(u),2)-Ly).*heaviside(rv(bndryVs(u),2)-Ly);
    
    %Get bond length change here.
    deltaR_dot = d0*rv_dot;
    rdot = dot(deltaR_dot,rb,2);
    
    %Flatten vertex movement.
    rv_dot = rv_dot(:);
    
    w = (nu*(D-T))./(m); %Walking component.  
    T_dot = rdot + nu*(D-T) - w.*T;
    
    %Get myosin dot here.
    m_dot = omega*(T - m);% omega*(T - m).*(T./m - 1);
    
    %Vertically concatenant all time derivatives.
    Xdot = vertcat(T_dot, m_dot, rv_dot);

end

