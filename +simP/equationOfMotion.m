function [ Xdot ] = equationOfMotion( t, X, d0, d0T, d1, l0, kappaL, K, cTopo, A0, Pb, bndryVs, bndryCs, Fbc )
%EQUATIONOFMOTION 

    nVerts = size(d0,2);    
    rv = reshape(X,nVerts,2); %Vertex positions sit on bottom.
    rb = mtimesx(d0,rv); 
    
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    T = kappaL.*(D - l0);
    
    % Calculate pressure
    p = zeros(size(d1,2),1);
    
    if (K > 0)
        r1 = squeeze(rv(:,1));
        r2 = squeeze(rv(:,2));
        for ii = 1:length(cTopo)
            z = size(cTopo(ii).order,2);
            G = diag(ones(z-1,1),1) - diag(ones(z-1,1),-1);
            G(1,z) = -1; G(z,1) = 1;
            
            rcv1 = r1(cTopo(ii).order')';
            rcv2 = r2(cTopo(ii).order')';
            cArea = sum((rcv1*G).*rcv2,2);
            if (length(A0) == 1)
                p(cTopo(ii).clist) = -K*(.5*abs(cArea) - A0);
            else
                p(cTopo(ii).clist) = -K*(.5*abs(cArea) - A0(cTopo(ii).clist));
            end
        end
        p(bndryCs) = Pb; %Put in boundary conditions.
    end

    deltaP = -d1*p;
    P = rb * [0,-1;1,0];
    P = bsxfun(@times,deltaP,P);

    rb = bsxfun(@rdivide,rb,D); %Norm lengths
    
    %Compute instantaneous tension vectors
    Tvec = bsxfun(@times,T,rb);
    
    %Get vertex movement here.
    rv_dot = -mtimesx(d0T,Tvec) + abs(d0T)*P;

    %Apply boundary conditions.
    rv_dot(bndryVs,:,:) = bsxfun(@plus,rv_dot(bndryVs,:,:),Fbc);
    
    %Flatten vertex movement.
    rv_dot = reshape(rv_dot,2*size(rv_dot,1),size(rv_dot,3));
        
    %Vertically concatenant all time derivatives.
    Xdot = rv_dot(:);
    
end

