function [ Xdot ] = equationOfMotion( t, X, d0, d0T, d1, K, cTopo, A0, tauL, alpha, Pb, bndryVs, bndryCs, omegaPulse, A, phiPulse,  Fbc )
%EQUATIONOFMOTION 

    nBonds = size(d0,1);
    nVerts = size(d0,2);
    
    T = X(1:nBonds,:);
    m = X(nBonds+1:2*nBonds,:); %Myosins sandwiched in the middle.
    
    % Add in `perturbed myosin' due to medial pulse
    mCell = A*sin((omegaPulse*t)+phiPulse);
    mCell(mCell<(A/1.05)) = 0;
    mCell = abs(d1)*mCell;
    m = m + mCell;
    
    rv = reshape(X(2*nBonds+1:end,:),nVerts,2,size(X,2)); %Vertex positions sit on bottom.
    rb = mtimesx(d0,rv); 
    
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    
    % Calculate pressure
    p = zeros(size(d1,2),1);
    p(bndryCs) = Pb; %Put in boundary conditions.
    
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
    
    %Get bond length change here.
    deltaR_dot = mtimesx(d0,rv_dot);
    rdot = squeeze(dot(deltaR_dot,rb,2));

    l = squeeze(D) - T;
    T_dot = rdot - l.*(T./m - 1)/tauL; 
    
    %Get myosin dot here.
    m_dot = alpha*(T - m)/tauL;
    
    %Flatten vertex movement.
    rv_dot = reshape(rv_dot,2*size(rv_dot,1),size(rv_dot,3));
        
    %Vertically concatenant all time derivatives.
    Xdot = vertcat(T_dot, m_dot, rv_dot);

    
end

