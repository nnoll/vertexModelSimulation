function [ J ] = jacobian( t, X, d0, d0T )
    % JACOBIAN 
    
    nBonds = size(d0,1);
    nVerts = size(d0,2);
    T = X(1:nBonds,:);
%     m = X(nBonds+1:2*nBonds); %Myosins sandwiched in the middle.
    rv = reshape(X(nBonds+1:end,:),nVerts,2,size(X,2));
    rb = mtimesx(d0,rv); 
%     Tvec = permute(bsxfun(@times,T,permute(rb,[1,3,2])),[1,3,2]);
%     rv_dot = -mtimesx(d0T,Tvec);
    
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    rb = bsxfun(@rdivide,rb,D);
    
    FxT = -bsxfun(@times,d0T,rb(:,1)');
    FyT = -bsxfun(@times,d0T,rb(:,2)');
    
    Fxx = -d0T*bsxfun(@times,d0,T./D - rb(:,1).*rb(:,1));
    Fxy = -d0T*bsxfun(@times,d0, -rb(:,1).*rb(:,2));
    Fyx = Fxy;
    Fyy = -d0T*bsxfun(@times,d0,T./D - rb(:,2).*rb(:,2));
    
    FTT = bsxfun(@times,d0*FxT,rb(:,1)) + bsxfun(@times,d0*FyT,rb(:,2));
    FTx = bsxfun(@times,d0*Fxx,rb(:,1)) + bsxfun(@times,d0*Fyx,rb(:,2));
    FTy = bsxfun(@times,d0*Fxy,rb(:,1)) + bsxfun(@times,d0*Fyy,rb(:,2));
    
    J = full([FTT,FTx,FTy;FxT,Fxx,Fxy;FyT,Fyx,Fyy]);
    
end

