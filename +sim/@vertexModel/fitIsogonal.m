function [ theta ] = fitIsogonal( this )
    % FIT ISOGONAL 
    
    deltaR = this.rV{length(this.rV)} - this.rV{1};
    deltaR = bsxfun(@minus,deltaR,mean(deltaR,1));
    
    Tri = this.Tri{1};
    
    % Construct inversion matrix.
    F = zeros(2*size(deltaR,1),max(Tri(:)));
    R = [0,-1;1,0];
    
    q21 = this.q0(Tri(:,2),:) - this.q0(Tri(:,1),:);
    q32 = this.q0(Tri(:,3),:) - this.q0(Tri(:,2),:);
    q13 = this.q0(Tri(:,1),:) - this.q0(Tri(:,3),:);
        
    triArea = .5*abs(q21(:,1).*q32(:,2) - q21(:,2).*q32(:,1));
    
    q21 = q21*R;
    q32 = q32*R;
    q13 = q13*R;
    
    q21 = bsxfun(@rdivide,q21,triArea);
    q32 = bsxfun(@rdivide,q32,triArea);
    q13 = bsxfun(@rdivide,q13,triArea);
    
    V = size(Tri,1);

    F((1:V)' + 2*V*(Tri(:,3)-1)) = q21(:,1); 
    F((1:V)' + 2*V*(Tri(:,1)-1)) = q32(:,1); 
    F((1:V)' + 2*V*(Tri(:,2)-1)) = q13(:,1); 
    
    F(((V+1):2*V)' + 2*V*(Tri(:,3)-1)) = q21(:,2); 
    F(((V+1):2*V)' + 2*V*(Tri(:,1)-1)) = q32(:,2); 
    F(((V+1):2*V)' + 2*V*(Tri(:,2)-1)) = q13(:,2);
    
%     F = [F;ones(1,size(F,2))];
    theta = pinv(F) * deltaR(:);

end

