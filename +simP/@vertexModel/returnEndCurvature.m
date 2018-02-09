function [ kappa ] = returnEndCurvature( this )
%RETURNENDCURVATURE 

    T = this.returnTension();
    P = this.returnPressure();
    
    R = T + 1;
    dP = this.d1{1}*P;
    dP = abs(dP);
    
    kappa = dP(:,size(dP,2)).*R(:,size(R,2)) ./ T(:,size(T,2));
    kappa = kappa(this.bulkE{1});
    
end

