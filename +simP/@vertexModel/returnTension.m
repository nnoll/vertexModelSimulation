function [ T ] = returnTension( this )
%PLOTEDGELENGTHS 

    T = zeros(size(this.d0{1},1),length(this.rV));
    for t = 1:length(this.rV)
       T(:,t) = sqrt(sum((this.d0{1}*this.rV{t}).^2,2));
    end
    
    T = T - 1;

end

