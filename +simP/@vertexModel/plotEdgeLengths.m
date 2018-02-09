function [ L ] = plotEdgeLengths( this )
%PLOTEDGELENGTHS 

    L = zeros(size(this.d0{1},1),length(this.rV));
    for t = 1:length(this.rV)
       L(:,t) = sqrt(sum((this.d0{1}*this.rV{t}).^2,2));
    end
    
    plot(L')


end

