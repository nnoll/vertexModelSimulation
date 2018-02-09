function [ ] = plotTension( this, delta )
%PLOTTENSION Summary of this function goes here
%   Detailed explanation goes here

    if (nargin == 1)
        delta = 1;
    end
    
    T = [this.T{:}];
    plot(T(:,1:delta:size(T,2))')

end

