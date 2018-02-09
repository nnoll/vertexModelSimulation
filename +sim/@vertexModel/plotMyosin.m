function [ ] = plotMyosin( this, delta )
%PLOTTENSION Summary of this function goes here
%   Detailed explanation goes here

    if (nargin == 1)
        delta = 1;
    end
    
    M = [this.m{:}];
    plot(M(:,1:delta:size(M,2))')

end

