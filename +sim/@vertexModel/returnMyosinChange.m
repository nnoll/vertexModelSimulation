function [ deltaM ] = returnMyosinChange( this )
%RETURNMYOSINCHANGE Summary of this function goes here
%   Detailed explanation goes here

    T = length(this.m);
    deltaM = this.m{T} - this.m{1};
    
end

