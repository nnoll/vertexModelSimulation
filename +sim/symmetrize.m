function [ Z ] = symmetrize( Z )
%SYMMETRIZE Summary of this function goes here
%   Detailed explanation goes here

    Z = [Z,Z(:,1)];

end

