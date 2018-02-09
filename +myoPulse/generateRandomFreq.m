function [ omega, phi ] = generateRandomFreq( T, sigma, ncells )
%GENERATERANDOMFREQ 

    omega = (2*pi/T) * ones(ncells,1);
    omega = omega + .5*pi/T*randn(ncells,1);
    
    phi = pi * sigma * randn(ncells,1);

end

