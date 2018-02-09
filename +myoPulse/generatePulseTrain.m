function [ Tpulse ] = generatePulseTrain( T, sigma, ncells, Tfinal )
    % GENERATE PULSE TRAIN 

    niter = round(3 * Tfinal/T);
    times = T + sigma*randn(ncells,niter);
    
    Tpulse = cumsum(times,2);
    
end

