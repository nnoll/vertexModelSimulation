function [ Fext ] = returnSinusoidExt( t, F, omega, nB, vp, vm, dir )
    % RETURN PULL EXT 

    Fext = zeros(nB,2);
  
    Fext(vm,dir) = -F*cos(omega*t);
    Fext(vp,dir) = F*cos(omega*t);

end

