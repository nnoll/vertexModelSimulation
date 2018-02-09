function [ Fext ] = returnEndingExt( t, ts, F, delta, nB, vp, vm, dir )
    % RETURN PULL EXT 

    Fext = zeros(nB,2);
    
    if ((t-ts) <= delta && delta > 0)
        Fmag = F - (F*(t-ts))/delta;
    else
        Fmag = 0; 
    end
    
    Fext(vm,dir) = -Fmag/length(vm);
    Fext(vp,dir) = Fmag/length(vp);

end

