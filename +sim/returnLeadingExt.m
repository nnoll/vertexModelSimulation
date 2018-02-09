function [ Fext ] = returnLeadingExt( t, F, delta, nB, vp, vm, dir )
    % RETURN PULL EXT 

    Fext = zeros(nB,2);
    
    if (t <= delta && delta > 0)
        Fmag = (F*t)/delta;
    else
        Fmag = F; 
    end
    
%     indM = 1:length(vm);
%     indP = 1:length(vp);

    Fext(vm,dir) = -Fmag/length(vm); %*(indM-1)/length(vm); %.*(length(vm)-indM)
    Fext(vp,dir) = Fmag/length(vp); %*(indP-1)/length(vp); %.*(length(vp)-indP)

%     Fext(vm(end-1:end),2) = 0;
%     Fext(vm(1:2),2) = 0;
    
%     Fext(vp(end-1:end),2) = 0;
%     Fext(vp(1:2),2) = 0;

end

