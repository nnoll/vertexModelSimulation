function [ phase, mag ] = computeFFT_rB( this, f0, deltaT )
%TAKEFFT_LENGTH Summary of this function goes here
%   Detailed explanation goes here

    
    timePts = 1:length(this.rV);
    RB = zeros(size(this.d0{1},1),length(this.rV));
      
    n = 1;
    for t = timePts
        
        if (n <= length(this.c2v) && this.TE{n} <= t)
            D = sparse(this.d0{n});
            n = n + 1;
        end
          
        r = D * this.rV{t};
        RB(:,t) = sqrt(sum(r.^2,2));
       
    end
    
    omega = 2*pi*f0;
    RB = bsxfun(@minus,RB,mean(RB,2));
    mag = max(abs(RB),[],2);
    
    phase = zeros(size(mag));
    timePts = deltaT*(timePts - 1);
    for ii = 1:length(mag)
        [fMax,xMax] = findpeaks(RB(ii,:)/mag(ii));
        [fMin,xMin] = findpeaks(-RB(ii,:)/mag(ii));
        mag(ii) = mag(ii)*median([fMax,fMin]);
        
        pks = sort(timePts([xMax,xMin]));
        z = pi*(0:(length(pks)-1));
        phase(ii) = mod(median(z-omega*pks)+pi,2*pi)-pi;
    end

end

