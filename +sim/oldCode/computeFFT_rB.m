function [ phase, mag ] = computeFFT_rB( this, f0 )
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
    
%     Ts = round(1/f0);
%     RB = RB(:,Ts+1:size(RB,2));

    L = size(RB,2);
    n = 2^nextpow2(L);
    RB = RB - 1;

    Ttwidle = fft(RB,n,2);
    Ttwidle = (Ttwidle/n);
    Ttwidle = Ttwidle(:,1:n/2+1);
    Ttwidle(:,2:end-1) = 2*Ttwidle(:,2:end-1);

    Fs = 1;
    f = Fs*(0:(n/2))/n;
%     plot(f,abs(Ttwidle))
%     pause
    [~,ind] = max(mean(abs(Ttwidle),1));
    Ttwidle = Ttwidle(:,ind);
    phase = angle(Ttwidle);
    mag = abs(Ttwidle);

end

