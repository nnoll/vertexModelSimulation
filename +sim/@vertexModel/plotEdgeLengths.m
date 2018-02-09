function [ ] = plotEdgeLengths( this, delta )

    if (nargin == 1)
        delta = 1;
    end
    
    timePts = 1:delta:length(this.rV);
    RB = zeros(size(this.d0{1},1),length(this.rV));
      
    n = 1;
    for t = timePts
        
        if (n <= length(this.c2v) && this.TE{n} <= t)
            D = this.d0{n};
            n = n + 1;
        end
          
        r = D * this.rV{t};
        RB(:,t) = sqrt(sum(r.^2,2));
       
    end
      
    plot(timePts,RB)
       

end

