function [ phase, mag, yUnique ] = produceAveragedFcn( this, phase, mag )
%PRODUCEAVERAGEDFCN 

%     phaseV = (abs(this.d0{1}')*phase)./sum(abs(this.d0{1}),1)';
%     magV = (abs(this.d0{1}')*mag)./sum(abs(this.d0{1}),1)';

    % Only store vertical edges.
    rB = abs(this.d0{1}*this.rV{1});
    rB = bsxfun(@rdivide,rB,sqrt(sum(rB.^2,2)));
    
    ind = abs(1-rB(:,2)) <= 1e-3;
    phaseV = phase(ind);
    magV = mag(ind);
  
    r = .5*abs(this.d0{1})*this.rV{1};
%     r = this.rV{1};
    rRound = round(r*10000)/10000;
    rRound = rRound(ind,:);
    yUnique = unique(rRound(:,2));
    
    phase = zeros(size(yUnique));
    mag = zeros(size(yUnique));
    n = 1;
    for y = yUnique'
        ind = (rRound(:,2) == y) & (abs(rRound(:,1)) < 3);
        phase(n) = mean(phaseV(ind));
        mag(n) = mean(magV(ind));
        n = n + 1;
    end

    [yUnique,ind] = sort(yUnique);
    mag = mag(ind);
    phase = phase(ind);

    mag = mag(2:end-1);
    phase = phase(2:end-1);
    yUnique = linspace(-.5,.5,length(yUnique));
    yUnique = yUnique(2:end-1);

end

