function [ ratioTraj ] = plotTensionMyoRatio( this, scale )
    % PLOT 
    
    
    if (nargin == 1)
        scale = 1;
    end
    
    for t = 1:length(this.rV)
        ratio(:,t) = (1 - this.T{t} ./ this.m{t}).^2;
    end
    
    clf
    timePts = scale * (1:size(ratio,2));
    cmap = jet(size(ratio,1));
    for ii = 1:size(ratio,1)
        plot(this.omega*timePts,ratio(ii,:),'Color',cmap(ii,:))
        hold on
    end
    
    hold all
    plot((this.omega*timePts),mean(ratio,1),'Color','k','LineWidth',5)
    set(gca,'yscale','log')
    ratioTraj = mean(ratio,1);

end

