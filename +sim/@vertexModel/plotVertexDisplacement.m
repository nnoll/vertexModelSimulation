function [ ] = plotVertexDisplacement( this, scale )
    % PLOT 
    
    if (nargin == 1)
        scale = 1;
    end
    
    deltaRF = this.rV{length(this.rV)} - this.rV{1};
    deltaRF = sqrt(sum(deltaRF.^2,2));
    
    rDiff = zeros(size(deltaRF,1),length(this.rV));
    timePts = scale * (1:size(rDiff,2));

    for t = 1:length(this.rV)
        deltaR = this.rV{t} - this.rV{1};
        deltaR = sqrt(sum(deltaR.^2,2));
        rDiff(:,t) = deltaR;% ./ deltaRF;
    end
    
    clf
    cmap = jet(size(rDiff,1));
    for ii = 1:size(rDiff,1)
        plot(timePts,rDiff(ii,:),'Color',cmap(ii,:))
        hold on
    end
    
%     hold all
%     plot(timePts,mean(rDiff,1),'Color','k','LineWidth',5)
%     rDiffTraj = mean(rDiff,1);

end

