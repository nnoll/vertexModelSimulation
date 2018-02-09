function [ r, deltaR ] = plotDisplacementVF( this, t )
%PLOTDISPLACEMENTVECTOR 

    if (nargin == 1)
        t = length(this.rV);
    end
    
    deltaR = this.rV{t} - this.rV{1};
    this.plotSnapshot(1)
    hold all
    quiver(this.rV{1}(:,1),this.rV{1}(:,2),deltaR(:,1),deltaR(:,2),2,'r','LineWidth',2)
    axis tight
    
    r = this.rV{1};
%     rB = .5*abs(this.d0{1})*r;
    [~,ind] = min(sum(r.^2,2));
    r = bsxfun(@minus,r,r(ind,:));
    
end

