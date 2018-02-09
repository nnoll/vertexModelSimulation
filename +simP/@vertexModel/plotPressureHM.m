function [ ] = plotPressureHM( this )
    % PLOT PRESSURE HM 

    [ p ] = this.returnPressure();
    fig = figure('units','normalized','outerposition',[0 0 1 1]);

    if (nargin == 1)
        delta = 1;
    end
    
    XLim = [100,-100];
    YLim = [100,-100];
    CLim = [inf,-inf];

    n = 1;
    for t = 1:delta:length(this.rV)
        
       if (n <= length(this.c2v) && this.TE{n} <= t)
           [ face ] = simP.returnFace( this.c2v{n}, this.rV{t} );
           intVerts = 1:size(this.c2v{n},2);
           intVerts(this.bV{1}) = [];
           bulkCells = sum(this.c2v{n}(:,intVerts),2) > 2;
           n = n + 1;
       end
       
       clf
       patch('Faces',face(bulkCells,:),'Vertices',this.rV{t},'FaceVertexCData',p(bulkCells,t),'EdgeColor','k','LineWidth',2);
       axis equal
       tmpX = get(gca,'XLim');
       
       if (tmpX(1) < XLim(1))
           XLim(1) = tmpX(1);
       end
       
       if (tmpX(2) > XLim(2))
           XLim(2) = tmpX(2);
       end
       
       tmpY = get(gca,'YLim');
       
       if (tmpY(1) < YLim(1))
           YLim(1) = tmpY(1);
       end
       
       if (tmpY(2) > YLim(2))
           YLim(2) = tmpY(2);
       end
       
       tmpC = get(gca,'CLim');
       
       if (tmpC(1) < CLim(1))
           CLim(1) = tmpC(1);
       end
       
       if (tmpC(2) > CLim(2))
           CLim(2) = tmpC(2);
       end
       
    end
    
    XLim(1) = XLim(1) - 1;
    XLim(2) = XLim(2) + 1;
    YLim(1) = YLim(1) - 2;
    YLim(2) = YLim(2) + 2;
    
    n = 1;
    nn = 1;
    
    for t = 1:delta:length(this.rV)
        
       clf
       if (nn <= length(this.c2v) && (this.TE{nn}+1) <= t)
           [ face ] = sim.returnFace( this.c2v{nn}, this.rV{t} );
           intVerts = 1:size(this.c2v{nn},2);
           intVerts(this.bV{1}) = [];
           bulkCells = sum(this.c2v{nn}(:,intVerts),2) > 2;
           nn = nn + 1;
       end
       
       patch('Faces',face(bulkCells,:),'Vertices',this.rV{t},'FaceVertexCData',p(bulkCells,t),'EdgeColor','k','LineWidth',2);
       set(gca,'XLim',XLim,'YLim',YLim)
       shading faceted
       colorbar()
       hold on
       scatter(this.rV{t}(:,1),this.rV{t}(:,2),'Fill')
       mov(n) = getframe(fig);
       n = n + 1;
       
    end
    
end

