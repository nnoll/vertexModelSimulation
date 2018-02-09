function [ ] = plotSnapshot( this, t, color )
%PLOTSNAPSHOT 
    
    flip = 0;
    for n = 1:length(this.c2v)
        if (t >= this.TE{n} && flip == 0)
            [ face ] = sim.returnFace( this.c2v{n}, this.rV{t} );
            intVerts = 1:size(this.c2v{n},2);
            intVerts(this.bV{1}) = [];
            bulkCells = sum(this.c2v{n}(:,intVerts),2) > 2;
        end
    end
    
    if (nargin == 2)
        patch('Faces',face(bulkCells,:),'Vertices',this.rV{t},'FaceColor','w','EdgeColor','k','LineWidth',2);
    else
        patch('Faces',face(bulkCells,:),'Vertices',this.rV{t},'FaceVertexCData',color,'EdgeColor','k','LineWidth',2);
        shading faceted
    end
    hold on
    scatter(this.rV{t}(:,1),this.rV{t}(:,2),'Fill','MarkerFaceColor','b');

end

