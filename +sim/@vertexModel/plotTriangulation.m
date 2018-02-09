function [ ] = plotTriangulation( this, color )
    % PLOT TRIANGULATION 
    
    if (nargin == 2)
        patch('Faces',this.Tri{1},'Vertices',this.q0,'FaceVertexCData',color);
        shading faceted
    else
        patch('Faces',this.Tri{1},'Vertices',this.q0,'FaceColor','w','LineWidth',3,'EdgeColor','b');
        hold on
        scatter(this.q0(:,1),this.q0(:,2),75,'Fill','MarkerFaceColor','b');
    end
    axis tight
    
end

