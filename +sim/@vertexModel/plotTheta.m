function [ theta ] = plotTheta( this )
    %PLOT THETA

    [ theta ] = this.fitIsogonal();
    [ face ] = sim.returnFace( this.c2v{1}, this.rV{1} );
    patch('Faces',face,'Vertices',this.rV{1},'FaceVertexCData',theta);
    shading faceted

end

