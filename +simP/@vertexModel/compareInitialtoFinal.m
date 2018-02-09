function [ ] = compareInitialtoFinal( this )
%COMPAREINITIALTOFINAL Summary of this function goes here
%   Detailed explanation goes here

    [ face0 ] = sim.returnFace( this.c2v{1}, this.rV{1} );
    [ faceF ] = sim.returnFace( this.c2v{end}, this.rV{end} );

    intVerts = 1:size(this.c2v{1},2);
    intVerts(this.bV{1}) = [];
    bulkCells0 = (sum(this.c2v{1}(:,intVerts),2) > 2);
    intVerts = 1:size(this.c2v{end},2);
    intVerts(this.bV{1}) = [];
    bulkCellsF = (sum(this.c2v{end}(:,intVerts),2) > 2);

    patch('Faces',face0(bulkCells0,:),'Vertices',this.rV{1},'FaceColor','none','EdgeColor','b','LineWidth',2);
    hold all
    patch('Faces',faceF(bulkCellsF,:),'Vertices',this.rV{length(this.rV)},'FaceColor','none','EdgeColor','r','LineWidth',2);

end

