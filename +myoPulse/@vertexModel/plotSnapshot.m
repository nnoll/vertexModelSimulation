function [ ] = plotSnapshot( this, t, color )
    % PLOTLATTICE Summary of this function goes here
    %   Detailed explanation goes here

   if (nargin < 3)
       color = 'k';
   end
   
   [ face ] = simP.returnFace( this.c2v{1}, this.rV{t} );
   intVerts = 1:size(this.c2v{1},2);
   intVerts(this.bV{1}) = [];
   bulkCells = sum(this.c2v{1}(:,intVerts),2) > 2;

   patch('Faces',face(bulkCells,:),'Vertices',this.rV{t},'FaceColor','none','EdgeColor',color,'LineWidth',2);
   hold on
   scatter(this.rV{t}(:,1),this.rV{t}(:,2),'MarkerFaceColor',color,'MarkerEdgeColor','k');
   axis equal

end

