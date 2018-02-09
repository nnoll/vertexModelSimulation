function [ ] = plotTensedLattice( this, t )
    % PLOTLATTICE Summary of this function goes here
    %   Detailed explanation goes here

    clf
    T = this.T{t};
    cmap = jet(256);
    T = (T - min(T))/(max(T)-min(T));
    T = round(255*T+1);
    for e = 1:size(this.d0{1})
       edgeVerts = (this.d0{1}(e,:) ~= 0);
       patch('Faces',[1,2],'Vertices',this.rV{t}(edgeVerts,:),'FaceColor','k','EdgeColor',cmap(T(e),:),'LineWidth',2)
    end


end

