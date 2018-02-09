function [ mov ] = plotTensedLattice( this )
    % PLOTLATTICE Summary of this function goes here
    %   Detailed explanation goes here

    Tt = this.returnTension();
    cmap = jet(256);
    fig = figure('units','normalized','outerposition',[0 0 1 1]);

    for t = 1:size(Tt,2)
        clf
        T = Tt(:,t);
        T = (T - min(T))/(max(T)-min(T));
        T = round(255*T+1);
        for e = 1:size(this.d0{1})
           edgeVerts = (this.d0{1}(e,:) ~= 0);
           patch('Faces',[1,2],'Vertices',this.rV{t}(edgeVerts,:),'FaceColor','k','EdgeColor',cmap(T(e),:),'LineWidth',2)
        end
        colorbar()
        mov(t) = getframe();
    end


end

