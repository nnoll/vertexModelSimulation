function [ cTopo ] = partition_cellTopo( adj, c2v, bndryCells )
    %PARTITION_CELL TOPO Takes in the list of vertices of each cell and returns
    % a data structure containing an ordered list for each vert partitioned into
    % cells of similar foldedness.

    Z = sum(c2v,2);
    [zlist,~,zcells] = unique(Z);
    mm = 1;
    for nn = 1:length(zlist)
        clist = find(zcells == nn);
        clist = clist(~ismember(clist,bndryCells));
        if (~isempty(clist))
            cTopo(mm).clist = clist;
            cTopo(mm).order = zeros(length(clist),zlist(nn));  
            for cc = 1:length(clist)
                c = clist(cc);
                cverts = find(c2v(c,:));
                oCverts = zeros(size(cverts));
                oCverts(1) = cverts(1);
                for ii = 2:zlist(nn)
                    conVerts = find(adj(:,oCverts(ii-1)));
                    conVerts = conVerts(ismember(conVerts,cverts));
                    conVerts = conVerts(~ismember(conVerts,oCverts));
                    oCverts(ii) = conVerts(1);
                end
                cTopo(mm).order(cc,:) = oCverts;
            end
            mm = mm + 1;
        end
    end

end

