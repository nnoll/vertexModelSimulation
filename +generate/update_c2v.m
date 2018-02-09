function [ c2v ] = update_c2v( c2v, adj, verts )
%UPDATE_C2V If a T1 process occurs, not only has vertex topology changed
%but also vertices identification with cells is changed. This must be taken
%into account.

%%Inputs 
% 1. c2v = List of vertices for each cell before the T1
% 2. adj = vertex adjacency after T1
% 3. verts = verts that underwent T1

%%Outputs
% 1. c2v = List of vertices for each cell after the T1

%Find involved cells. Partition into cells that will get new vertex and
%those that lost one.
vcell1 = find(c2v(:,verts(1)));
vcell2 = find(c2v(:,verts(2)));
involved_cells = union(vcell1,vcell2);
loser_cells = vcell1(ismember(vcell1,vcell2));
gainer_cells = involved_cells(~ismember(involved_cells,loser_cells));

%Gainer cells are easy. Just put a 1 in both columns.
for c = gainer_cells
    c2v(c,verts(1)) = 1;
    c2v(c,verts(2)) = 1;
end

%Loser cells require a bit more planning.
for ii=1:2
    v1 = verts(ii);
    v2 = verts(mod(ii,2)+1);
    conn_verts = find(adj(v1,:));
    conn_verts(conn_verts == v2) = []; 
    if (length(conn_verts) >= 1)
        ncells = find(c2v(:,conn_verts(1)));
    %     ncells2 = find(c2v(:,conn_verts(2)));
    %     ncells = ncells1(ismember(ncells1,ncells2));
        ncells = ncells(ismember(ncells,loser_cells));
        c2v(ncells,v2) = 0;
    else
        if (length(find(c2v(loser_cells(1),:))) <= 2)
            c2v(loser_cells(1),v2) = 0;
        else
            c2v(loser_cells(2),v2) = 0;
        end
    end
end

end

