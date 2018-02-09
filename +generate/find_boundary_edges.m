function [ bdry_edge, rv, adj ] = find_boundary_edges( adj, rv )
%FIND_BOUNDARY_EDGES Find and return distinction between boundary edges and
%bulk edges. Will prune the vertices to only fully connected.

N = size(adj,1);
bdry_edge = zeros(size(adj));

nnverts = sum(adj,1);
three_vert = find(nnverts==3);
bad_vert = find(nnverts==2);
ext_vert = [];
ext_nvert = [];

for n=1:N
    if (nnverts(n) == 3)
        conn_verts = find(adj(n,:));
        conn_2zverts = conn_verts(~ismember(conn_verts,three_vert));
        if (length(conn_2zverts)==1)
            ext_vert = horzcat(ext_vert,conn_2zverts);
            ext_nvert = horzcat(ext_nvert,n);
        elseif (length(conn_2zverts)==2)
            ext_vert = horzcat(ext_vert,n);
            ext_nvert = horzcat(ext_nvert,conn_verts(~ismember(conn_verts,conn_2zverts)));
        end
    end
end

for n=1:length(ext_vert)
    bdry_edge(ext_vert(n),ext_nvert(n)) = 1;
    bdry_edge(ext_nvert(n),ext_vert(n)) = 1;
end

bad_vert(ismember(bad_vert,ext_vert)) = [];

rv(:,bad_vert) = [];
adj(:,bad_vert) = [];
adj(bad_vert,:) = [];
bdry_edge(:,bad_vert) = [];
bdry_edge(bad_vert,:) = [];

end

