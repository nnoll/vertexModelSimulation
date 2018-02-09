function [ nadj ] = clean_adj( adj )
%CLEAN_ADJ Cleans up the boundaries of generated lattices.

nadj = adj;
nnverts = sum(adj,1);
two_verts = find(nnverts == 2);
for v = two_verts
    conn_verts = find(adj(v,:));
    two_conn_verts = conn_verts(ismember(conn_verts,two_verts));
    if (~isempty(two_conn_verts))
        for vv = two_conn_verts
            nadj(v,vv) = 0;
            nadj(vv,v) = 0;
        end
    end
end


end

