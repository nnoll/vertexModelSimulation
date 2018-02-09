function [ rot_sign, adj, LUT, adj_b ] = change_topology( adj, LUT, T1_bond, bverts, rv, D, adj_b )
%UPDATE_TOPOLOGY

rb = LUT*rv';
n = [0, -1; 1, 0]*rb(T1_bond,:)'/D(T1_bond);
cbonds = find(adj_b(T1_bond,:));
nverts = zeros(2,2);
r1 = zeros(2,2); r2 = zeros(2,2);

%Bulk T1s only.
v1 = bverts(1); v2 = bverts(2); 
tmp = find(adj(v1,:));
nverts(1,:) = tmp(~ismembc(tmp,v2));
tmp = find(adj(v2,:));
nverts(2,:) = tmp(~ismembc(tmp,v1));
for ii=1:2
    r1(ii,:) = rv(ii,nverts(1,:)) - rv(ii,v1)*ones(1,2);
    r2(ii,:) = rv(ii,nverts(2,:)) - rv(ii,v2)*ones(1,2);
end
[pair] = find_pair(r1, r2, n);

if (pair == 1)
    tmp = nverts;
    tmp(1,2) = nverts(2,1);
    tmp(2,1) = nverts(1,2);
    nverts = tmp;
else
    tmp = nverts;
    tmp(1,2) = nverts(2,2);
    tmp(2,2) = nverts(1,2);
    nverts = tmp;
end   

%Vertex adjacency.
n_v1 = find(adj(v1,:)); n_v2 = find(adj(v2,:));
adj(v1,:) = 0; adj(:,v1) = 0;
adj(v2,:) = 0; adj(:,v2) = 0;
adj(n_v1,v1) = 0; adj(n_v2,v2) = 0;
adj(v1,v2) = 1; adj(v2,v1) = 1;
adj(v1,nverts(1,:)) = 1; adj(nverts(1,:),v1) = 1;
adj(v2,nverts(2,:)) = 1; adj(nverts(2,:),v2) = 1;

if (LUT(T1_bond,v1) == 1)
    sn = sign(n'*(rv(:,nverts(1,1))+rv(:,nverts(1,2))-2*rv(:,v1)));
else
    sn = sign(n'*(rv(:,nverts(2,1))+rv(:,nverts(2,2))-2*rv(:,v2)));
end
rot_sign = sn;

%Bond/vertex LUT.
[LUT] = update_LUT(v1,v2,LUT,cbonds,adj);
[adj_b] = generate.create_adjB(LUT,size(LUT,2),size(LUT,1));

end

function [pairing] = find_pair(r1, r2, n)
    delta1 = r1(:,1) + r2(:,1); delta2 = r1(:,2) + r2(:,2);
    delta3 = r1(:,1) + r2(:,2); delta4 = r1(:,2) + r2(:,1);
    d1 = (n'*delta1)^2 + (n'*delta2)^2;
    d2 = (n'*delta3)^2 + (n'*delta4)^2;
    if (d1 > d2)
        pairing = 1;
    else
        pairing = 2;
    end
end

function [LUT] = update_LUT(v1,v2,LUT,cbonds,adj)
    for cb=cbonds
        old_verts = find(LUT(cb,:));
        still_connected = adj(old_verts(1),old_verts(2));
        if (~still_connected) %Then we switched one of the vertices.
            if (ismember(v1,old_verts))
                LUT(cb,v2) = LUT(cb,v1);
                LUT(cb,v1) = 0;
            else
                LUT(cb,v1) = LUT(cb,v2);
                LUT(cb,v2) = 0;
            end
        end
    end
end

