function [ Cdat, c2v, cadj ] = get_cells( rv )
%GET_CELLS Identifies cells within the simulated lattice. Outputs a list of
%vertices belonging to each cell.

%%Inputs 
%1. rv - simulated vertex positions

%%Outputs
%1. C - cell data structure
[C,V] = voronoin(rv'); 

%Build cell to vertex map.
c2v = zeros(size(C,1)-1,size(rv,2));
for v=1:length(V)
    for ii=1:length(V{v})
        if (V{v}(ii) > 1)
            c2v(V{v}(ii)-1,v) = 1;
        end
    end
end

%Build cell data structure.
Cdat = zeros(size(C,1)-1,3);
for c=1:size(Cdat,1)
    ind = find(c2v(c,:));
    Cdat(c,1) = mean(rv(1,ind));
    Cdat(c,2) = mean(rv(2,ind));
end

%Compute cell adjacency matrix.
cadj = (c2v*c2v')/2;
for n=1:size(cadj,1) %Remove diagonal pieces - just counts vertices cell has.
    cadj(n,n) = 0;
end

end

