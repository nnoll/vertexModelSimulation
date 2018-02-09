function [ rv, adj ] = seedLattice( N, L )
%GENERATE_LATTICE Generates a regular hexagonal lattice.

%Inputs:
%N = number of bulk cells to create
%L = length of edge.

%Outputs:
%rv = vertex positions
%adj = adjacency matrix

%Basis vectors for each cell centroid.
b1 = [0;sqrt(3)*L];
b2 = sqrt(3)*L/2*[sqrt(3);-1];

%Basis vectors for a cell's vertex positions.
R60 = .5*[1, -sqrt(3); sqrt(3), 1];
bv = zeros(2,6);
bv(:,1) = L*[1;0];
bv(:,2) = R60*bv(:,1);
bv(:,3) = R60*bv(:,2);
bv(:,4:6) = -bv(:,1:3);

%Generate cell centroid positions (not most elegant generation).
rc = zeros(2,7*N);
num_layers = floor(-.5 + (1/6)*sqrt(9+12*(N-1))) + 1;
delta = zeros(2,6);
delta(:,1) = b2; delta(:,2) = -b1; delta(:,3) = -(b1+b2);
delta(:,4:6) = -delta(:,1:3);

for nn=1:num_layers
    N_pc = 3*(nn-1)*(nn-2) + 1; %Number of past cells seen
    if (nn == 1)
        rc(:,1) = [0;0];
        rc(:,2) = b1;
        rc(:,3) = -b1;
        rc(:,4) = b2;
        rc(:,5) = -b2;
        rc(:,6) = +(b1+b2);
        rc(:,7) = -(b1+b2);
    else
        rr = 1;
        dd = 1;
        for cc=1:6*(nn-1)
            if (cc==1)
                rc(:,7*N_pc+7*(cc-1)+1) = (nn-1)*b1;
            else
                rc(:,7*N_pc+7*(cc-1)+1) = rc(:,7*N_pc+7*(cc-2)+1) + delta(:,dd);
                rr = rr + 1;
                if (rr > (nn-1))
                    rr = 1;
                    dd = dd + 1;
                end
            end
            rc(:,7*N_pc+7*(cc-1)+2) = rc(:,7*N_pc+7*(cc-1)+1) + b1;
            rc(:,7*N_pc+7*(cc-1)+3) = rc(:,7*N_pc+7*(cc-1)+1) - b1;
            rc(:,7*N_pc+7*(cc-1)+4) = rc(:,7*N_pc+7*(cc-1)+1) + b2;
            rc(:,7*N_pc+7*(cc-1)+5) = rc(:,7*N_pc+7*(cc-1)+1) - b2;
            rc(:,7*N_pc+7*(cc-1)+6) = rc(:,7*N_pc+7*(cc-1)+1) + (b1+b2);
            rc(:,7*N_pc+7*(cc-1)+7) = rc(:,7*N_pc+7*(cc-1)+1) - (b1+b2);
        end
    end
end

if (N > 3*(num_layers-1)*(num_layers) + 1) %Add cells at border.
    delta_N = N - (3*(num_layers-1)*(num_layers) + 1);
    nn = num_layers;
    N_pc = 3*nn*(nn-1) + 1;
    cells = randsample(6*nn, delta_N);
    rr = 1;
    dd = 1;
    iter = 0;
    for cc=1:6*nn
        if (cc==1)
            rc0 = nn*b1;
        else
            rc0 = rc0 + delta(:,dd);
            rr = rr + 1;
            if (rr > nn)
                rr = 1;
                dd = dd + 1;
            end
        end
        if (ismember(cc,cells))
            rc(:,7*N_pc+7*(iter)+1) = rc0;
            rc(:,7*N_pc+7*(iter)+2) = rc(:,7*N_pc+7*(iter)+1) + b1;
            rc(:,7*N_pc+7*(iter)+3) = rc(:,7*N_pc+7*(iter)+1) - b1;
            rc(:,7*N_pc+7*(iter)+4) = rc(:,7*N_pc+7*(iter)+1) + b2;
            rc(:,7*N_pc+7*(iter)+5) = rc(:,7*N_pc+7*(iter)+1) - b2;
            rc(:,7*N_pc+7*(iter)+6) = rc(:,7*N_pc+7*(iter)+1) + (b1+b2);
            rc(:,7*N_pc+7*(iter)+7) = rc(:,7*N_pc+7*(iter)+1) - (b1+b2);
            iter = iter + 1;
        end
    end
    num_layers = num_layers + 1;
end

%Remove redudancies in cell positions.
rc_tmp = round(rc*1e2)/1e2;
[~,ia,~] = unique(rc_tmp','rows');
rc = rc(:,ia);
%Generate vertex positions.
rv = zeros(2,6*size(rc,2));
for ii=1:size(rc,2)
    rv(:,6*(ii-1)+1:6*(ii)) = rc(:,ii)*ones(1,6) + bv;
end

%Remove redudancies in vertex positions.
rv_tmp = round(rv*1e2)/1e2;
[~,ia,~] = unique(rv_tmp','rows');
rv = rv(:,ia);

%Find adjacency matrix of our graph.
adj = zeros(size(rv,2));
for mm=1:size(rv,2)
    nnverts = 0;
    for nn=(mm+1):size(rv,2)
        deltav = rv(:,mm) - rv(:,nn);
        if (norm(deltav) <= L + .5 && norm(deltav) >= L - .5)
            adj(nn,mm) = 1;
            adj(mm,nn) = 1;
            nnverts = nnverts + 1;
        end
    end
end

% %Code for testing purposes
% rc = rc + [(2*(num_layers)+1)*L;(2*(num_layers)+1)*L] * ones(1,size(rc,2));
% rv = rv + [(2*(num_layers)+1)*L;(2*(num_layers)+1)*L] * ones(1,size(rv,2));
% tmp = zeros(2*(2*(num_layers)+1)*L);
% tmp2 = zeros(2*(2*(num_layers)+1)*L);
% 
% for ii=1:size(rc,2)
%     tmp(round(rc(2,ii)),round(rc(1,ii))) = 1;
% end
% 
% for ii=1:size(rv,2)
%     tmp2(round(rv(2,ii)),round(rv(1,ii))) = 1;
% end
% 
% rgb(:,:,1) = tmp;
% rgb(:,:,3) = tmp2;
% imshow(rgb)
% pause
    
end

