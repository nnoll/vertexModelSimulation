function [adj_b] = create_adjB(LUT,N,nbonds)
    adj_b = zeros(nbonds);
    for n=1:N
        shared_bonds = find(LUT(:,n));
        for sb=shared_bonds'
            for sbb=shared_bonds'
                if (sb == sbb)
                    adj_b(sb,sbb) = 1;
                else
                    if (LUT(sb,n)*LUT(sbb,n) == 1)
                        adj_b(sb,sbb) = -1;
                    else
                        adj_b(sb,sbb) = 1;
                    end
                end
            end
        end
    end
    adj_b = adj_b - eye(nbonds);
end
