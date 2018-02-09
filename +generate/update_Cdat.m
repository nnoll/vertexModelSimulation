function [ Cdat ] = update_Cdat( Cdat, c2v, rv )
%UPDATE_CDAT 

for c=1:size(Cdat,1)
    ind = find(c2v(c,:));
    Cdat(c,1) = mean(rv(1,ind));
    Cdat(c,2) = mean(rv(2,ind));
end

end

