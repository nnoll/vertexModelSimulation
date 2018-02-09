function [ extVerts ] = findExtVerts( this )
%FINDBOUNDARYVERTS 

    bV = this.bV{end};
    D = abs(this.d0{end});
    bE = sum(D(:,bV),2) >= 1;
    
    extVerts = find(sum(D(bE,:),1) >= 1);
    
    Adj = this.d0{end}'*this.d0{end};
    Adj = Adj - diag(diag(Adj));
    
    % Order external verts
    orderedEVerts = zeros(size(extVerts));
    orderedEVerts(1) = extVerts(1);
    for ii = 2:length(extVerts)
        nVerts = find(Adj(orderedEVerts(ii-1),:));
        nVerts = nVerts(ismember(nVerts,extVerts));
        nVerts = nVerts(~ismember(nVerts,orderedEVerts));
        [~,ind] = min(this.rV{end}(nVerts,1));
        orderedEVerts(ii) = nVerts(ind);
    end
    
    extVerts = orderedEVerts;

end

