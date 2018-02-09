function [ p ] = returnAreas( this, yMax )
%PLOTAREAS Summary of this function goes here
%   Detailed explanation goes here

    p = zeros(size(this.d1{1},2),length(this.rV));
    cTopo = this.cTopo{1};
    
    for t = 1:length(this.rV)
        rv = this.rV{t};
        
        r1 = squeeze(rv(:,1));
        r2 = squeeze(rv(:,2));
        for ii=1:length(cTopo)
            z = size(cTopo(ii).order,2);
            G = diag(ones(z-1,1),1) - diag(ones(z-1,1),-1);
            G(1,z) = -1; G(z,1) = 1;

            rcv1 = r1(cTopo(ii).order')';
            rcv2 = r2(cTopo(ii).order')';
            cArea = sum((rcv1*G).*rcv2,2);
            p(cTopo(ii).clist,t) = .5*(abs(cArea));
        end
        
    end
        
    rC = bsxfun(@rdivide,this.c2v{1}*this.rV{1},sum(this.c2v{1},2));
    ind = abs(rC(:,2))<5 & (~ismember(1:size(rC,1),this.bC{1}))';
    p = p(ind,:);
    
end

