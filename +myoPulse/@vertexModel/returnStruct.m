function [ Struct ] = returnStruct( this, alpha )
% RETURN STRUCT 

    d0 = this.d1{1};
    d1 = this.d0{1}';
    rV = this.rV{end};
    P = this.returnPressure();
    P = P(:,size(P,2));
    dP = this.d1{1}*P;
    
    T = this.returnTension();
    T = T(:,size(T,2));
    
    c2v = this.c2v{1};
    
    cellAdj = d0'*d0;
    cellAdj = cellAdj - diag(diag(cellAdj));
    cellAdj(cellAdj ~= 0) = 1;

    vertAdj = d1*d1';
    vertAdj = vertAdj - diag(diag(vertAdj));
    vertAdj(vertAdj ~= 0) = 1; 
    
    Struct = struct('Vdat',[],'Cdat',[],'Bdat',[]);
    
    for v = 1:size(rV,1)
        if (nargin == 1 || alpha == 0)
            Struct.Vdat(v).vertxcoord = rV(v,1);
            Struct.Vdat(v).vertycoord = rV(v,2);
        else
            Struct.Vdat(v).vertxcoord = round(alpha*rV(v,1));
            Struct.Vdat(v).vertycoord = round(alpha*rV(v,2));
        end
        Struct.Vdat(v).ncells = find(c2v(:,v));
        Struct.Vdat(v).nverts = find(vertAdj(v,:));
    end
        
    for c = 1:size(cellAdj,1)
        Struct.Cdat(c).nverts = find(c2v(c,:));
        Struct.Cdat(c).centroid.coord = mean(rV(Struct.Cdat(c).nverts,:),1);
        Struct.Cdat(c).numverts = length(Struct.Cdat(c).nverts);
        Struct.Cdat(c).ncells = find(cellAdj(c,:));
        Struct.Cdat(c).p = P(c);
    end
    
    for b = 1:size(d0,1)
        Struct.Bdat(b).cells = find(d0(b,:));
        Struct.Bdat(b).verts = find(d1(:,b));
        Struct.Bdat(b).tension = T(b);
        Struct.Bdat(b).curvature = dP(b)/T(b);
        if (length(Struct.Bdat(b).verts) == 2)
            Struct.Bdat(b).verts = [find(d1(:,b)==1);find(d1(:,b)==-1)];
            Struct.Bdat(b).length = sqrt( sum( (diff( rV(Struct.Bdat(b).verts,:),1,1)).^2, 2) );
        else
            Struct.Bdat(b).length = 0;
        end
    end

end

