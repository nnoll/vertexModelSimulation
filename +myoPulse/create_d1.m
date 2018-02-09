function [ d1 ] = create_d1( d0, rv, c2v )

    d1 = zeros(size(d0,1),size(c2v,1));
    rv = rv';
    for bb = 1:size(d0,1) 

        verts = find(d0(bb,:));
        v1c = c2v(:,verts(1));
        v2c = c2v(:,verts(2));
        vc = v1c .* v2c;
        cells = find(vc);
        rb = d0(bb,:)*rv';
        rb = [0,1;-1,0]*rb';

        Rc = mean(rv(:,c2v(cells(1),:)==1),2) - mean(rv(:,c2v(cells(2),:)==1),2);

        d1(bb,cells(1)) = sign(Rc'*rb);
        d1(bb,cells(2)) = -1*sign(Rc'*rb);
        
    end

end

