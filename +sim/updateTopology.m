function [ d0_n, c2v_n ] = updateTopology( t1Bonds, t1Verts, d0, c2v )
    % UPDATE TOPOLOGY
    
    d0_n = d0;
    c2v_n = c2v;
    
    for ii = 1:length(t1Bonds)
        t1B = t1Bonds(ii);
        v1 = t1Verts(ii,1);
        v2 = t1Verts(ii,2);
        
        b1 = find(d0(:,v1));
        b2 = find(d0(:,v2));
        
        b1 = b1(b1~=t1B); 
        b2 = b2(b2~=t1B); 
        
        if (length(b1) == 2 && length(b2) == 2) % T1 in the bulk.
            t1Cells = find(c2v(:,v1) .* c2v(:,v2) == 1);
            
            c1 = find(c2v(:,v1));
            c1 = c1(~ismember(c1,t1Cells));
            c2 = find(c2v(:,v2));
            c2 = c2(~ismember(c2,t1Cells));
            
            % v1 will be attached to t1Cells(1). v2 similarly. 
            b1_1C = find(prod(c2v(:,d0(b1(1),:)~=0),2));
            b1_2C = find(prod(c2v(:,d0(b1(2),:)~=0),2));
            
            b2_1C = find(prod(c2v(:,d0(b2(1),:)~=0),2));
            b2_2C = find(prod(c2v(:,d0(b2(2),:)~=0),2));
            
            if (ismember(t1Cells(2),b1_1C)) % This bond belongs to t1Cells(2)
                d0_n(b1(1),v2) = d0_n(b1(1),v1);
                d0_n(b1(1),v1) = 0;
            elseif (ismember(t1Cells(2),b1_2C))
                d0_n(b1(2),v2) = d0_n(b1(2),v1);
                d0_n(b1(2),v1) = 0;
            end
            
            if (ismember(t1Cells(1),b2_1C)) % This bond belongs to t1Cells(2)
                d0_n(b2(1),v1) = d0_n(b2(1),v2);
                d0_n(b2(1),v2) = 0;
            elseif (ismember(t1Cells(1),b2_2C))
                d0_n(b2(2),v1) = d0_n(b2(2),v2);
                d0_n(b2(2),v2) = 0;
            end

            % Update cell data structure
            c2v_n(c1,v2) = 1;
            c2v_n(c2,v1) = 1;
            c2v_n(t1Cells(1),v2) = 0;
            c2v_n(t1Cells(2),v1) = 0;
            
        elseif (length(b1) == 1 && length(b2) == 2) % Semi-boundary edge. V1 is boundary vert

        elseif (length(b1) == 2 && length(b2) == 1) % Semi-boundary edge. V2 is boundary vert
            t1Cells = find(c2v(:,v1) .* c2v(:,v2) == 1);
            
            c1 = find(c2v(:,v1));
            c1 = c1(~ismember(c1,t1Cells));
            c2 = find(c2v(:,v2));
            c2 = c2(~ismember(c2,t1Cells));
            
            b1_1C = find(prod(c2v(:,d0(b1(1),:)~=0),2));
            b1_2C = find(prod(c2v(:,d0(b1(2),:)~=0),2));
            
            b2_C = find(prod(c2v(:,d0(b2,:)~=0),2));
            
            if (all(~ismember(b1_1C,b2_C)))
                d0_n(b1(1),v2) = d0_n(b1(1),v1);
                d0_n(b1(1),v1) = 0;
            else
                d0_n(b1(2),v2) = d0_n(b1(2),v1);
                d0_n(b1(2),v1) = 0;
            end
            
            d0_n(b2,v1) = d0_n(b2,v2);
            d0_n(b2,v2) = 0;
            
            % Update cell data structure
            v2C = find(c2v(:,v2));

            c2v_n(c1,v2) = 1;
            c2v_n(c2,v1) = 1;
            
            b2C = v2C(~ismember(v2C,b2_C));
            
            c2v_n(b2C,v2) = 0;
            c2v_n(t1Cells(~ismember(t1Cells,b2C)),v1) = 0;
        end
    end

end