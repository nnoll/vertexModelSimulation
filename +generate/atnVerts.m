function [ rV ] = atnVerts( Tri, Q )
%COMPUTEVERTEXPOSITIONS 

    rV = zeros(size(Tri,1),2);
    R = [0,1;-1,0];
        
    for v = 1:size(Tri,1)

        q21 = Q(Tri(v,2),1:2) - Q(Tri(v,1),1:2);
        q32 = Q(Tri(v,3),1:2) - Q(Tri(v,2),1:2);
        q13 = Q(Tri(v,1),1:2) - Q(Tri(v,3),1:2);

        q21 = q21*R;
        q32 = q32*R;
        q13 = q13*R;

        S = abs(q21(:,1)*q32(:,2) - q21(:,2)*q32(:,1));
        pf = 1/(2*S);

        rV(v,:) = pf*(sum(Q(Tri(v,3),:).^2,2) * q21 + sum(Q(Tri(v,1),:).^2,2) * q32 + ...
                  sum(Q(Tri(v,2),:).^2,2) * q13);
              
    end

end

