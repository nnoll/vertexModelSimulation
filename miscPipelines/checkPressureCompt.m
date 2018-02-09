function [ chi, chiEst, P ] = checkPressureCompt( Struct )
    % BUILD CONTROL DISTRIBUTION 
            
    [ ~, bulkCells, ~, bulkVerts ] = fitDual.returnGraph( Struct, 0 ); 

    chi = zeros(length(bulkCells),1);
    chiEst = zeros(length(bulkCells),1);
    n = 1;

    bVerts = [Struct.Bdat.verts];
    kappa = [Struct.Bdat.curvature] .* [Struct.Bdat.length];
   
    for c = bulkCells
        cverts = Struct.Cdat(c).nverts;

        % Order verts
        rC = double([Struct.Vdat(cverts).vertxcoord;Struct.Vdat(cverts).vertycoord])';
        deltaRC = bsxfun(@minus,rC,mean(rC,1));

        theta = mod(atan2(deltaRC(:,2),deltaRC(:,1)),2*pi);
        [~,ind] = sort(theta);
        cverts = cverts(ind);
        rC = rC(ind,:);

        for ii = 1:length(cverts)

            p =  mod(ii,length(cverts)) + 1;
            m = mod(ii-2,length(cverts)) + 1;

            bondP = find( ((bVerts(1,:) == cverts(ii)) & (bVerts(2,:) == cverts(p))) | ...
                          ((bVerts(1,:) == cverts(p)) & (bVerts(2,:) == cverts(ii))) );
            if (bVerts(1,bondP) == cverts(ii))
                Sp = -1;
            else
                Sp = 1;
            end
            bondM = find( ((bVerts(1,:) == cverts(ii)) & (bVerts(2,:) == cverts(m))) | ...
                          ((bVerts(1,:) == cverts(m)) & (bVerts(2,:) == cverts(ii))) );
            if (bVerts(1,bondM) == cverts(ii))
                Sm = -1;
            else
                Sm = 1;
            end              
            delta1 = rC(p,:) - rC(ii,:); 
            delta1 = delta1/sqrt(sum(delta1.^2));
            delta2 = rC(m,:) - rC(ii,:); 
            delta2 = delta2/sqrt(sum(delta2.^2));

            vE = Struct.Vdat(cverts(ii)).nverts(~ismember(Struct.Vdat(cverts(ii)).nverts,cverts));

            bondE = find( ((bVerts(1,:) == cverts(ii)) & (bVerts(2,:) == vE)) | ...
                          ((bVerts(1,:) == vE) & (bVerts(2,:) == cverts(ii))) );
            if (bVerts(1,bondE) == cverts(ii))
                Se = -1;
            else
                Se = 1;
            end          
            deltaE = double([Struct.Vdat(vE).vertxcoord,Struct.Vdat(vE).vertycoord]);
            deltaE = deltaE - rC(ii,:); 
            deltaE = deltaE/sqrt(sum(deltaE.^2));

            angle1 = acos(dot(deltaE,delta1));
            angle2 = acos(dot(deltaE,delta2));
            
            chi(n) = chi(n) + log(sin(angle1)) - log(sin(angle2));
            if (~isempty(bondE) && ~isempty(bondP) && ~isempty(bondM))
                chiEst(n) = chiEst(n) + cot(angle1)*(Sp*kappa(bondP)-Se*kappa(bondE)) - cot(angle2)*(Se*kappa(bondE)-Sm*kappa(bondM));
            end
            
        end
        P(n) = Struct.Cdat(c).p;
        n = n + 1;

    end

end

