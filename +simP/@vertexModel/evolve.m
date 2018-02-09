function [ this ] = evolve( this, Tend, deltaT )
    % EVOLVE 
    
    if (nargin == 2)
        deltaT = 1;
    end
    
    Tfinal = (length(this.rV)+Tend-1);
    tspan = (length(this.rV)-1):deltaT:Tfinal;
        
    r0 = this.rV{end};
    X0 =  r0(:);
    Te = 1;
    while (~isempty(Te))
        
        if (this.Lc > 0)
            options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@(t,X) simP.T1event(t,X,0,sparse(this.d0{end}),this.Lc,this.bulkE{end}));
        else
            options = odeset('RelTol',1e-5,'AbsTol',1e-7);
        end

        EOM = @(t,X) simP.equationOfMotion( t, X, this.d0{end}, this.d0{end}', this.d1{end}, this.l0, this.kappaL, this.gamma,...
                                           this.cTopo{end}, this.a0, this.Pb, this.bV{end}, this.bC{end}, this.Fbc );
                                       
        if (this.Lc > 0)                               
            [t,X,Te,Xe] = ode15s(EOM,tspan,X0,options);
        else
            [t,X] = ode15s(EOM,tspan,X0,options);
            Te = [];
        end
        
        if (~isempty(Te))
            this.TE{end+1} = Te;
            
            nBonds = size(this.d0{end},1);
            nVerts = size(this.d0{end},2);
            
            % Store integration so far.
            X = X';
            Ntps = length(this.T);

            for ii = 2:(size(X,2)-1)
                this.rV{Ntps+ii-1} = reshape(X(:,ii),nVerts,2);
            end

            % Perform T1 and reboot integration.
            Xe = Xe';

            d0 = this.d0{end};
            c2v = this.c2v{end};

            TE = Xe(1:nBonds);
            mE = Xe((nBonds+1):2*nBonds);
            rVE = reshape(Xe((2*nBonds+1):length(Xe)),nVerts,2);
            rBE = d0*rVE; 
            rPB = .5*abs(d0)*rVE;
            lB = sqrt(sum(rBE.^2,2));

            rV0 = rVE;
            T0 = TE;
            M0 = mE;

            t1Bonds = find((lB < (this.Lc + .001)) & (sqrt(sum(rPB.^2,2)) <= 6) ); %Adjust for floating point error.
            t1Bonds = t1Bonds(ismember(t1Bonds,this.bulkE{end}));
            
            t1Verts = zeros(length(t1Bonds),2);
            for bb = 1:length(t1Bonds)
                t1Verts(bb,:) = find(d0(t1Bonds(bb),:));
            end
 
            % Update graph topology
            [ d0_n, c2v_n ] = sim.updateTopology(t1Bonds, t1Verts, d0, c2v);

            N = length(this.d0);
            this.d0{N+1} = d0_n;
            this.c2v{N+1} = c2v_n;
            this.bulkE{N+1} = sim.returnIntBonds(d0_n);
       
            % Update graph geometry
            for bb = 1:length(t1Bonds)
                t1Cells = find(c2v(:,t1Verts(bb,1)) .* c2v(:,t1Verts(bb,1)) == 1);
                R0 = .5*sum(rVE(t1Verts(bb,:),:),1);
                deltaR = (mean(rVE(c2v(t1Cells(1),:)==1,:),1) - mean(rVE(c2v(t1Cells(2),:)==1,:),1));
                deltaR = (this.Lc + .1) * (deltaR / sqrt(sum(deltaR.^2,2)));
                rV0(t1Verts(bb,1),:) = R0 + .5*deltaR;
                rV0(t1Verts(bb,2),:) = R0 - .5*deltaR;
            end

            [ this.Tri{N+1} ] = generate.triangulation( rV0, this.c2v{N+1}, find(sum(abs(this.d0{N+1}),1)==3) );

            r0 = rV0(:);

            X0 = r0(:);
            tspan = [Te,ceil(Te):1:Tfinal];
            
        else     
            X = X';
            nVerts = size(this.d0{end},2);
            Ntps = length(this.rV);
            for ii = 2:size(X,2)
                this.rV{Ntps+ii-1} = reshape(X(:,ii),nVerts,2);
            end
            
        end
        
    end
    
end

