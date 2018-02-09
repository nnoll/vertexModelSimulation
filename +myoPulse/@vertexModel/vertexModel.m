classdef vertexModel < handle
    %   Generate a primal structure to be used for discrete simulation of 
    %   our active vertex model.
    %
    
    %---------------------------------------------------------------------
    % Properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        rV     % Vertex Positions.
        T      % Tensions
        m      % Myosins
        
        bV     % Boundary verts
        bC     % Boundary cells
        bulkE  % Bulk Edges
        
        d0     % Adjacency matrix. 
        d1     % Bond matrix.
        c2v    % Cell/Vertex List.
        cTopo  % Used to compute areas.
        
        tauL   % Contractility time-scale
        alpha  % Myosin recruitment time-scale
        
        Pb     % Boundary pressure
        Kappa  % Area Elasticity Strength
        A0     % Desired size for each cell
        
        Lc     % Critical length under which you T1
        
%         tPulse % Times of pulses for each cell
        omegaPulse
        phiPulse
        A      % Amplitude of pulses
%         lambda % Decay of pulses
        
        TE     % Time of T1 event.
        Fbc    % Boundary conditions        
        
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = vertexModel( rV, T, m, bV, bC, d0, c2v, Pb, tauL, alpha, A, Tpulse, sigma, Kappa, A0, yMax, Fbc, Lc )

            this.rV{1} = rV;
            this.T{1} = T;
            this.m{1} = m;
            
            order = sim.orderVerts(rV(bV,:));
            this.bV{1} = bV(order);
            this.bC{1} = bC;
            
            this.Pb = Pb;
            this.tauL = tauL;
            this.alpha = alpha;
            this.Kappa = Kappa;
            this.A0 = A0;
            
            this.d0{1} = sparse(d0);
            this.c2v{1} = sparse(c2v);
            
            adj = d0'*d0;
            adj = adj - diag(diag(adj));
            adj = -adj;
            
            [ this.cTopo{1} ] = simP.partition_cellTopo( adj, c2v, this.bC{1} );
            
            this.A = A;
%             this.lambda = lambda;
%             Tfinal = 1000;
%             this.tPulse = myoPulse.generatePulseTrain(beta*lambda,gamma*lambda,size(c2v,1),Tfinal);
            [this.omegaPulse,this.phiPulse] = myoPulse.generateRandomFreq(Tpulse,sigma,size(c2v,1));
            
            % Silence cells about yMax
            rC = bsxfun(@rdivide,this.c2v{1}*this.rV{1},sum(this.c2v{1},2));
%             this.tPulse(abs(rC(:,2))>yMax,:) = Tfinal;
%             this.tPulse(bC,:) = Tfinal;
            this.omegaPulse(abs(rC(:,2))>yMax) = 0;
            this.phiPulse(abs(rC(:,2))>yMax) = 0;
            
            [ this.d1{1} ] = sparse(simP.create_d1( d0, rV, c2v ));
            
            this.bulkE{1} = sim.returnIntBonds(d0);

            this.Lc = Lc;
            this.TE{1} = 0;
            
            this.Fbc = Fbc(order,:); 
            
        end
        
    end
               
end


