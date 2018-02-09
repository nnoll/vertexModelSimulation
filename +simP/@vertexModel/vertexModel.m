classdef vertexModel < handle
    %   Generate a primal structure to be used for discrete simulation of 
    %   our active vertex model.
    %
    
    %---------------------------------------------------------------------
    % Properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        rV    % Vertex Positions.
        bV    % Boundary verts
        bC    % Boundary cells
        bulkE % Bulk Edges
        
        d0    % Adjacency matrix. 
        d1    % Bond matrix.
        c2v   % Cell/Vertex List.
        cTopo % Used to compute areas.
        
        gamma % Area stiffness
        a0    % Target Area
        Pb    % Boundary pressure
        Lc    % Critical length under which you T1
        l0    % Target Length
        kappaL
        
        TE    % Time of T1 event.
        
        Fbc   % Boundary conditions        
        
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = vertexModel( rV, bV, bC, d0, c2v, gamma, a0, Pb, Fbc, l0, kappaL, Lc )

            this.rV{1} = rV;
            order = sim.orderVerts(rV(bV,:));
            this.bV{1} = bV(order);
            this.bC{1} = bC;
            
            this.gamma = gamma;
            this.a0 = a0;
            this.Pb = Pb;
            this.l0 = l0;
            this.kappaL = kappaL;
            
            this.d0{1} = sparse(d0);
            this.c2v{1} = sparse(c2v);
            
            adj = d0'*d0;
            adj = adj - diag(diag(adj));
            adj = -adj;
            
            [ this.cTopo{1} ] = simP.partition_cellTopo( adj, c2v, this.bC{1} );
            [ this.d1{1} ] = sparse(simP.create_d1( d0, rV, c2v ));
            
            this.bulkE{1} = sim.returnIntBonds(d0);

            this.Lc = Lc;
            this.TE{1} = 0;
            
            this.Fbc = Fbc(order,:);      
        end
        
    end
               
end


