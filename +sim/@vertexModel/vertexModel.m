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
        bulkE    % Bulk Edges
        q0    % Initial triangulation points.
        
        d0    % Adjacency matrix. 
        c2v   % Cell/Vertex List.
        Tri   % Triangulation array
        
        T     % Tensions
        m     % Myosins
        
        nu    % Walking contractility time-scale
        omega % Myosin feedback time-scale.
        Lc    % Critical length under which you T1
        
        TE    % Time of T1 event.
        
        Fbc   % Boundary conditions
        Fext  % External Force.
        P     % Bulk Pressure inside tissue.
        
        mode  % Shear or Pull simulation.
        
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = vertexModel( rV, bV, q0, d0, c2v, Tri, T, m, nu, omega, Fbc, Fext, P, Lc, mode )

            this.rV{1} = rV;
            order = sim.orderVerts(rV(bV,:));
            this.bV{1} = bV(order);
            this.q0 = q0;
            
            this.d0{1} = sparse(d0);
            this.c2v{1} = sparse(c2v);
            this.Tri{1} = Tri;
            
            this.bulkE{1} = sim.returnIntBonds(d0);
            this.T{1} = T;
            this.m{1} = m;
            
            this.nu = nu;
            this.omega = omega;
            this.Lc = Lc;
            
            this.TE{1} = 0;
            
            this.Fbc = Fbc(order,:);
            this.Fext = Fext;
            this.P = P;
            
            this.mode = mode;
        end
        
    end
               
end


