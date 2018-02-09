function [ strain ] = returnStrain( this )
    % RETURN MODULUS 
    
    % Fit boundary shape change to a 2D transformation.
    [ l, r, u, d ] = sim.partitionBoundary( this.rV{1}, this.bV{1}, 1 );

    if (strcmp(this.mode,'pull'))
        Llong = mean(this.rV{1}(this.bV{1}(u),2)) - mean(this.rV{1}(this.bV{1}(d),2));
        LlongNew = mean(this.rV{length(this.rV)}(this.bV{1}(u),2)) - mean(this.rV{length(this.rV)}(this.bV{1}(d),2));
        Ltrans = mean(this.rV{1}(this.bV{1}(r),1)) - mean(this.rV{1}(this.bV{1}(l),1));
        LtransNew = mean(this.rV{length(this.rV)}(this.bV{1}(r),1)) - mean(this.rV{length(this.rV)}(this.bV{1}(l),1));
        strain(1) = (LlongNew - Llong)/Llong;
        strain(2) = (LtransNew - Ltrans)/Llong;
    elseif (strcmp(this.mode,'shear'))
        L = mean(this.rV{1}(this.bV{1}(u),2)) - mean(this.rV{1}(this.bV{1}(d),2));
        delta = .5*((max(this.rV{length(this.rV)}(this.bV{1}(l),1)) - min(this.rV{length(this.rV)}(this.bV{1}(l),1))) ...
                  + (max(this.rV{length(this.rV)}(this.bV{1}(r),1)) - min(this.rV{length(this.rV)}(this.bV{1}(r),1))) );
        strain = delta/L;
    end


end

