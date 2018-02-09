function [ XG, YG, ZG ] = interpolateGrid( x, y, Z )
%INTERPOLATEGRID 

    [X,Y] = ndgrid(x,y);
    F = griddedInterpolant(X,Y,Z,'cubic');
    [XG,YG] = ndgrid(linspace(min(x),max(x),100),linspace(min(y),max(y),100));

    ZG = F(XG,YG);

end

