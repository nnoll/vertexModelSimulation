function [ ZI, xG ] = interpolate1D( Z )
%INTERPOLATEGRID 

%     Z = [Z,Z(:,1)];

    x = linspace(-.5,.5,size(Z,2));
    xG = linspace(-.5,.5,500);
        
    ZI = interp1(x,Z',xG,'spline')';

end

