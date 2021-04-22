function [y, sy] = kinematics(x, sx, L, opt);
%
% forward kinematics and inverse kinematics for full Bayesian
% extension of Sober and Sabes model

if opt == 1, % forward kinematics: angle --> metric
    Rmu = [-sin(x(1))*sin(x(2)) -sin(x(1)); cos(x(1))*sin(x(2)) cos(x(1))];
    y = Rmu*L;
    [e1,e2] = size(sx);
    if e1 ~= e2, sx = diag(sx); end
    A = [-cos(x(1))*sin(x(2))*L(1)-cos(x(1))*L(2)  -sin(x(1))*cos(x(2))*L(1); ...
        -sin(x(1))*sin(x(2))*L(1)-sin(x(1))*L(2)  cos(x(1))*cos(x(2))*L(1)]+eps;
    sy = abs(det(A+eps))*A*(sx*inv(A+eps));
elseif opt == -1, % inverse kinematics: metric --> angle
    y(1) = atan(-x(1)/x(2));
    temp = (sqrt(x(1)^2+x(2)^2)-L(2))/L(1);
    if temp > 0, temp = rem(temp+1,2)-1; else temp = rem(temp-1,2)+1; end
    y(2) = asin(temp);
    [e1,e2] = size(sx);
    if e1 ~= e2, sx = diag(sx); end
    c = sqrt(1-((sqrt(x(1)^2+x(2)^2)-L(2))/L(1)).^2);
    xx = x(1)^2+x(2)^2;
    A = [-x(2)/xx x(1)/xx;...
        x(1)*c/(L(1)*sqrt(xx)) x(2)*c/(L(1)*sqrt(xx))];
    sy = abs(det(A+eps))*A*(sx*inv(A+eps));
else
    error('"opt" input inappropriate: should be either forward or inverse kinematics')
end

% sy = .5.*(sy+inv(sy).*det(sy)); % to make sy symmetric !!!