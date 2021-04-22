function [y, sy] = reftransfo(x, sx, b, H, sH, sT, opt)
%
% bayesian reference frame transformation

if opt == 1, % forward transformation: shoulder --> retinal
    bet = b*H;
elseif opt == -1, % inverse transformation: retinal --> shoulder
    bet = -b*H;
else
    error('"opt" input inappropriate: should be either forward or inverse transformation')
end

[e1,e2] = size(sx);
if e1 ~= e2, sx = diag(sx); end

T = [cos(bet) sin(bet); -sin(bet) cos(bet)];
y = T*x;

SH = abs(H).*sH.*[(sin(bet))^2 (cos(bet))^2; (cos(bet))^2 (sin(bet))^2];
sy = abs(det(T))*T*(sx*inv(T+eps)) + sT*diag([1 1]) + SH;

% sy = .5.*(sy+inv(sy).*det(sy)); % to make sy symmetric !!!