function [y, sy] = velocitycommand(x, sx, t, st, te, ste, L)
%
% computes the velocity command from current hand position (t) and
% an estimation of IHP (te)

a = L(1)*sin(t(2)) + L(2);
b = L(1)*cos(t(2))+eps;
JT = [-cos(t(1))*a (-b*sin(t(1))); ...
    -sin(t(1))*a  b*cos(t(1))];
ae = L(1)*sin(te(2)) + L(2);
be = L(1)*cos(te(2))+eps;
JTinv = [-cos(te(1))/ae  -sin(te(1))/ae; -sin(te(1))/be  cos(te(1))/be];
y = JT*(JTinv*x);

A11 = [(sin(te(1))/a)^2 (L(1)*sin(te(1))*cos(te(1))*cos(te(2))/a^3); (L(1)*sin(te(1))*cos(te(1))*cos(te(2))/a^3)  (L(1)*cos(te(1))*cos(te(2))/a^2)^2];
A12 = [(cos(te(1))/a)^2 (L(1)*sin(te(1))*cos(te(1))*cos(te(2))/a^3); (L(1)*sin(te(1))*cos(te(1))*cos(te(2))/a^3)  (L(1)*sin(te(1))*cos(te(2))/a^2)^2];
A21 = [(cos(te(1))/b)^2 (L(1)*sin(te(1))*cos(te(1))*sin(te(2))/b^3); (L(1)*sin(te(1))*cos(te(1))*sin(te(2))/b^3)  (L(1)*sin(te(1))*sin(te(2))/b^2)^2];
A22 = [(sin(te(1))/b)^2 (L(1)*sin(te(1))*cos(te(1))*sin(te(2))/b^3); (L(1)*sin(te(1))*cos(te(1))*sin(te(2))/b^3)  (L(1)*cos(te(1))*sin(te(2))/b^2)^2];
s11 = trace(A11*ste);
s12 = trace(A12*ste);
s21 = trace(A21*ste);
s22 = trace(A22*ste);
VJTinv = [s11 s12; s21 s22];
st = det(JTinv)*JTinv*(sx*inv(JTinv+eps));% + VJTinv;
sy = det(JT)*JT*(st*inv(JT+eps));

% sy = .5.*(sy+inv(sy).*det(sy)); % to make sy symmetric !!!