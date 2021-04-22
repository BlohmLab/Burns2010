function [Mdir, SMdir, alpVIS, alpPRO] = SoberSabesBayes(Hr, VHr, bHr, Xvis, VXvis, Tprop, VTprop, Xtar, VXtar, VT)
%
% uses a model for multi-sensory integration
% developed in Sober & Sabes 2003, 2005
% and extends it to be fully Bayesian
%
% Mdir: initial movement direction (velo based)
%
% Hr: head roll
% VHr: head roll variance
% bHr: roll compensation gain
% Xvis: visual IHP (visual coord)
% VXvis: variance of visual IHP: 2-vector
% Tprop: proprioceptive IHP (joint angles)
% VTprop: variance of proprioceptive IHP (joint angles): 2-vector
% Xtar: target position (vis coord)
% VXtar: variance of target: 2-vector
% VT: constant variance of ref. frame transformation

%% body parameters
XS = [250; -300]; % shoulder location
L1 = 300; % upper arm length (mm)
L2 = 450; % lower arm and hand length (mm)

%% proprioceptive IHP in visual coord.
[Xprop, VXprop] = kinematics(Tprop, VTprop, [L1; L2], 1); % forward kinematics
[XpropR, VXpropR] = reftransfo(Xprop+XS, VXprop, bHr, Hr, VHr, VT, 1); % shoulder-to-retinal transformation

%% visual IHP in proprioceptive coordinates
[XvisR, VXvisR] = reftransfo(Xvis+XS, VXvis, bHr, Hr, VHr, VT, -1); % retinal-to-shoulder transformation
[TvisR, VTvisR] = kinematics(XvisR-XS, VXvisR, [L1; L2], -1); % inverse kinematics

%% multi-sensory integration
VEXvis = inv(inv(diag(VXvis)+eps) + inv(VXpropR+eps) + eps); % visual coordinates
EXvis = VEXvis*inv(diag(VXvis+eps))*(Xvis+XS) + VEXvis*inv(VXpropR+eps)*XpropR;
alpVIS = 1/2*(VEXvis*inv(diag(VXvis+eps)) - VEXvis*inv(VXpropR+eps) + eye(2)); % visual weight matrix (visual coord.)

VETprop = inv(inv(diag(VTprop)+eps) + inv(VTvisR+eps)); % proprioceptive coordinates
ETprop = VETprop*inv(diag(VTprop)+eps)*Tprop' + VETprop*inv(VTvisR+eps)*TvisR';
alpPRO = 1/2*(VETprop*inv(diag(VTprop)+eps) - VETprop*inv(VTvisR+eps) + eye(2)); % visual weight matrix (proprio. coord.)

[XtarR, VXtarR] = reftransfo(Xtar+XS, VXtar, 1, Hr, 0, 0, 1); % shoulder-to-retinal transformation
% [XtarR, VXtarR] = reftransfo(Xtar+XS, VXtar, bHr, Hr, VHr, VT, 1); % shoulder-to-retinal transformation
xdotd = XtarR - EXvis; % desired movement vector
Vxdotd = VEXvis + VXtarR;

%% reference frame transformation of movement vector
[xdotdR, VxdotdR] = reftransfo(xdotd, Vxdotd, bHr, Hr, VHr, VT, -1); % retinal-to-shoulder transformation

%% velocity command model --> movement vector
[xdot, sxdot] = velocitycommand(xdotdR, VxdotdR, Tprop, VTprop, ETprop, VETprop, [L1; L2]);

%% movement direction
Mdir = atan2(xdot(2),xdot(1));

% variance orthogonal to movement direction
R = [cos(Mdir) -sin(Mdir);sin(Mdir) cos(Mdir)];
smdir = R*sxdot*inv(R);
SMdir = atan2(smdir(2,2),100);
 
% A0 = [xdot(1)*sqrt(xdot(1)^2+xdot(2)^2) xdot(2)*sqrt(xdot(1)^2+xdot(2)^2); -xdot(2) xdot(1)]./(xdot(1)^2+xdot(2)^2); % Taylor approximation of polar coordinate transformation
% smdir0 = det(A0)*A0*(sxdot*inv(A0+eps));
% smdir = .5.*(smdir0+inv(smdir0).*det(smdir0)); % to make smdir symmetric !!!
% SMdir = smdir(2,2);
