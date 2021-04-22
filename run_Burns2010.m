% run_Burns2010.m

shoulder = [250 -300]; % displacement between hand and shoulder in x,y (mm)
Hr = 10*pi/180; % Hr: head roll (red)
VHr = (1*pi/180)^2; % VHr: head roll variance (rad^2)
bHr = 1; % bHr: roll compensation gain
Xvis = [0 0] - shoulder;% Xvis: visual IHP (visual coord) (mm)
VXvis = [.75 .75].^2; % VXvis: variance of visual IHP: 2-vector (mm^2)
Tprop = [10 0]*pi/180; % Tprop: proprioceptive IHP (joint angles) (rad)
VTprop = ([.1 .1]*pi/180).^2; % VTprop: variance of proprioceptive IHP (joint angles): 2-vector (rad^2)
Xtar = [0 100] - shoulder; % Xtar: target position (vis coord)
VXtar = [.5 .5].^2; % VXtar: variance of target: 2-vector (mm^2)
VT = 1^2; % VT: constant variance of ref. frame transformation (mm^2)

[Mdir, SMdir, alpVIS, alpPRO] = SoberSabesBayes(Hr, VHr, bHr, Xvis, VXvis, Tprop, VTprop, Xtar, VXtar, VT)

%% Outputs:
% Mdir: initial movement direction (velocity based) (rad)
% SMdir: variance of movement direction (rad^2)
% alpVIS: visual weights (x,y) and covariances
% alpPRO: proprioceptive weights