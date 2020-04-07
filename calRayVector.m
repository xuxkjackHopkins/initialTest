function [pixel_opt, rayVector] = calRayVector(ydata, geoCaliPara, mirror_angle, pixel_size)
% This function is used to calculate the vector of projection ray from
% center of selected CCD pixel in pin-hole camera model.
% Input:
%   ydata (1x2) is pixel index in 2D optical image corresponding to P in CBCT
%       The direction of ydata(1): v, and ydata(2): u are consistent with
%       the 2D image direction, i.e., top to down;  left to right
%   geoCaliPara (1x12) is an array of the geometry calibration parameters, the elements are 
%       geoCaliPara(1) focal-to-detector length Lfd, 
%       geocaliPara(2) focal-to-object distance Lfo;
%   3D CBCT coord to 3D optical coord calibration: 3 rotations and 3 shifts
%       geoCaliPara(3)-(5) rotation around x-, y- and z-axis (degree); 
%       geoCaliPara(6) - (8): shift in x-, y- and z-axis
%       geoCaliPara(9)-(10) center of imaging plane (u0,v0); 
%       geoCaliPara(11)-(12) rotation center of raw images (RotCen_u, RotCen_v);
%   mirror_angle : the angle of the three-mirror system.
%   pixel_size :  the physical size of CCD pixel
%       physical size of CCD pixel, unit: mm/pixel. 13.0e-3 for MuriGlo and 13.5e-3 for Dual-use 

% Output:
%   pixel_opt : pixel position in the 3D optical coordinate.
%   rayVector: ray of the line connecting center of each Pixel and focal
%   point in 3D optical coordinate.
%  For detailed description of 3D optical coordinate, please refer:
%    Zhang, B., Iordachita, I., Wong, J. W., and Wang, K. K.-H., “Multi-projection bioluminescence tomography guided system for small animal radiation research platform ” Proc. SPIE 9701, Multimodal Biomedical Imaging XI, 97010J (2016).
%
%  Jack Xu, Johns Hopkins University, Dec 2019
Lfd = geoCaliPara(1);
Lfo = geoCaliPara(2);
% Lfo = x(1)*9;

theta_x = geoCaliPara(3)*3.1416/180; % rotation around x-axis
theta_y = geoCaliPara(4)*3.1416/180; % rotation around y-axis
theta_z = geoCaliPara(5)*3.1416/180; % rotation around z-axis

xot = geoCaliPara(6); % x-shift 
yot = geoCaliPara(7); % y-shift
zot = geoCaliPara(8);  % z-shift

u0 = geoCaliPara(9); % image position of origin of 3D optical coordinate
v0 = geoCaliPara(10);

% image rotation center;
% image rotation caused by the 3 mirrors
RotCen_u = geoCaliPara(11);
RotCen_v = geoCaliPara(12);

% mirror rotation angle, degree
theta0 = mirror_angle*3.1416/180; 
%% rotate 2D optical image
% in-plane image rotation caused by the rotation of the 3-mirror system
img_rot_deg = theta0; % mirror angle
img_Rot = [cos(img_rot_deg), sin(img_rot_deg); -sin(img_rot_deg), cos(img_rot_deg)];    
ydata_rot = img_Rot*[ydata(2) - RotCen_u; ydata(1) - RotCen_v] + [RotCen_u; RotCen_v]; % rotate image at true rotation center. (row, col) in matlab is (v,u) in 2D optical image plane.
%% ray connecting focal point and pixel position in 3D optical coordinate
% pixel position in 3D optical coordinate
opt_x = -(ydata_rot(1) - u0)*pixel_size;
opt_y = (ydata_rot(2) - v0)*pixel_size;
opt_z = Lfd + Lfo;
pinhole = [ 0 0 Lfo];
pixel_opt = [opt_x opt_y opt_z];
rayVector = pixel_opt - pinhole;
