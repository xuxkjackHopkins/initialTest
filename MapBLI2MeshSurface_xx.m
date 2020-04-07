function [Det_coor, Det_val] = MapBLI2MeshSurface_xx(img, pixel_size, mesh, time_factor, geometry_calibration_param, detector_threshold)
% This function is used to map 2D BLI data to 3D mesh surface with our previously published
% geometrical calibration method (Proc. SPIE 9701, Multimodal Biomedical Imaging XI, 97010J (2016)). 
% The center of each selected pixel was mapped to 3D mesh surface. The
% mapped position was used as the detector coordinate.
%   input:
%       img: filtered BL images, include data: img.Images, binning size:img.binnings_size, projection angle: img.projection_angle 
%       pixel_size: physical size of CCD pixel.
%       mesh: 3D mesh generated from CBCT image.
%       time_factor: time factor used to calibrate in vivo signal change.
%       geometry_calibration_param: parameters for geometrical calibration between 2D
%       BLI coordinate and 3D mesh coordinate.
%       detector_threshold.ratio : a lower bound threshold used to choose detector
%       values.
%       detector_threshold.angle: An upper bound angle relative to camera imaging axis.
%       detector_threshold.pixel_value: A lower bound value to choose pixels after background substraction for the image.
%   output:
%       Det_coor: Coordinates of detectors used for BLT reconstruction.
%       Det_val: Values of the detectors.
%
% Wirtten by Jack Xu, Johns Hopkins University, March 2020

Nfreq = length(img.Images); 
mirror_angle = img.projection_angle;         % imaging projection is equal to mirror rotating angle.
pixel_value_threshold = detector_threshold.pixel_value; % threshold used to choose pixels used for reconstruction
%pixel_index = [];
%pixel_number = zeros(Nfreq,1);
% Pick pixel based on pixel value threshold
for ii = 1:Nfreq
    image = img.Images{ii};
    pixel_value_threshold = pixel_value_threshold/img.exposure_time(ii);
    pixel_select_index = pickPixel(image, pixel_value_threshold);
    if ii > 1
        pixel_select_index = intersect(pixel_index_binning, pixel_select_index,'rows');
    end
    pixel_index_binning = pixel_select_index;
end
% Calculate the pixel coordinate at 1x1 binning pixel coordinate
BLI = imresize(img.Images{1}, img.binning_size(1), 'bilinear'); % Resize image to 1x1 binning, because geometrical calibration parameters were calculated based on 1x1 binning images.
if rem(log2(size(BLI,1)),1) ~= 0                           % if image size is not regular, padding 0 values to the boundary.
    image_size = size(BLI,1);
    log_num = round(log2(image_size));
    image_fullsize = 2^log_num;
    %         pad_size = (image_fullsize - image_size)/2;
    %         BLI{ii} = padarray(BLI{ii},[pad_size, pad_size],0,'both');
    pad_size = image_fullsize - image_size;
%    BLI{ii} = padarray(BLI{ii},[pad_size, pad_size],0,'pre'); % pad image before each array
    pixel_index_1x1_binning = pixel_index_binning*img.binning_size(1) - img.binning_size(1)/2 + pad_size;
else
    pixel_index_1x1_binning = pixel_index_binning*img.binning_size(1) - img.binning_size(1)/2; % use the center of each binned pixel as the pixel coordinate
end
pixel_opt = zeros(size(pixel_index_1x1_binning,1),3);
rayVector = zeros(size(pixel_index_1x1_binning,1),3);
parfor jj = 1 : size(pixel_index_1x1_binning,1)
    [pixel_opt(jj, :), rayVector(jj, :)] =calRayVector (pixel_index_1x1_binning(jj, :), geometry_calibration_param, mirror_angle, pixel_size); %pixel position in 3D optical coordinate
end
%% Select visible 3D mesh surface
surface_norm_visiable_threshold_deg = detector_threshold.angle;  % An upper bound angle between surface normal and camera imaging axis. 
rotation_angle = - mirror_angle;           % The mirror rotation angle is opposite to the CBCT rotation angle
% Calculate the surface nodes coordinate in 3D optical coordinate.
node_coord_in_3Doptical = Shift_and_rotate_CBCT_to_3D_optical(geometry_calibration_param, rotation_angle, mesh.nodes); 
% Calculate the center and area of each surface facet.
[surface_center_coord_CBCT, surface_area_CBCT] = ...
           mesh_face_center_area(mesh.face, mesh.nodes);
% Calculate surface normal of the 3D optical coordinates
surface_norm_3Doptical = - trisurfnorm(node_coord_in_3Doptical, mesh.face); % neagtive (-) means to convert all the surface normal toward to outside.
% Delete a range (1mm, empirical value) in y-axis(the body axis)to select
% surface facet.
y_min = min(surface_center_coord_CBCT(:, 2)) + 1;
y_max = max(surface_center_coord_CBCT(:, 2)) - 1;
surface_select_ind = surface_center_coord_CBCT(:, 2) > y_min ...
                    & surface_center_coord_CBCT(:, 2) < y_max;
% Find the upper surface facets which are visiable to the CCD camera
surface_norm_3Doptical_deg = acos(surface_norm_3Doptical(:,3))/pi*180;      % angle between surface normal and CCD camera.
surface_norm_visiable_index = surface_norm_3Doptical_deg < surface_norm_visiable_threshold_deg;       % angle less than the bound angle.      
surface_norm_visiable_selected_index = surface_norm_visiable_index & surface_select_ind;
surface_visible_triangle = mesh.face(surface_norm_visiable_selected_index, :);
%surface_visible_node_coord_in_3Doptical = node_coord_in_3Doptical(surface_norm_visiable_selected_index,:);
%% Find the intersection between the ray connection focal point and pixel and 3D mesh surface (start here)
count = zeros(size(pixel_index_1x1_binning,1),1);
triangle_its = zeros(size(pixel_index_1x1_binning,1),1); % index of trianlge where the intersection is in
u_its = zeros(size(pixel_index_1x1_binning,1),1); % baricenteric coordinate u
v_its = zeros(size(pixel_index_1x1_binning,1),1);% baricenteric coordinate v
t_its = zeros(size(pixel_index_1x1_binning,1),1);% distance from the origin to the intersection, see rayTriangleIntersection function below.
parfor m = 1:size(pixel_index_1x1_binning,1) % loop selected pixel number
    %origin = [0 0 geometry_calibration_param(2)]; %pinhole position
    origin = pixel_opt(m,:);
    ray = rayVector(m,:);
    direction = ray/norm(ray);
    for n = 1:size(surface_visible_triangle,1)
        v0 = node_coord_in_3Doptical(surface_visible_triangle(n,1),:);  % triangle vertice 1
        v1 = node_coord_in_3Doptical(surface_visible_triangle(n,2),:);  % triangle vertice 2
        v2 = node_coord_in_3Doptical(surface_visible_triangle(n,3),:);  % triangle vertice 3
        [flag, u, v, t] = rayTriangleIntersection(origin, direction, v0, v1, v2);   % Ray/Triangle intersection, flag 0 reject, 1 intersect
        if flag == 1
            if count(m) < 1
                count(m) = count(m) + 1;
                triangle_its(m) = n; 
                u_its(m) = u;
                v_its(m) = v;
                t_its(m) = t;
            elseif abs(t) < abs(t_its(m)) % choose the intersection that is close to the detector
                triangle_its(m) = n;
                u_its(m) = u;
                v_its(m) = v;
                t_its(m) = t;
            end                     
        end
    end
    
end

%% Delete selected pixel that has no intersection with the 3D mesh surface
index = count <1;
pixel_index_binning(index,:) = []; % delete selected pixels that are not mapped on the mesh surface.
triangle_its(index,:) = [];
u_its(index,:) = [];
v_its(index,:) = [];
det_coordinate = zeros(size(pixel_index_binning,1),3);
parfor det_num = 1 : size(pixel_index_binning,1)
    intersect_triangle = surface_visible_triangle(triangle_its(det_num),:)
    triangle_v1 = mesh.nodes(intersect_triangle(1),:);
    triangle_v2 = mesh.nodes(intersect_triangle(2),:);
    triangle_v3 = mesh.nodes(intersect_triangle(3),:);
    det_coordinate(det_num,:) = u_its(det_num)*triangle_v2 + v_its(det_num)*triangle_v3 + (1 - u_its(det_num) - v_its(det_num))*triangle_v1;
end
% %% Calculate the 2D optical coordinates cooresponding to the 3D CBCT surface facet center  
% xdata = zeros(size(surface_center_coord_ROI_visiable, 1), 4);
% parfor j = 1:size(surface_center_coord_ROI_visiable, 1)
%     xdata(j, :) = [rotation_angle, surface_center_coord_ROI_visiable(j, :)];    % xdata: [angle, surface center coordinate of CBCT]
% end
% ydata = zeros(size(xdata, 1), 2);
% parfor i = 1:size(xdata, 1)
%     ydata(i, :) = Cal_fwdCalc_RotCorr(geometry_calibration_param, xdata(i, :), pixel_size);     % ydata is the 2D optical coordinate of the visible surface facet center.
%     % ydata1 is column number and ydata2 is row number of 2D BL images.
% end
%% obtain the detector value.
Det_val_cell = cell(1, Nfreq);  % Detector value at each wavelength
ydata = pixel_index_binning;
parfor j = 1 : Nfreq
   for k = 1: size(ydata, 1)
      Det_val_cell{j}(k,1) =  img.Images{j}(round(ydata(k, 1)), round(ydata(k, 2)));
   end
end
Det_coor = det_coordinate; % Detector coordinate
Det_val = cell2mat(Det_val_cell);             % Detector value 
%% delete the hot noise spot.
% % % % % % eliminate hot spot noise
% % % % % % the hot spot noise is randomly appeared, it is unlikely to appear at the same position on all the BLIs
% % % % % % the strategy is to calculate the coefficient of variation(CV) for each detector, if the CV > threshold, then delete that point 
c_vairation_threshold = 1.0; % after previouse test, all of the normal detectors are <0.6, but abnormal points are > 1, so threshold is set to 1.
Det_val_mean = mean(Det_val, 2);
Det_val_std = std(Det_val, 1, 2);
c_vairation = Det_val_std./Det_val_mean;
c_vairation_ind = c_vairation>c_vairation_threshold;
Det_val(c_vairation_ind, :) = [];
Det_coor(c_vairation_ind, :) = [];
if isempty(Det_val)
    errordlg('No detectors' , 'Mapping error')
    return
end
%% Time factor correction: scale the detector value with time factor.
parfor n = 1:Nfreq
  Det_val(:,n)=Det_val(:,n)/time_factor(n);
end
end