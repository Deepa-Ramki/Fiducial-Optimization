% the fiducial detection is already performed and the coordinates used in
% the following code are obtained from the text file saved during fiducial
% detection (PD Code)

load x.mat     % x coordinates
load x_avg.mat % average of the x coordinates of the fiducials
load y.mat     % y coordinates
load y_avg.mat % average of the y coordinates of the fiducials
load names.mat % names of the files/images

% % % % % % Optimization using Fiducial Coordinates

res_weights = [];
for m = 1:74 % for all images
    A_x = x_coordinates(:, m);
    B_x = x_avg;
    cvx_begin
            variable x(17) 
            minimize(norm(A_x .* x - B_x)) % weights for x coordinates for an image
    cvx_end

    A_y = y_coordinates(:, m);
    B_y = y_avg;
    cvx_begin
            variable y(17)
            minimize(norm(A_y .* y - B_y)) % weights for y coordinates for an image
    cvx_end

    w = horzcat(x, y);
    weights = mean(w, 2); % weights of the fiducials
    res_weights = horzcat(res_weights, weights);
end

arr = [];

for i = 1:74
    [path_d, pos]=fileparts(string(join(names{1, i},'')));
    inpath = string(join([path_d, '\input', pos, '.json'], ''));
    otpath = string(join([path_d, '\', pos, '.dcm'], ''));
    % image_path = otpath;  
    input_SSR=fileread(inpath); % Read Input Image
    input_SSR=jsondecode(input_SSR);
    Version = '1.6';
    
    roi=input_SSR.CropRoi; % RoI
    position=input_SSR.Type; % Image Position 
    image_path = otpath; % ReadImage
    image=dicomread(image_path);
    image=image(:,:,1);
    
    World=input_SSR.CMM_WorldPoints; % CMM points based on 9 inch or 12inch
    Dist_pts=input_SSR.CMM_Dist_pts;

    C2DD=input_SSR.Marker_DD; % Tracker data
    C2R=input_SSR.Marker_Reference;
    
    %Resolution of the image
    [r,c,ch]=size(image);
    
    if max(image(:))<256
        image=uint8(image);
    else
        image=im2uint16(image);
    end
    
    % original = image;
    [imageT,xmin,ymin,xmax,ymax]=Crop_Image(image);
    [row,col] = size(imageT);
    aspRatio1 = row/col;
    if aspRatio1 > 0.9 && aspRatio1 < 1.1  
        image = imageT;    
    end
    output_SSR.CropRoi=[xmin;ymin;xmax;ymax];
    [r,c,ch]=size(image);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calibration
    rewarp = 0;
    txt_path = string(join([path_d, '\Output\PD\', pos, '_2D.txt'], ''));
    % Open the text file and read its contents
    fileID = fopen(txt_path, 'r');
    data = textscan(fileID, '%f %f %f');  
    fclose(fileID);
    
    dataArray = [data{1}, data{2}, data{3}];
    xy_bottom = dataArray(:, 1:3);
    xy_top = dataArray(:, 1:3);
    xy_orig = dataArray(:, 1:3);
    res_weights(:, i)
    weights = round(res_weights(:, i));
    weights = [transpose(1:17), weights];    
    
    sortedMatrix = sortrows(weights, 2); 
    bottom10Values = sortedMatrix(1:10, :);
    bottom10Values = sortrows(bottom10Values, 1);  
    valuesNotToRemove_bottom = bottom10Values(:, 1);
    top10Values = sortedMatrix(8:end, :);
    top10Values = sortrows(top10Values, 1);
    valuesNotToRemove_top = top10Values(:, 1);

    valuesNotToRemove_bottom = [valuesNotToRemove_bottom; 4; 7]; % including the 4th and 7th fiducials
    valuesNotToRemove_top = [valuesNotToRemove_top; 4; 7]; % including the 4th and 7th fiducials
    idx = ~ismember(xy_bottom(:, 3), valuesNotToRemove_bottom);
    xy_bottom(idx > 0, :) = [];
    idx1 = ~ismember(xy_top(:, 3), valuesNotToRemove_top);
    xy_top(idx1 > 0, :) = [];
    idx = ismember(xy_bottom(:, 3), [6, 9]);  % excluding the 6th and 9th fiducials
    xy_bottom(idx > 0, :) = [];
    idx1 = ismember(xy_top(:, 3), [6, 9]); % excluding the 6th and 9th fiducials
    xy_top(idx1 > 0, :) = [];

    %%%%% RPE when bottom 10 fiducials + Significant fiducials are chosen for calibration (Non significant fiducials are excluded)
    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=ref_conversion(position,World,Dist_pts,C2R,C2DD,path_d);
    xyBottom = xy_bottom(:,1:2);
    order_bottom = xy_bottom(:,3);
    xyBottom=[xyBottom(:,1),r-xyBottom(:,2)];
    Ref_CMM_pts_bottom = Ref_CMM_pts(order_bottom,:);
    [P0,reprojection_error0, x0, y0]=estimatecameramatrix(xyBottom,Ref_CMM_pts_bottom);
    P0;     [K,R_ct,Pc1, pp1, pv1] = decomposecamera(P0');    K_norm = K/K(3,3);    Rt3=[R_ct -R_ct*Pc1];    P_norm = K_norm*Rt3; 
    Error_bottom = registration_check(P0',Ref_dist_pts,pos,r,path_d,rewarp);

    %%%%% RPE when top 10 fiducials + Significant fiducials are chosen for calibration (Non significant fiducials are excluded)     
    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=ref_conversion(position,World,Dist_pts,C2R,C2DD,path_d);
    xyTop= xy_top(:,1:2);
    order_top = xy_top(:,3);
    xyTop=[xyTop(:,1),r-xyTop(:,2)];
    Ref_CMM_pts_top = Ref_CMM_pts(order_top,:);   
    [P1,reprojection_error1, x1, y1]=estimatecameramatrix(xyTop,Ref_CMM_pts_top);
    P1;     [K,R_ct,Pc1, pp1, pv1] = decomposecamera(P1');     K_norm = K/K(3,3);     Rt3=[R_ct -R_ct*Pc1];     P_norm = K_norm*Rt3;
    Error_top = registration_check(P1',Ref_dist_pts,pos,r,path_d,rewarp);

    %%%%% RPE when all the 17 fiducials are chosen for calibration    
    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=ref_conversion(position,World,Dist_pts,C2R,C2DD,path_d);
    xyOriginal= xy_orig(:,1:2);
    order_orig = xy_orig(:,3);
    xyOriginal = [xyOriginal(:,1),r-xyOriginal(:,2)];
    Ref_CMM_pts_orig = Ref_CMM_pts(order_orig,:);  
    [P2,reprojection_error2, x2, y2]=estimatecameramatrix(xyOriginal,Ref_CMM_pts_orig);
    P2;     [K,R_ct,Pc1, pp1, pv1] = decomposecamera(P1');     K_norm = K/K(3,3);     Rt3=[R_ct -R_ct*Pc1];     P_norm = K_norm*Rt3;
    Error_orig = registration_check(P2',Ref_dist_pts,pos,r,path_d,rewarp);

    Error_min = min([Error_orig, Error_bottom, Error_top]); % comparing the RPEs
    err = horzcat(Error_orig, Error_bottom, Error_top, Error_min)
    res_err = vertcat(res_err, err);    
end

    % % % PLOTS
    % if Error_bottom < Error_top
    %     figure;
    %     subplot(122)
    %     imshow(image); title("After Optimization")
    %     hold on;   
    %     xy_orig = sortrows(xy_orig, 3);
    %     plot(xy_bottom(:, 1), xy_bottom(:, 2), "o", 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'k')
    %     text(xy_bottom(:, 1), xy_bottom(:, 2), num2str(order_bottom), "Color", 'w')
    %     % text(xy_orig(:, 1), xy_orig(:, 2), num2str(weights(:, 2)), "Color", 'r', 'FontSize', 4)
    %     hold off;
    % else
    %     figure;
    %     imshow(image); title("Top")
    %     hold on;
    %     xy_orig = sortrows(xy_orig, 3);
    %     plot(xy_top(:, 1), xy_top(:, 2), "o", 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'k')
    %     text(xy_top(:, 1), xy_top(:, 2), num2str(order_top), "Color", 'w')
    %     % text(xy_orig(:, 1), xy_orig(:, 2), num2str(weights(:, 2)), "Color", 'r', 'FontSize', 4)
    %     hold off;
    % end
    % subplot(121);
    % imshow(image); title("Before Optimization")
    % hold on;   
    % % xy_orig = sortrows(xy_orig, 3);
    % plot(xy_orig(:, 1), xy_orig(:, 2), "o", 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'k')
    % text(xy_orig(:, 1), xy_orig(:, 2), num2str(sort(order_orig)), "Color", 'w')
    % % text(xy_orig(:, 1), xy_orig(:, 2), num2str(weights(:, 2)), "Color", 'r', 'FontSize', 4)
    % hold off;