res = [];
res_err = [];

for i = 1:74

    [path_d, pos]=fileparts(names(i)); % Read File
    inpath = string(join([path_d, '\input', pos, '.json'], ''));
    otpath = string(join([path_d, '\', pos, '.dcm'], ''));
    image_path = otpath;

    % Read Input Image
    input_SSR=fileread(inpath);
    input_SSR=jsondecode(input_SSR);
    
    Version = '1.6';
    roi=input_SSR.CropRoi; % RoI
    
    % Image Position 
    position=input_SSR.Type; 
    
    % ReadImage
    image_path = otpath;
    image=dicomread(image_path);
    image=image(:,:,1);
    
    % CMM points based on 9 inch or 12inch
    World=input_SSR.CMM_WorldPoints;
    Dist_pts=input_SSR.CMM_Dist_pts;
     
    % Tracker data
    C2DD=input_SSR.Marker_DD;
    C2R=input_SSR.Marker_Reference;
    
    % meta = dicominfo([path_d,'\meta\0.dcm']);
    
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
    
    rewarp = 0;
    txt_path = string(join([path_d, '\Output\PD\', pos, '_2D.txt'], ''));
    % Open the text file and read its contents
    fileID = fopen(txt_path, 'r');
    data = textscan(fileID, '%f %f %f');  
    fclose(fileID);
    dataArray = [data{1}, data{2}, data{3}];
    xy_o = dataArray(:, 1:3);           

    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=ref_conversion(position,World,Dist_pts,C2R,C2DD,path_d);

    xy= xy_o(:,1:2); % (x, y) coordinates of the fiducials
    order = xy_o(:,3); % fiducial numbers
    xy=[xy(:,1),r-xy(:,2)];
    Ref_CMM_pts_new=Ref_CMM_pts(order,:);
   
    [P0, reprojection_error0, x0, y0]=estimatecameramatrix(xy,Ref_CMM_pts_new);
    P0; % initial projection matrix
    imagePoints_x = xy(:, 1);
    imagePoints_y = xy(:, 2);
    px = x0; % projected x coordinates
    py = y0; % projected y coordinates
    rpe = reprojection_error0; 
    % epsilon = [-200:1:200]; % 870.1248
    % epsilon = [-200:10:200]; % 89.9660
    epsilon = [-200:10:200]; % 37.5924
    p_star = zeros(size(epsilon));
    d = [ones(1,17)]';
    bx = px - rpe;
    by = py - rpe;
    [K, R_ct, Pc1, pp1, pv1] = decomposecamera(P0');
    K_norm = K/K(3,3);
    Rt3=[R_ct -R_ct*Pc1];
    P_norm = K_norm*Rt3;
    Error0 = registration_check(P0',Ref_dist_pts,pos,r,path_d,rewarp); % rpe before optimization

    % % % % % % % % % % Optimization using Camera Projection Matrix
    for j=1:length(epsilon)
        cvx_begin quiet
            variable rx(17); % weights for x coordinates
            minimize ( norm( imagePoints_x.*rx-bx + epsilon(j)*d, 1) )
        cvx_end
        cvx_begin quiet
            variable ry(17); % weights for y coordinates
            minimize ( norm( imagePoints_y.*ry-by + epsilon(j)*d, 1) )
        cvx_end
        % p_star(j)= cvx_optval;
    end
    
    res = horzcat(rx, ry); 
    w = mean(res, 2); % weights for each fiducial
    weights = [transpose(1:17), w]
    
    sortedMatrix = sortrows(weights, 2); % sort the matrix based on the values in second column (ascending order)
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

    %%%%% PLOTS
    % img = imread("")
    % figure;
    % imshow(img); title("After Optimization")
    % hold on;   
    % xy_o = sortrows(xy_o, 3);
    % plot(xy2(:, 1), xy2(:, 2), "o", 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'k')
    % text(xy2(:, 1), xy2(:, 2), num2str(order1), "Color", 'w')
    % % text(xy_orig(:, 1), xy_orig(:, 2), num2str(weights(:, 2)), "Color", 'G', 'FontSize', 8)
    % hold off;

    % 
    % figure;
    % imshow(img); title("Before Optimization")
    % hold on;  
    % xy3 = sortrows(dataArray, 3)
    % % xy_orig = sortrows(xy_orig, 3);
    % plot(xy3(:, 1), xy3(:, 2), "o", 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'k')
    % text(xy3(:, 1), xy3(:, 2), num2str(sort(xy3(:, 3))), "Color", 'w')
    % % text(xy_orig(:, 1), xy_orig(:, 2), num2str(weights(:, 2)), "Color", 'r', 'FontSize', 4)
    % hold off;


