res = [];
res_err = [];
names = raw(1, 1:end);
for i = 1:74
    [path_d, pos]=fileparts(names(i));
    inpath = string(join([path_d, '\input', pos, '.json'], ''));
    otpath = string(join([path_d, '\', pos, '.dcm'], ''));
    image_path = otpath;
    % Read Input Image
    input_SSR=fileread(inpath);
    input_SSR=jsondecode(input_SSR);
    
    Version = '1.6';
    roi=input_SSR.CropRoi; %Roi
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calibration
    rewarp = 0;
    txt_path = string(join([path_d, '\Output\PD\', pos, '_2D.txt'], ''));
    % Open the text file and read its contents
    fileID = fopen(txt_path, 'r');
    data = textscan(fileID, '%f %f %f');  
    fclose(fileID);
    dataArray = [data{1}, data{2}, data{3}];
    xy_o = dataArray(:, 1:3);           

    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=Ref_conversion1(position,World,Dist_pts,C2R,C2DD,path_d);

    xy= xy_o(:,1:2); % (x, y) coordinates of the fiducials
    order = xy_o(:,3); % fiducial numbers
    xy=[xy(:,1),r-xy(:,2)];
    Ref_CMM_pts_new=Ref_CMM_pts(order,:);
    
    [P0, reprojection_error0, x0, y0]=estimateCameraMatrix1(xy,Ref_CMM_pts_new);
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
    Error0 = Registration_check1(P0',Ref_dist_pts,pos,r,path_d,rewarp); % rpe before optimization

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
    
    bottom10Values = sortedMatrix(1:10, :); % extract the bottom 10 values along with their corresponding serial numbers
    bottom10Values = sortrows(bottom10Values, 1);
    valuesNotToRemove = bottom10Values(:, 1);
    valuesNotToRemove = [valuesNotToRemove; 7; 4]; 
    idx = ~ismember(xy_o(:, 3), valuesNotToRemove);
    xy_o(idx > 0, :) = []; 
    
    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=Ref_conversion1(position,World,Dist_pts,C2R,C2DD,path_d);
    xy1 = xy_o(:,1:2);
    xy2 = xy_o(:,1:2);
    order1 = xy_o(:,3);
    xy1 = [xy1(:,1),r-xy1(:,2)]; 
    Ref_CMM_pts_new = Ref_CMM_pts(order1,:);
    [P1,reprojection_error1, x1, y1]=estimateCameraMatrix1(xy1,Ref_CMM_pts_new);
    P1; % projection matrix after optimization
    [K,R_ct,Pc1, pp1, pv1] = decomposecamera(P1');
    K_norm = K/K(3,3);
    Rt3=[R_ct -R_ct*Pc1];
    P_norm = K_norm*Rt3;
    Error1 = Registration_check1(P1',Ref_dist_pts,pos,r,path_d,rewarp);
    FinalError = min(Error0, Error1);
    err = horzcat(Error0, Error1, FinalError);
    res_err = vertcat(res_err, err);
end

        % similar to the method in which weights are obtained using
        % fiducial coordinates, the RPEs in various cases can be compared

    %%%%% PLOTS
    % img = imread("C:\Users\Srivibha\AppData\Local\Temp\293bb57f-a2b4-4528-b851-cc8de529edcf_all_images.zip.dcf\all_images\SurgeDataSet 61LP.png")
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


function ER=Registration_check1(P,W_dist_pts,pos,r,path_d,rewarp)
   if rewarp == 1
    folder = 'Rewarp\';
   else
       folder = '';
   end
    blob=imread(string(join([path_d,'\Output\PD\', pos, 'bw.png'], '')));
    stats=regionprops(blob,'Centroid');
    Centre=cat(1,stats.Centroid);
    
    projected_2d_pts=P*W_dist_pts;
    projected_2d_pts=projected_2d_pts./projected_2d_pts(3,:);
    
    projected_2d_pts= projected_2d_pts(1:2,:)';
    projected_2d_pts= [projected_2d_pts(:,1), r-projected_2d_pts(:,2)];
    
    distance_compute=pdist2(projected_2d_pts,Centre);
    
    actualSubscripts = find(distance_compute(:)<3);
    [row, col, z] = ind2sub(size(distance_compute), actualSubscripts);
    error_in_pixels = [];
    
    for i = 1:size(distance_compute,1)
        err = min(distance_compute(i,:));
        error_in_pixels = [error_in_pixels,err];
    end
    
    out = [row, col];
    val=distance_compute(row,col);
    ER=mean(diag(val));

end


function [Cmm_pts,ref_dist_pts,ref2cam]=Ref_conversion1(position,W,distpts,C2R,C2D,path_d)


 %Reference marker
    tx = C2R.tx;
    ty = C2R.ty;
    tz = C2R.tz;
    trans =[tx ty tz];
    q=C2R.Rotation;   
    rot=quat2rotm([-q(4) q(1) q(2) q(3)]);     
    Rt= [rot trans'; 0 0 0 1]; 
    cam2ref =Rt;
    
 %Detector Marker   
    tx2 = C2D.tx;
    ty2 = C2D.ty;
    tz2 = C2D.tz;   
    q=C2D.Rotation;   
    rot2=quat2rotm([-q(4) q(1) q(2) q(3)]);
    trans2 =[tx2 ty2 tz2];  
    Rt2 = [rot2 trans2'; 0 0 0 1]; 
    cam2DD =Rt2;

    ref2cam = inv(cam2ref);
    ref2DD = ref2cam*cam2DD;
    
    ref_DD= ref2DD*[W';ones(1,size(W,1))];
    Cmm_pts =ref_DD(1:3,:)';


    ref_dist_pts=ref2DD*[distpts';ones(1,size(distpts,1))];
%     ref_dist_pts =ref_dist_pts(1:3,:)';
end


function [camMatrix, reprojectionErrors, x, y] = estimateCameraMatrix1(imagePoints, worldPoints)

    X = worldPoints(:,1);
    Y = worldPoints(:,2);
    Z = worldPoints(:,3);
    
    M = numel(X);
    vec_1 = ones(M,1,'like',X);
    vec_0 = zeros(M,1,'like',X);
    
    % Normalize image points
    % [pts, ~, Tinv] = normalizePoints(imagePoints', 2, class(imagePoints));
    
    imagePoints=imagePoints';
    imagePoints=([imagePoints;vec_1']);
    [pts,T]=normalise2dpts(imagePoints);
    
    Tinv=inv(T);
    u = pts(1, :)';
    v = pts(2, :)';  
    A = [X      Y      Z      vec_1  vec_0  vec_0  vec_0  vec_0  -u.*X  -u.*Y  -u.*Z -u;
         vec_0  vec_0  vec_0  vec_0  X      Y      Z      vec_1  -v.*X  -v.*Y  -v.*Z -v];
    
    [V, D] = eig(A'*A, 'vector');
    [~, idx] = min(D);
    P = V(:, idx(1));
    
    camMatrix = reshape(P,4,3);
    
    % Set back the normalization transform
    camMatrix = camMatrix * Tinv';
    
    if camMatrix(end) < 0
        camMatrix = -camMatrix;
    end

    if nargout > 1
        p = [worldPoints, vec_1] * camMatrix;
        x = p(:, 1)./p(:, 3);
        y = p(:, 2)./p(:, 3);
        reprojectionErrors = sqrt((imagePoints(1,:)'-x).^2+(imagePoints(2,:)'-y).^2);
    end
end


