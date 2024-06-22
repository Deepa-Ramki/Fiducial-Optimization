% the fiducial detection is already performed and the coordinates used in
% the following code are obtained from the text file saved during fiducial
% detection (PD Code)

clc
clear all
close all
warning off;
res_order = [];
out_SSR.ErrorMessage =struct('ImagePath', {'ImagePath'}...
    , 'position', {'position'},'ErrorFunction', {'DD/PD/CB'}, 'ErrorMessage', {'Message'},...
    'PlateFiducials', {'Count'}, 'RPE_Error', {'Error'} );
out_SSR.Count = 1;
save out_SSR.mat out_SSR; 

res_err = [];
dataset_path ='C:\Users\Srivibha\Desktop\HTIC_FYP\DataSet - 17Fid - Copy';

files = dir(dataset_path);
dirFlags = [files.isdir];
Folders = files(dirFlags); 
FolderNames = {Folders(3:end).name} ;
resultx = [];
resulty = [];
for j = 1:length(FolderNames)
    dataset_path = "C:\Users\Srivibha\Desktop\HTIC_FYP\DataSet - 17Fid - Copy";
    dataset_path = string(join([dataset_path,'\',char(FolderNames(j))], ''));
    
    files = dir(dataset_path);
    dirFlags = [files.isdir];
    subFolders = files(dirFlags); 
    subFolderNames = {subFolders(3:end).name} ;
    for k = 1:length(subFolderNames)

        txtfiles = dir(fullfile(dataset_path,'\',char(subFolderNames(k)),'\Output\PD\', '*.txt'));
        numTextFiles = numel(txtfiles);
        
        for i = 1:numTextFiles
            txt_path = [dataset_path,'\',char(subFolderNames(k)),'\Output\PD\',txtfiles(i).name];
            fileID = fopen(string(join(txt_path, '')), 'r');
            data = textscan(fileID, '%f %f %f'); 
            fclose(fileID);
            dataArray = [data{1}, data{2}, data{3}];  
            
            if length(dataArray) == 17
                sortedArray = sortrows(dataArray, 3);
                x = sortedArray(:,1);
                y = sortedArray(:,2);

                newFirstRowx = [dataset_path,'\',char(subFolderNames(k)),'\', txtfiles(i).name(1:2)];
                newFirstRowy = [dataset_path,'\',char(subFolderNames(k)),'\', txtfiles(i).name(1:2)];

                originalCellArrayx = num2cell(x);
                originalCellArrayy = num2cell(y);

                insertedCellArrayx = vertcat({newFirstRowx}, originalCellArrayx);
                insertedCellArrayy = vertcat({newFirstRowy}, originalCellArrayy);
                
                resultx = horzcat(resultx, insertedCellArrayx);
                resulty = horzcat(resulty, insertedCellArrayy);

            else
                continue
            end
        end
    end
    x_coordinates = cell2mat(resultx(2:end,:));
    x_avg = [];
    for m = 1:17
        avg = mean(x_coordinates(m,:));
        x_avg(m) = avg;
    end
    x_avg = x_avg';
    
    y_coordinates = cell2mat(resulty(2:end,:));
    y_avg = [];
    for m = 1:17
        avg = mean(y_coordinates(m,:));
        y_avg(m) = avg;
    end
    y_avg = y_avg';
end   
names = resultx(1,:);

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
    % arr = horzcat(arr, sortrows(dataArray, 3), res_weights(:, i), round(res_weights(:, i)))
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
    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=Ref_conversion1(position,World,Dist_pts,C2R,C2DD,path_d);
    xyBottom = xy_bottom(:,1:2);
    order_bottom = xy_bottom(:,3);
    xyBottom=[xyBottom(:,1),r-xyBottom(:,2)];
    Ref_CMM_pts_bottom = Ref_CMM_pts(order_bottom,:);
    [P0,reprojection_error0, x0, y0]=estimateCameraMatrix1(xyBottom,Ref_CMM_pts_bottom);
    P0;     [K,R_ct,Pc1, pp1, pv1] = decomposecamera(P0');    K_norm = K/K(3,3);    Rt3=[R_ct -R_ct*Pc1];    P_norm = K_norm*Rt3; 
    Error_bottom = Registration_check1(P0',Ref_dist_pts,pos,r,path_d,rewarp);

            %%%%% RPE when top 10 fiducials + Significant fiducials are chosen for calibration (Non significant fiducials are excluded)     
    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=Ref_conversion1(position,World,Dist_pts,C2R,C2DD,path_d);
    xyTop= xy_top(:,1:2);
    order_top = xy_top(:,3);
    xyTop=[xyTop(:,1),r-xyTop(:,2)];
    Ref_CMM_pts_top = Ref_CMM_pts(order_top,:);   
    [P1,reprojection_error1, x1, y1]=estimateCameraMatrix1(xyTop,Ref_CMM_pts_top);
    P1;     [K,R_ct,Pc1, pp1, pv1] = decomposecamera(P1');     K_norm = K/K(3,3);     Rt3=[R_ct -R_ct*Pc1];     P_norm = K_norm*Rt3;
    Error_top = Registration_check1(P1',Ref_dist_pts,pos,r,path_d,rewarp);

            %%%%% RPE when all the 17 fiducials are chosen for calibration    
    [Ref_CMM_pts, Ref_dist_pts, ref2cam]=Ref_conversion1(position,World,Dist_pts,C2R,C2DD,path_d);
    xyOriginal= xy_orig(:,1:2);
    order_orig = xy_orig(:,3);
    xyOriginal = [xyOriginal(:,1),r-xyOriginal(:,2)];
    Ref_CMM_pts_orig = Ref_CMM_pts(order_orig,:);  
    [P2,reprojection_error2, x2, y2]=estimateCameraMatrix1(xyOriginal,Ref_CMM_pts_orig);
    P2;     [K,R_ct,Pc1, pp1, pv1] = decomposecamera(P1');     K_norm = K/K(3,3);     Rt3=[R_ct -R_ct*Pc1];     P_norm = K_norm*Rt3;
    Error_orig = Registration_check1(P2',Ref_dist_pts,pos,r,path_d,rewarp);

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

    Error_min = min([Error_orig, Error_bottom, Error_top]); % comparing the RPEs
    err = horzcat(Error_orig, Error_bottom, Error_top, Error_min)
    res_err = vertcat(res_err, err);    
end


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
    % ref_dist_pts = ref_dist_pts(1:3,:)';
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
