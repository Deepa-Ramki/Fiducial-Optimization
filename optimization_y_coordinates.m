clc
clear all
close all
warning off;

out_SSR.ErrorMessage =struct('ImagePath', {'ImagePath'}...
    , 'position', {'position'},'ErrorFunction', {'DD/PD/CB'}, 'ErrorMessage', {'Message'},...
    'PlateFiducials', {'Count'}, 'RPE_Error', {'Error'} );
out_SSR.Count = 1;
save out_SSR.mat out_SSR; 


dataset_path ='C:\Users\Srivibha\Desktop\HTIC_FYP\demo';

files = dir(dataset_path);
dirFlags = [files.isdir];
Folders = files(dirFlags); 
FolderNames = {Folders(3:end).name} ;
result = [];
for j = 1:length(FolderNames)
    dataset_path = 'C:\Users\Srivibha\Desktop\HTIC_FYP\demo';
    dataset_path = [dataset_path,'\',char(FolderNames(j))];
    
    files = dir(dataset_path);
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    subFolderNames = {subFolders(3:end).name} ;

    for k = 1:length(subFolderNames)
        txtfiles = dir(fullfile(dataset_path,'\',char(subFolderNames(k)),'\Output\PD\', '*.txt'));
        numTextFiles = numel(txtfiles);
        
        for i = 1:numTextFiles
            txt_path = [dataset_path,'\',char(subFolderNames(k)),'\Output\PD\',txtfiles(i).name];
            fileID = fopen(txt_path, 'r');
            data = textscan(fileID, '%f %f %f');  % The file contains three columns of numeric data; x coordinates of the fiducials, y coordinates of the fiducials and the fiducial number
            fclose(fileID);
            dataArray = [data{1}, data{2}, data{3}]; 
            
            if length(dataArray) == 17
                sortedArray = sortrows(dataArray, 3);
                y = sortedArray(:,2);

                newFirstRow = [dataset_path,'\',char(subFolderNames(k)),'\', txtfiles(i).name(1:2)];
                originalCellArray = num2cell(y);
                insertedCellArray = vertcat({newFirstRow}, originalCellArray);
                
                result = horzcat(result, insertedCellArray); % y-coordinates
            else
                continue
            end
        end
    end
    y_coordinates = cell2mat(result(2:end,:));
    y_avg = [];
    for m = 1:17
        avg = mean(y_coordinates(m,:));
        y_avg(m) = avg;
    end
    y_avg = y_avg';
end   

names = result(1,:);

res_y = [];
for m = 1:size(y_coordinates, 2) % code runs for all the images
    A = y_coordinates(:, m);
    B = y_avg;
    cvx_begin
            variable y(17)
            minimize(norm(A .* y - B)) % weights are obtained for x coordinates of the fiducials in a particular image
    cvx_end
    FirstRow = (names(m));

    originalcellarray = num2cell(y);
    insertedcellarray = vertcat({FirstRow}, originalcellarray);
    
    res_y = horzcat(res_y, insertedcellarray); % concatenated weights for all the images
end


