function ER=registration_check(P,W_dist_pts,pos,r,path_d,rewarp)
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


