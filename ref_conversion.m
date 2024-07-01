function [Cmm_pts,ref_dist_pts,ref2cam]=ref_conversion(position,W,distpts,C2R,C2D,path_d)


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


