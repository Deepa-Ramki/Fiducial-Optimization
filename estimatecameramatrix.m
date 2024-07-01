function [camMatrix, reprojectionErrors, x, y] = estimatecameramatrix(imagePoints, worldPoints)

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


