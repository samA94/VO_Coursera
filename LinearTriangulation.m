function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
%% LinearTriangulation
% Find 3D positions of the point correspondences using the relative
% position of one camera from another
% Inputs:
%     C1 - size (3 x 1) translation of the first camera pose
%     R1 - size (3 x 1) rotation of the first camera pose
%     C2 - size (3 x 1) translation of the second camera
%     R2 - size (3 x 1) rotation of the second camera pose
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Outputs: 
%     X - size (N x 3) matrix whose rows represent the 3D triangulated
%       points


%P1 = K * [eye(3), -C1];
%P2 = K * [eye(3), -C2];

P1 = K * [R1, -C1];
P2 = K * [R2, -C2];

A = zeros(6,4);
X = zeros(size(x1,1), 4);

for i=1:size(x1,1)
    skew1 = Vec2Skew([x1(i,:)'; 1]);
    skew2 = Vec2Skew([x2(i,:)'; 1]);
    
    A = [skew1*P1; skew2*P2];
    [U, D, V] = svd(A);
    
    X(i,:) = V(:,end)/V(end,end);
    
end

X = X(:,1:3);

end

