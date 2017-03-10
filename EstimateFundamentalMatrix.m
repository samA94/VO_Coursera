function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

u1 = x1(:,1);
u2 = x2(:,1);
v1 = x1(:,2);
v2 = x2(:,2);

A = zeros(size(u1,1), 9);

for i=1:size(x1,1)
    A(i,:) = [u1(i)*u1(i), u1(i)*v2(i), u1(i), v1(i)*u2(i), v1(i)*v2(i),...
        v1(i), u2(i), v2(i), 1];
end

[U, S, V] = svd(A);
size(S);
V = V';

x = V(:,end);

%x(9) = 0;

F = reshape(x, 3, 3);

[U, S, V] = svd(F);
S(3,3) = 0;

F = U*S*V';

%F = F./norm(F)

% S = [diag(x, 0); zeros(size(S, 1)-9, 9)];
% 
% A = U*S*V;
% 
% [U, S, V] = svd(A);
% 
% x = V(:,end)
% 
% F = reshape(x, 3, 3);
F = F./norm(F)';

end
