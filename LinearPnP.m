function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose translation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numerically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

C = zeros(3,1);
R = zeros(3,3);
x = [x, ones(size(x,1), 1)];

x = x';

% plot(xt(:,1))
% figure(2)
% plot(x(1,:))
% figure(3)
% xt(:,1);
% x(1,:);
% a = xt(:,1)' - x(1,:);
% plot(a)
% pause()

x = inv(K) * x;


x1 = X(:,1);
y1 = X(:,2);
z1 = X(:,3);

u = x(1,:);
v = x(2,:);


A = zeros(2*size(x,1), 12);
k = 1;

for i=1:size(x,1)
    
    A(k,:) = [0,0,0,0,-x1(i), -y1(i), -z1(i), -1, x1(i)*v(i), y1(i)*v(i),...
        z1(i)*v(i), v(i)];
    k = k + 1;
    
    A(k,:) = [-x1(i), -y1(i), -z1(i), -1, 0, 0, 0, 0, x1(i)*u(i), y1(i)*u(i),...
        z1(i)*u(i), u(i)];
    k = k + 1;
    
    
end

[U, S, V] = svd(A);

% V(:,end)
% P = reshape(V(:,end), [3,4])
% t = P(:,end);
% R = P(1:3, 1:3);

V = V(:,end);
P = [V(1:4)'; V(5:8)'; V(9:12)'];
t = P(:,end);
R = reshape(P(1:3,1:3), [3,3])';

[U, D, V] = svd(R);

if det(U*V')>.9 && det(U*V')<1.1
    R = U*V';
    t = t./D(1,1);
    
elseif det(U*V')>-1.1 && det(U*V')<-.9
    R = -U*V';
    t = -t./D(1,1);
    
else
    error=1
    
end
det(R)
%R = [R(1,1); R(2,2); R(3,3)]

C = -R'*t

% theta = acos(.5*(R(1,1)+R(2,2)+R(3,3)-1));
% e1 = (R(3,2)-R(2,3))/2/sin(theta);
% e2 = (R(1,3)-R(3,1))/(2*sin(theta));
% e3 = (R(2,1)-R(1,2))/(2*sin(theta));
% 
% R = [e1,e2,e3]

end
