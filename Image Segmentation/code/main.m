%% Part 1
clear all; clc;
load('../data/assignmentSegmentBrain.mat')

sigma = 1;
gauss2D = @(x, y, x_mean, y_mean) [exp(-0.5 * ((x - x_mean) .^ 2 + (y - y_mean) .^ 2) / (sigma*sigma))];
% The values of q vs the sum of absolute of the residual image
% 1.5 => 9.21, 1.7 => 309.9288, 1.8 => 303.9482, 2 => 303.6365, 2.5 => 381.3746
% 1.9 => 301.8408 1.85 => 302.2939
Y = imageData .* imageMask;
q = 1.95; K = 3; lim = 10; len = size(Y,1);
W = zeros(2*lim-1, 2*lim-1);
C = zeros(3,1);
for x = 1:size(W,1)
    for y = 1:size(W,2)
        W(x, y) = gauss2D(x,y,lim,lim);
    end
end
W = W ./ sum(sum(W));
B = ones(size(Y));
D = ones(len,len,K);
U = (ones(256, 256, K) / K);temp = U(:,:,1).*imageMask;temp(1,1)=1;
figure();imshow(temp,[]);title('Initial Membesship Image');
J = @(U,D) [sum(sum(sum(D .* (U .^ q))))];
[x,m] = kmeans(reshape(Y,256*256,1),K+1);
k=1; old_cost = inf;
for i = 1:K+1
    if m(i) > 0.01
       C(k) = m(i);
       k = k+1;
    end
end
new_cost = J(U,D);
Cost = [];

while new_cost < old_cost
    %Code to find the optimal values of memberships
    for k = 1:K
       D(:,:,k) = (Y .^ 2) - (2 * C(k) * Y .* conv2(B, W, 'same') .* imageMask) + ((C(k) ^ 2) * conv2(B.^2, W, 'same') .* imageMask);
    end
    temp = (1 ./ D) .^ (1 / (q-1));
    temp(isnan(temp))=0;
    U = temp ./ repmat(sum(temp,3),1,1,3);
    U(isnan(U))=0;
    % Code to find the optimal values of means
    for k = 1:K
       C(k) = sum(sum((U(:,:,k) .^ q) .* Y .* conv2(B, W, 'same') .* imageMask)) / sum(sum((U(:,:,k) .^ q) .* conv2(B.^2, W, 'same') .* imageMask));
    end
    % Code to find the optimal values of Bias
    num = 0; den = 0;
    for k = 1:K
       num = num + C(k) * (U(:,:,k) .^ q);
       den = den + C(k) * C(k) * (U(:,:,k) .^ q);
    end
    num = num .* Y .* imageMask;
    den = den .* imageMask;
    num(isnan(num)) = 0; den(isnan(den)) = 0;
    B = (conv2(num, W, 'same') .* imageMask) ./ (conv2(den, W, 'same') .* imageMask);
    B(isnan(B)) = 0;
    old_cost = new_cost;
    new_cost = J(U,D)
    Cost = [Cost new_cost];
%     A = U(:,:,1)*C(1)+U(:,:,2)*C(2)+U(:,:,3)*C(3);
%     diff = sum(sum(abs(imageData .* imageMask - (A .* B))))
end

figure();imshow(imageData,[]);title('Given Corrupted Data');
figure();plot(1:size(Cost,2), Cost);xlabel('Iteration');ylabel('Cost value');
title('The Value of Cost Function at each iteration');
A = U(:,:,1)*C(1)+U(:,:,2)*C(2)+U(:,:,3)*C(3);
R = imageData .* imageMask - (A .* B);
figure();imshow(abs(A),[]);title('The Bias Removedd Image - A');
figure();imshow(abs(R),[]);title('Residual Image - R');
figure();imshow(abs(U(:,:,1)),[]);title('Membership values for class 1');
figure();imshow(abs(U(:,:,2)),[]);title('Membership values for class 2');
figure();imshow(abs(U(:,:,3)),[]);title('Membership values for class 3');
figure();imshow(abs(B),[]);title('Bias Field as Image');

%% Part 2
p2