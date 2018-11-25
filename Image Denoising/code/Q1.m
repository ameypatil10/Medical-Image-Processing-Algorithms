clear all;
clc;

load('../data/assignmentImageDenoisingPhantom.mat');

imageNoiseless = abs(reshape(imageNoiseless,256*256,1));  % The Noiseless Data
X = abs(reshape(imageNoisy,256*256,1));
imageNoisy = X;                                           % The Observed Noisy Data

sigma = 1;   % Noise Level

likelihood_cost_function = @(X) [ (sum(X - imageNoisy) .^ 2) / (sigma .^ 2) ];
likelihood_grad_function = @(X) [ 2 * (X - imageNoisy) / (sigma .^ 2) ];

left_neighbor = @(X) [ reshape(circshift(reshape(X,256,256),[0,1]),256*256,1) ];
right_neighbor = @(X) [ reshape(circshift(reshape(X,256,256),[0,-1]),256*256,1) ];
top_neighbor = @(X) [ reshape(circshift(reshape(X,256,256),[1,0]),256*256,1) ];
bottom_neighbor = @(X) [ reshape(circshift(reshape(X,256,256),[-1,0]),256*256,1) ];

add_prior_cost_component = @(X,cost_component) [ cost_component(X,left_neighbor(X)) + cost_component(X,right_neighbor(X)) + cost_component(X,top_neighbor(X)) + cost_component(X,bottom_neighbor(X))];
add_prior_grad_component = @(X,grad_component) [ grad_component(X,left_neighbor(X)) + grad_component(X,right_neighbor(X)) + grad_component(X,top_neighbor(X)) + grad_component(X,bottom_neighbor(X))];

RRMSE = @(A,B) [sqrt(sum((A - B) .^ 2)) / sqrt(sum(A .^ 2))];
RRMSE0 = RRMSE(imageNoiseless,X)
%% Prior Number 1
    % alpha = 0.04 RRMSE = 0.2804

prior1_cost_component = @(X,S) [ sum((X - S) .^ 2) ];
prior1_grad_component = @(X,S) [ 2 * (X - S) ];

cost_function1 = @(X,alpha) [ (1 - alpha) * likelihood_cost_function(X) + alpha * add_prior_cost_component(X,prior1_cost_component) ];
grad_function1 = @(X,alpha) [ (1 - alpha) * likelihood_grad_function(X) + alpha * add_prior_grad_component(X,prior1_grad_component) ];

alpha1 = 0.04;
right_alpha1 = min(alpha1 * 1.2,0.9999);
left_alpha1 = max(0.0001,alpha1 * 0.8);

figure(1);hold on;title('Objective Function for Quadration Function Prior');xlabel('Iterations');ylabel('Cost Function');
X1 = gradient_descent(cost_function1, grad_function1, X, alpha1, 1);
X1_alpha_left = gradient_descent(cost_function1, grad_function1, X, left_alpha1, 0);
X1_alpha_right = gradient_descent(cost_function1, grad_function1, X, right_alpha1, 0);

RRMSE1 = RRMSE(imageNoiseless,X1)
RRMSE1_alpha_left = RRMSE(imageNoiseless,X1_alpha_left);
RRMSE1_alpha_right = RRMSE(imageNoiseless,X1_alpha_right);

%% Prior Number 2
    % alpha = 0.9800 gamma = 0.0016 RRMSE = 0.2358
               
gamma = 0.0016;
alpha2 = 0.9800;
right_alpha2 = min(alpha2 * 1.2,0.9999);
left_alpha2 = max(0.0001,alpha2 * 0.8);

prior2_cost_component = @(X,S) [sum((0.5 * ((X-S) .* (X-S)) .* (abs(X-S) <= gamma)) + (gamma * abs(X-S) - 0.5*gamma*gamma) .* (abs(X-S) > gamma))];
prior2_grad_component = @(X,S) [((X-S) .* (abs(X-S) <= gamma)) + (gamma * sign(X-S)) .* (abs(X-S) > gamma)];

cost_function2 = @(X,alpha) [ (1 - alpha) * likelihood_cost_function(X) + alpha * add_prior_cost_component(X,prior2_cost_component) ];
grad_function2 = @(X,alpha) [ (1 - alpha) * likelihood_grad_function(X) + alpha * add_prior_grad_component(X,prior2_grad_component) ];
figure(2);hold on;title('Objective Function for Huber Function Prior');xlabel('Iterations');ylabel('Cost Function');
X2 = gradient_descent(cost_function2, grad_function2, X, alpha2, 1);
X2_alpha_left = gradient_descent(cost_function2, grad_function2, X, left_alpha2, 0);
X2_alpha_right = gradient_descent(cost_function2, grad_function2, X, right_alpha2, 0);

RRMSE2 = RRMSE(imageNoiseless,X2)
RRMSE2_alpha_left = RRMSE(imageNoiseless,X2_alpha_left);
RRMSE2_alpha_right = RRMSE(imageNoiseless,X2_alpha_right);

%% Prior Number 3
    % alpha = 0.9900    gamma = 0.0008    RRMSE = 0.2362   
gamma = 0.0008;
alpha3 = 0.9900;
right_alpha3 = min(alpha3 * 1.2,0.9999);
left_alpha3 = max(0.0001,alpha3 * 0.8);

prior3_cost_component = @(X,S) [ sum(gamma * abs(X-S) - gamma*gamma * log(1 + (abs(X-S) / gamma)))];
prior3_grad_component = @(X,S) [ gamma * (X-S) ./ (abs(X-S) + gamma) ];

cost_function3 = @(X,alpha) [ (1 - alpha) * likelihood_cost_function(X) + alpha * add_prior_cost_component(X,prior3_cost_component) ];
grad_function3 = @(X,alpha) [ (1 - alpha) * likelihood_grad_function(X) + alpha * add_prior_grad_component(X,prior3_grad_component) ];
figure(3);hold on;title('Objective Function for 3rd Function Prior');xlabel('Iterations');ylabel('Cost Function');
X3 = gradient_descent(cost_function3, grad_function3, X, alpha3, 1);
X3_alpha_left = gradient_descent(cost_function3, grad_function3, X, left_alpha3, 0);
X3_alpha_right = gradient_descent(cost_function3, grad_function3, X, right_alpha3, 0);

RRMSE3 = RRMSE(imageNoiseless,X3)
RRMSE3_alpha_left = RRMSE(imageNoiseless,X3_alpha_left);
RRMSE3_alpha_right = RRMSE(imageNoiseless,X3_alpha_right);


%% Analysis
figure(10);imshow(reshape(imageNoisy,256,256));title('Noisy Image');
figure(11);imshow(reshape(X1,256,256));title('Prior1 Denoised Image');
figure(12);imshow(reshape(X2,256,256));title('Prior2 Denoised Image');
figure(13);imshow(reshape(X3,256,256));title('Prior3 Denoised Image');
figure(14);imshow(reshape(imageNoiseless,256,256));title('Noiseless Image');
