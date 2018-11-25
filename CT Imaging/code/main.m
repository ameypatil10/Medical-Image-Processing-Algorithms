%% Question Number 1:

clear all; clc;
phantom = phantom(128);
[X_img,Y_img] = meshgrid(1:128);
[X_interp,Y_interp] = meshgrid(1:0.711:128);
phantom180 = interp2(X_img,Y_img,phantom,X_interp,Y_interp);
dt = 1; dtheta = 1; 
t = floor(size(phantom180,1)/2) - mod(floor(size(phantom180,1)/2),dt);
Rf1 = myRadonTrans(phantom180,dt,dtheta,0.5);
Rf2 = myRadonTrans(phantom180,dt,dtheta,1);
Rf3 = myRadonTrans(phantom180,dt,dtheta,3);

figure(); show_colormap(Rf1,-t:dt:t,0:dtheta:179);title('Radon Transform for ds = 0.5 pixel width');hold off;
figure(); show_colormap(Rf2,-t:dt:t,0:dtheta:179);title('Radon Transform for ds = 1 pixel width');hold off;
figure(); show_colormap(Rf3,-t:dt:t,0:dtheta:179);title('Radon Transform for ds = 3 pixel width');hold off;

figure();
subplot(1,2,1);plot(-t:dt:t,Rf1(:,1));title('1D function plots for the Radon-transform for theta = 0 & ds = 0.5');
subplot(1,2,2);plot(-t:dt:t,Rf1(:,floor(size(Rf1,2)/2)));title('1D function plots for the Radon-transform for theta = 90 & ds = 0.5');
hold off;figure();
subplot(1,2,1);plot(-t:dt:t,Rf2(:,1));title('1D function plots for the Radon-transform for theta = 0 & ds = 1');
subplot(1,2,2);plot(-t:dt:t,Rf2(:,floor(size(Rf2,2)/2)));title('1D function plots for the Radon-transform for theta = 90 & ds = 1');
hold off;figure();
subplot(1,2,1);plot(-t:dt:t,Rf3(:,1));title('1D function plots for the Radon-transform for theta = 0 & ds = 3');
subplot(1,2,2);plot(-t:dt:t,Rf3(:,floor(size(Rf3,2)/2)));title('1D function plots for the Radon-transform for theta = 90 & ds = 3');
hold off;

%% Question Number 2

clear all; clc;
f = phantom(256);
my_theta = 0:3:177;
Rf = radon(f,my_theta);
BRf = iradon(Rf,my_theta,'linear','none',1,256);
w_max = 1;
RRMSE = @(A,B) [sqrt(sum(sum((A-B) .^ 2))) ./ sqrt(sum(sum(A .^ 2)))];

figure(); show_colormap(Rf);title('Radon Transform for phantom 256');hold off;
figure(); show_colormap(BRf);title('Back Projection of Radon Transform for phantom 256');hold off;
%%
I_Ram_Lak1 = myFilter(Rf, my_theta, 'Ram-Lak', w_max);
figure(); show_colormap(I_Ram_Lak1);title('Radon Transform for Ram-Lak & L = w_{max}');hold off;
I_Shepp_Logan1 = myFilter(Rf, my_theta, 'Shepp-Logan', w_max);
figure(); show_colormap(I_Shepp_Logan1);title('Radon Transform for Shepp-Logan & L = w_{max}');hold off;
I_Cosin1 = myFilter(Rf, my_theta, 'Cosine', w_max);
figure(); show_colormap(I_Cosin1);title('Radon Transform for Cosine & L = w_{max}');hold off;

I_Ram_Lak2 = myFilter(Rf, my_theta, 'Ram-Lak', w_max/2);
figure(); show_colormap(I_Ram_Lak2);title('Radon Transform for Ram-Lak & L = w_{max}/2');hold off;
I_Shepp_Logan2 = myFilter(Rf, my_theta, 'Shepp-Logan', w_max/2);
figure(); show_colormap(I_Shepp_Logan2);title('Radon Transform for Shepp-Logan & L = w_{max}/2');hold off;
I_Cosin2 = myFilter(Rf, my_theta, 'Cosine', w_max/2);
figure(); show_colormap(I_Cosin2);title('Radon Transform for Cosine & L = w_{max}/2');hold off;

mask1 = fspecial ('gaussian', 11, 1);
mask5 = fspecial ('gaussian', 51, 5);
S1 = conv2 (f, mask1, 'same');
S2 = conv2 (f, mask5, 'same');

Rf1 = radon(S1,my_theta);
Rf2 = radon(S2,my_theta);

R0 = myFilter(Rf,my_theta,'Ram-Lak',w_max);
R1 = myFilter(Rf1,my_theta,'Ram-Lak',w_max);
R2 = myFilter(Rf2,my_theta,'Ram-Lak',w_max);

RRMSE0 = RRMSE(f,R0);
RRMSE1 = RRMSE(S1,R1);
RRMSE2 = RRMSE(S2,R2);
figure();
subplot(1,2,1);show_colormap(f);title('Image without adding blur');
subplot(1,2,2);show_colormap(R0);title('Filtered Back projection of Blurless Image');hold off;
figure(); 
subplot(1,2,1);show_colormap(S1);title('Gaussain blurred image 1');
subplot(1,2,2);show_colormap(R1);title('Ram-Lak filtered Back Projection for Gaussian Blurred Image 1');hold off;
figure(); 
subplot(1,2,1);show_colormap(S2);title('Gaussain blurred image 2');
subplot(1,2,2);show_colormap(R2);title('Ram-Lak filtered Back Projection for Gaussian Blurred Image 2');hold off;

RRMSE0w = []; RRMSE1w = []; RRMSE5w = [];
freqs = linspace(-1, 1, size(Rf,1));
for l = (1/length(freqs)):(1/length(freqs)):1
    R0 = myFilter(Rf,my_theta,'Ram-Lak',l);
    R1 = myFilter(Rf1,my_theta,'Ram-Lak',l);
    R2 = myFilter(Rf2,my_theta,'Ram-Lak',l);
    RRMSE0w = [RRMSE0w,RRMSE(f,R0)];
    RRMSE1w = [RRMSE1w,RRMSE(S1,R1)];
    RRMSE5w = [RRMSE5w,RRMSE(S2,R2)];
end
figure(); plot(freqs,RRMSE0w); title('RRMSE values vs rectangle Length for blurless image'); xlabel('L');ylabel('RRMES');
figure(); plot(freqs,RRMSE1w); title('RRMSE values vs rectangle Length for gaussian blurred image 1'); xlabel('L');ylabel('RRMES');
figure(); plot(freqs,RRMSE5w); title('RRMSE values vs rectangle Length for gaussian blurred image 2'); xlabel('L');ylabel('RRMES');

%% Question Number 3

clear all; clc;
load('../data/CT_Chest.mat'); image_chest = imageAC;
load('../data/myPhantom.mat'); image_phantom = imageAC;

RRMSE = @(X,Y) [sqrt(sum(sum((X-Y) .^2))) / sqrt(sum(sum(X .^ 2)))];
RRMSE_chest = []; RRMSE_chest_best=1;
RRMSE_phantom = []; RRMSE_phantom_best=1;
for theta = 0:180
   my_theta = theta:1:theta+149;
   R_chest = radon(image_chest, my_theta);
   R_phantom = radon(image_phantom, my_theta);
   I_chest = iradon(R_chest, my_theta,'linear', 'Cosine', 1, size(image_chest,1));
   I_phantom = iradon(R_phantom, my_theta,'linear', 'Cosine', 1, size(image_phantom,1));
   r1 = RRMSE(image_chest,I_chest); r2 = RRMSE(image_phantom,I_phantom);
   if r1 < RRMSE_chest_best
       I_chest_best = I_chest; RRMSE_chest_best = r1;
   end
   if r2 < RRMSE_phantom_best
       I_phantom_best = I_phantom; RRMSE_phantom_best = r2;
   end
   RRMSE_chest = [RRMSE_chest,RRMSE(image_chest,I_chest)]; clear r1; clear r2;
   RRMSE_phantom = [RRMSE_phantom,RRMSE(image_phantom,I_phantom)];
end
figure(); plot(0:180, RRMSE_chest); title('RRMSE-chest');
figure(); plot(0:180, RRMSE_phantom); title('RRMSE-phantom');
figure(); subplot(1,2,1); show_colormap(I_chest_best); title('I chest best');
subplot(1,2,2); show_colormap(I_phantom_best); title('I phantom best');