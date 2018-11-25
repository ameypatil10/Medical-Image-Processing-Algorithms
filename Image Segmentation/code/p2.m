clc;
clear;
load('../data/assignmentSegmentBrainGmmEmMrf.mat');
[X, NU] = kmeans(imageData(:).*imageMask(:),4);
while(min(X)~=X(1)) 
    [X, NU] = kmeans(imageData(:).*imageMask(:),4);
end
Y = imageData(:).*imageMask(:);
NU = NU([2,3,4]);
index=find(NU==min(NU));
X(X==1)=index+1;
X = X - 1;
Beta = 1;
X = X.*imageMask(:);
VAR = [var(nonzeros((X==1).* imageData(:))),var(nonzeros((X==2).*imageData(:))),var(nonzeros((X==3).*imageData(:)))];


LOGPOSTERIOR = [];

for i = 1:100
    [GAM, X] = updateGAM(Y, NU, VAR, X, Beta, imageMask(:));
    NU = [updateNUk(GAM(:,1),Y,imageMask(:)),updateNUk(GAM(:,2),Y,imageMask(:)),updateNUk(GAM(:,3),Y,imageMask(:))];
    VAR = [updateVARk(GAM(:,1),Y,NU(1),imageMask(:)),updateVARk(GAM(:,2),Y,NU(2),imageMask(:)),updateVARk(GAM(:,3),Y,NU(3),imageMask(:))];
    
    LOGPOSTERIOR = [LOGPOSTERIOR,LOGP(X,Y,NU,VAR,Beta)];
    
end;

figure();
P = (X==1).*Y;
P = reshape(P,256,256);
imshow(P);
figure();
P = (X==2).*Y;
P = reshape(P,256,256);
imshow(P);
figure();
P = (X==3).*Y;
P = reshape(P,256,256);
imshow(P);

display(LOGPOSTERIOR);