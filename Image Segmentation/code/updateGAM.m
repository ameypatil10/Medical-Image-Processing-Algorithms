function [ GAM, Xnew ] = updateGAM( Y, NU, VAR, X, Beta, M )

    Xsq = reshape(X,256,256);
    XL = circshift(Xsq,1,2);
    XR = circshift(Xsq,-1,2);
    XU = circshift(Xsq,1,1);
    XD = circshift(Xsq,-1,1);
    
    F = @(y,u,v) -((y-u).^2) ./ (2*(v));
    GY1 = F(Y,NU(1),VAR(1));
    GY2 = F(Y,NU(2),VAR(2));
    GY3 = F(Y,NU(3),VAR(3));
    
    Prior1 = -Beta.*reshape((XL~=1).*(XL~=0)+(XR~=1).*(XR~=0)+(XU~=1).*(XU~=0)+(XD~=1).*(XD~=0),65536,1);
    Prior2 = -Beta.*reshape((XL~=2).*(XL~=0)+(XR==2).*(XR~=0)+(XU~=2).*(XU~=0)+(XD~=2).*(XD~=0),65536,1);
    Prior3 = -Beta.*reshape((XL~=3).*(XL~=0)+(XR==3).*(XR~=0)+(XU~=3).*(XU~=0)+(XD~=3).*(XD~=0),65536,1);
    
    Label1 = exp(GY1 + Prior1)./(sqrt(2*pi)*sqrt(VAR(1)));
    Label2 = exp(GY2 + Prior2)./(sqrt(2*pi)*sqrt(VAR(2)));
    Label3 = exp(GY3 + Prior3)./(sqrt(2*pi)*sqrt(VAR(3)));
    
    Total = Label1 + Label2 + Label3;

    Xnew = ((Label1 >= Label2).*(Label1 >= Label3).*1 + (Label2 >= Label3).*(Label2 > Label1).*2 + (Label3 > Label2).*(Label3 > Label1).*3);    
    Xnew(M==0)=0;
    Label1(M==0)=0;
    Label2(M==0)=0;
    Label3(M==0)=0;
    GAM = horzcat(Label1 ./ Total, Label2 ./ Total, Label3 ./ Total);
    
end

