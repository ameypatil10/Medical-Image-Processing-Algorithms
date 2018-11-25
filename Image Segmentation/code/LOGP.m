function [ ans ] = LOGP( X, Y, NU, VAR, Beta )

    Xsq = reshape(X,256,256);
    XL = circshift(Xsq,1,2);
    XR = circshift(Xsq,-1,2);
    XU = circshift(Xsq,1,1);
    XD = circshift(Xsq,-1,1);
    
    Prior1 = reshape((Xsq==1).*((XL~=1).*(XL~=0)+(XR~=1).*(XR~=0)+(XU~=1).*(XU~=0)+(XD~=1).*(XD~=0)),65536,1);
    Prior2 = reshape((Xsq==2).*((XL~=2).*(XL~=0)+(XR~=2).*(XR~=0)+(XU~=2).*(XU~=0)+(XD~=2).*(XD~=0)),65536,1);
    Prior3 = reshape((Xsq==3).*((XL~=3).*(XL~=0)+(XR~=3).*(XR~=0)+(XU~=3).*(XU~=0)+(XD~=3).*(XD~=0)),65536,1);
    Prior = -Beta.*(Prior1 + Prior2 + Prior3);
    
    F = @(y,u,v) -((y-u).^2) ./ (2*(v));
    
    Likelihood = (X==1).*( ( (Y-NU(1)).^2 )./( 2*VAR(1) ) ) + (X==2).*( ( (Y-NU(2)).^2 )./( 2*VAR(2) ) ) + (X==3).*( ( (Y-NU(3)).^2 )./( 2*VAR(3) ) );
    
    ans = sum( Prior ) + sum( Likelihood ) -0.5*sum(X==1)*log(VAR(1)) -0.5*sum(X==2)*log(VAR(2)) -0.5*sum(X==3)*log(VAR(3)); 
    

end

