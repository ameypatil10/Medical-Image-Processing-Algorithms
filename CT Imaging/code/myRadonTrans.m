function Rf = myRadonTrans(f,dt,dtheta,ds)
    if ~exist('dt','var')
          dt = 5;
    end
    if ~exist('dtheta','var')
          dtheta = 5;
    end
    if ~exist('ds','var')
          ds = 3;
    end
    t = floor(size(f,1)/2) - mod(floor(size(f,1)/2),dt);
    t = -t:dt:t;
    theta = 0:dtheta:179;
    Rf = zeros(size(t,2),size(theta,2));
    for i = 1:size(t,2)
        for j = 1:size(theta,2)
            Rf(i,j) = myIntegration(f,t(i),theta(j),ds);
        end
    end
end
