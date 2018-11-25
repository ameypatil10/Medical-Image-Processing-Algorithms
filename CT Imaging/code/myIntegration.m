function ans = myIntegration(I,t,theta,ds)
    boundary = ceil(size(I,1)/2); ans=0;
    f = @(x,y,boundary) [I(boundary+x,boundary-y)];
    s_lim = sqrt(2*(boundary^2) - t^2);
    for s = -s_lim:ds:s_lim
        x = floor(t * cosd(theta) + s * sind(theta));
        y = floor(t * sind(theta) - s * cosd(theta));
        if(abs(x) < boundary && abs(y) < boundary)
            ans = ans + f(x, y,boundary);
        end
    end
end
