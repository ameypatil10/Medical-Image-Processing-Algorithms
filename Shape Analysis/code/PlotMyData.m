function PlotMyData(pointSets, numOfPoints, numOfPointSets, TriangleIndex)
    
    switch nargin
        case 3
            for i = 1:numOfPointSets
                scatter(pointSets(1,:,i),pointSets(2,:,i),5,'filled');
                hold on;
            end
        case 4
            for i = 1:numOfPointSets
                scatter3(pointSets(1,:,i),pointSets(2,:,i),pointSets(3,:,i),5,'filled');
                hold on;
            end
    end
    hold off;
end
            