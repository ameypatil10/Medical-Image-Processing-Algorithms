function X = setCentroidToOrigin(pointSets, numOfPoints, numOfPointSets, dim)

    Mean = mean(pointSets,2);
    X = reshape(pointSets - repmat(Mean,1,numOfPoints), numOfPoints*dim, numOfPointSets);
    
end