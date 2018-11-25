function Xrotated = rotateToAlign(dim,numOfPoints,numOfPointSets,X,refranceSet)

    Xrotated = zeros(dim,numOfPoints,numOfPointSets);
    
    if dim == 2
        for i = 1:numOfPointSets
            x = reshape(X(:,i),dim,numOfPoints);
            [U,S,V] = svd(x * refranceSet');
            R = V * [1 0; 0 det(V * U')] * U';
            Xrotated(:,:,i) = R * x;
        end
    else 
        for i = 1:numOfPointSets
            x = reshape(X(:,i),dim,numOfPoints);
            [U,S,V] = svd(x * refranceSet');
            R = V * [1 0 0;0 1 0;0 0 det(V * U')] * U';
            Xrotated(:,:,i) = R * x;
        end
end