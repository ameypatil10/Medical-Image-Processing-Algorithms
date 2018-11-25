function X = normalizeData(Data,entries)
    
    temp = sqrt(sum(Data .^ 2,1));
    X = Data ./ repmat(temp,entries,1);
    
end