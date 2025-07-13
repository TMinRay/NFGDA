function groups =  cc2groups(asize,cc)
    groups = zeros(asize);
    for id = 1:cc.NumObjects
        groups(cc.PixelIdxList{id})=id;
    end
end