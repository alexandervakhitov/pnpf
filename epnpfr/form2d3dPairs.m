function [p3d, p2d] = form2d3dPairs(filePath, pts3d)
    fIn = fopen(filePath, 'r');
    dataRaw = fscanf(fIn, '%f ');
    data = reshape(dataRaw, 3, length(dataRaw)/3);
    p2d = data(2:3, :);
    p3d = [];
    cx = 576*2;
    cy = 384*2;
    for i = 1:size(data, 2)
        p3d = [p3d pts3d(data(1, i)+1, :)'];
        p2d(:, i) = p2d(:, i) - [cx; cy];
    end
end