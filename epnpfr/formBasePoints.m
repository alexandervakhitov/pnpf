function [BP] = formBasePoints(Ker, tarPtNum)
    BP = zeros(3, tarPtNum, size(Ker, 2));
    for colInd = 1:size(Ker, 2)
        BP(:, :, colInd) = reshape(Ker(:, colInd), 3, tarPtNum);
    end
end
