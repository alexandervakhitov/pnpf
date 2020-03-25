function BPans = formBPAns(BP, bvec, fest, tarPtNum)
    BPans = zeros(3, tarPtNum);    
    for kerInd = 1:size(BP, 3)
        BPans = BPans + BP(:, :, kerInd)*bvec(kerInd);
    end
    BPans(3, :) = BPans(3, :) * fest;
end