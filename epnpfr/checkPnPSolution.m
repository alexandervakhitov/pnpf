function ret = checkPnPSolution(currentSolution, pnpOpts)
    ret = 0;
    if (currentSolution.avgErr < 0)
        return;
    end
    if (currentSolution.avgErr > pnpOpts.errThr)
        return;
    end
    if (currentSolution.f > pnpOpts.fMax)
        return;
    end
    ret = 1;
    return;
end