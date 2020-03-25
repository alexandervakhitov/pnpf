function finalSolution = findBetterSolution(newSolution, currentSolution)
    finalSolution = currentSolution;
    if (newSolution.avgErr < 0)
        return;
    end    
    if (currentSolution.avgErr < 0)
        finalSolution = newSolution;
        return;
    end
    if (newSolution.avgErr < currentSolution.avgErr)
        finalSolution = newSolution;
    end
    return;
end