function err = loadErr(filePrefix)
    fileIn = fopen([filePrefix 'err.txt'], 'r');
    err = fscanf(fileIn, '%f ');
    fclose(fileIn);    
end