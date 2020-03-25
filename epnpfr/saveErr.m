function saveErr(filePrefix, err)
    fileOut = fopen([filePrefix 'err.txt'], 'w');
    fprintf(fileOut, '%f ', err);
    fclose(fileOut);
end
