function saveDistEstRes(filePref, R1, t1, f1, k1)
    fileOut = fopen([filePref 'res.txt'], 'w');
    fprintf(fileOut, '%f ', R1);
    fprintf(fileOut, '%f ', t1);
    fprintf(fileOut, '%f ', f1);
    fprintf(fileOut, '%f ', k1);
    fclose(fileOut);
end