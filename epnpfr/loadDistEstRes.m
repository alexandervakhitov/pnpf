function [R1, t1, f1, k1] = loadDistEstRes(filePref)
    fileIn = fopen([filePref 'res.txt'], 'r');
    dataRaw = fscanf(fileIn, '%f ');
    fclose(fileIn);
    R1 = reshape(dataRaw(1:9), 3, 3);
    pos = 10;
    t1 = dataRaw(pos:pos+2);
    pos = pos+3;
    f1 = dataRaw(pos);
    k1 = dataRaw(pos+1);
end