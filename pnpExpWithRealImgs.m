function pnpExpWithRealImgs
%sfm experiment
%we made 3 shots of the scene for 6 focal length settings (in folder
%sfmData) and built a 3D model from two 35 mm shots

%the script solve PnPf for all the images
%the matching results are saved in text files
%the script saves the reprojection errors and the focal lengths found by
%the EPnPfR and GPnPf+GN methods to the tab.txt file
%each file line has format: true focal length (from EXIF) & EPnPfR result
%from image 1 & GPnPf+GN result from image 1 & ... & EPnPfR result
%from image 3 & GPnPf+GN result from image 3 \\ \n
    addpath(genpath('toolbox\rpnp1.0'));
    addpath(genpath('toolbox\DLT'));
    addpath(genpath('toolbox\UPnP'))
    addpath(genpath('toolbox'));
    addpath(genpath('epnpfr'));
    addpath(genpath('levmar\matlab'))
    runMatching    
    prepareTableRealPnpExps
end

function runMatching
    debugFolder = 'sfmData\';
    ptsFileName = [debugFolder 'pts3d.mat'];
    load(ptsFileName);    
    showPts3dCoordsFileName = 'sfmData\showPts3d.txt';
    fileInShowPts = fopen(showPts3dCoordsFileName, 'r');
    showPtsRaw = fscanf(fileInShowPts, '%f ');
    fclose(fileInShowPts);
    showPts3d = reshape(showPtsRaw, 3, length(showPtsRaw)/3);
    fpix = 3491;
    cx = 576*2;
    cy = 384*2;    
    K = [fpix 0 cx;
        0 fpix cy;
        0 0 1];
    
    folderSuffs = {'18', '24', '35', '50', '75', '105'};
    folderSuffs = folderSuffs(1:end);
    for fInd = 1:length(folderSuffs)
        activeFolder = ['sfmData\' char(folderSuffs{fInd}) '\'];       
        fileLogIn = fopen([activeFolder 'fileLog.txt'], 'r');        
        thr = 3;
        for i = 0:2
            filePath = [activeFolder int2str(i) '.txt'];
            imgPath = fgetl(fileLogIn);
            [p3d, p2d] = form2d3dPairs(filePath, pts3d);
    %        [f1,R1,t1, k1est] = 
            [fFin RFin tFin goodInds] = pnpfmyRANSAC(p3d, p2d, 1e3, thr);                        
            
            p3df = p3d(:, goodInds);
            p2df = p2d(:, goodInds);
            
            [fG RG tG] = GPnP_f_GN(p3df, p2df);
            
            [err, projerr] = calcProjErr(RFin, tFin, p3df, fFin, 0, p2df); 
            [errG, projerrG] = calcProjErr(RG, tG, p3df, fG, 0, p2df); 
            filePrefix = [activeFolder int2str(i)];
            saveErr(filePrefix, err);    
            saveDistEstRes(filePrefix, RFin, tFin, fFin, 0);
            
            filePrefix = [activeFolder '_g_' int2str(i)];
            saveErr(filePrefix, errG);    
            saveDistEstRes(filePrefix, RG, tG, fG, 0);
            
            truefocal = str2num(char(folderSuffs{fInd})) / 35 * 3491
            fFin
            fG
            img1 = imread(imgPath);
            img2 = imresize(img1, 0.5);
            
            K(1,1) = fFin;
            K(2,2) = fFin;
            projs = K*(RFin*showPts3d+tFin*ones(1, size(showPts3d, 2)));
            projs(1, :) = projs(1, :) ./ projs(3, :);
            projs(2, :) = projs(2, :) ./ projs(3, :);
            projs(3, :) = projs(3, :) ./ projs(3, :);
            
            for ptInd = 1:size(projs, 2)
                cx = floor(projs(1, ptInd)+0.5);
                cy = floor(projs(2, ptInd)+0.5);
                img2(cy-6:cy+6, cx-6:cx+6, 1) = 255;
            end
            imwrite(img2, [activeFolder 'demo' int2str(i) '_' int2str(length(goodInds)) '.png']); 
            fileShowOut = fopen([activeFolder int2str(i) '_show.txt'], 'w');
            fprintf(fileShowOut, '%f ', projs(1:2, :));
            fclose(fileShowOut);
            
            if (size(RG, 1) > 0)
                K(1,1) = fG;
                K(2,2) = fG;
                projs = K*(RG*showPts3d+tG*ones(1, size(showPts3d, 2)));
                projs(1, :) = projs(1, :) ./ projs(3, :);
                projs(2, :) = projs(2, :) ./ projs(3, :);
                projs(3, :) = projs(3, :) ./ projs(3, :);                
                for ptInd = 1:size(projs, 2)
                    cx = floor(projs(1, ptInd)+0.5);
                    cy = floor(projs(2, ptInd)+0.5);
                    img2(cy-6:cy+6, cx-6:cx+6, 3) = 255;
                end
                imwrite(img2, [activeFolder 'demo' int2str(i) '_g_' int2str(length(goodInds)) '.png']);            
            end
        end
    end
end