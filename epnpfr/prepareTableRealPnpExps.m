function prepareTableRealPnpExps
    folderSuffs = {'18', '24', '35', '50', '75', '105'};
    errTab = [];
    fTab = [];
    for fInd = 1:length(folderSuffs)        
       truefocal = str2num(char(folderSuffs{fInd})) / 35 * 3491
       
       activeFolder = ['sfmData\' char(folderSuffs{fInd}) '\'];       
       errRow = [];
       errRowG = [];
       fRow = [truefocal];
       for i = 0:2
           filePrefix = [activeFolder int2str(i)];
           err = loadErr(filePrefix);                            
           [Re, te, fe, ke] = loadDistEstRes(filePrefix);
           filePrefixG = [activeFolder '_g_' int2str(i)];                      
           errG = loadErr(filePrefixG);
           errRow = [errRow err errG];
           if (errG > 0)
               [Rg, tg, fg, kg] = loadDistEstRes(filePrefixG);               
               fRow = [fRow fe fg];
           else
               fRow = [fRow fe -1];
           end
           
       end
       errTab = [errTab; errRow];
       fTab = [fTab; fRow];
    end
    errTab
    fTab
    fileOut = fopen('tab.txt', 'w');
    for i = 1:size(errTab, 1)        
        for j = 1:3
            fprintf(fileOut, '%.2f & ', fTab(i, 1));
            fprintf(fileOut, ' & %.2f (%.2f) & ', errTab(i, 2*(j-1)+1), fTab(i, 1 + 2*(j-1)+1));
            fprintf(fileOut, ' %.2f (%.2f) \\\\', errTab(i, 2*(j-1)+2), fTab(i, 1 + 2*(j-1)+2));
            fprintf(fileOut, '\n');
        end        
        
    end       
    fclose(fileOut);
end