function data = readData(fileName)

fileID = fopen(fileName,'r');
tmp = fscanf(fileID,'%g		%g',[2 Inf]);
fclose(fileID);

data = tmp(1,:) + 1i*tmp(2,:);