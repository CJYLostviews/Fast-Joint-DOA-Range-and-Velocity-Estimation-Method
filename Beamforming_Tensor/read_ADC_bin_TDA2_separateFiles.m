function [radar_data_Rxchain] = read_ADC_bin_TDA2_separateFiles(fileNameCascade,frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops,numRXPerDevice,numDevices)
dataFolder =fileNameCascade.dataFolderName;
fileFullPath_master = fullfile(dataFolder,fileNameCascade.master);
fileFullPath_slave1 = fullfile(dataFolder,fileNameCascade.slave1);
fileFullPath_slave2 = fullfile(dataFolder,fileNameCascade.slave2);
fileFullPath_slave3 = fullfile(dataFolder,fileNameCascade.slave3);
[radar_data_Rxchain_master] = readBinFile(fileFullPath_master, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
[radar_data_Rxchain_slave1] = readBinFile(fileFullPath_slave1, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
[radar_data_Rxchain_slave2] = readBinFile(fileFullPath_slave2, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
[radar_data_Rxchain_slave3] = readBinFile(fileFullPath_slave3, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
radar_data_Rxchain(:,:,1:4,:) = radar_data_Rxchain_master;%堆叠，按这个顺序在接收阵元上排列，4+4+4+4
radar_data_Rxchain(:,:,5:8,:) = radar_data_Rxchain_slave1;
radar_data_Rxchain(:,:,9:12,:) = radar_data_Rxchain_slave2;
radar_data_Rxchain(:,:,13:16,:) = radar_data_Rxchain_slave3;       
end