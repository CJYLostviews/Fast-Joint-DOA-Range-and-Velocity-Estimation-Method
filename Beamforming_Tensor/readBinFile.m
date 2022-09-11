function [adcData1Complex] = readBinFile(fileFullPath, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices)
Expected_Num_SamplesPerFrame = numSamplePerChirp*numChirpPerLoop*numLoops*numRXPerDevice*2;%乘以2是因为数据为虚实相间，如果纯实数这里不需要乘2
fp = fopen(fileFullPath, 'r');
fseek(fp,(frameIdx-1)*Expected_Num_SamplesPerFrame*2, 'bof');
adcData1 = fread(fp,Expected_Num_SamplesPerFrame,'uint16');
neg             = logical(bitget(adcData1, 16));
adcData1(neg)    = adcData1(neg) - 2^16;
adcData1 = adcData1(1:2:end) + sqrt(-1)*adcData1(2:2:end);%数据结构：实部，虚部，实部，虚部，依次排列，故需要合成在一起
adcData1Complex = reshape(adcData1, numRXPerDevice, numSamplePerChirp, numChirpPerLoop, numLoops);%重组数据，将二维数据重组
adcData1Complex = permute(adcData1Complex, [2 4 1 3]);%重新排列维度，重构后的维度为numSamplePerChirp*numLoops*numRXPerDevice*numChirpPerLoop
fclose(fp);
end