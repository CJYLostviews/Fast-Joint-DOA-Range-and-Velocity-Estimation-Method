close all;
clear all;
clc;

frameIdx=2;

DataFile.master='master_0000_data.bin';
DataFile.masterIdxFile='master_0000_idx.bin';
DataFile.slave1= 'slave1_0000_data.bin';
DataFile.slave1IdxFile='slave1_0000_idx.bin';
DataFile.slave2='slave2_0000_data.bin';
DataFile.slave2IdxFile='slave2_0000_idx.bin';
DataFile.slave3='slave3_0000_data.bin';
DataFile.slave3IdxFile='slave3_0000_idx.bin';
DataFile.dataFolderName='G:\Tensor_for_FMCW_MIMO_Ti\Paper_Indoor\adc_data_0909_02\';
fileName=DataFile;
load('G:\Tensor_for_FMCW_MIMO_Ti\Paper_Indoor\Beamforming_Tensor\calibrateResults_high_0904.mat');

adcSampleRate = 8e6;
numADCSample = 2.560000e+02;
chirpRampTime=numADCSample/adcSampleRate;
numSamplePerChirp   = round(chirpRampTime*adcSampleRate);
nchirp_loops = 64;
numChirpsInLoop = 12;
numChirpsPerFrame = nchirp_loops*numChirpsInLoop;
numChirpPerLoop = numChirpsPerFrame/nchirp_loops;
numLoops = nchirp_loops;     
numRXPerDevice = 4;
numDevices=1;
Data=read_ADC_bin_TDA2_separateFiles(fileName,frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops,numRXPerDevice,numDevices);
RxForMIMOProcess=[13 14 15 16 1 2 3 4 9 10 11 12 5 6 7 8];
TxForMIMOProcess=[12  11  10   9   8   7   6   5   4   3   2   1];
TxToEnable=TxForMIMOProcess;
numTX = length(TxToEnable);
TX_ref = TxToEnable(1);
RangeMat=calibResult.RangeMat;
PeakValMat=calibResult.PeakValMat;
Sampling_Rate_sps=params.Sampling_Rate_sps;
outData=Data;
chirpSlope = 78.986000061035156e12;
Slope_calib = 78986000000000; 
fs_calib = 8000000; 
calibrationInterp = 5;
for iTX = 1: numTX
    TXind = TxToEnable(iTX);   
    freq_calib = (RangeMat(TXind,:)-RangeMat(TX_ref,1))*fs_calib/Sampling_Rate_sps *chirpSlope/Slope_calib;       
    freq_calib = 2*pi*(freq_calib)/(numSamplePerChirp * calibrationInterp);
    correction_vec = (exp(1i*((0:numSamplePerChirp-1)'*freq_calib))');
    freq_correction_mat = repmat(correction_vec, 1, 1, nchirp_loops);
    freq_correction_mat = permute(freq_correction_mat, [2 3 1]);
    outData1TX = Data(:,:,:,iTX).*freq_correction_mat;
    phase_calib = PeakValMat(TX_ref,1)./PeakValMat(TXind,:);
    phase_correction_mat = repmat(phase_calib.', 1,numSamplePerChirp, nchirp_loops);
    phase_correction_mat = permute(phase_correction_mat, [2 3 1]);
    outData(:,:,:,iTX) = outData1TX.*phase_correction_mat;
end
Data=outData(:,:,RxForMIMOProcess,:);
TI_Cascade_TX_position_azi = [11 10 9 32 28 24 20 16 12 8 4 0 ];
TI_Cascade_TX_position_ele = [6 4 1 0 0 0 0 0 0 0 0 0];
TI_Cascade_RX_position_azi = [ 11:14 50:53 46:49 0:3  ];
TI_Cascade_RX_position_ele = zeros(1,16);
[IdTxForMIMOProcess ia ib] = intersect(TxForMIMOProcess, TxToEnable,'stable' );
D_TX = TI_Cascade_TX_position_azi(TxToEnable(ib)); 
D_TX_ele = TI_Cascade_TX_position_ele(TxToEnable(ib));
D_RX = TI_Cascade_RX_position_azi(RxForMIMOProcess);
D_RX_ele = TI_Cascade_RX_position_ele(RxForMIMOProcess);
RX_id_tot = [];
RX_id_tot_ele = [];
 for ii = 1:length(D_TX)
    RX_id_new = D_RX + sum(D_TX(ii));
    RX_id_tot = [RX_id_tot RX_id_new];
    RX_id_new_ele = D_RX_ele + D_TX_ele(ii);
    RX_id_tot_ele = [RX_id_tot_ele RX_id_new_ele];
end
D(:,1) = RX_id_tot;
D(:,2) = RX_id_tot_ele;
ind = find(D(:,2)==0);
[val ID_unique] = unique(D(ind,1));
antenna_azimuthonly = ind(ID_unique); 

Data=reshape(Data,size(Data,1), size(Data,2), size(Data,3)*size(Data,4));
Data=Data(:,:,antenna_azimuthonly);
DATA = permute(Data,[3 1 2]);
DATA = reshape(DATA,size(DATA,1),size(DATA,2)*size(DATA,3));

startFreqConst	=	77e9;
chirpSlope	=	Slope_calib;
adcSampleTime	=	1/adcSampleRate;
Tpr	=	numADCSample*adcSampleTime;	%	Pulse Duration
adcStartTimeConst	=	6e-06;	
carrierFrequency	=	startFreqConst+(adcStartTimeConst+Tpr/2)*chirpSlope;
speedOfLight	=	3e8;
lamda	=	speedOfLight/carrierFrequency;

j	=	sqrt(-1);
M	=	size(DATA,1);
K	=	2;

DivisionX	=	mat2cell(DATA,M,[numADCSample*ones(1,nchirp_loops)]);

index	=	1;
XX	=	zeros(M*nchirp_loops,numADCSample);
for	m	=	1:M
	for	s	=	1:nchirp_loops
		XX(index,:)	=	DivisionX{s}(m,:);
		index	=	index+1;
	end
end
X	=	XX;

for i	=	1:M
	Y(i,:,:)	=	X((i-1)*nchirp_loops+1:i*nchirp_loops,:);
end

[AR,AT,SSS] = comfac(Y,K);

eDOA	=	PHase_DOA(AR,[0:M-1])
eVelocity	=	PHase_Velocity(AT,[1:nchirp_loops]*(Tpr/lamda))
eRange	=	PHase_Range(SSS,eVelocity,lamda,speedOfLight,chirpSlope,[1:numADCSample]*adcSampleTime)

Xr	=	reshape(X,nchirp_loops,M,numADCSample);
for k	=	1:numADCSample
	Xr1(:,:,k) = Xr(:,:,k).';
end
Xr2	=	reshape(Xr1,M,nchirp_loops*numADCSample);

Searching_DOA   =   -30:0.01:30;

R_DOA	=	Xr2*Xr2'/nchirp_loops*numADCSample;
[EV,D]	=	eig(R_DOA);
EVA	=	diag(D)';
[EVA,I]	=	sort(EVA);
EV	=	fliplr(EV(:,I));
En	=	EV(:,K+1:end);
for	i	=	1:length(Searching_DOA)
	a	=	exp(j*pi*(0:M-1).'*sind(Searching_DOA(i)));
	SP(i)	=	real(1./det(a'*En*En'*a));
end
SP	=	(SP-min(SP)*ones(1,size(Searching_DOA,2)))/(max(SP)-min(SP));	
[Peaks_Temp,DOA_Temp]	=	findpeaks(SP,Searching_DOA);
[~,Index]	=	sort(abs(Peaks_Temp),'descend');
eDOA_MUSIC	=	sort(DOA_Temp(Index(1:K)))

Searching_Velocity  =   -10:0.01:10;

Xv	=	zeros(nchirp_loops,M*numADCSample);
for	i	=	1:numADCSample
	Temp	=	zeros(nchirp_loops,M);
	for	ii	=	1:M
		Temp(:,ii)	=	X((ii-1)*nchirp_loops+1:ii*nchirp_loops,i);
	end
	Xv(:,(i-1)*M+1:i*M)	=	Temp;
end

R_v	=	Xv*Xv'/M*numADCSample;
[EV,D]	=	eig(R_v);
EVA	=	diag(D)';
[EVA,I]	=	sort(EVA);
EV	=	fliplr(EV(:,I));
En	=	EV(:,K+1:end);
for	i	=	1:length(Searching_Velocity)
	b	=	(exp(j*2*pi*(Tpr/lamda)*(Searching_Velocity(i))))*exp(j*2*pi*(0:nchirp_loops-1).'*(Tpr/lamda)*Searching_Velocity(i));
	SPv(i)	=	real(1./det(b'*En*En'*b));
end
SPv	=	(SPv-min(SPv)*ones(1,size(Searching_Velocity,2)))/(max(SPv)-min(SPv));	
[Peaks_Temp,Velocity_Temp]	=	findpeaks(SPv,Searching_Velocity);
[~,Index]	=	sort(abs(Peaks_Temp),'descend');
Velocity_MUSIC	=	sort(Velocity_Temp(Index(1:K)))

eRange_2DMUSIC_FFT	=	zeros(1,K);
for kk	=	1:K
	a	=	exp(j*pi*(0:M-1).'*sind(eDOA_MUSIC(kk)));
	b	=	(exp(j*2*pi*(Tpr/lamda)*(Velocity_MUSIC(kk))))*exp(j*2*pi*(0:nchirp_loops-1).'*(Tpr/lamda)*Velocity_MUSIC(kk));	
	cc	=	kron(a,b);
	YY	=	cc'*X;
	F	=	fft(YY,numADCSample);
	f_beat	=	(0:length(F)-1)'*adcSampleRate/length(F);
	[max_num,max_indx]	=	max(F);
	range	=	(speedOfLight/(2*chirpSlope))*(f_beat-(Velocity_MUSIC(kk)/lamda));
	eRange_2DMUSIC_FFT(kk)	=	(speedOfLight/(2*chirpSlope))*(f_beat(max_indx)-(Velocity_MUSIC(kk)/lamda));
end
eRange_2DMUSIC_FFT











