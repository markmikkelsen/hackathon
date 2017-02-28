function [AllFramesFTrealign, MRS_struct] = Spectral_Registration2(MRS_struct)
%Spectral Registration is a time-domain frequency-and-phse correction as
%per Near et al. 2014 [under review].
% OnWhat = 0 for spectro data OnWhat=1 for water data

ii = MRS_struct.ii; % MM
MRS_struct.fids.data_align = zeros(size(MRS_struct.fids.data));
CorrPars = zeros(size(MRS_struct.fids.data,2),2);
MRS_struct.out.f_results(ii,:) = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.out.ph_results(ii,:) = zeros(1,size(MRS_struct.fids.data,2));

nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts, 'MaxIter', 1e5);

%We want to run this code four times - once for A, once for B, ...
% once for C, once for D.
SpecRegLoop = 3;
SubspecToAlign = repmat([3 2 1 0], [1 size(MRS_struct.fids.data,2)/4]);
tLim = size(MRS_struct.fids.data,1)/2; % use first half of time-domain data

while(SpecRegLoop>(-1)) % -2
    % read spectro data
    if any(SpecRegLoop == 0:3)
        flatdata(:,1,:)=real(MRS_struct.fids.data(1:tLim,(SubspecToAlign==SpecRegLoop)));
        flatdata(:,2,:)=imag(MRS_struct.fids.data(1:tLim,(SubspecToAlign==SpecRegLoop)));
        
        % For HERMES subexperiments, correct to a point 10% into the file
        AlignRow=ceil(size(flatdata,3)/10);
        MRS_struct.fids.flattarget=squeeze(flatdata(:,:,AlignRow));
    else
        clear flatdata
        flatdata(:,1,:)=real(MRS_struct.fids.data_align(1:tLim,:));
        flatdata(:,2,:)=imag(MRS_struct.fids.data_align(1:tLim,:));
        
        % For fifth run through correct to first transient in OFF/OFF
        % experiment
        AlignRow=ceil(size(flatdata,3)/4)*3+1;
        MRS_struct.fids.flattarget=squeeze(flatdata(:,:,AlignRow));
    end
    
    %Time domain Frequency and Phase Correction
    %Preliminary to fitting:
    parsGuess=[0 0]; %initial freq and phase guess
    parsFit = zeros([size(flatdata,3) 2]);
    input.dwelltime=1/MRS_struct.p.sw;
    time=((0:1:(MRS_struct.p.npoints-1)).'/MRS_struct.p.sw);
    %Fitting to determine frequency and phase corrections.
    reverseStr = ''; % MM    
    for corrloop=1:size(flatdata,3)
        target=MRS_struct.fids.flattarget(:);
        transient=squeeze(flatdata(:,:,corrloop));
        input.data=transient(:);
        parsFit(corrloop,:)=nlinfit(input,target,@FreqPhaseShiftNest,parsGuess,nlinopts);
        parsGuess = parsFit(corrloop,:); %Carry parameters from point to point
        
        % MM
        msg = sprintf('\nSpectral Registration - Fitting transient: %d', corrloop);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    %Applying frequency and phase corrections.
    if any(SpecRegLoop == 0:3)
        corrloop_d = find(SubspecToAlign==SpecRegLoop);
    end
    
    for corrloop=1:size(flatdata,3)
        if any(SpecRegLoop == 0:3)
            %Need to get the slot right to put data back into
            MRS_struct.fids.data_align(:,corrloop_d(corrloop))=MRS_struct.fids.data(:,corrloop_d(corrloop)).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
        else
            corrloop_d=corrloop;
            MRS_struct.fids.data_align(:,corrloop_d)=MRS_struct.fids.data_align(:,corrloop_d).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
            CorrPars(corrloop_d,:)=parsFit(corrloop,:);
        end
    end
    
    if any(SpecRegLoop == 0:3)
        CorrPars(corrloop_d,:) = parsFit;
        MRS_struct.out.f_results(ii,corrloop_d) = -parsFit(:,1)';
        MRS_struct.out.ph_results(ii,corrloop_d) = parsFit(:,2)';
    end
    
    if(SpecRegLoop==0) % -1
        MRS_struct.out.SpecRegParsFit(:,:,ii) = CorrPars;

        FullData = MRS_struct.fids.data_align;
        FullData = FullData.* repmat( (exp(-(time)*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data,2)]);
        AllFramesFTrealign = fftshift(fft(FullData,MRS_struct.p.ZeroFillTo,1),1);
        
        %In FD, move Cr to 3.02 and get phase 'right' as opposed to 'consistent'.   
        freqrange = MRS_struct.spec.freq(ii,:) >= 2.925 & MRS_struct.spec.freq(ii,:) <= 3.125;
        
        FrameMax = max(real(AllFramesFTrealign(freqrange,:)),[],1);
        FrameMaxPos = zeros(1,size(AllFramesFTrealign,2));
        for jj = 1:size(AllFramesFTrealign,2)
            FrameMaxPos(jj) = find(real(AllFramesFTrealign(:,jj)) == FrameMax(jj));
        end
        CrFreqShift = MRS_struct.spec.freq(ii,FrameMaxPos);
        if ~strcmp(MRS_struct.gabafile{ii}((end-3):end),'.mat')
            CrFreqShift = CrFreqShift - 3.02;
        else
            CrFreqShift = CrFreqShift - 3.027;
        end
        CrFreqShift_pts = CrFreqShift / abs(MRS_struct.spec.freq(ii,1) - MRS_struct.spec.freq(ii,2));
        CrFreqShift_pts = round(CrFreqShift_pts);
        
        if strcmp(MRS_struct.gabafile{ii}((end-3):end),'.mat')
            CrFreqShift = -CrFreqShift;
            CrFreqShift_pts = -CrFreqShift_pts;
        end
        
        for jj = 1:size(AllFramesFTrealign,2)
            AllFramesFTrealign(:,jj) = circshift(AllFramesFTrealign(:,jj), CrFreqShift_pts(jj));
        end
        
        MRS_struct.out.CrFreqShift(ii,:) = CrFreqShift;
        MRS_struct.out.f_results(ii,:) = MRS_struct.out.f_results(ii,:) - CrFreqShift * MRS_struct.p.LarmorFreq;
        
        freq = MRS_struct.spec.freq(ii,freqrange);
%         for jj = 1:4
%             SubspecToAlign_ChoCr = SubspecToAlign==abs(jj-4);
            SubspecToAlign_ChoCr = 1:size(MRS_struct.fids.data,2);
            %Do some detective work to figure out the initial parameters
%             ChoCrMeanSpec = mean(AllFramesFTrealign(freqrange,SubspecToAlign_ChoCr),2);
%             Baseline_offset = real(ChoCrMeanSpec(1)+ChoCrMeanSpec(end))/2;
%             Area_estimate = (max(real(ChoCrMeanSpec))-min(real(ChoCrMeanSpec)))*Width_estimate*4;
%             ChoCr_initx = [Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1].*[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 1];

%             ChoCrMeanSpecFit = FitChoCr(freqrange, ChoCrMeanSpec, ChoCr_initx, MRS_struct.p.LarmorFreq);
            
            CrMeanSpec = mean(AllFramesFTrealign(freqrange,SubspecToAlign_ChoCr),2);
            Baseline_offset = real(CrMeanSpec(1)+CrMeanSpec(end))/2;
            Width_estimate = 0.05;%ppm
            Area_estimate = (max(real(CrMeanSpec))-min(real(CrMeanSpec)))*Width_estimate*4;
            Cr_initx = [Area_estimate Width_estimate 3.02 0 Baseline_offset 0 ].*[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1];
            
            CrMeanSpecFit = FitCr(freq, CrMeanSpec, Cr_initx);
            
            MRS_struct.out.CrMeanSpecFit(jj,:,ii) = CrMeanSpecFit./[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1];
%             AllFramesFTrealign(:,SubspecToAlign_ChoCr) = AllFramesFTrealign(:,SubspecToAlign_ChoCr)*exp(1i*pi/180*CrMeanSpecFit(4));%phase

%             MRS_struct.out.ph_results(ii,SubspecToAlign_ChoCr) = MRS_struct.out.ph_results(ii,SubspecToAlign_ChoCr) + CrMeanSpecFit(4); % phase estimates (deg)
%         end
%         AllFramesFTrealign = AllFramesFTrealign2;
                 
        %figure(5)
        %plot(freqrange,ChoCrMeanSpec,freqrange,TwoLorentzModel(ChoCrMeanSpecFit./[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 1],freqrange))
        
        %Fit just the Cr in the aligned mean spectrum to get CrFWHMHz
        freqrange = MRS_struct.spec.freq(ii,:) >= 2.6 & MRS_struct.spec.freq(ii,:) <= 3.11;
        freq = MRS_struct.spec.freq(ii,freqrange);
        
        Cr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 ].*[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 ];
        CrMeanSpec = mean(AllFramesFTrealign(freqrange,:),2);
        CrMeanSpecFit = FitCr(freq, CrMeanSpec, Cr_initx);
        
        %Some Output
        MRS_struct.out.FreqStdevHz(MRS_struct.ii)=std(parsFit(:,1),1);
        % GO 01/29/16: output frequency correction in Hz for every
        % average
        MRS_struct.out.FreqHz(MRS_struct.ii,:)=parsFit(:,1).';
        MRS_struct.out.CrFWHMHz(MRS_struct.ii)=CrMeanSpecFit(2);
        
        %Decide which rows to reject based on 3-sigma
        % Reject any point where the fit params - phase
        %  or freq are > 3stdev away from the mean
        MeanFrameParams = mean(CorrPars, 1);
        UpperLim = repmat(MeanFrameParams + 3*std(CorrPars,1), [size(AllFramesFTrealign,2) 1]);
        LowerLim = repmat(MeanFrameParams - 3*std(CorrPars,1), [size(AllFramesFTrealign,2) 1]);
        %but don't reject on linear, const baseline fit vals
        rejectframe = gt(CorrPars, UpperLim);
        rejectframe = rejectframe + lt(CorrPars, LowerLim);
        MRS_struct.out.reject(:,MRS_struct.ii) = max(rejectframe,[],2);
        %Balance up rejects
        IsRejectOnOrOff=SubspecToAlign(MRS_struct.out.reject(:,MRS_struct.ii)==1);
        stepsize=find(SubspecToAlign~=(SubspecToAlign(1)),1)-1;
        for kk=1:size(MRS_struct.out.reject(:,MRS_struct.ii),1)
            %first find whether reject is ON or OFF
            if(MRS_struct.out.reject(kk,MRS_struct.ii)==1)
                IsRejectOnOrOff=SubspecToAlign(kk);
                if IsRejectOnOrOff==SubspecToAlign(1)
                    MRS_struct.out.reject(kk+stepsize,MRS_struct.ii)=1;
                else
                    MRS_struct.out.reject(kk-stepsize,MRS_struct.ii)=1;
                end
            end
        end
        
    end
    
    SpecRegLoop=SpecRegLoop-1;
    
end

end



