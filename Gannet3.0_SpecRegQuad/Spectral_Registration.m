function [AllFramesFTrealign, MRS_struct] = Spectral_Registration(MRS_struct, OnWhat, Dual)
%Spectral Registration is a time-domain frequency-and-phse correction as
%per Near et al. 2014 [under review].
% OnWhat = 0 for spectro data OnWhat=1 for water data




%MRS_struct.p.parsFit=[]; % MM
ii=MRS_struct.ii; % MM

%Dual-channel option only applies registration separately to ONs and OFFs
SpecRegLoop=0;
if(nargin==3)
    if(Dual==1)
        %We want to run this code twice, once for ONs, once for OFFs.
        SpecRegLoop=1;
    end
    if(Dual==2)
        %We want to run this code four times - once for A, once for B, ...
        % once for C, once for D.
        SpecRegLoop=3;
        MRS_struct.fids.ON_OFF = repmat([1 0 0 0], 1, (size(MRS_struct.fids.data,2)/4)); 
        MRS_struct.fids.ON_OFF = circshift(MRS_struct.fids.ON_OFF, SpecRegLoop, 2); 

    end

end

c = 1;
while(SpecRegLoop>(-1))
    
    
    
    if(OnWhat) %Read water data
        %First, take the complex data and turn it into a real matrix
        flatdata(:,1,:)=real(MRS_struct.fids.data_water);
        flatdata(:,2,:)=imag(MRS_struct.fids.data_water);
        
        
    else % read spectro data
        
        if(nargin==3)
            if(Dual==1)
                %This code runs twice, first for ONs, second for OFFs.
                %SpecRegLoop;
                size(real(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF==SpecRegLoop))));
                
                flatdata(:,1,:)=real(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF==SpecRegLoop)));
                flatdata(:,2,:)=imag(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF==SpecRegLoop)));
            end
            if(Dual==2) % ADH 
                %This code runs 4 times - once for each of A B C D.
                %SpecRegLoop;
                
                flatdata(:,1,:)=real(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF ==1)));
                flatdata(:,2,:)=imag(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF ==1)));

            end

        else
            %First, take the complex data and turn it into a real matrix
            flatdata(:,1,:)=real(MRS_struct.fids.data);
            flatdata(:,2,:)=imag(MRS_struct.fids.data);
        end
    end
    
    %Correct to a point 10% into the file (seems better that the actual
    %beginning) - this might not actually be a great idea - might be better
    %to use a window as Jamie does but something to come back to
    AlignRow=ceil(size(flatdata,3)/10); 
    MRS_struct.fids.flattarget=squeeze(flatdata(:,:,AlignRow));
    
    
    %Time domain Frequency and Phase Correction
    %Preliminary to fitting:
    parsGuess=[0 0]; %initial freq and phase guess
    parsFit = zeros([size(flatdata,3) 2]);
    input.dwelltime=1/MRS_struct.p.sw;
    time=((0:1:(MRS_struct.p.npoints-1)).'/MRS_struct.p.sw);
    %Fitting to determine frequency and phase corrections.
    reverseStr = ''; % MM
    for corrloop=1:size(flatdata,3)
        
        %corrloop
        target=MRS_struct.fids.flattarget(:);
        transient=squeeze(flatdata(:,:,corrloop));
        input.data=transient(:);
        parsFit(corrloop,:)=nlinfit(input,target,@FreqPhaseShiftNest,parsGuess);
        parsGuess = parsFit(corrloop,:); %Carry parameters from point to point
        
        % MM
        msg = sprintf('Spectral Registration - Fitting transient: %d', corrloop);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    % MM
    %MRS_struct.p.parsFit=[MRS_struct.p.parsFit parsFit];
    MRS_struct.p.parsFit(:,:,c,ii)=parsFit;
    f_results(:,:,c,ii) = -MRS_struct.p.parsFit(:,1,c,ii)'; % freq estimates (Hz)
    ph_results(:,:,c,ii) = MRS_struct.p.parsFit(:,2,c,ii)'; % phase estimates (deg)
    
    if(OnWhat)
        %Applyng frequency and phase corrections.
        for corrloop=1:size(flatdata,3)
            MRS_struct.fids.data_water(:,corrloop)=MRS_struct.fids.data_water(:,corrloop).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
        end
        FullData = MRS_struct.fids.data_water;
        FullData = FullData.* repmat((exp(-(time)*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data_water,2)]);
        AllFramesFTrealign=fftshift(fft(FullData,MRS_struct.p.ZeroFillTo,1),1);
        
    else
                
        %Applying frequency and phase corrections.
        for corrloop=1:size(flatdata,3)
            
            if(nargin==3)
                if(Dual==1)
                    %Need to get the slot right to put data back into
                    averages_per_dynamic=find(MRS_struct.fids.ON_OFF~=(MRS_struct.fids.ON_OFF(1)),1)-1;
                    dyn=floor((corrloop-1)/averages_per_dynamic); %number of cycles in
                    ind=mod((corrloop-1),averages_per_dynamic)+1; %number in current cycle
                    
                    if(SpecRegLoop==1)
                        if(MRS_struct.fids.ON_OFF(1)==1)
                            corrloop_d = dyn*averages_per_dynamic*2+ind;
                        else
                            corrloop_d = dyn*averages_per_dynamic*2+averages_per_dynamic+ind;
                        end
                    else
                        if(MRS_struct.fids.ON_OFF(1)==1)
                            corrloop_d = dyn*averages_per_dynamic*2+averages_per_dynamic+ind;
                        else
                            corrloop_d = dyn*averages_per_dynamic*2+ind;
                        end
                    end
                    
                    MRS_struct.fids.data_align(:,corrloop_d)=MRS_struct.fids.data(:,corrloop_d).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
                end
                
                
                 if(Dual==2) %ADH
                    %Need to get the slot right to put data back into

                    
                    corrloop_d=(corrloop-1)*4+SpecRegLoop+1;
                    
                    MRS_struct.fids.data_align(:,corrloop_d)=MRS_struct.fids.data(:,corrloop_d).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
           
                 end

                
                
                
                CorrPars(corrloop_d,:)=parsFit(corrloop,:);
            else
                corrloop_d=corrloop;
                MRS_struct.fids.data_align(:,corrloop_d)=MRS_struct.fids.data(:,corrloop_d).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
                CorrPars(corrloop_d,:)=parsFit(corrloop,:);
            end
            
        end
        
        if(SpecRegLoop==0)
            FullData = MRS_struct.fids.data_align;
            FullData = FullData.* repmat( (exp(-(time)*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data,2)]);
            AllFramesFTrealign=fftshift(fft(FullData,MRS_struct.p.ZeroFillTo,1),1);
            
            
            %In FD, move Cr to 3.02 and get phase 'right' as opposed to 'consistent'.
            ChoCrFitLimLow=2.6;
            ChoCrFitLimHigh=3.6;
            %Still need ranges for Creatine align plot
            z=abs(MRS_struct.spec.freq(ii,:)-ChoCrFitLimHigh);
            cclb=find(min(z)==z);
            z=abs(MRS_struct.spec.freq(ii,:)-ChoCrFitLimLow);
            ccub=find(min(z)==z);
            freqrange=MRS_struct.spec.freq(ii,cclb:ccub);
            if isempty(freqrange) % MM (for simulated data)
                [ccub,cclb] = deal(cclb,ccub);
                freqrange=MRS_struct.spec.freq(ii,cclb:ccub);
            end
            %Do some detective work to figure out the initial parameters
            ChoCrMeanSpec = mean(AllFramesFTrealign(cclb:ccub,:),2);
            Baseline_offset=real(ChoCrMeanSpec(1)+ChoCrMeanSpec(end))/2;
            Width_estimate=0.05;%ppm
            Area_estimate=(max(real(ChoCrMeanSpec))-min(real(ChoCrMeanSpec)))*Width_estimate*4;
            ChoCr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1].*[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 1];
            
            
            if(nargin==3)
                if(Dual==1)
                    %This bit is silly - we don't want to do OFF-to-ON based on the CR signal
                    ChoCrMeanSpecON = mean(AllFramesFTrealign(cclb:ccub,(MRS_struct.fids.ON_OFF==1)),2);
                    ChoCrMeanSpecOFF = mean(AllFramesFTrealign(cclb:ccub,(MRS_struct.fids.ON_OFF==0)),2);
                    ChoCrMeanSpecFitON = FitChoCr(freqrange, ChoCrMeanSpecON, ChoCr_initx,MRS_struct.p.LarmorFreq);
                    ChoCrMeanSpecFitOFF = FitChoCr(freqrange, ChoCrMeanSpecOFF, ChoCr_initx,MRS_struct.p.LarmorFreq);
                    AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1))=AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1))*exp(1i*pi/180*(ChoCrMeanSpecFitON(4)));%phase
                    AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0))=AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0))*exp(1i*pi/180*(ChoCrMeanSpecFitOFF(4)));%phase
                    
                    ChoCrFreqShiftON = ChoCrMeanSpecFitON(3);
                    ChoCrFreqShiftON = ChoCrFreqShiftON - 3.02*MRS_struct.p.LarmorFreq;
                    ChoCrFreqShiftON = ChoCrFreqShiftON ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(ii,2) - MRS_struct.spec.freq(ii,1) ));
                    ChoCrFreqShift_pointsON = round(ChoCrFreqShiftON);
                    AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1))=circshift(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1)), [-ChoCrFreqShift_pointsON 0]);%freq
                    ChoCrFreqShiftOFF = ChoCrMeanSpecFitOFF(3);
                    ChoCrFreqShiftOFF = ChoCrFreqShiftOFF - 3.02*MRS_struct.p.LarmorFreq;
                    ChoCrFreqShiftOFF = ChoCrFreqShiftOFF ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(ii,2) - MRS_struct.spec.freq(ii,1) ));
                    ChoCrFreqShift_pointsOFF = round(ChoCrFreqShiftOFF);
                    AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0))=circshift(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0)), [-ChoCrFreqShift_pointsOFF 0]);%freq
                    
                end
                
%                 MRS_struct.out.FreqStdevHz(MRS_struct.ii)=std(parsFit(:,1),1);
%                 MRS_struct.out.CrFWHMHz(MRS_struct.ii)=mean([ChoCrMeanSpecFitON(2) ChoCrMeanSpecFitOFF(2)]);
            if(Dual == 2)
%                     ChoCrMeanSpec_1 = mean(AllFramesFTrealign(cclb:ccub,1:4:320),2);
%                     ChoCrMeanSpec_2 = mean(AllFramesFTrealign(cclb:ccub,2:4:320),2);
%                     ChoCrMeanSpec_3 = mean(AllFramesFTrealign(cclb:ccub,3:4:320),2);
%                     ChoCrMeanSpec_4 = mean(AllFramesFTrealign(cclb:ccub,4:4:320),2);
%                     
%                     ChoCrMeanSpecFit_1 = FitChoCr(freqrange, ChoCrMeanSpec_1, ChoCr_initx,MRS_struct.p.LarmorFreq);
%                     ChoCrMeanSpecFit_2 = FitChoCr(freqrange, ChoCrMeanSpec_2, ChoCr_initx,MRS_struct.p.LarmorFreq);
%                     ChoCrMeanSpecFit_3 = FitChoCr(freqrange, ChoCrMeanSpec_3, ChoCr_initx,MRS_struct.p.LarmorFreq);
%                     ChoCrMeanSpecFit_4 = FitChoCr(freqrange, ChoCrMeanSpec_4, ChoCr_initx,MRS_struct.p.LarmorFreq);
%                    
%                     AllFramesFTrealign(:,1:4:320)=AllFramesFTrealign(:,1:4:320)*exp(1i*pi/180*(ChoCrMeanSpecFit_1(4)));%phase
%                     AllFramesFTrealign(:,2:4:320)=AllFramesFTrealign(:,2:4:320)*exp(1i*pi/180*(ChoCrMeanSpecFit_1(4)));%phase
%                     AllFramesFTrealign(:,3:4:320)=AllFramesFTrealign(:,3:4:320)*exp(1i*pi/180*(ChoCrMeanSpecFit_1(4)));%phase
%                     AllFramesFTrealign(:,4:4:320)=AllFramesFTrealign(:,4:4:320)*exp(1i*pi/180*(ChoCrMeanSpecFit_1(4)));%phase
%                   
%                     ChoCrFreqShift_1 = ChoCrMeanSpecFit_1(3);
%                     ChoCrFreqShift_1 = ChoCrFreqShift_1 - 3.02*MRS_struct.p.LarmorFreq;
%                     ChoCrFreqShift_1 = ChoCrFreqShift_1 ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(ii,2) - MRS_struct.spec.freq(ii,1) ));
%                     ChoCrFreqShift_points_1 = round(ChoCrFreqShift_1);
%                     AllFramesFTrealign(:,1:4:320)=circshift(AllFramesFTrealign(:,1:4:320), [-ChoCrFreqShift_points_1 0]);%freq
%                     
%                     ChoCrFreqShift_2 = ChoCrMeanSpecFit_2(3);
%                     ChoCrFreqShift_2 = ChoCrFreqShift_2 - 3.02*MRS_struct.p.LarmorFreq;
%                     ChoCrFreqShift_2 = ChoCrFreqShift_2 ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(ii,2) - MRS_struct.spec.freq(ii,1) ));
%                     ChoCrFreqShift_points_2 = round(ChoCrFreqShift_2);
%                     AllFramesFTrealign(:,2:4:320)=circshift(AllFramesFTrealign(:,2:4:320), [-ChoCrFreqShift_points_2 0]);%freq
%                     
%                     ChoCrFreqShift_3 = ChoCrMeanSpecFit_3(3);
%                     ChoCrFreqShift_3 = ChoCrFreqShift_3 - 3.02*MRS_struct.p.LarmorFreq;
%                     ChoCrFreqShift_3 = ChoCrFreqShift_3 ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(ii,2) - MRS_struct.spec.freq(ii,1) ));
%                     ChoCrFreqShift_points_3 = round(ChoCrFreqShift_3);
%                     AllFramesFTrealign(:,3:4:320)=circshift(AllFramesFTrealign(:,3:4:320), [-ChoCrFreqShift_points_3 0]);%freq
%                     
%                     ChoCrFreqShift_4 = ChoCrMeanSpecFit_4(3);
%                     ChoCrFreqShift_4 = ChoCrFreqShift_4 - 3.02*MRS_struct.p.LarmorFreq;
%                     ChoCrFreqShift_4 = ChoCrFreqShift_4 ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(ii,2) - MRS_struct.spec.freq(ii,1) ));
%                     ChoCrFreqShift_points_4 = round(ChoCrFreqShift_4);
%                     AllFramesFTrealign(:,4:4:320)=circshift(AllFramesFTrealign(:,4:4:320), [-ChoCrFreqShift_points_4 0]);%freq
                    
             else
                ChoCrMeanSpecFit = FitChoCr(freqrange, ChoCrMeanSpec, ChoCr_initx,MRS_struct.p.LarmorFreq);
                
                % MM
%                 MRS_struct.out.f_results(ii,:) = -(MRS_struct.p.parsFit(:,1,ii)' + -(ChoCrMeanSpecFit(3) - 3.02*MRS_struct.p.LarmorFreq)); % freq estimates (Hz)
%                 MRS_struct.out.ph_results(ii,:) = MRS_struct.p.parsFit(:,2,ii)' + ChoCrMeanSpecFit(4); % phase estimates (deg)
                
                MRS_struct.out.ChoCrMeanSpecFit = ChoCrMeanSpecFit./[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 1];
                AllFramesFTrealign=AllFramesFTrealign*exp(1i*pi/180*(ChoCrMeanSpecFit(4)));%phase
                ChoCrFreqShift = ChoCrMeanSpecFit(3);
                ChoCrFreqShift = ChoCrFreqShift - 3.02*MRS_struct.p.LarmorFreq;
                ChoCrFreqShift = ChoCrFreqShift ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(ii,2) - MRS_struct.spec.freq(ii,1) ));
                ChoCrFreqShift_points = round(ChoCrFreqShift);
                AllFramesFTrealign=circshift(AllFramesFTrealign, [-ChoCrFreqShift_points 0]);%freq
            end
            %figure(5)
            %plot(freqrange,ChoCrMeanSpec,freqrange,TwoLorentzModel(ChoCrMeanSpecFit./[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 1],freqrange))
            
            
            %Fit just the Cr in the aligned mean spectrum to get CrFWHMHz
            CrFitLimLow=2.6;
            CrFitLimHigh=3.11;
            %Still need ranges for Creatine align plot
            z=abs(MRS_struct.spec.freq(ii,:)-CrFitLimHigh);
            clb=find(min(z)==z);
            z=abs(MRS_struct.spec.freq(ii,:)-CrFitLimLow);
            cub=find(min(z)==z);
            freqrange=MRS_struct.spec.freq(ii,clb:cub);
            if isempty(freqrange) % MM (for simulated data)
                [cub,clb] = deal(clb,cub);
                freqrange=MRS_struct.spec.freq(ii,clb:cub);
            end
            Cr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 ].*[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 ];
            CrMeanSpec = mean(AllFramesFTrealign(clb:cub,:),2);
            CrMeanSpecFit = FitCr(freqrange, CrMeanSpec, Cr_initx);
            
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
            IsRejectOnOrOff=MRS_struct.fids.ON_OFF(MRS_struct.out.reject(:,MRS_struct.ii)==1);
            stepsize=find(MRS_struct.fids.ON_OFF~=(MRS_struct.fids.ON_OFF(1)),1)-1;
            for kk=1:size(MRS_struct.out.reject(:,MRS_struct.ii),1)
                %first find whether reject is ON or OFF
                if(MRS_struct.out.reject(kk,MRS_struct.ii)==1)
                    IsRejectOnOrOff=MRS_struct.fids.ON_OFF(kk);
                    if IsRejectOnOrOff==MRS_struct.fids.ON_OFF(1)
                        MRS_struct.out.reject(kk+stepsize,MRS_struct.ii)=1;
                    else
                        MRS_struct.out.reject(kk-stepsize,MRS_struct.ii)=1;
                    end
                end
            end
            
            
            
            
            
            
        end
        end
    SpecRegLoop
    SpecRegLoop=SpecRegLoop-1;
    c = c + 1;
 
end



end

MRS_struct.out.f_results(ii,:) = zeros(320,1);
MRS_struct.out.f_results(ii,1:4:end) = f_results(:,:,1,ii);
MRS_struct.out.f_results(ii,2:4:end) = f_results(:,:,2,ii);
MRS_struct.out.f_results(ii,3:4:end) = f_results(:,:,3,ii);
MRS_struct.out.f_results(ii,4:4:end) = f_results(:,:,4,ii);

MRS_struct.out.ph_results(ii,:) = zeros(320,1);
MRS_struct.out.ph_results(ii,1:4:end) = ph_results(:,:,1,ii);
MRS_struct.out.ph_results(ii,2:4:end) = ph_results(:,:,2,ii);
MRS_struct.out.ph_results(ii,3:4:end) = ph_results(:,:,3,ii);
MRS_struct.out.ph_results(ii,4:4:end) = ph_results(:,:,4,ii);

end



