function [AllFramesFTrealign_MCMC, MRS_struct] = MCMC_Realign(AllFramesFTrealign, MRS_struct)

ii = MRS_struct.ii;

% Hadamard reconstruction of subexperiments
orig.ON_ON = mean(AllFramesFTrealign(:,1:4:end),2);
orig.OFF_ON = mean(AllFramesFTrealign(:,2:4:end),2);
orig.ON_OFF = mean(AllFramesFTrealign(:,3:4:end),2);
orig.OFF_OFF = mean(AllFramesFTrealign(:,4:4:end),2);

orig.GABA =  orig.ON_ON + orig.OFF_ON - orig.ON_OFF - orig.OFF_OFF;
orig.GSH =   orig.ON_ON - orig.OFF_ON + orig.ON_OFF - orig.OFF_OFF;
orig.sum =   orig.ON_ON + orig.OFF_ON + orig.ON_OFF + orig.OFF_OFF;
orig.resid = orig.ON_ON - orig.OFF_ON - orig.ON_OFF + orig.OFF_OFF;

% % Set some starting conditions
% initPower = rms(orig.resid);
% 
f0 = mean(MRS_struct.out.f_results(ii,:));
f0 = f0 / MRS_struct.p.LarmorFreq / abs(MRS_struct.spec.freq(ii,1)-MRS_struct.spec.freq(ii,2));
% % sigma = std(MRS_struct.out.f_results(ii,:));
% % pd = makedist('Normal',mu,sigma);
% % x = min(MRS_struct.out.f_results(ii,:)):0.01:max(MRS_struct.out.f_results(ii,:));
% % prior_f = pdf(pd,x);
% foo = fitdist(MRS_struct.out.ph_results(ii,:)', 'normal');
% prior_f = random(foo, 1000,1);
 
ph0 = mean(MRS_struct.out.ph_results(ii,:));
% % sigma = std(MRS_struct.out.ph_results(ii,:));
% % pd = makedist('Normal',mu,sigma);
% % x = min(MRS_struct.out.ph_results(ii,:)):0.01:max(MRS_struct.out.ph_results(ii,:));
% % prior_ph = pdf(pd,x);
% foo = fitdist(MRS_struct.out.ph_results(ii,:)', 'normal');
% prior_ph = random(foo, 1000,1);

% jj = 1;
% MCMCpower = Inf;
% MCMCpower = zeros(size(prior_ph));
% nPts = 1:length(orig.resid);

lb = find(MRS_struct.spec.freq(ii,:) >= 4.68-0.5);
ub = find(MRS_struct.spec.freq(ii,:) <= 4.68+0.5);
waterRange = intersect(lb,ub);
lb = find(MRS_struct.spec.freq(ii,:) >= 3.185-0.05+0.02);
ub = find(MRS_struct.spec.freq(ii,:) <= 3.185+0.05+0.02);
ChoRange = intersect(lb,ub);
lb = find(MRS_struct.spec.freq(ii,:) >= 4.25);
ub = find(MRS_struct.spec.freq(ii,:) <= 0.5);
metabRange = intersect(lb,ub);

% X0 = [f0*ones(1,4) ph0*ones(1,4)];
X0 = f0*ones(1,4);
lb = min(MRS_struct.out.f_results(ii,:));
lb = lb / MRS_struct.p.LarmorFreq / abs(MRS_struct.spec.freq(ii,1)-MRS_struct.spec.freq(ii,2));
lb = lb*ones(1,4);
ub = max(MRS_struct.out.f_results(ii,:));
ub = ub / MRS_struct.p.LarmorFreq / abs(MRS_struct.spec.freq(ii,1)-MRS_struct.spec.freq(ii,2));
ub = ub*ones(1,4);

callFunc = @(c) RMS_Power1(orig.ON_ON, orig.OFF_ON, orig.ON_OFF, orig.OFF_OFF, metabRange, c);
% parsFit = fminsearch(callFunc, X0, optimset('TolX',1e-8,'MaxIter',1e5,'MaxFunEvals',1e4));
parsFit = fmincon(callFunc,X0,[],[],[],[],lb,ub,[],optimset('TolX',1e-8,'MaxIter',1e5,'MaxFunEvals',1e4));

% MCMC.ON_ON = circshift(orig.ON_ON * exp(1i*pi/180*parsFit(5)), round(parsFit(1)));
% MCMC.OFF_ON = circshift(orig.OFF_ON * exp(1i*pi/180*parsFit(5)), round(parsFit(2)));
% MCMC.ON_OFF = circshift(orig.ON_OFF * exp(1i*pi/180*parsFit(6)), round(parsFit(3)));
% MCMC.OFF_OFF = circshift(orig.OFF_OFF * exp(1i*pi/180*parsFit(7)), round(parsFit(4)));

MCMC.ON_ON = circshift(orig.ON_ON, round(parsFit(1)));
MCMC.OFF_ON = circshift(orig.OFF_ON, round(parsFit(2)));
MCMC.ON_OFF = circshift(orig.ON_OFF, round(parsFit(3)));
MCMC.OFF_OFF = circshift(orig.OFF_OFF, round(parsFit(4)));

X0 = ph0*ones(1,4);
lb = min(MRS_struct.out.ph_results(ii,:));
lb = lb*ones(1,4);
ub = max(MRS_struct.out.ph_results(ii,:));
ub = ub*ones(1,4);

% MCMC.ON_ON2(:,1) = real(MCMC.ON_ON);
% MCMC.ON_ON2(:,2) = imag(MCMC.ON_ON);
% MCMC.OFF_ON2(:,1) = real(MCMC.OFF_ON);
% MCMC.OFF_ON2(:,2) = imag(MCMC.OFF_ON);
% MCMC.ON_OFF2(:,1) = real(MCMC.ON_OFF);
% MCMC.ON_OFF2(:,2) = imag(MCMC.ON_OFF);
% MCMC.OFF_OFF2(:,1) = real(MCMC.OFF_OFF);
% MCMC.OFF_OFF2(:,2) = imag(MCMC.OFF_OFF);

callFunc = @(c) RMS_Power2(real(MCMC.ON_ON), real(MCMC.OFF_ON), real(MCMC.ON_OFF), real(MCMC.OFF_OFF), metabRange, c);
% parsFit = fminsearch(callFunc, X0, optimset('TolX',1e-8,'MaxIter',1e5,'MaxFunEvals',1e4));
parsFit = fmincon(callFunc,X0,[],[],[],[],lb,ub,[],optimset('TolX',1e-8,'MaxIter',1e5,'MaxFunEvals',1e4));

MCMC.ON_ON = orig.ON_ON * exp(1i*pi/180*parsFit(1));
MCMC.OFF_ON = orig.OFF_ON * exp(1i*pi/180*parsFit(2));
MCMC.ON_OFF = orig.ON_OFF * exp(1i*pi/180*parsFit(3));
MCMC.OFF_OFF = orig.OFF_OFF * exp(1i*pi/180*parsFit(4));


MCMC.GABA  = MCMC.ON_ON + MCMC.OFF_ON - MCMC.ON_OFF - MCMC.OFF_OFF;
MCMC.GSH   = MCMC.ON_ON - MCMC.OFF_ON + MCMC.ON_OFF - MCMC.OFF_OFF;
MCMC.sum   = MCMC.ON_ON + MCMC.OFF_ON + MCMC.ON_OFF + MCMC.OFF_OFF;
MCMC.resid = MCMC.ON_ON - MCMC.OFF_ON - MCMC.ON_OFF + MCMC.OFF_OFF;


function power = RMS_Power1(A,B,C,D,freqrange,c)

freq1 = c(1); freq2 = c(2); freq3 = c(3); freq4 = c(4);

A = circshift(A, round(freq1));
B = circshift(B, round(freq2));
C = circshift(C, round(freq3));
D = circshift(D, round(freq4));

power = rms(real(A - B - C + D));
% GABA = A + B - C - D;
% power = rms(GABA(freqrange));

function power = RMS_Power2(A,B,C,D,freqrange,c)

% freq1 = c(1); freq2 = c(2); freq3 = c(3); freq4 = c(4);
% phi1 = c(5); phi2 = c(6); phi3 = c(7); phi4 = c(8);
phi1 = c(1); phi2 = c(2); phi3 = c(3); phi4 = c(4);

% A = circshift(A, round(freq1));
% B = circshift(B, round(freq2));
% C = circshift(C, round(freq3));
% D = circshift(D, round(freq4));

% A2 = A(:,1)+1i*A(:,2);
% B2 = B(:,1)+1i*B(:,2);
% C2 = C(:,1)+1i*C(:,2);
% D2 = D(:,1)+1i*D(:,2);

A = A * exp(1i*pi/180*phi1);
B = B * exp(1i*pi/180*phi2);
C = C * exp(1i*pi/180*phi3);
D = D * exp(1i*pi/180*phi4);

% power = rms(A - B - C + D);
resid = A - B - C + D;
power = rms(real(resid(freqrange)));






