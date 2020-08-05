clear; close all; clc;

addpath plot_hht_modified/


%% load signal and load true labels
SUBJECT = 'S4';  %'S4' or 'X11'
% We are ignoring O3VR subject, because half of the trials are not used 
% accroding to true_labels_O3VR.txt. We want to avoid complications and 
% focus on the idea.

load(['../data/BCICompetitionIII/IIIb/', SUBJECT,'b_mat/', SUBJECT, 'b.mat'])

C3a_C3p = s(:,1)';
C4a_C4p = s(:,2)';
Fs = HDR.SampleRate;
dt = 1/Fs;
N = numel(C3a_C3p);
t = 0:dt:(N-1)*dt;

trigLabels = csvread(['../data/BCICompetitionIII/IIIb/true_labels_', SUBJECT,'.txt']);

% fill NaN labels (test trials) with correct labels
m = zeros(1,N);
m(HDR.TRIG) = trigLabels; 

%% SELECT left or right channel for further analysis
x = C3a_C3p; % left hemisphere
% x = C4a_C4p; % right hemisphere


%% see
m1 = m; m1(find(m==1)) = 1; m1(find(m~=1)) = 0; 
m2 = m; m2(find(m==2)) = 2; m2(find(m~=2)) = 0; 
plot(t,x,'k')
hold on
stem(t,m1*12,'r') % times 12 for visibility
stem(t,m2*12,'g')
hold off
legend({'sig','1','2'})
xlabel('Time [s]')
ylabel Amplitude
title('Signal and trial triggers indicated for each class')

%% store each trial separately for further processing
beg = HDR.TRIG;
nTrials = numel(beg);
nanTrialsIdx = [];
for iTrial = 1:nTrials
    label = trigLabels(iTrial);
    
    if(iTrial+1 <= nTrials)
        seg = x(beg(iTrial):beg(iTrial+1));
    else
        seg = x(beg(iTrial):end); % last trial
    end
    
    allSegments{iTrial} = seg;
    
    % mark trials containing nans
    if(any(isnan(seg)))
        nanTrialsIdx = [nanTrialsIdx, iTrial];
    end

end

%% EMD on selected trial
selectedTrial = 47; 

s = allSegments{selectedTrial};
N = numel(s);
t = 0:dt:(N-1)*dt; % individual trial axis


% uncomment if you want to apply bandpass filter 7-13 Hz before EMD
fa = 7; fb = 13;
[b,a] = butter(10,[fa fb]/(Fs/2));
x = filter(b,a,s);
% cwt(s,Fs) % check T-F if you like

[IMF, ~, ~] = emd1(s);

%% plot all IMFs from selected trial
nIMFs = size(IMF,1);
nPlots = nIMFs + 1;

h_all = figure;
for iPlot = 1:nPlots
    iIMF = iPlot-1;
    if iPlot == 1
        subplot(nPlots,1,iPlot)
        plot(t,s,'r');
        title(['signal (trial = ', num2str(selectedTrial),')'])
        grid on
    else
        subplot(nPlots,1,iPlot)
        plot(t, IMF(iIMF,:))
        title(['IMF no: ', num2str(iIMF)])
        grid on
    end % if iPlot == 1
end % iPlot

%% plot HHT from selected trial
figure
for SELECTED_IMF = 1:(size(IMF,1)-1) % we drop last residual because it is overpowering according to Huang
    I = IMF(SELECTED_IMF,:);
    % analytic signal
    z = hilbert(I); % WATCH OUT ! hilbert do not calculate hilbert transform but already an analytic signal !
    % instantaneous amplitude
    a = abs(z);
    % instantaeous energy
    e = a.^2;
    % instantaneuos phase
    th = angle(z);
    % instantaneous frequency
    f = 1/(2*pi)*diff(th)/dt;
    
    timeAxis = linspace(0,(N-2)*dt,N-1); % one sample lost during diff
    
    % color by component
%     plot(timeAxis, f, '.') 

    % color by energy
      scatter(timeAxis, f, 10, e(1:(N-1)), 'filled')
      set(gca,'color',[0 0 0])

    hold on
end
hold off
cmap = hot;
colormap(cmap)
colorbar
xlabel('Time [s]');
ylabel('Frequency');
xlim([0 timeAxis(end)])
%--- linear scale of Y
ylim([0 0.5*Fs]);
%--- log scale of Y
% set(gca, 'YScale','log')
% ylim([1 0.5*Fs])
%---
title(['Hilbert Spectrum (trial = ', num2str(selectedTrial),')'])


%% look how selected IMF of all signal correlates with class
nFirstTrials = nTrials; 

disp('Indicators ...')
% indicator signal for class 1
I_1 = [];
% indicator signal for class 2
I_2 = [];

for iTrial = 1:nFirstTrials
    label = trigLabels(iTrial);
    L = numel(allSegments{iTrial});
    
    % decide which class
    switch label
        case 1
            I_1 = [I_1 ones(1,L)];
            I_2 = [I_2 zeros(1,L)];
        case 2
            I_1 = [I_1 zeros(1,L)];
            I_2 = [I_2 ones(1,L)];

        otherwise
            error('unrecognized trial label!')
    end

end
figure
area(1:numel(I_1),I_1*100,'FaceColor','r','FaceAlpha',0.3)
hold on
area(1:numel(I_2),I_2*100,'FaceColor','g','FaceAlpha',0.3)
hold off
xlabel('Sample')
legend({'1','2'})
title('Timeline with indicated class')

%% calc energy envelope of selected IMF 

SELECTED_IMF = 1;
nFirstTrials = nTrials; % limit number of trials for faster exploration or take all (nTrials)

disp('Energy ...')
E = [];
F = [];
IMF_all = [];
count1 = 0;
count2 = 0;
for iTrial = 1:nFirstTrials
    disp([num2str(iTrial), ' / ', num2str(nFirstTrials)]);
    s = allSegments{iTrial};
    N = numel(s);

     if(any(isnan(s)))
         E = [E, nan(1,N)];
     else                 
         % uncomment if you want to apply bandpass filter 7-13 Hz before EMD
         fa = 7; fb = 13;
         [b,a] = butter(10,[fa fb]/(Fs/2));
         s = filter(b,a,s);
         
         [IMF, ~, ~] = emd1(s);
         imf = IMF(SELECTED_IMF,:);
%          imf = IMF(selectedIMF,:)+IMF(selectedIMF+1,:); % maybe take two neighbouring IMFs?
         
        % separate imfs from different class for further averaging
        label = trigLabels(iTrial);
        switch label
            case 1
                count1 = count1+1;
                class_1_IMFs_nFirstTrials{count1} = imf;
            case 2
                count2 = count2+1;
                class_2_IMFs_nFirstTrials{count2} = imf;          
            otherwise
                error('unrecognized trial label!')
        end

         E = [E, imf.^2];
         
         % extract also frequency to verify that you are looking on right oscillations       
         f = 1/(2*pi)*diff(angle(hilbert(imf)))/dt;
         F = [F, f];
     end
end

disp('OK')

%% plot IMF envelope overlapped with the class label indicators
% upper energy envelope
p = findpeaks2(E);
emax = spline([0 p numel(E)+1],[0 E(p) 0], 1:numel(E));

% plot to check
figure
plot(1:numel(emax),emax,'k')
hold on
area(1:numel(I_1),I_1*100,'FaceColor','r','FaceAlpha',0.3)
area(1:numel(I_2),I_2*100,'FaceColor','g','FaceAlpha',0.3)
hold off
legend({'env','1','2'})
title(['Energy envelope of IMF no', num2str(SELECTED_IMF), ' and class indicators (nFirstTrials)'])
xlabel('Sample')
%% check IMF frequency
figure
plot(1:numel(F),F)
ylabel('Frequency [Hz]')
xlabel Sample
title(['Instantaneuous frequency of IMF no',  num2str(SELECTED_IMF), ' (nFirstTrials)'])

%% correlations 

r1 = corr(E',I_1', 'Rows', 'complete') %  'Rows', 'complete' to omit NaNs
r2 = corr(E',I_2', 'Rows', 'complete')

%% correlation-like coefficients
% Negative correlation make it difficult to interpret, 
% so omit substracting mean when calculating correlation.
% Thus the lower q the more mu desynchronization is correlated to class label
q1 = corr_no_mean(E, I_1)
q2 = corr_no_mean(E, I_2)



%% divide for two classes
beg = HDR.TRIG;
nTrials = numel(beg);
count1 = 0;
count2 = 0;
nanTrialsIdx = [];
for iTrial = 1:nTrials
    label = trigLabels(iTrial);
    
    seg = allSegments{iTrial};

    % STORE SEPARATELY for further averaging
    % uncomment to remove trials with nan values (emd wont process them) 
    if(any(isnan(seg)))
        continue;
    end
    
    % decide to which cell store
    switch label
        case 1
            count1 = count1+1;
            class_1{count1} = seg;
        case 2
            count2 = count2+1;
            class_2{count2} = seg;
        otherwise
            error('unrecognized trial label!')
    end
end

%% look on averaged signal energy from all trials

numberOfSeconds = 8;
av_class_1 = zeros(1, numberOfSeconds*Fs);
nSeg = numel(class_1)
for iSeg = 1:nSeg
    av_class_1 = av_class_1 + class_1{iSeg}(1:(numberOfSeconds*Fs));
end
av_class_1_E = (av_class_1/nSeg).^2;

av_class_2 = zeros(1, numberOfSeconds*Fs);
nSeg = numel(class_2)
for iSeg = 1:nSeg
    av_class_2 = av_class_2 + class_2{iSeg}(1:(numberOfSeconds*Fs));
end
av_class_2_E = (av_class_2/nSeg).^2;


figure
t = 0:dt:(numberOfSeconds-dt);
plot(t, av_class_1_E,'r')
hold on
plot(t, av_class_2_E,'g')
hold off
legend({'av energy 1', 'av energy 2'})
xlabel('Time from trial onset [s]')
title('Averaged trials signal energy')
ylabel('Energy')

%% check signal average power after 7-13 Hz bandpass filtration

fa = 7; fb = 13;
[b,a] = butter(10,[fa fb]/(Fs/2));

numberOfSeconds = 8;
av_class_1 = zeros(1, numberOfSeconds*Fs);
nSeg = numel(class_1)
for iSeg = 1:nSeg
    x = class_1{iSeg}(1:(numberOfSeconds*Fs));
    y = filter(b,a,x);
    av_class_1 = av_class_1 + y;
end
av_class_1_E = (av_class_1/nSeg).^2;

av_class_2 = zeros(1, numberOfSeconds*Fs);
nSeg = numel(class_2)
for iSeg = 1:nSeg
    x = class_2{iSeg}(1:(numberOfSeconds*Fs));
    y = filter(b,a,x);
    av_class_2 = av_class_2 + y;
end
av_class_2_E = (av_class_2/nSeg).^2;

figure
t = 0:dt:(numberOfSeconds-dt);
plot(t, av_class_1_E,'r')
hold on
plot(t, av_class_2_E,'g')
hold off
legend({'av energy 1', 'av energy 2'})
xlabel('Time from trial onset [s]')
title('Averaged FILTERED trials signal energy')
ylabel('Energy')

%% look on averaged selected IMF energy from all trials
numberOfSeconds = 8;
av_IMF_class_1 = zeros(1, numberOfSeconds*Fs);

nSeg = numel(class_1_IMFs_nFirstTrials);
av_IMF_class_1 = zeros(1, numberOfSeconds*Fs);
for iSeg = 1:nSeg
    av_IMF_class_1 = av_IMF_class_1 + class_1_IMFs_nFirstTrials{iSeg}(1:(numberOfSeconds*Fs));
end
av_IMF_class_1_E = (av_IMF_class_1/nSeg).^2;

nSeg = numel(class_2_IMFs_nFirstTrials);
av_IMF_class_2 = zeros(1, numberOfSeconds*Fs);
for iSeg = 1:nSeg
    av_IMF_class_2 = av_IMF_class_2 + class_2_IMFs_nFirstTrials{iSeg}(1:(numberOfSeconds*Fs));
end
av_IMF_class_2_E = (av_IMF_class_2/nSeg).^2;


% upper energy envelopes of envelopes
p = findpeaks2(av_IMF_class_1_E);
emax_1 = spline([0 p numel(av_IMF_class_1_E)+1],[0 av_IMF_class_1_E(p) 0], 1:numel(av_IMF_class_1_E));
p = findpeaks2(av_IMF_class_2_E);
emax_2 = spline([0 p numel(av_IMF_class_2_E)+1],[0 av_IMF_class_2_E(p) 0], 1:numel(av_IMF_class_2_E));

figure
t = 0:dt:(numberOfSeconds-dt);
plot(t, av_IMF_class_1_E,'r')
hold on
plot(t, av_IMF_class_2_E,'g')
plot(t, emax_1, 'r','LineWidth',5);
plot(t, emax_2, 'g','LineWidth',5);
hold off
legend({'av energy 1', 'av energy 2'})
xlabel('Time from trial onset [s]')
ylabel('Energy')
title('Averaged trials selected IMF energy')



%% use CWT to examine signals T-F plane (animation)

% figure
% class_1_trials = find(trigLabels==1)';
% class_2_trials = find(trigLabels==2)';
% for ii = 1:numel(class_1_trials)
%     iTrialA = class_1_trials(ii);
%     iTrialB = class_2_trials(ii);
%     s1 = allSegments{iTrialA};
%     s2 = allSegments{iTrialB};
%     if(any(isnan(s1)) || any(isnan(s2)) )
%         continue;
%     else
%         subplot(2,1,1)
%         cwt(s1,Fs)
%         subplot(2,1,2)
%         cwt(s2,Fs)
%         pause(1)
%     end
% end

