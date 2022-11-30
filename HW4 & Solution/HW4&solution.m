load('BarrelSpikes.mat');
T=10;
sample_interval = .001;
t = sample_interval:sample_interval:T;
spike_inds1 = find(spike1);
spike_inds2 = find(spike2);
spike_times1 = t(spike_inds1);
spike_times2 = t(spike_inds2);
%Q1

t1 = (0:length(spike1)-1);
t2 = (0:length(spike2)-1);

subplot(2,1,1)
plot(t1,spike1)
title('spike1')

subplot(2,1,2)
plot(t2,spike2)
title('spike2')
xlabel('Time (s)')

figure()
subplot(2, 1, 1)
hold on
plot(spike_times1, ones(1, length(spike_times1)), '.k', 'MarkerSize', 10)
plot(spike_times2, 2*ones(1, length(spike_times2)), '.r', 'MarkerSize', 10)
hold off
xlabel('Time (seconds)')
ylim([0 3])
ylabel('Spike Train')
title('Rat Cortex Experiment')
legend({'Spike 1', 'Spike 2'})
set(gca, 'FontSize', 14)
subplot(2, 1, 2)
hold on
plot(spike_times1, ones(1, length(spike_times1)), '.k', 'MarkerSize', 10)
plot(spike_times2, 2*ones(1, length(spike_times2)), '.r', 'MarkerSize', 10)
hold off
xlabel('Time (seconds)')
ylim([0 3])
xlim([0 1])
ylabel('Spike Train')
title('Rat Cortex Experiment')
legend({'Spike 1', 'Spike 2'})
set(gca, 'FontSize', 14)

%Q2




    [auto_corr1, lags] = xcorr(spike1 - mean(spike1), 200, 'coeff');
    [auto_corr2, ~] = xcorr(spike2 - mean(spike2), 200, 'coeff');
    figure()
    subplot(2, 1, 1)
    plot(lags, auto_corr1, 'k', 'LineWidth', 2)
    xlabel('Lag (milliseconds)')
    ylabel('Autocorrelation')
    title('Spike 1')
    set(gca, 'FontSize', 14)
    subplot(2, 1, 2)
    plot(lags, auto_corr2, 'k', 'LineWidth', 2)
    xlabel('Lag (milliseconds)')
    ylabel('Autocorrelation')
    title('Spike 2')
    set(gca, 'FontSize', 14)

%Q3
f = 1000*(0:(5000))/5000;
test = conj(fft(spike1)).*fft(spike1);
length(test);
psd1b=conj(fft(spike1)).*fft(spike1);
psd1b(1:10)
psd1a(1:10)
psd1a=ffr(xcorr(spike1));
psd1a(1:10)
figure();
plot(f,psd1a(1:length(f)))
plot(f,psd1b(1:length(f)))




%Q4

    %TW = 4;
    %ntapers = 2*TW - 1;
    %params.Fs = 1000;
    %params.tapers = [TW ntapers];
    %params.fpass = [0 500];
    %[spike1_S, f] = mtspectrumpb(spike1, params);
    %[spike2_S, ~] = mtspectrumpb(spike2, params);

%Q5 

%Create Model Parameters
    spike_burn_in = 140; %need at least 140 lags to fit model
    hist_inds = [3:7, 120:140];
    n_hist_params = length(hist_inds);
    y1 = spike1(spike_burn_in+1:end);
    y2 = spike2(spike_burn_in+1:end);
    spike_inds1 = find(y1);
    spike_inds2 = find(y2);
    N1 = length(spike_inds1);
    N2 = length(spike_inds2);
    xHist1 = zeros(length(y1), n_hist_params);
    xHist2 = zeros(length(y1), n_hist_params);
    


    %Populate History Matrix
    i = 1;
    for j = (spike_burn_in + 1):length(spike1)
        xHist1(i, :) = spike1(j - hist_inds);
        xHist2(i, :) = spike2(j - hist_inds);
        i = i + 1;
    end
    
    %Fit Model
    [b1, dev1, stats1] = glmfit(xHist1, y1, 'poisson', 'log');
    [b2, dev2, stats2] = glmfit(xHist2, y2, 'poisson', 'log');
    lambda1 = exp(b1(1) + xHist1*b1(2:end));
    lambda2 = exp(b2(1) + xHist2*b2(2:end));
    
    %Time-Rescale
    Z1(1) = sum(lambda1(1:spike_inds1(1)));
    Z2(1) = sum(lambda2(1:spike_inds2(1)));
    for i = 2:N1
        Z1(i) = sum(lambda1(spike_inds1(i-1) + 1:spike_inds1(i)));
    end
    for i = 2:N2
        Z2(i) = sum(lambda2(spike_inds2(i-1) + 1:spike_inds2(i)));
    end
    [eCDF1, zvals1] = ecdf(Z1);
    [eCDF2, zvals2] = ecdf(Z2);
    mCDF1 = 1 - exp(-zvals1);
    mCDF2 = 1 - exp(-zvals2);
    
    %Model 1: KS
    ci1 = 1.36/sqrt(N1);
    ci2 = 1.36/sqrt(N2);
    figure()
    subplot(1, 2, 1)
    hold on
    plot(mCDF1, eCDF1, 'k', 'LineWidth', 2)
    plot([0 1], [0 1] + ci1, '--k', 'LineWidth', 1)
    plot([0 1], [0 1] - ci1, '--k', 'LineWidth', 1)
    hold off
    xlabel('Model CDF')
    ylabel('Emirical CDF')
    xlim([0 1])
    ylim([0 1])
    title('Spike 1 KS Plot')
    set(gca, 'FontSize', 14)
    subplot(1, 2, 2)
    hold on
    plot(mCDF2, eCDF2, 'k', 'LineWidth', 2)
    plot([0 1], [0 1] + ci2, '--k', 'LineWidth', 1)
    plot([0 1], [0 1] - ci2, '--k', 'LineWidth', 1)
    hold off
    xlabel('Model CDF')
    ylabel('Emirical CDF')
    xlim([0 1])
    ylim([0 1])
    title('Spike 2 KS Plot')
    set(gca, 'FontSize', 14)


