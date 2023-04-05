clc;
clear all ;
close all;
%%
DOF_data = importdata('DOF_measurements.txt');
I_data = importdata('I_measurements.txt')';
%%
DOF_data.data=DOF_data.data';
I_data.data = I_data.data';
%%
t_anim = linspace(1,120,width(DOF_data.data));
dt_anim = t_anim(2) - t_anim(1);
figure(100)
plot(t_anim , DOF_data.data(1:20,:))
title('Measured signals as a function of time')
xlabel("Time[sec]")
ylabel('Amplitude')
%%
% animate_mov_data(DOF_data.data)

%%
%P_neutral
%calculations
P_time_avg = calc_time_avg(DOF_data.data);
P_neutral = P_time_avg(:,end) ;

dp = DOF_data.data - P_neutral;
dp_dpt = dp*dp';
cov_mat_K_dp = calc_time_avg(dp_dpt);
[eig_vec , eig_val] = eig(dp_dpt);

%normalized Eig valuss calculation
eig_val_vec = zeros(1,width(eig_val));
for i = 1:width(eig_val)
    eig_val_vec(i) = eig_val(i,i);
end
%%

counter = 1:1:width(eig_val);
figure(101)
plot(counter,sort(eig_val_vec/norm(eig_val_vec),'descend'),'*')
title('normlized eig vals')
xlabel("index")
ylabel('normalized magnitude')
%signal matrix calulation
u = eig_vec \ dp;
dominant_control_signals = u(31:45,:);
%%
figure(102)
% plot(t_anim,dominant_control_signals);
plot(t_anim , u(45,:)')
hold on
% plot(t_anim , u(42,:)')
plot(t_anim , u(35,:)')
% plot(t_anim , u(36,:)')
plot(t_anim , u(30,:)')
hold off
legend('most significant signal', '10th' ,'15th ')
xlabel("Time[sec]")
ylabel('Amplitude')
%%
figure(103)
subplot(2,2,1)
title('hi')
plot(t_anim , I_data.data(1,:))
hold on
plot(t_anim , u(45,:)')
hold off
xlabel("Time[sec]")
ylabel('Amplitude')
legend('yaw rate','ctrl sig 1')
subplot(2,2,2)
plot(t_anim , I_data.data(1,:))
hold on
plot(t_anim , u(40,:)')
hold off
xlabel("Time[sec]")
ylabel('Amplitude')
legend('yaw rate','ctrl sig 5')

subplot(2,2,3)
plot(t_anim , I_data.data(1,:))
hold on
plot(t_anim , u(35,:)')
hold off
xlabel("Time[sec]")
ylabel('Amplitude')
legend('yaw rate','ctrl sig 10')

subplot(2,2,4)
plot(t_anim , I_data.data(1,:))
hold on
plot(t_anim , u(31,:)')
hold off
xlabel("Time[sec]")
ylabel('Amplitude')
sgtitle('Yaw rate VS control signal')
legend('yaw rate','ctrl sig 15')

%%
%4
four_ctrl_signals = u(42:45,:);
tau_vec = 0.1:0.1:10;
p_corr = calc_corr(dt_anim,tau_vec , four_ctrl_signals, I_data.data);
%%
sz = size(p_corr);
p_val= zeros(height(four_ctrl_signals),sz(3));
for l = 1:1:height(I_data.data)-3
    for j = 1:1:height(four_ctrl_signals)
        for i = 1:1:sz(3)
            p_val(j+(l-1)*height(four_ctrl_signals),i) = p_corr(j,l,i);
        end
    end
figure(111)
subplot(3,3,l)
plot(tau_vec ,p_val(1+(l-1)*height(four_ctrl_signals),:))
hold on
plot(tau_vec ,p_val(2+(l-1)*height(four_ctrl_signals),:))
plot(tau_vec ,p_val(3+(l-1)*height(four_ctrl_signals),:))
plot(tau_vec ,p_val(4+(l-1)*height(four_ctrl_signals),:))
hold off
legend("sig 1","sig 2","sig 3","sig 4")
xlabel('\tau [sec]')
ylabel('Amplitude')
title(strrep(I_data.textdata{l},"_"," "))
end
sgtitle('correlation coefficient as a function of time difference for each inertial motion')%KEEP YAW RATE V VERTICAL VBODY UP VBODY LEFT

%%
% Define the time axis
rel_rows = 1000;
alpha = 0.5;
beta = 100;
t_final = 100;
step_size= 1e-2;
standard_deviation = 1;
x_0 = 0;
t_0 = 0;
[x,t] = create_Langevin(x_0,t_0,alpha,t_final,step_size,standard_deviation,0);
len_rel = length(x);
data_matrix =zeros(rel_rows,len_rel) ;
z_arr = zeros(1,len_rel) ;
%%
figure(1)
    hold on
    for i = 1:1:rel_rows
        [x,t,z] = create_Langevin(x_0,t_0,alpha,t_final,step_size,standard_deviation,0);
        data_matrix(i,:) = x;
        z_arr(i,1:end-1) = z;
        plot(t,x)
        title("Realizations of the process")
        xlabel("Time[sec]")
        ylabel('Amplitude')
    end
hold off

%%
%Ensemble average calculation
X_ensemble_avg = mean(data_matrix(:,:));
figure(2)
plot(t,X_ensemble_avg)
title('Ensemble average')
xlabel('Time[sec]')
ylabel('amplitude')
rms_X_ensemble = rms(mean(data_matrix(:,:)));
%%
figure(3)
title('histogram')
histfit(data_matrix(:,end),100)
xlabel('Final value')
ylabel('# of realizations')
%%
X_time_avg = calc_time_avg(data_matrix);

%%
figure(4)
plot(t, X_time_avg)
title('time avrage of the process X')
xlabel("Time[sec]")
ylabel('Amplitude')
%%
rms_X_time = rms(mean(X_time_avg(:,end))); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 1.1.c
%numeric Autocorrelation calculatuon
% 1
theoretical_R_X_tt = (standard_deviation^2/(2*alpha))*(exp(-alpha*abs(0))-exp(-alpha*(2*t)));
R_X = mean(data_matrix.^2)*step_size;
%%
figure(5)
plot(t,theoretical_R_X_tt)
hold on
plot(t,R_X)
title('Ensemble average for the autocorrelation function of X')
xlim([0 20])
hold off
xlabel("Time[sec]")
ylabel('Amplitude')
legend('theroretical','analytic')
%%
%2
tau_sample_vec = 0:0.2:10;
theoretical_R_X_tau = (standard_deviation^2/(2*alpha))*(exp(-alpha*abs(tau_sample_vec)));
ts = 3/step_size;
time_mean_X_1 = calc_time_mean(tau_sample_vec , step_size , t , ts , data_matrix , 1);
time_mean_X_5 = calc_time_mean(tau_sample_vec , step_size , t , ts , data_matrix , 5);
time_mean_X_10 = calc_time_mean(tau_sample_vec , step_size , t , ts , data_matrix , 10);
time_mean_X_100 = calc_time_mean(tau_sample_vec , step_size , t , ts , data_matrix , 100);
%%
figure(6)
plot(tau_sample_vec, theoretical_R_X_tau , '--')
hold on

plot(tau_sample_vec ,mean(time_mean_X_1(:,ts:end),2))
plot(tau_sample_vec ,mean(time_mean_X_5(:,ts:end),2))
plot(tau_sample_vec ,mean(time_mean_X_10(:,ts:end),2))
plot(tau_sample_vec ,mean(time_mean_X_100(:,ts:end),2))
legend("theory" , '1 realizations ','5 realizations','10 realizations', '100 realizations')


xlim([0 10])
% legend("calculation - mean of all realizations" , 'calculation - a few realizations')
hold off
xlabel("\tau [sec]")
ylabel('Amplitude')
title('autocorrelation of X as a function of Tau')
%%
% [Pxx,~,~,f_new] = dft_periodogram(data_matrix(1:10,:),z_arr(1:10,:),1);
dft_matrix = create_dft_matrix(width(data_matrix)-ts+1);
%%
%1.1.D
%numerical sol
[f,peri_ensemble_X , peri_ensemble_Z ,peri_ensemble_X_Z, peri_X ,peri_Z,peri_X_Z] = calc_periodogram(step_size,t,ts,data_matrix,z_arr,dft_matrix);

%%
%analytical sol
figure(7)
semilogx(f , 20*log10(step_size*peri_X(1:4,:)))
title('realizations of PSD X ')
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')
S_X_theory = (standard_deviation/sqrt(step_size))^2*((alpha^2+(2*pi*f).^2)).^-1; 

figure(8)
semilogx(f , 20*log10(step_size*peri_ensemble_X))
hold on
semilogx(f , 20*log10(S_X_theory))
xlim([10^-2 10])
hold off
title('ensemble avg of PSD X')
legend('calculation','theory')
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')
%%
% 1.1.E
%analytic sol
S_X_Z_theory = (standard_deviation/sqrt(step_size))^2*(alpha+(2j*pi*f)).^-1;
figure(9)
semilogx(f , 20*log10(step_size*abs(peri_X_Z(1:4,:))))
title('realizations of cross PSD X,Z ')
xlim([10^-2 10])
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')

figure(10)
semilogx(f , 20*log10(abs(peri_ensemble_X_Z)))
hold on
semilogx(f , 20*log10(abs(S_X_Z_theory)))
xlim([10^-2 10])
hold off
title('ensemble average of Cross PSD XZ ')
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')
legend('calculation','theory')

%%
%coherence calc
coherence_X_Z_theo = abs(S_X_Z_theory).^2./(S_X_theory*standard_deviation^2); 
coherence_X_Z = abs(peri_ensemble_X_Z).^2./abs(peri_ensemble_X.*peri_ensemble_Z);
figure(11)
semilogx(f ,coherence_X_Z)
hold on
semilogx(f ,step_size*coherence_X_Z_theo)
hold off
legend('analytical','theoretical')
title('Calculated coherence as a function of frequency')
xlabel("Frequency [Hz]")
ylabel('Amplitude')
xlim([10^-2 10])
ylim([0.9 1.1])


%%
%------------------------------------------------------------------------
%1.2.A
%1.1.A.1
Y_data_matrix = zeros(height(data_matrix),width(data_matrix));
for j = 1:1:height(data_matrix)
    Y_data_matrix(j,:) = create_Langevin(x_0,t_0,beta,t_final,step_size,standard_deviation,data_matrix(j,:));
end

%%
figure(12)
plot(t,Y_data_matrix(1:round(rel_rows*0.5),:))
title("Realizations of the process")
xlabel("Time [sec]")
ylabel('Amplitude')
%%
figure(13)
title('histogram')
histfit(Y_data_matrix(:,end),100)
xlabel('Final value')
ylabel('# of realizations')
%%
%1.1.B.2
rms_Y_ensemble = rms(mean(Y_data_matrix(:,:)));
%%
figure(14)
plot(t,mean(Y_data_matrix))
hold on
plot(t,mean(Y_data_matrix(1:100,:)))
plot(t,mean(Y_data_matrix(1:10,:)))

hold off
title('Ensemble average')
legend('1000 realizations','100 realizations','10 realizations')
xlabel("Time [sec]")
ylabel('Amplitude')
%%
figure(15)
plot(t, calc_time_avg(Y_data_matrix))
title('time average as a function of time')
xlabel("Time [sec]")
ylabel('Amplitude')
%%
rms_Y_time = rms(mean(Y_time_avg(:,end)));
%%
%1.1.C
% 1
% R_Y_theo = -(standard_deviation^2)*(beta*exp(-alpha*(2*t*0))-alpha*exp(-beta*(2*t))) / ((alpha^2-beta^2)*(2*alpha*beta));
R_Y_theo = -(standard_deviation^2)*(beta*exp(-alpha*(2*t*0))-alpha*exp(-beta*(2*t))) / ((alpha^2-beta^2)*(2*alpha*beta));

R_Y = mean(Y_data_matrix.^2)*step_size;
figure(16)
plot(t,R_Y)
hold on 
plot(t,R_Y_theo)
title('Ensemble average for the autocorrelation function of Y')
xlim([0 40])
hold off
xlabel("Time [sec]")
ylabel('Amplitude')
%%
% 2
ts = 3/step_size;
tau_sample_vec_Y = 0:0.5:20;
time_mean_Y_1 = calc_time_mean(tau_sample_vec_Y , step_size , t , ts , Y_data_matrix , 1);
time_mean_Y_5 = calc_time_mean(tau_sample_vec_Y , step_size , t , ts , Y_data_matrix , 5);
time_mean_Y_10 = calc_time_mean(tau_sample_vec_Y , step_size , t , ts , Y_data_matrix , 10);
time_mean_Y_100 = calc_time_mean(tau_sample_vec_Y , step_size , t , ts , Y_data_matrix , 100);
% time_mean_Y_1000 = calc_time_mean(tau_sample_vec_Y , step_size , t , ts , Y_data_matrix , 1000);

%%
% R_Y_theo_tau = -(standard_deviation^2)*(beta*exp(-alpha*(tau_sample_vec_Y))-alpha*exp(-beta*(2*t+tau_sample_vec_Y))) / ((alpha^2-beta^2)*(4*alpha*beta));
% R_Y_theo_tau = -(standard_deviation^2)*(beta*exp(-alpha*(tau_sample_vec_Y))) / ((alpha^2-beta^2)*(4*alpha*beta));
R_Y_theo_tau = -(standard_deviation^2)*(beta*exp(-alpha*(tau_sample_vec_Y))-alpha*exp(-beta*(tau_sample_vec_Y))) / ((alpha^2-beta^2)*(2*alpha*beta));

figure(17)
plot(tau_sample_vec_Y,R_Y_theo_tau ,'--')

hold on
plot(tau_sample_vec_Y ,mean(time_mean_Y_1(:,ts:end),2))
plot(tau_sample_vec_Y ,mean(time_mean_Y_5(:,ts:end),2))
plot(tau_sample_vec_Y ,mean(time_mean_Y_10(:,ts:end),2))
plot(tau_sample_vec_Y ,mean(time_mean_Y_100(:,ts:end),2))
% plot(tau_sample_vec_Y ,mean(time_mean_Y_1000(:,ts:end),2))
legend("theory" , '1 realizations ','5 realizations','10 realizations', '100 realizations')

xlim([0 20])
% ylim([-1 2])
hold off
xlabel("\tau [sec]")
ylabel('Amplitude')
title('autocorrelation of Y as a function of Tau')
%%
%1.1.D
[f,peri_ensemble_X,peri_ensemble_Y ,peri_ensemble_X_Y,peri_X,peri_Y,peri_X_Y] = calc_periodogram(step_size,t,ts,data_matrix,Y_data_matrix,dft_matrix);
figure(18)
semilogx(f , 20*log10(step_size*peri_Y(1:4,:)))
title('realizations of PSD Y ')
xlim([10^-2 10])
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')

figure(19)
semilogx(f , 20*log10(step_size*peri_ensemble_Y))
xlim([10^-2 10])
title('ensamble average of PSD Y')
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')

%%
%1.1.E
figure(20)
semilogx(f , 20*log10(abs(step_size*peri_X_Y(1:4,:))))
title('realizations of cross PSD X,Y ')
xlim([10^-2 10])
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')
figure(21)
semilogx(f , 20*log10(abs(peri_ensemble_X_Y)))
xlim([10^-2 10])
title('ensemble average of PSD YX')
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')
%%
%coherence
coherence_X_Y = calc_coherence(peri_ensemble_X,peri_ensemble_Y,peri_ensemble_X_Y);
figure(22)
semilogx(f ,coherence_X_Y)
title('Calculated coherence of PSD XY as a function of frequency')
xlim([10^-2 10])
xlabel("Frequency [Hz]")
ylabel('Amplitude [dB]')
%%
R_X_Y = mean(Y_data_matrix.*data_matrix)*step_size;
R_X_Y_theo_tt =(standard_deviation^2/(2*alpha*(beta-alpha))).*((exp(2*alpha.*t)-1).*(exp(beta.*t)-exp(alpha.*t)).*exp(-beta.*t-2*alpha.*t)) ;
figure(122)
plot(t,R_X_Y)
title('Ensemble average for the corss correlation function of X,Y')
hold on
plot(t,R_X_Y_theo_tt)
hold off
xlim([0 40])
xlabel("Time [sec]")
ylabel('Amplitude')
%%
%------------------------------------------------------------------------
%1.2.B
ln_data_matrix = zeros(height(data_matrix),width(data_matrix));
m = 1;
v = 1;
vd = v/step_size;
log_normal_mu = log(m^2/sqrt(vd+m^2));
log_normal_sigma = sqrt(log(1+(vd/m^2)));
%%
pd=makedist('lognormal', 'mu', log_normal_mu, 'sigma', log_normal_sigma);
ln_Z=random(pd,height(ln_data_matrix),length(t));

%%
%1.1.A
for j = 1:1:height(ln_data_matrix)
    ln_data_matrix(j,:) = create_Langevin(x_0,t_0,beta,t_final,step_size,standard_deviation,ln_Z(j,:));
end
%%
figure(23)
        plot(t,log(ln_data_matrix(5:round(rel_rows*0.02),:)))
        title("Realizations of the process")
        xlabel("Time")
        ylabel('Amplitude')
        %%
figure(24)
title('histogram')
histfit(log(ln_data_matrix(:,end)),100)
xlabel('Final value')
ylabel('# of realizations')
dist_log_normal = fitdist(log(ln_data_matrix(:,end)),'Normal');
%%
%1.1.B
ln_ensemble_avg = mean(ln_data_matrix);
%%
figure(25)
plot(t,ln_ensemble_avg)
title('Ensemble average with log normal')
xlabel("Time [sec]")
ylabel('Amplitude')
xlim([0 30])
%%
rms_ln_ensemble = rms((ln_ensemble_avg(:,2:end)));
%%
ln_time_avg = calc_time_avg(ln_data_matrix);
%%
figure(26)
plot(t, ln_time_avg)
title('time average as a function of time')
xlabel("Time [sec]")
ylabel('Amplitude')
%%
rms_ln_time = rms(mean(ln_time_avg(:,end)));
%%
% 1.1.c
%numeric Autocorrelation calculatuon
% 1
tau_sample_vec_ln = 0:0.1:15;
R_X_ln = mean(ln_data_matrix.^2);
%%
R_X_theo_ln = log_normal_mu^2/alpha^2 + (log_normal_sigma/(2*alpha))*exp(-alpha*abs(tau_sample_vec_ln)) ;
theoretical_R_X_tt_ln = (m^2/alpha^2)*(1-2*exp(-alpha.*t)+exp(-2*alpha.*t))+ (vd/(2*alpha))*(1-exp(-2*alpha.*t));
R_X_theo_ln_2 = m^2/alpha^2 + (v/(2*alpha))*exp(-alpha*abs(tau_sample_vec_ln)) ;
% theoretical_R_X_tt_ln_2 = exp(-alpha*2.*t).*(((vd/(2*alpha)/sqrt(step_size)).*(exp(alpha.*t)-1))+(m^2/alpha^2).*(exp(alpha.*t)-1).^2);

%%
figure(27)
plot(t,R_X_ln)
hold on
plot(t,theoretical_R_X_tt_ln)
hold off
title('Ensemble average for the autocorrelation function of X')
xlim([0 150])
xlabel("Time [sec]")
ylabel('Amplitude')
% 2
%%
ts = 1.5/step_size;
time_mean_ln_X = calc_time_mean(tau_sample_vec_ln , step_size , t , ts , ln_data_matrix , rel_rows);
%%
figure(28)
% plot(tau_sample_vec_ln , mean(time_mean_ln_X(:,200:200:1000),2),"-*",'LineWidth',1.5)
hold on
% plot(tau_sample_vec_ln , time_mean_ln_X(1,:))
plot(tau_sample_vec_ln , R_X_theo_ln_2*step_size,"o",'LineWidth',1.5)
plot(tau_sample_vec_ln , (time_mean_ln_X(:,round(length(time_mean_ln_X)/20):round(length(time_mean_ln_X)/20):round(length(time_mean_ln_X)/2))))

% xlim([0 10])
title('autocorrelation as a function of time difference')
ylabel('amplitude')
xlabel('\tau [sec]')
legend('theory', 'calculation - a few realizations')
hold off
%%
%1.1.D
%%
[f,peri_ensemble_ln_X,peri_ensemble_ln_Z ,peri_ensemble_ln_X_Z,peri_ln_X,peri_ln_Z,peri_ln_X_Z] = calc_periodogram(step_size,t,ts,ln_data_matrix,ln_Z);
%%
S_X_theory_ln = (vd*(alpha^2+4*pi^2.*f.^2).^-1)+((m^2/alpha^2)); 

%%
figure(29)
semilogx(f , 20*log10(step_size*peri_ln_X(1:4,:)))
title('realizations of PSD X ')
xlim([10^-2 10])
xlabel("Frequency [Hz]")
ylabel('Amplitude')
figure(30)
semilogx(f , 20*log10(step_size*peri_ensemble_ln_X))
hold on
semilogx(f , 20*log10(step_size*S_X_theory_ln))
hold off
xlim([10^-2 10])
title('ensamble average of PSD X')
xlabel("Frequency [Hz]")
ylabel('Amplitude')
legend('calculation','theory')
%%
%1.1.E
S_X_Z_theory_ln = (vd*(alpha+2*pi.*f*1j).^-1)+((m^2/alpha));
%%
figure(31)
semilogx(f , 20*log10(abs(step_size*peri_ln_X_Z(1:4,:))))
title('realizations of cross PSD X,Z ')
xlim([10^-2 10])
xlabel("Frequency [Hz]")
ylabel('Amplitude')
figure(32)
semilogx(f , 20*log10(abs(peri_ensemble_ln_X_Z)))
hold on
semilogx(f , 20*log10(abs(S_X_Z_theory_ln)))
hold off
xlim([10^-2 10])
title('ensemble average of CPSD XZ')
legend('calculation','theory')
xlabel("Frequency [Hz]")
ylabel('Amplitude')
%%
%coherence
coherence_ln_X_Z = calc_coherence(peri_ensemble_ln_X,peri_ensemble_ln_Z,peri_ensemble_ln_X_Z);
%%
S_Z_ln = m^2+v;
coherence_theo_ln = abs(S_X_Z_theory_ln).^2./(S_Z_ln.*S_X_theory_ln);
%%
figure(33)
semilogx(f ,coherence_ln_X_Z)
hold on
% semilogx(f ,coherence_theo_ln)
fplot(1)
hold off
title('Calculated coherence of PSD XZ as a function of frequency')
xlim([10^-2 10])
ylim([0.95 1.05])
legend('calculation','theory')
xlabel("Frequency [Hz]")
ylabel('Amplitude')
%%
%this function solves the langevin equation.
%it takes initial conditions , alpha ,t final, step size , standard
%deviation and an input, for white gaussian noise with mean of zero and std
%of the input standard deviation enter zero.
function[x,t,z] =create_Langevin(x_0,t_0,alpha,t_fin,del_t,sig,input)

    t = t_0:del_t:t_fin;
    sig_d = sig/sqrt(del_t);

    x = zeros(1,length(t));
    z =sig_d*sqrt(del_t)*randn(1,length(t)-1);
    x(1) = x_0;
    t(1) = t_0;    

    if mean(input) ~=0
        for i = 1:length(t)-1
            x(i+1) = x(i) * (1-alpha*del_t) + del_t*input(i);
        end
    else
        for i = 1:length(t)-1
            x(i+1) = x(i) * (1-alpha*del_t) + z(i);
        end
    end
end
function[time_avg] = calc_time_avg(data_matrix)
    [height , width] = size(data_matrix);
    p = zeros(height, width);
    time_avg = zeros(1,height);
    for i = 1:1:height
        p(i,:) = cumsum(data_matrix(i,:));
        for j = 1:1:length(p)
            time_avg(i,j) = p(i,j)/j;
        end
    end
end
function[f,peri_ensemble_data_1 ,peri_ensemble_data_2 , peri_ensemble_data_1_2 , peri_data_1,peri_data_2,peri_data_1_2] = calc_periodogram(step_size,t,ts,data_1,data_2,dft_matrix)
    N = length(t(ts:end));%samples in time
    K = 1:1:N;%pos in list
    R = 1/step_size;%sampling rate
%     K = 1:1:(10*N/R);%pos in list
    f = (R*K)/N;
%     f = 0.01:0.1:10;
%     % Compute the frequency bins
%     N = size(data_1, 2);
%     Fs = 1; % Sampling frequency (Hz)
%     f = Fs*(0:(N/2))/N;
    
    data_1_after_noise = data_1(:,ts:end);
%     transformed_data_1_arr = fft(data_1_after_noise');
%     [transformed_data_1_arr,~] = dft(data_1_after_noise,t,f,step_size,ts);
%     transformed_data_1_arr = transformed_data_1_arr';
    [transformed_data_1_arr] = (dft_matrix*(data_1_after_noise(:,:)'))';

    % preallocation 
    peri_data_1 = (1/width(transformed_data_1_arr)) * abs(transformed_data_1_arr).^2;
    peri_ensemble_data_1 = mean(peri_data_1);

    if mean(mean(data_2),2) ~= 0
        data_2_after_noise = data_2(:,ts:end);
%         transformed_data_2_arr = fft(data_2_after_noise');
%         transformed_data_2_arr = transformed_data_2_arr';
%         [transformed_data_2_arr,~] = dft(data_2_after_noise,t,f,step_size,ts);
%         [transformed_data_2_arr] = comp_dft(data_2_after_noise(1:10,:));
        [transformed_data_2_arr] = (dft_matrix*(data_2_after_noise(:,:)'))';


        
        peri_data_1_2 = (1/width(data_1_after_noise)) * (transformed_data_1_arr.*conj(transformed_data_2_arr));
        peri_data_2 = (1/width(transformed_data_2_arr)) * abs(transformed_data_2_arr).^2;
        peri_ensemble_data_1_2 = mean(peri_data_1_2);
        peri_ensemble_data_2 = mean(peri_data_2);
    else
        peri_ensemble_data_1_2 = 0;
        peri_ensemble_data_2 = 0;
        peri_data_1_2 = 0;
        peri_data_2 = 0;
    end
end
function[time_mean] = calc_time_mean(tau_samples,step_size ,t,ts,data_matrix ,realizations)
    tau_vec_normlized = round(tau_samples/step_size);
    time_mean = zeros(length(tau_vec_normlized),length(t));
    for i = 1:1:length(tau_vec_normlized)
        for j = ts:length(t)
            if(j+tau_vec_normlized(i)<width(data_matrix))
                time_mean(i,j) = step_size* mean(data_matrix(1:realizations,j).*data_matrix(1:realizations,j+tau_vec_normlized(i)));
            else
                continue
            end
        end
    end


end
function[coherence_1_2] = calc_coherence(S_1,S_2,S_1_2)
        coherence_1_2 = abs(S_1_2).^2./abs(S_1.*conj(S_2));
end
function[corr_coeff] = calc_corr(step_size,tau_samples,u,z)   
    tau_vec_normlized = round(tau_samples/step_size);
    t = width(z);
    cov_x_z = zeros(height(u),height(z),length(tau_samples));
    var_x =  zeros(1,height(u));
    var_z =  zeros(1,height(z));
    for i = 1:1:length(tau_vec_normlized)
        for j = 1:1:height(z)
            var_z(j) = mean(z(j,:).*z(j,:)) - mean(z(j,:)).*mean(z(j,:));
    
            for k = 1:1:height(u)
                var_x(k) =  mean(u(k,:).*u(k,:)) - mean(u(k,:)).*mean(u(k,:));
    
                if(j+tau_vec_normlized(i)<t)
                    cov_x_z(k,j,i) =mean(u(k,1:end-tau_vec_normlized(i)).*z(j,tau_vec_normlized(i):end-1)) - mean(u(k,1:end-tau_vec_normlized(i)))*mean(z(j,tau_vec_normlized(i):end-1));
                else
                    continue
                end
            end
        end
    end
    c_x_c_z = sqrt(var_x'.*var_z);
    corr_coeff = cov_x_z./c_x_c_z;
end

function[transformed_data,transformed_data_tweaked] = dft(data,t,f,step_size,ts)
    transformed_data = zeros(height(data),length(f)) ;
    dft_mat = zeros(data);
    n = 0:1:length(f)-1;
    for i = 1:1:10
        for j = 1:1:width(f)
            transformed_data(i,j) = sum(data(i,:).* exp( (-2j * pi .* j .* n(j)) ./ length(f) ));
        end
        disp(i)
    end
    transformed_data_tweaked = sum(transformed_data,2)*step_size;

%     transformed_data(i,j) = data(:,i).* exp( (-2j * pi .* 1/length(f) .* t(i)) ./ length(f) );

end


% Compute the DFT coefficients for each signal
function[dft_coeffs,f] = comp_dft(data)
    dft_coeffs = zeros(size(data));
    for i = 1:size(data,1)
        for k = 1:size(data,2)
            dft_coeffs(i,k) = sum(data(i,:).*exp(-2*pi*1i*(k-1)*(0:size(data,2)-1)/(size(data,2))));
        end
    end

end
function X = dft_single(x)
% DFT calculation without using built-in function
% Input: x - Input signal
% Output: X - DFT of the input signal

N = length(x);  % Length of the signal
X = zeros(1, N); % Allocate memory for the output signal

    for k = 1:N % For each frequency component
        for n = 1:N % For each time sample
            X(k) = X(k) + x(n)*exp(-1j*2*pi*(k-1)*(n-1)/N); % DFT equation
        end
    end
end



function [Pxx,Pyy,Pxy, f] = dft_periodogram(x, y, fs)
% DFT calculation and periodogram without using built-in function
% Inputs:
%   x - Input signal matrix of size M by N
%   y - Input signal matrix of size M by N
%   fs - Sampling rate in Hz
% Outputs:
%   Pxx - Power spectral density estimate (periodogram)
%   f - Frequency vector in Hz
    
    [M, N] = size(x);  % Size of the input signal matrix
    Pxx = zeros(M, N); % Allocate memory for the output periodogram

    for m = 1:M % For each row of the input signal matrices
        for n = 1:N % For each time sample
            X = 0; Y = 0;
            for k = 1:N % For each frequency component
                X = X + x(m, k)*exp(-1j*2*pi*(k-1)*(n-1)/N); % DFT equation for signal x
                Y = Y + y(m, k)*exp(-1j*2*pi*(k-1)*(n-1)/N); % DFT equation for signal y
            end
            Pxx_yy(m, :) = Pxx_yy(m, :) + (abs(X)^2 + abs(Y)^2)/(N*fs); % Periodogram calculation
            Pxx(m, :) = Pxx(m, :) + (abs(X)^2)/(N*fs); % Periodogram calculation
            Pyy(m, :) = Pyy(m, :) + (abs(X)^2 + abs(Y)^2)/(N*fs); % Periodogram calculation
            Pxy(m, :) = Pxy(m, :) + (abs(X * Y))/(N*fs); % Periodogram calculation
        end
        f = (0:N/2) * fs / N; % Frequency vector in Hz
    end
end

function DFT_matrix = create_dft_matrix(N)
    % Create a vector of indices from 0 to N-1
    n = 0:N-1;
    
    % Calculate the DFT matrix using the formula
    W = exp(-1i*2*pi/N);
    DFT_matrix = zeros(N,N);
    for k=0:N-1
        for j=0:N-1
            DFT_matrix(k+1,j+1) = W^(k*j);
        end
    end
end