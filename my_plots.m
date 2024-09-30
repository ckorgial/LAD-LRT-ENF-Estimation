
load ladreg_sota_comparison_recordings_all_results.mat

%% 1.

[~, order] = sort(MSE_LAD_E_single,'descend');
figure(1);
pf = plot(1:130, MSE_LAD_E_single(order),'r',...
    1:130, MSE_E_single(order),'b',...
    1:130, MSE_single(order),'g');
set(gca, 'YScale', 'log');
grid on; 
pf(1).LineWidth=2;
axis([0 131 10^(-6) 10^(-1)]);
xlabel('Recording Index');ylabel('MSE');
leg = legend('${\rm{LAD-E-single}}$','${\rm{E-single}}$','${\rm{single}}$');
set(leg,'Interpreter','latex');
%% 2.

[~, order] = sort(MSE_LAD_P_MLE,'descend');
figure(2);
pf = plot(1:130, MSE_LAD_P_MLE(order),'r',...
    1:130, MSE_P_MLE(order),'b',...
    1:130, MSE_MLE(order),'g');
set(gca, 'YScale', 'log');
grid on; 
pf(1).LineWidth=2;
axis([0 131 10^(-6) 10^(-1)]);
xlabel('Recording Index');ylabel('MSE');
leg = legend('${\rm{LAD-P-MLE}}$','${\rm{P-MLE}}$','${\rm{MLE}}$');
set(leg,'Interpreter','latex');
%% 3. 

[~, order] = sort(MSE_LAD_P_WMLE,'descend');
figure(3);
pf = plot(1:130, MSE_LAD_P_WMLE(order),'r',...
    1:130, MSE_P_WMLE(order),'b',...
    1:130, MSE_WMLE(order),'g');
set(gca, 'YScale', 'log');
grid on; 
pf(1).LineWidth=2;
axis([0 131 10^(-6) 10^(-1)]);
xlabel('Recording Index');ylabel('MSE');
leg = legend('${\rm{LAD-P-WMLE}}$','${\rm{P-WMLE}}$','${\rm{WMLE}}$');
set(leg,'Interpreter','latex');

%% 4.

index = 41;

estimatedFreq1 = f_LAD_P_MLE{index};
estimatedFreq2 = f_P_MLE{index};
referenceFreq  = f_ref{index};

offset_red = 0.01; 
offset_blue = 0.07; 
estimatedFreq1 = estimatedFreq1 + offset_red;
estimatedFreq2 = estimatedFreq2 + offset_blue;
referenceFreq = referenceFreq - offset_blue;

figure;
hold on;
plot(1:length(estimatedFreq1), estimatedFreq1, 'r', 'LineWidth', 1);
plot(1:length(estimatedFreq2), estimatedFreq2, 'g', 'LineWidth', 1);
plot(1:length(referenceFreq), referenceFreq, 'b', 'LineWidth', 1);
xlabel('Frame Index');
ylabel('Frequency (Hz)');
leg = legend('${\rm{LAD-P-MLE}}$', '${\rm{P-MLE}}$', 'Ground Truth ENF');
grid on;
hold off;
set(leg,'Interpreter','latex');