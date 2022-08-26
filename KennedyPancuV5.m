% Content: converting real/img(TRF)to data to aquire natural freq || -3dB Method for damping
%
% Date: 04/08/2022
%
% Author: Rohit Avadhani
%
% Version 6: this method has improved the identification of the middle
% point by interpolating directly the middle point of the maximum distance
% with the circle by using a line passing from the center of the circle. 
% It has been added the capacity to compute the damping/freq/psi with the
% mobility which is indicated as the better one. Moreover the code can also
% add some correction to center the first and the second natural
% frequencies
% In this version I have also added the -3db Method to check the data
% obtained with the Kennedy-Pancu Method in terms of damping.



%% Loading the experimental Data
clear
close all
clc
%Estimator H1
% (tip acceleration)
H1_12_r = table2array(readtable('H1_1,3_sv002.csv')); % freq - real - imag
% (base_acceleration)
H1_13_r = table2array(readtable('H1_1,2_sv002.csv')); % freq - real - imag

% This flag change the direction of the switch angle which is affected by
% 180 degree when the inertance is activated. Set equal to 0 for mobility
% an receptance. If we activate the mobility there is a shift of only 90
% deg.
inertance = 0;
mobility = 1;

% Delta frequency - change in freq sweep
dF = H1_12_r(2,1);

% Plotting H1
figure()
subplot(2,1,1)
plot(H1_12_r(:,1),angle(H1_12_r(:,2)+1i*H1_12_r(:,3))*180/pi,'-k','linewidth',2) % freq vs phase
hold on
plot(H1_13_r(:,1),angle(H1_13_r(:,2)+1i*H1_13_r(:,3))*180/pi,'-r','linewidth',2)
% hold on
% plot(H1_14_r(:,1),angle(H1_14_r(:,2)+1i*H1_14_r(:,3))*180/pi,'-b','linewidth',2)
grid on
box on
xlim([10 20])
xlabel('Frequency [Hz]', 'Interpreter','latex')
ylabel('Amplitude [deg]', 'Interpreter','latex')
title('Phase-H1', 'Interpreter','latex')
set(gca,'FontSize',10,'FontName','Times New Roman')
legend('Acc44-FM','Acc46-NFM')
subplot(2,1,2)
semilogy(H1_12_r(:,1),((H1_12_r(:,2)).^2+(H1_12_r(:,3)).^2).^0.5,'-k','linewidth',2)
hold on
semilogy(H1_13_r(:,1),((H1_13_r(:,2)).^2+(H1_13_r(:,3)).^2).^0.5,'-r','linewidth',2)
% hold on
% semilogy(H1_14_r(:,1),((H1_14_r(:,2)).^2+(H1_14_r(:,3)).^2).^0.5,'-b','linewidth',2)
grid on
box on
xlim([10 20])
xlabel('Frequency [Hz]', 'Interpreter','latex')
ylabel('Amplitude [m/s/N]', 'Interpreter','latex')
title('Amplitude-H1', 'Interpreter','latex')
set(gca,'FontSize',10,'FontName','Times New Roman')
legend('Acc44-FM','Acc46-NFM')

% Nyquist Plot H1
f1 = 8;
f2 = 15;
figure()
plot(H1_12_r(f1/dF:f2/dF,2),H1_12_r(f1/dF:f2/dF,3),'-k','linewidth',2)
hold on
plot(H1_13_r(f1/dF:f2/dF,2),H1_13_r(f1/dF:f2/dF,3),'-r','linewidth',2)
% hold on
% plot(H1_14_r(f1/dF:f2/dF,2),H1_14_r(f1/dF:f2/dF,3),'-b','linewidth',2)
box on
grid on
xlabel('Real Part', 'Interpreter','latex')
ylabel('Imaginary Part', 'Interpreter','latex')
title('Nyquist Plot - H1', 'Interpreter','latex')
set(gca,'FontSize',10,'FontName','Times New Roman')
legend('Acc44-FM','Acc46-NFM')



%% Extraction of parameter for the Kennedy-Pancu Method

f_ = H1_12_r(f1/dF:f2/dF,1);
FRF_H22 = H1_13_r(f1/dF:f2/dF,2)+1i*H1_13_r(f1/dF:f2/dF,3);
FRF_H12 = H1_12_r(f1/dF:f2/dF,2)+1i*H1_12_r(f1/dF:f2/dF,3);

%% NATURAL FREQUENCIES
% This block allows to get the natural frequencies in the Nyquist plot by
% looking at the point point where the there is the maximum distance
% between two consecutive points

% Correction - use only in the case the frequency is not precise
corr_1 = 1; % number of position to go forward or backward (first peak)
corr_2 = 0; % number of position to go forward or backward (second peak)

% Distance between points in the complex plane
distH12 = ((diff(real(FRF_H12))).^2+(diff(imag(FRF_H12))).^2).^0.5;


% Scaled frequency
frq_s = (f_(1:end-1)+dF/2);

% Position and value of the Antiresonace
[val_antR,pos_antR] = min(abs(FRF_H12));

% Plotting the distance
figure()
semilogy(frq_s,distH12)

xlabel('Frequency [Hz]')
ylabel('Distance')

% Getting the position and the value of the maximum distance points
[max_H12(1),idxH12(1)] = max(distH12(1:pos_antR));
[max_H12(2),idxH12(2)] = max(distH12(pos_antR:end));
idxH12(2) = numel(distH12(1:pos_antR-1))+idxH12(2);


% Correction - to be modified on the base of the case
idxH12(1) = idxH12(1)+corr_1;

idxH12(2) = idxH12(2)+corr_2;


% Getting the natural frequencies looking for the maximum distance points
fn_H12 = sort(frq_s(idxH12));


% Plot FRF in the Amplitude plot with the natural frequencies
figure()
semilogy(f_,mag2db(abs(FRF_H12)))

grid on;box on
xline(fn_H12,'-k')

xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
legend('FRF of Damper')

disp('-------------------------------------------------------------------')
disp('Natural Frequencies - Kennedy Pancu Method')
fprintf('H12-1    H12-2    \n')
fprintf('%.2f   %.2f', [fn_H12(1),fn_H12(2)])
fprintf('\n')
disp('-------------------------------------------------------------------')
%% receptance 
Rm = (H1_12_r(f1/dF:f2/dF,2));
Im = (H1_12_r(f1/dF:f2/dF,3));

A = [sum(Rm.^2)  sum(Rm.*Im) sum(Rm);
    sum(Rm.*Im) sum(Im.^2)  sum(Im);
    sum(Rm)     sum(Im)     numel(Rm)];
b = -[sum((Rm.^2+Im.^2).*Rm);
    sum((Rm.^2+Im.^2).*Im);
    sum((Rm.^2+Im.^2))];
xc = A\b;
X0_22A = xc(1)/(-2);
Y0_22A = xc(2)/(-2);
R_22A = (-xc(3)+X0_22A^2+Y0_22A^2)^0.5;

% Plot receptance in the Nyquist plot for the first circle
figure()
plot((H1_12_r(f1/dF:f2/dF,2)),(H1_12_r(f1/dF:f2/dF,3)),'-k','LineWidth',2)

hold on
th = 0:pi/50:2*pi;
xunit = R_22A * cos(th) + X0_22A;
yunit = R_22A * sin(th) + Y0_22A;
h = plot(xunit, yunit,'-o','LineWidth',2);
xline(800,'linewidth',3)
yline(0,'linewidth',3)
xlabel('Real Axis')
ylabel('Immaginary Axis')
grid on
title('plotting receptance of damper')
hold off




%% -3dB Method
% This block is used to double-check the damping ratio obtained with the
% -3dB Method

% H12 - First Peak
%maxdB = max(mag2db(abs(FRF_H12(1:pos_antR))));
maxdB = 1/2*(mag2db(abs(FRF_H12(idxH12(2))))+mag2db(abs(FRF_H12(idxH12(2)+1))));
min3dB = maxdB-3;
frA = interp1(mag2db(abs(FRF_H12(1:idxH12(2)))),f_(1:idxH12(2)),min3dB);
frB = interp1(mag2db(abs(FRF_H12(idxH12(2):281))),f_(idxH12(2):281),min3dB);

figure()
plot(f_,mag2db(abs(FRF_H12)),'LineWidth',1)
hold on
yline(maxdB,'LineWidth',2)
yline(min3dB,'LineWidth',2)
xline(fn_H12(2),'LineWidth',2)
xline(frA,'LineWidth',2)
xline(frB,'LineWidth',2)
grid on;box on
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
title('H12 - First mode/peak')
% 
zita_1_H12_HP = (frB-frA)/2/fn_H12(2);

disp('-------------------------------------------------------------------')
disp('Damping Ratio - -3dB Method')
fprintf('H12-1\n')
fprintf('%.4f', [zita_1_H12_HP])
fprintf('\n')
disp('-------------------------------------------------------------------')

%% Mode-Shape

% we define the sign on the base of the ratio of the center. If the
% inertance and the receptance are used Y centre is the one to be used
% while for the mobility we need the Y
% if mobility == 1
%     if X0_12A/X0_22A > 0
%         sign1 = 1;
%     else
%         sign1 = -1;
%     end
% 
%     if X0_12B/X0_22B > 0
%         sign2 = 1;
%     else
%         sign2 = -1;
%     end
% else
%     if Y0_12A/Y0_22A > 0
%         sign1 = 1;
%     else
%         sign1 = -1;
%     end
% 
%     if Y0_12B/Y0_22B > 0
%         sign2 = 1;
%     else
%         sign2 = -1;
%     end
% 
% end
% 
% PSI = [R_12A/R_12A R_12B/R_12B;
%        R_22A/R_12A*sign1 R_22B/R_12B*sign2]
% 
% 
% 
% %% Assumed Mass Matrix Results
% % We assume the mass matrix and we obtain the values of the matrices
% 
% % Real Mass matrix
% M = [189.6/1000 0;
%      0 189.6/1000*1]; % Obtained from measurement (it does not count the 
%                     % the dynamic effect and the effect of the shaker)
% 
% % Relative Damping Matrix
% Z = [2*zita_1_H22*fn_H12(1)*2*pi 0;
%      0 2*zita_2_H22*fn_H12(2)*2*pi];
% 
% % Wn Matrix
% Wn = [(fn_H22(1)*2*pi)^2 0;
%      0 (fn_H22(2)*2*pi)^2];
% 
% % Cm - Modal Damping Matrix
% Cm = PSI.'*M*PSI*Z;
% 
% % Km - Modal Stiffness Matrix
% Km = PSI.'*M*PSI*Wn;
% 
% % C - Real Damping Matrix
% C = inv(PSI.')*Cm*inv(PSI)
% 
% % K - Real Stiffness Matrix
% K = inv(PSI.')*Km*inv(PSI)
% 
% 
% k2 = -(K(1,2)+K(2,1))/2;
% k1 = K(1,1)-k2;
% k3 = K(2,2)-k2;
% 
% c2 = -(C(1,2)+C(2,1))/2;
% c1 = C(1,1)-c2;
% c3 = C(2,2)-c2;
% 
% function circ(x,y,r)
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% plot(xunit, yunit,'--');
% end
% 
% function [max_v,ind] = get_max(x)
% [max_v(1), ind(1)] = max(x);
% x(ind(1))      = -Inf;
% [max_v(2), ind(2)] = max(x);
% end