%% Modelling & Control of a 3DOF Helicopter
% Erfan Saremi
%%
clear
close all
clc
tic
%% Defining variables
% Rotational angles & Thrusts
syms eps tow rho real
syms deps dtow drho real
syms ddeps ddtow ddrho real
syms F_b F_f omega_f omega_b real
%%
% Masses and gravitional acceleration
g = 9.81;
M_h = 3.15;
M_w = 4.533;
%%
% Lenghs
L_w = 0.74;
L_m = 1.17;
L_h = 0.20;
L_rho = 0.078;
%%
% Inertias
I_rho = 0.11;
I_eps = 6.78;
syms I_tow real
Itow = @(eps) -6.37 * ( eps ^ 2 ) + 1.032 * eps + 7.196;
%%
% Grvaitaional momentum
syms M_g real
Mg = @(eps) -2.8834 * ( eps ^ 2 ) + 5.2111 * eps + 7.1683;
%%
% Friction & air resistance coefficients
k_d_rho = 0.053;
k_s_rho = 0.04;
k_d_eps = 1.78;
k_s_eps = 0;
k_d_tow = 0.5;
k_s_tow = 0;
%% 
% Friction & air resistance momentums
M_f_rho = k_d_rho * drho + sign(drho) * k_s_rho;
M_f_eps = k_d_eps * deps + sign(deps) * k_s_eps;
M_f_tow = k_d_tow * deps + sign(dtow) * k_s_tow;
%%  Dynamic Model
F_f = ((7.26 * (10 ^ (-8))) * (omega_f ^ 2)) - (5.93 * (10 ^ (-5)) * omega_f);
F_b = ((7.26 * (10 ^ (-8))) * (omega_b ^ 2)) - (5.93 * (10 ^ (-5)) * omega_b);
% u_sum = F_f + F_b;
% u_diff = F_b - F_f;
syms u_sum u_diff real
ddrho = ( u_diff * L_h - M_h * g * L_rho * sin(rho) * cos(eps) - M_f_rho) / I_rho;
ddeps = ((u_sum * L_m * cos(rho) - M_g - M_f_eps) - ((I_tow - I_rho) * (dtow ^ 2) * cos(eps) * sin(eps))) / I_eps;
ddtow = ( - u_sum * L_m * cos(eps) * sin(rho) - M_f_tow) / I_tow;
%% Linearized State Space Model
x_1 = rho;
x_2 = drho;
x_3 = eps;
x_4 = deps;
x_5 = tow;
x_6 = dtow;
dx_1 = x_2;
dx_2 = ddrho;
dx_3 = x_4;
dx_4 = ddeps;
dx_5 = x_6;
dx_6 = ddtow;
dx = [dx_1;
    dx_2;
    dx_3;
    dx_4;
    dx_5;
    dx_6];
x = [x_1;
    x_2;
    x_3;
    x_4;
    x_5;
    x_6];
u = [u_sum;
    u_diff];

A = [0, 1, 0, 0, 0, 0;
    -(M_h * g * L_rho)/I_rho, -(k_d_rho/I_rho), 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0;
    0, 0, -(Mg(0)/I_eps), -(k_d_eps/I_eps), 0, 0;
    0, 0, 0, 0, 0, 1;
    -Mg(0) / Itow(0), 0, 0, 0, 0, -k_d_tow / Itow(0)];
B = jacobian(dx, u);
B = subs(B, [rho, tow, eps], [0, 0, 0]);
B = double(B);
C = [0 , 0, 1, 0, 0, 0;
    0 , 0, 0, 0, 1, 0];
D = 0;
%% State Space Model and Transfer Function
state_space = ss(A, B, C, D);
transfer_function = tf(state_space);
[num, den] = tfdata(transfer_function, 'v');
figure, step(state_space, 30);
title("\fontsize{18}\fontname{Times}Step Response of the Whole System");
xlabel("\fontsize{16}\fontname{Times}Time");
ylabel("\fontsize{16}\fontname{Times}Amplitude");
figure, impulse(state_space)
%% Giving Controling Inputs
% Custom Inputs
tf_1 = transfer_function(1,1);
tf_4 = transfer_function(2,2);

[num_1,dem_1] = tfdata(tf_1, 'v');
[num_2,dem_2] = tfdata(tf_4, 'v');

s = sym('s');

nH = poly2sym(num_1, s);
dH = poly2sym(dem_1, s);

H = simplify(nH / dH);

U_1 = 1/s;
U_2 = 23/s^2 + 6*s -10;

Y_1 = U_1 * H;
Y_2 = U_2 * H;

y_1 = ilaplace(Y_1);
y_1 = simplify(y_1);
y_2 = ilaplace(Y_2);
y_2 = simplify(y_2);

t = linspace(0, 100, 1000);
% figure,
% p = ezplot(y_1,t), ylim([0, 0.3]), xlim([0, 40]);
% p.Color = 'black';
% p.LineWidth = 2.5;
% p.LineStyle = '-.';
% hold on ;
% step(tf_1, 40);

f = matlabFunction(y_1);
y_1_2 = f(t);

plot(t,y_1_2, 's',"MarkerSize",3)
legend('\fontsize{16}\fontname{Times}Giving Step Function Manualy to TF in Frequense Space','\fontsize{16}\fontname{Times}From Step Function', '\fontsize{16}\fontname{Times}Giving Step Function Manualy to TF in Time Space');
title('');
ylabel("\fontsize{16}\fontname{Times}Amplitude");
xlabel("\fontsize{16}\fontname{Times}Time");
grid on;
hold off;
figure, 
ezplot(y_2,t), ylim([-5, 70]), xlim([0, 20]);
ylabel("\fontsize{16}\fontname{Times}Amplitude");
xlabel("\fontsize{16}\fontname{Times}Time(Seconds)");
title("\fontsize{18}\fontname{Times}Custom Input");
grid on
%% Separating the systems (Decoupling)
tf_travel = transfer_function(2,2); % Input: u_diff
tf_elevation = transfer_function(1,1); % Input: u_sum
ss_travel = ss(tf_travel);  % Input: u_diff
ss_elevation = ss(tf_elevation); % Input: u_sum

figure, 
subplot(1,2,1)
hold on
step(ss_elevation, 60);
impulse(ss_elevation);
title("\fontsize{18}\fontname{Times}Elevation Response to Step and Impulse Inputs");
legend("\fontsize{16}\fontname{Times}Step Response", "\fontsize{16}\fontname{Times}Impulse Response");
ylabel("\fontsize{16}\fontname{Times}Amplitude");
xlabel("\fontsize{16}\fontname{Times}Time");
grid on
hold off
subplot(1,2,2)
hold on
step(ss_travel, 200);
impulse(ss_travel);
title("\fontsize{18}\fontname{Times}Travel Response to Step and Impulse Inputs");
legend("\fontsize{16}\fontname{Times}Step Response", "\fontsize{16}\fontname{Times}Impulse Response");
ylim([-250 100]);
ylabel("\fontsize{16}\fontname{Times}Amplitude");
xlabel("\fontsize{16}\fontname{Times}Time");
grid on
hold off
stepinfo(ss_travel)
stepinfo(ss_elevation)
%% Linear simulation result 
tt = linspace(0,30,1000);
uu = exp(-tt);
figure, 
subplot(1,2,1)
lsim(tf_elevation,uu,tt);
title("\fontsize{18}\fontname{Times}Linear Simulation Result of Elevation");
ylabel("\fontsize{16}\fontname{Times}Amplitude");
xlabel("\fontsize{16}\fontname{Times}Time");
grid on
subplot(1,2,2)
lsim(tf_travel,uu,tt);
title("\fontsize{18}\fontname{Times}Linear Simulation Result of Travel");
ylabel("\fontsize{16}\fontname{Times}Amplitude");
xlabel("\fontsize{16}\fontname{Times}Time");
grid on
%% Initials
x0_elevation = [1, 0]';
x0_travel = [1, 0, 1 ,0]';
figure, 
subplot(1,2,1)
initial(ss_elevation, x0_elevation)
title("\fontsize{18}\fontname{Times}Initial Result of Elevation");
ylabel("\fontsize{16}\fontname{Times}Amplitude");
xlabel("\fontsize{16}\fontname{Times}Time");
grid on
subplot(1,2,2)
initial(ss_travel, x0_travel);
title("\fontsize{18}\fontname{Times}Initial Result of Travel");
ylabel("\fontsize{16}\fontname{Times}Amplitude");
xlabel("\fontsize{16}\fontname{Times}Time");
grid on
%% Characteristic Equetion Roots
[num_t, den_t] = tfdata(tf_travel, 'v');
[num_e, den_e] = tfdata(tf_elevation, 'v');
[zeros_travel, poles_travel, gain_travel] = tf2zpk(num_t, den_t);
[zeros_elevation, poles_elevation, gain_elevation] = tf2zpk(num_e, den_e);
%% Routh-Hurwitz Stability
Routh_Hurwitz_Stability(den_e);
Routh_Hurwitz_Stability(den_t);
%% Root Locus and Bode Diagram (Control Designer)
s = tf('s');

% Elevation
C_CSD_elevation =  (1.2358e12 * (s+1120)) / (s+1.126e06);
figure,
subplot(2,2,1), rlocus(feedback(series(C_CSD_elevation,tf_elevation),1));
grid on;ylabel("\fontsize{16}\fontname{Times}Imaginary Axis");xlabel("\fontsize{16}\fontname{Times}Real Axis");
title('\fontsize{18}\fontname{Times}Roots Locus Diagram for Eelevation');
subplot(2,2,2), step(feedback(C_CSD_elevation*tf_elevation, 1));
title('\fontsize{18}\fontname{Times}Step Input for Eelevation');
grid on;ylabel("\fontsize{16}\fontname{Times}Amplitude");xlabel("\fontsize{16}\fontname{Times}Time");
subplot(2,2,3), step(feedback(C_CSD_elevation,tf_elevation));
title('\fontsize{18}\fontname{Times}Control Signal for Eelevation');
grid on;ylabel("\fontsize{16}\fontname{Times}Amplitude");xlabel("\fontsize{16}\fontname{Times}Time");
subplot(2,2,4), impulse(feedback(C_CSD_elevation*tf_elevation, 1));
title('\fontsize{18}\fontname{Times}Impulse input for Eelevation');
grid on;ylabel("\fontsize{16}\fontname{Times}Amplitude");xlabel("\fontsize{16}\fontname{Times}Time");

% Travel
C_CSD_travel =  (-14.991 * (s+0.05942) * (s+0.03009)) /(s * (s+7.973));
figure,
subplot(1,2,1), bode(feedback(series(C_CSD_travel,tf_travel),1));
grid on;xlabel("\fontsize{16}\fontname{Times}Frequency");
title('\fontsize{18}\fontname{Times}Roots Locus Diagram for Travel');
subplot(3,2,2), step(feedback(C_CSD_travel*tf_travel, 1));
title('\fontsize{18}\fontname{Times}Step Input for Travel');
grid on;ylabel("\fontsize{16}\fontname{Times}Amplitude");xlabel("\fontsize{14}\fontname{Times}Time");
subplot(3,2,4), step(feedback(C_CSD_travel,tf_travel));
title('\fontsize{18}\fontname{Times}Control Signal for Travel');
grid on;ylabel("\fontsize{16}\fontname{Times}Amplitude");xlabel("\fontsize{14}\fontname{Times}Time");
subplot(3,2,6), impulse(feedback(C_CSD_travel*tf_travel, 1));
title('\fontsize{18}\fontname{Times}Impulse input for Travel');
grid on;ylabel("\fontsize{16}\fontname{Times}Amplitude");xlabel("\fontsize{14}\fontname{Times}Time");

%% PID Controller Designing
% Time Space (Elevation)
C_p_elevation = pid(8.3008);
C_i_elevation = pid(0,1.1772);
C_pi_elevation = pid(1.16, 1.09);
C_pd_elevation = pid(1.45e+03, 0, 250);
C_pdf_elevation = pid(986, 0, 229, 0.00022);
C_pid_elevation = pid(163, 66.3, 100);
C_pidf_elevation = pid(168, 103, 67, 0.0153);
% Frequence Space (Travel)
C_p_travel = pid(-0.043);
C_i_travel = pid(0, -0.0252);
C_pi_travel = pid(-0.0468, -3.8e-05);
C_pd_travel = pid(-0.151, 0, -0.82);
C_pdf_travel = pid(-0.198, 0, -3, 1.25);
C_pid_travel = pid( -0.197 -0.00288, -1.96);
C_pidf_travel = pid(-0.177, -0.00364, -1.77, 1.14);
% P
ans_travel_p = feedback(series(ss_travel, C_p_travel), 1);
cs_travel_p = feedback(C_p_travel,ss_travel);
indis_travel_p = feedback(ss_travel, C_p_travel);
outdis_travel_p = feedback(1, series(ss_travel, C_p_travel));
ans_elevation_p = feedback(series(ss_elevation, C_p_elevation), 1);
cs_elevation_p = feedback(C_p_elevation,ss_elevation);
indis_elevation_p = feedback(ss_elevation, C_p_elevation);
outdis_elevation_p = feedback(1, series(ss_elevation, C_p_elevation));
% I
ans_travel_i = feedback(series(ss_travel, C_i_travel), 1);
cs_travel_i = feedback(C_i_travel,ss_travel);
indis_travel_i = feedback(ss_travel, C_i_travel);
outdis_travel_i = feedback(1, series(ss_travel, C_i_travel));
ans_elevation_i = feedback(series(ss_elevation, C_i_elevation), 1);
cs_elevation_i = feedback(C_i_elevation,ss_elevation);
indis_elevation_i = feedback(ss_elevation, C_i_elevation);
outdis_elevation_i = feedback(1, series(ss_elevation, C_i_elevation));
% PD
ans_travel_pd = feedback(series(ss_travel, C_pd_travel), 1);
cs_travel_pd = feedback(C_pd_travel,ss_travel);
indis_travel_pd = feedback(ss_travel, C_pd_travel);
outdis_travel_pd = feedback(1, series(ss_travel, C_pd_travel));
ans_elevation_pd = feedback(series(ss_elevation, C_pd_elevation), 1);
cs_elevation_pd = feedback(C_pd_elevation,ss_elevation);
indis_elevation_pd = feedback(ss_elevation, C_pd_elevation);
outdis_elevation_pd = feedback(1, series(ss_elevation, C_pd_elevation));
% PI
ans_travel_pi = feedback(series(ss_travel, C_pi_travel), 1);
cs_travel_pi = feedback(C_pi_travel,ss_travel);
indis_travel_pi = feedback(ss_travel, C_pi_travel);
outdis_travel_pi = feedback(1, series(ss_travel, C_pi_travel));
ans_elevation_pi = feedback(series(ss_elevation, C_pi_elevation), 1);
cs_elevation_pi = feedback(C_pi_elevation,ss_elevation);
indis_elevation_pi = feedback(ss_elevation, C_pi_elevation);
outdis_elevation_pi = feedback(1, series(ss_elevation, C_pi_elevation));
% PDF
ans_travel_pdf = feedback(series(ss_travel, C_pdf_travel), 1);
cs_travel_pdf = feedback(C_pdf_travel,ss_travel);
indis_travel_pdf = feedback(ss_travel, C_pdf_travel);
outdis_travel_pdf = feedback(1, series(ss_travel, C_pdf_travel));
ans_elevation_pdf = feedback(series(ss_elevation, C_pdf_elevation), 1);
cs_elevation_pdf = feedback(C_pdf_elevation,ss_elevation);
indis_elevation_pdf = feedback(ss_elevation, C_pdf_elevation);
outdis_elevation_pdf = feedback(1, series(ss_elevation, C_pdf_elevation));
% PID
ans_travel_pid = feedback(series(ss_travel, C_pid_travel), 1);
cs_travel_pid = feedback(C_pid_travel,ss_travel);
indis_travel_pid = feedback(ss_travel, C_pid_travel);
outdis_travel_pid = feedback(1, series(ss_travel, C_pid_travel));
ans_elevation_pid = feedback(series(ss_elevation, C_pid_elevation), 1);
cs_elevation_pid = feedback(C_pid_elevation,ss_elevation);
indis_elevation_pid = feedback(ss_elevation, C_pid_elevation);
outdis_elevation_pid = feedback(1, series(ss_elevation, C_pid_elevation));
% PIDF
ans_travel_pidf = feedback(series(ss_travel, C_pidf_travel), 1);
ans_elevation_pidf = feedback(series(ss_elevation, C_pidf_elevation), 1);
cs_travel_pidf = feedback(C_pidf_travel,ss_travel);
cs_elevation_pidf = feedback(C_pidf_elevation,ss_elevation);
indis_elevation_pidf = feedback(ss_elevation, C_pidf_elevation);
outdis_elevation_pidf = feedback(1, series(ss_elevation, C_pidf_elevation));
indis_travel_pidf = feedback(ss_travel, C_pidf_travel);
outdis_travel_pidf = feedback(1, series(ss_travel, C_pidf_travel));
%% Elevation P
figure,
subplot(1,2,1)
hold on
step(ans_elevation_p, 30);
impulse(ans_elevation_p);
title("\fontsize{16}\fontname{Times}Elevation Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_elevation_p)
hold off

subplot(1,2,2)
hold on
step(cs_elevation_p);
impulse(cs_elevation_p);
title("\fontsize{16}\fontname{Times}Elevation Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_elevation_p)
hold off

figure,
subplot(2,2,1); grid on
step(indis_elevation_p);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_elevation_p);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_elevation_p)
subplot(2,2,3); grid on
step(outdis_elevation_p);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_elevation_p);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_elevation_p)
%% Elevation I
figure,
subplot(1,2,1)
hold on
step(ans_elevation_i, 100);
impulse(ans_elevation_i);
title("\fontsize{16}\fontname{Times}Elevation Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_elevation_i)
hold off

subplot(1,2,2)
hold on
step(cs_elevation_i);
impulse(cs_elevation_i);
title("\fontsize{16}\fontname{Times}Elevation Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_elevation_i)
hold off

figure,
subplot(2,2,1); grid on
step(indis_elevation_i);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_elevation_i);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_elevation_i)
subplot(2,2,3); grid on
step(outdis_elevation_i);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_elevation_i);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_elevation_i)
%% Elevation PI
figure,
subplot(1,2,1)
hold on
step(ans_elevation_pi, 70);
impulse(ans_elevation_pi);
title("\fontsize{16}\fontname{Times}Elevation Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_elevation_pi)
hold off

subplot(1,2,2)
hold on
step(cs_elevation_pi);
impulse(cs_elevation_pi);
title("\fontsize{16}\fontname{Times}Elevation Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_elevation_pi)
hold off

figure,
subplot(2,2,1); grid on
step(indis_elevation_pi);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_elevation_pi);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_elevation_pi)
subplot(2,2,3); grid on
step(outdis_elevation_pi);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_elevation_pi);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_elevation_pi)
%% Elevation PD
figure,
hold on
step(ans_elevation_pd);
impulse(ans_elevation_pd);
title("\fontsize{16}\fontname{Times}Elevation Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_elevation_pd)
hold off

figure,
subplot(2,2,1); grid on
step(indis_elevation_pd);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_elevation_pd);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_elevation_pd)
subplot(2,2,3); grid on
step(outdis_elevation_pd);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_elevation_pd);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_elevation_pd)
%% Elevation PDF
figure,
subplot(1,2,1)
hold on
step(ans_elevation_pdf);
impulse(ans_elevation_pdf);
title("\fontsize{16}\fontname{Times}Elevation Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_elevation_pdf)
hold off

subplot(1,2,2)
hold on
step(cs_elevation_pdf);
impulse(cs_elevation_pdf);
title("\fontsize{16}\fontname{Times}Elevation Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_elevation_pdf)
hold off

figure,
subplot(2,2,1); grid on
step(indis_elevation_pdf);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_elevation_pdf);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_elevation_pdf)
subplot(2,2,3); grid on
step(outdis_elevation_pdf);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_elevation_pdf);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_elevation_pdf)
%% Elevation PID
figure,
hold on
step(ans_elevation_pid);
impulse(ans_elevation_pid);
title("\fontsize{16}\fontname{Times}Elevation Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_elevation_pid)
hold off

figure,
subplot(2,2,1); grid on
step(indis_elevation_pid);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_elevation_pid);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_elevation_pid)
subplot(2,2,3); grid on
step(outdis_elevation_pid);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_elevation_pid);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_elevation_pid)
%% Elevation PIDF
figure,
subplot(1,2,1)
hold on
step(ans_elevation_pidf);
impulse(ans_elevation_pidf);
title("\fontsize{16}\fontname{Times}Elevation Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_elevation_pidf)
hold off

subplot(1,2,2)
hold on
step(cs_elevation_pidf);
impulse(cs_elevation_pidf);
title("\fontsize{16}\fontname{Times}Elevation Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_elevation_pidf)
hold off

figure,
subplot(2,2,1); grid on
step(indis_elevation_pidf);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_elevation_pidf);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_elevation_pidf)
subplot(2,2,3); grid on
step(outdis_elevation_pidf);
title("\fontsize{16}\fontname{Times}Elevation Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_elevation_pidf);
title("\fontsize{16}\fontname{Times}Elevation Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_elevation_pidf)
%% Travel P
figure,
subplot(1,2,1)
hold on
step(ans_travel_p);
impulse(ans_travel_p);
title("\fontsize{16}\fontname{Times}Travel Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_travel_p)
hold off

subplot(1,2,2)
hold on
step(cs_travel_p);
impulse(cs_travel_p);
title("\fontsize{16}\fontname{Times}Travel Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_travel_p)
hold off

figure,
subplot(2,2,1); grid on
step(indis_travel_p);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_travel_p);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_travel_p)
subplot(2,2,3); grid on
step(outdis_travel_p);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_travel_p);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_travel_p)
%% travel I
figure,
subplot(1,2,1)
hold on
step(ans_travel_i, 100);
impulse(ans_travel_i);
title("\fontsize{16}\fontname{Times}Elevation Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_travel_i)
hold off

subplot(1,2,2)
hold on
step(cs_travel_i);
impulse(cs_travel_i);
title("\fontsize{16}\fontname{Times}Travel Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_travel_i)
hold off

figure,
subplot(2,2,1); grid on
step(indis_travel_i);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_travel_i);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_travel_i)
subplot(2,2,3); grid on
step(outdis_travel_i);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_travel_i);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_travel_i)
%% Travel PI
figure,
subplot(1,2,1)
hold on
step(ans_travel_pi, 70);
impulse(ans_travel_pi);
title("\fontsize{16}\fontname{Times}Travel Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_travel_pi)
hold off

subplot(1,2,2)
hold on
step(cs_travel_pi);
impulse(cs_travel_pi);
title("\fontsize{16}\fontname{Times}Travel Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_travel_pi)
hold off

figure,
subplot(2,2,1); grid on
step(indis_travel_pi);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_travel_pi);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_travel_pi)
subplot(2,2,3); grid on
step(outdis_travel_pi);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_travel_pi);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_travel_pi)
%% Travel PD
figure,
subplot(1,2,1)
hold on
step(ans_travel_pd);
impulse(ans_travel_pd);
title("\fontsize{16}\fontname{Times}Travel Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_travel_pd)
hold off

subplot(1,2,2)
hold on
step(cs_travel_pdf);
impulse(cs_travel_pdf);
title("\fontsize{16}\fontname{Times}Travel Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_travel_pdf)
hold off

figure,
subplot(2,2,1); grid on
step(indis_travel_pd);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_travel_pd);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_travel_pd)
subplot(2,2,3); grid on
step(outdis_travel_pd);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_travel_pd);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_travel_pd)
%% Elevation PDF
figure,
subplot(1,2,1)
hold on
step(ans_travel_pdf);
impulse(ans_travel_pdf);
title("\fontsize{16}\fontname{Times}Travel Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_travel_pdf)
hold off

subplot(1,2,2)
hold on
step(cs_travel_pdf);
impulse(cs_travel_pdf);
title("\fontsize{16}\fontname{Times}Travel Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_travel_pdf)
hold off

figure,
subplot(2,2,1); grid on
step(indis_travel_pdf);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_travel_pdf);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_travel_pdf)
subplot(2,2,3); grid on
step(outdis_travel_pdf);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_travel_pdf);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_travel_pdf)
%% Travel PID
figure,
subplot(1,2,1)
hold on
step(ans_travel_pid);
impulse(ans_travel_pid);
title("\fontsize{16}\fontname{Times}Travel Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_travel_pid)
hold off

subplot(1,2,2)
hold on
step(cs_travel_pid);
impulse(cs_travel_pid);
title("\fontsize{16}\fontname{Times}Travel Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_travel_pid)
hold off

figure,
subplot(2,2,1); grid on
step(indis_travel_pid);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_travel_pid);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_travel_pid)
subplot(2,2,3); grid on
step(outdis_travel_pid);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_travel_pid);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_travel_pid)
%% Travel PIDF
figure,
subplot(1,2,1)
hold on
step(ans_travel_pidf);
impulse(ans_travel_pidf);
title("\fontsize{16}\fontname{Times}Travel Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(ans_travel_pidf)
hold off

subplot(1,2,2)
hold on
step(cs_travel_pidf);
impulse(cs_travel_pidf);
title("\fontsize{16}\fontname{Times}Travel Control Signal Response to Step and Impulse Inputs");
grid on
legend('\fontsize{16}\fontname{Times}Step', '\fontsize{16}\fontname{Times}Impulse');
ylabel("\fontsize{18}\fontname{Times}Amplitude");
xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(cs_travel_pidf)
hold off

figure,
subplot(2,2,1); grid on
step(indis_travel_pidf);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,2); grid on
impulse(indis_travel_pidf);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone to Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(indis_travel_pidf)
subplot(2,2,3); grid on
step(outdis_travel_pidf);
title("\fontsize{16}\fontname{Times}Travel Input Disturbance Rejection Respone to Step Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
subplot(2,2,4); grid on
impulse(outdis_travel_pidf);
title("\fontsize{16}\fontname{Times}Travel Output Disturbance Rejection Respone Impulse Input");
ylabel("\fontsize{18}\fontname{Times}Amplitude"); xlabel("\fontsize{18}\fontname{Times}Time");
stepinfo(outdis_travel_pidf)
%% 
toc