%% SSY285 Linear Control System Design: Assignment - 3
%% Group 11 - Fikri Farhan Witjaksono,Chintalapudi Adhitya Reddy and Manikanta Venkatesh 

clc
clear all
close all

%% Assignment 2 values

R = 1;
Ke = 0.1;
Kt = 0.1;
J1 = 1e-5;
J2 = 4e-5;
Bf = 2e-3;
D1 = 20;
D2 = 2;

AP = [ 0 0 0 1 0 ;
     0 0 0 0  1;
     0 (D2/Bf) (-D2/Bf) 0 0 ;
     (-D1/J1) (D1/J1) 0 (-(Kt*Ke)/(R*J1)) 0;
     (D1/J2) ((-D1-D2)/J2) (D2/J2) 0 0];
 
B1= [0 0;
    0 0;
    0 (1/Bf);
   (Kt/(R*J1)) 0 ;
     0 0];
 
C1 = [0 1 0 0 0;
     0 0 0 0 1];
 
D1 = zeros(2);

S = ss(AP,B1,C1,D1);
p = c2d(S,0.001);
Add=p.A;
Bdd=p.B;
Cdd=p.C;
Ddd=p.D;

%% Assignment 3
%% Task (a)
mu_1=0;
sigma_1=0.3/3;
sigma_2=0.1*1/3;

a1 = sigma_1^2;
a2 = sigma_2^2;
R1 = [a1 0;
      0 a2];

N = Bdd; % x_dot1 = Ax + Bu + Nw, 
         % x_dot2 = Ax + Bu + Bw, 
         % Hence B = N

%% Task (b)
mu_2=0;
sigma_3=0.02/3;
sigma_4=0.01/3;

a3 = sigma_3^2;
a4 = sigma_4^2;
R2 = [a3 0;
      0 a4];


%% Task (c)
h = 0.001;
Qn = R1;
Rn = R2;
Nn = zeros(2);

G_kalman = ss(Add,[Bdd Bdd],Cdd,[Ddd Ddd],h);
[KEST,K,P] = kalman(G_kalman,Qn,Rn,Nn)

% Observer eigenvalue
Ob = Add-(K*Cdd) ;
Obsv_eig = eig(Ob)

%% Task (d)

x = [1 1 1 1 1 1];
Qx = diag(x);

u = [1 1];
Qu = 5*diag(u);

A_lq = [Add [0;0;0;0;0]; 0,0,0,0,-1,1];

B_lq = [Bdd;[0,0]];

C_lq = [Cdd, [0;0]];

% LQ solution
[L_LQ,S,P] = dlqr(A_lq,B_lq,Qx,Qu);
L_P = L_LQ(:,1:5);
L_I = L_LQ(:,6);

sim ('Group_11_Assignment_3_Task_D_Simulink_model');


