%% SSY285 Linear Control System Design:  Assignment - 2
%% Mechanical Group 11 - Manikanta Venkatesh, Chintalapudi Adhitya Reddy and Fikri Farhan Witjaksono

clc;
clear all;
close all;
%% Task (a)
% Ke = 0.1;
% Kt = 0.1;
% J1 = 1e-5;
% J2 = 4e-5;
% Bf = 2e-3;
% D2 = 2;

R = sym('R','positive');
syms D1 D2 Bf Kt Ke J1 J2 

A = [ 0 0 0 1 0 ;
     0 0 0 0  1;
     0 (D2/Bf) (-D2/Bf) 0 0 ;
     (-D1/J1) (D1/J1) 0 (-(Kt*Ke)/(R*J1)) 0;
     (D1/J2) ((-D1-D2)/J2) (D2/J2) 0 0];

B= [0 0;
    0 0;
    0 (1/Bf);
   (Kt/(R*J1)) 0 ;
     0 0];
 
C1 = [0 1 0 0 0;
     0 0 0 0 1];
 
 C2 = [ 0 0 0 ((-Ke)/R) 0;
     0 (D2/Bf) (-D2/Bf) 0 0];
 
D = zeros(2);

W_r = [ B A*B (A)^2*B (A)^3*B (A)^4*B];

Rank_w_r = rank(W_r)
%% Controllability check
if Rank_w_r == size(A,2)
    controllability_a1 = 1
    disp('system is controllable')
else controllability_a1 = 0
    disp('System is not controllable')
end 
%% Observeability of the system

%case 1
W_o_1 = [C1;C1*A;C1*(A^2);C1*(A^3);C1*(A^4)];

%case 2
W_o_2 = [C2;C2*A;C2*(A^2);C2*(A^3);C2*(A^4)];

Rank_w_o_1 = rank(W_o_1)
Rank_w_o_2 = rank(W_o_2)

%obsereveability check
if Rank_w_o_1 == size(A,2)
    observeability_a1 = 1
    disp('system with C1 inputs is observable')
else observeability_a1 = 0
    disp('system with C1 inputs  is  not observable')
end
    
if Rank_w_o_2 == size(A,2)
    observeability_a2 = 1
    disp ('system with C2 inputs  is observable')
else observeability_a2 = 0
    disp ('system with C2 inputs  is  not observable')
end 


%% Task (b)

Ke = 0.1;
Kt = 0.1;
J1 = 1e-5;
J2 = 4e-5;
Bf = 2e-3;
D2 = 2;
D1 = 20;

A = [ 0 0 0 1 0 ;
     0 0 0 0  1;
     0 (D2/Bf) (-D2/Bf) 0 0 ;
     (-D1/J1) (D1/J1) 0 (-(Kt*Ke)/(R*J1)) 0;
     (D1/J2) ((-D1-D2)/J2) (D2/J2) 0 0];

eig_A =eig(A)

 Matrix= [0*eye(5)-A;C2] 
 
 Matrix_1 = [0*eye(5)-A;C1] 
 Rank_M =rank(Matrix)

 %subsystem 2
if Rank_M == size(A,2)
    disp('The system is detectable')
else disp('The system is not detectable')
end 

%subsystem 1
 Rank_M1 =rank(Matrix_1)

if Rank_M1 == size(A,2)
    disp('The system is detectable')
else disp('The system is not detectable')
end 

%% Task (c)

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

% Define the 1st Controllability and Observability Matrix
% and count for number of uncontrollable states

O1 = obsv(AP,C1);
K1 = ctrb(AP,B1);

con_obsv_1 = cond(O1);
con_ctrb_1 = cond(K1);

if rank(K1)== length(AP)
    controlability =1 
    disp('system is controllable')
else controlability = 0
    disp('system is not controllable')
end
    
if rank(O1)== length(AP)
    controlability_c1 = 1
     disp('system is controllable')
else controlability_c1 = 0
    disp('system is not controllable')
end 



% Define the 2nd Controllability and Observability Matrix
% and count for number of uncontrollable states

B2= [0 0;
    0 0;
    0 (1/Bf);
   (Kt/(R*J1)) 0 ;
     0 0];

C2 = [ 0 0 0 ((-Ke)/R) 0;
       0 (D2/Bf) (-D2/Bf) 0 0];
   
D2 = [(1/R) 0;
    0 (1/Bf)];
 
O2 = obsv(AP,C2);
K2 = ctrb(AP,B2);

con_obsv_2 = cond(O2);
con_ctrb_2 = cond(K2);

if rank(K2)== length(AP)
    controlability_c2 =1 
    disp('system is controllable')
else controlability_c2 = 0
    disp('system is not controllable')
end
    
if rank(O2)== length(AP)
    controlability_c2 = 1
     disp('system is controllable')
else controlability_c2 = 0
    disp('system is not controllable')
end 


%% Task (d),(e)

% Direct method using zero order hold equation 
Ts=0.001
R = 1;
Ke = 0.1;
Kt = 0.1;
J1 = 1e-5;
J2 = 4e-5;
Bf = 2e-3;
D1 = 20;
D2 = 2;

A =[ 0 0 0 1 0 ;
     0 0 0 0  1;
     0 (D2/Bf) (-D2/Bf) 0 0 ;
     (-D1/J1) (D1/J1) 0 (-(Kt*Ke)/(R*J1)) 0;
     (D1/J2) ((-D1-D2)/J2) (D2/J2) 0 0];
 
 B= [0 0;
    0 0;
    0 (1/Bf);
   (Kt/(R*J1)) 0 ;
     0 0];
 
Ad = expm(A*Ts)
Bd = double(inv(A)*(Ad-eye(5))*B)

Ts1 = linspace(0,0.001,1000);
for i = 0:0.00002:0.001
    Bd1 = expm(A*i)*B;
end

%Alternate method using MATLAB functions
S = ss(AP,B1,C1,D);
p = c2d(S,0.001)
Add=p.A
Bdd=p.B
%% Task (f)
eigenvalues_Ad = eig(Ad)

% Check if it is minimal order (realization)
% minimal order if it is both controllable and observable
% A non-minimal system has either uncontrollable or
% unobservable modes or both.


O1d = obsv(Ad,C1);
K1d = ctrb(Ad,B1);

rankk_1d = rank(K1d)
ranko_1d = rank(O1d)

if rank(O1d)== length(Ad)
    Observability_O1d =1 
    disp('system is observable')
else Observability_O1d = 0
    disp('system is not observable')
end

if rank(K1d)== length(Ad)
    controlability_K1d =1 
    disp('system is controllable')
else controlability_K1d = 0
    disp('system is not controllable')
end

O2d = obsv(Ad,C2);
K2d = ctrb(Ad,B2);

rankk_2d =  rank(K2d)
ranko_2d = rank(O2d)

if rank(K2d)== length(Ad)
    controlability_K2d = 1
     disp('system is controllable')
else controlability_K2d = 0
    disp('system is not controllable')
end 

if rank(O2d)== length(Ad)
   Observability_O2d =1 
    disp('system is observable')
else Observability_O2d = 0
    disp('system is not observable')
end