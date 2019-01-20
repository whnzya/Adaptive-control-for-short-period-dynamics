% linearize the 2 degree model using matlab function 
% with q and w as the states, need further modification to this state

%% save the data generated from this linearization
load 'SP_linsys_TM.mat';
A=A2_TM;
B=B2_TM;
C=C2_TM;
%Cqw=eye(2);
D=D2_TM;

%lqr design for non_normalized system
P_sp=ss(A,B,C,D);
zpk(P_sp)
Q=diag(1,1);
R=0.01;
K=lqr(P_sp,Q,R);
eig(A-B*K)
V=400;

%% lqr design for modified states system
%Aqw(2,1)=Aqw(2,1)*op.States(3).x(1)/V^2;
%Aqw(2,2)=Aqw(2,2)*(op.States(3).x(1))^2/V^2;
%Bqw(1,1)=Bqw(1,1)*op.States(3).x(1)/V^2;
%Cqw=[0 1];

% P_sp=ss(Aqw,Bqw,Cqw,Dqw);
% q=[1 0];
% Q=diag(q);
% R=100;
% K=lqr(P_sp,Q,R);
% eig(Aqw-Bqw*K)
%% property of the system

%controlability of the system
co=ctrb(A,B);
m=rank(co);
X='rank of the controlability matrix=';
disp(X);
disp(m)

%observability of the system
to=obsv(A,C);
X='rank of the observability matrix=';
m=rank(to);
disp(X);
disp(m);

%% design controller for qw system

  [n,~]=size(A);
  [~,m]=size(B); %number of inputs
  [rr,~]=size(C);

% formulate the tracking system for analysis
 A_bar=[zeros(rr, rr) C ;zeros(n,rr) A];
 B_bar=[0 ;B];
 C_bar=[0 C];
 D_bar=0;
 P_bar=ss(A_bar,B_bar,C_bar,D_bar);
 
 %[Gm,Pm,Wcg,Wcp]
 dim=1000;
 Gm=zeros(dim,1); 
 Pm=zeros(dim,1);
 Wcg=zeros(dim,1);
 Wcp=zeros(dim,1);
 RT=zeros(dim,1);
 ST=zeros(dim,1);
 Overshoot=zeros(dim,1);
 Undershoot=zeros(dim,1);
 QQ=zeros(dim,1);
 RDD=zeros(dim,1); %sigma(1+L)
 RDD2=zeros(dim,1); %sigma(1+L-1)
 i=1;
%% tracking system design
for q1=1:1:dim
 %q1=0.1;
 QQ(i,1)=q1;
 q_bar=[q1 0 0];
 Q_bar=diag(q_bar);
 R_bar=1;
 K_bar=lqr(P_bar,Q_bar,R_bar);
 
 Ke_sp=K_bar(:,1);
 Kx_sp=K_bar(:,2:end);
 
 %analysis the system to find the appropriate p1
 
 %time domain analysis
 
 %define the closed loop system
   Acl=[0 C;-B*Ke_sp A-B*Kx_sp];
   Bcl=[-1;0;0];
   Ccl=[0 1 0];
   Dcl=0;
   Pcl=ss(Acl,Bcl,Ccl,Dcl);
   
 %generate the closed loop performance
   S = stepinfo(Pcl);
   ST(i,1)=S.SettlingTime;
   RT(i,1)=S.RiseTime;
   Overshoot(i,1)=S.Overshoot;
   Undershoot(i,1)=S.Undershoot;
   
 %frequency domain analysis 
 %margin(Pcl)
  [Gm(i,1),Pm(i,1),Wcg(i,1),Wcp(i,1)] = margin(Pcl);
 %define the open loop system
   Aol=A_bar;
   Bol=B_bar;
   Col=K_bar;
   Dol=0;
   L=ss(Aol,Bol,Col,Dol);
 %analysis the singular value of the return difference   
   RD=1+L;
   fun=sigma(RD);
   RDD(i,1)=min(fun);
   
   RD2=1+inv(L);
   fun=sigma(RD2);
   RDD2(i,1)=min(fun);
   i=i+1;
end
save('LQR.mat','Gm','Pm','Wcg','Wcp','RT','ST','Overshoot','Undershoot','QQ','RDD','RDD2');
%%
figure 
plot(QQ,Wcg);
title('crossover frequency change with penalty parameter');

figure
subplot(2,3,1)
plot(Wcg,Gm);
title('Gain margin change with crossover frequency');
xlabel('crossover frequency/Hz');
ylabel('Gain margin in dB')

subplot(2,3,2)
plot(Wcg,ST);
title('Settling time change with crossover frequency');
xlabel('crossover frequency/Hz');
ylabel('Settling time in s')

subplot(2,3,3)
plot(Wcg,RT);
title('Rise time change with crossover frequency');
xlabel('crossover frequency/Hz');
ylabel('Rise time in s')

subplot(2,3,4)
plot(Wcg,Overshoot);
title('Overshoot change with crossover frequency');
xlabel('crossover frequency/Hz');
ylabel('Percentage overshoot')

subplot(2,3,5)
plot(Wcg,RDD);
title('Sigma of return reference change with crossover frequency');
xlabel('crossover frequency/Hz');
ylabel('sigma(I+L)')

subplot(2,3,6)
plot(Wcg,RDD2);
title('Sigma of (I+inv(L)) change with crossover frequency');
xlabel('crossover frequency/Hz');
ylabel('sigma(I+inv(L))')
%% analysis the controller

%choose q1=200 for the result generated from the previous section;
 q1=200;
 q_bar=[q1 0 0];
 Q_bar=diag(q_bar);
 R_bar=1;
 K_bar=lqr(P_bar,Q_bar,R_bar);
 Ke_sp=K_bar(:,1);
 Kx_sp=K_bar(:,2:end);
 
   %define the closed loop system
   Acl=[0 C;-B*Ke_sp A-B*Kx_sp];
   Bcl=[-1;0;0];
   Ccl=[zeros(2,1) eye(2)];
   Dcl=0;
   Pcl=ss(Acl,Bcl,Ccl,Dcl);
   
   figure (1)
   step(Pcl)
   grid
   
   %% define the open loop system (verified)
   Aol=[0 C;zeros(2,1) A];
   Bol=[0;B];
   Col=[Ke_sp Kx_sp];
   Dol=0;
   L1=-ss(Aol,Bol,Col,Dol);
   L=-tf(L1);

% %% define the controller state space
%  Ac=0;
%  Bc=Cqw;
%  Cc=Ke_sp;
%  Dc=Kx_sp;
%  Pcc=ss(Ac,Bc,Cc,Dc)
 %%
 %Lc=-Pcc*P_sp

 %% time domain analysis
   figure (2)
   step(Pcl)
   grid
   title('step response of the system');
% open loop steady state representation
   figure (3)
   subplot(1,2,1)
   sigma(L,[],2);
   title('sigma(I+L)');
      grid
   
   subplot(1,2,2)
   sigma(L,[],3)
   title('sigma(I+inv(L))');
   grid
   %%
   T=inv(L+eye(1))*(L) %complementary sensitivity
   SS=inv(1+L);%sensitity
   figure (4)
   bode(SS,T);
   legend('S','T');
   grid
   %figure (5)
   %[Gm,Pm,Wcg,Wcp]=margin(L1)
   % margin(L1)
   %n = norm(sys,Inf)
   %grid
   figure (6)
   nyquist(L1)
  % axis([-2,0.1,-2,2])
      grid
   %% LQG design
   rho=10000;
   Q0=[10^(-7) 0;0 10^(-7)]+1/rho*(B*B');
   R0=0.1;
   Ak=A';
   Bk=C';
   Ck=B';
   Dk=D';
   Pk=ss(Ak,Bk,Ck,Dk);
   
   Kf=lqr(Pk,Q0,R0);
   eig(Ak-Bk*Kf)
   
   %% LQG controller analysis
   Aclqg=[zeros(1,3);-B*Ke_sp A-B*Kx_sp-Kf'*C];
   Bclqg=[C;Kf'*C];
   Cclqg=K_bar;
   Dclqg=zeros(1,2);
   Pclqg=ss(Aclqg,Bclqg,Cclqg,Dclqg);
   
   Alqg=A2_TM;
   Blqg=B2_TM;
   Clqg=eye(2);
   Dlqg=D2_TM;
   Plqg=ss(Alqg,Blqg,Clqg,Dlqg);
  
   Llqg=Pclqg*Plqg;
   
   figure
   subplot(1,3,1)
   sigma(Llqg,[],2);
   title('sigma(I+L)');
   
   subplot(1,3,2)
   sigma(Llqg,[],3)
   title('sigma(I+inv(L))');
   
   subplot(1,3,3)
   nyquist(-Llqg)
   title('nyquist L');
   
   %%
   T=inv(Llqg+eye(1))*(Llqg); %complementary sensitivity
   SS=inv(1+L);%sensitity
  
   figure
   bode(SS,T);
   legend('S','T');
   
   %%
   save('controller_parameter.mat','Ke_sp','Kx_sp','Kf');