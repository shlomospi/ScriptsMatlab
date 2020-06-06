clc;
clear all;
close all;
k= 4;  %number of layers
tk = 0.25e-3; %m
t = tk*4;
E1 = (138e3); %(MPa)
E2 = 9e3; %(MPa)
G12 =6.9e3; %(MPa)
nu12 =0.3;
nu21 = (nu12*E2)/E1;
alphak = [0.88e-6 ;31e-6; 0];  %1/C
Q = zeros(3);
Q(1,1) = E1/(1-nu12*nu21);
Q(1,2) =(nu12*E2)/(1-nu12*nu21);
Q(2,1) = (nu12*E2)/(1-nu12*nu21);
Q(2,2) =  E2/(1-nu12*nu21);
Q(3,3) = G12;
% [ 45, -45, 90,0]
theta =[pi/4, -pi/4, pi/2, 0];

 %Q for the entire plate
 Qtot = 0;
 alphatot=0;
 for K =1:k
     theta_k = theta(1,K);
     c = cos(theta_k);s= sin(theta_k);
     Tsig = [c^2 s^2 2*s*c ; s^2 c^2 -2*s*c;...
         -s*c s*c c^2-s^2];
     Teps = [c^2 s^2 s*c ; s^2 c^2 -s*c;...
    -2*s*c 2*s*c c^2-s^2];
     Qtot =  Qtot+(1/t)*(tk.*inv(Tsig)*Q*Teps)
     alphatot=alphatot +(tk.*(inv(Tsig)*Q)*alphak)
 end
 alphatot = (1/t).*inv(Qtot)*alphatot;
 
 %strain and stress

%  sigt = [sigxx; sigyy; sigxy];
 epsmK=zeros(3,12);epstrK=zeros(3,4);
 sigmK = zeros(3,12);sigtrK=zeros(3,4);
 i=1;
 for K= 1:k
     theta_k = theta(1,K);
      c = cos(theta_k);s= sin(theta_k);
      Tsig = [c^2 s^2 2*s*c ; s^2 c^2 -2*s*c;...
          -s*c s*c c^2-s^2];
      Teps = [c^2 s^2 s*c ; s^2 c^2 -s*c;...
          -2*s*c 2*s*c c^2-s^2];
      epsmK(:,i:i+2) = Teps*inv(Qtot);
      sigmK(:,i:i+2)=Q*Teps*inv(Qtot);
      i=i+3;
      epstrK(:,K) = Teps*alphatot;
      sigtrK(:,K) = Q*(Teps*alphatot-alphak);
 end
 % 3 - DT = 0, sig = {0,0,50}
 epsT=zeros(3,4);epsbarT=zeros(3,4);
sigT=zeros(3,4);sigbarT=zeros(3,4);
sig0=[0;0;50];DT=0;
i=1;
x1=2*[k,k];
x2 = [1,1]
figure(1);
for K= 1:k
    theta_k = theta(1,K);
    c = cos(theta_k);s= sin(theta_k);
    Tsig = [c^2 s^2 2*s*c ; s^2 c^2 -2*s*c;...
        -s*c s*c c^2-s^2];
    Teps = [c^2 s^2 s*c ; s^2 c^2 -s*c;...
        -2*s*c 2*s*c c^2-s^2];
    epsT(:,K)= epsmK(:,i:i+2) *sig0+DT*epstrK(:,K);
    epsbarT(:,K) = inv(Teps)*epsT(:,K);
    sigT (:,K)=  sigmK(:,i:i+2)*sig0+DT*sigtrK(:,K);
    sigbarT(:,K) =inv(Tsig)*sigT(:,K);
    i=i+3;
    
    %plot
  
    sig11 = [0,sigbarT(1,K)];
    sig22 = [0,sigbarT(2,K)];
    sig12 = [0,sigbarT(3,K)];
    
    plot(sig11,x1, 'b', 'Linewidth',20);hold on;
    plot(sig22,x1, 'r', 'Linewidth', 15); hold on;
    plot(sig12,x1, 'g', 'Linewidth', 10); hold on;
    legend('\sigma_{11}', '\sigma_{22}', '\sigma_{12}')
    plot(sig11,x2, 'b', 'Linewidth',20);hold on;
    plot(sig22,x2, 'r', 'Linewidth', 15); hold on;
    plot(sig12,x2, 'g', 'Linewidth', 10); hold on;
    x1=x1-1;
    x2=x2+1;

end
grid on;title('Laminate Stress Analysis')
ylim([0,9]);xlim([-80,130]);xlabel('Stress [MPa]'); 
ylabel('Lamina Layer Position');

