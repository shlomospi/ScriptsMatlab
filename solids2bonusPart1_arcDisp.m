clc
n=2; %# of beams
R = 1; %Radius

P = 1; %Force
E = 1;
I = 1;
x = sym('x');
delta_arc=pi*R^3*P/(2*E*I); %displacement for arc case
for e=0:60 %temp loop for a test
alpha=pi/n; %for every n, compute alpha
d=2*R*cos((pi-alpha)/2); %for every n, compute beam length

if n == 0
    disp('choose a whole positive number')
    break
end
if n<0
    disp('choose a whole positive number')
    break
end
if mod(n,1)~=0
    disp('choose a whole positive number')
    break
elseif n == 1
    D=0;
    disp('no bending no delta')
    
elseif n == 2 %is here so we could use a loop for even n that always adds the previous torque
  
    M=P*x*sin((pi-alpha)/2);
    Mtag=x*sin((pi-alpha)/2);
    fun = M*Mtag;
    Delta_B=2/(E*I)*(int(fun,0,d));
    D=double(Delta_B);
    plot(n,D/delta_arc,'b--o');
    
elseif mod(n,2)==0
    
    dUdP=2/(E*I)*int((P*x*sin((pi-alpha)/2))*(x*sin((pi-alpha)/2)),0,d); % skipping the first step, we add the energy, torgue and dm/dp as starting values
    M=P*d*sin((pi-alpha)/2);
    Mtag=d*sin((pi-alpha)/2);

    for i=2:(n/2) %for every rod on the left but the first and the last:
        
    beta=(pi-((i*2)-1)*alpha)/2; %find beta of i
    Mi=M+P*x*sin(beta);      %find the torque function for this rod using torque in the end of the previous rod.
    Mtagi=Mtag+x*sin(beta);       %find dm/dp
    fun = Mi*Mtagi;
    dUdP=dUdP+2/(E*I)*int(fun,0,d);    %compute energy derivative and add it
    M=M+P*d*sin(beta);           %find the torgue at the end of the rod for next iteration.
    Mtag=Mtag+d*sin(beta);
    end
    D=double(dUdP);
    
elseif mod(n,2)==1
    
    dUdP=2/(E*I)*int((P*x*sin((pi-alpha)/2))*(x*sin((pi-alpha)/2)),0,d); % skipping the first step, we add the energy, torgue and dm/dp as starting values
    M=P*d*sin((pi-alpha)/2);
    Mtag=d*sin((pi-alpha)/2);
    for i=2:((n-1)/2)           %for every rod on the left but the first and the last:
        beta=(pi-((i*2)-1)*alpha)/2;
        Mi=M+P*x*sin(beta);      %find the torque function for this rod using torque in the end of the previous rod.
        Mtagi=Mtag+x*sin(beta);       %find dm/dp
        fun = Mi*Mtagi;
        dUdP=dUdP+2/(E*I)*int(fun,0,d);    %compute energy derivative and add it
        M=M+P*d*sin(beta);           %find the torgue at the end of the rod for next iteration.
        Mtag=Mtag+d*sin(beta);
    end
    dUdP=dUdP+2/(E*I)*(M*Mtag*d/2); %add the middle rod energy change
    D=double(dUdP);
    
end
plot(n,D/delta_arc,'b--o')
n=n+1;%%%%%%%%%%%part of test only
hold on
end