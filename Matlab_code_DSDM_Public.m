%% Matlab_code_DSDM_for_decomposing_period-4_bifurcation_signal
% Figure 1 is the frequency response of the period-4_bifurcation_signal.
% Figure 2-4 plot the period-4 bifurcation signal, typical vibration component, and period-4 bifurcation component.
% Figure 5-7 illustrate the decomposition and reconstruction errors for the period-4 bifurcation signal.
% For the paper titled ¡®Digital synchronous decomposition and period-N
% bifurcation size identification in dynamic systems: application to a milling process¡¯
clc;
clear
close all
   
M = 400;% Number of data samples collected per revolution
dt = 1/100/M;% The period of the period-4 bifurction signal is 1/100 s.
tt = dt:dt:0.04;
x1 = 0.7*sin(50*pi*tt)+0.6*sin(50*2*pi*tt)+0.45*sin(150*pi*tt)+ 0.4*sin(250*pi*tt)+0.25*sin(350*pi*tt) ; % Period-4 bifurcation component
x2 =  0.8*sin(2*50*2*pi*tt)+ 0.55*sin(4*50*2*pi*tt)+ 0.25*sin(6*50*2*pi*tt);                             % Typical vibration component without constant term
% x2 =  0.8*sin(2*50*2*pi*tt)+ 0.55*sin(4*50*2*pi*tt)+ 0.25*sin(6*50*2*pi*tt)+1;                         % Typical vibration component with constant term
x = x1+x2;
fs= 1/dt;

s_n = size(tt);
n=0:s_n(2)-1;
figure(1)
x_x = x;
y=fft(x_x(1:s_n(2)),s_n(2)); % Fast Fourier transform
mag=abs(y); 
f=n*fs/s_n(2); %Frequency
figure(1)
plot(f(1:s_n(2)/2),mag(1:s_n(2)/2)*2/s_n(2)); 

axis([0,s_n(2)/4,0,1.5])
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title('Frequency response')

%% 
Dis = zeros(M,4);% 4 periods
M_Dis = zeros(M,4);
bifucation_component = zeros(M,4);% Decomposed bifucation component
T_component = zeros(M,4);% Original typical vibration component
B_component = zeros(M,4);% Original bifucation component
tic;
for i = 1:M

 Dis(i,:) = x( i: M : 0+3*M+i );
 M_Dis(i,1) = 1/4*(Dis(i,1)+Dis(i,2)+Dis(i,3)+Dis(i,4));
 M_Dis(i,2) = 1/4*(Dis(i,1)+Dis(i,2)+Dis(i,3)+Dis(i,4));
 M_Dis(i,3) = 1/4*(Dis(i,1)+Dis(i,2)+Dis(i,3)+Dis(i,4));
 M_Dis(i,4) = 1/4*(Dis(i,1)+Dis(i,2)+Dis(i,3)+Dis(i,4));

 bifucation_component(i,:) = Dis(i,:) - M_Dis(i,:); % Decomposed period-4 bifucation component

end
toc;
box on

figure(2) 
plot(x(1:M*4),'.b')
xlabel('Sampling point {\itn}')
ylabel('Amplidute')
title('Period-4 bifurcation signal')
axis([0,4*M,-4,4])
box on
figure (3)
xlabel('Sampling point {\itn}')
ylabel('Amplidute')
axis([0,4*M,-4,4])
for  i = 1:M
    hold on
    plot(i:M:3*M+i,M_Dis(i,:),'.b')
end
title('Typical vibration component ')
box on

figure(4)
xlabel('Sampling point {\itn}')
ylabel('Amplidute')
axis([0,4*M,-2,2])
for  i = 1:M
    
    hold on
    plot(i:M:3*M+i,bifucation_component(i,:),'.b')
    T_component(i,:) = x2( i: M : 0+3*M+i );%  
    B_component(i,:) = x1( i: M : 3*M+i );
end
box on
title('Period-4 bifurcation component ')

figure (5)
xlabel('Sampling point {\itn}')
ylabel('Amplidute')

for  i = 1:M
    
    hold on
    plot(i:M:3*M+i,T_component(i,:)-M_Dis(i,:),'b.')
 
end
title('Decomposition errors of typical vibration component')
box on
axis([0,4*M,-8*10^-15,8*10^-15])
figure (6)
xlabel('Sampling point {\itn}')
ylabel('Amplidute')
% title('Error')
for  i = 1:M
    
    hold on
    plot(i:M:3*M+i,B_component(i,:)-bifucation_component(i,:),'b.')
    
end
title('Decomposition errors of period-4 bifurcation component')
box on
axis([0,4*M,-8*10^-15,8*10^-15])
figure (7)
xlabel('Sampling point {\itn}')
ylabel('Amplidute')

for  i = 1:M
    
    hold on
    plot(i:M:3*M+i,B_component(i,:)+T_component(i,:)-(M_Dis(i,:)+bifucation_component(i,:)),'.b')
    
end
title('Reconstruction errors for period-4 bifurcation signal')
box on
axis([0,4*M,-3*10^-16,3*10^-16])
%%
