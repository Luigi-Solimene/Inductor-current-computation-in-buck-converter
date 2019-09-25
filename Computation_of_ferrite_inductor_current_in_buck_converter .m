clear all
close all
clc

%%%% Initialization of the L(i) parameters 
Lh=9.1e-5;
Ll=3.2968e-7;
I30=0.2588;
I70=0.4549;
%%%%

%%%% Initialization of the simulation parameters
U_i=15;         % Input voltage of the converter  
U_d=0.5;        % Diode threshold voltage
R_on=0.3;       % On resistance of the converter
R_off=0.05;     % Off resistance of the converter
U_o=3.3;        % Load voltage
f_s=465e3;      % Switching frequency
T_s=1/f_s;      % Switching period  
I_o=0.45;       % Average load current 
N_points=5000;  % Number of points of the simulation arrays
%%%%

%%%% Definition of the L(i) curve
G30=(0.7*Lh-Ll)/(Lh-Ll);
G70=(0.3*Lh-Ll)/(Lh-Ll);
IL=(I70*cot(pi*G30)-I30*cot(pi*G70))/(cot(pi*G30)-cot(pi*G70));
sigma=(cot(pi*G30)-cot(pi*G70))/(I30-I70);
di=0.001;
I=0:di:10;
L=Ll+(Lh-Ll)/2*(1-2/pi*atan(sigma*(I-IL)));
Lfp=0.5*(max(L)+min(L));
phi_tab=0.5*(cumsum(L(1:end-1)*di)+cumsum(L(2:end)*di));
I1_tab=I(1:end-1)+di/2;
%%%%

%%%% Execution of the fixed point algorithm
delta=(U_o+U_d+R_off*I_o)/(U_i+U_d-(R_on-R_off)*I_o);
t=0:T_s/(N_points-1):T_s;
u_l=(U_i-U_o)*max(-sign(t-delta*T_s),0)-(U_o+U_d)*max(sign(t-delta*T_s),0);
u_l=u_l-mean(u_l)+(U_i-U_o)*delta-(U_o+U_d)*(1-delta);
n0=round((N_points-1)*delta);       
tau1=Lfp/R_on;
tau2=Lfp/R_off;
h1=exp(-t(1:n0)/tau1);
h2=exp(-(t(n0+1:end)-delta*T_s)/tau2);
dt=T_s/(N_points-1);
aus1=u_l(1:n0).*(fliplr(h1));
a0=dt*sum(aus1);
aus2=u_l(n0+1:end).*(fliplr(h2));
a1=dt*sum(aus2);
y0=(a0*exp(-(1-delta)*T_s/tau2)+a1)/(1-exp(-delta*T_s/tau1-(1-delta)*T_s/tau2));
for ii=1:n0
    phi0(ii)=dt*sum(u_l(1:ii).*fliplr(h1(1:ii)));
end
phi0(1:n0)=phi0(1:n0)+y0*exp(-t(1:n0)/tau1);
for ii=(n0+1):length(t)
    phi0(ii)=dt*sum(u_l(n0+1:ii).*fliplr(h2(1:ii-n0)));
end
phi0(n0+1:end)=phi0(n0+1:end)+phi0(n0)*exp(-(t(n0+1:end)-delta*T_s)/tau2);
phi=phi0;
i=phi/Lfp;
ir=i-interp1(I1_tab,phi_tab,i,'linear','extrap')/Lfp;
for j=1:100 %iterazioni Punto Fisso
   aus1=R_on*ir(1:n0).*(fliplr(h1));
   a0=-dt*sum(aus1);
   aus2=R_off*ir(n0+1:end).*(fliplr(h2));
   a1=-dt*sum(aus2); 
   y0=(a0*exp(-(1-delta)*T_s/tau2)+a1)/(1-exp(-delta*T_s/tau1-(1-delta)*T_s/tau2));
   for ii=1:n0
    phi(ii)=-dt*sum(R_on*ir(1:ii).*fliplr(h1(1:ii)));
   end
   phi(1:n0)=phi(1:n0)+y0*exp(-t(1:n0)/tau1);
   for ii=n0+1:length(t);
     phi(ii)=-dt*sum(R_off*ir(n0+1:ii).*fliplr(h2(1:ii-n0)));
   end
   phi(n0+1:end)=phi(n0+1:end)+phi(n0)*exp(-(t(n0+1:end)-delta*T_s)/tau2);
   phi=phi+phi0;
   i=phi/Lfp+ir;
   irold=ir;
   ir=i-interp1(I1_tab,phi_tab,i,'linear','extrap')/Lfp;
   max(abs(ir-irold));
   mean(i);
end
%%%%

t_u=t(1:end-1)+diff(t)/2;
u_l=diff(phi)./diff(t);

figure(1)
plot(t_u,u_l)
axis square

figure(2)
plot(t,i)
axis square