clear; 
clc;
%set value
Z12=0+0.12i; 
Zsh=0+0.10i; 
y12=1/Z12; 
Y11=y12; 
Y12=-y12; 
Y21=Y12; 
Y22=Y11;
%Conductance and Susceptance Values
B11=imag(Y11); 
G11=real(Y11); 
B12=imag(Y12); 
G12=real(Y12); 
B21=imag(Y21); 
G21=real(Y21); 
B22=imag(Y22); 
G22=real(Y22); 
Gsh=real(1/Zsh); 
Bsh=imag(1/Zsh); 
V1=1.05+0i;
%Initialization
V22=1.0; 
Vsh=1.0; 
L2=0; 
Lsh=0;
syms V2 l vsh lsh;
%The power flow equations of Bus 2
P2=V2*(V1*(G21*cos(l-0)+B21*sin(l-0))+V2*(G22*cos(l-l)+B22*sin(l-l))); 
Q2=V2*(V1*(G21*sin(l-0)-B21*cos(l-0))+V2*(G22*sin(l-l)-B22*cos(l-l)));
%The power flow equations of STATCOM
Psh=(V2^2)*Gsh-V2*vsh*(Gsh*cos(l-lsh)+Bsh*sin(l-lsh)); %Psh 
Qsh=-(V2^2)*Bsh-V2*vsh*(Gsh*sin(l-lsh)-Bsh*cos(l-lsh)); %Qsh 
%Operating constraints of STATCOM
PE=(vsh^2)*Gsh-V2*vsh*(Gsh*cos(l-lsh)-Bsh*sin(l-lsh)); %PE=RE 
F=V2-1.03;
%The power mismatch equations
P=-0.8-P2-Psh; 
Q=-0.25-Q2-Qsh;
%The Jacobian Matrix
J=jacobian([PE;F;P;Q],[lsh,vsh,l,V2]); 
n=0;
put=[];
DELTA=ones(4,1);%Correcting the values of variables 
while max(abs(DELTA))>1e-6
%Substitute the updated values into the variables
pe=double(subs(PE,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh])); 
f=double(subs(F,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh])); 
p=double(subs(P,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh])); 
q=double(subs(Q,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh])); 
d=[-pe;-f;-p;-q];
VJ=double(subs(J,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh])); 
DELTA=VJ\d;
n=n+1;
%put(n)=max(abs(d(3)-d(4)));
put(n)=max(abs(DELTA));
V22=V22+DELTA(4);%Bus 2 voltage magnitude 
L2=L2+DELTA(3); %Bus 2 voltage angle
Vsh=Vsh+DELTA(2); %STATCOM voltage angle 
Lsh=Lsh+DELTA(1); %STATCOM voltage magnitude
fprintf('Iteration %i  Voltage Magnitude at bus 2:%2.5e  Voltage Angle at bus 2:%2.5e\n',n,V22,L2/pi*180);
end
%ToleranceV=sqrt((fbus21+deltaf)^2+(ebus21+deltae)^2)-sqrt(fbus21^2+ebus21^2);
%ToleranceV=double(abs(ToleranceV));
%A1(k)=ToleranceV;
%plot
figure('NumberTitle','off','Name','-Max⁡|V_i^(k+1)-V_i^k|'); 
plot(-put,'Linewidth', 4, 'MarkerSize', 4);
xlabel('Number of iterations'); 
ylabel('-Max⁡|V_i^(k+1)-V_i^k|'); 
title('Convergence characteristics');
%calculate P1 Q1 P2 Q2 Psh Qsh
P1=V1*(V1*(G11)+V22*(G12*cos(L2)-B12*sin(L2)));
Q1=V1*(V1*(-B11)+V22*(-G12*sin(L2)-B12*cos(L2)));
P2=double(subs(P2,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh]));
Q2=double(subs(Q2,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh]));
Psh=double(subs(Psh,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh])); 
Qsh=double(subs(Qsh,[V2,l,vsh,lsh],[V22,L2,Vsh,Lsh]));
fprintf('Bus1 active power output: %4i Bus1 reactive power output: %4i\n',P1,Q1);
fprintf('STATCOM branch active power injection: %4i STATCOM branch reactive power injection: %4i\n',Psh,Qsh);
fprintf('Bus2 active power output: %4i Bus2 reactive power output: %4i\n',P2,Q2);