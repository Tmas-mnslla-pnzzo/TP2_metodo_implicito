graphics_toolkit("gnuplot")

%constantes:

T_ext=40;
T_con=20;

l1=3;
l2=5;
l3=5;
l4=2;
%W / (m K)
k1=1;
k2=0.400;
k3=0.030;
k4=0.250;


%https://www.studocu.com/es/document/universidade-de-vigo/termodinamica-e-transmision-de-calor/tablas-de-calor/17847611
%https://ingemecanica.com/tutoriales/pesos.html#otros
%https://www.arquimaster.com.ar/descargas/articulo410.pdf

c1=840*274.15 ; % J/K.Kg;
c2=1330; % J/°C.Kg;
c3=1673.6; % J/°C.Kg; 
c4=1000*274.15; % J/K.Kg; 

%kg / m3
P1=1800;
P2=1000;
P3=20;
P4=600;

%parametros:

t_final=50;

Dt=0.5;
Dx=0.5;
h=1;

beta1=k1/(c1*P1);
beta2=k2/(c2*P2);
beta3=k3/(c3*P3);
beta4=k4/(c4*P4);

n1=l1/Dx;
n2=l2/Dx;
n3=l3/Dx;
n4=l4/Dx;

gamma1=((Dt/(Dx*Dx))*beta1);
gamma2=((Dt/(Dx*Dx))*beta2);
gamma3=((Dt/(Dx*Dx))*beta3);
gamma4=((Dt/(Dx*Dx))*beta4);

%funciones

function mi=matriz_i(ni,a)
  N=ni+2;
  a1=a;
  a2=(1-2*a);
  a3=a;
  M0=[1 zeros(1,N-1)];
  MN=[zeros(1,N-3) a (1+a) (-2*a)];
  M1=diag(a1*ones(1,N-1),1)+diag(a2*ones(1,N))  + diag(a3*ones(1,N-1),-1);
  M1=[M0;M1];
  M1(N+1,:)=MN;
  M1(2,:)=[];
  mi=M1;
endfunction

function map=mapa(n,gam,U,T_con,t_final,Dt)
  M=matriz_i(n-1,gam);
  F=[ones(size(M)(2)-1,1)*T_con];
  V=[U(1); F];
  for i=1:1:(1/Dt)*t_final
    F=[U(i); F];
    F=M*F;
    V=[V F];
    F(1,:)=[];
  end
  map=V';
endfunction

M1=matriz_i(n1-1,0.5);
U1=[T_ext;ones(size(M1)(2)-1,1)*T_con];
V=[U1];
for j=1:1:(1/Dt)*t_final
  U1=M1*U1;
  V=[V U1];
end

Z1=V';
X1=[0: Dx/h: l1/h];
Y1=[0: Dt: t_final];

U2=Z1(:,n1+1);
Z2=mapa(n2,0.5,U2,T_con,t_final,Dt);
X2=[0: Dx/h: l2/h];

U3=Z2(:,n2+1);
Z3=mapa(n3,0.5,U3,T_con,t_final,Dt);
X3=[0: Dx/h: l3/h];

U4=Z3(:,n3+1);
Z4=mapa(n4,0.5,U4,T_con,t_final,Dt);
X4=[0: Dx/h: l4/h];

Z_f1=horzcat(Z1,Z2);
X_f1=horzcat(X1,X2+l1/h);

Z_f2=horzcat(Z_f1,Z3);
X_f2=horzcat(X_f1,X3+l2/h+l1/h);

Z_f3=horzcat(Z_f2,Z4);
X_f3=horzcat(X_f2,X4+l3/h+l2/h+l1/h);

surf(X_f3,Y1,Z_f3)
