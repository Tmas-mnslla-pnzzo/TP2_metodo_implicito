graphics_toolkit qt

%constantes:
T_ext=40;
T_con=20;

%cm
l1=3;
l2=12;
l3=5;
l4=1;

%W / (m K)
k1=1;
k2=0.400;
k3=0.030;
k4=0.250;

% J/K.Kg;
c1=840; 
c2=790; 
c3=1600; 
c4=1000;  

%kg / m3
P1=1800;
P2=1000;
P3=20;
P4=600;

%parametros:
h=2;
t_final=20000;
r1=1;
r2=5;

Dt=(h/r1);
Dx=(h/r2);

beta1=k1/(c1*P1)*100*100;
beta2=k2/(c2*P2)*100*100;
beta3=k3/(c3*P3)*100*100;
beta4=k4/(c4*P4)*100*100;

n1=(l1/Dx)*h;
n2=(l2/Dx)*h;
n3=(l3/Dx)*h;
n4=(l4/Dx)*h;

gamma1=((Dt/(Dx*Dx))*beta1);
gamma2=((Dt/(Dx*Dx))*beta2);
gamma3=((Dt/(Dx*Dx))*beta3);
gamma4=((Dt/(Dx*Dx))*beta4);

gamma1<0.5;
gamma2<0.5;
gamma3<0.5;
gamma4<0.5;

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
  V=V';
  V=[U V];
  V(:,end)=[];
  map=V;
endfunction

M1=matriz_i(n1-1,gamma1);
U1=[T_ext;ones(size(M1)(2)-1,1)*T_con];
V=[U1];

for j=1:1:(1/Dt)*t_final
  U1=M1*U1;
  V=[V U1];
end

Gf=6;

Z1=V'; X1=[0: Dx/h: l1]; Y=[0: Dt*(2**Gf): t_final];

U2=Z1(:,n1+1); Z2=mapa(n2,gamma2,U2,T_con,t_final,Dt); X2=[0: Dx/h: l2];
 
U3=Z2(:,n2+1); Z3=mapa(n3,gamma3,U3,T_con,t_final,Dt); X3=[0: Dx/h: l3];

U4=Z3(:,n3+1); Z4=mapa(n4,gamma4,U4,T_con,t_final,Dt); X4=[0: Dx/h: l4];

Z_f1=horzcat(Z1,Z2); X_f1=horzcat(X1,X2+l1);

Z_f2=horzcat(Z_f1,Z3); X_f2=horzcat(X_f1,X3+l2+l1);

Z=horzcat(Z_f2,Z4); X=horzcat(X_f2,X4+l3+l2+l1);

Zf=Z'(:,end);

for gf=1:1:Gf
  Z(2:2:end,:)=[];
end

subplot (2, 2, 3);
  plot(X,Zf); title('Temperatura en el tiempo final. U(x,t_f)' ); xlabel ("Distancia [cm]"); ylabel ("Temperatura [°C]");
subplot (2, 2, 1);
  waterfall(X,Y,Z); title('Temperatura U(x,t)' ); xlabel ("Distancia [cm]"); ylabel ("Tiempo [s]"); zlabel("Temperatura [°C]");

