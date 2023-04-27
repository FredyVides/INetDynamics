% function x=SpSolver(A,y,L,tol,sl,delta)
%%
% Example:
% t=0:2/34:12;
% y=sin(pi*t)+.5*sin(3*pi*t)+.25*sin(7*pi*t)+.15*cos(7*pi*t)-.1*sin(9*pi*t);
% N=30;
% M=N-1;
% L=8;
% [H,H0,H1]=HankelSMData(y,N,1);
% D=@(j,n)(n==j);
% H00=H0(:,1:M).';
% H11=H1(:,1:M).';
% A=diag(ones(1,N-1),1);
% A(N,:)=SpSolver(H00,H11(:,N),L,1e-5).';
% A1=SpSolver(H00,H11,1+D(1:N,N)*(L-1),1e-5).';
% y0=H(:,1);
% z0=y0;
% for k=1:(length(y)-1), y0=[y0 A*y0(:,k)];z0=[z0 A1*z0(:,k)];end
% subplot(211);hold on;plot(t,y,'b'),plot(t,y0(1,:),'r.-');hold off;
% subplot(212);hold on;plot(t,y,'b'),plot(t,z0(1,:),'r.-');hold off;
% s=0:1/100:1;
% S1=exp(2*pi*i*s);
% figure;
% subplot(121);hold on;plot(S1,'b'),plot(eig(full(A)),'r.','markersize',9);hold off;
% axis tight;axis square;
% subplot(122);hold on;plot(S1,'b'),plot(eig(full(A1)),'r.','markersize',9);hold off;
% axis tight;axis square;
% z=y+1e-4*randn(1,length(y));
% [H,H0,H1]=HankelSMData(z,N,1);
% H00=H0(:,1:M).';
% H11=H1(:,1:M).';
% A=diag(ones(1,N-1),1);
% A(N,:)=SpSolver(H00,H11(:,N),L,1e-2).';
% A1=SpSolver(H00,H11,1+D(1:N,N)*(L-1),1e-2).';
% y0=H(:,1);
% z0=y0;
% for k=1:(length(y)-1), y0=[y0 A*y0(:,k)];z0=[z0 A1*z0(:,k)];end
% figure;
% subplot(211);hold on;plot(t,z,'b'),plot(t,y0(1,:),'r.-');hold off;
% subplot(212);hold on;plot(t,z,'b'),plot(t,z0(1,:),'r.-');hold off;
% s=0:1/100:1;
% S1=exp(2*pi*i*s);
% figure;
% subplot(121);hold on;plot(S1,'b'),plot(eig(full(A)),'r.','markersize',9);hold off;
% axis tight;axis square;
% subplot(122);hold on;plot(S1,'b'),plot(eig(full(A1)),'r.','markersize',9);hold off;
% axis tight;axis square;
% tau = 1e-5; mu=1/2; MaxIt = 1e5; tol= 1e-12;
% sigma=1e-3;
% [u,s]=svd([H0;H1(end,:)]',0);
% rk=sum(diag(s)>sigma);
% u=u(:,1:rk);
% H00=(H0*u).';
% H11=(H1*u).';
% tic, for k=1:N, C(k,:)=DouglasRachford(H00,H11(:,k),5e-3,tau,mu,MaxIt,tol);end;toc
% tic,A1=SpSolver(H00,H11,N,5e-3).';toc
% y0=H(:,1);
% z0=y0;
% for k=1:(length(z)-1), y0=[y0 A1*y0(:,k)];z0=[z0 C*z0(:,k)];end
% figure;
% subplot(211);hold on;plot(t,z,'b'),plot(t,y0(1,:),'r.-');hold off;
% subplot(212);hold on;plot(t,z,'b'),plot(t,z0(1,:),'r.-');hold off;
% 
%%
function X=SpSolver(A,Y,L,tol,sl,delta)
if nargin<=4 
	sl=1;delta=tol;
elseif nargin<=5 
	delta=tol;
end
N=size(Y,2);
%if length(L)==1, L=L*ones(1,N);end
[u,s,~]=svd(A,0);
rk=sum(diag(s)>tol);
A=u(:,1:rk)'*A;
Y=u(:,1:rk)'*Y;
for k=1:N
	X(:,k)=SpLSSolver(A,Y(:,k),L,tol,sl,delta);
end
end
