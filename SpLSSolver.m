% function x=SpLSSolver(A,y,L,tol,sl,delta)
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
% A(N,:)=SpLSSolver(H00,H11(:,N),L,1e-5).';
% for k=1:N, A1(k,:)=SpLSSolver(H00,H11(:,k),1+D(k,N)*L,1e-4).';end
% y0=H(:,1);
% z0=y0;
% for k=1:(length(y)-1), y0=[y0 A*y0(:,k)];z0=[z0 A1*z0(:,k)];end
% subplot(211);plot(t,y,'b',t,y0(1,:),'r.-');
% subplot(212);plot(t,y,'b',t,z0(1,:),'r.-');
% s=0:1/100:1;
% S1=exp(2*pi*i*s);
% figure;
% subplot(121);plot(S1,'b',eig(A),'r.','markersize',9)
% axis tight;axis square;
% subplot(122);plot(S1,'b',eig(A1),'r.','markersize',9)
% axis tight;axis square;
%%
function x=SpLSSolver(A,y,L,tol,sl,delta)
	if nargin<=4 
		sl=1;delta=tol;
	elseif nargin<=5 
		delta=tol;
	end
N=size(A,2);
K=1;
Error = abs(1+tol);
w=sparse(N,1);
x=w;
[u,s,v]=svd(A,0);
rk=sum(diag(s)>tol);
u=u(:,1:rk);
s=s(1:rk,1:rk);
v=v(:,1:rk);
c=v*(s\(u'*y));
ac=abs(c);
[~,f]=sort(-ac);
L=min(max(sum(ac(f)>delta),1),L);
while K<=L && Error>tol
	ff=sort(f(1:K));
	x=w;
	[u,s,v]=svd(A(:,ff),0);
	rk=sum(diag(s)>tol);
	u=u(:,1:rk);
	s=s(1:rk,1:rk);
	v=v(:,1:rk);
	c=v*(s\((u'*y)));
	x(ff,1)=c;
	Error = norm(A*x-y);
	K=K+sl;
end
end
