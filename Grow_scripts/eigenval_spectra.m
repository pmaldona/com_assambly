function [stab,LL]=eigenval_spectra(C,S,n,min_mort)
max_r=1;
%min_mort=0.5;
figure
contF=0;
%
r_hist=zeros(C*n,1);
x_hist=zeros(C*n,1);
stab=zeros(n,1);
LL=zeros(n,1);
%
for i=1:n
M=comm_matrix(C,S,-1);
%M=comm_matrix_rnd_trophic(C,S);
%TPsp=TP_shortestpath(M)

[F,X,R,L] = is_feasible(M,max_r,min_mort);
stab(i)=-max(real(eig(M)));
LL(i)=L;
if F
    plot(eig(M),'bo')
  
%   subplot(3,1,2)=hist(R)
%   subplot(3,1,3)=hist(X)
    contF=contF+1;
else
  	plot(eig(M),'ro') 
%   subplot(3,1,2)=hist(R)
%   subplot(3,1,3)=hist(X)
end
r_hist((i-1)*S+1:i*S,1)=R;
x_hist((i-1)*S+1:i*S,1)=X;
hold on
end
hold off
% figure
% hist(r_hist)
% title('R')
% figure
% hist(x_hist)
% title('X')
figure
plot(stab,LL,'ko')
xlabel('Stability (-max(real(eig))')
ylabel('Feasibility (min(X_eq))')
LL
contF
end
