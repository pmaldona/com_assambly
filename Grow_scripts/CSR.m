function CSR(n)
max_r=1;
min_mort=0.0;
stab=zeros(n,1);
LL=zeros(n,1);
M_stab=zeros(n,n);
M_fact=zeros(n,n);
C=linspace(0.1,0.3,n);
S=round(linspace(10,100,n));
for x=1:n
    for y=1:n
        M=comm_matrix(C(x),S(y)) ;
        [F,X,R,L] = is_feasible(M,max_r,min_mort);
        M_stab(x,y)=-max(real(eig(M)));
        M_fact(x,y)=L;
    end
end
figure
imagesc(flipud(M_stab))
figure
imagesc(flipud(M_fact))
end
