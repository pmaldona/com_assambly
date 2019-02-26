max_r=1;
min_mort=0;
A=[-1 0 -1 0 -1;
    0 -1 -1 -1 0;
    1 1 -1 0 -1;
    0 1 0 -1 0;
    1 0 1 0 -1];
A=A.*unifrnd(.75,1.25,5,5).*(1-eye(5))-eye(5);
[F,X,R,L] = is_feasible(A,max_r,min_mort);
T1=[1 5 3 2 4];
T=[2 4 1 3 5];
v_stab=zeros(length(T),1);
v_fact=zeros(length(T),1);
for i=1:5
    max_r=R(sort(T(1:i)));
    min_mort=max_r;

    M=A(sort(T(1:i)),:);
    N=M(:,sort(T(1:i)))
    %[F,X,R,L] = is_feasible(N,max_r,min_mort);
    eq=-N\R(sort(T(1:i)));
    LL=min(eq);
    J=diag(eq)*N ;
    v_stab(i)=-max(real(eig(J)));
    v_fact(i)=LL;
    %sort(real(eig(N)));
    R(sort(T(1:i)))
    figure
    if LL>0
        plot(eig(J),'bo')
    else
        plot(eig(J),'ro')
    end
end
% v_stab;
% figure
% subplot(2,1,1)=plot(1:length(T),v_stab,'bo');
% figure
% subplot(2,1,2)=plot(1:length(T),v_fact,'ro');