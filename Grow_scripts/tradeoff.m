function tradeoff(C,S,n_valores_diag)
number_mortalities=3;
replicas=10;
mm=linspace(0,0.005,number_mortalities);
diagM=linspace(-5,-.1,n_valores_diag);
%
for i=1:length(mm)
    i
stab =zeros(n_valores_diag*replicas,1);
LL   =zeros(n_valores_diag*replicas,1);
    for j=1:length(diagM)
        for k=1:replicas
            M=comm_matrix(C,S,diagM(j));
            [F,X,R,L] = is_feasible(M,1,mm(i));
            stab( (j-1)*replicas+k)=-max(real(eig(M)));
            LL(   (j-1)*replicas+k)=L;
        end   
    end
        plot((stab+1),(LL+1),'o')
        xlabel('Stability')
        ylabel('Factibility')
        lsline
        hold on
end
    hold off
end
