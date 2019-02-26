function M=comm_matrix_rnd_trophic(C,S)
%random trophic matrix
diagM=-.1;
A=rand(S)<C;
A=triu(A,1);
M=A.*(2*(rand(S)<=0.5)-1);
M=M-M';
M=M.*rand(S);
M=M.*(1-eye(S))+eye(S)*diagM;
end