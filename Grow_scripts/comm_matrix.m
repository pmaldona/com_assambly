function M=comm_matrix(C,S,diagM)
CON=0.5;
%diagM=-10;
A=TheGenNicheModel(C,S,CON);
f=triu(A+A')==2;
A(f)=0;
%
M=-A;
%M(f')=1
M=M-M';
M=M.*rand(length(M)).*1;
M=M.*(1-eye(length(M)))+eye(length(M))*diagM;
end
