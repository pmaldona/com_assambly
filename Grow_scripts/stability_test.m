function out=stability_test(M)
%assess local stability of matrix M
le=length(M);
e=eig(M);
s=sum(sign(real(e)));
if s==-le
    out=1;
else
    out=0;
end
end
  
    

