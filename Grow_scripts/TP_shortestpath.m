function TPsp=TP_shortestpath(A)
%LGN
%10 JULIO 2015
S=size(A,1);
%Calculamos las basales sin considerar canibalismo porque si hay basales
%can?bales quedan fuera en la pre-selecci?n de las matrices de adyascencia
basal=find(sum(A)==0);
n=length(basal);
aux=zeros(n,S);
%Deleting self effects (Loops) para evitar warning de "biograph"
Na=A;
V=diag(Na);
M=diag(V,0);
Mat=Na-M;

for i=1:n
[dist, path, pred] = shortestpath(biograph(Mat),basal(i));
aux(i,:)=dist;
end

TPsp=min(aux,[],1);
TPsp(basal)=0;
TPsp=TPsp';