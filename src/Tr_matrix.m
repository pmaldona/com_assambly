function [M,b] = Tr_matrix(con, n)
% Generates an niche model trophic matrix with desire connectance "con" 
% and "n" species, that is fully connected.
M=TheGenNicheModel(con,n,0.8);
M=M-transpose(M);
M(M==-1)=0;
G=digraph(M);
while(length(unique(conncomp(G,'Type','weak')))>1 || length(find(M))==0)
    M=TheGenNicheModel(con,n,0.8);
    M=M-transpose(M);
    M(M==-1)=0;
    G=digraph(M);
end

v=TroPos(M);
b=length(v(v==0));

end
