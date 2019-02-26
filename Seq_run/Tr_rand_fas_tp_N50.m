addpath('../src/');

con=[0.1,0.2,0.3];
cp=0:1/35:0.35;
fas=1;
df={@d_c};
ind=[];

k=0;

for i=1:length(cp)
    for ii=1:length(con)
        for iii=1:length(fas)
            k=k+1;
            ind=[ind;cp(i),0,0,con(ii),fas(iii),k];
        end
    end
end


for i=1:length(cp)
    for iii=1:length(fas)
            k=k+1;
            ind=[ind;cp(i),0,0,con(ii),fas(iii),k];
        end
end


for i=1:length(cp)
        for iii=1:length(fas)
            k=k+1;
            ind=[ind;cp(i),0,0,con(ii),fas(iii),k];
        end
end


mps=1;
tp=0;
cd ../save/N_50/mps1
N=50;
addpath('../../..');
parfor i=1:length(ind(:,1))
    [Tr,b]=Tr_matrix(ind(i,4),N);
    com_scale(Tr,ind(i,1),1,ind(i,2),1,0,1,0,1,@d_c,ind(i,3),1,0.2,0,ind(i,5),1,1,0.2,0.2,tp,mps,1,0.0,500,ind(i,6),false);
end

mps=0;
tp=0;
cd ../mps0

parfor i=1:length(ind(:,1))
    [Tr,b]=Tr_matrix(ind(i,4),N);
    com_scale(Tr,ind(i,1),1,ind(i,2),1,0,1,0,1,@d_c,ind(i,3),1,0.2,0,ind(i,5),1,1,0.2,0.2,tp,mps,1,0.0,500,ind(i,6),false);  
end

mps=1;
tp=1;
cd ../tp1

parfor i=1:length(ind(:,1))
    [Tr,b]=Tr_matrix(ind(i,4),N);
    com_scale(Tr,ind(i,1),1,ind(i,2),1,0,1,0,1,@d_c,ind(i,3),1,0.2,0,ind(i,5),1,1,0.2,0.2,tp,mps,1,0.0,500,ind(i,6),false);
end
