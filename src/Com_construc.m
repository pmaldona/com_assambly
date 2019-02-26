function [A,B,R,R_tr,L,L_tr,bas_e,con,con_tr,cp,mu,am,cm,an,me,me_tr,gama,gama_tr,var,var_tr,eta,eta_tr,ict,Nv,dist,x,r] = Com_construc(Tr,cpi,fcp,mui,fmu,cmi,fcm,ami,fam,dnt,ani,fan,rit,mit,fas,mnti,nb,sg,sgd,tsp,msp,max_r,min_mort,it,bc)

N=length(Tr);


R=zeros(it,N);
R_tr=zeros(it,1);
bas_e=zeros(it,N);
con=zeros(it,N);
con_tr=zeros(it,1);
cp=zeros(it,N);
mu=zeros(it,N);
am=zeros(it,N);
an=zeros(it,N);
cm=zeros(it,N);
me=zeros(it,N);
me_tr=zeros(it,1);
gama=zeros(it,N);
gama_tr=zeros(it,1);
var=zeros(it,N);
var_tr=zeros(it,1);
eta=zeros(it,N);
eta_tr=zeros(it,1);
ict=zeros(5,it,N);
L=zeros(it,N);
L_tr=zeros(it,1);
dist=NaN(it,N,N);
Nv=NaN(it,N,N);
x=NaN(it,N,N);
r=NaN(it,N,N);




for k=1:it
    
    A=NT_Community(Tr,cpi,fcp,mui,fmu,cmi,fcm,ami,fam,dnt,ani,fan,rit,mit,fas,mnti,nb,sg,sgd,tsp,msp,max_r,min_mort,false,0.02,0.025,bc);
    
    R(k,1)=A.res;
    R_tr(k)=A.res_tr;
    L(k,1)=A.L;
    L_tr(k)=A.L_tr;
    bas_e(k,1)=A.bas;
    con(k,1)=A.con;
    con_tr(k)=A.con_tr;
    cp(k,1)=A.cp;
    mu(k,1)=A.mu;
    am(k,1)=A.am;
    cm(k,1)=A.cm;
    an(k,1)=A.an;
    me(k,1)=A.me;
    me_tr(k)=A.me_tr;
    gama(k,1)=A.gama;
    gama_tr(k)=A.gama_tr;
    var(k,1)=A.var;
    var_tr(k)=A.var_tr;
    eta(k,1)=A.eta;
    eta_tr(k)=A.eta_tr;
    dist(k,1,1:length(A.k))=A.k;
    Nv(k,1,1:length(A.Nv))=A.Nv;
    x(k,1,1:length(A.Nv))=A.X;
    r(k,1,1:length(A.Nv))=A.R;
    
    
    B=A.NT_SubCommunity_init_Links(2);
    R(k,2)=B.res;
    L(k,2)=B.L;
    bas_e(k,2)=B.bas;
    con(k,2)=B.con;
    cp(k,2)=B.cp;
    mu(k,2)=B.mu;
    am(k,2)=B.am;
    cm(k,2)=B.cm;
    an(k,2)=B.an;
    me(k,2)=B.me;
    gama(k,2)=B.gama;
    var(k,2)=B.var;
    eta(k,2)=B.eta;
    ict(:,k,1)=B.ict;
    dist(k,2,1:length(B.k))=B.k;
    Nv(k,2,1:length(B.Nv))=B.Nv;
    x(k,2,1:length(A.Nv))=A.X;
    r(k,2,1:length(A.Nv))=A.R;
    
    
    S(k)=2;
    
    for j=3:N
        B=B.NT_SubCommunity_add_Links(A,1);
        
        R(k,j)=B.res;
        L(k,j)=B.L;
        bas_e(k,j)=B.bas;
        con(k,j)=B.con;
        cp(k,j)=B.cp;
        mu(k,j)=B.mu;
        am(k,j)=B.am;
        cm(k,j)=B.cm;
        an(k,j)=B.an;
        me(k,j)=B.me;
        gama(k,j)=B.gama;
        var(k,j)=B.var;
        eta(k,j)=B.eta;
        ict(:,k,j-1)=B.ict;
        dist(k,j,1:length(B.k))=B.k;
        Nv(k,j,1:length(B.Nv))=B.Nv;
        x(k,j,1:length(A.Nv))=A.X;
        r(k,j,1:length(A.Nv))=A.R;
    end
end

end