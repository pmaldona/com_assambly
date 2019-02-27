classdef NT_Community
    %Class object for Ecological Community
    properties(Access = public)
        %Matrices
        adj
        com
        tr
        tr_i
        tp
        coef_C
        coef_b
        %Vectors
        
        Nv
        esp
        ict
        X
        R
        k
        k_p
        k_m
        m
        m_p
        m_m
        a
        a_p
        a_m
        
        %input Dobles
        
        rit
        mit
        fcp
        fmu
        fcm
        fam
        fan
        fas
        mnti
        sg
        sgd
        tsp
        msp
        max_r
        min_mort
        nbf
        
        %calculated doubles
        cp
        mu
        am
        cm
        an
        res
        res_tr
        con
        con_tr
        me_tr
        me
        gama_tr
        gama
        var_tr
        var
        eta_tr
        eta
        bas
        bas_tr
        L
        L_tr
        
        
        %bool
        mst
        est
        F
        bc
        
        %Probability function
        di
        
        
    end
    
    properties(Access=private)
        %Index vector of interacions
        idx
        
    end
    
    methods (Access=public)
        
        
        function C = NT_Community(Tr,cp,fcp,mu,fmu,cm,fcm,am,fam,dnt,an,fan,rit,mit,fas,mnti,nbf,sg,sgd,tsp,msp,max_r,min_mort,tr_cal,res,dr,bc)
            % Master Community constructor
            if nargin < 13
            elseif mnti<sg*sqrt(3)
                error('Mean should be grather that square tree times std');
            elseif 1<sgd*sqrt(3)
                error('std should be less that one over square tree');
            else
                C.est=false;
                C.mst=true;
                
                C.fcp=fcp;
                C.fmu=fmu;
                C.fam=fam;
                C.fcm=fcm;
                C.fan=fan;
                C.fas=fas;
                C.sgd=sgd;
                C.mnti=mnti;
                C.nbf=nbf;
                C.sg=sg;
                C.tsp=tsp;
                C.msp=msp;
                C.max_r=max_r;
                C.min_mort=min_mort;
                C.tr=Tr;
                C.adj=-C.tr+transpose(C.tr);
                C.bc=bc;
                C.mit=mit;
                
                %trophic level calculation
                C = C.Tro_Pos(C.tr);
                %species index calculation
                C.esp=1:length(C.Nv);
                C.di=dnt;
                C.rit=rit;
                
                % Trophic skeleton community matrix, resilen an feasibility
                % calculation
                n=length(C.Nv);
                
                for i = 1:n
                    d(i)=(i-1)*n+i;
                end
                
                C=C.bc_assing();
                C=C.type_matrix();
                C.res=1e5;
              
                B=zeros(n);
                ted=(1/sgd^2-1)/2;
                te=(mnti/sg^2-1)/2;
                
                if(tr_cal)
                    
                    while(~(C.res <= res+dr && C.res >= res-dr))
                        B_d=betarnd(ted*ones(n,1),ted*ones(n,1)).*(-2);
                        B_d(C.Nv~=0)=nbf*B_d(C.Nv~=0);
                        B(logical(eye(size(B))))=B_d;
                        C.coef_C=C.adj.*(betarnd(te*ones(n),te*ones(n)))*2*mnti+B;
                        C.coef_b=abs(C.adj).*(betarnd(te*ones(n),te*ones(n)))*2*mnti;
                        C.com=zeros(n);
                        for i=1:length(C.esp)
                            for ii=1:length(C.esp)
                                tm=C.tp(i,:);
                                if(ii==i)
                                    C.com(i,ii)=C.coef_C(i,ii);
                                elseif(tm(ii)==1)
                                    C.com(i,ii)=C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                                elseif(tm(ii)==-1)
                                    C.com(i,ii)=-C.coef_C(ii,i)/(msp+(1-msp)*(sum(abs(C.coef_C(ii,C.tp(ii,:)==C.tp(ii,i))))+tsp*(sum(abs(C.coef_C(ii,C.tp(ii,:)~=C.tp(ii,i)&C.tp(ii,:)~=0))))));
                                elseif(tm(ii)==4)
                                    C.com(i,ii)=fcp*C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                                end
                            end
                        end
                        C=C.prop();
                        
                    end
                else
                    B_d=betarnd(ted*ones(n,1),ted*ones(n,1)).*(-2);
                    B_d(C.Nv~=0)=nbf*B_d(C.Nv~=0);
                    B(logical(eye(size(B))))=B_d;
                    C.coef_C=C.adj.*(betarnd(te*ones(n),te*ones(n)))*2*mnti+B;
                    C.coef_b=abs(C.adj).*(betarnd(te*ones(n),te*ones(n)))*2*mnti;
                    C.com=zeros(n);
                    for i=1:length(C.esp)
                        for ii=1:length(C.esp)
                            tm=C.tp(i,:);
                            if(ii==i)
                                C.com(i,ii)=C.coef_C(i,ii);
                            elseif(tm(ii)==1)
                                C.com(i,ii)=C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                            elseif(tm(ii)==-1)
                                C.com(i,ii)=-C.coef_C(ii,i)/(msp+(1-msp)*(sum(abs(C.coef_C(ii,C.tp(ii,:)==C.tp(ii,i))))+tsp*(sum(abs(C.coef_C(ii,C.tp(ii,:)~=C.tp(ii,i)&C.tp(ii,:)~=0))))));
                            elseif(tm(ii)==4)
                                C.com(i,ii)=fcp*C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                            end
                        end
                    end
                    C=C.prop();
                end
                
                C.tr_i=C.tr;
                C.res_tr=C.res;
                C.con_tr=C.con;
                C.bas_tr=C.bas;
                C.gama_tr=C.gama;
                C.var_tr=C.var;
                C.eta_tr=C.eta;
                C.me_tr=C.me;
                C.L_tr=C.L;
                
                %Added index calculation
                C.idx=randperm(n*n);
                C.idx(ismember(C.idx,d))=[];
                C.idx(ismember(C.idx,find(C.adj)))=[];
                
                if(cp>0 || mu>0 || cm>0 || am>0 || an>0)
                    C = C.perm_gen();
                end
                
                
                pool=1:5;
                
                cp=round((n-1)*n*cp);
                mu=round((n-1)*n*mu);
                cm=round((n-1)*n*cm);
                am=round((n-1)*n*am);
                an=round((n-1)*n*an);
                
                if(cp>0 || mu>0 || cm>0 || am>0 || an>0)
                    
                    while cp>0 || mu>0 || cm>0 || am>0 || an>0
                        
                        type = pool(randi(numel(pool), 1, 1));
                        
                        if type==1 && cp>0;
                            C=C.assing_cp();
                            cp=cp-2;
                        elseif type==2 && mu>0;
                            C=C.assing_mu();
                            mu=mu-2;
                        elseif type==3 && cm>0;
                            C=C.assing_cm();
                            cm=cm-1;
                        elseif type==4 && am>0;
                            C=C.assing_am();
                            am=am-1;
                        elseif type==5 && an>0;
                            C=C.assing_an();
                            an=an-2;
                        end
                        
                    end
                end
                
                
                
               %Community matrix, feasibility an stability calculation.
                C=C.type_matrix();
                C=C.re_assiang_an();
                C=C.type_matrix();
                
                C_tmp_C=C.adj.*(betarnd(te*ones(n),te*ones(n)))*2*mnti;
                C_tmp_b=abs(C.adj).*(betarnd(te*ones(n),te*ones(n)))*2*mnti;
                C.coef_C(C.coef_C==0)=C_tmp_C(C.coef_C==0);
                C.coef_b(C.coef_b==0)=C_tmp_b(C.coef_b==0);
                for i=1:length(C.esp)
                    for ii=1:length(C.esp)
                        tm=C.tp(i,:);
                        if(ii==i)
                            C.com(i,ii)=C.coef_C(i,ii);
                        elseif(tm(ii)==2)
                            C.com(i,ii)=fmu*C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                        elseif(tm(ii)==3)
                            C.com(i,ii)=fcm*C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                        elseif(tm(ii)==4)
                            C.com(i,ii)=fcp*C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                        elseif(tm(ii)==5)
                            C.com(i,ii)=fam*C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                        elseif(Tr(i,ii)==0 && Tr(ii,i)==0)
                            if(tm(ii)==1)
                                C.com(i,ii)=fan*C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                            elseif(tm(ii)==-1)
                                C.com(i,ii)=-fan*C.coef_C(ii,i)/(msp+(1-msp)*(sum(abs(C.coef_C(ii,C.tp(ii,:)==C.tp(ii,i))))+tsp*(sum(abs(C.coef_C(ii,C.tp(ii,:)~=C.tp(ii,i)&C.tp(ii,:)~=0))))));
                            end
                        elseif(Tr(i,ii)~=0 || Tr(ii,i)~=0)
                            if(tm(ii)==1)
                                C.com(i,ii)=C.coef_b(i,ii)*C.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(C.coef_C(i,C.tp(i,:)==C.tp(i,ii))))+tsp*(sum(abs(C.coef_C(i,C.tp(i,:)~=C.tp(i,ii)&C.tp(i,:)~=0))))));
                            elseif(tm(ii)==-1)
                                C.com(i,ii)=-C.coef_C(ii,i)/(msp+(1-msp)*(sum(abs(C.coef_C(ii,C.tp(ii,:)==C.tp(ii,i))))+tsp*(sum(abs(C.coef_C(ii,C.tp(ii,:)~=C.tp(ii,i)&C.tp(ii,:)~=0))))));
                            end
                        end
                    end
                end
                
                % Asimetric coeficient assignment
                for i=1:(length(C.esp)-1)
                    for ii=(i+1):length(C.esp)
                        tm=C.tp(i,:);
                        if(tm(ii)~=4 && tm(ii)~=5)
                            C.com(ii,i)=fas*C.com(ii,i);
                        end
                    end
                end
                C=C.Tro_Pos(C.tp==-1);
                C=C.prop();
                
            end
        end
        
        function Sc=NT_SubCommunity_init_Links(C,S)
            
            if(C.mst)
                Sc=NT_Community();
                
                Sc.est=false;
                Sc.mst=false;
                
                Sc.fcp=C.fcp;
                Sc.fmu=C.fmu;
                Sc.fam=C.fam;
                Sc.fcm=C.fcm;
                Sc.fan=C.fan;
                Sc.fas=C.fas;
                Sc.sgd=C.sgd;
                Sc.mnti=C.mnti;
                Sc.sg=C.sg;
                Sc.tsp=C.tsp;
                Sc.nbf=C.nbf;
                Sc.rit=C.rit;
                Sc.ict=zeros(5,1);
                Sc.di=C.di;
                Sc.bc=C.bc;
                
                fcp=C.fcp;
                fmu=C.fmu;
                fam=C.fam;
                fcm=C.fcm;
                fan=C.fan;
                fas=C.fas;
                sgd=C.sgd;
                mnti=C.mnti;
                sg=C.sg;
                tsp=C.tsp;
                msp=C.msp;
                te=(mnti/sg^2-1)/2;
                
                %Species selection
                
                %basal are added to the pool
                pool_b=find(sum(C.tr)==0);
                %first basal added
                fs=pool_b(randi(numel(pool_b), 1, 1));
                Sc.esp(1)=fs;
                
                
                for i=1:(S-1)
                    for ii=Sc.esp
                        %trophic or basa added
                        pool_b=[pool_b,find(C.tr(ii,:)==1)];
                    end
                    
                    pool_b=unique(pool_b);
                    pool_b(find(ismember(pool_b,Sc.esp)))=[];
                    
                    if isempty(pool_b)
                        
                    else
                        Sc.esp=[Sc.esp,pool_b(randi(numel(pool_b), 1, 1))];
                    end
                end
                
                Sc.esp=sort(Sc.esp);
                Sc.Nv=C.Nv(Sc.esp);
                
                for i=1:length(Sc.esp)
                    for ii=1:length(Sc.esp)
                        Sc.adj(i,ii)=C.adj(Sc.esp(i),Sc.esp(ii));
                        Sc.tr(i,ii)=C.tr(Sc.esp(i),Sc.esp(ii));
                        Sc.tp(i,ii)=C.tp(Sc.esp(i),Sc.esp(ii));
                        Sc.coef_C(i,ii)=C.coef_C(Sc.esp(i),Sc.esp(ii));
                        Sc.coef_b(i,ii)=C.coef_b(Sc.esp(i),Sc.esp(ii));
                        
                    end
                end
                
                for i=1:length(Sc.esp)
                    for ii=1:length(Sc.esp)
                        tm=Sc.tp(i,:);
                        if(ii==i)
                            Sc.com(i,ii)=Sc.coef_C(i,ii);
                        elseif(tm(ii)==2)
                            Sc.com(i,ii)=fmu*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==3)
                            Sc.com(i,ii)=fcm*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==4)
                            Sc.com(i,ii)=fcp*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==5)
                            Sc.com(i,ii)=fam*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==1)
                            Sc.com(i,ii)=fan*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==-1)
                            Sc.com(i,ii)=-fan*Sc.coef_C(ii,i)/(msp+(1-msp)*(sum(abs(Sc.coef_C(ii,Sc.tp(ii,:)==Sc.tp(ii,i))))+tsp*(sum(abs(Sc.coef_C(ii,Sc.tp(ii,:)~=Sc.tp(ii,i)&Sc.tp(ii,:)~=0))))));
                        end
                        
                    end
                end
                
                for i=1:(length(Sc.esp)-1)
                    for ii=(i+1):length(Sc.esp)
                        tm=Sc.tp(i,:);
                        if(tm(ii)~=4 && tm(ii)~=5)
                            Sc.com(ii,i)=fas*Sc.com(ii,i);
                        end
                    end
                end
                
                X=C.tp-C.tp';
                X(X==-2)=0;
                X(X==2)=0;
                X=abs(X);
                X=X+C.tp;
                X(X==10)=5;
                X(X==6)=3;
                
                aux=[];
                for i=Sc.esp
                    aux=[aux,i];
                    for ii=setdiff(Sc.esp,aux)
                        
                        if(X(i,ii)==1)
                            Sc.ict(1)=Sc.ict(1)+1;
                        elseif(X(i,ii)==-1)
                            Sc.ict(1)=Sc.ict(1)+1;
                        elseif(X(i,ii)==2)
                            Sc.ict(2)=Sc.ict(2)+1;
                        elseif(X(i,ii)==3)
                            Sc.ict(3)=Sc.ict(3)+1;
                        elseif(X(i,ii)==4)
                            Sc.ict(4)=Sc.ict(4)+1;
                        elseif(X(i,ii)==5)
                            Sc.ict(5)=Sc.ict(5)+1;
                        end
                    end
                end
                Sc.R=C.R(Sc.esp);
                Sc=Sc.prop();
            end
        end
        
        function Sc=NT_SubCommunity_add_Links(Si,C,S)
            if(C.mst)
                Sc=NT_Community();
                
                Sc.est=false;
                Sc.mst=false;
                
                Sc.est=false;
                Sc.mst=false;
                
                
                Sc.fcp=C.fcp;
                Sc.fmu=C.fmu;
                Sc.fam=C.fam;
                Sc.fcm=C.fcm;
                Sc.fan=C.fan;
                Sc.fas=C.fas;
                Sc.sgd=C.sgd;
                Sc.mnti=C.mnti;
                Sc.sg=C.sg;
                Sc.tsp=C.tsp;
                Sc.nbf=C.nbf;
                Sc.rit=C.rit;
                Sc.di=C.di;
                Sc.bc=C.bc;
                
                fcp=C.fcp;
                fmu=C.fmu;
                fam=C.fam;
                fcm=C.fcm;
                fan=C.fan;
                fas=C.fas;
                sgd=C.sgd;
                mnti=C.mnti;
                sg=C.sg;
                tsp=C.tsp;
                msp=C.msp;
                te=(mnti/sg^2-1)/2;
                
                Sc.esp=Si.esp;
                
                add_esp=[];
                Sc.ict=zeros(5,1);
                
                pool_b=[Sc.esp,find(sum(C.tr)==0)];
                pool_b=unique(pool_b);
                
                for i=1:S
                    for ii=[Sc.esp,add_esp]
                        %Trophic or basal species added
                        pool_b=[pool_b,find(C.tr(ii,:)==1)];
                    end
                    
                    pool_b=unique(pool_b);
                    pool_b(find(ismember(pool_b,[Sc.esp,add_esp])))=[];
                    
                    if isempty(pool_b)
                        
                    else
                        add_esp=[add_esp,pool_b(randi(numel(pool_b), 1, 1))];
                    end
                end
                
                X=C.tp-C.tp';
                X(X==-2)=0;
                X(X==2)=0;
                X=abs(X);
                X=X+C.tp;
                X(X==10)=5;
                X(X==6)=3;
                
                aux=[];
                
                add_esp=sort(add_esp);
                for i=sort([Sc.esp,add_esp])
                    aux=[aux,i];
                    for ii=setdiff(add_esp,aux)
                        
                        if(X(i,ii)==1)
                            Sc.ict(1)=Sc.ict(1)+1;
                        elseif(X(i,ii)==-1)
                            Sc.ict(1)=Sc.ict(1)+1;
                        elseif(X(i,ii)==2)
                            Sc.ict(2)=Sc.ict(2)+1;
                        elseif(X(i,ii)==3)
                            Sc.ict(3)=Sc.ict(3)+1;
                        elseif(X(i,ii)==4)
                            Sc.ict(4)=Sc.ict(4)+1;
                        elseif(X(i,ii)==5)
                            Sc.ict(5)=Sc.ict(5)+1;
                        end
                    end
                end
                
                Sc.esp=sort([Sc.esp,add_esp]);
                Sc.Nv=C.Nv(Sc.esp);
                
                for i=1:length(Sc.esp)
                    for ii=1:length(Sc.esp)
                        Sc.adj(i,ii)=C.adj(Sc.esp(i),Sc.esp(ii));
                        Sc.tr(i,ii)=C.tr(Sc.esp(i),Sc.esp(ii));
                        Sc.tp(i,ii)=C.tp(Sc.esp(i),Sc.esp(ii));
                        Sc.coef_C(i,ii)=C.coef_C(Sc.esp(i),Sc.esp(ii));
                        Sc.coef_b(i,ii)=C.coef_b(Sc.esp(i),Sc.esp(ii));
                        
                    end
                end
                
                for i=1:length(Sc.esp)
                    for ii=1:length(Sc.esp)
                        tm=Sc.tp(i,:);
                        if(ii==i)
                            Sc.com(i,ii)=Sc.coef_C(i,ii);
                        elseif(tm(ii)==2)
                            Sc.com(i,ii)=fmu*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==3)
                            Sc.com(i,ii)=fcm*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==4)
                            Sc.com(i,ii)=fcp*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==5)
                            Sc.com(i,ii)=fam*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==1)
                            Sc.com(i,ii)=fan*Sc.coef_b(i,ii)*Sc.coef_C(i,ii)/(msp+(1-msp)*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)==Sc.tp(i,ii))))+tsp*(sum(abs(Sc.coef_C(i,Sc.tp(i,:)~=Sc.tp(i,ii)&Sc.tp(i,:)~=0))))));
                        elseif(tm(ii)==-1)
                            Sc.com(i,ii)=-fan*Sc.coef_C(ii,i)/(msp+(1-msp)*(sum(abs(Sc.coef_C(ii,Sc.tp(ii,:)==Sc.tp(ii,i))))+tsp*(sum(abs(Sc.coef_C(ii,Sc.tp(ii,:)~=Sc.tp(ii,i)&Sc.tp(ii,:)~=0))))));
                        end
                        
                    end
                end
                for i=1:(length(Sc.esp)-1)
                    for ii=(i+1):length(Sc.esp)
                        tm=Sc.tp(i,:);
                        if(tm(ii)~=4 && tm(ii)~=5)
                            Sc.com(ii,i)=fas*Sc.com(ii,i);
                        end
                    end
                end
                Sc.R=C.R(Sc.esp);
                Sc=Sc.prop();
            end
        end
    end
    
    
    
    methods (Access=private)
        
        function C=Tro_Pos(C,Na)
            % Leslie Garay-Narvaez & Jose D. Flores.
            % February 28, 2011
            % Algorithm S.Levine(1980)
            
            %Deleting self effects (Loops)
            V=diag(Na);
            M=diag(V,0);
            Mat=Na-M;
            
            %Vector TP
            TP=zeros(size(Mat,1),1);%Vector de ceros (vector x en algoritmo de Levine)
            
            
            ss=sum(Mat);%Suma columnas de matriz de adyacencia, entregando cantidad de
            %presas por nodo =In_degree.
            
            s=find(ss);%Identifica las especies no basales
            
            T1=Mat(:,s);%Saca la submatriz de interacciones entre
            %no basales y basales (cuadrante superior, equivale a matriz R
            %de Levine)+ la submatriz de las interacciones
            %entre especies no basales (cuadrante inferior, equivale a
            %matriz Q de Levine)
            
            D=diag(1./ss(s));%Divide la energ?a por la cantidad de presas de cada
            %especie y as? da peso a los links
            
            T2=T1*D;%Con este calculo de peso a los links dependiendo de la cantidad de
            %presas
            Q=T2(s,:)';%Se recupera la matriz Q de Levine que contiene las interacciones
            %entre las especies no basales como matriz de
            %transici?n(=adyacencia transpuesta)
            L=ones(max(size(Q)),1);%Vector de unos que viene del supuesto de Levine que
            %dice que, se espera que la posici?n tr?fica
            %promedio sea 1+ que la posici?n tr?fica de sus
            %recursos, recordar que se asume que las basales
            %tienen posici?n tr?fica=0
            y=pinv(eye(size(Q))-Q)*L;%Se calcula las posici?n tr?fica seg?n ecuaci?n 4-5 en
            %Levine
            TP(s)=y;%Se junta posici?n tr?fica para basales y no basales
            C.Nv=TP;
            
        end
        
        %         function [mcp,mmu,mcm,mam] = non_trophic_m(C)
        %         %Entrega el la medida de "no-troficidad" de la matriz de adyacencia ingresada C.adj
        %
        %         C.adj=sign(C.adj);
        %
        %         A=C.adj+transpose(C.adj);
        %
        %         A(logical(eye(size(A))))=zeros(1,length(A));
        %
        %         A=triu(A);
        %
        %         M= A==-2;
        %
        %         mcp=2*nnz(M)/length(A)^2;
        %
        %         L= A==2;
        %
        %         mmu=2*nnz(L)/length(A)^2;
        %
        %         J= A==1;
        %
        %         mcm=nnz(J)/length(A)^2;
        %
        %
        %         K= A==-1;
        %
        %         mam=nnz(K)/length(A)^2;
        %
        %         end
        
        function C = inter_porp(C)
            % calculates percentage of interactions of the community matrix
            
            if(isempty(C.tp))
                C.type_matrix(C);
            end
            C.an=0;
            C.cm=0;
            C.am=0;
            C.cp=0;
            C.mu=0;
            for i=1:length(C.esp)
                for j=C.tp(:,i)'
                    if(j==1)
                        C.an=C.an+1;
                    elseif(j==-1)
                        C.an=C.an+1;
                    elseif(j==2)
                        C.mu=C.mu+1;
                    elseif(j==3)
                        C.cm=C.cm+1;
                    elseif(j==4)
                        C.cp=C.cp+1;
                    elseif(j==5)
                        C.am=C.am+1;
                    end
                end
            end
            
            C.an=C.an/(length(C.esp))^2;
            C.mu=C.mu/(length(C.esp))^2;
            C.cm=C.cm/(length(C.esp))^2;
            C.cp=C.cp/(length(C.esp))^2;
            C.am=C.am/(length(C.esp))^2;
        end
        
        
        function C=perm_gen(C)
            %Monte-Carlo of permutations as a function of a distribution
            %that favors the position (favors means putting them in the
            %first places of the indices vector) of the community matrix 
            %respect to distribution C.dnt determined by the attribute C.Nv 
            %and the distance between the species.
            
            %tophic level normalization
            amax=max(C.Nv);
            amin=min(C.Nv);
            Ntl=(C.Nv-amin)/(amax-amin);
            psg=C.rit;
            l=length(C.idx);
            la=length(C.esp);
            
            % Function that calculates the score on each position
            function pri=por_idx(k)
                e=k-la*floor((k-1)/la);
                b=floor((k-1)/la)+1;
                %pri=C.di(Ntl(e))/integral(C.di,0,1,'ArrayValued',true)*(normpdf(abs(Ntl(e)-Ntl(b)),C.mit,psg)/(normcdf(1,C.mit,psg)-normcdf(0,C.mit,psg)));
                pri=C.di(Ntl(e))/max(C.di(Ntl))*(normpdf(abs(Ntl(e)-Ntl(b)),C.mit,psg)/normpdf(1,C.mit,psg));
            end
            
            % score vector generation
            pri=arrayfun(@por_idx,C.idx);
            
            
            %init output
            idt=[];
            %Montecarlo of indexes
            for m=1:(l-1)

                pm=max(pri/sum(pri));
                pre=pm*rand();
                prr=pri(1)/sum(pri);
                
                while(pre>prr)
                    
                    perm=randperm(length(C.idx));
                    C.idx=C.idx(perm);
                    pri=pri(perm);
                    
                    
                    pre=pm*rand();
                    prr=pri(1)/sum(pri);
                end

                idt=[idt,C.idx(1)];
                C.idx(1)=[];
                pri(1)=[];
            end
            C.idx=[idt,C.idx];
            
        end
        
        function C=assing_cp(C)
            %competitive added interacion function
            
            n=length(C.esp);
            B=zeros(n);
            idx=C.idx;
            idx(ismember(idx,find(C.adj)))=[];
            idx(ismember(idx,find(transpose(C.adj))))=[];
            
            if(length(idx)>1)
                
                cpi=[idx(1),n*(idx(1)-1)-(n^2-1)*floor((idx(1)-1)/n)+1];
                B(cpi)=-1;
                C.adj=C.adj+B;
            end
            C.idx(ismember(C.idx,find(C.adj)))=[];
            
        end
        
        function C=assing_an(C)
            %Trophic added interacion
            
            n=length(C.esp);
            B=zeros(n);
            idx=C.idx;
            idx(ismember(idx,find(C.adj)))=[];
            idx(ismember(idx,find(transpose(C.adj))))=[];
            
            if(length(idx)>1)

                cpi=idx(1);
                cpi_tr=n*(idx(1)-1)-(n^2-1)*floor((idx(1)-1)/n)+1;
                B(cpi)=-1;
                B(cpi_tr)=1;
                C.adj=C.adj+B;
            end
            C.idx(ismember(C.idx,find(C.adj)))=[];
            
        end
        
        function C=assing_mu(C)
            % mutla added interaction
            
            n=length(C.esp);
            B=zeros(n);
            idx=C.idx;
            idx(ismember(idx,find(C.adj)))=[];
            idx(ismember(idx,find(transpose(C.adj))))=[];
            
            if(length(idx)>1)

                cpi=[idx(1),n*(idx(1)-1)-(n^2-1)*floor((idx(1)-1)/n)+1];
                B(cpi)=1;
                C.adj=C.adj+B;
            end
            C.idx(ismember(C.idx,find(C.adj)))=[];
        end
        
        function C=assing_cm(C)
            %commensal added interactions
            
            B=zeros(length(C.esp));
            idx=C.idx;
            idx(ismember(idx,find(C.adj)))=[];
            idx(ismember(idx,find(transpose(C.adj))))=[];
            
            if(length(idx)>0)
                B(idx(1))=1;
                C.adj=C.adj+B;
            end
            C.idx(ismember(C.idx,find(C.adj)))=[];
        end
        
        function C=assing_am(C)
            %Ammensal added interaction
            
            B=zeros(length(C.esp));
            idx=C.idx;
            idx(ismember(idx,find(C.adj)))=[];
            idx(ismember(idx,find(transpose(C.adj))))=[];
            %verificación de indices
            
            if(length(idx)>0)
                B(idx(1))=-1;
                % se modifica la matriz
                C.adj=C.adj+B;
            end
            C.idx(ismember(C.idx,find(C.adj)))=[];
        end
        
        function C=re_assiang_an(C)
            %Trophic interation reassingment to preferve basal species
            
            function ji=jd(x,y)
                ji = 1 - sum(x & y)/sum(x | y);
            end
            
            
            ver=true;
            C.tr(C.tp==-1)=1;
            B=C.tr;
            C.bas=length(find(sum(C.tr)==0));
            C.bas_tr=length(find(sum(C.tr_i)==0));
            %length(find(sum(C.tr)==0))
            if(C.bas==C.bas_tr)
                ver=false;
            end
            
            while(ver)
                
                % possible basal candidates that no modify the exisiting
                % ones
                pbc=find(sum(C.tr)~=0 & sum(C.tr_i)==0);
                
                n=sum(C.tr(:,pbc));
                n=min(n(n~=0));
                
                % Basal candidates with minimal preys to change
                pbc=find(sum(C.tr)==n & sum(C.tr_i)==0);
         
                %sampled basal candidate
                can=datasample(pbc,1);
                
                % for loop over preys to remplace
                A=C.tr;
                
                % Possible assignment places where didn't exsist previos
                % ones.
                A_tp=C.tp;
                
                % preys to remplace without modifying trophic skeleton
                % matrix
                j=find(C.tr(:,can))';
                for i=j
                    %Replacement candidates that do not contain i as a prey
                    
                    r_can=find(~A_tp(i,:) & sum(A)~=0);
                    r_can(r_can==can)=[];
                    r_can(r_can==i)=[];
                    if(isempty(r_can))
                      % if there are no replacement candidates, the 
                      % interaction is eliminated
                        A(i,can)=0;
                        continue;
                    end
                    
                    %vector of distance of similarity of the candidates
                    sim=zeros(length(r_can),1);
                    l=length(r_can);
                    %vector of similarity distance calculation 
                    for ii=1:l
                        
                        if(ii==can)
                            sim(ii)=1;
                        else
                            sim(ii)=jd(A(:,can),A(:,ii));
                        end
                        
                    end
                    %vector of tophic level distance calculation
                    r_tn=abs(C.Nv(r_can)-C.Nv(can)+1)/max(abs(C.Nv(r_can)-C.Nv(can)+1));
                    % finally  election vector is the product of similarity
                    % vector by the tophic distance vector.
                    r_tn=r_tn.*sim;
                    % choosed index is index the minimun of electon vector
                    [scr,ind]=min(r_tn);
                    
                    %trophic link remplacement
                    A(i,r_can(ind))=1;
                    A(i,can)=0;
                    
                end
                C.tr=A;
                C.bas=length(find(sum(A)==0));
                if(C.bas==C.bas_tr)
                    ver=false;
                    %length(find(sum(C.tr)==0))
                    C.adj(B~=0)=0;
                    C.adj(B'~=0)=0;
                    C.adj(A~=0)=-1;
                    C.adj(A'~=0)=1;
                end
            end
            
        end
        
        function C=Feas_r_bt(C)
            %Fasibility domain calculation via a ployhedra intersection
            %(not in use)
    
            addpath('/home/pmaldona/Documents/MATLAB/Add-ons/bt-1.1')
            S1=size(C.com);
            
            A=C.com;
            
            b=zeros(1,S1(1));
            c=zeros(1,S1(1));
            for i=1:S1(1)
                if(~any(C.tr(i,:)==1) || all(C.tr(i,:)==0))
                    b(i)=1;
                else
                    c(i)=1;
                end
            end
            
            B=diag(b);
            D=diag(c);
            
            
            A(S1(1)+1,:)=zeros(1,S1(1));
            P.V=normr(A');
            
            %P.D=normc(A');
            %P.V=zeros(S1(1),1);
            
            P=polyh(P,'v');
            
            Q.V=zeros(S1(1),1);
            Q.D=B;
            Q.L=D;
            Q=polyh(Q,'v');
            
            R=Q&P;
            clear Q;
            clear P;
            
            if(isempty(R) || dim(R)<(S1(2)+1))
                C.F = false;
                C.R=zeros(S1(1),1);
                C.X=C.R;
                return;
            end
            
            V=vrep(R);
            %r=mean(normc(V.D)');
            C.R=mean(normc(V.V)');
            C.R=C.R/norm(r,Inf);
            
            C.F=true;
            
        end
        
        function C=Feas_cen(C)
            %Feasibility and grow rate calulation via centroid calculation
            %(not in use)
            
            C.R=mean(normc(-C.com)');
            C.R=C.R'/norm(C.R,Inf);
            
            if(cond(C.com)<1e5)
                C.X=-C.com\C.R;
                if(all(C.X>0))
                    C.F=1;
                else
                    C.F=0;
                end   
            else
                C.X=-pinv(C.com)*C.R;
                C.F=0;
            end
            C.L=min(C.X);
        end
            

        
        function C = is_feasible(C)
            % LP-feasibility calulation.
            
            max_r=C.max_r;
            min_mort=C.min_mort;
            
            M=C.com;
            sign_of_M = sign(M) ;
            tmp = sign_of_M + sign_of_M' ;  % tmp[i,j]~=0 iff there is a non-trophic relation between i and j
            Trophic_M = sign_of_M ;
            Trophic_M(tmp~=0) = 0 ; % qualitative community matrix with only the signs of trophic interactions
            
            C.bas = length(find(sum(Trophic_M==-1)==0));  % number of basal species
            non_basals = find(sum(Trophic_M==-1)>0) ; % detects non-basal species
            
            n=length(M) ;
            Aeq=[eye(n) M zeros(n,1)] ;
            beq = zeros(n,1) ;
            
            % tmp = ones(n,1) ;
            % tmp(basals) = -1 ;
            % A = [ zeros(n) -eye(n) ones(n,1) ; diag(tmp) zeros(n) ones(n,1)] ;
            % b=zeros(2*n,1) ;
            
            A = [ zeros(n) -eye(n) ones(n,1)] ;
            b=zeros(n,1) ;
            
            ub = [max_r*ones(n,1) ; inf(n+1,1)] ;
            lb = -ub ;
            ub(non_basals) = -min_mort ;
            
            f = [zeros(2*n,1) ; -1] ;
            
            [x,~,exit_flag]=linprog(f,A,b,Aeq,beq,lb,ub,optimoptions('linprog','Display','off'));
            F = exit_flag >0 && x(end) > 0 ;
            
            if  F>0
                
                C.F=true;
                C.X=x(n+1:end-1);
                C.R=x(1:n);
                C.L=x(end);
            else
                C.F=false;
                C.X=zeros(n,1) ;
                C.R=C.X ;
                C.L=0 ;
            end
        end
        
        function C=esta(C)
            %assess local stability of matrix M
            le=length(C.com);
            e=eig(C.com);
            s=sum(sign(real(e)));
            if s==-le
                C.est=1;
            else
                C.est=0;
            end
            
        end
        
        function C=bc_assing(C)
            %Basal competitive interacion addind
            if(C.bc)
                bas=find(sum(C.tr)==0);
                
                for i=bas
                    for j=bas
                        if(j~=i)
                            C.adj(i,j)=-1;
                        end
                    end
                end
            end
        end
        
        
        function C= type_matrix(C)
            %type matrix calulation
            n=length(C.esp);
            C.tp=zeros(n);
            for i=1:n
                tv=find(C.adj(i,:));
                for ii=tv
                    if(C.adj(i,ii)==1)
                        if(C.adj(ii,i)==-1)
                            %tp(i,ii)='tr';
                            C.tp(i,ii)=1;
                        elseif(C.adj(ii,i)==1)
                            %tp(i,ii)='mu';
                            C.tp(i,ii)=2;
                        else
                            %tp(i,ii)='cm';
                            C.tp(i,ii)=3;
                        end
                    else
                        if(C.adj(ii,i)==1)
                            %tp(i,ii)='tr';
                            C.tp(i,ii)=-1;
                        elseif(C.adj(ii,i)==-1)
                            %tp(i,ii)='cp';
                            C.tp(i,ii)=4;
                        else
                            %tp(i,ii)='am';
                            C.tp(i,ii)=5;
                        end
                    end
                end
            end
            
            
        end
        
        function C=prop(C)
            %Summary of properties calulation
            
            if(C.mst)
                %                 C=C.Feas_r_bt();
                C.L=NaN;
                C=C.is_feasible();
                %             C=C.Feas_cen();
                C=C.esta();
                n=length(C.R);
            else
            
            C.X=-C.com\C.R;
            C.L=min(C.X);
            n=length(C.R);
            C.bas=n-length(find(C.Nv));
            
        end
            
            z=length(find(C.com))/(2*n);
            C.k=sum(C.com'~=0)/(2*z);
            C.k_p=sum(C.com'>0)/z;
            C.k_m=sum(C.com'<0)/z;
            
            M=diag(C.X)*C.com;
            C.m=mean(abs(M'));
            
            M(M>0)=0;
            C.m_m=mean(abs(M'));
            
            M=diag(C.X)*C.com;
            M(M<0)=0;
            C.m_p=mean(abs(M'));
            
            M=C.com;
            C.a=mean(abs(M'));
            
            M(M>0)=0;
            C.a_m=mean(abs(M'));
            
            M=C.com;
            M(M<0)=0;
            C.a_p=mean(abs(M'));
            %         if(cond(C.com)>1e5)
            %             C.com
            %         end
            
            C.res=-max(real(eig(diag(C.X)*C.com)));
            
            C.con=length(find(C.com))/n^2;
            s=n*C.con;
            
            comR=-C.com;
            dia=diag(comR);
            comR(logical(eye(size(comR))))=zeros(n,1);
            comR=comR./dia;
            
            C.me=s*mean2(comR);
            C.var=s*std2(comR)^2;
            C.gama=s/C.var*(mean2(transpose(comR).*comR)-C.me^2/s);
            C.eta=std(C.R./dia);
            %         [C.cp,C.mu,C.cm,C.am] = non_trophic_m(C);
            C=C.inter_porp();
            C.bas=length(find(sum(C.tr)==0));
        end
    end
    
end

