function com_scale(Tr,cpi,fcp,mui,fmu,cmi,fcm,ami,fam,df,ani,fan,rit,mit,fas,mnti,nb,sg,sgd,tsp,msp,max_r,min_mort,rep,ind,bc)

n=length(Tr);
b=find(sum(Tr)==0);
C=(length(find(Tr))*2+n+length(b)^2)/n^2;
pm=cpi+mui+cmi+ami+C;
dn=func2str(df);

if(pm>1)
    warning('La suma porcental las interacciones no mutualistas no puede ser mayor que el 100%');
    return;
end





name=sprintf('com_scale_NM_ind(%d)',ind);
fin=strcat('begin_',name);
disp(fin);
if(size(dir(strcat(name,'.mat')),1)>0)
    warning('existe archivo')
    return;
else
    [A,B,R,R_tr,L,L_tr,bas_e,con,con_tr,cp,mu,am,cm,an,me,me_tr,gama,gama_tr,var,var_tr,eta,eta_tr,ict,Nv,dist,x,r] = Com_construc(Tr,cpi,fcp,mui,fmu,cmi,fcm,ami,fam,df,ani,fan,rit,mit,fas,mnti,nb,sg,sgd,tsp,msp,max_r,min_mort,rep,bc);
    save(strcat(name,'.mat'),'R','R_tr','L','L_tr','bas_e','con','con_tr','cp','mu','am','cm','an','me','me_tr','gama','gama_tr','var','var_tr','eta','eta_tr','ict','Nv','dist','x','r','cpi','fcp','mui','fmu','cmi','fcm','ami','fam','ani','fan','rit','mit','fas','mnti','nb','sg','sgd','tsp','msp','max_r','min_mort','rep','ind','dn');
end

fin=strcat('end_',name);
disp(fin);

end
