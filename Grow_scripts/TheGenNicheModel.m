function   A=TheGenNicheModel(C,S,CON)
%LGN
%ACTUALIZADO 13/AGOSTO/2013
%Modificado 28/NOVIEMBRE/2013 == Ahora tiene como output el valor de nicho
%aleatorio y ordenado
%MODIFICADO DICIEMBRE/2013 PARA TENER COMO OUTPUTS VALORES DE NICHO
%Adapted from Williams & Mart?nez 2000; stouffer et al 2006
%Modification Jan 03, 2011 (change the algorith to work with theb General Niche Model (GNM))

%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%
%C=Conectancia
%S=Riqueza de especies
%CON=Valor de contigudad de la dieta
%%%%%%%%%%%%%%%%%%%%%%%%%%OUTPUTS%%%%%%%%%%%%%%%%%%%%%%
%A=Matriz de adyascencia
%Nsort=Valor de nicho ordenado a partir del que se genera en definitiva la
%topolog?a
%Nrand=Valor de nicho aleatorio



%NICHE VALUE
%assign a "niche value" parameter to each i species drawn uniformly from the
%interval [0,1] that has mean=1/2*(0+1)=0.5
N=unifrnd(0,1,S,1);
Nrand=N;
%Assume that species are ordered by niche values n1<ni
Nsort=sortrows(N);


%RANGE VALUE
%assign the range value to each i species (beta no puede ser <0)
Beta=(1/(2*C))-1;
X=betarnd(1,Beta,S,1);
RNM=Nsort.*X;%cambie R por RNM (rango niche model), ya que ahora el rango es un rango reducido para generalized niche model

%GNM
%GENERALIZED NICHE MODEL
%for the generalized niche model a new parameter is incorporated, CON = {0,1}, and the
%range, R is reduced
R=CON*RNM;


%AT LEAST ONE BASAL SPECIES
%by settig the range of species with the lowest ni =zero
R(1,1)=0;
RR=R;

%CENTRE OF THE RANGE (adapted from )
%species i consumes all species falling in a range ri, that is placed by
%uniformly drawing the centre of the range Ci
Ce=zeros(S,1);
RRhalf=RR/2;%dividing the range in two parts

    for i=1:S  
        Ce(i)=unifrnd(RRhalf(i),Nsort(i));
        
    end
%GNM
%Sin embargo las especies siguen estando distribuidas uniformemente en el
%eje de los recursos, entonces qu? pasa con las especies que quedan sin
%depredador luego de la reduccion del rango? las presas esperadas sin
%ubicacion luego de la reduccion del rango son

deltak=round((1-CON)*RNM*S);%esta cantidad de presas son escogidas al azar de las i especies con valor de nicho i menor o igual que el valor de nicho de la especie j

    
    
    
%ADJACENCY MATRIX
Ceminus=Ce-RRhalf;
Ceplus=Ce+RRhalf;
A=zeros(S,S);
    for i=1:S
    A(:,i)=Ceplus(i)>Nsort & Nsort>Ceminus(i);
    end
    
 pj=find(deltak>0);%encuentra lo depredadores que carecen de algunas presas por la disminucion del rango
 nj=deltak(pj);%indica cuantas presas por depredador hay que agregar
 subA=A(:,pj);%rescata los vectores columna de cada uno de ellos de la matriz A (en estos estan contenidas sus presas)
    
 for k=1:size(pj,1)%este for es para ir depredador a depredador
    
     sp=nj(k);
     
     predpreys=subA(:,k);%este es el vector columna que contiene las presas del depredador k 
     
    %for l=1:sp%este for es para ir presa a presa dependiendo de cuantas presas haya que agregar
    
    while sum(predpreys)<sum(subA(:,k))+sp 
        
         potprey=find(subA(:,k)==0);%identidad de potenciales presas
    
         potprey2=Nsort(potprey)<=Nsort(pj(k,1));
        
         if sum(potprey2) == 0 || sum(potprey2)<sp == 1
            break    
         end   
         
         efprey=round(unifrnd(1,sum(potprey2),1,1));%escoge aleatoriamente cual o cuales de las posibles presas sera presa efectiva
        
         predpreys(potprey(efprey),1)=1;%da la identidad la nueva presa    
    end

    if sum(potprey2) == 0 || sum(potprey2)<sp == 1
         break    
    end

    A(:,pj(k))=predpreys;
 
end

end
 
