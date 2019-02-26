function ensamblaje_full_free(M, max_r,min_mort)
% Ensambla una comunidad a partir de la matriz MASTER A
% A es matriz comunitaria, a_ij es el efecto de j sobre i
% Criterios de ensamblaje: (a) se parte de una basal. (b) la sp que llage
% debe tener al menos un recurso. (c) el ensamblaje prosigue en t+1 si la123
% matrix en t es factible y estable.
% los "r" se dejan cambiar en cada t.

sign_of_M = sign(M) ;
tmp = sign_of_M + sign_of_M' ;  % tmp[i,j]~=0 iff there is a non-trophic relation between i and j
Trophic_M = sign_of_M ;
Trophic_M(tmp~=0) = 0 ; % qualitative community matrix with only the signs of trophic interactions

basals = find(sum(Trophic_M==-1)==0);   % detects basal species
%non_basals = find(sum(Trophic_M==-1)>0) ; ; % detects non-basal species

[F,X,R,L] = is_feasible(M,max_r,min_mort); % Calculamos los "r" consistentes con equilibrios positivos, si los hay

Factibility = zeros (length(M),1) ;
Resilience = zeros(length(M),1) ;
%% assemblying the community E and sequence S
%
figure
for i=1:length(M)
    if i==1
    	init_basal=basals(randi([1 length(basals)])) ; % choses a basal at random
    	S=init_basal ;
        %E = M(init_basal,init_basal) ;
    else
    	Candidates_1 = setdiff(1:length(M), S) ; % selects elements of M that are not in E
        Candidates_2 = find(sum(M(S,:)<0)) ; % selects those species that have at least one resource in E
        Candidates_3 = intersect(Candidates_1,Candidates_2) ; % intersection
        Candidates   = union( Candidates_3, setdiff(basals, S) ) ; % adding basals
        assert(length(Candidates)>=1) ;
        Next_sp = Candidates( randi([1 length(Candidates) ])) ; % Choose randomly among candidates
        
    	S = [S Next_sp] ;
    end
    E = M(S,S) ;
    [F,X,RR,L] = is_feasible(E,max_r,min_mort);
    
    Eq = -E\RR ;
    
    Factibility(i) = min(Eq) ;
    Resilience(i)  = -max(real(eig(  diag(Eq)*E  ))) ;
    %aut=[real(eig(  diag(Eq)*E  )) real(eig(E))];
    J =  diag(Eq)*E  ;
    %
    subplot(6,6,i)
    plot(eig(J),'b*')
    grid
    %axis([-2 2 -.5 .5])    
end
% R
% Factibility
% Resilience
figure
plot(1:length(M),Factibility,'ro')
refline(0,0)
xlabel('Time')
ylabel('Factibility level')
figure
plot(1:length(M),Resilience,'bo')
refline(0,0)
xlabel('Time')
ylabel('Resilience')
figure
plot(eig(J),'ko')
grid
end




