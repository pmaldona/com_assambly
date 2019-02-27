# Ecological Community Assembly

A *Matlab* class that generates an ecological community assembly by adding species in function of their trophic relation. 
## Contents

- [Class](https://github.com/pmaldona/com_assembly#class)
  - [Constructor](https://github.com/pmaldona/com_assambly#constructor-nt_community)
    - [Inputs](https://github.com/pmaldona/com_assambly#inputs)
    - [Model](https://github.com/pmaldona/com_assambly#model)
    - [Outputs](https://github.com/pmaldona/com_assambly#outputs)
      - [Matices](https://github.com/pmaldona/com_assambly#matrices)
      - [Vectors](https://github.com/pmaldona/com_assambly#vectors)
      - [Doubles](https://github.com/pmaldona/com_assambly#doubles)
      - [Bools](https://github.com/pmaldona/com_assambly#bools)
  - [Init assembly](https://github.com/pmaldona/com_assambly#init-assambly-nt_subcommunity_init_links)
    - [Inputs](https://github.com/pmaldona/com_assambly#inputs-1)
    - [Model](https://github.com/pmaldona/com_assambly#model-1)
    - [Outputs](https://github.com/pmaldona/com_assambly#outputs-1)
  - [Add to assembly](https://github.com/pmaldona/com_assambly/blob/master/README.md#add-to-assambly-nt_subcomunity_add_links)
    - [Inputs](https://github.com/pmaldona/com_assambly#inputs-2)
    - [Model](https://github.com/pmaldona/com_assambly#model-2)
    - [Outputs](https://github.com/pmaldona/com_assambly#outputs-2)



## Class

The class in defined in the file `NTCommunity.m` which contain the Community constructor (master community). This community has a set of parameter that describe the interaction model of the Community.  

### Constructor `NT_Community`

Call:
`NT_community(Tr,cp,fcp,mu,fmu,cm,fcm,am,fam,dnt,an,fan,rit,mit,fas,mnti,nbf,sg,sgd,tsp,msp,max_r,tr_cal,res,dr,bc)`

#### Inputs
The inputs defines model and interaction type of the community matrix, and conditions of the grow rates:
- `Tr` Input trophic matrix, can be define by the Generalized Niche Model [Stouffer et. al](https://www.pnas.org/content/103/50/19015) `Tr_init.m`
- `cp` Fraction of added competitive interactions of the Community
- `fcp` Force of mutual interactions of the Community
- `mu` Fraction of added mutual interactions of the Community
- `fmu` Force of competitive interaction of the Community
- `cm` Fraction of added comensal interactions of the Community
- `fcm` Force of comensal interactions of the Community
- `am` Fraction added of anmensal interactions of the Community
- `fam` Force of anmensal interactions of the Community
- `an` Fraction of added trophic interaction of the Community (that aren't in `Tr`)
- `fan` Force of added trophic interaction of the Community
- `dnt` Function for preferential added interaction position (normalized tropic level domain)
- `rit` Range (standard deviation) of the distance between interact species of added interactions <img src="https://latex.codecogs.com/svg.latex?\sigma_r"> , in units of normalized trophic level
- `mit` Mean of the distance between interact species of added interactions <img src="https://latex.codecogs.com/svg.latex?;\mu_r"> , in units of normalized trophic level
- `fas` Asymmetric force
- `mnti`Coefficients mean of added interactions
- `nbf` Non basal diagonal factor
- `sg`  Non diagonal values standard deviation
- `sgd` Diagonal values standard deviation
- `tsp` Non-same interaction spending time factor (0 only same type spending , 1 all type spendig)
- `msp` Overall spendig time factor (1 non-spending time model, 0 Holling type-1 model)
- `max_r` Maximum value of grow-rate (imposition for LP-optimization)
- `min_mort` Minimum value of non basal grow-rate (imposition for LP-optimization)
- `tr_cal` Bool value to force resilience of the trophic initial community (tophic skeleton community)
- `res` Resilience for the trophic initial community (if `tr_cal` is `true`)
- `dr` Range of acceptance for resilience(`res` <img src="https://latex.codecogs.com/svg.latex?\pm"> `dr`, if `tr_cal` is `true`)   
- `bc` basal Competition (bool value, that if it's true, add competition between basal species)
 #### Model
Our goal is analize the stablitiy and feasibility of matrix a bases dynamics ecology system model. Where the vector biomases <img src="https://latex.codecogs.com/svg.latex?\bold{x}"> a  determined by:

![img](https://latex.codecogs.com/svg.latex?\dot{\bold{x}}=\bold{x}\left(\bold{r}+\bold{A}\bold{x}\right))

Here <img src="https://latex.codecogs.com/svg.latex?\bold{A}"> is known as Community Matrix, and  <img src="https://latex.codecogs.com/svg.latex?\bold{r}"> as grow rate. 

In first instance, initial tropic interactions are obtained form the `Tr` matrix, this will be defined as trophic skeleton. The number <img src="https://latex.codecogs.com/svg.latex?n"> of species correspond to the `Tr` matrix dimension. A requirement for `Tr` is to be fully connected. Interactions are added by type, proportional to values `cp`, `mu`, `cm`, `am` ,`an` according to the total connectance <img src="https://latex.codecogs.com/svg.latex?n(n-1)/2">. The position of interactions are randomly choose proportional to an function depends on the trophic levels if the interacting species. This function correspond to:

![img](https://latex.codecogs.com/svg.latex?\mathbb{P}_{ij}\propto%20f_d(tl_n(i))\exp\left({\frac{|tl_n(i)-tl_n(j)|-\mu_{it}}{\sigma_r}}\right))

where <img src="https://latex.codecogs.com/svg.latex?;tl_n(i)"> is the normalized trophic level ([Livine](https://www.sciencedirect.com/science/article/pii/002251938090288X)) of species <img src="https://latex.codecogs.com/svg.latex?i">, <img src="https://latex.codecogs.com/svg.latex?f_d"> correspond to the `dfn` function (matlab `@` call), and <img src="https://latex.codecogs.com/svg.latex?\mu_{ti}"> and <img src="https://latex.codecogs.com/svg.latex?\sigma_r"> correspond to `mit`and `rit` respectively. In the adding interactions process, the basal species of `Tr` are preserved. Finally the bool value `bc` add competition between basal species.

The non-diagonal community matrix coefficients <img src="https://latex.codecogs.com/svg.latex?A_{ij}"> are defined by holling type I model:

![img](https://latex.codecogs.com/svg.latex?A_{ij}=f_{t}\frac{b_{ij}C_{ij}}{msp+(1-msp)\left(\sum_{k\in\mathcal{T}_i}C_{ik}+tp\sum_{k\in%20\mathcal{A}_i\setminus\mathcal{T}_i}C_{ik}\right)})


<img src="https://latex.codecogs.com/svg.latex?\mathcal{A}_i"> is the set of all connected species with <img src="https://latex.codecogs.com/svg.latex?i"> (that interact with <img src="https://latex.codecogs.com/svg.latex?i">), and <img src="https://latex.codecogs.com/svg.latex?\mathcal{T}_i"> is the set of species that have same type interaction of <img src="https://latex.codecogs.com/svg.latex?i"> with <img src="https://latex.codecogs.com/svg.latex?\j"> considering this, <img src="https://latex.codecogs.com/svg.latex?tp="> `tp` and <img src="https://latex.codecogs.com/svg.latex?msp="> `msp`. Coefficients <img src="https://latex.codecogs.com/svg.latex?C_{ij}"> and <img src="https://latex.codecogs.com/svg.latex?b_{ij}"> are random generated from an beta distribution multiplied by an median factor: <img src="https://latex.codecogs.com/svg.latex?\mu\mathcal{B}(\alpha,\alpha)"> with <img src="https://latex.codecogs.com/svg.latex?\alpha=\mu/(2\sigma^2)">, here <img src="https://latex.codecogs.com/svg.latex?\mu="> `mnti` and <img src="https://latex.codecogs.com/svg.latex?\sigma="> `sg`.  In case that <img src="https://latex.codecogs.com/svg.latex?ij"> is a predation interacion, <img src="https://latex.codecogs.com/svg.latex?b_{ij}=1">. <img src="https://latex.codecogs.com/svg.latex?f_t"> are force of interacion, this values are taken form `fcp`,`fmu`, `fcm`, `fam` and `fan`. Interactions that belongs to the `Tr` matrix don't have force factor. 
The diagonal elements of community matrix are randomly genereted by normal distribution with mean one and standard deviation <img src="https://latex.codecogs.com/svg.latex?\sigma_d"> = `sgd`: <img src="https://latex.codecogs.com/svg.latex?A_{ii}\sim\mathcal{N}(1,\sigma_d)">. We add an bool variable `nbf` that makes diagonal values for basal species zero.

Our aim is to have local stability of the community. For this it's necessary that the biomasses change rate be cero, *i.e.*:

<img src="https://latex.codecogs.com/svg.latex?\dot{\bold{x}^*}=\bold{x}^*\left(\bold{r}+\bold{A}\bold{x}^*\right)=\bold{0}">

so we have:

<img src="https://latex.codecogs.com/svg.latex?\bold{r}=-\bold{A}\bold{x^*}">

<img src="https://latex.codecogs.com/svg.latex?\bold{x}^*"> represents stability point of biomasses. At this point we have freedom of choise for the grow rates <img src="https://latex.codecogs.com/svg.latex?\bold{r}"> and so <img src="https://latex.codecogs.com/svg.latex?\bold{x}^*"> (biomasses), due there are related by lineal transformation. To solve this, we search the feasibility via an LP-optimization that maximize the minimal of biomasses <img src="https://latex.codecogs.com/svg.latex?\text{min}(\{x_i\})">, holding the related equality, boundig the grow rates via a Chebyshov norm, <img src="https://latex.codecogs.com/svg.latex?\|\bold{r}\|_{\infty}<m_r"> and lower bound non-basal grow rates <img src="https://latex.codecogs.com/svg.latex?r_i>m_m">.   Here <img src="https://latex.codecogs.com/svg.latex?m_m="> `min_mort` and <img src="https://latex.codecogs.com/svg.latex?m_r="> `max_r`.

Resilience is equal to minus the maximum eigenvalue of the Jacobean evaluated at the stability point <img src="https://latex.codecogs.com/svg.latex?\bold{x}^*">:

<img src="https://latex.codecogs.com/svg.latex?R=-\lambda^{\uparrow}_1(J(x^*))=-\lambda^{\uparrow}_1(\text{diag}(\bold{x}^*)\bold{A})">

So, if <img src="https://latex.codecogs.com/svg.latex?\Large&space;R>0"> then the community is locally stable, if <img src="https://latex.codecogs.com/svg.latex?\Large&space;R<0"> then is locally unstable.

We first calculate the resilience and feasibility of the trophic skeleton community matrix, associated to their interactions (without considering added interactions). If `tr_cal` = `true` the program it will calculate resilience until <img src="https://latex.codecogs.com/svg.latex?R="> `res`<img src="https://latex.codecogs.com/svg.latex?\pm"> `dr` is reached (this step it's a little tricky because it enters in a while loop, and the program maybe can't exit the loop. It's better use `tr_cal`=`false`). After this, the rest interactions are adeed and feasibility and resilience calculated.



#### Outputs

The output constitute a structure that contain a set of properties and variables that builds the community, The strucutre is sparatre by type data:

##### Matrices
- `adj` Adjoint matrix of the community
- `com` Community matrix (<img src="https://latex.codecogs.com/svg.latex?\bold{A}">)
- `tr` Trophic matrix
- `tr_i` Trophic skeleton
- `tp` Type interaction matrix
- `coef_C` <img src="https://latex.codecogs.com/svg.latex?\Large&space;C_{ij}"> matrix coeficients
- `coef_b` <img src="https://latex.codecogs.com/svg.latex?\Large&space;b_{ij}"> matrix coeficients

##### Vectors
- `Nv` Trophic levels vector (<img src="https://latex.codecogs.com/svg.latex?\bold{tl}">) [Livine](https://www.sciencedirect.com/science/article/pii/002251938090288X)
- `esp` Species index vector
- `X` Biomasses vector at stability (<img src="https://latex.codecogs.com/svg.latex?\bold{x}">)
- `R` Grow rate vector (<img src="https://latex.codecogs.com/svg.latex?\bold{r}">)
- `k` Vector of percentage of interaction [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `k_p` Vector of percentage of positive interaction [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `k_m` Vector of percentage of negative interaction [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `m` Vector of mean population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `m_p` Vector of mean positive population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `m_m` Vector of mean negative population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `a` Vector of mean percapita interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `a_p` Vector of mean percapita population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `a_m` Vector of mean percapita population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)

##### Doubles
- `cp` Porcentage of competitive interactions
- `mu` Porcentage of mutual interactions
- `am` Porcentage of amensal interactions
- `cm` Porcentage of comensal interactions
- `an` Porcentage of trophic interactions
- `res` Community resilience (<img src="https://latex.codecogs.com/svg.latex?R">)
- `res_tr` Trophic skeleton community resilence (<img src="https://latex.codecogs.com/svg.latex?R">)
- `con` Community connectance
- `con_tr` Trophic skeleton community connectance
- `me` Average total interaction faced for one species [Barbier & Arnoldi](http://dx.doi.org/10.1101/147728)
- `me_tr`Average total interaction faced for one species of trophic skeleton community [Barbier & Arnoldi](http://dx.doi.org/10.1101/147728)
- `gama` Carring capacity deviation [Barbier & Arnoldi](http://dx.doi.org/10.1101/147728)
- `var` Total interaction deviation faced for one species [Barbier & Arnoldi](http://dx.doi.org/10.1101/147728)
- `var_tr` Total interaction deviation faced for one species of trophic skeleton community [Barbier & Arnoldi](http://dx.doi.org/10.1101/147728)
- `eta` reciprocal interaction effect [Barbier & Arnoldi](http://dx.doi.org/10.1101/147728)
- `eta_tr` reciprocal interaction effect of trophic skeleton community [Barbier & Arnoldi](http://dx.doi.org/10.1101/147728)
- `bas` number of basal species of the community
- `L` minimal of biomasses
- `L_tr` minimal of biomasses of the trophic skeleton community

##### Bools
- `mst` True if community is a aster one (first constructed)
- `est` True if community is locally stable
- `F` True if community is feasible
- `bc` True if there are competition between basal species
- 
### Init assambly `NT_SubCommunity_init_Links`
Call: `Sc=NT_SubComunity_init_Links(C,S)`
#### Inputs
- `S` Number of desire species on sub-community
- `C` Master community Class object
#### Model
Constructions of sub-community starts from the master community (first initialized community). The process of assembly begins by taking an basal species. Once achieved this, candidates to add are direct predators of present species or new basal ones. This algorithm iterates until the `S` species are reached. The grow rates of the sub-community are directly taken from the master community `C` grow rates, as a subset of the present species. Sub-Community matrix is constructed in the same way as the master community, but only considering the susbet species in the sub-community. Then biomasses is calculated via the known lineal transformation (<img src="https://latex.codecogs.com/svg.latex?\Large&space;\bold{r}=-\bold{A}\bold{x}^*">) and finally all other properties (resilience, min of biomasses, etc.).
#### Outputs 
Returns a class object as constructor `NT_community`

### Add to assambly `NT_SubComunity_add_Links`
Call: `Sc=NT_SubComunity_add_Links(Si,C,S)`
#### Inputs
- `S` Number of desire species to add
- `C` Master community Class object
- `Si` Sub-community Class object
 #### Model
Constructor of sub-community starting of the sub-community `Si` adding `S` desired species from master community `C`.
Has the same procedure as ``NT_SubCommunity_init_Links``
#### Outputs
Returns a class object as constructor `NT_community`

