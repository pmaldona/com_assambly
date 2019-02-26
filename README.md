# Ecological Community Assambly

A *Matlab* class that generates an ecological community assambly by adding species in funtion of their trophic relation. 
## Contents

- [Class](https://github.com/pmaldona/com_assembly#class)
-- [Constructor](https://github.com/pmaldona/com_assembly##constructor)


## Class

The class in defined in the file [NTCommunity.m](https://github.com/pmaldona/com_assembly/blob/master/src/NT_Comunity.m) which contain the Community constructor (master community). This comunity has a set of paramter that describe the intaraction model of the Community.  

### Constructor `NT_Community`

Call:
`NT_community(Tr,cp,fcp,mu,fmu,cm,fcm,am,fam,dnt,an,fan,rit,mit,fas,mnti,nbf,sg,sgd,tsp,msp,max_r,bc)`

#### Inputs
The inputs defines model and interaction type of the community matrix, and conditions of the grow rates:
- `Tr` Input trophic matrix, can be define by the Generlized Niche Model [Stouffer et. al](https://www.pnas.org/content/103/50/19015) `TheGenNicheModel.m`
- `cp` Fraction of added competitive interaction of the Community
- `fcp` Force of mutual interaction of the Commmedian
- `mu` Fraction of added mutual interaction of the Community
- `fmu` Force of competitive interaction of the Commmedian
- `cm` Fraction of added comensal interaction of the Community
- `fcm` Force of comensal interaction of the Community
- `am` Fraction added of anmensal interaction of the Community
- `fam` Force of anmensal interaction of the Community
- `an` Fraction of added trophic interaction of the Community (that aren't in `Tr`)
- `fan` Force of added trophic interaction of the Community
- `dnt` funtion for preferencial added interaction position (normalized tropic level domain)
- `rit` Range (standard deviation) of the distance between interact species of added interactions $$\sigma_r$$ , in units of normalized trophic level
- `mit` Median of the distance between interact species of added interactions $$\mu_r$$ , in units of normalized trophic level
- `fas` Asmietric force
- `mnti` median of coefficents for added interactions
- `nbf` non basal diagonal factor
- `sg` standard deviaton of non diagonal values
- `sgd` standard deviation of diagonal values
- `tsp` non-same interaction spending time factor (0 only same type spending , 1 all type spendig)
- `msp` overall spendig time factor (1 non-spendig time model, 0 Holling type-1 model)
- `max_r` maximum value of grow-rate (imposition for LP-optimization)
- `min_mort` minimum value for de minimal of biomasses (imposition for LP-optimization)
- `bc` basal competition (bool value, that if it's true, add competition between basal species)
 #### Model
Our goal is analize the stablitiy and feasibility of matrix a bases dynamics ecology system model. Where the vector biomases <img src="https://latex.codecogs.com/svg.latex?\Large&space;\bold{x}"> a  determined by:

$$\dot{\bold{x}}=\bold{x}\left(\bold{r} + \bold{A} \bold{x}\right)$$

Here $$\bold{A}$$ is known as Community Matrix, and  $$\bold{r}$$ as grow rate. 

In first instace, initial tropic interactions are obtainend form the `Tr` matrix, this will be defined as trophic skeleton. The number $$n$$ of species correspond to the `Tr`matrix dimension. A requirment for `Tr` is to be fully conected. Interactions are added by type proportional to values `cp`, `mu`, `cm`, `am` ,`an` according to the total connectance $$n(n-1)/2$$. The position of the interactions are randomly choosen proportional to an funtion depends on the trophic levels if the interacting species. This function correspond to:

$$\mathbb{P}_{ij} \propto f_d(tl_n(i))\exp\left({\frac{|tl_n(i)-tl_n(j)|-\mu_{it}}{\sigma_r}}\right)$$

where $$tl_n(i)$$ is the normalized trophic level ([Livine](https://www.sciencedirect.com/science/article/pii/002251938090288X)) of species $$i$$, $$f_d$$ correspond to the `dfn` function (matlab `@` call), and $$\mu_{ti}$$ and $$\sigma_r$$ correspondo to `mit`and `rit` respectivly. In the adding interaction process, the basal species of `Tr` are preserved. Finally the bool value `bc` add competition between basal species.

The non-diagonal community matrix coeficients $$A_{ij}$$ are defined by holling type I model:

$$A_{ij}=f_{t} \frac{b_{ij}C_{ij}}{msp+(1-msp)\left(\sum_{k\in \mathcal{T}_i} C_{ik}+tp\sum_{k\in0 \mathcal{A}_i\setminus\mathcal{T}_i}C_{ik}\right)}$$


$$\mathcal{A}_i$$ is the set of all connected species with $$i$$ (that interact with $$i$$), and $$\mathcal{T}_i$$ is the set of species that have same type intraction of $$i$$ with $$j$$ considering this, $$tp=$$ `tp` and $$msp=$$ `msp`. Coefficients $$C_{ij}$$ and $$b_{ij}$$ are random generated from an beta distribution multiplied by an median factor: $$\mu \mathcal{B}(\alpha,\alpha)$$ with $$\alpha=\mu/(2\sigma^2)$$, here $$\mu=$$ `mnti` and $$\sigma=$$ `sg`.  In case that $$ij$$ is a predation interacion, $$b_{ij}=1$$. $$f_t$$ are force of interacion, this values are taken form `fcp`,`fmu`, `fcm`, `fam` and `fan`. Interactions that belongs to the `Tr` matriz don't have force factor. 
The diagonal elements of community matrix are randomly genereted by normal distribution with mean 1 and standard deviation $$\sigma_d$$ = `sgd`: $$A_{ii} \sim \mathcal{N}(1,\sigma_d)$$. We add an bool variable `nbf` that makes diagonal values for basal species zero.

Our aim is to have local stability of the community. For this it's necessary that the biomasses change rate be cero, *i.e.*:

$$\dot{\bold{x}^*}=\bold{x}^*\left(\bold{r} + \bold{A} \bold{x}^*\right)=\bold{0}$$

so we have:

$$\bold{r}=-\bold{A}\bold{x^*}$$

$$\bold x^*$$ represents satbility point for the biomases. At this point we have freedom of choise for the grow rate $$\bold r$$ and so $$\bold x^*$$ due there are releted by lineal transformation. To solve this, we search the feasibility via an LP-optimization that maximize the minimal of bioases $$\text{min}(\{x_i\})$$ holding the related equality, boundig the grow rates via a Chebyshov norm, $$\|\bold{r}\|_{\infty}< m_r$$ and lower bound non-basal grow rates $$r_i>m_m$$.   Here $$m_m=$$ `min_mort` and $$m_r=$$ `max_r`. 

#### Outputs

The output constitute a structure where can be obtained a set of properties and variables that builds the community, The strucutre is sparatre by type data:

##### Matrices
- `adj` Adjoint matrix of the community
- `com` Community matrix ($$\bold A$$)
- `tr` Trophic matrix
- `tr_i` Trophic skeleton
- `tp` Type interaction matrix
- `coef_C` $$C_{ij}$$ matrix coeficients
- `coef_b` $$b_{ij}$$ matrix coeficients

##### Vectors
- `Nv` Trophic levels vector ($$\bold{tl}$$) [Livine](https://www.sciencedirect.com/science/article/pii/002251938090288X)
- `esp` Species index vector
- `X` Biomasses vector at stability ($$\bold x$$)
- `R` Grow rate vector ($$\bold r$$)
- `k` Vector of porcentage of interaction [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `k_p` Vector of porcentage of positive interaction [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `k_m` Vector of porcentage of negative interaction [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `m` Vector of mean population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `m_p` Vector of mean positive population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `m_m` Vector of mean negative population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `a` Vector of mean percapita interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `a_p` Vector of mean percapita population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)
- `a_m` Vector of mean percapita population interaction strength [Wootton & Stouffer](https://link.springer.com/article/10.1007/s12080-015-0279-3)

##### Doubles
- `cp` Porcentage of cometitive interactions
- `mu` Porcentage of mutual interactions
- `am` Porcentage of amensal interactions
- `cm` Porcentage of comensal interactions
- `an` Porcentage of trophic interactions
- `res` Community resilence ($$R$$)
- `res_tr` Trophic seleton community resilence ($$R$$)
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

Resilence is equal to minus the maximun eiqgenvalue of the Jacobiean evaluated at the stablity point $$\bold x^*$$:

$$R=-\lambda^{\uparrow}_1(J(x^*))=-\lambda^{\uparrow}_1(\text{diag}(\bold x^*) \bold A)$$

So, if $$R > 0$$ then the community is localy satble, if $$R < 0$$ then is localy instable.

##### Bools
- `mst` True if community is a aster one (first constructed)
- `est` True if comunity is locally stable
- `F` True if community is feasible
- `bc` True if there are competition between basal species
- 
### Init assambly `NT_SubCommunity_init_Links`
Call: `Sc=NT_SubComunity_init_Links(C,S)`
#### Inputs
- `S` Number of desire species on sub-community
- `C` Master community Class object
#### Model
Constructions of sub-community starts from the master community (first initalized community). The process of assambly begins by taking an basal species. Once achived this, candidates to add are direct predetors of present species or new basal ones. This algorithm iterates until the `S` species are reached. The grow rates of the sub-community are directly taked from the master community `C` grow rates, as a subset of the present species. Sub-Community matrix is constructed in the same way as the master community, but only considering the susbet species in the sub-community. Then biomasses is calculeted via the known lineal trasnformation ($$\bold{r}=-\bold{A}\bold{x^*}$$) and finally all other propertires (resilence, min of biomasses, etc.).
### Outputs 
Returns a class object as constructor `NT_community`

### Add assambly `NT_SubComunity_add_Links`
Call: `Sc=NT_SubComunity_add_Links(Si,C,S)`
#### Inputs
- `S` Number of desire species to add
- `C` Master community Class object
- `Si` Sub-community Class object
 #### Model
Constructor of sub-community starting of the sub-community `Si` adding `S` desierd species from master community `C`.
Has the same porcedure as ``NT_SubCommunity_init_Links``
#### Outputs
Returns a class object as constructor `NT_community`
