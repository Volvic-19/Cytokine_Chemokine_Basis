# Introduction


## 1. Cytokine and chemokine 

**1. general introduction to cytokine superfamily:**

 * Biochemical and molecular property: 
   * What are they: Cytokines are small (typically 5~25 kd) regulatory proteins/peptides with pleiotropic functions.
   * Where are they syntehsized and secreted
   * Who responds to them, how do they regulate downstream pathways: Examples of typical receptors that are targeted by cytokines

 * Classification of cytokine super family:
   * Cytokine includes chemokines, interferon, interleukins, lymphokines, tumor necrosis factor etc.
   * Brief description of those different classes, with a few examples of key members in each class.
   * What is special about chemokine: specially adapted for chemotaxis of cells. Crucial in directing the migration of immune cells to a target sie (e.g. Direct the migration of white blood cells to infected or damaged tissues)


**2. their important roles in IMD**

* Broadly involved in a variety of biological processes (inflmmation, hematopoiesis, development etc), but most importnantly, immunity
* Involvement of cytokine in immunity: humoral and cell-mediated immunity
* Some examples of best studied members whose dysfunctions are reported to be causal in IMD:
  * (Following papers can be good reference. But prefer covering diseases/cytokines that are used in our basis/projection.)
  * Inflammatory bowel disease (IBD)[Friedrich, 2019] (review)
  * Immune-mediated inflammatory diseases [Burkhard Möller 2006] (review)
  * Inflammatory disease [Shinwan Kany 2019](review)

**3. Transition paragraph to part 2:**
* why it is tempting to study genetic architecture of cytokines; what is the method to choose.


## 2. GWAS

**1. introduction to GWAS summary statistics**

* What is GWAS: What does it take as input; How does it work; What is its output.
* What is GWAS summary statistics
* Interpretation of key variables: p value, effect size, MAF, etc.

**2. Previous GWAS studies on cytokines**

**3. Challenges to apply GWAS to study cytokines in the context of IMD**

* rare disease does not have enough sample size
* high dimension data (100,000s of SNPs); most SNPs are not causal driver SNPs (LD); hard to study shared genetic architectures across multiple variants and between multiple diseases simultaneously

**4. Transition paragraph to part 3:**
* reduction of dimesion with shrinkage + PCA as a solution to address above problems


## 3. Bayesian shrinkage approach + PCA 

**1. origin, definition** 
* We simplify model by the assumption that maximum one SNP is causal in each region.
* Show only the most general/important equations (Detailed math for quantative traits should be in method part)
* Show graphic abstract of shrinkage, decomposition (basis) and projection (testing dataset).

**2. strength of it: Why it helps us address the above problem in part 2.**
* shriknage weight are naturally adjusted for LD
* by keeping maximun one calsual SNP in one region, aovid double counting, decrease computational complexity
* mitigate overfitting: effect from technical noise and sampling variation are weakened in final metric, aim to capture biological signal, which is the true genetic articheture that we are interested in

**3. briefly restate the aim of this project:**
* Creating a low-dimension cytokine/chemokine basis that summarises the shared genetic risks across a set of IMD.
* Use that basis to project independent datasets, get insights in the underlying genetic factors contributing to traits (by interpreting principle components)








