# Introduction


## 1. Cytokines and chemokines 

**1. General introduction to the cytokine superfamily:**

 * Biochemical and molecular property: 
   * What are they: Cytokines are small (typically 5~25 kd) regulatory proteins/peptides with pleiotropic functions.
   * Where are they synthesized and secreted
   * Which cells respond to them and how do they regulate downstream pathways: examples of typical receptors that are targeted by cytokines

 * Classification of the cytokine superfamily:
   * Cytokines include chemokines, interferons, interleukins, lymphokines, tumor necrosis factors etc.
   * Brief description of those different classes, with a few examples of key members for each class.
   * What is special about chemokines: molecules specially adapted for chemotaxis of cells. Crucial in directing the migration of immune cells to a target site (e.g. white blood cells to infected or damaged tissues)


**2. Their important roles in IMD**

* Broadly involved in a variety of biological processes (inflammation, haematopoiesis, development etc), but most importantly, immunity
* Involvement of cytokines in immunity: humoral and cell-mediated immunity
* Some examples of best studied members whose dysfunctions are reported to be causative in IMD:
  * (Following papers can be good reference. But we will mainly cover diseases/cytokines that are eventually used in our basis/projection.)
  * Inflammatory bowel disease (IBD)[Friedrich, 2019] (review)
  * Immune-mediated inflammatory diseases [Burkhard Möller 2006] (review)
  * Inflammatory disease [Shinwan Kany 2019] (review)

**3. Transition paragraph to part 2:**
* Why is it interesting to study the genetic architecture of cytokines; what is the method of choice.


## 2. GWAS

**1. Introduction to GWAS summary statistics**

* What is GWAS: What is its input; How does it work; What is its output.
* What is GWAS summary statistics
* Interpretation of key variables: p value, effect size, MAF, etc.

**2. Previous GWAS studies on cytokines**

**3. Challenges in applying GWAS to study cytokines in the context of IMD**

* Rare diseases do not provide sufficiently large sample size
* High dimension data (100,000s of SNPs); most SNPs are not causal driver SNPs (LD); hard to study shared genetic architectures across multiple variants and between multiple diseases simultaneously

**4. Transition paragraph to part 3:**
* Reduction of data dimesion with shrinkage + PCA as a solution to address the above problems


## 3. Bayesian shrinkage approach + PCA 

**1. Origin and definition** 
* Simplify the model by assuming that only one SNP has a causative role in each region.
* Show only the most general/important equations (Detailed maths for quantative traits should be in the method part)
* Show the graphic abstract of shrinkage, decomposition (basis) and projection (testing dataset).

**2. Strength of it: Why it helps us to address the problem stated in part 2.**
* Shrinkage weights are naturally adjusted for LD
* By keeping only one causal SNP for each region, avoid double counting, decrease computational complexity
* Mitigate overfitting: the effect from technical noise as well as the sampling variation are weakened in the final metric; capture biological signal, which is the true genetic architecture that we are interested in

**3. Briefly state the aim of this project:**
* Creating a low-dimension cytokine/chemokine basis that summarises the shared genetic risks across a set of IMD.
* Use that basis to project independent datasets, get insights in the underlying genetic factors contributing to traits (by interpreting principle components)








