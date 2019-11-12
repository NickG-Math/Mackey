\mainpage General Information
\tableofcontents
 \section intro Introduction
This is a C++ header only library devoted to numerically computing the \f$RO(G)\f$ homology of a point. It is hosted <a href="https://github.com/NickG-Math/Mackey">here</a>

For a quick demonstration in the case of \f$G=C_4\f$ you can use one of the available binaries <a href="https://github.com/NickG-Math/Mackey/tree/master/bin">here</a>.


\section req Requirements
 * C++17 and the standard library.
 * <a href=" http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</a>, a header only library for matrix manipulation. I've tested this with Eigen 3.3.7
 * Optional: For improved performance you can use the Intel MKL with Eigen and further combine with OpenMP for multithreading.
 * Optional: To draw the multiplication graphs you will need Graphviz.
 
\section install Installation

* To install simply clone/download the folder <a href="https://github.com/NickG-Math/Mackey">repository</a> and include it in your path. You will also need to do the same with Eigen.

* See the page \ref use for a tutorial on using the library.

* As for compiler support, I have tested the code with the following C++ compilers: GCC 9.2 (Linux), Clang 10 (Linux and MacOS), Intel Compiler 19 (Linux and Windows), MSVC 19 (Windows). Remember to use the option <CODE>-std=c++17</CODE>. For more information on compiler options, see the \ref perf page.


\section status Current Status

* The project is effectively complete for \f$G\f$ a cyclic group of prime power order. The only input that's needed are the equivariant chains at the bottom level for the spheres corresponding to nonnegative linear combinations of irreducible representations; we call these "standard chains". The standard chains can be easily computed from geometric equivariant decompositions by hand, and then fed into the program as explained in \ref how. It might be worth it to automate this process as well; after all, the differentials of the standard chains can all be obtained using the fact that the homology at the bottom level has to be trivial apart from top dimension.

* For general cyclic groups a few aspects that involve transferring need reworking. The problem is that non prime-power cyclic groups the diagram of subgroups is not a vertical tower but a somewhat more complicated diagram, so care has to be taken to account for all these extra transfers and restrictions. Ultimately this is the only part that needs changing.

* The bulletpoint above also applies to general finite abelian groups. We also need to specify the order of the elements of the group and how they relate to the subgroup diagram to form our equivariant bases.

* For non abelian groups we have the added complication of needing the real representation theory of our group. And of course we need the standard chains for these groups as well. 

* For coefficients other than \f$\mathbb Z\f$ a lot more things start to break, as transferring becomes more complicated as non cyclic modules are involved in the free Mackey functors.

\section doc Documentation

This documentation is organized in pages:
 
* The  <a href="namespaces.html">Related pages</a> consisting of \ref index, \ref math, \ref use, \ref algo, \ref perf. These explain how the program works, starting from the math and moving to more technical territory regarding the actual implementation. I recommend starting from these.

* The pages <a href="namespaces.html">Namespaces</a>, <a href="annotated.html">Classes</a>  and <a href="files.html">Files</a> that are automatically generated pages by doxygen from the source code. These offer a much more indepth look into all classes and functions of this project. Note that only public and protected members and named namespaces are documented.

* I recommend starting with the <a href="namespaces.html">related pages</a> before moving and then moving to the other tabs. If you just want to use this library for computations, you only really have to go over the \ref how section.


\page math From Math to Code
\tableofcontents
\section Briefly

There are three fundamental ideas that make our code work:

* Homology of chain complexes of free \f$\mathbb Z\f$ modules can be algorithmically computed by turning the differentials into matrices and then diagonalizing them (Smith Normal Form). 

* Free Mackey functors are determined by their bottom level groups. The higher levels are obtained by taking fixed points, and this can be done algorithmically on our matrices as long as equivariant bases are used (this is an ordered basis in which elements in the same orbit are written consecutively).

* Box products of Mackey functors are tensor products on the bottom level. This is not true for the higher levels, however using the bulletpoint above, higher levels are obtained by transferring. This is how box products of chain complexes of free Mackey functors are computed.

These three ideas dictate our approach: We work primarily on the bottom level, before transferring to get the higher ones. Homology is computed in each level separately and transfers/restrictions/Weyl group actions on the generators can be computed on the chains. But the chains are always free, so we can also algorithmically compute the effect of these operations.

By converting all our differentials to matrices, using equivariant bases throughout, we can reduce our computations to pure linear algebra (over the integers), avoiding any symbolic math.


\section imd In more detail

Here's how the code works in more detail (for simplicity we specialize to the \f$G=C_4\f$ case, although everything works equally well with general prime powers):

\subsection add The additive structure

* The input are the bottom levels of the chains of the spheres \f$S^{n\sigma+m\lambda}\f$ for \f$n,m\ge 0\f$. We could do away with only \f$S^{\sigma},S^{\lambda}\f$, but this would result in taking arbitrarily many box products and devastate run-time performance. Instead, if we use the spheres  \f$S^{n\sigma+m\lambda}\f$ for \f$n,m\ge 0\f$ we only have to take double box products at worst, and that's only for part of the multiplicative structure.

* The data of a Chain complex are the ranks and differentials. The differentials are stored as matrices, but the ranks are stored as integer arrays and not integers. This is crucial as for example \f$\mathbb Z[C_4]\f$ transfers completely differently from \f$\mathbb Z[C_2]\oplus \mathbb Z[C_2]\f$ even though they both have rank \f$4=2+2\f$ over \f$\mathbb Z\f$. 
With our conventions, \f$\mathbb Z[C_2]\oplus \mathbb Z[C_2]\f$ has rank \f$[2,2]\f$ while \f$\mathbb Z[C_4]\f$ has rank \f$4\f$.

* We transfer both ranks and differentials to higher levels. While transferring ranks is straightforward, transferring differentials is quite a bit more complicated and requires to have already transferred the ranks of the domain and range of the differentials.

* Using the classical homology algorithm we compute the groups of the Mackey functor at every level. We also compute their generators (as elements in the Chain complex)

* We transfer/restrict and compute the Weyl group action on the group generators. This concludes the Mackey functor computation for the  \f$S^{n\sigma+m\lambda}\f$, \f$n,m\ge 0\f$.

* To obtain the chains of the rest of the spheres, we box the Chains we already have. Boxing is more complicated compared to just taking tensor products, as we have to use equivariant bases throughout to transfer properly. However the most convenient bases for tensoring are not equivariant, and in the end we have to change bases through permutation matrices.

* We then perform the same procedure with transferring to get the entrire \f$RO(G)\f$ homology.

\subsection mult The multiplicative structure

Once we have the additive structure, we can work on multiplying the additive generators. 

* First we restrict the generators to the bottom level.

* We then take the product of their restrictions as an element of the box product of chain complexes at bottom level.

* The product of restrictions is a restriction and as we are working with free Mackey functors, restriction is an injection. By inverting it we can get the product of generators at a higher level as element an of the box product.

* We finally take homology and write that product in terms of the generators.

Once we know how to multiply any two additive generators, we have in effect determined the multiplicative structure (see below for a catch).

\subsection factor Factorization

Even if we can multiply any two generators, that doesn't mean we can automatically write any element as a product of our preferred generators. If we know the expression then we can easily check it, however factorizing is a lot more complicated:

* First we form a multiplication table, where all generators (in a range of course) are multiplied with the "basic irreducibles". These can be the Euler and orientation classes. 

* Once we have that we can draw a directed colored graph by connecting \f$a\f$ with \f$ab\f$ for \f$b\f$ a basic irreducible; we color this edge red. If multiplication by b is an isomorphism i.e. \f$a=(ab)/b\f$ then we also connect \f$ab\f$ with \f$a\f$; we color this edge blue.

* Since the product \f$ab\f$ may not be a generator, but rather a multiple of it, we need to allow multiples of generators as distinct nodes. We never allow trivial (0) multiples of generators.

* To obtain a factorization, we simply need to connect 1 with any node in the graph. For the most efficient factorizations, we want to minimize the alternation of blue and red edges. This is done by a generalized Dikjstra algorithm.

* For the generators not connected to 1 (eg \f$s_3\f$) we perform the same process using different sources for our graph.

\section caveat A caveat
\subsection cyclic Cyclic Generators

* The way we prove that say a transfer map is multiplication by \f$2\f$, is by computing the generators at the domain and target, compute the transfer of the domain generator and compare with the target. Of course, there are usually multiple choices of generators, but up to isomorphism we get the same Mackey functor. 

* There is a caveat however that appears when computing the multiplicative structure: If we prove that \f$ab\f$ and \f$cd\f$ are both generators of the same cyclic group, then we can't conclude that they are equal. Eg if the group is \f$\mathbb Z/4\f$ or \f$\mathbb Z\f$ then they differ by a sign. Still, since we are interested in generating the \f$RO(G)\f$ homology, as opposed to finding exact relations, we don't have to distinguish between cyclic generators and we don't need to stress over this detail.

* If we are interested in exact relations, then we can resolve the multiple generator ambiguity as follows: \f$ab\f$ and \f$cd\f$ are produced by tensoring different complexes, and if we have an explicit chain homotopy between them then we can compare directly. For example if they are obtained by tensoring \f$C_*(S^{n\sigma+m\lambda})\otimes C_*S^{\lambda}\f$ and \f$C_*(S^{n\sigma+(m-1)\lambda})\otimes C_*S^{2\lambda}\f$ then we can compare them by using \f$C_*(S^{n\sigma+(m-1)\lambda})\otimes C_*S^{\lambda}\otimes C_*S^{\lambda}\f$ as a stepping stone. 

* The problem with the above approach is that we need to take more and more box products, which is the most costly operation for runtime.


\subsection noncycl Non cyclic generators

* There is a situtation where this caveat cannot be sidestepped and that's when we have non cyclic groups. Here's an example where this problem comes up: If we have \f$\mathbb Z\oplus \mathbb Z/2\f$ with generators \f$x,y\f$ respectively then we can't automatically distinguish \f$x\f$ from \f$x+y\f$ as there is an automorphism of \f$\mathbb Z\oplus \mathbb Z/2\f$ exchanging them. In that case the difference between \f$ab\f$ and \f$cd\f$ generating the same group can be much more severe than multiplication with an integer coprime to the group's order (or a sign).

* One way out of this is to apply the approach for cyclic generators, breaking down our box products further, until they can be compared.

* Alternatively (and this is the approach we take in practice) is to ignore these products and make no statement as to the equality of \f$ab\f$ and \f$cd\f$. This gives us less data to work with, but at least in the \f$C_4\f$ case this is enough to write the factorization of any element. 



\page use How to Use
\tableofcontents
\section how Step 0: Setting the Group Parameters

* For every group there are certain mandatory parameters that need to be set for the library to work. We have included an example (for \f$G=C_4\f$) on how to set them in the file Implementation.h available in the <a href="https://github.com/NickG-Math/Mackey/tree/master/Demo">Demo</a> folder. These parameters all live in the \ref GroupSpecific "GroupSpecific" namespace and we will go over them in more detail below.

* There are also some optional ones that allow you to identify and print the names of the computed Mackey functors. This functionality is disabled by default, but can easily be turned on by defining the macro ```MACKEY_NAMES``` and setting the optional parameters. The way this is done is explained in the file Optional_Implementation.h in <a href="https://github.com/NickG-Math/Mackey/tree/master/Demo">Demo</a>. The implementation there is for \f$G=C_4\f$ and the Mackey functors in the \f$RO(C_4)\f$ homology.

\subsection var Global variables

The global variables that need to be set are:

* \ref GroupSpecific::Variables::prime "prime" : the \f$p\f$ in \f$G=C_{p^n}\f$.
* \ref GroupSpecific::Variables::power "power": the \f$n\f$ in \f$G=C_{p^n}\f$.
* \ref GroupSpecific::Variables::reps "reps" : the number of nontrivial irreducible real representations of \f$G\f$.
* \ref GroupSpecific::Variables::sphere_dimensions "sphere_dimensions" : the array consisting of the dimensions of those representations (so we must fix an order for them beforehand).

\subsection fun The standard differentials

* There is one function that needs to be manually defined, and that's \ref GroupSpecific::Function::StandardDiff "StandardDiff" creating the differential for the standard chains. In practice what this means is setting the correct matrix given the input sphere. Apart from the case work that comes from math, I have made the construction of the differential as painless as possible. 

* This construction is done through the \ref Mackey::altmatrix "altmatrix" function, that creates alternating matrices of the desired size and the desired "pattern". This pattern is repeated cyclically in the columns of the matrix. 

* An example: The matrix of size 4x4 with pattern \f$a,b\f$ is <br>
\f$\begin{matrix} a&b&a&b\\ b&a&b&a\\ a&b&a&b \\  b&a&b&a \end{matrix}\f$ <br> 
If we use the pattern \f$a,b,c,d\f$ instead we get <br>
\f$\begin{matrix} a&b&c&d\\ b&c&d&a\\ c&d&a&b \\  d&a&b&c \end{matrix}\f$ <br> 

* The matrices appearing in the Standard Chains all look like this due to equivariance.

\section next Calling the library

Once Step 0 is complete, you can include ```<Mackey/Compute.h>``` to access the methods relating to the additive and multiplicative structure, and ```<Mackey/Factorization.h>``` to access the factorization methods. For a demonstration you can use the files included in the <a href="https://github.com/NickG-Math/Mackey/Demo">Demo</a> folder.


\subsection step1add The additive structure

The file ```<Mackey/Compute.h>``` exposes the method \ref Mackey::ROHomology "ROHomology" that computes the homology of a given sphere as a Mackey functor. Example: The code

<CODE> auto M= ROHomology<rank_t,diff_t>({2,-2}); </CODE>

sets 

\f$ M=H_*(S^{2\sigma-2\lambda})\f$

Here the typenames ```rank_t,diff_t``` can be set to ```Eigen::Matrix<char,1,-1>``` and ```Eigen::Matrix<char,-1,-1>``` respectively for maximum performance, as long as the order of the group is \f$ <127 \f$. 

Here ```M``` is an object of class \ref Mackey::MackeyFunctor "MackeyFunctor" so you should read the documentation of that on how to extract that information. If the optional parameters are set then it's also possible to extract the name of the Mackey functor as seen in the "C4Verify.h" file in <a href="https://github.com/NickG-Math/Mackey/tree/master/Demo">Demo</a>


\subsection step1mult The multiplicative structure

The file ```<Mackey/Compute.h>``` also exposes the method \ref Mackey::ROGreen "ROGreen" that multiplies two generators in the Green functor \f$H_{\star}(S)\f$. Example: The code

```auto linear_combination= ROGreen<rank_t,diff_t>(2,{0,2,-2},{1,3,-4},0,0);```

multiplies the generators of

\f$ H_0^{C_4}(S^{2\sigma-2\lambda}) \otimes H_1^{C_4}(S^{3\sigma-4\lambda}) \to H_1^{C_4}(S^{5\sigma-6\lambda})  \f$

writing the answer as a linear combination of the generators in the box product. The first argument of \ref Mackey::ROGreen "ROGreen" indicates the level the generators live in (level 0=bottom, level 1= one higher etc.) so for \f$C_4\f$, level=2 is the top level. The second and third entries are the degrees of the two generators, while the last two are needed if we have noncyclic groups. Then \f$1,2\f$ selects the second and third generators of these noncyclic groups respectively (remember that in C++ counting starts from \f$0\f$.

The result of the computation ```linear_combination``` is an Eigen array (```rank_t```) that contains the coefficients eg if it's ```[2,1]``` then the product of generators is 2 times the first generator plus 1 times the second. For convenience we omit any signs etc. and identify generators of the same cyclic groups (see \ref caveat).


\subsection step1fact Factorization

The file ```<Mackey/Factorization.h>``` also exposes the class \ref Mackey::Factorization "Factorization" whose constructor creates the multiplication graph and factorizes all generators using the given sources. First construct as:

<CODE>auto F= Factorization<rank_t, diff_t> F({ -5,-5 }, { 5,5 }, { {0,1,0},{2,2,0},{0,0,1},{2,0,1} }, { "asigma", "u2sigma", "alambda", "ulambda" });</CODE>

This will work to factorize within the range of \f$S^{-5\sigma-5\lambda}\f$ to \f$S^{5\sigma+5\lambda}\f$ by multiplying with the basic irreducibles \f$ a_{\sigma}, u_{2\sigma}, a_{\lambda}, u_{\lambda}\f$ that live in degrees \f$[0,1,0],[2,2,0],[0,0,1],[2,0,1]\f$ respectively.

After that, to actually get the factorizations use

<CODE>F.compute_with_sources({[0,0,0]}, {"1"});</CODE>

and <CODE>std::cout<< F.getname(i) </CODE>

to print the name of the ```i```-th generator. This name will be nonempty as long as it's created by multiplying/dividing the basic irreducibles with ```1```. As such, it will fail for say ```s_3```. To improve it use instead

<CODE>F.compute_with_sources({[0,0,0],[-3,0,-2]}, {"1","s3"});</CODE>

For more details see the code in TestFactorization.cpp of the <a href="https://github.com/NickG-Math/Mackey/tree/master/Demo">Demo</a> folder


\page algo Algorithm Details
\tableofcontents
\section smith Smith Normal Form

For the Smith Normal Form we use a variant of the classical row-column elimination algorithm. Interestingly, for our matrices, an entry divides or is divided by any other. We optimize for this by not finding the minimum elements of the matrix, and instead working with the first nonzero element in each row and column. This is much faster for our matrices, but slower and much more dangerous for random matrices: the Smith normal form coefficient matrices \f$P,Q\f$ can easily overflow in that case (while the usual row-column elimination algorithm can give accurate results).

\section cob Change of Basis

If we have bases for modules \f$A,B\f$ then \f$A\otimes B\f$ can be canonically given two lexicographical bases, that we call left and right convenient bases. 
But if \f$A,B\f$ are the bottom levels of free Mackey functors then they have equivariant bases and the tensor product also gets an equivariant basis, called the canonical one. 
The left and right convenient bases are used to write the left and right differentials in a simple manner (hence their designation as convenient). The canonical bases are used to transfer.
\section box Box product

Computing the tensor product of Chain complexes breaks down to computing the left and right differentials \f$L(x\otimes y)=dx\otimes y\f$ and \f$R(x\otimes y)=(-1)^{|x|}x\otimes dy\f$ respectively. If we use the convenient bases explained in the previous section, these are just block diagonal matrices with the blocks being the differential from the original chains \f$d\f$. To get \f$L,R\f$ w.r.t. the canonical bases, we need to apply the change of basis matrices explained above. Once we do that the total differential of the tensor product is just \f$L+R\f$. To be more accurate, \f$L\f$ is not a single differential, but rather a sequence of them, one for each summand of the Box product. That is, \f$(C\otimes D)_n\to (C\otimes D)_{n-1}\f$ is a map \f$C_n\otimes D_0\oplus\cdots \oplus C_0\otimes D_n\to C_0\otimes D_{n-1}\oplus\cdots\oplus C_{n-1}\otimes D_0\f$. Each summand of this map is computed separately into an \f$L\f$ and an \f$R\f$, and then these are mixed together to form the total differential. The mixing specifies that we start with a block \f$ L_0\f$, then place \f$ R_0\f$ directly below it, then \f$L_1\f$ adjecent to the right of \f$ R_0\f$ etc.

\section graph Graphs

* For weighted graphs we want the shortest path from a given source to all other points. I use a straightforward implementation of Dikjstra's algorithm using std::priority_queue.

* For graphs of two colors, we are interested in the paths from the source to all points with the minimum alternations of colors (an alternation of colors means switching from division to multiplication and vice-versa). 
This problem can be easily reduced to the previous bullet, by using a sort of "dual" graph where now the nodes are colored and edges between same colored nodes have weight 0, while for different colored nodes we get weight 1.
To get the new graph simply duplicate the nodes of the original, color the originals by red and the new ones by blue, and quadruple the edges (so we using all combinations of colored nodes).
After that we can find the red and blue paths starting and ending from a red/blue source and compare them in length, choosing the shortest one.


\page perf Performance
\tableofcontents

First, two observations:

* The Linux binary runs measurably faster than the Windows one.
* Out of all compilers, Clang seems to produce marginally faster code. 

\section compoptions Compiler Options


* I recommend the following compiler options (GCC, Clang): <CODE>-Ofast --funroll-loops -march=native </CODE>
* With the Intel and Microsoft compilers Eigen recommends <CODE>-inline-forceinline</CODE> option.
* Note: ```-Ofast``` doesn't actually reduce the accuracy of our results, since we only use integers even if the data-type is sometimes a floating point (see \ref intvsfloat for an explanation as to why we do that). 
\section thread Multithreading

* While the <a href="https://github.com/NickG-Math/Mackey/tree/master/bin">binaries</a> are all single-threaded, the most (by far) computationally intensive calculations can be multithreaded extemelly easily and efficiently. This is as simple as adding a ```#pragma omp parallel for``` before certain loops that compute the additive/multiplicative structure in a range. There is no need to lock anything.

* There is one caveat: While the loop iterations are independent, they are not all equally intensive. A sphere like \f$S^{2\sigma+\lambda}\f$ is cheaper to compute compared to \f$S^{6\sigma+8\lambda}\f$ which is in turn much cheaper compared to \f$S^{6\sigma-8\lambda}\f$ as the latter one involves a box product. In the multiplicative structure we may have to take double box products, and these are by comparison much more expensive in run-time as they involve arbitrarily large permutation matrices.

* So it's important to equally divide the work amongst the thread. At this point, this has to be done manually.

\section intvsfloat Integers vs Floats

I use integers (or indeed ```char``` and ```short```) for the majority of the computations; that's usually the fastest method and makes the most sense (as all numbers appearing are actually integers).
There is one important exception: Matrix multiplication. Eigen is much slower with integer matrix multiplication compared to floating points, and the Intel MKL does not even support integer matrix multiplication.
So when we need to multiply matrices we cast them to floats. This is only needed for the Homology algorithm, which is at the very end of the pipeline (together with the Smith Normal Form) so we can benefit from smaller integer types before casting.

\section memo Memoizing ChangeBasis

To form the Box product of Chains we need the change of basis matrices. These matrices only depend on the ranks of the given Chains, call them rank1 and rank2. 
The ranks that actually come up in our computations always look like [?,order,...,order,?] where order is the order of the group and ?\f$\le\f$ order.
This means that we very effectively memoize this function for better performance.


\section bottle Bottlenecks

* The biggest performance bottleneck is found in the multiplicative structure. That's when we apply some large change of basis matrices through Eigen's permutation matrix product. This is a memory bottleneck.
* The second biggest bottleneck lies in the transfering very large differentials. To transfer we need to delete certain rows of the matrix, and this is done by copying the remaining rows into a new matrix.  This is a memory bottleneck.
* A more minor bottleneck is the Smith normal form computation. The problem is that it involves going through our matrix both by rows and by columns, which is not ideal for cache locality. The SNF is both memory and core bound.
