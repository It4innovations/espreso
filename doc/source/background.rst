



.. .. math::
.. 
..     \begin{array}{rl}
..       CP				& \text{Coarse Problem }\\
..       DD       			& \text{Domain Decomposition }\\
..       DOF      		& \text{Degrees of Freedom }\\
..       FETI     			&\text{ Finite Tearing and Interconnecting DD method }\\
..       TFETI    		&\text{ Total FETI method }\\
..       HTFETI   		&\text{ Hybrid Total FETI  method}\\ 
..       SpMV			& \text{Sparse Matrix Vector Multiplication }\\
..       RHS			&\text{ Right-Hand Side }\\
..       RBM        &\text{ Rigid Body Motions }\\
..       SpDS			&\text{ Sparse Direct Solver }\\
..       \mathbf{u}_i 		&\text{ subdomain nodal displacements }\\ 
..       \boldsymbol{\lambda} 		&\text{ Lagrange multipliers }\\ 
..       \mathbf{K}_{i}          	& \text{non-regularized subdomain stiffness  matrix }\\
..       \hat{\mathbf{K}}_{i}     	& \text{regularized subdomain stiffness  matrix }\\
..       \mathbf{T}_{i}           	&\text{ subdomain regularization  matrix }\\
..       \mathbf{M}_{i}           	& \text{subdomain preconditioner matrix }\\
..       \mathbf{R}_{i}           	& \text{subdomain kernel matrix }\\
..       \mathbf{B}_{1,i}    	& \text{subdomain gluing matrix - constraint matrix  }\\
..       \mathbf{B}_{0,i}       	& \text{subdomain corner matrix }\\
..       \mathbf{F}_0			& \text{cluster FETI operator }\\
..       \mathbf{G}_0			& \text{cluster matrix }\\
..       \mathbf{S}_0			& \text{cluster coarse problem matrix }\\
..       \mathbf{GG}^\top		& \text{coarse problem matrix of the TFETI method }\\
..     \end{array}


Domain Decomposition
====================

For more than two decades, the Finite Element Tearing and
Interconnecting method (FETI, for more detail see, e.g., in [5]_, [7]_) has been
successfully used in the engineering community, primarily for very large
problems arising from the discretization of partial differential
equations. In such an approach the original structure is artificially
decomposed into several non-overlapping subdomains. Mutual continuity of
primal variables between neighboring subdomains is enforced afterwards
by dual variables, i.e., Lagrange multipliers (LM). They are usually
obtained iteratively by some of the Krylov subspace methods, then the
primal solution is evaluated locally on each subdomain. The FETI method
is also applicable for contact problems with inequality constrains (see,
e.g., [3]_, [4]_)

In 2006 Dostál et al. [2]_ introduced a new variant of an algorithm called
Total FETI (or TFETI) in which Dirichlet boundary condition is enforced
also by LM. Please note, there is no difference, if the hybrid variant
stems from FETI or TFETI. The impact is identical in both cases. This
paper addresses the modification of the TFETI variant, and is therefore
called HTFETI.

The HTFETI method is a variant of hybrid FETI methods introduced by
Klawonn and Rheinbach (see, e.g., [8]_, [9]_) for FETI and FETI-DP. In the
original approach, the authors gather a number of subdomains in clusters
which provides a three-level domain decomposition approach. Each cluster
consists of a certain number of subdomains and for these, a FETI-DP
system is set up. The clusters are then solved by a traditional FETI
approach using projections to treat the non trivial kernels. In
contrast, in HTFETI, a TFETI approach is used for the subdomains in each
cluster and the FETI approach using projections for clusters.

The main advantage of HTFETI is the ability to solve problems decomposed
into a very large number of subdomains. 

.. _label_fig1:
.. figure:: _images/FETI-motiv.*
   :width: 70%
   :align: center
   
   The two approaches how to use HTFETI method and their comparison with TFETI for problems of size up to 3 billion DOF.




To solve extremely large problems, the configuration of the HTFETI has
to be different (*Approach A: 1 cluster per compute node with large
number of subdomains - up to several thousand*) than for solving
moderate problems (*Approach B: 1 cluster per CPU core with small number
of subdomains - up to a few hundred*). In particular, Approach B can use
subdomains of approximately 1,000 DOF to solve a 1 billion DOF problem
on approximately 200 compute nodes with 64 GB of RAM. Solving the same
problem with TFETI requires subdomains of approximately 90,000 DOF.

The comparison of HTFETI and TFETI methods is shown in :numref:`label_fig1`. 
Using TFETI with large number of DOF per
subdomain is an example of how to solve large problems using this
technique. The figure shows that its scalability is not ideal, but is
acceptable up to 3 billion DOF. However, the solution time in this range
is still shorter than it would be with the use of the HTFETI technique.
The figure clearly shows that for problems greater than 3 billion DOF,
the HTFETI Approach A is better.

However, the HTFETI can be configured to the other extreme in which
small number of subdomains per cluster (here 216) and very small number
of DOF per subdomain (here 1536 DOF) are used. This approach is
significantly faster than TFETI with large number of DOF per subdomain.
Figure also shows that TFETI can be configured with smaller number of
DOF per subdomain (here 20,577) to reduce processing time. However, this
option is suitable only for very small problems as its scalability is
unsatisfactory.

In sum, the figure shows that HTFETI provides better solution time for
“smaller problems” up to 3 billion DOF and better scalability for large
problems.

In the next section, the theory behind the HTFETI is described, followed
by the description of the parallel implementation in the ESPRESO
library.

Hybrid Total FETI Method 
========================

The FETI method, widely known for more than two deca-des, is an
efficient tool for solving large-scale problems in structural mechanics
via Krylov subspace methods adapted to parallel machines. Although this
paper is focused on the HTFETI, this section will first introduce the
original FETI method on a simple cantilever beam 

.. _label_fig2:   
.. figure:: _images/problem_discretization.*
   :width: 70%
   :align: center

   `Cantilever beam, FEM discretization`

followed by an introduction of its hybrid variant.


.. _label_FETI:

FETI
----

In the simple engineering problem depicted in 
:numref:`label_fig2`, the aim is to get deformation of
the structure. A cantilever beam, fixed on the left side and loaded by
the pressure :math:`p_y` at the top, is discretizated with FEM [11]_. The
number of nodes is :math:`n_n^0`, and number of elements is :math:`n_e`.
Above the mesh global stiffness matrix :math:`{\mathbf{K}}_g` and global
RHS :math:`{\mathbf{f}}_g` (nodal forces) are assembled, whereas both
objects are included into nodal equilibrium equation

.. math:: 
  :label: label_Ku_f_undecomp
  
  {\mathbf{K}}_g{\mathbf{u}}_g={\mathbf{f}}_g.

To get the vector :math:`{\mathbf{u}}_g` (nodal displacements),
Dirichlet boundary condition has to be taken into account due to a
singular matrix :math:`{\mathbf{K}}_g` since it has a non-empty
null-space. Then the linear system can be solved by an iterative or
direct solver.

Clearly, this undecomposed approach has its own limitations. The size of
:math:`{\mathbf{K}}_g` can overload the computer memory which is one of
the reasons one might employ domain decomposition methods. In the FETI
case, the mesh is decomposed into four smaller submeshes (see :numref:`label_fig3`-a) 

.. _label_fig3:
.. figure:: _images/decomp_feti_hfeti.*
   :width: 70%
   :align: center
  
   `Domain decomposition`

.. :ref:`label_fig2`
to avoid assembling the global objects
:math:`\mathbf{K}_g`, :math:`\mathbf{f}_g`. Decomposition generally
causes an increase in the number of nodes in cuts between subdomains,
thus the global size of the unknowns of the decomposed problem is always
bigger than the original one. In fact, it is not a problem, because any
object in a FETI algorithm stored at one computational node does not
have such dimension. Equilibrium equation for :math:`i`-th subdomain is

.. math::
  :label:  label_eq_1

   \mathbf{K}_i \mathbf{u}_i = \mathbf{f}_i-\mathbf{B}^T_i \boldsymbol{\lambda},~ i\in (1,2,3,4),
   

where :math:`\mathbf{K}_i` and :math:`\mathbf{f}_i` with the same
meaning as before are assembled for all subdomains separately. On the
RHS beside the vector :math:`\mathbf{f}_i` the product
:math:`\mathbf{B}^T_i \boldsymbol{\lambda}` appears. It is actually an
additive nodal force vector which acts on the interface between
subdomains and it substitutes influence transmitting from neighbors.
Those four systems in :eq:`label_eq_1` are supplemented by the global constraint
condition

.. math::
  :label: label_eq_Bu

   \sum_{i=1}^{4} \mathbf{B}_i\mathbf{u}_i=\mathbf{o}. 
     \label{eq:sum_Biui}

Equilibrium equations :eq:`label_eq_1` together with condition :eq:`label_eq_Bu` can be written globally as

.. math:: 
  :label: label_eq:KKTprimal_detailed

     \left(
       \begin{array}{lllll}
       \mathbf{K}_1 & \mathbf{O}         & \mathbf{O}         &
       \mathbf{O}   & \mathbf{B}_{1}^T  \\
       \mathbf{O}   & \mathbf{K}_2       & \mathbf{O}         & 
       \mathbf{O}   & \mathbf{B}_{2}^T \\
       \mathbf{O}   & \mathbf{O}         & \mathbf{K}_3       & 
       \mathbf{O}   & \mathbf{B}_{3}^T \\
       \mathbf{O}   & \mathbf{O}         & \mathbf{O}         & 
       \mathbf{K}_4 & \mathbf{B}_{4}^T\\
      \mathbf{B}_{1}& \mathbf{B}_{2}&  
      \mathbf{B}_{3} & \mathbf{B}_{4}& \mathbf{O}
     \end{array}
   \right)
   \left(\begin{array}{c}
       \mathbf{u}_1 \\ \mathbf{u}_2 \\ \mathbf{u}_3 \\ \mathbf{u}_4 \\ 
       \boldsymbol{\lambda}
   \end{array}\right)=
   \left(\begin{array}{c}
       \mathbf{f}_1 \\ \mathbf{f}_2 \\ \mathbf{f}_3 \\ \mathbf{f}_4 \\
       \mathbf{c}
   \end{array}\right)

or shortly

.. math::
  :label:  label_eq_KKT1

   \begin{array}{l}
       \mathbf{K}\mathbf{u} + \mathbf{B}^T\boldsymbol{\lambda} = \mathbf{f}, \\
       \mathbf{B}\mathbf{u} = \mathbf{c}.
     \end{array}

The meaning of symbols in :eq:`label_eq:KKTprimal_detailed` is obvious to detailed expression in
:eq:`label_eq_KKT1`. 
Generally, vector :math:`\mathbf{c}` contains zero entries. Vector of LM

.. math::

   \boldsymbol{\lambda} = \left(
           (\boldsymbol{\lambda}^{d,1})^T,~
           (\boldsymbol{\lambda}^{g,12})^T,~
           (\boldsymbol{\lambda}^{g,23})^T,~
           (\boldsymbol{\lambda}^{g,34})^T 
                     \right)^T

consists of four parts where the first subvector
:math:`\boldsymbol{\lambda}^{d,1}` enforces Dirichlet boundary
condition, the second one :math:`\boldsymbol{\lambda}^{g,12}` enforces
connectivity between subdomains 1 and 2, etc. The division of global
:math:`\boldsymbol{\lambda}` into four parts does not relate to the
number of subdomains. These four sets correspond to four cuts caused by
the decomposition. The splitting into four subvectors shows how
:math:`\boldsymbol{\lambda}` can be stored according to numeric
superindexes. For instance, the first part
:math:`\boldsymbol{\lambda}^{d,1}` is related to the first subdomain
only, and it will never touch another one, therefore it is stored only
on the first subdomain. The second part
:math:`\boldsymbol{\lambda}^{g,12}` connects subdomain 1 and 2, thus
they are stored only there. Evidently, there is no need to assemble the
global :math:`\boldsymbol{\lambda}` vector on a single node.

In the next step, the primal variables will be eliminated. From 
Eq. :eq:`label_eq_KKT1` the vector of unknown displacements is

.. math::
  :label: label_eq:u_prim

   \mathbf{u} = \mathbf{K}^{+} (\mathbf{f}-\mathbf{B}^T \boldsymbol{\lambda})
      + \mathbf{R} \boldsymbol{\alpha},
..      \label{eq:u_prim}

where
:math:`\mathbf{R}=\text{diag}\left(\mathbf{R}_1,~\mathbf{R}_2,~\mathbf{R}_3,~\mathbf{R}_4\right)`
is block-diagonal matrix, the columns of which define the basis of null
space of :math:`\mathbf{K}`, and the vector :math:`\boldsymbol{\alpha}`
contains its amplitudes. Mutual relation between them is

.. math::
  :label: label_eq:KR_O

   \mathbf{K}\mathbf{R}=\mathbf{O}~\text{or}~\mathbf{R}^T\mathbf{K}=\mathbf{O}.
     \label{eq:KR_O}

The term :eq:`label_eq:u_prim` with the second condition in :eq:`label_eq_KKT1` eliminates the primal variables

.. math::

   \mathbf{B}\mathbf{K}^{+} \mathbf{f}-\mathbf{B}\mathbf{K}^{+} \mathbf{B}^T \boldsymbol{\lambda}
     + \mathbf{B}\mathbf{R} \boldsymbol{\alpha}=\mathbf{c}
      \label{eq:KKTdual_a}

similarly, as a combination of the first equation from :eq:`label_eq_KKT1` and Eq :eq:`label_eq:KR_O`.

.. math::

   \mathbf{R}^T\mathbf{Ku}+\mathbf{R}^T\mathbf{B}^T\boldsymbol{\lambda} = \mathbf{R}^T\mathbf{f}.
      \label{eq:KKTdual_b}

The last two terms together can be expressed as

.. math::
  :label: label_eq:KKTdual

   \begin{array}{l}
     \mathbf{F}\boldsymbol{\lambda} + \mathbf{G}^T\boldsymbol{\alpha} = \mathbf{d}, \\
     \mathbf{G}\boldsymbol{\lambda} = \mathbf{e}
   \end{array}

in which the following substitutions are used

.. math::

   \begin{array}{ll}
       \mathbf{F}=\mathbf{BK}^+\mathbf{B}^T-\mathbf{c}, & \mathbf{G}=-\mathbf{BR}^T, \\
       \mathbf{d}=\mathbf{BK}^+\mathbf{f}, & \mathbf{e}=-\mathbf{R}^T \mathbf{f}.
     \end{array}


.. _label_pcgm:

Projected Conjugate Gradient Method
-----------------------------------

In the previous steps, the primal variables :math:`\mathbf{u}` were
eliminated. The newly derived system of linear equations :eq:`label_eq:KKTdual` can be
favorably solved by iterative methods, e.g., the conjugate gradient
method (CGM).

.. math::

     \begin{array}{rl} 
       1. &  \text{set:} ~ \varepsilon>0,~i_{max}>0, ~\boldsymbol{\lambda}_0,\\
       2. & \overline{\mathbf{g}_0} = \mathbf{g}_0 +\mathbf{G}^T \boldsymbol{\alpha}_0 ~~ \text{where} ~~ 
      \mathbf{g}_0 = \mathbf{F}\boldsymbol{\lambda}_0 -\mathbf{d}\\
      3. & \mathbf{w}_0 = \overline{\mathbf{g}_0} \\
      4. & \mathbf{do}~i=0,1,...,i_{max} \\
      5. &  \hspace{5mm} \rho_i = -(\overline{\mathbf{g}_{i}},\overline{\mathbf{g}_{i}})/
       (\mathbf{w}_i,\mathbf{Fw}_i) \\
      6. & \hspace{5mm}\boldsymbol{\lambda}_{i+1} = \boldsymbol{\lambda}_{i+1} + \mathbf{w}_i \rho_i \\
      7. & \hspace{5mm}\overline{\mathbf{g}_{i+1}}=\mathbf{g}_{i+1}+\mathbf{G}^T\boldsymbol{\alpha}_i,~~\mathbf{g}_{i+1} = \mathbf{g}_{i} + \mathbf{Fw}_i \rho_i\\
      8. & \hspace{5mm}\mathbf{if}~\|\overline{\mathbf{g}_{i+1}}\|<\varepsilon  \\
      9. & \hspace{10mm} \mathbf{break} \\
      10. & \hspace{5mm} \mathbf{end if} \\
      11.&\hspace{5mm} \boldsymbol{\gamma} = (\overline{\mathbf{g}_{i+1}},\overline{\mathbf{g}_{i+1}})
       / (\overline{\mathbf{g}}_{i},\overline{\mathbf{g}_{i}})\\ 
      11. & \hspace{5mm} \mathbf{w}_{i+1} = \overline{\mathbf{g}_{i+1}} + \mathbf{w}_i \gamma_i\\
      12. & \mathbf{end}~\mathbf{do} \\
       & \\
       &  \hspace{15mm}\mathbf{Algorithm}~1.
     \end{array}

An approximation of LM used in Algorithm :math:`1` in :math:`i`-th
iteration is considered in the form

.. math::
  :label: label_eq:lambda_iter

   \boldsymbol{\lambda}_i = \boldsymbol{\lambda}_o + 
   \sum_{j=0}^{i}\mathbf{w}_j \rho_j.
   \label{eq:lambda_iter}

If an initial guess is chosen as a linear combination of basis vectors
of matrix :math:`\mathbf{G}^T` in the form of
:math:`\boldsymbol{\lambda}_0=\mathbf{G}^T \mathbf{y}`, the vector
:math:`\mathbf{y}` can be simply determined using the second equation in
:eq:`label_eq:KKTdual`
as

.. math:: \boldsymbol{\lambda}_0 = (\mathbf{G}\mathbf{G}^T)^{-1}\mathbf{e}.

Such initial guess fully accomplishes the second equality in :eq:`label_eq:KKTdual`,
therefore the rest of the approximation in :eq:`label_eq:lambda_iter` lies in the kernel of
:math:`\mathbf{G}`.

During the iterative process, just the partial gradient
:math:`\mathbf{g}_i` is kept in the memory, and before using the
complete form :math:`\overline{\mathbf{g}_{i+1}}` it is modified by
:math:`\boldsymbol{\alpha}_i` to satisfy condition
:math:`\mathbf{G}\overline{\mathbf{g}_{i+1}}=\mathbf{o}`. Firstly, it
appears if the initial conjugate vector :math:`\mathbf{w}_0` is used as
a contribution to the :math:`\boldsymbol{\lambda}_0`. To keep
:math:`\mathbf{w}_0` in the kernel of :math:`\mathbf{G}`, vector

.. math::

   \boldsymbol{\alpha}_0 = -(\mathbf{G}\mathbf{G}^T)^{-1}
   \mathbf{G}(\mathbf{F}\boldsymbol{\lambda}_0 -\mathbf{d})

has to be evaluated to satisfy the condition

.. math::

   \mathbf{G}\mathbf{w}_0 = \mathbf{G}\overline{\mathbf{g}_0} = \mathbf{G}(\mathbf{F}\boldsymbol{\lambda}_0 -\mathbf{d}) + 
     \mathbf{G}\mathbf{G}^T \boldsymbol{\alpha}_0=\mathbf{o}.

The fulfillment of the condition above is equivalent to the projection
of the part of the gradient

.. math::

   \mathbf{w}_0=\overline{\mathbf{g}}_0= \mathbf{P}
   (\mathbf{F}\boldsymbol{\lambda}_0 -\mathbf{d})=
   \mathbf{P}\mathbf{g}_0

by the orthogonal projector

.. math:: \mathbf{P} =  \mathbf{I}-\mathbf{G}^T(\mathbf{G}\mathbf{G}^T)^{-1} \mathbf{G}.

Due to the projection applied to each iteration, Algorithm 1 is called
*Projected Conjugate Gradient Method* (PCGM) (for more details see [6]_).
The iterative process can be accelerated by Lumped or Dirichlet
preconditioner, in this case the preconditioned gradient has to be
projected again (for more detail see, e.g., [10]_).

From FETI to HTFETI
-------------------

FETI algorithm described in subsection :ref:`label_FETI` is efficient
(besides other effects related to the condition of the system etc.)
until CP

.. math:: \mathbf{GG}^T\in\mathbb{R}^{n_{cp}\times n_{cp}}

can be effectively factorized from the time point of view. Its dimension
:math:`n_{cp}` in linear elasticity depends on the number of subdomains
:math:`n_{sub}` (in 2D: :math:`n_{cp}=3\cdot n_{sub}`, in 3D:
:math:`n_{cp}=6\cdot n_{sub}`), and it is equal to all RBM in a
decomposed structure. In :numref:`label_fig3`, as
illustrated using the FETI technique, each subdomain of the decomposed
beam has three RBM (in 2D), therefore :math:`n_{cp}=3 \cdot 4 = 12`. CP
dimension reduction by HTFETI can be demonstrated even in this small and
simple benchmark. The permutation and splitting of matrix
:math:`\mathbf{B}` according to :numref:`label_fig3`-c reads

.. math::

     \left(
       \begin{array}{llll|ll|l}
       \mathbf{K}_1 & \mathbf{O}         & \mathbf{O}         &
       \mathbf{O}   & \mathbf{B}_{0,1}^T & \mathbf{O}         & \mathbf{B}_{1,1}^T  \\
       \mathbf{O}   & \mathbf{K}_2       & \mathbf{O}         & 
       \mathbf{O}   & \mathbf{B}_{0,2}^T & \mathbf{O}         & \mathbf{B}_{1,2}^T \\
       \mathbf{O}   & \mathbf{O}         & \mathbf{K}_3       & 
       \mathbf{O}   & \mathbf{O}         & \mathbf{B}_{0,3}^T & \mathbf{B}_{1,3}^T \\
       \mathbf{O}   & \mathbf{O}         & \mathbf{O}   & 
       \mathbf{K}_4 & \mathbf{O}         & \mathbf{B}_{0,4}^T & \mathbf{B}_{1,4}^T\\
       \hline
       \mathbf{B}_{0,1} & \mathbf{B}_{0,2}& 
       \mathbf{O}& \mathbf{O}& \mathbf{O}& \mathbf{O}& \mathbf{O}\\
       \mathbf{O}& \mathbf{O}&  
      \mathbf{B}_{0,3} & \mathbf{B}_{0,4}& \mathbf{O}& \mathbf{O}& \mathbf{O}\\
       \hline
      \mathbf{B}_{1,1}& \mathbf{B}_{1,2}&  
      \mathbf{B}_{1,3} & \mathbf{B}_{1,4}& \mathbf{O}& \mathbf{O}& \mathbf{O}
     \end{array}
   \right)
   \left(\begin{array}{c}
       \mathbf{u}_1 \\ \mathbf{u}_2 \\ \mathbf{u}_3 \\ \mathbf{u}_4 \\ \hline
       \boldsymbol{\lambda}_{0,1} \\ \boldsymbol{\lambda}_{0,2} \\\hline 
       \boldsymbol{\lambda}_{1}
   \end{array}\right)=
   \left(\begin{array}{c}
       \mathbf{f}_1 \\ \mathbf{f}_2 \\ \mathbf{f}_3 \\ \mathbf{f}_4 \\ \hline 
       \mathbf{c}_{0,1} \\ \mathbf{c}_{0,2} \\\hline 
       \mathbf{c}_{1}
   \end{array}\right)

\ where matrix

.. math::

   \mathbf{B}_0=
     \left(\begin{array}{llll}
       \mathbf{B}_{0,1} & \mathbf{B}_{0,2}& \mathbf{O}& \mathbf{O} \\
       \mathbf{O}& \mathbf{O} & \mathbf{B}_{0,3} & \mathbf{B}_{0,4} 
     \end{array}\right)

consists of constraints gluing subdomains into 2 clusters, and

.. math::

   \mathbf{B}_1=
     \left(\begin{array}{llll}
       \mathbf{B}_{1,1} & \mathbf{B}_{1,2}& \mathbf{B}_{1,3} & \mathbf{B}_{1,4} 
     \end{array}\right)

contains the rest of the equality constraints. This permuted system of
the linear equation still has the structure for a FETI algorithm as in 
:eq:`label_eq:KKTprimal_detailed`.
The next equation shows the system after the permutation

.. math::
   :label: label_eq:KKTprimal_reordered

     \left(
     \begin{array}{lll|lll|l}
       \mathbf{K}_1 & \mathbf{O}         & \mathbf{B}_{0,1}^T &
       \mathbf{O}   & \mathbf{O}         & \mathbf{O}         & \mathbf{B}_{1,1}^T  \\
       \mathbf{O}   & \mathbf{K}_2       & \mathbf{B}_{0,2}^T & 
       \mathbf{O}   & \mathbf{O}         & \mathbf{O}         & \mathbf{B}_{1,2}^T \\
       \mathbf{B}_{0,1}   & \mathbf{B}_{0,2}& 
       \mathbf{O}& \mathbf{O}& \mathbf{O}& \mathbf{O}& \mathbf{O}\\
       \hline
       \mathbf{O}   & \mathbf{O}         & \mathbf{O}          & 
       \mathbf{K}_3  & \mathbf{O}        & \mathbf{B}_{0,3}^T & \mathbf{B}_{1,3}^T \\
       \mathbf{O}   & \mathbf{O}         & \mathbf{O}   & 
       \mathbf{O}   & \mathbf{K}_4       & \mathbf{B}_{0,4}^T & \mathbf{B}_{1,4}^T\\
       \mathbf{O}& \mathbf{O}&  \mathbf{O}&
      \mathbf{B}_{0,3} & \mathbf{B}_{0,4}& \mathbf{O}&  \mathbf{O}\\
       \hline
      \mathbf{B}_{1,1}& \mathbf{B}_{1,2}& \mathbf{O}& 
      \mathbf{B}_{1,3}& \mathbf{B}_{1,4}&  \mathbf{O}& \mathbf{O}
     \end{array}
   \right)
   \left(\begin{array}{c}
       \mathbf{u}_1 \\ \mathbf{u}_2 \\  \boldsymbol{\lambda}_{0,1} \\ \hline
       \mathbf{u}_3 \\ \mathbf{u}_4 \\  \boldsymbol{\lambda}_{0,2} \\ \hline
       \boldsymbol{\lambda}_{1}
   \end{array}\right) =
   \left(\begin{array}{c}
       \mathbf{f}_1 \\ \mathbf{f}_2 \\ \mathbf{c}_{0,1} \\ \hline
       \mathbf{f}_3 \\ \mathbf{f}_4 \\ \mathbf{c}_{0,2} \\ \hline
       \mathbf{c}_{1}
   \end{array}\right).

This system can be rewritten in a simplified form

.. math::
  :label: label_eq:KKThprimal

     \left(
     \begin{array}{l|l|l}
       \tilde{\mathbf{K}}_{1} & \mathbf{O} & \tilde{\mathbf{B}}^T_{1} \\
       \hline
       \mathbf{O} & \tilde{\mathbf{K}}_{2} & \tilde{\mathbf{B}}^T_{2} \\
       \hline
       \tilde{\mathbf{B}}_{1} & \tilde{\mathbf{B}}_{2} & \mathbf{O}  
     \end{array}
   \right)
   \left(\begin{array}{l} \tilde{\mathbf{u}}_1 \\ \hline \tilde{\mathbf{u}}_2 \\
       \hline
           \tilde{\boldsymbol{\lambda}} 
       \end{array} \right) = 
   \left(\begin{array}{l} \tilde{\mathbf{f}}_1 \\ \hline \tilde{\mathbf{f}}_2 \\
       \hline
           \tilde{\mathbf{c}} 
       \end{array} \right)

which can be solved by HTFETI method. The classical FETI method belongs
to the group of dual Schur complement methods in which the primal
variables :math:`\mathbf{u}` are eliminated and satisfied exactly in
each iteration. Dual variables :math:`\boldsymbol{\lambda}`, reaction
forces between subdomains, are obtained by an iterative solver. In the
case of HTFETI, not only the primal variables but also a subset of dual
variables :math:`\boldsymbol{\lambda}_0` are eliminated. This
modification causes subdomains 1 and 2 to be partially glued onto
interface into cluster 1, and subdomains 3 and 4 into cluster 2. After
that, the structure behaves like a problem decomposed into two parts,
and the dimension of CP is two times smaller accordingly
(:math:`n_{cp}= 3 \cdot n_c = 3 \cdot 2=6` where :math:`n_c` is the
number of clusters). The important thing is that the linear system :eq:`label_eq:KKThprimal`
can be handled by the same methodology applied onto in subsections
:ref:`label_FETI` and :ref:`label_pcgm`, so long as each object is relabeled by
a tilde.

HTFETI - Optimal implementation
-------------------------------

In the following text the methodology will be described on the first
cluster and its subdomains 1, 2.

 To keep optimal properties of the HTFETI code, the cluster stiffness
matrix :math:`\tilde{\mathbf{K}}_1` can not be assembled, but used
implicitly. It requires modified routines for factorization and kernel
detection.
The factorization is obtained by solving
:math:`\tilde{\mathbf{K}}_1\tilde{\mathbf{x}}=\tilde{\mathbf{b}}` more
detailed written as

.. math::
   :label: eq:KKTprimal_kplus

     \left(
     \begin{array}{ll}
       \mathbf{K}_{1:2} & \mathbf{B}^T_{0,1:2} \\
       \mathbf{B}_{0,1:2} &  \mathbf{O}  
     \end{array}
   \right)
   \left(\begin{array}{l}
     \mathbf{x} \\ \boldsymbol{\mu}  
   \end{array}\right) =
   \left(\begin{array}{l}
       \mathbf{b} \\ \mathbf{z}  
   \end{array}\right)

with :math:`\mathbf{K}_{1:2}=\text{diag}(\mathbf{K}_1,~\mathbf{K}_2)`
and :math:`\mathbf{B}_{0,1:2}=(\mathbf{B}_{0,1},~\mathbf{B}_{0,2})`.
Notation :math:`1:2` points to the first and last ordinal number of the
subdomains on the cluster. So for the second cluster it would be
:math:`3:4`. The system :eq:`eq:KKTprimal_kplus` can be taken as a FETI problem to be solved not
iteratively but by a direct solver. After the dualization it reads

.. math::

     \left(
     \begin{array}{ll}
       \mathbf{F}_{0,1:2} & \mathbf{G}^T_{0,1:2} \\
       \mathbf{G}_{0,1:2} & \mathbf{O}  
     \end{array}
   \right)
   \left(\begin{array}{l}
    \boldsymbol{\mu}  \\ \boldsymbol{\beta}  
   \end{array}\right) =
   \left(\begin{array}{l}
       \mathbf{d}_{0,1:2} \\ \mathbf{e}_{0,1:2}  
   \end{array}\right)
   \label{eq:KKTdual_kplus}

with substitutions

.. math::

     \begin{array}{lr}
       \mathbf{F}_{0,1:2}=\mathbf{B}_{0,1:2}\mathbf{K}_{1:2}^+\mathbf{B}_{0,1:2}^T-\mathbf{z},
       & \mathbf{G}_{0,1:2}=-\mathbf{R}_{1:2}^T\mathbf{B}_{0,1:2}^T, \\
       \mathbf{d}_{0,1:2}=\mathbf{B}_{0,1:2}\mathbf{K}_{1:2}^+\mathbf{b}, &
       \mathbf{e}_{0,1:2}=-\mathbf{R}_{1:2}^T \mathbf{b}.
     \end{array}

To obtain the vector
:math:`\tilde{\mathbf{x}}=(\mathbf{x}^T,~\boldsymbol{\mu}^T)^T`, both
systems , are subsequently solved in three steps

.. math::

     \begin{array}{l}
       \boldsymbol{\beta} = \mathbf{S}_{0,1:2}^+
     \left(\mathbf{G}_{0,1:2}\mathbf{F}_{0,1:2}^{-1}\mathbf{g}_{0,1:2}-\mathbf{e}_{0,1:2}\right)\\
     \boldsymbol{\mu} = \mathbf{F}_{0,1:2}^{-1}
     \left(\mathbf{g}_{0,1:2}-\mathbf{G}_{0,1:2}^T\boldsymbol{\beta}\right) \\
     \mathbf{x} = \mathbf{K}_{1:2}^+
     \left(\mathbf{b}-\mathbf{B}_{0,1:2}^T\boldsymbol{\mu}\right)+
     \mathbf{R}_{1:2}\boldsymbol{\beta}
   \end{array}

in which singular Schur complement

.. math:: \mathbf{S}_{0,1:2}=\mathbf{G}_{0,1:2}\mathbf{F}_{0,1:2}^{-1}\mathbf{G}_{0,1:2}^T

appears.

In the HTFETI method, not only do the subdomain matrices
:math:`\mathbf{K}_i` have to be factorized, but also
:math:`\mathbf{F}_{0,1:2}` and :math:`\mathbf{S}_{0,1:2}` - one pair on
each cluster. The dimension of :math:`\mathbf{F}_{0,1:2}` is controlled
by the number of LM which glue subdomains of the cluster, and the
dimension of :math:`\mathbf{S}_{0,1:2}` is given by the number of
subdomains per cluster multiplied by the defect of :math:`\mathbf{K}_i`.

To get the kernel :math:`\tilde{\mathbf{R}}_1` of the first cluster, the
following term is written

.. math::

     \left(
     \begin{array}{ll|l}
       \mathbf{K}_{1} & \mathbf{O} & \mathbf{B}^T_{0,1} \\
       \mathbf{O} & \mathbf{K}_{2} & \mathbf{B}^T_{0,2} \\
       \hline
       \mathbf{B}_{0,1} & \mathbf{B}_{0,2} & \mathbf{O}  
     \end{array}
   \right)
   \left(\begin{array}{l} \mathbf{R}_1~ \mathbf{O} \\
                           \mathbf{O}~~ \mathbf{R}_2 \\
                           \hline
                           \mathbf{O}~~ \mathbf{O} 
       \end{array} \right) 
   \left(\begin{array}{l} \mathbf{H}_1 \\ \mathbf{H}_2
       \end{array} \right)  = 
       \left(\begin{array}{l} \mathbf{O} \\ \mathbf{O} \\\hline \mathbf{O}
       \end{array} \right) 
     \label{eq:K1_clust}

or shortly

.. math::

   \left(
     \begin{array}{l|l}
       \mathbf{K}_{1:2} & \mathbf{B}^T_{0,1:2} \\
       \hline
       \mathbf{B}_{0,1:2} &  \mathbf{O}  
     \end{array}
   \right)
   \left(\begin{array}{c} \mathbf{R}_{1:2} \\
       \hline
                           \mathbf{O}
       \end{array} \right) 
       \mathbf{H}_{1:2} = 
       \left(\begin{array}{c} \mathbf{O} \\\hline \mathbf{O} 
   \end{array} \right).

..  \label{eq:K1_clust}

where the kernel is given by

.. math:: \tilde{\mathbf{R}}_1=(\mathbf{R}_{1:2}^T,~\mathbf{O})^T\mathbf{H}_{1:2}.

Assuming that the subdomian kernels :math:`\mathbf{R}_1` and
:math:`\mathbf{R}_2` are already known, the determination of
:math:`\mathbf{H}_{1:2}` remains. The first equation

.. math:: \mathbf{K}_{1:2}\mathbf{R}_{1:2}\mathbf{H}_{1:2}=\mathbf{O}

from does not impose any special conditions onto
:math:`\mathbf{H}_{1:2}`, because the product
:math:`\mathbf{K}_{1:2} \mathbf{R}_{1:2}` is already a zero matrix. The
second equation reads

.. math::

   \mathbf{B}_{0,1:2} \mathbf{R}_{1:2} \mathbf{H}_{1:2} = -\mathbf{G}_{0,1:2}^T \mathbf{H}_{1:2} = \mathbf{O}
     \label{eq:G_012_H12_0}

and it says, :math:`\mathbf{H}_{1:2}` is kernel of
:math:`\mathbf{G}_{0,1:2}^T`. It can be beneficial if kernels of
subdomains stored in :math:`\mathbf{R}_{1,2}` are obtained analytically
and from one coordinate system (see, e.g., eq. (3.7) in ). Then the
searched kernel is

.. math::

   \mathbf{H}_{1:2}=
   \left(\begin{array}{cc} \mathbf{I}_{3,3}& \mathbf{I}_{3,3}\end{array} \right)^T

where :math:`\mathbf{I}_{3,3}\in \mathbb{R}^{3\times3}` is an identity
matrix. Matrix :math:`\mathbf{H}_{1:2}` is also the kernel of
:math:`\mathbf{S}_{0,1:2}`, which can be regularized by it and then
easily factorized.


.. rubric:: Reference


.. 1 : ESPRESO

.. [1] ESPRESO - Exascale Parallel FETI Solver, http://espreso.it4i.cz.

.. 2 : Dostl2006

.. [2] Z. Dostal; D. Horak; and R. Kucera: Total FETI-an easier implementable variant of the FETI method for numerical solution of elliptic PDE.   Communications in Numerical Methods in Engineering, 22(12):1155--1162, jun 2006.

.. 3 : Dostl2012

.. [3] Z. Dostal; T. Kozubek; T. Brzobohaty; A. Markopoulos; and O. Vlach.: Scalable TFETI with optional preconditioning by conjugate projector: for transient frictionless contact problems of elasticity.   Computer Methods in Applied Mechanics and Engineering, 247-248:37--50, nov 2012.

.. 4 : Dostl2009

.. [4] Z. Dostal, T. Kozubek, V. Vondrak, T. Brzobohaty, and A. Markopoulos.: Scalable TFETI algorithm for the solution of multibody contact problems of elasticity.  International Journal for Numerical Methods in Engineering,
  pages n/a--n/a, 2009.

.. 5 : FarManRou-CMAME-94

.. [5] C. Farhat, J. Mandel, and F.-X. Roux.: Optimal convergence properties of the FETI domain decomposition method.   Computer Methods in Applied Mechanics and Engineering, 115:365--385, 1994.

.. 6 : Farhat1991

.. [6] C. Farhat and F.-X. Roux.: A method of finite element tearing and interconnecting and its
  parallel solution algorithm.   International Journal for Numerical Methods in Engineering, 32(6):1205--1227, oct 1991.

.. 7 : GosRey-ACME-06

.. [7] P. Gosselet and C. Rey.: Non-overlapping domain decomposition methods in structural mechanics.   Archives of Computational Methods in Engineering, 13(4):515--572, 2006.

.. 8 : Klawonn2010
.. [8] A. Klawonn and O. Rheinbach.: Highly scalable parallel domain decomposition methods with an application to biomechanics.  ZAMM, 90(1):5--32, jan 2010.

.. 9 : Klawonn2015

.. [9] A. Klawonn, M. Lanser, and O. Rheinbach: Toward extremely scalable nonlinear domain decomposition methods for elliptic partial differential equations.  SIAM J. Sci. Comput., 37(6):C667--C696, jan 2015.

.. 10 : TosWid2005

.. [10] A. Toselli and O. B. Widlund.: Domain Decomposition Methods {\textemdash} Algorithms and Theory.  Springer-Verlag, 2005.

.. 11 : OCZIENK

.. [11] O. C. Zienkiewicz and R. L. Taylor.: The finite element method. fifth edition.   Bautechnik, 79(2):122--123, feb 2002.



.. |The two approaches how to use HTFETI method and their comparison with TFETI for problems of size up to 3 billion DOF.| image:: FETI-motiv.pdf
