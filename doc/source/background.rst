

Background of the library
=========================
   
DDM & Krylov methods
====================

Finite Element Tearing and Interconnecting
------------------------------------------

Let us consider a system of linear equations in the form

.. math::

   \begin{aligned}
   \boldsymbol{A}\boldsymbol{x} & = & \boldsymbol{b}+\boldsymbol{B}^{\top}\boldsymbol{\lambda}\\
   \boldsymbol{B}\boldsymbol{x} & = & \boldsymbol{o}.\end{aligned}

where

-  :math:`\boldsymbol{A}` is block diagonal stiffness matrix (blocks can
   be singular),

-  :math:`\boldsymbol{b}` is right-hand side vector contains nodal
   forces,

-  :math:`\boldsymbol{\lambda}` are Lagrange multipliers,

-  :math:`\boldsymbol{B}` is equality constraints matrix.

The solution - nodal displacements - is

.. math:: \boldsymbol{x}=\boldsymbol{A}^{+}\left(\boldsymbol{b}+\boldsymbol{B}^\top\boldsymbol{\lambda}\right)+\boldsymbol{N}\boldsymbol{\alpha}.

The columns of the matrix :math:`\boldsymbol{N}` create a basis of the
null space of the matrix :math:`\boldsymbol{A}` (in elasticity so-called
rigid body modes), :math:`\boldsymbol{\alpha}` are its amplitudes. Using
previous terms and condition

.. math:: \boldsymbol{o}=\boldsymbol{N}^\top\left(\boldsymbol{b}+\boldsymbol{B}^\top\boldsymbol{\lambda}\right)

original problem can be rewritten following way

.. math::

   \begin{aligned}
   \boldsymbol{B}\boldsymbol{A}^{+}\boldsymbol{B}^\top\boldsymbol{\lambda}+\boldsymbol{B}\boldsymbol{N}\boldsymbol{\alpha} & = & -\boldsymbol{B}\boldsymbol{A}^{+}\boldsymbol{b}\\
   \boldsymbol{N}^\top\boldsymbol{B}^\top\boldsymbol{\lambda} & = & -\boldsymbol{N}^\top\boldsymbol{b}.\end{aligned}

The same equations in the standard FETI notation are

.. math::

   \begin{aligned}
   \boldsymbol{F}\boldsymbol{\lambda}+\boldsymbol{G}\boldsymbol{\alpha} & = & \boldsymbol{d}_{0}\\
   \boldsymbol{G}^\top\boldsymbol{\lambda} & = & \boldsymbol{e}.\end{aligned}

Orthogonal projection
~~~~~~~~~~~~~~~~~~~~~

Gradient in :math:`i-`\ th iteration can be expressed as

.. math:: \boldsymbol{g}_{i}=\boldsymbol{F}\boldsymbol{\lambda}_{i}-\boldsymbol{d}_{0}+\boldsymbol{G}\boldsymbol{\alpha}_{i},

and iterative process can be accelerated using Petrov-Galerking
condition

.. math:: \boldsymbol{G}^\top\boldsymbol{g}_{i,opt}=\boldsymbol{o},\

the optimal coefficients :math:`\boldsymbol{\alpha}_{i}` are given by

.. math:: \boldsymbol{\alpha}_{i,opt}=-\left(\boldsymbol{G}^\top\boldsymbol{G}\right)^{-1}\boldsymbol{G}^\top\left(\boldsymbol{F}\boldsymbol{\lambda}_{i}-\boldsymbol{d}_{0}\right).

Optimized gradient is

.. math:: \boldsymbol{g}_{i,opt}={\boldsymbol{F}\boldsymbol{\lambda}_{i}-\boldsymbol{d}_{0}}-\boldsymbol{G}\left(\boldsymbol{G}^\top\boldsymbol{G}\right)^{-1}\boldsymbol{G}^\top\left({\boldsymbol{F}\boldsymbol{\lambda}_{i}-\boldsymbol{d}_{0}}\right),

so to having an optimal gradient it is just the same, if the projection
is applied onto
:math:`\boldsymbol{g}_{\lambda,i}={\boldsymbol{F}\boldsymbol{\lambda}_{i}-\boldsymbol{d}_{0}}`
by projector

.. math:: \boldsymbol{P}=\boldsymbol{I}-\boldsymbol{Q},\ \mbox{where}\boldsymbol{\ Q}=\boldsymbol{G}\left(\boldsymbol{G}^\top\boldsymbol{G}\right)^{-1}\boldsymbol{G}^\top.

.. math:: \boldsymbol{g}_{i,opt}=\bar{\boldsymbol{g}}_{i}=\boldsymbol{P}\boldsymbol{g}_{\lambda,i}.

Original problem is now in the form

.. math:: \boldsymbol{P}\boldsymbol{F}\boldsymbol{\lambda}=\boldsymbol{P}\boldsymbol{d}_{0}.

The solution :math:`\boldsymbol{\lambda}` can be split in two parts
with the same projectors used before

.. math:: \boldsymbol{\lambda}=\boldsymbol{P}\boldsymbol{\lambda}+\boldsymbol{Q}\boldsymbol{\lambda}=\bar{\boldsymbol{\lambda}}+\hat{\boldsymbol{\lambda}},

:math:`\bar{\boldsymbol{\lambda}}\in\mbox{Ker}\boldsymbol{G}^\top` will
be obtained from iterative process (bar over
:math:`\boldsymbol{\lambda}` points out, vector is projected), and
:math:`\hat{\boldsymbol{\lambda}}\in\mbox{Im}\boldsymbol{G}` is obtained
directly from

.. math:: \hat{\boldsymbol{\lambda}}=\boldsymbol{Q}\boldsymbol{\lambda}=\boldsymbol{G}\left(\boldsymbol{G}^\top\boldsymbol{G}\right)^{-1}{\boldsymbol{G}^\top\boldsymbol{\lambda}}=\boldsymbol{G}\left(\boldsymbol{G}^\top\boldsymbol{G}\right)^{-1}{\boldsymbol{e}}.

Finally the problem becomes

.. math:: \boldsymbol{P}\boldsymbol{F}\bar{\boldsymbol{\lambda}}=\boldsymbol{P}\left(\boldsymbol{d}_{0}-\boldsymbol{F}\hat{\boldsymbol{\lambda}}\right),

.. math:: \boldsymbol{P}\boldsymbol{F}\bar{\boldsymbol{\lambda}}=\boldsymbol{P}\boldsymbol{d},

or using bar over matrix and vectors (pointing out to the projection)
is expressed as

.. math:: \bar{\boldsymbol{F}}\bar{\boldsymbol{\lambda}}=\bar{\boldsymbol{d}}.

Linear system
-------------

Linear system of equations is written in the form

.. math::

   \boldsymbol{Ax}=\boldsymbol{b}.

Using incomplete Choleski factorization

.. math::

   \boldsymbol{A}^{-1 } \approx{\boldsymbol{M}} = \boldsymbol{L}^{-T} \boldsymbol{L}^{-1}

.. math::

   \boldsymbol{PAP}^{-1 } \approx{\boldsymbol{M}} = \boldsymbol{L}^{-T} \boldsymbol{L}^{-1}

can be preconditioned and original problem is rewritten into

.. math:: \boldsymbol{L}^{-1}\boldsymbol{A}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\boldsymbol{x}=\boldsymbol{L}^{-1}\boldsymbol{b}\\

or with new notation

.. math:: \hat{\boldsymbol{A}}\hat{\boldsymbol{x}}=\hat{\boldsymbol{b}}.

Conjugate gradient methods
--------------------------

Essential algorithm of conjugate gradient method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

System can be solved using following algorithm:

| [H] initialization; :math:`i=0,~\boldsymbol{x}_{0},
      ~ \boldsymbol{g}_{0}=\boldsymbol{Ax}_{0} - \boldsymbol{b},
      ~ \boldsymbol{w}_{0}=\boldsymbol{g}_{0}, ~n > 0,~\varepsilon>0`
|  [alg:cg]

Preconditioned conjugate gradient method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Transformation of the system preconditioned
by factor :math:`\boldsymbol{L}` to the system with left preconditioner
is derived in following table.

.. math::
   \begin{array}{|c|c|c|}
   \hline
   \boldsymbol{1.} & \boldsymbol{2.} & \boldsymbol{3.} \\
   \hat{\boldsymbol{g}}_{i}=\hat{\boldsymbol{A}}\hat{\boldsymbol{x}}_i-\hat{\boldsymbol{b}} &
   \boldsymbol{L}^{-1}\boldsymbol{g}_{i}=\boldsymbol{L}^{-1} \left(\boldsymbol{A}\boldsymbol{x}_i-\boldsymbol{b}\right) &
   \boldsymbol{g}_{i}=\boldsymbol{A}\boldsymbol{x}_i-\boldsymbol{b}\\
   
   \hat{\boldsymbol{x}}_{i+1}=\hat{\boldsymbol{x}}_{i} + \hat{\boldsymbol{w}}_i\rho_i   &  
   \boldsymbol{L}^\top{\boldsymbol{x}}_{i+1}=\boldsymbol{L}^\top \left({\boldsymbol{x}}_{i} + {\boldsymbol{w}}_i\rho_i\right) &
   \boldsymbol{x}_{i+1}=\boldsymbol{x}_{i} + \boldsymbol{w}_i\rho_i \\
   \hat{\boldsymbol{g}}_{i+1}=\hat{\boldsymbol{g}}_{i} + \hat{\boldsymbol{A}}\hat{\boldsymbol{w}}_i\rho_i&  
   \boldsymbol{L}^{-1}\boldsymbol{g}_{i+1}=\boldsymbol{L}^{-1}\left(\boldsymbol{g}_{i} + \boldsymbol{Aw}_i\rho \right) &
   \boldsymbol{g}_{i+1}=\boldsymbol{g}_{i} + \boldsymbol{Aw}_i\rho\\
   \hat{\boldsymbol{w}}_{i+1}=\hat{\boldsymbol{g}}_{i+1} + \hat{\boldsymbol{w}}_i\gamma_i&  
   \boldsymbol{L}^\top\boldsymbol{w}_{i+1}=\boldsymbol{L}^{-1}\boldsymbol{g}_{i+1} + \boldsymbol{L}^\top\boldsymbol{w}_i\gamma_i &
   \boldsymbol{w}_{i+1}=\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}\boldsymbol{g}_{i+1} + \boldsymbol{w}_i\gamma_i\\
   
   \rho_{i}= -
   {{\left(\hat{\boldsymbol{g}}_{i},\hat{\boldsymbol{g}}_{i}\right)}
     \over 
   {\left(\hat{\boldsymbol{w}}_{i},\hat{\boldsymbol{A}}\hat{\boldsymbol{w}}_{i}\right)}}&  
   \rho_{i}= -
   {{\left(\boldsymbol{L}^{-1}\boldsymbol{g}_{i},\boldsymbol{L}^{-1}\boldsymbol{g}_{i}\right)}
     \over 
   {\left(\boldsymbol{L}^\top\boldsymbol{w}_{i},\boldsymbol{L}^{-1}\boldsymbol{A}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\boldsymbol{w}_{i}\right)}} &
     
   \rho_{i}= -
   {{\left(\boldsymbol{g}_{i},\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}\boldsymbol{g}_{i}\right)}
     \over 
   {\left(\boldsymbol{w}_{i},\boldsymbol{A}\boldsymbol{w}_{i}\right)}} \\
   \gamma_{i}= 
   {{\left(\hat{\boldsymbol{g}}_{i+1},\hat{\boldsymbol{g}}_{i+1}\right)}
     \over 
   {\left(\hat{\boldsymbol{g}}_{i},\hat{\boldsymbol{g}}_{i}\right)}}&  
     \gamma_{i}= 
   {{\left(\boldsymbol{L}^{-1}\boldsymbol{g}_{i+1},\boldsymbol{L}^{-1}\boldsymbol{g}_{i+1}\right)}
     \over 
   {\left(\boldsymbol{L}^{-1}\boldsymbol{g}_{i},\boldsymbol{L}^{-1}\boldsymbol{g}_{i}\right)}} & 
     \gamma_{i}= 
   {{\left(\boldsymbol{g}_{i+1},\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}\boldsymbol{g}_{i+1}\right)}
     \over 
   {\left(\boldsymbol{g}_{i},\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}\boldsymbol{g}_{i}\right)}}\\
   \hline
   \end{array}

First column in table above contains standard terms derived in previous
section, only the differ is, solved system is marked by symbol 
:math:`\hat{}`. Second column is formally same as the previous one, but
original matrix :math:`\boldsymbol{A}` and right-hand-side vector
:math:`\boldsymbol{b}` is used. Finally the last column is modification
of second column so that an approximation of solution
:math:`\boldsymbol{x}_i` is generated directly (not
:math:`\hat{\boldsymbol{x}}_i`) and moreover, instead of factor
:math:`\boldsymbol{L}` undecomposed preconditioned
:math:`\boldsymbol{M}}` is applied
(:math:`\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}=\boldsymbol{M}`).

| [H] initialization; :math:`i=0,~\boldsymbol{x}_{0},
      ~ \boldsymbol{g}_{0}=\boldsymbol{Ax}_{0} - \boldsymbol{b},
      ~ \boldsymbol{w}_{0}={\boldsymbol{M}}\boldsymbol{g}_{0}, ~n > 0,~\varepsilon>0`

Projected conjugate gradient method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

System to be solved

.. math::

   \begin{array}{l}
     \boldsymbol{Ax}=\boldsymbol{b} + \boldsymbol{B}^\top\boldsymbol{\lambda}\\ 
     \boldsymbol{Bx}=\boldsymbol{c}
     \end{array}

differs from the previous case by additional part including Lagrange
multipliers. An approximation of the solution :math:`\boldsymbol{x}` can
be written in the form

.. math:: \boldsymbol{x}_i = \boldsymbol{x}_0 + \sum \limits_{j=0}^{i-1}\boldsymbol{w}_j \rho_j.

If initial guess satisfies

.. math:: \boldsymbol{B}\boldsymbol{x}_0 = \boldsymbol{c},

due to second equality in ([eq:Ax:sub:`bB`\ x\ :sub:`c`]) each
:math:`\boldsymbol{w}_j` has to be orthogonal to basis
:math:`\boldsymbol{B}^\top`

.. math::

   \boldsymbol{B}\boldsymbol{w}_j = \boldsymbol{o},~\forall j,~ j \in (0,1,...,i-1)

Due to Lagrange multipliers current gradient is

.. math::

   \boldsymbol{g}_{i}=\boldsymbol{g}_{x,i} - \boldsymbol{B}^\top\boldsymbol{\lambda}_i 
               ~\mbox{ where } ~
               \boldsymbol{g}_{x,i} = \boldsymbol{Ax}_{i} - \boldsymbol{b}.

and next approximation of the solution is

.. math:: \boldsymbol{x}_{i+1}=\boldsymbol{x}_{i} + \boldsymbol{w}_i\rho_i.

 Next gradient is

.. math:: \boldsymbol{g}_{i+1}=\boldsymbol{g}_{i+1,x} - \boldsymbol{B}^\top\boldsymbol{\lambda}_{i+1}

 and finally conjugate direction in iteration :math:`i+1` is

.. math:: \boldsymbol{w}_{i+1}=\boldsymbol{g}_{i+1} + \boldsymbol{w}_i\gamma_i.

To describe treatment with :math:`\boldsymbol{w}_j` (orthogonal to
:math:`\boldsymbol{B}^\top`) properly we start with initial conjugate
direction which is also initial gradient

.. math:: \boldsymbol{w}_0=\boldsymbol{g}_0=\boldsymbol{g}_{x,0}-\boldsymbol{B}^\top\boldsymbol{\lambda_0}.

Orthogonality condition ([eq:Bw:sub:`0`]) is enforced by choosing an
optimal LM

.. math:: \boldsymbol{\lambda}_0 =\left(\boldsymbol{B}\boldsymbol{B}^\top \right)^{-1}\boldsymbol{B}\boldsymbol{g}_{0,x}

.. math:: \boldsymbol{\lambda}_i =\left(\boldsymbol{B}\boldsymbol{B}^\top \right)^{-1}\boldsymbol{B}\boldsymbol{g}_{i,x}

obtained from

.. math::

   \boldsymbol{B}\boldsymbol{g}_{o}=
     \boldsymbol{B}\boldsymbol{g}_{o,x}-\boldsymbol{B}\boldsymbol{B}^\top\boldsymbol{\lambda}_0
     = \boldsymbol{o}.

This optimized gradient can be written in the form

.. math::

   \boldsymbol{g}_{0} = 
       \boldsymbol{g}_{x,0}- 
         \boldsymbol{B}^\top 
           \left(  
               \boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{B} 
         \boldsymbol{g}_{x,0}

or introducing projector

.. math:: \boldsymbol{P}=\boldsymbol{I}-\boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{B}

briefly

.. math:: \boldsymbol{g}_{0}=\boldsymbol{P}\boldsymbol{g}_{x,0}.

Orthogonality condition ([eq:Bw:sub:`0`]) in next iteration is applied
into

.. math::

   \boldsymbol{w}_1=\boldsymbol{g}_{x,1}-\boldsymbol{B}^\top\boldsymbol{\lambda_1}+
     \boldsymbol{w}_0\gamma_0.

where :math:`\boldsymbol{w}_0` is already orthogonal to
:math:`\boldsymbol{B}^\top` from previous step. Applying former approach
it can be written

.. math:: \boldsymbol{w}_{1}=\boldsymbol{P}\boldsymbol{g}_{x,1} + \boldsymbol{w}_0\gamma_0.

or generally

.. math:: \boldsymbol{w}_{i+1}=\boldsymbol{P}\boldsymbol{g}_{x,i+1} + \boldsymbol{w}_i\gamma_i.

To get constants :math:`\rho_i` and :math:`\gamma_i` projected gradient
is used

.. math::

   \rho_{i}= -{{\left(\boldsymbol{Pg}_{x,i},\boldsymbol{Pg}_{x,i}\right)}
         \over 
       {\left(\boldsymbol{w}_{i},\boldsymbol{Aw}_{i}\right)}},~
     \gamma_{i}= 
     {{\left(\boldsymbol{Pg}_{x,i+1},\boldsymbol{Pg}_{x,i+1}\right)}
         \over 
         {\left(\boldsymbol{Pg}_{x,i},\boldsymbol{Pg}_{x,i}\right)}}

and because of the fundamental property of projector

.. math:: \boldsymbol{P}^n =  \boldsymbol{P},~n=1,2,3,...

they can be written

.. math::

   \rho_{i}= -{{\left(\boldsymbol{g}_{x,i},\boldsymbol{Pg}_{x,i}\right)}
         \over 
       {\left(\boldsymbol{w}_{i},\boldsymbol{Aw}_{i}\right)}},~
     \gamma_{i}= 
     {{\left(\boldsymbol{g}_{x,i+1},\boldsymbol{Pg}_{x,i+1}\right)}
         \over 
         {\left(\boldsymbol{g}_{x,i},\boldsymbol{Pg}_{x,i}\right)}}.

Complete algorithm is written bellow

| [H] initialization;
:math:`i=0,~\boldsymbol{x}_{0}=\boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{c},
~\boldsymbol{g}_{x,0}=\boldsymbol{Ax}_{0} - \boldsymbol{b},
~\boldsymbol{w}_{0}=\boldsymbol{P}\boldsymbol{g}_{x,0}, ~n > 0,~\varepsilon>0`
|  [alg:PCG:sub:`v`\ ar\ :sub:`1`]

Projected preconditioned conjugate gradient method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

System to be solved is

.. math::

   \begin{array}{l}
       \hat{\boldsymbol{A}}\hat{\boldsymbol{x}}=\hat{\boldsymbol{b}} + \hat{\boldsymbol{B}}^\top{\boldsymbol{\lambda}}\\ 
       \hat{\boldsymbol{B}}\hat{\boldsymbol{x}}=\hat{\boldsymbol{c}} 
     \end{array}

which is preconditioned by factor :math:`\boldsymbol{L}`. Written with
original objects and factor of preconditioner

.. math::

   \begin{array}{l}
       \boldsymbol{L}^{-1}\boldsymbol{A}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\boldsymbol{x}=
       \boldsymbol{L}^{-1}\left(\boldsymbol{b} + \boldsymbol{B}^\top{\boldsymbol{\lambda}}\right)\\ 
       \boldsymbol{B}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\boldsymbol{x}=\boldsymbol{c} 
     \end{array}.

Similarly as before we will try use precondtitioner
:math:`{\boldsymbol{M}}` instead of factor
:math:`\boldsymbol{L}`. Orthogonality condition ([eq:Bw:sub:`0`]) is
influenced by modified matrix :math:`\boldsymbol{B}`

.. math::

   \hat{\boldsymbol{B}}\hat{\boldsymbol{w}}_j = 
     \boldsymbol{B}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\boldsymbol{w}_j = \boldsymbol{o}

and very first conjugate direction is

.. math::

   \begin{array}{l}
       \hat{\boldsymbol{w}}_0=\hat{\boldsymbol{g}}_{0}=\hat{\boldsymbol{g}}_{x,0}-\hat{\boldsymbol{B}}^\top{\boldsymbol{\lambda_0}}
           \\
           \boldsymbol{L}^\top\boldsymbol{w}_0=\boldsymbol{L}^{-1}\boldsymbol{g}_0=\boldsymbol{L}^{-1}\left(\boldsymbol{g}_{x,0}
           -\boldsymbol{B}^\top{\boldsymbol{\lambda}}_0\right)
     \end{array}.

We need find such :math:`\boldsymbol\lambda_0` for which gradient
:math:`\hat{\boldsymbol{g}}_0` will be orthogonal to
:math:`\hat{\boldsymbol{B}}^\top`. Therefore

.. math::

   \hat{\boldsymbol{B}}\hat{\boldsymbol{w}}_0= \boldsymbol{B}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\boldsymbol{w}_0=
       \boldsymbol{B}\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}\left(\boldsymbol{g}_{x,0}
       -\boldsymbol{B}^\top{\boldsymbol{\lambda}}_0\right)=\boldsymbol{o}.

Because
:math:`{\boldsymbol{M}}=\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}`,
optimized vector :math:`\boldsymbol\lambda_0` is

.. math::

   {\boldsymbol{\lambda}}_0
         =   
         \left(\boldsymbol{B}{\boldsymbol{M}}\boldsymbol{B}^\top \right)^{-1}
         \boldsymbol{B}^\top{\boldsymbol{M}}\boldsymbol{g}_{x,0}.

Projected gradient is

.. math:: \hat{\boldsymbol{g}}_0=\boldsymbol{L}^{-1}\boldsymbol{P}_{\!o}~\boldsymbol{g}_{x,0}

where

.. math::

   \boldsymbol{P}_{\!o} = 
     \boldsymbol{I}-\boldsymbol{B}^\top\left(\boldsymbol{B}{\boldsymbol{M}}\boldsymbol{B}^\top \right)^{-1}
     \boldsymbol{B}{\boldsymbol{M}}

is oblique projector augmented by preconditioner
:math:`{\boldsymbol{M}}`. The norm of initial gradient is

.. math::

   \sqrt{\left(\hat{\boldsymbol{g}}_0,\hat{\boldsymbol{g}}_0\right)}=
     \sqrt{\left(\boldsymbol{g}_{x,0},\boldsymbol{P}_{\!o}^\top{\boldsymbol{M}}\boldsymbol{P}_{\!o}~\boldsymbol{g}_{x,0}\right)}.

Any other conjugate vector is

.. math:: \hat{\boldsymbol{w}}_{i+1}=\hat{\boldsymbol{g}}_{i+1}+\hat{\boldsymbol{w}}_{i}\gamma_i

which still has to be orthogonal to :math:`\boldsymbol{B}^\top`. Second
part with :math:`\boldsymbol{w}_i` already fulfils this conditions from
previous iteration thus we need specify ’only’ optimal
:math:`\hat{\boldsymbol{g}}_{i+1}`. But this step is equivalent to
procedure applied on initial gradient described before, therefore
:math:`(i+1)`-th optimized gradient is

.. math:: \hat{\boldsymbol{g}}_{i+1}=\boldsymbol{L}^{-1}\boldsymbol{P}_{\!o}~\boldsymbol{g}_{x,i+1}.

Finally the transformation from preconditioned system by factor
:math:`\boldsymbol{L}` to the new one preconditioned by unfactorized
:math:`{\boldsymbol{M}}` is shown on gradient and conjugate
direction in the table below. Vectors in first column are generated
implicitly. Via step in column :math:`2` in last column there are
vectors which are performed explicitly.

.. math::

   \begin{array}{|c|c|c|}
     \hline
       {\boldsymbol{1.}} & {\boldsymbol{2.}} & {\boldsymbol{3.}} \\
     \hat{\boldsymbol{w}}_{i+1}=\hat{\boldsymbol{g}}_{i+1}+\hat{\boldsymbol{w}}_{i}\gamma_i & 
           \boldsymbol{L}^{T^{}}\boldsymbol{w}_{i+1}=
           \boldsymbol{L}^{-1}\boldsymbol{P}\!_o\boldsymbol{g}_{i+1} + \boldsymbol{L}^\top \boldsymbol{w}_i\gamma_i
                &
           \boldsymbol{w}_{i+1}=
           {\boldsymbol{M}}\boldsymbol{P}_{\!o}\boldsymbol{g}_{i+1} + \boldsymbol{w}_i\gamma_i\\
                \hat{\boldsymbol{g}}_{i+1}=\boldsymbol{L}^{-1}\boldsymbol{P}_{\!o}~\boldsymbol{g}_{x,i+1}& 
                \boldsymbol{L}^{-1}\boldsymbol{g}_{i+1}=\boldsymbol{L}^{-1}\boldsymbol{P}_{\!o}~\boldsymbol{g}_{x,i+1}
                &
           \boldsymbol{g}_{i+1}=\boldsymbol{P}_{\!o}~\boldsymbol{g}_{x,i+1}. \\
     \hline
     \end{array}

Itself algorithm follows. Stopping criteria is measured by projected
gradient (without preconditioner :math:`{\boldsymbol{M}}`)
and thus two projectors are introduced.

| [H] initialization;
:math:`i=0,~\boldsymbol{x}_{0}=\boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{c},
~ \boldsymbol{g}_{x,0}=\boldsymbol{Ax}_{0} - \boldsymbol{b},
~ \boldsymbol{w}_{0}=\boldsymbol{P}_{\!o}\boldsymbol{g}_{x,0}, ~n > 0,~\varepsilon>0`

Conjugate gradient methods - orthogonal and oblique projectors
--------------------------------------------------------------

Projected conjugate gradient method - var. II.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Projection method to get the solution of KKT system can be introduced
different way compared to section [sec:ProjConGraMet].

.. math::

   \begin{array}{l}
     \boldsymbol{Ax}=\boldsymbol{b} + \boldsymbol{B}^\top\boldsymbol{\lambda}\\ 
     \boldsymbol{Bx}=\boldsymbol{c}.
     \end{array}

Solution :math:`\boldsymbol{x}` can be split into two parts

.. math:: \boldsymbol{x} = \boldsymbol{x}_{{Im B}} + \boldsymbol{x}_{{Ker B}^\top}.

which using orthogonal projectors

.. math::

   \begin{array}{l}
       \boldsymbol{Q} = \boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{B} \\
       \boldsymbol{P} = \boldsymbol{I}-\boldsymbol{Q}
     \end{array}

are

.. math::

   \begin{array}{l}
         \boldsymbol{x}_{{Ker B}^\top} =   \boldsymbol{Px} \\
         \boldsymbol{x}_{{Im B}} = \boldsymbol{x}_{0} = \boldsymbol{Qx} =
         \boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{Bx} =
         \boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{c}.
     \end{array}

Now it is very natural to use projector :math:`\boldsymbol{P}`,
multiply first equation and get modified system in the form

.. math::

   \begin{array}{l}
     \boldsymbol{PAP}  \boldsymbol{x}_{{Ker B}^\top} =
            \boldsymbol{P}\left(\boldsymbol{b}-\boldsymbol{Ax}_{0}\right) \\
     \boldsymbol{Bx}_{{Ker B}^\top}=\boldsymbol{o},
     \end{array}

where moreover the condition

.. math:: \boldsymbol{x}_{{KerB}^\top} =  \boldsymbol{Px}_{{KerB}^\top}

was used. Finally the system to be solved using substitution
:math:`\bar{\boldsymbol{x}} = \boldsymbol{x}_{{KerB}^\top}` is

.. math:: \boldsymbol{PAP}  \bar{\boldsymbol{x}} = \boldsymbol{P}\left(\boldsymbol{b}-\boldsymbol{Ax}_0\right).

Formally we can use algorithm [alg:cg] with replacing matrix of the
linear system and rhs by their projected versions, and choosing
appropriate initial guess :math:`\boldsymbol{x}_0`.

| [H] initialization;
:math:`i=0,~\boldsymbol{x}_{0}=\boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{c},~
 \bar{\boldsymbol{x}}_{0}=\boldsymbol{0}`
:math:`\boldsymbol{g}_{x,0}=\boldsymbol{Ax}_{0} - \boldsymbol{b},~
      ~ \boldsymbol{w}_{0}=\boldsymbol{Pg}_{x,0}, ~n > 0,~\varepsilon>0`

Algorithm contains many appearances of projector :math:`\boldsymbol{P}`.
Thanks to the property
:math:`\boldsymbol{P}=\boldsymbol{P}^k=\boldsymbol{P}\boldsymbol{P}\cdots \boldsymbol{P}`
we can write simplified product

.. math::

   {\left(\boldsymbol{Pg}_{x,i},\boldsymbol{Pg}_{x,i}\right)}=
   {\left(\boldsymbol{g}_{x,i},\boldsymbol{Pg}_{x,i}\right)}.

Another product

.. math::

   \left(\boldsymbol{w}_{i},\boldsymbol{PAPw}_{i}\right)=
   \left(\boldsymbol{w}_{i},\boldsymbol{Aw}_{i}\right)

is simplified because conjugate direction

.. math::

   \boldsymbol{w}_i = \mbox{span} 
     \left(
       \boldsymbol{Pg}_{x,0},
       ~\boldsymbol{Pg}_{x,1},
       \cdots,
       ~\boldsymbol{Pg}_{x,i}
     \right)

is actually given by linear combination of projected gradients,
therefore it is not necessary redundantly perform product
:math:`\boldsymbol{PAP}` in noticed dot product. Accepting those facts
we can modify the algorithm following way:

| [H] initialization;
:math:`i=0,~\boldsymbol{x}_{0}=\boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{c},~
 \bar{\boldsymbol{x}}_{0}=\boldsymbol{0}`
:math:`\boldsymbol{g}_{x,0}=\boldsymbol{Ax}_{0} - \boldsymbol{b},
      ~ \boldsymbol{w}_{0}=\boldsymbol{Pg}_{x,0}, ~n > 0,~\varepsilon>0`
| 
:math:`\boldsymbol{x} \approx \boldsymbol{x}_0 + \bar{\boldsymbol{x}}_i`

Scheme of this algorithm is identical to version of algorithm
[alg:PCG:sub:`v`\ ar\ :sub:`1`].

Projected preconditioned conjugate gradient method - var. II.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Applying of preconditioning on the projected problem

.. math:: \boldsymbol{PAP}  \bar{\boldsymbol{x}} = \boldsymbol{P}\left(\boldsymbol{b}-\boldsymbol{Ax}_0\right).

leads to

.. math:: \boldsymbol{L}^{-1}\boldsymbol{PAP}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\bar{\boldsymbol{x}} = \boldsymbol{L}^{-1}\boldsymbol{P}\left(\boldsymbol{b}-\boldsymbol{Ax}_0\right)

and shortly

.. math:: \hat{\boldsymbol{A}}\hat{\boldsymbol{x}}=\hat{\boldsymbol{b}}.

For such problem the fundamental terms are again those below

.. math::

   \begin{array}{l}
   \hat{\boldsymbol{g}}_{x,i}=\hat{\boldsymbol{A}}\hat{\boldsymbol{x}}_i-\hat{\boldsymbol{b}} \\
   \hat{\boldsymbol{x}}_{i+1}=\hat{\boldsymbol{x}}_{i} + \hat{\boldsymbol{w}}_i\rho_i  \\ 
   \hat{\boldsymbol{g}}_{x,i+1}=\hat{\boldsymbol{g}}_{x,i} + \hat{\boldsymbol{A}}\hat{\boldsymbol{w}}_i\rho_i\\
   \hat{\boldsymbol{w}}_{i+1}=\hat{\boldsymbol{g}}_{i+1} + \hat{\boldsymbol{w}}_i\gamma_i\\
       \rho_{i}= -
       {{\left(\hat{\boldsymbol{g}}_{x,i},\hat{\boldsymbol{g}}_{x,i}\right)}
           \over 
         {\left(\hat{\boldsymbol{w}}_{i},\hat{\boldsymbol{A}}\hat{\boldsymbol{w}}_{i}\right)}}\\
       \gamma_{i}= 
       {{\left(\hat{\boldsymbol{g}}_{x,i+1},\hat{\boldsymbol{g}}_{x,i+1}\right)}
           \over 
         {\left(\hat{\boldsymbol{g}}_{x,i},\hat{\boldsymbol{g}}_{i}\right)}}
       \end{array}

and the transformation from the variant of the algorithm with factor
:math:`\boldsymbol{L}` to the case with preconditioner
:math:`{\boldsymbol{M}}` is in the next table.

.. math::

   \begin{array}{|c|c|}
       \hline
       \boldsymbol{2.} & \boldsymbol{3.} \\
       \boldsymbol{L}^{-1}\boldsymbol{P}\boldsymbol{g}_{x,i}=\boldsymbol{L}^{-1}\left(\boldsymbol{P}\boldsymbol{A}\boldsymbol{P}\bar{\boldsymbol{x}}_i+\boldsymbol{P}\boldsymbol{A}\boldsymbol{x}_{0}-\boldsymbol{P}\boldsymbol{b}\right) &
       \boldsymbol{Pg}_{x,i}=\boldsymbol{PAP}\bar{\boldsymbol{x}}_i+\boldsymbol{PA}\boldsymbol{x}_0-\boldsymbol{Pb}\\
       \boldsymbol{L}^\top{\bar{\boldsymbol{x}}}_{i+1}=\boldsymbol{L}^\top \left({\bar{\boldsymbol{x}}}_{i} + {\boldsymbol{w}}_i\rho_i\right) &
       \bar{\boldsymbol{x}}_{i+1}=\bar{\boldsymbol{x}}_{i} + \boldsymbol{w}_i\rho_i \\
       \boldsymbol{L}^{-1}\boldsymbol{P}\boldsymbol{g}_{x,i+1}=\boldsymbol{L}^{-1}\left(\boldsymbol{P}\boldsymbol{g}_{x,i} + \boldsymbol{P}\boldsymbol{APw}_i\rho \right) &
       \boldsymbol{Pg}_{x,i+1}=\boldsymbol{Pg}_{x,i} + \boldsymbol{PAPw}_i\rho\\
       \boldsymbol{L}^\top\boldsymbol{w}_{i+1}=\boldsymbol{L}^{-1}\boldsymbol{g}_{i+1} + \boldsymbol{L}^\top\boldsymbol{w}_i\gamma_i &
       \boldsymbol{w}_{i+1}={\boldsymbol{M}}\boldsymbol{Pg}_{x,i+1} + \boldsymbol{w}_i\gamma_i\\
       \rho_{i}= -
       {{\left(\boldsymbol{L}^{-1}\boldsymbol{Pg}_{x,i},\boldsymbol{L}^{-1}\boldsymbol{Pg}_{i}\right)}
           \over 
         {\left(\boldsymbol{L}^\top\boldsymbol{w}_{i},\boldsymbol{L}^{-1}\boldsymbol{PAP}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\boldsymbol{w}_{i}\right)}} &
       \rho_{i}= -
       {{\left(\boldsymbol{Pg}_{i},{\boldsymbol{M}}\boldsymbol{Pg}_{x,i}\right)}
           \over 
       {\left(\boldsymbol{w}_{i},\boldsymbol{PAP}\boldsymbol{w}_{i}\right)}} \\
           \gamma_{i}= 
       {{\left(\boldsymbol{L}^{-1}\boldsymbol{Pg}_{i+1},\boldsymbol{L}^{-1}\boldsymbol{Pg}_{i+1}\right)}
           \over 
         {\left(\boldsymbol{L}^{-1}\boldsymbol{Pg}_{i},\boldsymbol{L}^{-1}\boldsymbol{Pg}_{i}\right)}} & 
           \gamma_{i}= 
           {{\left(\boldsymbol{Pg}_{x,i+1},{\boldsymbol{M}}\boldsymbol{Pg}_{x,i+1}\right)}
           \over 
         {\left(\boldsymbol{Pg}_{x,i},{\boldsymbol{M}}\boldsymbol{Pg}_{x,i}\right)}}\\
       \hline
     \end{array}

Algorithm in second column can be simplified more. Important is the
term for conjugate direction

.. math:: \boldsymbol{w}_{i+1}={\boldsymbol{M}}\boldsymbol{Pg}_{x,i+1} + \boldsymbol{w}_i\gamma_i.

which has to be orthogonal to basis :math:`\boldsymbol{B}^\top`. By the
former knowledge the part :math:`\boldsymbol{w}_{i}` already fulfils
orthogonality condition, so there is remaining

.. math:: {\boldsymbol{M}}\boldsymbol{Pg}_{x,i+1}.

Another projection applied onto :math:`\boldsymbol{w}` is an easy way
how to satisfy
:math:`\boldsymbol{BP}\boldsymbol{w}_{x,i+1}=\boldsymbol{o}`.

| [H] initialization;
:math:`i=0,~\boldsymbol{x}_{0}=\boldsymbol{B}^\top\left(\boldsymbol{B}\boldsymbol{B}^\top\right)^{-1}\boldsymbol{c},~
 \bar{\boldsymbol{x}}_{0}=\boldsymbol{0}`
| :math:`{\boldsymbol{M}}\approx\left(\boldsymbol{PAP}\right)^{-1}`
:math:`\boldsymbol{g}_{x,0}=\boldsymbol{Ax}_{0} - \boldsymbol{b},
      ~ \boldsymbol{w}_{0}=\boldsymbol{PMPg}_{x,0}, ~n > 0,~\varepsilon>0`
| 
:math:`\boldsymbol{x} \approx \boldsymbol{x}_0 + \bar{\boldsymbol{x}}_i`

Conjugate gradient methods based on Lanczos method
--------------------------------------------------

ordinary
~~~~~~~~

System of linear equations ([eq:Ax:sub:`b`]) will be solved using an
approximation

.. math::

   \boldsymbol{x}_{i+1} = \boldsymbol{x}_0 + \sum\limits_{j=0}^{i}\boldsymbol{v}_j y_{j} = 
    \boldsymbol{x}_0 + \boldsymbol{V}_i\boldsymbol{y}_i.

Related gradient in :math:`i+1` iteration is

.. math::

   \boldsymbol{g}_{i+1} = 
    \boldsymbol{g}_0 + \boldsymbol{AV}_i\boldsymbol{y}_i

which is orthogonal to current basis :math:`\boldsymbol{V}_i`

.. math::

   \boldsymbol{V}_i^\top\boldsymbol{g}_{i+1} =\boldsymbol{o}=
    \boldsymbol{V}_i^\top\boldsymbol{g}_0 +\boldsymbol{V}_i^\top\boldsymbol{AV}_i\boldsymbol{y}_i.

and thus

.. math:: \boldsymbol{y}_i=-\left(\boldsymbol{V}_i^\top\boldsymbol{AV}_i\right)^{-1}\boldsymbol{V}_i^\top\boldsymbol{g}_0.

Finally an approximation is

.. math::

   \boldsymbol{x}_{i+1} = \boldsymbol{x}_0 - \boldsymbol{V}_i\left(\boldsymbol{V}_i^\top\boldsymbol{AV}_i\right)^{-1}\boldsymbol{V}_i^\top\boldsymbol{g}_0 .

During iterative process norm of gradient is compared with given
tolerance :math:`\varepsilon`, and it is stopped when
:math:`\|\boldsymbol{g}_{i+1}\| < \varepsilon`. Basis
:math:`\boldsymbol{V}_i` is generated simultaneously during iterations
by Lanczos algorithm, and it can be written as following equality

.. math::

   \boldsymbol{AV}_i= \boldsymbol{V}_{i+1}\boldsymbol{H}_{i+1,i}.

Matrix :math:`\boldsymbol{V}_{i}` contains vectors of orthogonal basis,
and :math:`\boldsymbol{H}_{i+1,i}` is tridiagonal and rectangular matrix

.. math::

   \boldsymbol{H}_{i+1,i}=
     \left(
       \begin{array}{ccccc}
         h_{0,0} & h_{0,1}&   0    &         &            \\
         h_{1,0} & h_{1,1}& h_{1,2}&   0     &            \\
            0    & h_{2,1}& h_{2,2}& \ddots  &      0     \\
                 &   0    & h_{3,2}& \ddots  &  h_{i-1,i} \\
                 &        &    0   & \ddots  &  h_{i,i}   \\
                 &        &        &    0    &  h_{i+1,i}
       \end{array}
     \right).

Deleting last row we get matrix :math:`\boldsymbol{H}_{i,i}` which is
SPD. Such matrix can be obtained from the equality
([eq:AVi:sub:`V`\ ip1Hip1i]) multiplying by
:math:`\boldsymbol{V}_i^\top` from the left

.. math:: \boldsymbol{V}_i^\top\boldsymbol{A}\boldsymbol{V}_i=\boldsymbol{H}_{i,i}.

In the variant of Lanczos algorithm with full orthogonalization in
:math:`i-th` iteration it generates equality

.. math::

   \boldsymbol{Av}_i= \boldsymbol{v}_{0}h_{0,i} +
     \boldsymbol{v}_{1}h_{1,i} + \cdots
   \boldsymbol{v}_{i+1}h_{i+1,i}.

First basis vector :math:`\boldsymbol{v}_0` is just normed initial
gradient

.. math:: \boldsymbol{v}_0= {\boldsymbol{g}_{0} \over \|\boldsymbol{g}_{0}\|}

therefore the process can continue with

.. math:: \boldsymbol{Av}_0= \boldsymbol{v}_{0}h_{0,0} + \boldsymbol{v}_{1}h_{1,0},

and first entry of matrix :math:`\boldsymbol{H}` is

.. math:: \boldsymbol{v}^\top_0\boldsymbol{Av}_0=h_{0,0}.

Next direction :math:`\boldsymbol{v}_{1}` is already given

.. math:: \boldsymbol{v}_{1}h_{1,0} = \boldsymbol{Av}_0-\boldsymbol{v}_{0}h_{0,0},

and the constant :math:`h_{1,0}` is just

.. math:: h_{1,0} = \|\boldsymbol{Av}_0-\boldsymbol{v}_{0}h_{0,0}\|.

Generally in :math:`i`-th iteration

.. math::

   \boldsymbol{Av}_i= \boldsymbol{v}_{0}h_{0,i} + \boldsymbol{v}_{1}h_{1,i} + \cdots +
                    \boldsymbol{v}_{i}h_{i,i} + \boldsymbol{v}_{i+1}h_{i+1,i}

and

.. math:: h_{j,i} = \boldsymbol{v}_j^\top\boldsymbol{Av}_i, ~j = 1,2,\cdots,i

.. math:: \boldsymbol{v}_{i+1}h_{i+1,i} =  \boldsymbol{Av}_i-\sum\limits_{j=0}^{i}\boldsymbol{v}_{j}h_{j,i},

| [H] initialization; :math:`i=0,~\boldsymbol{x}_{0}`
| :math:`\boldsymbol{g}_{0}=\boldsymbol{Ax}_{0} - \boldsymbol{b} `
| :math:`\boldsymbol{v}_{0}={\boldsymbol{g}_{0} \over {\|\boldsymbol{g}_{0}\|}} `
|  :math:`\boldsymbol{t}=\boldsymbol{Av}_i `
 :math:`h_{i+1,i} = \|\boldsymbol{t}\|`
:math:`\boldsymbol{v}_{i+1} = {\boldsymbol{t} / h_{i+1,i}}`
:math:`\boldsymbol{y}_i=-\boldsymbol{H}_{i,i}^{-1}\boldsymbol{V}_i^\top\boldsymbol{g}_0`
:math:` \boldsymbol{g}_{i+1}=\boldsymbol{g}_0 + \boldsymbol{AV}_i\boldsymbol{y}_i`
:math:`i=i+1`
:math:`\boldsymbol{x} \approx \boldsymbol{x}_0 + \boldsymbol{x}_i`

Preconditioned variant
~~~~~~~~~~~~~~~~~~~~~~

Main terms of Lanczos algorithm used for solving of preconditioned
system ([eq:Ax:sub:`bp`\ recond]) are in table below.

.. math::

   \begin{array}{|c|c|c|}
       \hline
       \boldsymbol{1.} & \boldsymbol{2.} & \boldsymbol{3.} \\
       \hat{\boldsymbol{g}}_{i}=\hat{\boldsymbol{A}}\hat{\boldsymbol{x}}_i-\hat{\boldsymbol{b}} &
       \boldsymbol{L}^{-1}\boldsymbol{g}_{i}=\boldsymbol{L}^{-1} \left(\boldsymbol{A}\boldsymbol{x}_i-\boldsymbol{b}\right) &
       \boldsymbol{g}_{i}=\boldsymbol{A}\boldsymbol{x}_i-\boldsymbol{b}\\
       \hat{\boldsymbol{x}}_{i+1}=\hat{\boldsymbol{x}}_{0} + \hat{\boldsymbol{V}}_i\hat{\boldsymbol{y}}_i   &  
       \boldsymbol{L}^\top{\boldsymbol{x}}_{i+1}=\boldsymbol{L}^\top \left({\boldsymbol{x}}_{0} +
               {\boldsymbol{V}}_i\hat{\boldsymbol{y}}_i \right) &
       \boldsymbol{x}_{i+1}=\boldsymbol{x}_{i} + \boldsymbol{V}_i\boldsymbol{y}_i \\
       \hat{\boldsymbol{A}} \hat{\boldsymbol{V}}_i = \hat{\boldsymbol{V}}_{i+1} \boldsymbol{H}_{i+1,i} &
       \boldsymbol{L}^{-1}\boldsymbol{A} \boldsymbol{V}_i = \boldsymbol{L}^\top\boldsymbol{V}_{i+1} \boldsymbol{H}_{i+1,i} &
       \boldsymbol{A} \boldsymbol{V}_i =   {\boldsymbol{M}}^{-1}\boldsymbol{V}_{i+1} \boldsymbol{H}_{i+1,i}\\
       \hline
     \end{array}

Preconditioner :math:`{\boldsymbol{M}}` is unconventionally
written in inverted form. In the variant of Lanczos algorithm with full
orthogonalization in :math:`i-th` iteration it generates equality. First
basis vector, in case a preconditioner is used, is
:math:`\hat{\boldsymbol{v}}_0`

.. math::

   \begin{array}{rll}
       \hat{\boldsymbol{v}}_0&=& {\hat{\boldsymbol{g}}_{0} \over \|\hat{\boldsymbol{g}}_{0}\|}
   \\
       \boldsymbol{L}^\top\boldsymbol{v}_0&=& 
       {
         {\boldsymbol{L}^{-1}\boldsymbol{g}_{0}} 
           \over 
         {\|\sqrt{\left(\boldsymbol{L}^{-1}\boldsymbol{g}_{0},\boldsymbol{L}^{-1}\boldsymbol{g}_{0}\right) \|}} 
       }
   \\
       \boldsymbol{v}_0&=& 
       {
         {\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}\boldsymbol{g}_{0}} 
           \over 
         {\|\sqrt{\left(\boldsymbol{g}_{0},\boldsymbol{L}^{-T}\boldsymbol{L}^{-1}\boldsymbol{g}_{0}\right) \|}} 
       }
   \\
       \boldsymbol{v}_0&=& 
       {
         {{\boldsymbol{M}}\boldsymbol{g}_{0}} 
           \over 
           {\|\sqrt{\left(\boldsymbol{g}_{0},{\boldsymbol{M}}\boldsymbol{g}_{0}\right) \|}} 
       }
   \\
       \boldsymbol{v}_0&=& 
   {\boldsymbol{M}}\boldsymbol{z}_{0} 
     \end{array}

.. math::

   \begin{array}{c}
         \boldsymbol{z}_{0}= 
         {
         {\boldsymbol{g}_{0}} 
           \over 
           {\|\sqrt{\left(\boldsymbol{g}_{0},{\boldsymbol{M}}\boldsymbol{g}_{0}\right) \|}} 
       }
       \end{array}

.. math::

   \hat{\boldsymbol{A}}\hat{\boldsymbol{v}}_i =
     \hat{\boldsymbol{v}}_{0}h_{0,i} +
     \hat{\boldsymbol{v}}_{1}h_{1,i} + \cdots + 
     \hat{\boldsymbol{v}}_{i+1}h_{i+1,i}

.. math::

   \begin{array}{rll}
     h_{k,i}&=&\hat{\boldsymbol{v}}_{k}^\top\hat{\boldsymbol{A}}\hat{\boldsymbol{v}}_i,~ k\in(0,1,\cdots,i) \\
            &=&\boldsymbol{v}_{k}^\top\boldsymbol{A}\boldsymbol{v}_i
   \end{array}

.. math::

   \begin{array}{rll}
       \hat{\boldsymbol{v}}_{i+1}h_{i+1,i}&=&
       \hat{\boldsymbol{A}}\hat{\boldsymbol{v}}_i
       -\hat{\boldsymbol{v}}_{0}h_{0,i} 
       -\hat{\boldsymbol{v}}_{1}h_{1,i} - \dots  
       -\hat{\boldsymbol{v}}_{i}h_{i,i}  
       \\
       \boldsymbol{L}^\top\boldsymbol{v}_{i+1}h_{i+1,i}&=&
       \boldsymbol{L}^{-1}\boldsymbol{A}\hat{\boldsymbol{v}}_i
      -\boldsymbol{L}^\top\left( \boldsymbol{v}_{0}h_{0,i} 
       -\boldsymbol{v}_{1}h_{1,i} - \dots  
     -\boldsymbol{v}_{i}h_{i,i} \right)
     \\
       \boldsymbol{LL}^\top\boldsymbol{v}_{i+1}h_{i+1,i}&=&
     \boldsymbol{A}\hat{\boldsymbol{v}}_i
       -\boldsymbol{LL}^\top\left(\boldsymbol{v}_{0}h_{0,i} 
       -\boldsymbol{v}_{1}h_{1,i} - \dots  
     -\boldsymbol{v}_{i}h_{i,i} \right)
     \\
       \boldsymbol{M}^{-1}\boldsymbol{v}_{i+1}h_{i+1,i}&=&
     \boldsymbol{A}\hat{\boldsymbol{v}}_i
       -\boldsymbol{M}^{-1}\left(\boldsymbol{v}_{0}h_{0,i} 
       -\boldsymbol{v}_{1}h_{1,i} - \dots  
     -\boldsymbol{v}_{i}h_{i,i} \right)
     \\
       \boldsymbol{z}_{i+1}h_{i+1,i}&=&
     \boldsymbol{A}\hat{\boldsymbol{v}}_i
       -\boldsymbol{z}_{0}h_{0,i} 
       -\boldsymbol{z}_{1}h_{1,i} - \dots  
     -\boldsymbol{z}_{i}h_{i,i} 
     \end{array}

.. math::

   \begin{array}{rll}
       h_{i+1,i} &=&\sqrt {
             \left(
               h_{i+1,i}\hat{\boldsymbol{v}}_{i+1}^\top\hat{\boldsymbol{v}}_{i+1} h_{i+1,i}
            \right)}
            \\
             &=&\sqrt {
             \left(
               h_{i+1,i}\boldsymbol{v}_{i+1}^\top\boldsymbol{L}\boldsymbol{L}^\top\boldsymbol{v}_{i+1} h_{i+1,i}
            \right)}
            \\
             &=&\sqrt {
             \left(
               h_{i+1,i}\boldsymbol{v}_{i+1}^\top\boldsymbol{M}^{-1}\boldsymbol{v}_{i+1} h_{i+1,i}
            \right)}
            \\
             &=&\sqrt {
             \left(
               h_{i+1,i}\boldsymbol{v}_{i+1}^\top\boldsymbol{z}_{i+1} h_{i+1,i}
            \right)}
          \end{array}

.. math::

   \begin{array}{rll}
       \hat{\boldsymbol{A}}\hat{\boldsymbol{v}}_i &=&
       \hat{\boldsymbol{v}}_{0}h_{0,i} +
       \hat{\boldsymbol{v}}_{1}h_{1,i} + \cdots + 
       \hat{\boldsymbol{v}}_{i+1}h_{i+1,i}\\
     \boldsymbol{L}^{-1}\boldsymbol{A}\boldsymbol{L}^{-T}\boldsymbol{L}^\top\boldsymbol{v}_i &=&
     \boldsymbol{L}^\top\boldsymbol{v}_{0}h_{0,i} +
     \boldsymbol{L}^\top\boldsymbol{v}_{1}h_{1,i} + \cdots +
     \boldsymbol{L}^\top\boldsymbol{v}_{i+1}h_{i+1,i}\\
     \boldsymbol{L}^{-1}\boldsymbol{A}\boldsymbol{v}_i &=&
     \boldsymbol{L}^\top\boldsymbol{v}_{0}h_{0,i} +
     \boldsymbol{L}^\top\boldsymbol{v}_{1}h_{1,i} + \cdots +
     \boldsymbol{L}^\top\boldsymbol{v}_{i+1}h_{i+1,i}\\
     \boldsymbol{A}\boldsymbol{v}_i &=&
     \boldsymbol{M}^{-1}\left(\boldsymbol{v}_{0}h_{0,i} +
     \boldsymbol{v}_{1}h_{1,i} + \cdots +
     \boldsymbol{v}_{i+1}h_{i+1,i}\right).
   \end{array}
   