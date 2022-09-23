GradedDimKLR
============

SageMath_ code for computing the dimensions of the weight spaces of the cyclotomic KLR algebras attached to symmetrisable quivers using the amazing formula

    .. math:

        \dim_q e(i)R^\Lambda_\alpha e(j)
          = \sum_{w\in\Sym(i,j)}
              \prod_{t=1}^n[N^\Lambda(w,i,t)}_{i_t}] q_{i_i}^{N^{\Lambda,i,t)-1}}

from the paper HuShi_, "Graded dimensions and monomial bases for the cyclotomic quiver Hecke algebras" of Hu and Shi .

The attached code provides a sigle command ``klr_cyclotomic_dimension`` that
computes the ses the Hu-Shi formula to compute the dimension. The synaxt of
this command is:

.. code-bHuShi_lock:: python

   klr_cyclotomic_dimension(C, L, i, [j, verbose])

where:

* ``C`` specifies the *Cartan type* of the quiver. This can either be a a list
  of the form ``['A', 3]`` (finite type :math:`A_3`), ``['A', 3, 1]`` (affine type
  :math:`A_3^{(1)}`), ``['B',4]`` (finite type :math:`B_4`), etc.

* ``L`` is a *list* that specifies the dominant weight. For example,
  ``[0,2,2,3]`` represents the dominant weight
  :math:`\Lambda_0+2\Lambda_2+\Lambda_3`.

* ``i`` is an `n`-tuple of vertices of the quiver

* ``j`` is an `n`-tuple of vertices of the quiver. If ``j`` is omitted then
  ``j`` is set equal to ``i``.

* ``verbose`` an optional parameter that wehn set to ``True`` causes additional
  information to from the HuShi_ formula to be printed.

Examples
--------

.. code-block:: python

    sage: klr_cyclotomic_dimension(['D',4],[2], [2,3,4,1], verbose=True)
    Subgroup of permutations = < () >
    N(1,t)-1: 0  0  0  0
    N(w,t):   1  1  1  1
    X(w): 1
    1
    sage: klr_cyclotomic_dimension(['D',4],[2], [2,3,4,1], [2,3,4,1], verbose=True)
    Subgroup of permutations = < () >
    N(1,t)-1: 0  0  0  0
    N(w,t):   1  1  1  1
    X(w): 1
    1
    sage: klr_cyclotomic_dimension(['D',4],[2], [2,3,4,1], [2,4,3,1], verbose=True)
    Subgroup of permutations = < (1,2) >
    N(1,t)-1: 0  0  0  0
    N(w,t):   1  1  1  1
    X(w): 1
    1
    sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,1,2], verbose=True)
    Subgroup of permutations = < () >
    N(1,t)-1: 0  0  1
    N(w,t):   1  1  2
    X(w): q^2 + 1
    q^2 + 1
    sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,2,1], verbose=True)
    Subgroup of permutations = < () >
    N(1,t)-1: 0  0  1
    N(w,t):   1  1  2
    X(w): q^2 + 1
    q^2 + 1
    sage:
    sage: klr_cyclotomic_dimension(['A',2,1],[0], [0,1,2], [0,2,1], verbose=True)
    Subgroup of permutations = < (1,2) >
    N(1,t)-1: 0  0  1
    N(w,t):   1  1  1
    X(w): q
    q
    sage: klr_cyclotomic_dimension(['A',1],[1,1],[1],[1], verbose=True)
    Subgroup of permutations = < () >
    N(1,t)-1: 1
    N(w,t):   2
    X(w): q^2 + 1
    q^2 + 1
    sage: klr_cyclotomic_dimension(['A',1],[1,1],[1],[1], verbose=True)
    Subgroup of permutations = < () >
    N(1,t)-1: 1
    N(w,t):   2
    X(w): q^2 + 1
    q^2 + 1
    sage: klr_cyclotomic_dimension(['B',3],[2], [2,3,3,2,1], verbose=True)
    Subgroup of permutations = < (), (1,2), (0,3), (0,3)(1,2) >
    N(1,t)-1: 0  1 -1  0  1
    N(w,t):   1  2  0  1  2
    X(w): 0
    N(w,t):   1  2  2  1  2
    X(w): (q^4 + 1)*(q^2 + 1)^2/q^2
    N(w,t):   1  0 -2  1  2
    X(w): 0
    N(w,t):   1  0  0  1  2
    X(w): 0
    (q^4 + 1)*(q^2 + 1)^2/q^2

Usage
-----

* With a local installation of SageMath_, start  SageMath_ and attach the file ``graded_dim_klr`` using:

.. code-block:: python

   sage: %attach graded_dim_klr
   sage: klr_cyclotomic_dimension(['A',3],[2], [2,3,3,2,1], [2,3,2,3,1])
   sage: klr_cyclotomic_dimension(['B',3],[2], [2,3,3,2,1], verbose=True)


* Using SageMathCell_, cut-and-paste the code from ``graded_dim_klr`` into a cell
  and then type the ``klr_cyclotomic_dimension`` commands into the bottom of the
  cell.

.. code-block:: python

   klr_cyclotomic_dimension(['A',3],[2], [2,3,3,2,1], [2,3,2,3,1])
   klr_cyclotomic_dimension(['B',3],[2], [2,3,3,2,1], verbose=True)


Andrew Mathas
Copyright (C) 2022

------------

GNU General Public License, Version 3, 29 June 2007

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License (GPL_) as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

.. _GPL: http://www.gnu.org/licenses/gpl.html
.. _HuShi: https://arxiv.org/abs/2108.05508G
.. _SageMath: https://www.sagemath.org/
.. _SageMathCell: https://sagecell.sagemath.org/
