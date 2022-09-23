GradedDimKLR
============

SageMath_ code for computing the dimensions of the weight spaces of the cyclotomic KLR algebras attached to symmetrisable quivers using the amazing formula

    .. math:

        \dim_q e(i)R^\Lambda_\alpha e(j)
          = \sum_{w\in\Sym(i,j)}
              \prod_{t=1}^n[N^\Lambda(w,i,t)}_{i_t}] q_{i_i}^{N^{\Lambda,i,t)-1}}

from the paper HuShi_, "Graded dimensions and monomial bases for the cyclotomic quiver Hecke algebras" of Hu and Shi .

The attached code provides a sigle command `klr_cyclotomic_dimension` that
computes the ses the Hu-Shi formula to compute the dimension

Usage
-----

* With a local installation of SageMath_, start  SageMath_ and attach the file `graded_dim_klr` using:

.. code-block:: python

   sage: %attach graded_dim_klr
   sage: klr_cyclotomic_dimension(['A',3],[2], [2,3,3,2,1], [2,3,2,3,1])
   sage: klr_cyclotomic_dimension(['B',3],[2], [2,3,3,2,1], verbose=True)


* Using SageMathCell_, cut-and-paste the code from `graded_dim_klr` into a cell
  and then type the `klr_cyclotomic_dimension` commands into the bottom of the
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
