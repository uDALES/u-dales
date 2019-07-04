Exmaple information
================

`NAMSTRESS` has been created in
the `namoptions.exp` file which include a switch (lstress) to calculate
all the terms of the stress tensor. It's important to notice that the
subroutine calculates and writes, in ASCII or NetCDF, the file `stressbudget.exp` the budget for the elements $u'^2$, $u'v'$, $u'w'$,
$v'^2$, $v'w'$, and $w'^2$ of the stress tensor. Finally, at the end of
the file the sum of the elements of the diagonal divided by 2 is
written. This corresponds to the TKE budget calculated in subroutine
*modbudget.f90*.

Momentum budget (exmaple)
===============

The prognostic equation for the turbulent fluctuations, $u'_i$, reads
[@stull:introbl chapter 4, eq. 4.1.1]:

$$\frac{\partial u'_i}{\partial t} +  \overline{u}_k\frac{\partial u'_i}{\partial{x_k}} +
                                      u'_k\frac{\partial \overline{u}_i}{\partial{x_k}} +
                                      u'_k \frac{\partial u'_i}{\partial{x_k}} =
  \delta_{i3} \left( \frac{\theta'_v}{\overline{\theta}_v} \right) g +
   f_c \epsilon_{ik3} u'_k -
  \frac{1}{\overline{\rho}} \frac{\partial p'}{\partial x_i} + 
  \left(\nu \frac{\partial^2 u'_i}{\partial x^2_k}\right)' +
  \frac{\partial(\overline{u'_i u'_k})}{\partial x_k}.$$

By using the condition of incompressibility,
$\partial u'_k/ \partial x_k=0$, this equation can be written:

$$\frac{\partial u'_i}{\partial t} + \frac{\partial  \overline{u}_k u'_i}{\partial{x_k}} +
  \frac{\partial u'_k \overline{u}_i}{\partial{x_k}} +
  \frac{\partial u'_k u'_i}{\partial{x_k}} =
  \delta_{i3} \left( \frac{\theta'_v}{\overline{\theta}_v} \right) g +
  f_c \epsilon_{ik3} u'_k -
  \frac{1}{\overline{\rho}} \frac{\partial p'}{\partial x_i} + 
  \left(\nu \frac{\partial^2 u'_i}{\partial x^2_k}\right)'+
  \frac{\partial(\overline{u'_i u'_k})}{\partial x_k}. 
  \label{EQ_TEND_FLUCT}$$


Viscous dissipation
------------------------

In the term
$\overline{u'_j \nu \frac{\partial^2 u'_i}{\partial x^2_k}}$,
$\nu \frac{\partial^2 u'_i}{\partial x^2_k}$ is defined at the
$u_i$-point and is interpolated to the correct edge by averaging in the
$j$-direction. $u'_j$ is available at the $u_j$-point and is
interpolated to the $ij$ edge by averaging in the $i$-direction. Note
that in the LES model, $\nu$ includes (or in fact: *only* includes) the
subgrid viscosity.

The calculation of the `u_term, v_term and w_term` is as follows:


    call diffu(u_term)
    call diffv(v_term)
    call diffw(w_term)

    call cyclicx(u_term)
    call cyclicx(v_term)
    call cyclicx(w_term)

    do k=1,k1
       uterm_avl(k) = sum(u_term(2:i1,2:j1,k))/rslabs
       vterm_avl(k) = sum(v_term(2:i1,2:j1,k))/rslabs
       wterm_avl(k) = sum(w_term(2:i1,2:j1,k))/rslabs
    enddo

    call MPI_ALLREDUCE(uterm_avl, uterm_av, k1,  MY_REAL, &
                           MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(vterm_avl, vterm_av, k1,  MY_REAL, &
                           MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(wterm_avl, wterm_av, k1,  MY_REAL, &
                           MPI_SUM, comm3d, mpierr)

    do k=1,k1
        u_term(:,:,k) = u_term(:,:,k) - uterm_av(k)
        v_term(:,:,k) = v_term(:,:,k) - vterm_av(k)
        w_term(:,:,k) = w_term(:,:,k) - wterm_av(k)
    enddo
