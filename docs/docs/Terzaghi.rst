Terzaghi benchmark
===================

Introduction
-------------------
This is the simplest test case to get in touch with most of the features 
currently available in GReS.
Terzaghi problem is a well problem benchmark which describes the 1D consolidation process of a fully saturated porous media 
column subjected to a sudden constant load.
The governing equations are the Biot consolidation equations, which couple the equilibrium equations for the 
displacement field with the continuity equation for the pore pressure.

Problem Description
--------------------

Consider a homogeneous, isotropic, fully saturated soil layer of thickness \( H \) subjected to a sudden uniform load :math:`p_T` applied at the top surface. 

The domain is a simple 1x1x10 m sample discretized with regular cubes of size :math:`0.5\times0.5\times0.5.`.

.. figure:: ../../Tests_Release/Terzaghi_Biot/Images/Setting_Terzaghi.png
   :align: center
   :width: 500
   :figclass: align-center


- **Governing equations**

1. Equilibrium equation
   
   .. math::
      \nabla \cdot \sigma = - \mathbf{F}

   where :math:`\sigma` is the total stress tensor and :math:`\mathbf F` is a volume force (e.g gravity).

   The constitutive relationship

   .. math:: 
      \sigma = C : \varepsilon - \alpha p \mathbf{I}

   relates the stress tensor to strain :math:`\varepsilon` and pressure :math:`p` according to Biot's theory. 
   Coefficient :math:`\alpha` is Biot coefficient and is typically set equal to :math:`1`.


2. Flow mass conservation
   
   Single phase fully saturated flow is governed by a mass balance equation

   .. math:: 
          \nabla \cdot \bigg( \dfrac{K}{\mu} \nabla p \bigg) - \alpha \dfrac{\partial \epsilon}{\partial t} - \dfrac{1}{M} \dfrac{\partial p}{\partial t} = q

   where :math:`K` is the permeability tensor, :math:`M` is Biot modulus, :math:`\epsilon = tr(\varepsilon)` is the volumetric strain and 
   :math:`q` is the volumetric flux. 



- **Boundary conditions**

   The following boundary conditions are applied: 

   .. math:: 
      \begin{align}
      &\mathbf u = 0      \ \ \  \text{at} z = 0 \\
      &\mathbf u \cdot n = 0   \ \ \     \text{at} \ x = [0,1] \ || \ y = [0,1] \\
      &p = 0 \ \ \ \text{at} \ z = H \\
      &q = 0 \ \ \ \text{at} \ x = [0,1] \ || \ y = [0,1] \ || \ z = 0; 
      \end{align}




- **Initial conditions**
  
   Since the load is assumed to be applied suddenly, inhomogeneous initial conditions apply 
   for both pressure and vertical displacements.

     .. math:: 
      \begin{align}
      & p_0 = \dfrac{\alpha M}{K_u + 4 \mu / 3}p_T \\
      &    u_0(z) = \dfrac{1}{K_u+4\mu/3}p_T(H-z)
      \end{align}


- **Analytical Solution**

   The analytical solution for the Terzaghi problem in terms of pore pressure is given by:

   .. math::
      p(z, t) = \frac{4p_0}{\pi} \sum_{n=0}^{\infty} \frac{1}{2n+1} \exp \left( -\frac{(2n+1)^2 \pi^2 c_v t}{4H^2} \right) \cos \left( \frac{(2n+1) \pi z}{2H} \right)

   While vertical displecement are given by:

   .. math::
      \begin{split}
       u(z,t) = c_M p_0 & \bigg\{ (H-z) - \dfrac{8H}{\pi^2} \sum_{m=0}^\infty \dfrac{1}{(2m+1)^2} \\ 
       &\exp \bigg( \bigg[ \dfrac{-2(m+1)^2 \pi^2 c \,t}{4H^2}\bigg]\bigg) \, \cos \bigg[\dfrac{(2m+1)\pi z}{2H}\bigg]\bigg\}+u_0
      \end{split}




   .. mat:currentmodule:: read




Creating the Mesh
----------------------
The mesh can be read in `vtk` or `gmsh` creating an object
of :class:`Mesh` class:

.. mat:class:: Mesh(obj,fileName)