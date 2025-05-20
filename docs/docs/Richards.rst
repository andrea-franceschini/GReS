Richards benchmark
===================

Introduction
-------------------
This is the simplest test case to get in touch with how to describe
unsaturated flow in a porous media using Richards modulo currently available
in GReS.
This benchmark has been organized as short overview of the model,
followed by two experiments.





- **Governing equations**

The mass balance equation for a single phase flow on a unsaturated porous
media is given by

.. math:: 
    \sigma(S) \dfrac{\partial p}{\partial t} + \nabla \cdot \bigg( \lambda(S) K \nabla ( p + \gamma z)\bigg)= q

where :math:`\sigma(S)` is the store coefficient, :math:`\lambda(S)`
is the fluid mobility, :math:`K` is the permeability tensor and :math:`q`
is the volumetric flux.
The store coefficient is given by

.. math:: 
    \sigma(S) = S\bigg(\alpha+\beta\phi\bigg) + \phi\dfrac{dS}{dp}

* :math:`S` is the fluid saturation. 
* :math:`\alpha` is the rock compressibility.
* :math:`\beta` is the fluid compressibility.
* :math:`\phi` is the porosity of the media.
  
Meanwhile, the saturation is described as

.. math:: 
    S_{e} = \dfrac{S-S_{r}}{S_{s}-S_{r}}

* :math:`S_{e}` is the effective saturation.
* :math:`S_{r}` is the residual saturation.
* :math:`S_{s}` is the maximum saturation for the fluid.

The fluid mobility is given by

.. math::
    \lambda(S) = \dfrac{k_{r}(S)}{\mu}

* :math:`k_{r}(S)` is relative permeability of the fluid.
* :math:`\mu` is the dynamic viscosity of the fluid.


The variables :math:`S_{e}` and :math:`k_{r}` depend on pressure, and the
code is designed to handle their relationships in two different ways: through
tabulated data or via a function, the van Genuchten model [VanG80]_,
[VanG85]_, [MRST]_.





- **Tabular Data**

In this case, the variables :math:`S_{e}` and :math:`k_{r}` are store
in a table that contains the information about the variable and it's correspondent
pressure (see `fig:Tabular Curves`_). 

.. _fig:Tabular Curves:

    .. container:: figures

        .. image:: figs/Richards/Tab_CapCurve200.png
            :width: 49%
            :name: CapilaryCurve

        .. image:: figs/Richards/Tab_RelPerm200.png
            :width: 49%
            :name: RelativeCurve


    **Tabular Curves**: Example of tabular curves for the effective saturation and the
    relative permeability.



.. class:: TabularCurve

   Class to represent the :math:`S_{e}` or :math:`k_{r}` variable and
   it's dependence on pressure.

   :ivar tabW: a list of points with the values of :math:`S_{e}` or :math:`k_{r}` and it's respective pressure associated.
   :ivar derivW: a list of first derivative, computed using foward differenciating, for the tabW values.
   :ivar derivW2: a list of second derivative, computed using foward differenciating, for the derivW values.
   :ivar nPoints: number of points descibring tabW.

   In this object class, the :math:`S_{e}` or :math:`k_{r}` value and it's derivatives
   at a specific pressure is obtained through linear interpolations of the point
   lists (tabW, derivW, derivW2).

.. note:: 
 - The data stored in this class assumes that the pressure is strictly positive.
 - The forward differentiation process creates a list with one fewer point than the original.
 - The pressure position of this derivative is at the midpoint between the points used in thi process.
 - Two extrapolation regions are always necessary due to this differentiation process.
 - One region covers pressures greater than the maximum defined in the list.
 - The other covers values between zero and the first defined pressure value.
 - It is recommended to be cautious with these extrapolation regions, as they can interfere with the convergence of the non-linear solver.





- **Van Genuchten Model**

The relationship between capillary pressure and water saturation is

.. math:: 
    S_{e} = \bigg(1+(\epsilon p)^{n}\bigg)^{-m}

* :math:`\epsilon` is related to the average size of pores
* :math:`n` and :math:`m` are experimental parameters.

The relative permeability be given by the van Genuchten–Mualem model,

.. math:: 
    k_{r} = S_{e}^{k}\bigg[1-\bigg(1 - S_{e}^{1/m}\bigg)^{m}\bigg]^{2},
    \quad m=1-1/n

or by the van Genuchten-Burdine

.. math:: 
    k_{r} = S_{e}^{2}\bigg[1-\bigg(1 - S_{e}^{1/m}\bigg)^{m}\bigg],
    \quad m=1-2/n

and the :math:`k` is a connectivity factor, usually 0.5.




Experiment 1
-------------------

This experiment consist in a simulation of a pressure dropping inside a column. 
The boundary and initial conditions are describe in `fig:Exp 1 - Conditions`_

.. _fig:Exp 1 - Conditions:

    .. image:: figs/Richards/Exp1_Boundary.png
        :align: center
        :scale: 50 %

    **Conditions**: Initial and Boundary condition for this experiment.

The proprieties for this model are

+------------------------+-----------+----------+
|                        | Value     | unit     |
+========================+===========+==========+
| Dimensions             | 1x1x10    | m        |
+------------------------+-----------+----------+
| Partitions             | 4x4x40    | cells    |
+------------------------+-----------+----------+
| Permeability           | 1e-13     | m²       |
+------------------------+-----------+----------+
| Porosity               | 0.3       | adm      |
+------------------------+-----------+----------+
| Rock Specific Weight   | 21        | KPa      |
+------------------------+-----------+----------+
| Water Specific Weight  | 9.81      | KPa      |
+------------------------+-----------+----------+
| Water Compressibility  | 0         | KPa      |
+------------------------+-----------+----------+
| Viscosity              | 1.157e-11 | KPa m    |
+------------------------+-----------+----------+
| Gravity                | 9.81      | m/s²     |
+------------------------+-----------+----------+
| Residual Saturation    | 0         | adm      |
+------------------------+-----------+----------+
| Maximum Saturation     | 1         | adm      |
+------------------------+-----------+----------+

and :math:`S_{e}` and :math:`k_{r}` are given and two different forms, by a 
table with the files "pcCurve*.dat" and "KrCurve*.dat", and by the proposed
function, where :math:`\epsilon = 0.3592\ (1/KPa)`, :math:`k = 0.5`,
:math:`n = 3.1769` and :math:`m = 0.6852`.

The profile for capillary pressure and water saturation for the initial condition is
describe in the `fig:Exp 1 - Profile`_.

.. _fig:Exp 1 - Profile:
    
    .. figure:: figs/Richards/Exp1_Perfil.png
        :align: center
        :scale: 75 %

    **Profile for the initial condition**: Illustrated the profile for the capillary 
    pressure and the saturation in the initial condition.


With this settings, the script "Main.m" run simulation of the Richards equations
for this model.

.. _fig:Exp 1 - Results:

    .. container:: figures

        .. image:: figs/Richards/Exp1_pressure.png
            :width: 49%

        .. image:: figs/Richards/Exp1_saturation.png
            :width: 49%

    **Results**: admensional pressure and saturation fields for three different time 
    instance.

In `fig:Exp 1 - Results`_, we observe a pressure drop throughout the entire domain
as time progresses.
This leads to a decrease in saturation, as also illustrated in the figure.


To validate the two approaches for representing :math:`S_{e}` and :math:`k_{r}`, using
either a table or an function, we executed the script "Main.m" for both cases under
identical conditions.
The simulation results were saved and can be compared by running the script
"Validation.m", which provides a direct visualization of the differences
between the two methods.

.. _fig:Exp 1 - Validation:

    .. container:: figures

        .. image:: figs/Richards/Exp1_pressure_validation.png
            :width: 49%

        .. image:: figs/Richards/Exp1_saturation_validation.png
            :width: 49%

    **Validation**: Comparison between the solution using a tabular curve and a
    analytical function for the variables :math:`S_{e}` and :math:`k_{r}`.

In `fig:Exp 1 - Validation`_, we observe an agreement between the results, which
was expected and validate both methods.

.. note:: 
  - The number of points is a important factor for the performance of the non-linear solver, as demonstrated in this experiment. When using 200 points to describe the variables :math:`S_{e}` and :math:`k_{r}` under the conditions defined for this study, it becomes necessary to reduce the time step as the simulation progresses. However, this adjustment is not required when using tables with 2000 points.
  - To compare the solution obtained using table and function, it's necessary to impose the same conditions.  However, achieving a meaningful comparison is necessary to use the same time steps in the non-linear solver.




Experiment 2
-------------------

This second experiment compares methodology code in GReS with one present
in [Varela]_ for the MRST code.
The simulation comprises of one-dimensional soil block with an
residual moisture content and a total length of 1 m.
A pressure head is imposed at the top boundary, initiating downward
water infiltration.
As a result, an infiltration front develops and progresses through the
block over time.

The proprieties for this model are describe as

+-----------------------+------------+-------+
|                       | Value      | unit  |
+=======================+============+=======+
| Dimensions            | 1x1x1      | m     |
+-----------------------+------------+-------+
| Partitions            | 1x1x30     | cells |
+-----------------------+------------+-------+
| Permeability          | 9.4018e-12 | m²    |
+-----------------------+------------+-------+
| Porosity              | 1.0        | adm   |
+-----------------------+------------+-------+
| Rock Specific Weight  | 2.1e4      | N/m³  |
+-----------------------+------------+-------+
| Water Specific Weight | 9.8066e3   | N/m³  |
+-----------------------+------------+-------+
| Water Compressibility | 0          | N/m³  |
+-----------------------+------------+-------+
| Viscosity             | 1e-3       | mPa s |
+-----------------------+------------+-------+
| Gravity               | 9.8066     | m/s²  |
+-----------------------+------------+-------+
| Residual Saturation   | 0.102      | adm   |
+-----------------------+------------+-------+
| Maximum Saturation    | 0.368      | adm   |
+-----------------------+------------+-------+

where the parameters to Van Genuchten model are: :math:`\epsilon = 0.3592\ (1/KPa)`,
:math:`k = 0.5`, :math:`n = 3.1769` and :math:`m = 0.6852`.

No-flow (Neumann) boundary conditions are imposed along the lateral sides of
the block and Dirichlet boundary conditions are applied at the top and bottom
boundaries, as detailed in Table

+-----------------------------+-------+------+
|                             | Value | unit |
+=============================+=======+======+
| Boundary condition (Top)    | -0.75 | m    |
+-----------------------------+-------+------+
| Boundary condition (Bottom) | -10.0 | m    |
+-----------------------------+-------+------+
| Initial condition (Top)     | -10.0 | m    |
+-----------------------------+-------+------+


The simulations are run using identical conditions for both cases, and the results
of a 3 days simulation is show in `fig:Exp 2 - Validation`_.

.. _fig:Exp 2 - Validation:

    .. container:: figures

        .. image:: figs/Richards/Exp2_GReSxMRST_head_pressure.png
            :width: 49%

        .. image:: figs/Richards/Exp2_GReSxMRST_saturation.png
            :width: 49%

    **Validation**: Comparison between the solution using GReS an MRST

The curves illustrate the evolution of pressure head, capturing the progression
of the infiltration process and the associated changes in saturation within the
block.
A high degree of similarity is observed between the results, indicating that
the GReS model produces outcomes comparable to those of the MRST model.

Reference
-------------------

.. [VanG80] 
    M. Th. van Genuchten.
    A Closed-form Equation for Predicting the Hydraulic Conductivity
    of Unsaturated Soils.
    In: Soil Science Society of America Journal 44.5 (1980), pp. 892–898.
    doi: https://doi.org/10.2136/sssaj1980.03615995004400050002x

.. [VanG85]
    Martinus Van Genuchten and D.R. Nielsen.
    On Describing and Predicting the Hydraulic Properties of Unsaturated Soils.
    In: Annales Geophysicae 3 (Jan. 1985), pp. 615–628.

.. [MRST]
    Knut-Andreas Lie.
    An Introduction to Reservoir Simulation Using MATLAB/GNU Octave:
    User Guide for the MATLAB Reservoir Simulation Toolbox (MRST).
    Cambridge University Press, 2019

.. [Varela]
    Varela, Jhabriel.
    Implementation of an MPFA/MPSA-FV solver for the unsaturated flow in deformable porous media.  
    (2018).
