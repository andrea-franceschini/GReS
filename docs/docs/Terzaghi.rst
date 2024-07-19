Terzaghi benchmark
===================

Introduction
-------------------
This is the simplest test case to get in touch with most of the features 
currently available in GReS.

.. figure:: ../../Tests_Release/Terzaghi_Biot/Images/Setting_Terzaghi.png
   :align: center
   :width: 500
   :figclass: align-center

.. mat:currentmodule:: read

Mesh
----------------------
The mesh can be read in `vtk` or `gmsh` creating an object
of :class:`Mesh` class:

.. mat:class:: Mesh(obj,fileName)