Terzaghi Benchmark
==================

Introduction
------------

This tutorial describes how to set up a simple model for Terzaghi Benchmark.

This benchamark solves coupled poromechanics equations.

Mesh
------------
GReS is able to read .VTK and .msh files. First create an object of the ``Mesh()`` class.
Then call either the method ``readGMSHmesh`` or ``readVTKmesh`` depending on your mesh file format.  

