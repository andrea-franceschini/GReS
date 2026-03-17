%% GReS - Input Files
%% Introduction
% This tutorial explains how to pass input in GReS and describes all the
% objects required to configure a single-domain simulation. Details on
% multi-domain simulations are provided in a separate tutorial.

dirpath = fullfile(gres_root,'Tutorial','inputFiles');
cd(dirpath)

%% Input format
% GReS accepts input in three equivalent forms. The user can pass a MATLAB
% |struct|, a path to an XML file, or a key-value cell array:
%
%   s    = struct('key1',val1,'key2',val2,...);   % struct
%   file = 'myFile.xml';                          % XML file path
%   kv   = {'key1',val1,'key2',val2,...};         % key-value cell array
%
% All three forms carry the same information and are interchangeable. In
% particular, a struct and a key-value array are related by
%
%   s == struct(kv{:});
%
% while a struct and an XML file are related by
%
%   s == readstruct(file, AttributeSuffix="");
%
% The conversion from any of the three forms into a struct is handled
% internally by |readInput|, which is called by every GReS constructor:

help readInput

%%
% In what follows we show how to create each GReS object using the struct
% or key-value form. The equivalent XML representation is shown alongside
% each object so the reader can recognise the correspondence between the
% two.

%% Simulation parameters
% The |SimulationParameters| object stores global time-stepping and solver
% settings that are shared by all domains and interfaces in a simulation.
% It is typically the first object to be created, and its content does not
% depend on the spatial discretisation.
%
% The most common parameters are passed directly as struct fields:

input = struct('Start',0.0,'End',10.01,'DtInit',1e-1,'DtMax',1e1,'DtMin',1e-2,'incrementFactor',1.0);
simparams = SimulationParameters(input)

%%
% Any field that is not explicitly provided is automatically set to a
% default value by the constructor. The full set of parameters is:
%
% * |Start| / |End|: initial and final simulation time.
% * |DtInit|: initial time step size.
% * |DtMax| / |DtMin|: upper and lower bounds on the adaptive time step.
% * |incrementFactor|: factor by which the time step is multiplied after a
%   successful nonlinear solve (should be > 1).
% * |choppingFactor|: factor by which the time step is reduced after a
%   failed step and backstep (should be < 1).
% * |RelativeTolerance|: nonlinear convergence criterion measured relative
%   to the norm of the initial residual.
% * |AbsoluteTolerance|: nonlinear convergence criterion on the absolute
%   residual norm.
% * |MaxNLIteration|: maximum number of nonlinear iterations allowed before
%   the solver declares failure and performs a backstep.
% * |MaxConfigurationIteration|: maximum iterations of an outer
%   configuration loop, such as the active-set iteration in contact
%   mechanics problems.
%
% The XML equivalent uses the same names as attributes of a single
% self-closing element. Omitted attributes still fall back to defaults:
%
%   <SimulationParameters
%       Start            = "0.0"
%       End              = "10.01"
%       DtInit           = "1e-1"
%       DtMax            = "1e1"
%       DtMin            = "1e-2"
%       incrementFactor  = "1.0"
%   />

simparamsXML = SimulationParameters('simparam.xml')

%% Materials
% The |Materials| object collects all material definitions for a domain.
% Because a simulation may involve multiple solids with different mechanical
% or hydrogeological properties — and at most one fluid phase — the object
% is built up incrementally rather than in a single call.
%
% We start by creating an empty container:

mat = Materials();

%%
% *Solid material.* Each solid is identified by a unique |name| and is
% associated with one or more mesh regions through |cellTags|. The cell
% tags must be disjoint across different solids, since each mesh cell can
% belong to only one material.

mat.addSolid('name',"sand",'cellTags',1);

%%
% *Constitutive law.* The mechanical response of a solid is specified
% separately via |addConstitutiveLaw|. The second argument is the name of
% the constitutive class to instantiate; all subsequent arguments are
% passed to that class constructor as key-value pairs. This design makes it
% straightforward to swap one law for another without touching the rest of
% the input.

mat.addConstitutiveLaw("sand","Elastic",'youngModulus',1e3,'poissonRatio',0.25);

%%
% *Porous-rock properties.* Hydrogeological properties are added through
% |addPorousRock|. The |permeability| field is flexible: a scalar implies
% isotropic permeability, a 3-vector sets the three diagonal entries of the
% permeability tensor, and a 6-vector provides all independent entries in
% Voigt notation for a full symmetric tensor.

mat.addPorousRock("sand","specificWeight",21.0,"permeability",1e-12,"porosity",0.375);

%%
% *Capillary curves.* Retention curves are attached to a solid with
% |addCapillaryCurves|. The |type| argument selects the analytical model
% (here the Mualem-van Genuchten model); the remaining arguments are the
% model parameters.

mat.addCapillaryCurves("sand","type","mualem","beta",2.0,"n",1.0,"kappa",1.0);

%%
% *Fluid.* GReS currently supports a single fluid phase per simulation,
% which is automatically applied to every cell where a flow model is
% active. The fluid is characterised by its dynamic viscosity, specific
% weight, and compressibility.

mat.addFluid('dynamicViscosity',1e-3,'specificWeight',9.81,'compressibility',4.4e-7);

%%
% *XML equivalent.* The XML structure mirrors the method-call hierarchy
% above. Each |Solid| block may contain nested |Constitutive|, |PorousRock|
% and |Curves| sub-blocks; attribute names are identical to the key-value
% arguments used in the MATLAB calls. Note that |Fluid| must appear before
% any |Solid| block, because porous-rock properties are computed using
% fluid data at construction time.
%
%   <Materials>
%     <Fluid
%         dynamicViscosity = "1e-3"
%         specificWeight   = "9.81"
%         compressibility  = "4.4e-7"/>
%     <Solid name="sand" cellTags="1">
%       <Constitutive>
%         <Elastic youngModulus="1e3" poissonRatio="0.25"/>
%       </Constitutive>
%       <PorousRock
%           specificWeight = "21.0"
%           permeability   = "1e-12"
%           porosity       = "0.375"/>
%       <Curves type="mualem" beta="2.0" n="1.0" kappa="1.0"/>
%     </Solid>
%   </Materials>

matXML = Materials('materials.xml');

%% Boundary Conditions
% The |Boundaries| object manages all boundary conditions assigned to a
% single domain. It must be initialised with the domain grid, which is
% needed to resolve entity lists and build influence maps.

Boundaries(grid);

%%
% *Adding a BC.* Individual boundary conditions are registered one at a
% time with |addBC|. Each BC requires the following fields:
%
% * |name|: a unique string identifier used to reference the BC elsewhere.
% * |type|: the mathematical nature of the constraint. |"dirichlet"|
%   enforces an essential (strong) condition; |"neumann"| applies a natural
%   (weak) condition. Custom types are also supported and are recognised at
%   the physics-solver level.
% * |field|: the type of mesh entity that carries the BC value — the
%   so-called _source_ entity (e.g., |"surface"| or |"cell"|).
% * |variable|: the name of the physical variable being constrained (e.g.,
%   |"displacements"| or |"pressure"|). This name is not validated at
%   construction time; it is checked later by the physics solver.
% * |entityListType|: how the source entities are identified. Use |"tags"|
%   for mesh region tags, |"bcList"| for an explicit array of entity IDs,
%   or |"bcListFile"| for a path to a text file listing the IDs.
% * |entityList|: the actual tags or IDs, interpreted according to
%   |entityListType|.
% * |components| _(optional)_: for vector variables, the subset of
%   components to constrain — e.g., |[1 2 3]| or equivalently |"x,y,z"|.
%   If omitted, the BC applies to all components.

bc.addBC('name',"firstBC",...
        'type',"dirichlet",...
        'field',"surface",...
        'entityListType',"tags",...
        'entityList',1,...
        'variable',"displacements",...
        'components',[1 2 3])

%%
% *Source vs. target entities.* GReS makes a distinction between the
% _source_ entities, where the BC value is prescribed (determined by
% |field| and |entityList|), and the _target_ entities, where the
% constrained degree of freedom actually lives. The two sets can differ —
% for instance, a BC defined on a surface may constrain DOFs at nodes.
% GReS builds a geometric influence map between the two sets, but this
% requires knowing where the DOFs are located. The map is therefore
% computed in a separate step, once the DOF type is known:

bc.computeTargetEntities("firstBC", entityField.node);
M = bc.getEntitiesInfluence("firstBC");

%%
% When source and target entities coincide, the influence map reduces to
% the identity matrix.

%%
% *BC events.* A BC event associates a value with a specific time instant.
% When multiple events are registered for the same BC, GReS constructs a
% piecewise-linear time history and interpolates between events during the
% simulation. The two fields of a BC event are:
%
% * |time|: the time instant at which the event occurs. Alternatively, a
%   function handle |f(t)| can be provided; it will be evaluated at each
%   time step and used as a scalar multiplier on |value|.
% * |value|: the BC value at that instant. It can be a scalar (uniform
%   across all entities), a per-entity array, a path to a file of values,
%   or a spatial function handle |f(x,y,z)|.

bc.addBCEvent("firstBC",'time',1,'value',0.0);
bc.addBCEvent("firstBC",'time',3,'value',1.0);
bc.addBCEvent("firstBC",'time',6,'value',-2.0);

%%
% *XML equivalent.* Each |BC| block corresponds to one |addBC| call; each
% |BCevent| child element corresponds to one |addBCEvent| call. Attribute
% names are identical to the key-value arguments:
%
%   <BoundaryConditions>
%     <BC name="firstBC" type="Dirichlet" field="surface"
%         variable="displacements"
%         entityListType="tags" entityList="1"
%         components="1,2,3">
%       <BCevent time="1" value="0.0"/>
%       <BCevent time="3" value="1.0"/>
%       <BCevent time="6" value="-2.0"/>
%     </BC>
%   </BoundaryConditions>
%
% Further examples:
%
% _Dirichlet constraint loaded from file, all components fixed:_
%
%   <BC name="lateralFix" type="Dirichlet" field="surface"
%       variable="displacements"
%       entityListType="bcListFile" entityList="Input/listLateralFix">
%     <BCevent time="0" value="0"/>
%   </BC>
%
% _Time-varying Neumann flux on an explicit list of cells:_
%
%   <BC name="FluxCell" type="Neumann" field="cell"
%       variable="pressure"
%       entityListType="bcList" entityList="1,2,3,4">
%     <BCevent time="0" value="0"/>
%     <BCevent time="1" value="-0.360"/>
%     <BCevent time="6" value="-0.360"/>
%     <BCevent time="7" value="0"/>
%   </BC>

%% Outputs
% Output settings control what GReS writes to disk and when. The
% |printTimes| field defines the simulation times at which output is
% produced; it accepts either a comma-separated list of values or a path
% to a text file with one time per line. If a print time falls inside a
% time step, GReS linearly interpolates between the solutions at the
% bracketing time levels.
%
% * |outputFile|: base path for ParaView-compatible output. GReS writes a
%   |.pvd| index file and a subfolder of |.vtk| files, one per print time.
%   If this field is omitted, no visualisation output is produced.
% * |matFileName|: base path for a |.mat| file that accumulates history
%   quantities — time series of selected variables collected during the run.
% * |saveHistory|: set to |"1"| to enable history output or |"0"| to
%   disable it.
%
%   <Output
%       printTimes  = "15,30,60,90,120,180"
%       outputFile  = "Output/results"
%       matFileName = "Output/results"
%       saveHistory = "1"
%   />
%
% The output folder is created automatically if it does not already exist.

%% Discretizer
% The |Discretizer| is the central object of a single-domain simulation.
% It acts as a container that brings together the grid, the materials, the
% boundary conditions, and one or more physics solvers into a coherent
% unit. Having all domain data in one place allows the physics solvers to
% access geometry, material properties and constraints through a single
% handle, rather than through a collection of separate arguments.
%
% A Discretizer is instantiated with key-value input:

domain = Discretizer('Boundaries', bound, ...
                     'OutState',   printUtils, ...
                     'Materials',  mat, ...
                     'Grid',       grid);

%%
% Note that |SimulationParameters| and |OutState| are not passed to the
% Discretizer directly. Because they govern the global simulation loop
% rather than a specific domain, they are shared across all domains and
% interfaces and are managed at the top-level simulation object.

%% PhysicsSolver
% Up to this point the model is physics-agnostic: we have defined a
% geometry, material properties and boundary constraints, but we have not
% yet specified which equations GReS should solve. This is done by adding
% one or more |PhysicsSolver| objects to the Discretizer.
%
% Each PhysicsSolver is responsible for three things: registering the
% variable fields it owns in the |DoFManager|, assembling its contribution
% to the global Jacobian and residual, and writing its fields to the output
% at each time step.
%
% To add a solver to the domain, call |addPhysicsSolver| with the class
% name of the solver followed by any solver-specific parameters:

domain.addPhysicsSolver("Poromechanics",'targetRegions',1);

%%
% The same operation in XML is expressed through a |Solver| block, where
% the tag name of each child element must match the class name of the
% corresponding solver:
%
%   <Solver>
%       <Poromechanics targetRegions="1"/>
%   </Solver>
%
% *Important:* a variable field can only be registered once per domain,
% regardless of how many solvers or target regions are involved. The
% following block, for example, would raise an error because
% |"displacements"| would be registered by |Poromechanics| twice:
%
%   <Solver>
%     <Poromechanics targetRegions="1"/>
%     <Poromechanics targetRegions="3,4"/>
%   </Solver>
%
% To see all physics solvers currently available in GReS:

listPhysicsSolvers

%% Working with a single XML file
% In practice it is convenient to collect all the blocks described above
% into a single XML file, and to use the |fileName| attribute to delegate
% verbose sections (such as materials or boundary conditions) to separate
% files. This keeps the main input file readable while still allowing
% fine-grained control over individual components.
%
%   <Problem>
%     <SimulationParameters fileName="Input/simParam.xml"/>
%     <Output
%         outputFile  = "Output/Results"
%         matFileName = "Output/results"
%         printTimes  = "1,4,8"/>
%     <Domain>
%       <Geometry meshFile="Input/ReservoirTest_Hexa.msh"/>
%       <Gauss nGP="2"/>
%       <Materials          fileName="Input/materials.xml"/>
%       <BoundaryConditions fileName="Input/boundaryConditions.xml"/>
%       <Solver>
%         <BiotFullySaturated>
%           <Poromechanics      targetRegions="1,2,3"/>
%           <SinglePhaseFlowFEM targetRegions="2,3"/>
%         </BiotFullySaturated>
%       </Solver>
%     </Domain>
%   </Problem>
%
% Given a single input file, there are two ways to assemble a GReS
% simulation.
%
% The first approach, better suited to prototyping, is to parse the file
% into a struct and then instantiate each object individually. This gives
% full control over how each object is created and makes it easy to modify
% or inspect intermediate results:
%
%   params = readInput(xmlFileName);
%   bc = Boundaries(grid, params.Domain(1).BoundaryConditions);
%
% The second approach is to call |buildModel|, which automates the above
% steps and returns a ready-to-use array of |Discretizer| objects, one per
% |Domain| block. If interface solvers are present, the function also
% returns a cell array of |InterfaceSolver| objects:
%
%   domains = buildModel(xmlFileName);
