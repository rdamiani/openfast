.. _UsrPtfm-theory:

User-Platform Theory
++++++==============

This document discusses the theory behind the UsrPtfmLoad module to include Craig-Bampton's reduction capabilities. The
theoretical foundation, numerical tools, and some special handling in
the implementation will be introduced. References will be provided in
each section detailing the theories and numerical tools.

Matrix notation is used to denote vectorial or
vectorial-like quantities. For example, an underline denotes a vector
:math:`\underline{u}`, an over bar denotes unit vector :math:`\bar{n}`,
and a double underline denotes a tensor
:math:`\underline{\underline{\Delta}}`. Note that sometimes the
underlines only denote the dimension of the corresponding matrix.

Overview
------------------

This section focuses on the theory behind the UsrPtfmLoad module.

UsrPtfmLoad relies on a dynamics system reduction
via the Craig-Bampton (C-B) method together with a static-improvement
method (SIM), greatly reducing the number of modes needed to obtain an
accurate solution. Most of the theory is transfered from SubDyn's theory, and the user is referred to that module's documentation :cite:`damianixxx`..

It is common practice in the offshore wind industry to perform sequentially-coupled laods analyses for fixed-bottom wind turbines.
The substructure and foundation designers provide C-B reduced mass, damping, and stiffness matrices together with C-B reduced hydrodynamic and gravitational forces at the interface with the tower. 

The wind turbine designers can then analyze aeroelastic loads on the tower and rotor-nacelle-assembly (RNA) through an aero-servo-elastic code capable of simultaneously processing the C-B information at tower base.
In this fashion, the respective IPs (for the turbine and substructure) do not need to be shared, while still accounting for all the actions on the system (aero-hydro-servo-elastic).
The C-B reduced set of matrices and loads is an expedient toward reducing the large number of DOFs (~ :math:`{10^3}`) associated with a standard
finite-element analysis of a typical multimember structure, thus leading to an increase in computational efficiency during wind turbine system dynamic simulations.
The level of fidelity in the overall system response can be maintained by guaranteeing that the substructure reduced-DOF model retains the fundamental low-frequency response modes of the original finite-element model.

UsrPtfmLoad expects the user to provide (6+m)x(6+m) equivalent mass, damping, and stiffness matrices as obtained via a C-B dyanmic model reduction, where :math:`m` is the number of retained modes.
Additionally, the C-B equivalent forces (:math:`(6+m)x1`) associated with the :math:`m` are also required for every time step to be analyzed.


The following sections discuss the integration of UsrPtfmLoad within the FAST
framework, the theory pertaining to this C-B implementation, but more details can be found in :cite:`damianixxx`. The
state-space formulation to be used in the time-domain simulation is
also presented. 

Integration with the FAST Modularization Framework
--------------------------------------------------

Based on a new modularization framework :cite:`jonkman2013`, FAST joins an
aerodynamics module, a hydrodynamics module, a control and electrical
system (servo) module, and structural-dynamics (elastic) modules to
enable coupled nonlinear aero-hydro-servo-elastic analysis of land-based
and offshore wind turbines in the time domain.  :numref:`flow-chart2` shows the basic
layout of the SubDyn module within the FAST modularization framework.

.. _flow-chart2:

.. figure:: figs/flowchart2.png
   :width: 70%
           
   SubDyn layout within the modularization framework


In the existing loosely coupled time-integration scheme, the glue-code
transfers data at each time step. Such data includes hydrodynamic loads,
substructure response, loads transmitted to the TP, and TP response
among SubDyn, HydroDyn, and ElastoDyn. At the interface nodes, the TP
displacement, rotation, velocity, and acceleration are inputs to SubDyn
from ElastoDyn, and the reaction forces at the TP are outputs of SubDyn
for input to ElastoDyn. SubDyn also outputs the substructure
displacements, velocities, and accelerations for input to HydroDyn to
calculate the hydrodynamic loads that become inputs for SubDyn. In
addition, SubDyn can calculate the member forces, as requested by the
user. Within this scheme, SubDyn tracks its states and integrates its
equations through its own solver.

In a tightly coupled time-integration scheme (yet to be implemented),
SubDyn sets up its own equations, but its states and those of other
modules are tracked and integrated by a solver within the glue-code that
is common to all of the modules.

SubDyn is implemented in a state-space formulation that forms the
equation of motion of the substructure system with physical DOFs at the
boundaries and modal DOFs representing all interior motions. At each
time step, loads and motions are exchanged between modules through the
driver code; the modal responses are calculated inside SubDyn’s
state-space model; and the next time-step responses are calculated by
the SubDyn integrator for loose coupling and the global system
integrator for tight coupling.


.. _UsrPtfmCB:

Dynamic System of Equations and C-B Reduction 
---------------------------------------------

The main equations of motion for the substructure can be found in  cite:`damianixxx` together with the C-B reduction methodology.
Skipping the derivation, the reduced equations of motion can be written as:

.. math:: :label: main4b

        \begin{bmatrix} 
        	\tilde{M}_{BB} & \tilde{M}_{Bm} \\
                \tilde{M}_{mB} & I 
        \end{bmatrix} 
        \begin{bmatrix} 
        	\ddot{U_{TP}} \\ 
                \ddot{q_m} 
        \end{bmatrix} +
        \begin{bmatrix} 
	         	0 & 0 \\
	                0 & 2\zeta \Omega_m 
        \end{bmatrix}
         \begin{bmatrix} 
	        	\dot{U_{TP}} \\ 
	                \dot{q_m} 
        \end{bmatrix} +
        \begin{bmatrix} \tilde{K}_{BB} & 0 \\
			0      & \Omega_m^2
        \end{bmatrix} 
        \begin{bmatrix} 
        	U_{TP} \\ 
                q_m 
        \end{bmatrix} =
        \begin{bmatrix} \tilde{F}_{TP} \\
                        \tilde{F}_m  
                        \end{bmatrix}  
   
with

.. math:: :label: tilde_partitionsb

	\tilde{M}_{BB} = T_I^T \bar{M}_{BB} T_I
	
	\tilde{M}_{Bm} = T_I^T \bar{M}_{Bm}
	
	\tilde{M}_{mB} = \tilde{M}_{Bm}^T 
	
	\tilde{K}_{BB} = T_I^T \bar{K}_{BB} T_I 

	\tilde{F}_{TP} = F_{TP} + T_I^T \bar{F}_{HDR} + T_I^T \bar{F}_{Rg} + T_I^T \bar{\Phi}_{R}^T \left( F_L + F_{Lg} \right)

	\tilde{F}_{m} = \Phi_m^T \left( F_L + F_{Lg} \right)

Equation :eq:`main4b` represents the equations of motion of the substructure after
the C-B reduction. The total DOFs of the substructure are reduced from
(6 x total number of nodes) to (6 + *m*). 

The main objective of the UsrPtfmLoad module is to calculate the output :math:`{-F_{TP}}`,  i.e., the force applied to the TP by the substructure.
      
The user will need to provide the C-B equivalent mass :math:`{M_{usr}}`, damping :math:`{C_{usr}}`, and stiffness :math:`{K_{usr}}` 
represented by the first three matrices in :eq:`main4b`.

Note that :math:`{\tilde{M}_{BB}}` is a (*6*\ ×\ *6*) matrix, :math:`{\tilde{M}_{Bm}}` is a (*6*\ ×\ *m*) matrix, 
:math:`{\tilde{M}_{mB}}` is a (*m*\ ×\ *6*) matrix, :math:`I` is the (*m*\ ×\ *m*) identity matrix; :math:`\Omega_m^2` is a (*m*\ ×\ *m*)  matrix; 
:math:`2\zeta \Omega_m` is a (*m*\ ×\ *m*)  matrix; and :math:`{\tilde{K}_{BB}}` is a (*6*\ ×\ *6*) matrix.  
Therefore :math:`{M_{usr}}`, :math:`{C_{usr}}`, and :math:`{K_{usr}}`  are ((*6+m*)\ ×\ (*6+m*)) matrices.

Additionally, the user will provide ((*6+m*)\ ×\ *1*) forces, :math:`{F_{usr}}`, for every time step. With obvious meaning of the new symbols, :eq:`main4b` can be written as in :eq:`main5b`: 

.. math:: :label: main5b

        \begin{bmatrix} 
        	M_{usr11} & M_{usr12} \\
                M_{usr21} & I 
        \end{bmatrix} 
        \begin{bmatrix} 
        	\ddot{U_{TP}} \\ 
                \ddot{q_m} 
        \end{bmatrix} +
        \begin{bmatrix} 
	         	0 & 0 \\
	                0 & C_{usr22} 
        \end{bmatrix}
         \begin{bmatrix} 
	        	\dot{U_{TP}} \\ 
	                \dot{q_m} 
        \end{bmatrix} +
        \begin{bmatrix} K_{usr11} & 0 \\
			0          & K_{usr22}
        \end{bmatrix} 
        \begin{bmatrix} 
        	U_{TP} \\ 
                q_m 
        \end{bmatrix} =
        \begin{bmatrix} F_{TP} + F_{usr11} \\
                        F_{usr21}  
        \end{bmatrix}  


:math:`F_{usr11}` is a (*6*\ ×\ *1*) vector, and :math:`F_{usr21}` is a (*m*\ ×\ *1*) vector.

:math:`q_m` are a reduced set of generalized modal DOFs, which are chosen as the first few (*m*) eigenvectors that
are arranged by increasing eigenfrequencies, and that :math:`{\Phi_m}` (*L*\ ×\ *m* matrix) represents the retained, mass-normalized, internal eigenmodes, 
and  :math:`{\Omega_m}` is the diagonal (*m*\ ×\ *m*) matrix
containing the corresponding eigenfrequencies. As done in SubDyn, the user
decides how many modes to retain, including possibly zero or all modes.
Retaining zero modes corresponds to a Guyan (static) reduction;
retaining all modes corresponds to keeping the full finite-element
model.

Also to underline is that in SubDyn, the only damping matrix term retained is the one
associated with internal DOF damping. This assumption has implications
on the damping at the interface with the turbine system, as discussed in
Section :ref:`_TowerTurbineCpling`. The diagonal (*m*\ ×\ *m*) :math:`\zeta` matrix contains the modal
damping ratios corresponding to each retained internal mode. In SubDyn,
the user provides damping ratios (in percent of critical damping
coefficients) for the retained modes. Perhaps this can be changed here, as :math:`{C_{usr}}` is provided by the user.
                        
During initialization, UsrPtfmLoad reads in :math:`{M_{usr}}`, :math:`{C_{usr}}`, and :math:`{K_{usr}}` and the :math:`{F_{usr}}`
The substructure response at each time step can then be obtained by using
the state-space formulation discussed in the next section.


.. _UPSS:

State-Space Formulation    
~~~~~~~~~~~~~~~~~~~~~~~~~~

A state-space formulation of the substructure structural dynamics
problem was devised to integrate UsrPtfmLoad within the FAST modularization
framework. The state-space formulation was developed in terms of inputs,
outputs, states, and parameters. The notations highlighted here are
consistent with those used in :cite:{jonkman2013}. Inputs (identified by *u*)
are a set of values supplied to UsrPtfmLoad that, along with the states, are
needed to calculate future states and the system’s output. Outputs (*y*)
are a set of values calculated by and returned from UsrPtfmLoad that depend
on the states, inputs, and/or parameters through output equations (with
functions *Y*). States are a set of internal values of UsrPtfmLoad that are
influenced by the inputs and used to calculate future state values and
the output. As was done in SubDyn, only continuous states are employed in UsrPtfmLoad. Continuous
states (*x*) are states that are differentiable in time and
characterized by continuous time differential equations (with functions
*X*). Parameters (*p*) are a set of internal system values that are
independent of the states and inputs. Furthermore, parameters can be
fully defined at initialization and characterize the system’s state
equations and output equations.

In UsrPtfmLoad, the inputs are defined as:

.. math:: :label: UPinputs

	u = \begin{bmatrix}
		u1 \\ 
		u2 \\
		u3 \\
	     \end{bmatrix} = \begin{bmatrix}
	     			U_{TP} \\
	     			\dot{U}_{TP}  \\
	     			\ddot{U}_{TP}  \\
	     		     \end{bmatrix}
	     			

where :math:`{ U_{TP},\dot{U}_{TP}, \quad \textrm{and} \quad \ddot{U}_{TP}}` are TP deflections (6 DOFs), velocities, and
accelerations, respectively. Note that, in contrast to what is done in SubDyn, :math:`F_L` (i.e., the hydrodynamic forces on every interior node of the
substructure from HydroDyn, and :math:`{F_{HDR}}` (the analogous forces at the boundary nodes) are not included, as they are effectively included in :math:`{F_{usr}}`.

In first-order form, the states are defined as:

.. math:: :label: UPstates

	x = \begin{bmatrix}
		x1 \\ 
		x2 \\
 	     \end{bmatrix} = \begin{bmatrix}
	     			q_m  \\
	     			\dot{q}_m  \\
	     		     \end{bmatrix}
	     		     
	     		     
From the system equation of motion, the state equation corresponding to
Eq. :eq:`main4b` can be written as a standard linear system state equation:

.. math:: :label: UPstate_eq

	\dot{x} = X = A_{UP} x +B_{UP} u + F_{UP}

where

.. math:: :label: UP_ABFx

	A_{UP} = \begin{bmatrix}
		0 & I \\ 
		-\Omega_m^2 & -2 \zeta \Omega_m
            \end{bmatrix} = \begin{bmatrix}
				0 & I \\ 
				-K_{usr22} & -C_{usr22}
            \end{bmatrix}
            
	B_{UP} = \begin{bmatrix}
		0 & 0  & 0  \\ 
		0 & 0  & -\tilde{M}_{mB} 
            \end{bmatrix} = \begin{bmatrix}
		           0 & 0  & 0  \\ 
		           0 & 0  & -M_{usr21} 
            \end{bmatrix}
            
	F_{UP} = \begin{bmatrix}
		0 \\ 
		\Phi_m^T (F_{L}+F_{Lg}) 
            \end{bmatrix} = \begin{bmatrix}
				0 \\ 
				F_{usr21} 
                            \end{bmatrix}


Note that the dimensions of the matrix partitions in :eq:`UP_ABFx` are:

.. math:: :label: UP_ABFxdims

	-K_{usr22} , -C_{usr22} \rightarrow (m \times m) \quad \textrm{bottom right partitions of} \quad K_{usr} , C_{usr}
	
        -M_{usr21}              \rightarrow (m \times 6) \quad \textrm{bottom left partition of} \quad M_{usr}

	F_{usr21}               \rightarrow (m \times 1) \quad \textrm{bottom partition of} \quad F_{usr}
                            


In UsrPtfmLoad, the outputs to the ElastoDyn module are the reaction forces at the transition piece :math:`-F_{TP}`:

.. math:: :label: UPy1

	y1 = Y_1 =-F_{TP}

By examining Eq. :eq:`main5b` , the output equation can be found after substituting for :math:`{\ddot{q}_{m}}` as:

.. math:: :label: UPY1
	
	 \ddot{q}_{m} = F_{usr21} - \begin{bmatrix}  
	                             \Omega_m^2 &  2\zeta \Omega_m  
	                            \end{bmatrix} x - \begin{bmatrix}
	                                                    0 & 0 &  M_{usr21} 
	                                               \end{bmatrix} u
	 
	 Y_1 =C_{UP} x + D_{UP} \bar{u} + F_{UP}


where

.. math:: :label: CupDupFup

	C_{UP} = \begin{bmatrix}  M_{usr12} K_{usr22} & + M_{usr12} C_{usr22} \end{bmatrix}
	
	D_{UP} = \begin{bmatrix}  -K_{usr11} & 0 & M_{usr12}M_{usr21} - M_{usr11}  \end{bmatrix}
	
	F_{UP} = \begin{bmatrix}  F_{usr11} -M_{usr12} F_{usr21} \end{bmatrix}
	

Note that the dimensions of the matrix partitions in :eq:`CupDupFup` are:

.. math:: :label: CupDupFupdims

	K_{usr11}              \rightarrow (6 \times 6) \quad \textrm{top left partition of} \quad K_{usr} 
	
	K_{usr22}              \rightarrow (m \times m) \quad \textrm{bottom right partition of} \quad K_{usr} 
	
	C_{usr22}              \rightarrow (m \times m) \quad \textrm{bottom right partition of} \quad C_{usr} 
        
        M_{usr12}              \rightarrow (6 \times m) \quad \textrm{top right partition of} \quad M_{usr}
        
        M_{usr21}              \rightarrow (6 \times m) \quad \textrm{bottom left partition of} \quad M_{usr}

	F_{usr11}               \rightarrow (6 \times 1) \quad \textrm{top partition of} \quad F_{usr}
	
	F_{usr21}               \rightarrow (m \times 1) \quad \textrm{bottom partition of} \quad F_{usr}
                            



.. _UPTimeIntegration:

Time Integration  
~~~~~~~~~~~~~~~~~

At time :math:`{t=0}`, the initial states are specified as initial conditions (all
assumed to be zero in UsrPtfmLoad) and the initial inputs are supplied to
UsrPtfmLoad. During each subsequent time step, the inputs and states are
known values, with the inputs :math:`u(t)` coming from ElastoDyn, and
the states :math:`x(t)` known from the previous time-step integration. All of the
parameter matrices are read-in and calculated in the UsrPtfmLoad initiation module. With
known :math:`u(t)` and :math:`x(t)`, :math:`{\dot{x}(t)}` can be calculated using the state equation :math:`{\dot{x}(t)=X(u,x,t)}` (see Eq. :eq:`UPstate_eq`), and
the outputs :math:`y_1(t)` can be calculated solving Eq. :eq:`UPY1`.  The next time-step states :math:`{x(t + \Delta t)}` are
obtained by integration:

.. math:: :label: UPintegration

	\left [ u(t), \dot{x}(t), x(t) \right ] \xrightarrow[]{\text{Integrate}}  x(t + \Delta t)
	
	
For loose coupling, UsrPtfmLoad uses its own integrator, whereas for tight
coupling, the states from all the modules will be integrated
simultaneously using an integrator in the glue-code. UsrPtfmLoad’s built-in
time integrator options for loose coupling are:

-  Fourth-order explicit Runge-Kutta

-  Fourth-order explicit Adams-Bashforth predictor

-  Fourth-order explicit Adams-Bashforth-Moulton predictor-corrector

-  Implicit second-order Adams-Moulton.

For more information, consult any numerical methods reference, e.g.,
:cite:`chapra2010`.

.. _UPsim:

Static-Improvement Method
~~~~~~~~~~~~~~~~~~~~~~~~~
To account for the effects of static gravity (member self-weight) and
buoyancy forces, one would have to include all of the structural axial
modes in the C-B reduction. This inclusion often translates into
hundreds of modes to be retained for practical problems. An alternative
method is thus promoted to reduce this limitation and speed up SubDyn.
This method is denoted as SIM, and computes two static solutions at each
time step: one based on the full system stiffness matrix and one based
on the reduced stiffness matrix. The dynamic solution then proceeds as
described in the previous sections, and at each time step the
time-varying dynamic solution is superimposed on the difference between
the two static solutions, which amounts to quasi-statically accounting
for the contribution of those modes not directly included within the
dynamic solution.

Recalling the previous C-B formulation :eq:`CB3`, and adding the total static
deflection of all the internal DOFs (:math:`U_{L0}`), and subtracting the static
deflection associated with C-B modes (:math:`U_{L0m}`), the SIM formulation is cast as
in :eq:`SIM`:

.. math::   :label: SIM

	U_L = \hat{U}_L + U_{L0} - U_{L0m} = \underbrace{\Phi_R U_R + \Phi_m q_m}_{\hat{U}_L}  +  U_{L0} - U_{L0m} 
 
 
Eq. :eq:`SIM` can be rewritten as:

.. math::  :label: SIM2

        \begin{bmatrix} 
        	U_R \\ 
                U_L 
        \end{bmatrix} =
	  \begin{bmatrix} 
        	I & 0 & 0 & 0 \\
           \Phi_R & \Phi_m & \Phi_L & -\Phi_m 
        \end{bmatrix} 
        \begin{bmatrix} 
        	U_R \\ 
                q_m \\
                q_{L0} \\
                q_{m0}
        \end{bmatrix}

with:

.. math::  :label: UL0

	U_{L0} = \Phi_L q_{L0}
	
.. math::  :label: UL0m

	U_{L0m} = \Phi_m q_{m0}
	

where :math:`{q_{m0}}` and :math:`{q_{L0}}` are the *m* and *L* modal coefficients that are assumed to be
operating in a static fashion. For Eqs. :eq:`SIM2` and :eq:`UL0` to be valid, and are
calculated under the C-B hypothesis that the boundary nodes are fixed.

The static displacement vectors can also be calculated as follows:


.. math::  :label: SIM3
	
	K_{LL} U_{L0} = F_L + F_{Lg}

By making use of :eq:`UL0`, and by pre-multiplying both sides times , Eq. :eq:`SIM3` can be
rewritten as: :math:`{\Phi_L^T K_{LL} \Phi_L q_{L0} = \Phi_L^T  \left( F_L + F_{Lg} \right) = \tilde{F}_L }` or, recalling that :math:`{\Phi_L^T K_{LL} \Phi_L = \Omega_L^2}`, as: :math:`{\Omega_L^2 q_{L0} =\tilde{F}_L }`, or equivalently in terms of :math:`U_{L0}`:

.. math::  :label: UL02

	U_{L0} = \Phi_L \left[ \Omega_L^2 \right]^{-1} \tilde{F}_L 

Similarly:

.. math::  :label: UL0m2

	U_{L0m} = \Phi_m \left[ \Omega_m^2 \right]^{-1} \tilde{F}_m 

Note that: :math:`{ \dot{U}_{L0} = \dot{q}_{L0} = \dot{U}_{L0m} = \dot{q}_{m0} =0 }` and :math:`{ \ddot{U}_{L0} = \ddot{q}_{L0} = \ddot{U}_{L0m} = \ddot{q}_{m0} =0 }`.

The dynamic component :math:`{ \hat{U} = \begin{bmatrix} \hat{U}_R \\ \hat{U}_R \end{bmatrix} }` is calculated following the usual procedure
described in Sections :ref:`_UPSS` -- :ref:`_UPTimeIntegration`. For example, states are still
calculated and integrated as in Eq. :eq:`UPstate_eq`, and the output to ElastoDyn, i.e.,
the reaction provided by the substructure at the TP interface, is also
calculated as it was done previously in Eqs. :eq:`y1` and :eq:`Y1`.

However, the state-space formulation is slightly modified to allow for
the calculation of the outputs to HydroDyn as:

.. math:: :label: y2sim

	y_2= = \begin{bmatrix}
        	\bar{U}_R \\
                     U_L  \\
           	\dot{\bar{U}}_R  \\
           	\dot{U}_L \\
           	\ddot{\bar{U}}_R  \\
           	\ddot{U}_L \\
	     \end{bmatrix} = \begin{bmatrix}  
	     	\bar{U}_R \\
	     	\hat{U}_L + U_{L0} - U_{L0m} \\
	     	\dot{\bar{U}}_R  \\
		\dot{U}_L \\
		\ddot{\bar{U}}_R  \\
           	\ddot{U}_L \\
	     \end{bmatrix}

.. math:: :label: Y2sim

  Y_2 = C_2 x + D_2 u + F_{Y2}

where the matrices now have the following meaning:

.. math:: :label:  C2D2FY2sim

	C_2 = \begin{bmatrix}
	       0 & 0 \\
	       \Phi_m & 0 \\
	       0 & 0 \\
	       0 & \Phi_m \\
	       0 & 0 \\
	   -\Phi_m \Omega_m^2 & -2 \Phi_m \zeta \Omega_m \\
	      \end{bmatrix}
	      
	D_2 = \begin{bmatrix}
	       T_I & 0 & 0 & 0 & 0 \\
	       \bar{\Phi}_R T_I & 0 & 0 & 0 & 0 \\
	       0 & T_I  & 0 & 0 & 0 \\
	       0 & \bar{\Phi}_R T_I & 0 & 0 & 0 \\
	       0 & 0 & T_I  & 0 & 0  \\
	       0 & 0 & \bar{\Phi}_R T_I - \Phi_m \tilde{M}_{mB} &  \Phi_m \Phi_m^T & 0 
	      \end{bmatrix}

	F_{Y2} = \begin{bmatrix}
	       0 \\
	       U_{L0} - U_{L0m} \\
	       0 \\
	       0 \\
	       0 \\
	       \Phi_m \Phi_m^T F_{Lg} 
	      \end{bmatrix}


Finally, the element forces can be calculated as:

.. math:: :label: el_loads_sim

	\text{Element Inertia load:} ~~ F_I^e = [m] \ddot{U}_e 
	
	\text{Element Static load:} ~~ F_S^e = [k] U_e = [k] \left[ \hat{U}_e + U_{L0,e} - U_{L0m,e} \right] 
	
with the element node DOFs expressed as:

.. math::  :label: Uesim

	U_e = \hat{U}_e + U_{L0,e} - U_{L0m,e}

where the SIM decomposition is still used with :math:`\hat{U}_e` denoting the
time-varying components of the elements nodes’ displacements, and :math:`U_{L0,e}` and :math:`U_{L0m,e}` are
derived from the parent :math:`U_{L0}` and :math:`U_{L0m}` arrays of displacements, respectively.








