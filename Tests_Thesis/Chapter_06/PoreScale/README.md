# PoreScale

Thesis section: 6.3.2.

Represents pore-scale fluid-structure interaction using a hexahedral void-space grid and independently discretized solid grains. A finite-volume flow solve supplies pressure that is transferred to grain surfaces through mortar interpolation, followed by a mechanical solve with grain-contact constraints. The two stages form a one-way flow-to-mechanics coupling.

