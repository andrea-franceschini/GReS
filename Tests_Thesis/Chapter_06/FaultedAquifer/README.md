# Faulted aquifer

Thesis section: 6.3.3; mesh reference: Figure 6.10.

Models ten years of withdrawal from two wells in sandy layers separated from a stiff rock formation by a fault inclined 40 degrees from vertical. Flow is solved first; its pressure change then loads a mechanical model with stabilized P0 mortar contact, friction angle 30 degrees, and zero cohesion. Reservoir contraction induces fault slip against the comparatively rigid rock.

Run `runFaultedTestCase` after `initGReS`, from this folder. The driver calls `runFlowSimulation` followed by `runContactMechanicsSimulation`. Native GReS grids reproduce the 2,000 x 1,000 x 100 m setting; precise horizon undulations and fault position remain reconstruction parameters because their formulas are not specified in the thesis. Additional functions are in `Utils/`, with materials in `InputFlow/` and `InputMech/`.

`Outputs/` is overwritten on each run and contains flow and mechanical histories, pressure-transfer data, and mesh plots. Positive flow values represent drawdown; mechanics converts them to a negative pressure change in MPa. MATLAB numerical reproduction has not been verified in the preparation environment.
