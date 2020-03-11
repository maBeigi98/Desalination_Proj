# Desalination Model

## Introduction

The goal of this project is ultimately to present a model of desalination based on the forces of osmosis and hydrostatic pressure.

One common approach to desalination of seawater is filtration across a semipermeable membrane just porous enough to permit water molecules to pass through but not other ions like NaCl.

Osmosis is a chemical process by which a solution separated by semipermeable membrane breaks its intuitive equilibrium in such a way that one side of the membrane experiences a swell in volume from the other. The new equilibrium of the tank occurs with equal solute concentrations. My first set of simulations demonstrate this process from even water levels and disparate concentrations to even concentration with disparate water levels. In the context of desalination, osmosis is relevant because it pulls freshwater on one side of the membrane towards the salt water on the other side, which effectively undoes the filtration. Osmotic pressure is the primary force in this system working against desalination.

Hydrostatic pressure is the second force that acts against desalination. This pressure represents the gravitational pressure on water in a tank and since the goal of our model is to move an entire quantity of water across a membrane, the disparate heights of the columns of water result in a strong pressure to return the water back across the filter to towards an even height. My second set of simulations demonstrate how this force naturally works against and limits osmosis until the force from the hydrostatic pressure which works to even the heights of the water columns and the force of the osmotic pressure which works to even the solute concentrations of the water columns balance each other out.

The final step of my project introduces a new force of a piston which introduces a third pressure which counteracts both and ultimately filters the water.

## Variables

| Variable in Code | Explanation |
| --- | --- |
| dVdt | change in volume |
| Flux | movement of water |
| hydraulicCond | set to .0001 for a semi-permeable membrane |
| membraneSurfaceArea | 10cm |
| LDepth, RDepth | depth of water computed from volume, length, and width  |
| LVolume, RVolume | initially 1000cm |
| tankLength | set to 20cm |
| tankWidth | set to 10cm |
| LOsmoticPressure, ROsmoticPressure | osmotic pressure |
| constVantHoffFactor | set to 2 for ionic composition of NaCl |
| constOsmoticCoefficient | set to .93 for NaCl |
| LMolarity, RMolarity | molarity of solute |
| LMolsSolute, RMolsSolute | mols of solute in each side |
| LGramSolute, RGramSolute | initial grams of solute in each side |
| molarWeight | molar weight of solute |
| constIdealGas | 8.314 – ideal gas constant |
| temperature | set to 290K average surface temperature of sea water |
| LHydrostaticPressure, RHydrostaticPressure |  hydrostatic pressure |
| constLiquidDens | 997 – the density of water |
| constGravity | 9.87 – force of gravity |

## Simulations—All simulations are done over 100s with a deltaTime of .0001

1. Osmosis Model
   1. Simulation 1: 35 grams of NaCl in Left Tank, 0 in Right Tank
   2. Simulation 2: 35 grams of NaCl in Left Tank, 20 in Right Tank
   3. Simulation 3: 20 grams of NaCl in Left Tank, 35 in Right Tank
   4. Simulation 4: 0 grams of NaCl in Left Tank, 35 in Right Tank

Each of these simulations are produced by the script and demonstrate as expected the movement of the volume of liquid from the Tank with the greater concentration to the tank with the lesser concentration. Furthermore, where both tanks have some amount of salt the osmosis happens at a less extreme pace over time.

2. Hydrostatic Pressure Model
   1. Simulation 1: 35 grams of NaCl in Left Tank, 0 in Right Tank
   2. Simulation 2: 35 grams of NaCl in Left Tank, 20 in Right Tank
   3. Simulation 3: 20 grams of NaCl in Left Tank, 35 in Right Tank
   4. Simulation 4: 0 grams of NaCl in Left Tank, 35 in Right Tank

These four simulations used the exact same conditions except also involved the hydrostatic pressure shifting with the changing depths and slowing the pull of osmosis.

3. Desalination
   1. Simulation 1: 35 grams of NaCl in Left Tank, 0 in Right Tank, 500 kPa in Piston
   2. Simulation 2: 35 grams of NaCl in Left Tank, 0 in Right Tank, 1000 kPa

This section dabbled into the desalination and introduced a piston to help drive the movement of the water.

## Assumptions and Limitations

The biggest problem with my model proved to be the fact that the volume is driven into the negatives. There are a few possible fixes to this but I do not they resolve the root issue of the surface area of the membrane not being proportional to the depths of the water. However, I could find no resource to guide me on my calculations there.

My model works on some basic assumptions about osmosis and diffusion. I am assuming that the solute is always perfectly distilled in the solvent. In practice, this cannot be the case and the pressures pulling and pushing across the membrane would realistically affect the dissolution. With more time and experience with the subject material, I could have included a model of how salt interacts, builds up, and may clog the membrane over time. I have neglected other smaller details related to diffusion and osmosis like the potential for the system having more salt than can be dissolved in the water, the flow of water from the higher column back across the membrane, and neglected the reality that seawater in fact contains in abundance a few more ions than just NaCl.

The model does not explore the more specific physics of filtration. For example, I am also assuming that the friction within the membrane and tank is nonexistent and is not impacting the temperature of the system which is absolute and controlled. Furthermore, The equations I am using are themselves models that are still being developed and refined. I made choices to use specific equations (e.g. some do not accept the vant Hoff factor) according to what seemed most achievable and appropriate for solutions composed of ionic solutes with a water solvent.

## Conclusions and Future Direction

My model serves as an effective basis for simple filtration and osmosis. It also can be used to test how a constant pressure perhaps applied from a piston may offset the chemical forces. In addition to addressing some of the oversights described in the Assumptions and Limitations section, some future endeavors may explore how drawing the water from different depths of seawater may impact the desalination process or how different membranes may speed up the process. While my simulations do not treat time as itself a variable worthy of optimization, this model can definitely serve as a basis for examining more closely how the speed of desalination is affected under various conditions.
