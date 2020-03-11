% Modeling Osmosis and Desalination 

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 1.1
% ---------------------

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;

% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 35; % g -- the amount of solute in left side of tank
RGramSolute = .0001; % g

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
membraneSurfaceArea = 100; % cm^2


% Anonymous Functions
% --------------------------

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP)((dO-dP)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial volumes
LVolume(1) = 1000; % cm^3
RVolume(1) = 1000; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,0);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;

    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
 
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end

hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 1.1: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 1.2
% ---------------------

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;

% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 35; % g -- the amount of solute in left side of tank
RGramSolute = 20; % g

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
membraneSurfaceArea = 100; % cm^2


% Anonymous Functions
% --------------------------

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP)((dO-dP)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial volumes
LVolume(1) = 1000; % cm^3
RVolume(1) = 1000; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,0);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;

    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
 
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end

hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 1.2: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 1.3
% ---------------------

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;

% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 20; % g -- the amount of solute in left side of tank
RGramSolute = 35; % g

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
membraneSurfaceArea = 100; % cm^2


% Anonymous Functions
% --------------------------

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP)((dO-dP)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial volumes
LVolume(1) = 1000; % cm^3
RVolume(1) = 1000; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,0);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;

    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
 
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end

hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 1.3: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 1.4
% ---------------------

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;

% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 0; % g -- the amount of solute in left side of tank
RGramSolute = 35; % g

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
membraneSurfaceArea = 100; % cm^2


% Anonymous Functions
% --------------------------

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP)((dO-dP)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial volumes
LVolume(1) = 1000; % cm^3
RVolume(1) = 1000; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,0);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;

    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
 
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end

hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 1.4: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 2.1
% ---------------------

% Final Project: Modeling Osmosis and Desalination 

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LDepth = 1:1:numIterations; % Arrays for water depth in each side of tank
RDepth = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LHydrostaticPressure = 1:1:numIterations; % Arrays for hydrostatic pressure
RHydrostaticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;




% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 35; % g -- the amount of solute in left side of tank
RGramSolute = 0; % g

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl
constLiquidDens = .001; % kg/m -- liquid density of water helps compute 
                       % hydrostatic pressure in each side of tank
constGravity = 9.8; % m/s^2

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
tankLength = 20; % cm
tankWidth = 10; % cm
membraneSurfaceArea = 100; % cm^2
halfBase = .5 * tankLength * tankWidth; % cm^2

% Anonymous Functions
% --------------------------
% Computing the Hydrostatic Pressure according to water depth
HydrostaticPressure = @(Depth)(constLiquidDens*constGravity*Depth);

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP)((dO-dP)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial depth of water
LDepth(1) = 10; % cm
RDepth(1) = 10; % cm

% Initial volumes
LVolume(1) = LDepth(1)*halfBase; % cm^3
RVolume(1) = RDepth(1)*halfBase; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 


LHydrostaticPressure(1) = HydrostaticPressure(LDepth(1));
RHydrostaticPressure(1) = HydrostaticPressure(RDepth(1));

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations     
    deltaHP = RHydrostaticPressure(i) - LHydrostaticPressure(i); % Flows L->R Pos, R->L Neg
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,deltaHP);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;
     
    LDepth(i+1) = LVolume(i)/halfBase;
    RDepth(i+1) = RVolume(i)/halfBase;
    
    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
    
    LHydrostaticPressure(i+1) = HydrostaticPressure(LDepth(i+1));
    RHydrostaticPressure(i+1) = HydrostaticPressure(RDepth(i+1));
     
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end
 
hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 2.1: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 2.2
% ---------------------

% Final Project: Modeling Osmosis and Desalination 

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LDepth = 1:1:numIterations; % Arrays for water depth in each side of tank
RDepth = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LHydrostaticPressure = 1:1:numIterations; % Arrays for hydrostatic pressure
RHydrostaticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;




% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 35; % g -- the amount of solute in left side of tank
RGramSolute = 20; % g

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl
constLiquidDens = .001; % kg/m -- liquid density of water helps compute 
                       % hydrostatic pressure in each side of tank
constGravity = 9.8; % m/s^2

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
tankLength = 20; % cm
tankWidth = 10; % cm
membraneSurfaceArea = 100; % cm^2
halfBase = .5 * tankLength * tankWidth; % cm^2

% Anonymous Functions
% --------------------------
% Computing the Hydrostatic Pressure according to water depth
HydrostaticPressure = @(Depth)(constLiquidDens*constGravity*Depth);

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP)((dO-dP)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial depth of water
LDepth(1) = 10; % cm
RDepth(1) = 10; % cm

% Initial volumes
LVolume(1) = LDepth(1)*halfBase; % cm^3
RVolume(1) = RDepth(1)*halfBase; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 


LHydrostaticPressure(1) = HydrostaticPressure(LDepth(1));
RHydrostaticPressure(1) = HydrostaticPressure(RDepth(1));

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations     
    deltaHP = RHydrostaticPressure(i) - LHydrostaticPressure(i); % Flows L->R Pos, R->L Neg
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,deltaHP);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;
     
    LDepth(i+1) = LVolume(i)/halfBase;
    RDepth(i+1) = RVolume(i)/halfBase;
    
    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
    
    LHydrostaticPressure(i+1) = HydrostaticPressure(LDepth(i+1));
    RHydrostaticPressure(i+1) = HydrostaticPressure(RDepth(i+1));
     
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end
 
hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 2.2: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 2.3
% ---------------------

% Final Project: Modeling Osmosis and Desalination 

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LDepth = 1:1:numIterations; % Arrays for water depth in each side of tank
RDepth = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LHydrostaticPressure = 1:1:numIterations; % Arrays for hydrostatic pressure
RHydrostaticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;




% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 20; % g -- the amount of solute in left side of tank
RGramSolute = 35; % g

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl
constLiquidDens = .001; % kg/m -- liquid density of water helps compute 
                       % hydrostatic pressure in each side of tank
constGravity = 9.8; % m/s^2

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
tankLength = 20; % cm
tankWidth = 10; % cm
membraneSurfaceArea = 100; % cm^2
halfBase = .5 * tankLength * tankWidth; % cm^2

% Anonymous Functions
% --------------------------
% Computing the Hydrostatic Pressure according to water depth
HydrostaticPressure = @(Depth)(constLiquidDens*constGravity*Depth);

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP)((dO-dP)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial depth of water
LDepth(1) = 10; % cm
RDepth(1) = 10; % cm

% Initial volumes
LVolume(1) = LDepth(1)*halfBase; % cm^3
RVolume(1) = RDepth(1)*halfBase; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 


LHydrostaticPressure(1) = HydrostaticPressure(LDepth(1));
RHydrostaticPressure(1) = HydrostaticPressure(RDepth(1));

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations     
    deltaHP = RHydrostaticPressure(i) - LHydrostaticPressure(i); % Flows L->R Pos, R->L Neg
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,deltaHP);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;
     
    LDepth(i+1) = LVolume(i)/halfBase;
    RDepth(i+1) = RVolume(i)/halfBase;
    
    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
    
    LHydrostaticPressure(i+1) = HydrostaticPressure(LDepth(i+1));
    RHydrostaticPressure(i+1) = HydrostaticPressure(RDepth(i+1));
     
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end
 
hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 2.3: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 2.4
% ---------------------

% Final Project: Modeling Osmosis and Desalination 

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LDepth = 1:1:numIterations; % Arrays for water depth in each side of tank
RDepth = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LHydrostaticPressure = 1:1:numIterations; % Arrays for hydrostatic pressure
RHydrostaticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;




% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 0; % g -- the amount of solute in left side of tank
RGramSolute = 35; % g

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl
constLiquidDens = .001; % kg/m -- liquid density of water helps compute 
                       % hydrostatic pressure in each side of tank
constGravity = 9.8; % m/s^2

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
tankLength = 20; % cm
tankWidth = 10; % cm
membraneSurfaceArea = 100; % cm^2
halfBase = .5 * tankLength * tankWidth; % cm^2

% Anonymous Functions
% --------------------------
% Computing the Hydrostatic Pressure according to water depth
HydrostaticPressure = @(Depth)(constLiquidDens*constGravity*Depth);

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP)((dO-dP)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial depth of water
LDepth(1) = 10; % cm
RDepth(1) = 10; % cm

% Initial volumes
LVolume(1) = LDepth(1)*halfBase; % cm^3
RVolume(1) = RDepth(1)*halfBase; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 


LHydrostaticPressure(1) = HydrostaticPressure(LDepth(1));
RHydrostaticPressure(1) = HydrostaticPressure(RDepth(1));

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations     
    deltaHP = RHydrostaticPressure(i) - LHydrostaticPressure(i); % Flows L->R Pos, R->L Neg
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,deltaHP);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;
     
    LDepth(i+1) = LVolume(i)/halfBase;
    RDepth(i+1) = RVolume(i)/halfBase;
    
    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
    
    LHydrostaticPressure(i+1) = HydrostaticPressure(LDepth(i+1));
    RHydrostaticPressure(i+1) = HydrostaticPressure(RDepth(i+1));
     
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end
 
hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 2.4: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 3.1
% ---------------------

% Final Project: Modeling Osmosis and Desalination 

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LDepth = 1:1:numIterations; % Arrays for water depth in each side of tank
RDepth = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LHydrostaticPressure = 1:1:numIterations; % Arrays for hydrostatic pressure
RHydrostaticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;




% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 35; % g -- the amount of solute in left side of tank
RGramSolute = 0; % g
Piston = 500;

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl
constLiquidDens = .001; % kg/m -- liquid density of water helps compute 
                       % hydrostatic pressure in each side of tank
constGravity = 9.8; % m/s^2

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
tankLength = 20; % cm
tankWidth = 10; % cm
membraneSurfaceArea = 100; % cm^2
halfBase = .5 * tankLength * tankWidth; % cm^2

% Anonymous Functions
% --------------------------
% Computing the Hydrostatic Pressure according to water depth
HydrostaticPressure = @(Depth)(constLiquidDens*constGravity*Depth);

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP,Piston)((dO-dP-Piston)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial depth of water
LDepth(1) = 10; % cm
RDepth(1) = 10; % cm

% Initial volumes
LVolume(1) = LDepth(1)*halfBase; % cm^3
RVolume(1) = RDepth(1)*halfBase; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 


LHydrostaticPressure(1) = HydrostaticPressure(LDepth(1));
RHydrostaticPressure(1) = HydrostaticPressure(RDepth(1));

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations     
    deltaHP = RHydrostaticPressure(i) - LHydrostaticPressure(i); % Flows L->R Pos, R->L Neg
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,deltaHP,Piston);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;
     
    LDepth(i+1) = LVolume(i)/halfBase;
    RDepth(i+1) = RVolume(i)/halfBase;
    
    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
    
    LHydrostaticPressure(i+1) = HydrostaticPressure(LDepth(i+1));
    RHydrostaticPressure(i+1) = HydrostaticPressure(RDepth(i+1));
     
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end
 
hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 3.1: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;

%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION 3.2
% ---------------------

% Final Project: Modeling Osmosis and Desalination 

% Initialization
% --------------------------

% Simulation Variables
simLength = 100; % Length of simulation in seconds
deltaTime = .0001; % Timestep
numIterations = simLength/deltaTime; % Number of iterations to complete
                                     % simulation
times = 0:deltaTime:simLength; % An array to keep track of time                                     

LVolume = 1:1:numIterations; % Arrays for water in each side of tank
RVolume = 1:1:numIterations;

LDepth = 1:1:numIterations; % Arrays for water depth in each side of tank
RDepth = 1:1:numIterations;

LOsmoticPressure = 1:1:numIterations; % Arrays for osmotic pressure
ROsmoticPressure = 1:1:numIterations;

LHydrostaticPressure = 1:1:numIterations; % Arrays for hydrostatic pressure
RHydrostaticPressure = 1:1:numIterations;

LMolarity = 1:1:numIterations; % Arrays for molar concentration
RMolarity = 1:1:numIterations;




% Variables
% -------------------------

hydraulicCond = .0001; % m/s -- a property of membranes that can be used
                       % to measure permeability. This value represents
                       % semi-permeable
temperature = 290; % K -- the average surface temperature of seawater
LGramSolute = 35; % g -- the amount of solute in left side of tank
RGramSolute = 0; % g
Piston = 1000;

% Constants
% --------------------------

molarWeight = 58; % g/mol -- molar weight of NaCl

constOsmoticCoefficient = .93;

constIdealGas = 8.314; % J/(K*mol)
constVantHoffFactor = 2; % -- a scalar that helps estimation of osmotic
                         % pressure and changes based on chemical breakdown
                         % of solute. This is set to 2 for NaCl
constLiquidDens = .001; % kg/m -- liquid density of water helps compute 
                       % hydrostatic pressure in each side of tank
constGravity = 9.8; % m/s^2

LMolsSolute = LGramSolute/molarWeight;
RMolsSolute = RGramSolute/molarWeight;

% Tank and Membrane Dimensions
% --------------------------
tankLength = 20; % cm
tankWidth = 10; % cm
membraneSurfaceArea = 100; % cm^2
halfBase = .5 * tankLength * tankWidth; % cm^2

% Anonymous Functions
% --------------------------
% Computing the Hydrostatic Pressure according to water depth
HydrostaticPressure = @(Depth)(constLiquidDens*constGravity*Depth);

% Computing the Osmotic Pressure according to molar concentration of solute
OsmoticPressure = @(Molarity)(Molarity*constIdealGas*temperature*constVantHoffFactor*constOsmoticCoefficient);

% Computing flux according to net force in system
flux = @(dO,dP,Piston)((dO-dP-Piston)*hydraulicCond);

% Computing change in volume
dVolumedt = @(J)(J*membraneSurfaceArea*.0001);

% Initial Conditions
% --------------------------

% Initial depth of water
LDepth(1) = 10; % cm
RDepth(1) = 10; % cm

% Initial volumes
LVolume(1) = LDepth(1)*halfBase; % cm^3
RVolume(1) = RDepth(1)*halfBase; % cm^3

% Initial Concentrations
LMolarity(1) = (LMolsSolute/LVolume(1))*1000; % mol/L
RMolarity(1) = (RMolsSolute/RVolume(1))*1000; 


LHydrostaticPressure(1) = HydrostaticPressure(LDepth(1));
RHydrostaticPressure(1) = HydrostaticPressure(RDepth(1));

LOsmoticPressure(1) = OsmoticPressure(LMolarity(1));
ROsmoticPressure(1) = OsmoticPressure(RMolarity(1));

for i=1:numIterations     
    deltaHP = RHydrostaticPressure(i) - LHydrostaticPressure(i); % Flows L->R Pos, R->L Neg
    deltaOP = LOsmoticPressure(i) - ROsmoticPressure(i);
     
    Flux = flux(deltaOP,deltaHP,Piston);
     
    dVdt = dVolumedt(Flux)*1000;
     
    LVolume(i+1) = LVolume(i) + dVdt;
    RVolume(i+1) = RVolume(i) - dVdt;
     
    LDepth(i+1) = LVolume(i)/halfBase;
    RDepth(i+1) = RVolume(i)/halfBase;
    
    LMolarity(i+1) = LMolsSolute/LVolume(i+1);
    RMolarity(i+1) = RMolsSolute/RVolume(i+1);
    
    LHydrostaticPressure(i+1) = HydrostaticPressure(LDepth(i+1));
    RHydrostaticPressure(i+1) = HydrostaticPressure(RDepth(i+1));
     
    LOsmoticPressure(i+1) = OsmoticPressure(LMolarity(i+1));
    ROsmoticPressure(i+1) = OsmoticPressure(RMolarity(i+1));
end
 
hold on;
plot(times,LVolume);
plot(times,RVolume);
title("Sim 3.2: Volume in Tanks");
xlabel("Time (s)");
ylabel("Volume (cubic centimeters)");
legend("Left Side Volume", "Right Side Volume");
hold off;
figure;