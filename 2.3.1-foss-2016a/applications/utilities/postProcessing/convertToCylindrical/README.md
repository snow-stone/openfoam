# convertToCylindrical

## source

## Functionality
convert velocity field (x,y,z) to cylindrical components (r,theta,z)

## Input
U   
UMean

## Output
Ucyl
UMeanCyl

## How
convertToCylindrical -time '150.4' -fields '(U)'   
convertToCylindrical -time '150.4' -fields '(UMean)'

## limit
The model must be oriented with the x-y plan at the r-theta plane and the z-axis must be the center axis of rotation. If not check `Usage` of the code.

# Cope with limit

# Dig the mine of OpenFOAM
