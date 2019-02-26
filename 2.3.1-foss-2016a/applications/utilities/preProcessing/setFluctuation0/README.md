# setFluctuation0

## source
https://bitbucket.org/lesituu/computewallquantities

## Functionality
add `cylindrical symmetrical` mean flow and fluctuations

## Line
mean flow is imposed in direction `vector n_ (0, 0, 0)` : no mean flow yet   
fluctuation enveloppes : `envelopeUzRMS` and `envelopeUrRMS`.   
For the fact that the randoms are generated initially in `[0,1]`   
=> `([0,1]-0.5)*2` => `[-1,1]` : uniform distribution [-1,1] (variable name `randomField`) has a   
variance = 2^2 * 1/12 = 1/3 meaning std = 1/sqrt(3)   
=> if we want a std equaling `envelopeUrRMS`, we need to do   
` (sqrt(3) * randomField) * envelopeUrRMS ` to make the std of `(sqrt(3) * randomField)` to be unity !

## Developpment Approacheseses

## Input

## Output
r   
U1

## How
setFluctuation0 -time '5'

## Limit

# Cope with limit

## Dig the mine of OpenFOAM
