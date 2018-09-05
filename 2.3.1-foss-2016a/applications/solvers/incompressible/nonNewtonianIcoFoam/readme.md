# Solver
nonNewtonianIcoFoamPS_profiling1

## Version
2.3.1-foss-2016a

## Prototype
原版是这一版的prototype

## Functionality
- add profiling : 将runTime.write()的所有时间步输出到文件[支持并行运行solver，并输出到文件]
- add passive scalar : 被动标量T (扩散率DT)[求解T方程放在U修正后]

## How To Use
1. nonNewtonianIcoFoamPS_profiling1
2. mpirun -n 4 nonNewtonianIcoFoamPS_profiling1 -parallel

## Limit
- 对于profiling : 目前把write()执行完毕后的"当地时间"输出注释掉，需要时可去掉注释

## Cope with limit

## Note
