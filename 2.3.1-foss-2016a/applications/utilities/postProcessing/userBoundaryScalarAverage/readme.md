# Utility
userBoundaryScalarReadWrite2Txt

## Version
2.3.1-foss-2016a

## Prototype
range

## Functionality
- 对boundaryField上的某一个patch, 对一个scalarField (volScalarField or surfaceScalarField) 求个max min mean
- 可选择输出为纯txt格式方便python后处理
- 求和的话，就能算出那个面上的通量(几乎等同于patchIntegrate phi)

## How To Use
userBoundaryScalarReadWrite2Txt phi movingWall -time '0.5'

## Limit

## Cope with limit

## Note
