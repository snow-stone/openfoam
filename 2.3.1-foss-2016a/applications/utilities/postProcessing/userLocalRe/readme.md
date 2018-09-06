# Utility
userLocalRe

## Version
2.3.1-foss-2016a

## Prototype
utility range

## Functionality
- 计算LocalRe: velocity * D / nu
- 并输出`string("LocalRe_")+string("mag")+velocityFieldName`

## How To Use
1. 保证时间步文件中有速度场和粘性场(命名为`nu`)
2. userLocalRe UMean -time '0.5' 或 userLocalRe U -time '0.5' 等，第一个argument是流场即可

## Limit
LocalRe有量纲，因为常量D有量纲

## Cope with limit

## Note
