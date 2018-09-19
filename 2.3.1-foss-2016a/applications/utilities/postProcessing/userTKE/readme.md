# Utility
userTKE

## Version
2.3.1-foss-2016a

## Prototype
utility range

## Functionality
- 读取U_mean
- 由U_mean计算湍动能
- option : noWriteField 不写TKE 到时间步
- option : noWriteLog 不写TKE体积平均后的值到userDefinedLog

## How To Use
userTKE 0.5 -time '0.1:0.5', 第一个argument为U_mean所在时间步，可选择option有二:  
userTKE 0.5 -time '0.1:0.5' -noWriteField (不写TKE)

## Limit

## Cope with limit

## Note
