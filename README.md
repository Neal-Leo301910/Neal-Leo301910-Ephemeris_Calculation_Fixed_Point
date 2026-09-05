# Neal-Leo301910-Ephemeris_Calculation_Fixed_Point

用于太阳星历计算的 C99 定点数实现。所有公开数值均为 Q16.16：
角度单位为度，距离单位为天文单位（AU），时间为距 J2000.0 的天数。

```sh
make test
```

`ephemeris_days_since_j2000()` 将公历 UTC 时间转换为定点天数；
`ephemeris_sun_position()` 返回太阳的地心黄经、黄纬和日地距离。实现仅使用整数
运算（包括 CORDIC 正弦），不依赖浮点运行时。