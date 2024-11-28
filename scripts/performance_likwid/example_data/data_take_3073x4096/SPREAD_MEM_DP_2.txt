--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+----------+
|    Region Info    |  Core 0  |  Core 14 |
+-------------------+----------+----------+
| RDTSC Runtime [s] | 2.977116 | 2.977079 |
|     call count    |        2 |        2 |
+-------------------+----------+----------+
+------------------------------------------+---------+------------+------------+
|                   Event                  | Counter |   Core 0   |   Core 14  |
+------------------------------------------+---------+------------+------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 4944445000 | 4195855000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 2891364000 | 2455573000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 3957772000 | 3360262000 |
|              PWR_PKG_ENERGY              |   PWR0  |    93.7183 |    81.5170 |
|              PWR_DRAM_ENERGY             |   PWR3  |    13.6040 |    11.8306 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |     122920 |          0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 1110075000 | 1074023000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      81488 |          0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |    1080242 |          0 |
|               CAS_COUNT_RD               | MBOX0C0 |    4286775 |    1103856 |
|               CAS_COUNT_WR               | MBOX0C1 |    6168885 |    1439594 |
|               CAS_COUNT_RD               | MBOX1C0 |    4285411 |    1116524 |
|               CAS_COUNT_WR               | MBOX1C1 |    6180062 |    1455031 |
|               CAS_COUNT_RD               | MBOX2C0 |    4251751 |    1103678 |
|               CAS_COUNT_WR               | MBOX2C1 |    6115684 |    1439969 |
|               CAS_COUNT_RD               | MBOX3C0 |    4037798 |    1005571 |
|               CAS_COUNT_WR               | MBOX3C1 |    6188546 |    1412693 |
|               CAS_COUNT_RD               | MBOX4C0 |    4423230 |    1003218 |
|               CAS_COUNT_WR               | MBOX4C1 |    6567766 |    1409143 |
|               CAS_COUNT_RD               | MBOX5C0 |    4025507 |    1003526 |
|               CAS_COUNT_WR               | MBOX5C1 |    6174429 |    1409948 |
+------------------------------------------+---------+------------+------------+
+-----------------------------------------------+---------+------------+------------+------------+--------------+
|                     Event                     | Counter |     Sum    |     Min    |     Max    |      Avg     |
+-----------------------------------------------+---------+------------+------------+------------+--------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 9140300000 | 4195855000 | 4944445000 |   4570150000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 5346937000 | 2455573000 | 2891364000 |   2673468500 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 7318034000 | 3360262000 | 3957772000 |   3659017000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |   175.2353 |    81.5170 |    93.7183 |      87.6176 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |    25.4346 |    11.8306 |    13.6040 |      12.7173 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |     122920 |          0 |     122920 |        61460 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  | 2184098000 | 1074023000 | 1110075000 |   1092049000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |      81488 |          0 |      81488 |        40744 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |    1080242 |          0 |    1080242 |       540121 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |    5390631 |    1103856 |    4286775 | 2.695316e+06 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |    7608479 |    1439594 |    6168885 | 3.804240e+06 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |    5401935 |    1116524 |    4285411 | 2.700968e+06 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |    7635093 |    1455031 |    6180062 | 3.817546e+06 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |    5355429 |    1103678 |    4251751 | 2.677714e+06 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |    7555653 |    1439969 |    6115684 | 3.777826e+06 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |    5043369 |    1005571 |    4037798 | 2.521684e+06 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |    7601239 |    1412693 |    6188546 | 3.800620e+06 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |    5426448 |    1003218 |    4423230 |      2713224 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |    7976909 |    1409143 |    6567766 | 3.988454e+06 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |    5029033 |    1003526 |    4025507 | 2.514516e+06 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |    7584377 |    1409948 |    6174429 | 3.792188e+06 |
+-----------------------------------------------+---------+------------+------------+------------+--------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |    2.9771 |    2.9771 |
|        Runtime unhalted [s]       |    1.1147 |    0.9467 |
|            Clock [MHz]            | 1894.8925 | 1895.4498 |
|                CPI                |    0.5848 |    0.5852 |
|             Energy [J]            |   93.7183 |   81.5170 |
|             Power [W]             |   31.4796 |   27.3815 |
|          Energy DRAM [J]          |   13.6040 |   11.8306 |
|           Power DRAM [W]          |    4.5695 |    3.9739 |
|            DP [MFLOP/s]           |  375.9641 |  360.7640 |
|          AVX DP [MFLOP/s]         |    3.0123 |         0 |
|          Packed [MUOPS/s]         |    0.4315 |         0 |
|          Scalar [MUOPS/s]         |  372.8692 |  360.7640 |
|  Memory read bandwidth [MBytes/s] |  544.1072 |  136.2167 |
|  Memory read data volume [GBytes] |    1.6199 |    0.4055 |
| Memory write bandwidth [MBytes/s] |  803.9001 |  184.1564 |
| Memory write data volume [GBytes] |    2.3933 |    0.5482 |
|    Memory bandwidth [MBytes/s]    | 1348.0073 |  320.3731 |
|    Memory data volume [GBytes]    |    4.0132 |    0.9538 |
|       Operational intensity       |    0.2789 |    1.1261 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |    5.9542 |    2.9771 |    2.9771 |    2.9771 |
|        Runtime unhalted [s] STAT       |    2.0614 |    0.9467 |    1.1147 |    1.0307 |
|            Clock [MHz] STAT            | 3790.3423 | 1894.8925 | 1895.4498 | 1895.1712 |
|                CPI STAT                |    1.1700 |    0.5848 |    0.5852 |    0.5850 |
|             Energy [J] STAT            |  175.2353 |   81.5170 |   93.7183 |   87.6176 |
|             Power [W] STAT             |   58.8611 |   27.3815 |   31.4796 |   29.4306 |
|          Energy DRAM [J] STAT          |   25.4346 |   11.8306 |   13.6040 |   12.7173 |
|           Power DRAM [W] STAT          |    8.5434 |    3.9739 |    4.5695 |    4.2717 |
|            DP [MFLOP/s] STAT           |  736.7281 |  360.7640 |  375.9641 |  368.3641 |
|          AVX DP [MFLOP/s] STAT         |    3.0123 |         0 |    3.0123 |    1.5062 |
|          Packed [MUOPS/s] STAT         |    0.4315 |         0 |    0.4315 |    0.2157 |
|          Scalar [MUOPS/s] STAT         |  733.6332 |  360.7640 |  372.8692 |  366.8166 |
|  Memory read bandwidth [MBytes/s] STAT |  680.3239 |  136.2167 |  544.1072 |  340.1620 |
|  Memory read data volume [GBytes] STAT |    2.0254 |    0.4055 |    1.6199 |    1.0127 |
| Memory write bandwidth [MBytes/s] STAT |  988.0565 |  184.1564 |  803.9001 |  494.0282 |
| Memory write data volume [GBytes] STAT |    2.9415 |    0.5482 |    2.3933 |    1.4708 |
|    Memory bandwidth [MBytes/s] STAT    | 1668.3804 |  320.3731 | 1348.0073 |  834.1902 |
|    Memory data volume [GBytes] STAT    |    4.9670 |    0.9538 |    4.0132 |    2.4835 |
|       Operational intensity STAT       |    1.4050 |    0.2789 |    1.1261 |    0.7025 |
+----------------------------------------+-----------+-----------+-----------+-----------+
--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Solve, Group 1: MEM_DP
+-------------------+-----------+-----------+
|    Region Info    |   Core 0  |  Core 14  |
+-------------------+-----------+-----------+
| RDTSC Runtime [s] | 99.145200 | 99.145190 |
|     call count    |         1 |         1 |
+-------------------+-----------+-----------+
+------------------------------------------+---------+--------------+--------------+
|                   Event                  | Counter |    Core 0    |    Core 14   |
+------------------------------------------+---------+--------------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 455020800000 | 450753700000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 179437500000 | 185200400000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 245549700000 | 253449200000 |
|              PWR_PKG_ENERGY              |   PWR0  |    3197.4800 |    2840.3330 |
|              PWR_DRAM_ENERGY             |   PWR3  |     467.7027 |     485.7474 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |       218232 |            0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 122060500000 | 122165300000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       302358 |            0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |      2178750 |            0 |
|               CAS_COUNT_RD               | MBOX0C0 |    316436800 |    327618200 |
|               CAS_COUNT_WR               | MBOX0C1 |     87128320 |     82200630 |
|               CAS_COUNT_RD               | MBOX1C0 |    316326600 |    327494800 |
|               CAS_COUNT_WR               | MBOX1C1 |     87225280 |     82267840 |
|               CAS_COUNT_RD               | MBOX2C0 |    316340200 |    327451800 |
|               CAS_COUNT_WR               | MBOX2C1 |     87165400 |     82083780 |
|               CAS_COUNT_RD               | MBOX3C0 |    315234300 |    325484700 |
|               CAS_COUNT_WR               | MBOX3C1 |     86399020 |     81645190 |
|               CAS_COUNT_RD               | MBOX4C0 |    316630100 |    325638500 |
|               CAS_COUNT_WR               | MBOX4C1 |     87622440 |     81765940 |
|               CAS_COUNT_RD               | MBOX5C0 |    315523300 |    325459400 |
|               CAS_COUNT_WR               | MBOX5C1 |     86653320 |     81738080 |
+------------------------------------------+---------+--------------+--------------+
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
|                     Event                     | Counter |      Sum     |      Min     |      Max     |      Avg     |
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 905774500000 | 450753700000 | 455020800000 | 452887250000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 364637900000 | 179437500000 | 185200400000 | 182318950000 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 498998900000 | 245549700000 | 253449200000 | 249499450000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    6037.8130 |    2840.3330 |    3197.4800 |    3018.9065 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     953.4501 |     467.7027 |     485.7474 |     476.7251 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |       218232 |            0 |       218232 |       109116 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  | 244225800000 | 122060500000 | 122165300000 | 122112900000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |       302358 |            0 |       302358 |       151179 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |      2178750 |            0 |      2178750 |      1089375 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |    644055000 |    316436800 |    327618200 |    322027500 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |    169328950 |     82200630 |     87128320 |     84664475 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |    643821400 |    316326600 |    327494800 |    321910700 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |    169493120 |     82267840 |     87225280 |     84746560 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |    643792000 |    316340200 |    327451800 |    321896000 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |    169249180 |     82083780 |     87165400 |     84624590 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |    640719000 |    315234300 |    325484700 |    320359500 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |    168044210 |     81645190 |     86399020 |     84022105 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |    642268600 |    316630100 |    325638500 |    321134300 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |    169388380 |     81765940 |     87622440 |     84694190 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |    640982700 |    315523300 |    325459400 |    320491350 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |    168391400 |     81738080 |     86653320 |     84195700 |
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |   99.1452 |   99.1452 |
|        Runtime unhalted [s]       |   69.1800 |   71.4018 |
|            Clock [MHz]            | 1895.4240 | 1895.3245 |
|                CPI                |    0.3944 |    0.4109 |
|             Energy [J]            | 3197.4800 | 2840.3330 |
|             Power [W]             |   32.2505 |   28.6482 |
|          Energy DRAM [J]          |  467.7027 |  485.7474 |
|           Power DRAM [W]          |    4.7174 |    4.8994 |
|            DP [MFLOP/s]           | 1231.3211 | 1232.1858 |
|          AVX DP [MFLOP/s]         |    0.1880 |         0 |
|          Packed [MUOPS/s]         |    0.0272 |         0 |
|          Scalar [MUOPS/s]         | 1231.1287 | 1232.1858 |
|  Memory read bandwidth [MBytes/s] | 1224.2191 | 1264.6648 |
|  Memory read data volume [GBytes] |  121.3754 |  125.3854 |
| Memory write bandwidth [MBytes/s] |  337.0854 |  317.4021 |
| Memory write data volume [GBytes] |   33.4204 |   31.4689 |
|    Memory bandwidth [MBytes/s]    | 1561.3045 | 1582.0669 |
|    Memory data volume [GBytes]    |  154.7958 |  156.8543 |
|       Operational intensity       |    0.7886 |    0.7788 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |  198.2904 |   99.1452 |   99.1452 |   99.1452 |
|        Runtime unhalted [s] STAT       |  140.5818 |   69.1800 |   71.4018 |   70.2909 |
|            Clock [MHz] STAT            | 3790.7485 | 1895.3245 | 1895.4240 | 1895.3742 |
|                CPI STAT                |    0.8053 |    0.3944 |    0.4109 |    0.4026 |
|             Energy [J] STAT            | 6037.8130 | 2840.3330 | 3197.4800 | 3018.9065 |
|             Power [W] STAT             |   60.8987 |   28.6482 |   32.2505 |   30.4494 |
|          Energy DRAM [J] STAT          |  953.4501 |  467.7027 |  485.7474 |  476.7251 |
|           Power DRAM [W] STAT          |    9.6168 |    4.7174 |    4.8994 |    4.8084 |
|            DP [MFLOP/s] STAT           | 2463.5069 | 1231.3211 | 1232.1858 | 1231.7534 |
|          AVX DP [MFLOP/s] STAT         |    0.1880 |         0 |    0.1880 |    0.0940 |
|          Packed [MUOPS/s] STAT         |    0.0272 |         0 |    0.0272 |    0.0136 |
|          Scalar [MUOPS/s] STAT         | 2463.3145 | 1231.1287 | 1232.1858 | 1231.6572 |
|  Memory read bandwidth [MBytes/s] STAT | 2488.8839 | 1224.2191 | 1264.6648 | 1244.4419 |
|  Memory read data volume [GBytes] STAT |  246.7608 |  121.3754 |  125.3854 |  123.3804 |
| Memory write bandwidth [MBytes/s] STAT |  654.4875 |  317.4021 |  337.0854 |  327.2437 |
| Memory write data volume [GBytes] STAT |   64.8893 |   31.4689 |   33.4204 |   32.4447 |
|    Memory bandwidth [MBytes/s] STAT    | 3143.3714 | 1561.3045 | 1582.0669 | 1571.6857 |
|    Memory data volume [GBytes] STAT    |  311.6501 |  154.7958 |  156.8543 |  155.8251 |
|       Operational intensity STAT       |    1.5674 |    0.7788 |    0.7886 |    0.7837 |
+----------------------------------------+-----------+-----------+-----------+-----------+
