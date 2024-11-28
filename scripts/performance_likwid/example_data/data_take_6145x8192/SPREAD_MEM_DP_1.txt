--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+-----------+
|    Region Info    |   Core 0  |
+-------------------+-----------+
| RDTSC Runtime [s] | 13.771290 |
|     call count    |         2 |
+-------------------+-----------+
+------------------------------------------+---------+-------------+
|                   Event                  | Counter |    Core 0   |
+------------------------------------------+---------+-------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 36331100000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 17484480000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 23928180000 |
|              PWR_PKG_ENERGY              |   PWR0  |    435.1104 |
|              PWR_DRAM_ENERGY             |   PWR3  |     63.4649 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      761098 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  8735507000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      402100 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     8148557 |
|               CAS_COUNT_RD               | MBOX0C0 |    21114160 |
|               CAS_COUNT_WR               | MBOX0C1 |    22220920 |
|               CAS_COUNT_RD               | MBOX1C0 |    21112350 |
|               CAS_COUNT_WR               | MBOX1C1 |    22223270 |
|               CAS_COUNT_RD               | MBOX2C0 |    21118530 |
|               CAS_COUNT_WR               | MBOX2C1 |    22227100 |
|               CAS_COUNT_RD               | MBOX3C0 |    19788240 |
|               CAS_COUNT_WR               | MBOX3C1 |    21813560 |
|               CAS_COUNT_RD               | MBOX4C0 |    19798310 |
|               CAS_COUNT_WR               | MBOX4C1 |    21813520 |
|               CAS_COUNT_RD               | MBOX5C0 |    19786820 |
|               CAS_COUNT_WR               | MBOX5C1 |    21809990 |
+------------------------------------------+---------+-------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |   13.7713 |
|        Runtime unhalted [s]       |    6.7410 |
|            Clock [MHz]            | 1895.2823 |
|                CPI                |    0.4813 |
|             Energy [J]            |  435.1104 |
|             Power [W]             |   31.5955 |
|          Energy DRAM [J]          |   63.4649 |
|           Power DRAM [W]          |    4.6085 |
|            DP [MFLOP/s]           |  639.2884 |
|          AVX DP [MFLOP/s]         |    4.8504 |
|          Packed [MUOPS/s]         |    0.6762 |
|          Scalar [MUOPS/s]         |  634.3274 |
|  Memory read bandwidth [MBytes/s] |  570.3154 |
|  Memory read data volume [GBytes] |    7.8540 |
| Memory write bandwidth [MBytes/s] |  613.9537 |
| Memory write data volume [GBytes] |    8.4549 |
|    Memory bandwidth [MBytes/s]    | 1184.2691 |
|    Memory data volume [GBytes]    |   16.3089 |
|       Operational intensity       |    0.5398 |
+-----------------------------------+-----------+
--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Solve, Group 1: MEM_DP
+-------------------+------------+
|    Region Info    |   Core 0   |
+-------------------+------------+
| RDTSC Runtime [s] | 756.475800 |
|     call count    |          1 |
+-------------------+------------+
+------------------------------------------+---------+---------------+
|                   Event                  | Counter |     Core 0    |
+------------------------------------------+---------+---------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 3526386000000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 1405611000000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 1923470000000 |
|              PWR_PKG_ENERGY              |   PWR0  |    24139.6600 |
|              PWR_DRAM_ENERGY             |   PWR3  |     3677.4780 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        714138 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  954016300000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       1138324 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |      10452010 |
|               CAS_COUNT_RD               | MBOX0C0 |    2711245000 |
|               CAS_COUNT_WR               | MBOX0C1 |     621698200 |
|               CAS_COUNT_RD               | MBOX1C0 |    2711730000 |
|               CAS_COUNT_WR               | MBOX1C1 |     621910600 |
|               CAS_COUNT_RD               | MBOX2C0 |    2711082000 |
|               CAS_COUNT_WR               | MBOX2C1 |     621775600 |
|               CAS_COUNT_RD               | MBOX3C0 |    2700666000 |
|               CAS_COUNT_WR               | MBOX3C1 |     620968800 |
|               CAS_COUNT_RD               | MBOX4C0 |    2700815000 |
|               CAS_COUNT_WR               | MBOX4C1 |     620967500 |
|               CAS_COUNT_RD               | MBOX5C0 |    2700904000 |
|               CAS_COUNT_WR               | MBOX5C1 |     620937500 |
+------------------------------------------+---------+---------------+
+-----------------------------------+------------+
|               Metric              |   Core 0   |
+-----------------------------------+------------+
|        Runtime (RDTSC) [s]        |   756.4758 |
|        Runtime unhalted [s]       |   541.9189 |
|            Clock [MHz]            |  1895.4423 |
|                CPI                |     0.3986 |
|             Energy [J]            | 24139.6600 |
|             Power [W]             |    31.9107 |
|          Energy DRAM [J]          |  3677.4780 |
|           Power DRAM [W]          |     4.8613 |
|            DP [MFLOP/s]           |  1261.2511 |
|          AVX DP [MFLOP/s]         |     0.1166 |
|          Packed [MUOPS/s]         |     0.0163 |
|          Scalar [MUOPS/s]         |  1261.1326 |
|  Memory read bandwidth [MBytes/s] |  1373.6491 |
|  Memory read data volume [GBytes] |  1039.1323 |
| Memory write bandwidth [MBytes/s] |   315.4212 |
| Memory write data volume [GBytes] |   238.6085 |
|    Memory bandwidth [MBytes/s]    |  1689.0703 |
|    Memory data volume [GBytes]    |  1277.7408 |
|       Operational intensity       |     0.7467 |
+-----------------------------------+------------+
