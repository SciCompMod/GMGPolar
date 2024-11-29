--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+
|    Region Info    |  Core 0  |
+-------------------+----------+
| RDTSC Runtime [s] | 0.868439 |
|     call count    |        2 |
+-------------------+----------+
+------------------------------------------+---------+------------+
|                   Event                  | Counter |   Core 0   |
+------------------------------------------+---------+------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 2293625000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 1108274000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 1516729000 |
|              PWR_PKG_ENERGY              |   PWR0  |    27.4759 |
|              PWR_DRAM_ENERGY             |   PWR3  |     3.8530 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      23000 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  546170800 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      15090 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     141851 |
|               CAS_COUNT_RD               | MBOX0C0 |    1124935 |
|               CAS_COUNT_WR               | MBOX0C1 |    1223141 |
|               CAS_COUNT_RD               | MBOX1C0 |    1125863 |
|               CAS_COUNT_WR               | MBOX1C1 |    1224573 |
|               CAS_COUNT_RD               | MBOX2C0 |    1125711 |
|               CAS_COUNT_WR               | MBOX2C1 |    1223495 |
|               CAS_COUNT_RD               | MBOX3C0 |    1033513 |
|               CAS_COUNT_WR               | MBOX3C1 |    1173982 |
|               CAS_COUNT_RD               | MBOX4C0 |    1033731 |
|               CAS_COUNT_WR               | MBOX4C1 |    1173919 |
|               CAS_COUNT_RD               | MBOX5C0 |    1033310 |
|               CAS_COUNT_WR               | MBOX5C1 |    1173299 |
+------------------------------------------+---------+------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |    0.8684 |
|        Runtime unhalted [s]       |    0.4273 |
|            Clock [MHz]            | 1895.2334 |
|                CPI                |    0.4832 |
|             Energy [J]            |   27.4759 |
|             Power [W]             |   31.6382 |
|          Energy DRAM [J]          |    3.8530 |
|           Power DRAM [W]          |    4.4367 |
|            DP [MFLOP/s]           |  630.3401 |
|          AVX DP [MFLOP/s]         |    1.3762 |
|          Packed [MUOPS/s]         |    0.2072 |
|          Scalar [MUOPS/s]         |  628.9110 |
|  Memory read bandwidth [MBytes/s] |  477.3301 |
|  Memory read data volume [GBytes] |    0.4145 |
| Memory write bandwidth [MBytes/s] |  530.0478 |
| Memory write data volume [GBytes] |    0.4603 |
|    Memory bandwidth [MBytes/s]    | 1007.3778 |
|    Memory data volume [GBytes]    |    0.8748 |
|       Operational intensity       |    0.6257 |
+-----------------------------------+-----------+
--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Solve, Group 1: MEM_DP
+-------------------+-----------+
|    Region Info    |   Core 0  |
+-------------------+-----------+
| RDTSC Runtime [s] | 50.175300 |
|     call count    |         1 |
+-------------------+-----------+
+------------------------------------------+---------+--------------+
|                   Event                  | Counter |    Core 0    |
+------------------------------------------+---------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 238321600000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  91380540000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 125047700000 |
|              PWR_PKG_ENERGY              |   PWR0  |    1596.5970 |
|              PWR_DRAM_ENERGY             |   PWR3  |     240.0646 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        65560 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  63917010000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |        83556 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       417868 |
|               CAS_COUNT_RD               | MBOX0C0 |    161758300 |
|               CAS_COUNT_WR               | MBOX0C1 |     35933340 |
|               CAS_COUNT_RD               | MBOX1C0 |    161730500 |
|               CAS_COUNT_WR               | MBOX1C1 |     35898040 |
|               CAS_COUNT_RD               | MBOX2C0 |    161691500 |
|               CAS_COUNT_WR               | MBOX2C1 |     35959670 |
|               CAS_COUNT_RD               | MBOX3C0 |    160728700 |
|               CAS_COUNT_WR               | MBOX3C1 |     35872620 |
|               CAS_COUNT_RD               | MBOX4C0 |    160764700 |
|               CAS_COUNT_WR               | MBOX4C1 |     35837650 |
|               CAS_COUNT_RD               | MBOX5C0 |    160705400 |
|               CAS_COUNT_WR               | MBOX5C1 |     35781980 |
+------------------------------------------+---------+--------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |   50.1753 |
|        Runtime unhalted [s]       |   35.2314 |
|            Clock [MHz]            | 1895.4030 |
|                CPI                |    0.3834 |
|             Energy [J]            | 1596.5970 |
|             Power [W]             |   31.8204 |
|          Energy DRAM [J]          |  240.0646 |
|           Power DRAM [W]          |    4.7845 |
|            DP [MFLOP/s]           | 1273.9499 |
|          AVX DP [MFLOP/s]         |    0.0733 |
|          Packed [MUOPS/s]         |    0.0113 |
|          Scalar [MUOPS/s]         | 1273.8740 |
|  Memory read bandwidth [MBytes/s] | 1233.9191 |
|  Memory read data volume [GBytes] |   61.9123 |
| Memory write bandwidth [MBytes/s] |  274.5999 |
| Memory write data volume [GBytes] |   13.7781 |
|    Memory bandwidth [MBytes/s]    | 1508.5190 |
|    Memory data volume [GBytes]    |   75.6904 |
|       Operational intensity       |    0.8445 |
+-----------------------------------+-----------+
