--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+-----------+
|    Region Info    |   Core 0  |
+-------------------+-----------+
| RDTSC Runtime [s] | 12.819910 |
|     call count    |         2 |
+-------------------+-----------+
+------------------------------------------+---------+-------------+
|                   Event                  | Counter |    Core 0   |
+------------------------------------------+---------+-------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 40325030000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 19384130000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 26526550000 |
|              PWR_PKG_ENERGY              |   PWR0  |    401.1357 |
|              PWR_DRAM_ENERGY             |   PWR3  |     51.1964 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      761098 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 11401930000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      402100 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     8148557 |
|               CAS_COUNT_RD               | MBOX0C0 |     6608816 |
|               CAS_COUNT_WR               | MBOX0C1 |    10888620 |
|               CAS_COUNT_RD               | MBOX1C0 |     6600158 |
|               CAS_COUNT_WR               | MBOX1C1 |    10887390 |
|               CAS_COUNT_RD               | MBOX2C0 |     6604833 |
|               CAS_COUNT_WR               | MBOX2C1 |    10888900 |
|               CAS_COUNT_RD               | MBOX3C0 |     5844602 |
|               CAS_COUNT_WR               | MBOX3C1 |    10661730 |
|               CAS_COUNT_RD               | MBOX4C0 |     5851240 |
|               CAS_COUNT_WR               | MBOX4C1 |    10660920 |
|               CAS_COUNT_RD               | MBOX5C0 |     5847198 |
|               CAS_COUNT_WR               | MBOX5C1 |    10657240 |
+------------------------------------------+---------+-------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |   12.8199 |
|        Runtime unhalted [s]       |    7.4734 |
|            Clock [MHz]            | 1895.3762 |
|                CPI                |    0.4807 |
|             Energy [J]            |  401.1357 |
|             Power [W]             |   31.2901 |
|          Energy DRAM [J]          |   51.1964 |
|           Power DRAM [W]          |    3.9935 |
|            DP [MFLOP/s]           |  894.7215 |
|          AVX DP [MFLOP/s]         |    5.2104 |
|          Packed [MUOPS/s]         |    0.7264 |
|          Scalar [MUOPS/s]         |  889.3924 |
|  Memory read bandwidth [MBytes/s] |  186.4941 |
|  Memory read data volume [GBytes] |    2.3908 |
| Memory write bandwidth [MBytes/s] |  322.7220 |
| Memory write data volume [GBytes] |    4.1373 |
|    Memory bandwidth [MBytes/s]    |  509.2162 |
|    Memory data volume [GBytes]    |    6.5281 |
|       Operational intensity       |    1.7571 |
+-----------------------------------+-----------+
--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Solve, Group 1: MEM_DP
+-------------------+-------------+
|    Region Info    |    Core 0   |
+-------------------+-------------+
| RDTSC Runtime [s] | 2044.964000 |
|     call count    |           1 |
+-------------------+-------------+
+------------------------------------------+---------+---------------+
|                   Event                  | Counter |     Core 0    |
+------------------------------------------+---------+---------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 8454334000000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 3807578000000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 5210383000000 |
|              PWR_PKG_ENERGY              |   PWR0  |    64403.1500 |
|              PWR_DRAM_ENERGY             |   PWR3  |     8056.4770 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        714138 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 2518901000000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       1138324 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |      10452010 |
|               CAS_COUNT_RD               | MBOX0C0 |    2237158000 |
|               CAS_COUNT_WR               | MBOX0C1 |     973524800 |
|               CAS_COUNT_RD               | MBOX1C0 |    2236350000 |
|               CAS_COUNT_WR               | MBOX1C1 |     973510700 |
|               CAS_COUNT_RD               | MBOX2C0 |    2236891000 |
|               CAS_COUNT_WR               | MBOX2C1 |     973822300 |
|               CAS_COUNT_RD               | MBOX3C0 |    2210220000 |
|               CAS_COUNT_WR               | MBOX3C1 |     972488600 |
|               CAS_COUNT_RD               | MBOX4C0 |    2211140000 |
|               CAS_COUNT_WR               | MBOX4C1 |     972300100 |
|               CAS_COUNT_RD               | MBOX5C0 |    2209971000 |
|               CAS_COUNT_WR               | MBOX5C1 |     971560900 |
+------------------------------------------+---------+---------------+
+-----------------------------------+------------+
|               Metric              |   Core 0   |
+-----------------------------------+------------+
|        Runtime (RDTSC) [s]        |  2044.9640 |
|        Runtime unhalted [s]       |  1467.9759 |
|            Clock [MHz]            |  1895.4357 |
|                CPI                |     0.4504 |
|             Energy [J]            | 64403.1500 |
|             Power [W]             |    31.4935 |
|          Energy DRAM [J]          |  8056.4770 |
|           Power DRAM [W]          |     3.9397 |
|            DP [MFLOP/s]           |  1231.8019 |
|          AVX DP [MFLOP/s]         |     0.0431 |
|          Packed [MUOPS/s]         |     0.0060 |
|          Scalar [MUOPS/s]         |  1231.7581 |
|  Memory read bandwidth [MBytes/s] |   417.5480 |
|  Memory read data volume [GBytes] |   853.8707 |
| Memory write bandwidth [MBytes/s] |   182.6835 |
| Memory write data volume [GBytes] |   373.5813 |
|    Memory bandwidth [MBytes/s]    |   600.2316 |
|    Memory data volume [GBytes]    |  1227.4520 |
|       Operational intensity       |     2.0522 |
+-----------------------------------+------------+
