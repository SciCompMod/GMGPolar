--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+
|    Region Info    |  Core 0  |
+-------------------+----------+
| RDTSC Runtime [s] | 3.215063 |
|     call count    |        2 |
+-------------------+----------+
+------------------------------------------+---------+-------------+
|                   Event                  | Counter |    Core 0   |
+------------------------------------------+---------+-------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 10103740000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  4863090000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  |  6654335000 |
|              PWR_PKG_ENERGY              |   PWR0  |     99.4171 |
|              PWR_DRAM_ENERGY             |   PWR3  |     12.8004 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      122920 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  2850725000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       81488 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     1080242 |
|               CAS_COUNT_RD               | MBOX0C0 |     1428678 |
|               CAS_COUNT_WR               | MBOX0C1 |     2669220 |
|               CAS_COUNT_RD               | MBOX1C0 |     1429962 |
|               CAS_COUNT_WR               | MBOX1C1 |     2669035 |
|               CAS_COUNT_RD               | MBOX2C0 |     1430908 |
|               CAS_COUNT_WR               | MBOX2C1 |     2668359 |
|               CAS_COUNT_RD               | MBOX3C0 |     1271884 |
|               CAS_COUNT_WR               | MBOX3C1 |     2567472 |
|               CAS_COUNT_RD               | MBOX4C0 |     1273473 |
|               CAS_COUNT_WR               | MBOX4C1 |     2567783 |
|               CAS_COUNT_RD               | MBOX5C0 |     1269357 |
|               CAS_COUNT_WR               | MBOX5C1 |     2568106 |
+------------------------------------------+---------+-------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |    3.2151 |
|        Runtime unhalted [s]       |    1.8749 |
|            Clock [MHz]            | 1895.5601 |
|                CPI                |    0.4813 |
|             Energy [J]            |   99.4171 |
|             Power [W]             |   30.9223 |
|          Energy DRAM [J]          |   12.8004 |
|           Power DRAM [W]          |    3.9814 |
|            DP [MFLOP/s]           |  889.5436 |
|          AVX DP [MFLOP/s]         |    2.7893 |
|          Packed [MUOPS/s]         |    0.3996 |
|          Scalar [MUOPS/s]         |  886.6778 |
|  Memory read bandwidth [MBytes/s] |  161.3258 |
|  Memory read data volume [GBytes] |    0.5187 |
| Memory write bandwidth [MBytes/s] |  312.7274 |
| Memory write data volume [GBytes] |    1.0054 |
|    Memory bandwidth [MBytes/s]    |  474.0533 |
|    Memory data volume [GBytes]    |    1.5241 |
|       Operational intensity       |    1.8765 |
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
| RDTSC Runtime [s] | 520.092800 |
|     call count    |          1 |
+-------------------+------------+
+------------------------------------------+---------+---------------+
|                   Event                  | Counter |     Core 0    |
+------------------------------------------+---------+---------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 2167689000000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  969320900000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 1326443000000 |
|              PWR_PKG_ENERGY              |   PWR0  |    16199.7800 |
|              PWR_DRAM_ENERGY             |   PWR3  |     2029.8250 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        218232 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  645109700000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |        302358 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       2178750 |
|               CAS_COUNT_RD               | MBOX0C0 |     527319600 |
|               CAS_COUNT_WR               | MBOX0C1 |     241276900 |
|               CAS_COUNT_RD               | MBOX1C0 |     527615600 |
|               CAS_COUNT_WR               | MBOX1C1 |     241484700 |
|               CAS_COUNT_RD               | MBOX2C0 |     527729400 |
|               CAS_COUNT_WR               | MBOX2C1 |     241321200 |
|               CAS_COUNT_RD               | MBOX3C0 |     524868500 |
|               CAS_COUNT_WR               | MBOX3C1 |     240993200 |
|               CAS_COUNT_RD               | MBOX4C0 |     525016700 |
|               CAS_COUNT_WR               | MBOX4C1 |     240625800 |
|               CAS_COUNT_RD               | MBOX5C0 |     524542400 |
|               CAS_COUNT_WR               | MBOX5C1 |     241058700 |
+------------------------------------------+---------+---------------+
+-----------------------------------+------------+
|               Metric              |   Core 0   |
+-----------------------------------+------------+
|        Runtime (RDTSC) [s]        |   520.0928 |
|        Runtime unhalted [s]       |   373.7125 |
|            Clock [MHz]            |  1895.4350 |
|                CPI                |     0.4472 |
|             Energy [J]            | 16199.7800 |
|             Power [W]             |    31.1479 |
|          Energy DRAM [J]          |  2029.8250 |
|           Power DRAM [W]          |     3.9028 |
|            DP [MFLOP/s]           |  1240.4109 |
|          AVX DP [MFLOP/s]         |     0.0358 |
|          Packed [MUOPS/s]         |     0.0052 |
|          Scalar [MUOPS/s]         |  1240.3742 |
|  Memory read bandwidth [MBytes/s] |   388.4959 |
|  Memory read data volume [GBytes] |   202.0539 |
| Memory write bandwidth [MBytes/s] |   178.0311 |
| Memory write data volume [GBytes] |    92.5927 |
|    Memory bandwidth [MBytes/s]    |   566.5269 |
|    Memory data volume [GBytes]    |   294.6466 |
|       Operational intensity       |     2.1895 |
+-----------------------------------+------------+
