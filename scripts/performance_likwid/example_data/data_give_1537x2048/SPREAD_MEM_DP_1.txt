--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+
|    Region Info    |  Core 0  |
+-------------------+----------+
| RDTSC Runtime [s] | 0.806196 |
|     call count    |        2 |
+-------------------+----------+
+------------------------------------------+---------+------------+
|                   Event                  | Counter |   Core 0   |
+------------------------------------------+---------+------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 2543534000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 1225751000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 1677433000 |
|              PWR_PKG_ENERGY              |   PWR0  |    25.3612 |
|              PWR_DRAM_ENERGY             |   PWR3  |     3.1573 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      23000 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  712833400 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      15090 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     141851 |
|               CAS_COUNT_RD               | MBOX0C0 |     281705 |
|               CAS_COUNT_WR               | MBOX0C1 |     562894 |
|               CAS_COUNT_RD               | MBOX1C0 |     282426 |
|               CAS_COUNT_WR               | MBOX1C1 |     563732 |
|               CAS_COUNT_RD               | MBOX2C0 |     281820 |
|               CAS_COUNT_WR               | MBOX2C1 |     562802 |
|               CAS_COUNT_RD               | MBOX3C0 |     234409 |
|               CAS_COUNT_WR               | MBOX3C1 |     551552 |
|               CAS_COUNT_RD               | MBOX4C0 |     234236 |
|               CAS_COUNT_WR               | MBOX4C1 |     552392 |
|               CAS_COUNT_RD               | MBOX5C0 |     234076 |
|               CAS_COUNT_WR               | MBOX5C1 |     551400 |
+------------------------------------------+---------+------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |    0.8062 |
|        Runtime unhalted [s]       |    0.4726 |
|            Clock [MHz]            | 1895.3472 |
|                CPI                |    0.4819 |
|             Energy [J]            |   25.3612 |
|             Power [W]             |   31.4579 |
|          Energy DRAM [J]          |    3.1573 |
|           Power DRAM [W]          |    3.9162 |
|            DP [MFLOP/s]           |  885.7330 |
|          AVX DP [MFLOP/s]         |    1.4825 |
|          Packed [MUOPS/s]         |    0.2232 |
|          Scalar [MUOPS/s]         |  884.1935 |
|  Memory read bandwidth [MBytes/s] |  122.9415 |
|  Memory read data volume [GBytes] |    0.0991 |
| Memory write bandwidth [MBytes/s] |  265.5252 |
| Memory write data volume [GBytes] |    0.2141 |
|    Memory bandwidth [MBytes/s]    |  388.4667 |
|    Memory data volume [GBytes]    |    0.3132 |
|       Operational intensity       |    2.2801 |
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
| RDTSC Runtime [s] | 134.726000 |
|     call count    |          1 |
+-------------------+------------+
+------------------------------------------+---------+--------------+
|                   Event                  | Counter |    Core 0    |
+------------------------------------------+---------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 569206200000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 252580000000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 345636600000 |
|              PWR_PKG_ENERGY              |   PWR0  |    4218.3210 |
|              PWR_DRAM_ENERGY             |   PWR3  |     526.1740 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        65560 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 168975100000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |        83556 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       417868 |
|               CAS_COUNT_RD               | MBOX0C0 |    108219700 |
|               CAS_COUNT_WR               | MBOX0C1 |     49033000 |
|               CAS_COUNT_RD               | MBOX1C0 |    108067400 |
|               CAS_COUNT_WR               | MBOX1C1 |     48914270 |
|               CAS_COUNT_RD               | MBOX2C0 |    108089900 |
|               CAS_COUNT_WR               | MBOX2C1 |     48994020 |
|               CAS_COUNT_RD               | MBOX3C0 |    107149300 |
|               CAS_COUNT_WR               | MBOX3C1 |     48842420 |
|               CAS_COUNT_RD               | MBOX4C0 |    107309500 |
|               CAS_COUNT_WR               | MBOX4C1 |     48938470 |
|               CAS_COUNT_RD               | MBOX5C0 |    107001100 |
|               CAS_COUNT_WR               | MBOX5C1 |     48721600 |
+------------------------------------------+---------+--------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |  134.7260 |
|        Runtime unhalted [s]       |   97.3794 |
|            Clock [MHz]            | 1895.4440 |
|                CPI                |    0.4437 |
|             Energy [J]            | 4218.3210 |
|             Power [W]             |   31.3104 |
|          Energy DRAM [J]          |  526.1740 |
|           Power DRAM [W]          |    3.9055 |
|            DP [MFLOP/s]           | 1254.2413 |
|          AVX DP [MFLOP/s]         |    0.0273 |
|          Packed [MUOPS/s]         |    0.0042 |
|          Scalar [MUOPS/s]         | 1254.2130 |
|  Memory read bandwidth [MBytes/s] |  306.7972 |
|  Memory read data volume [GBytes] |   41.3336 |
| Memory write bandwidth [MBytes/s] |  139.3970 |
| Memory write data volume [GBytes] |   18.7804 |
|    Memory bandwidth [MBytes/s]    |  446.1942 |
|    Memory data volume [GBytes]    |   60.1140 |
|       Operational intensity       |    2.8110 |
+-----------------------------------+-----------+
