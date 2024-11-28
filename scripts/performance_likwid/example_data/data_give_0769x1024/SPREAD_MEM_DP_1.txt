--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+
|    Region Info    |  Core 0  |
+-------------------+----------+
| RDTSC Runtime [s] | 0.208327 |
|     call count    |        2 |
+-------------------+----------+
+------------------------------------------+---------+-----------+
|                   Event                  | Counter |   Core 0  |
+------------------------------------------+---------+-----------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 646879000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 313140400 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 428515500 |
|              PWR_PKG_ENERGY              |   PWR0  |    6.7346 |
|              PWR_DRAM_ENERGY             |   PWR3  |    0.8150 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      4800 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 178310800 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      2800 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     11152 |
|               CAS_COUNT_RD               | MBOX0C0 |     50979 |
|               CAS_COUNT_WR               | MBOX0C1 |    111230 |
|               CAS_COUNT_RD               | MBOX1C0 |     50777 |
|               CAS_COUNT_WR               | MBOX1C1 |    110292 |
|               CAS_COUNT_RD               | MBOX2C0 |     50076 |
|               CAS_COUNT_WR               | MBOX2C1 |    110059 |
|               CAS_COUNT_RD               | MBOX3C0 |     37817 |
|               CAS_COUNT_WR               | MBOX3C1 |    102235 |
|               CAS_COUNT_RD               | MBOX4C0 |     38226 |
|               CAS_COUNT_WR               | MBOX4C1 |    102884 |
|               CAS_COUNT_RD               | MBOX5C0 |     37961 |
|               CAS_COUNT_WR               | MBOX5C1 |    101932 |
+------------------------------------------+---------+-----------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |    0.2083 |
|        Runtime unhalted [s]       |    0.1207 |
|            Clock [MHz]            | 1895.4163 |
|                CPI                |    0.4841 |
|             Energy [J]            |    6.7346 |
|             Power [W]             |   32.3271 |
|          Energy DRAM [J]          |    0.8150 |
|           Power DRAM [W]          |    3.9123 |
|            DP [MFLOP/s]           |  856.4439 |
|          AVX DP [MFLOP/s]         |    0.4820 |
|          Packed [MUOPS/s]         |    0.0900 |
|          Scalar [MUOPS/s]         |  855.9158 |
|  Memory read bandwidth [MBytes/s] |   81.6671 |
|  Memory read data volume [GBytes] |    0.0170 |
| Memory write bandwidth [MBytes/s] |  196.1932 |
| Memory write data volume [GBytes] |    0.0409 |
|    Memory bandwidth [MBytes/s]    |  277.8603 |
|    Memory data volume [GBytes]    |    0.0579 |
|       Operational intensity       |    3.0823 |
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
| RDTSC Runtime [s] | 35.383850 |
|     call count    |         1 |
+-------------------+-----------+
+------------------------------------------+---------+--------------+
|                   Event                  | Counter |    Core 0    |
+------------------------------------------+---------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 149678300000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  66030420000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  |  90357870000 |
|              PWR_PKG_ENERGY              |   PWR0  |    1106.0830 |
|              PWR_DRAM_ENERGY             |   PWR3  |     131.3437 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        22816 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  44190880000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |        26128 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |        58742 |
|               CAS_COUNT_RD               | MBOX0C0 |     11649560 |
|               CAS_COUNT_WR               | MBOX0C1 |      3704309 |
|               CAS_COUNT_RD               | MBOX1C0 |     11611320 |
|               CAS_COUNT_WR               | MBOX1C1 |      3703609 |
|               CAS_COUNT_RD               | MBOX2C0 |     11570020 |
|               CAS_COUNT_WR               | MBOX2C1 |      3689521 |
|               CAS_COUNT_RD               | MBOX3C0 |     11484930 |
|               CAS_COUNT_WR               | MBOX3C1 |      3685958 |
|               CAS_COUNT_RD               | MBOX4C0 |     11570040 |
|               CAS_COUNT_WR               | MBOX4C1 |      3700925 |
|               CAS_COUNT_RD               | MBOX5C0 |     11490690 |
|               CAS_COUNT_WR               | MBOX5C1 |      3676280 |
+------------------------------------------+---------+--------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |   35.3839 |
|        Runtime unhalted [s]       |   25.4573 |
|            Clock [MHz]            | 1895.4404 |
|                CPI                |    0.4411 |
|             Energy [J]            | 1106.0830 |
|             Power [W]             |   31.2595 |
|          Energy DRAM [J]          |  131.3437 |
|           Power DRAM [W]          |    3.7120 |
|            DP [MFLOP/s]           | 1248.9172 |
|          AVX DP [MFLOP/s]         |    0.0162 |
|          Packed [MUOPS/s]         |    0.0030 |
|          Scalar [MUOPS/s]         | 1248.8997 |
|  Memory read bandwidth [MBytes/s] |  125.4838 |
|  Memory read data volume [GBytes] |    4.4401 |
| Memory write bandwidth [MBytes/s] |   40.0827 |
| Memory write data volume [GBytes] |    1.4183 |
|    Memory bandwidth [MBytes/s]    |  165.5664 |
|    Memory data volume [GBytes]    |    5.8584 |
|       Operational intensity       |    7.5433 |
+-----------------------------------+-----------+
