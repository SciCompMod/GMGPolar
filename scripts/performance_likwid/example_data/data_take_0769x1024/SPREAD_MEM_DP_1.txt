--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+
|    Region Info    |  Core 0  |
+-------------------+----------+
| RDTSC Runtime [s] | 0.223871 |
|     call count    |        2 |
+-------------------+----------+
+------------------------------------------+---------+-----------+
|                   Event                  | Counter |   Core 0  |
+------------------------------------------+---------+-----------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 584409900 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 282742800 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 386937000 |
|              PWR_PKG_ENERGY              |   PWR0  |    7.2479 |
|              PWR_DRAM_ENERGY             |   PWR3  |    0.9509 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      4800 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 136642200 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      2800 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     11152 |
|               CAS_COUNT_RD               | MBOX0C0 |    211634 |
|               CAS_COUNT_WR               | MBOX0C1 |    250371 |
|               CAS_COUNT_RD               | MBOX1C0 |    213095 |
|               CAS_COUNT_WR               | MBOX1C1 |    251079 |
|               CAS_COUNT_RD               | MBOX2C0 |    212852 |
|               CAS_COUNT_WR               | MBOX2C1 |    251197 |
|               CAS_COUNT_RD               | MBOX3C0 |    186577 |
|               CAS_COUNT_WR               | MBOX3C1 |    228706 |
|               CAS_COUNT_RD               | MBOX4C0 |    186175 |
|               CAS_COUNT_WR               | MBOX4C1 |    228104 |
|               CAS_COUNT_RD               | MBOX5C0 |    186305 |
|               CAS_COUNT_WR               | MBOX5C1 |    227980 |
+------------------------------------------+---------+-----------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |    0.2239 |
|        Runtime unhalted [s]       |    0.1090 |
|            Clock [MHz]            | 1895.2887 |
|                CPI                |    0.4838 |
|             Energy [J]            |    7.2479 |
|             Power [W]             |   32.3754 |
|          Energy DRAM [J]          |    0.9509 |
|           Power DRAM [W]          |    4.2477 |
|            DP [MFLOP/s]           |  610.8517 |
|          AVX DP [MFLOP/s]         |    0.4485 |
|          Packed [MUOPS/s]         |    0.0838 |
|          Scalar [MUOPS/s]         |  610.3602 |
|  Memory read bandwidth [MBytes/s] |  342.0930 |
|  Memory read data volume [GBytes] |    0.0766 |
| Memory write bandwidth [MBytes/s] |  410.9322 |
| Memory write data volume [GBytes] |    0.0920 |
|    Memory bandwidth [MBytes/s]    |  753.0252 |
|    Memory data volume [GBytes]    |    0.1686 |
|       Operational intensity       |    0.8112 |
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
| RDTSC Runtime [s] | 13.625510 |
|     call count    |         1 |
+-------------------+-----------+
+------------------------------------------+---------+-------------+
|                   Event                  | Counter |    Core 0   |
+------------------------------------------+---------+-------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 63087340000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 24110480000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 32993420000 |
|              PWR_PKG_ENERGY              |   PWR0  |    430.9313 |
|              PWR_DRAM_ENERGY             |   PWR3  |     58.5331 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |       22816 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 16693060000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       26128 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       58742 |
|               CAS_COUNT_RD               | MBOX0C0 |    29167960 |
|               CAS_COUNT_WR               | MBOX0C1 |     4692538 |
|               CAS_COUNT_RD               | MBOX1C0 |    29215580 |
|               CAS_COUNT_WR               | MBOX1C1 |     4703591 |
|               CAS_COUNT_RD               | MBOX2C0 |    29233670 |
|               CAS_COUNT_WR               | MBOX2C1 |     4735096 |
|               CAS_COUNT_RD               | MBOX3C0 |    29050570 |
|               CAS_COUNT_WR               | MBOX3C1 |     4723825 |
|               CAS_COUNT_RD               | MBOX4C0 |    29046210 |
|               CAS_COUNT_WR               | MBOX4C1 |     4712213 |
|               CAS_COUNT_RD               | MBOX5C0 |    29060190 |
|               CAS_COUNT_WR               | MBOX5C1 |     4688654 |
+------------------------------------------+---------+-------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |   13.6255 |
|        Runtime unhalted [s]       |    9.2957 |
|            Clock [MHz]            | 1895.4075 |
|                CPI                |    0.3822 |
|             Energy [J]            |  430.9313 |
|             Power [W]             |   31.6268 |
|          Energy DRAM [J]          |   58.5331 |
|           Power DRAM [W]          |    4.2958 |
|            DP [MFLOP/s]           | 1225.1784 |
|          AVX DP [MFLOP/s]         |    0.0422 |
|          Packed [MUOPS/s]         |    0.0079 |
|          Scalar [MUOPS/s]         | 1225.1329 |
|  Memory read bandwidth [MBytes/s] |  820.9269 |
|  Memory read data volume [GBytes] |   11.1855 |
| Memory write bandwidth [MBytes/s] |  132.7201 |
| Memory write data volume [GBytes] |    1.8084 |
|    Memory bandwidth [MBytes/s]    |  953.6470 |
|    Memory data volume [GBytes]    |   12.9939 |
|       Operational intensity       |    1.2847 |
+-----------------------------------+-----------+
