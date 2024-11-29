--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+
|    Region Info    |  Core 0  |
+-------------------+----------+
| RDTSC Runtime [s] | 3.435907 |
|     call count    |        2 |
+-------------------+----------+
+------------------------------------------+---------+------------+
|                   Event                  | Counter |   Core 0   |
+------------------------------------------+---------+------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 9104843000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 4382045000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 5997069000 |
|              PWR_PKG_ENERGY              |   PWR0  |   108.4174 |
|              PWR_DRAM_ENERGY             |   PWR3  |    15.4438 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |     122920 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 2184098000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      81488 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |    1080242 |
|               CAS_COUNT_RD               | MBOX0C0 |    5028529 |
|               CAS_COUNT_WR               | MBOX0C1 |    5330393 |
|               CAS_COUNT_RD               | MBOX1C0 |    5027551 |
|               CAS_COUNT_WR               | MBOX1C1 |    5330151 |
|               CAS_COUNT_RD               | MBOX2C0 |    5030219 |
|               CAS_COUNT_WR               | MBOX2C1 |    5331873 |
|               CAS_COUNT_RD               | MBOX3C0 |    4704726 |
|               CAS_COUNT_WR               | MBOX3C1 |    5274794 |
|               CAS_COUNT_RD               | MBOX4C0 |    4706134 |
|               CAS_COUNT_WR               | MBOX4C1 |    5275915 |
|               CAS_COUNT_RD               | MBOX5C0 |    4702566 |
|               CAS_COUNT_WR               | MBOX5C1 |    5272807 |
+------------------------------------------+---------+------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |    3.4359 |
|        Runtime unhalted [s]       |    1.6895 |
|            Clock [MHz]            | 1895.2601 |
|                CPI                |    0.4813 |
|             Energy [J]            |  108.4174 |
|             Power [W]             |   31.5542 |
|          Energy DRAM [J]          |   15.4438 |
|           Power DRAM [W]          |    4.4948 |
|            DP [MFLOP/s]           |  638.3501 |
|          AVX DP [MFLOP/s]         |    2.6100 |
|          Packed [MUOPS/s]         |    0.3739 |
|          Scalar [MUOPS/s]         |  635.6685 |
|  Memory read bandwidth [MBytes/s] |  543.8978 |
|  Memory read data volume [GBytes] |    1.8688 |
| Memory write bandwidth [MBytes/s] |  592.6295 |
| Memory write data volume [GBytes] |    2.0362 |
|    Memory bandwidth [MBytes/s]    | 1136.5273 |
|    Memory data volume [GBytes]    |    3.9050 |
|       Operational intensity       |    0.5617 |
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
| RDTSC Runtime [s] | 192.119300 |
|     call count    |          1 |
+-------------------+------------+
+------------------------------------------+---------+--------------+
|                   Event                  | Counter |    Core 0    |
+------------------------------------------+---------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 905250500000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 352023700000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 481718600000 |
|              PWR_PKG_ENERGY              |   PWR0  |    6135.8010 |
|              PWR_DRAM_ENERGY             |   PWR3  |     923.2806 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |       218232 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 244225800000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       302358 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |      2178750 |
|               CAS_COUNT_RD               | MBOX0C0 |    675091900 |
|               CAS_COUNT_WR               | MBOX0C1 |    156282300 |
|               CAS_COUNT_RD               | MBOX1C0 |    674837900 |
|               CAS_COUNT_WR               | MBOX1C1 |    156198300 |
|               CAS_COUNT_RD               | MBOX2C0 |    675018600 |
|               CAS_COUNT_WR               | MBOX2C1 |    156169200 |
|               CAS_COUNT_RD               | MBOX3C0 |    671036200 |
|               CAS_COUNT_WR               | MBOX3C1 |    155837500 |
|               CAS_COUNT_RD               | MBOX4C0 |    671340400 |
|               CAS_COUNT_WR               | MBOX4C1 |    155975100 |
|               CAS_COUNT_RD               | MBOX5C0 |    670925500 |
|               CAS_COUNT_WR               | MBOX5C1 |    155805700 |
+------------------------------------------+---------+--------------+
+-----------------------------------+-----------+
|               Metric              |   Core 0  |
+-----------------------------------+-----------+
|        Runtime (RDTSC) [s]        |  192.1193 |
|        Runtime unhalted [s]       |  135.7191 |
|            Clock [MHz]            | 1895.4377 |
|                CPI                |    0.3889 |
|             Energy [J]            | 6135.8010 |
|             Power [W]             |   31.9375 |
|          Energy DRAM [J]          |  923.2806 |
|           Power DRAM [W]          |    4.8058 |
|            DP [MFLOP/s]           | 1271.3188 |
|          AVX DP [MFLOP/s]         |    0.0970 |
|          Packed [MUOPS/s]         |    0.0141 |
|          Scalar [MUOPS/s]         | 1271.2195 |
|  Memory read bandwidth [MBytes/s] | 1345.2476 |
|  Memory read data volume [GBytes] |  258.4480 |
| Memory write bandwidth [MBytes/s] |  311.8956 |
| Memory write data volume [GBytes] |   59.9212 |
|    Memory bandwidth [MBytes/s]    | 1657.1432 |
|    Memory data volume [GBytes]    |  318.3692 |
|       Operational intensity       |    0.7672 |
+-----------------------------------+-----------+
