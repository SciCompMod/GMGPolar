--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+----------+
|    Region Info    |  Core 0  |  Core 14 |
+-------------------+----------+----------+
| RDTSC Runtime [s] | 0.195575 | 0.195520 |
|     call count    |        2 |        2 |
+-------------------+----------+----------+
+------------------------------------------+---------+-----------+-----------+
|                   Event                  | Counter |   Core 0  |  Core 14  |
+------------------------------------------+---------+-----------+-----------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 325476600 | 268330000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 221000700 | 266572200 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 302459900 | 365060200 |
|              PWR_PKG_ENERGY              |   PWR0  |    6.3079 |    5.5280 |
|              PWR_DRAM_ENERGY             |   PWR3  |    0.8639 |    0.7761 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      4800 |         0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  69489420 |  67152790 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      2800 |         0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     11152 |         0 |
|               CAS_COUNT_RD               | MBOX0C0 |    241688 |     41040 |
|               CAS_COUNT_WR               | MBOX0C1 |    343129 |     60902 |
|               CAS_COUNT_RD               | MBOX1C0 |    209337 |     40871 |
|               CAS_COUNT_WR               | MBOX1C1 |    304239 |     61167 |
|               CAS_COUNT_RD               | MBOX2C0 |    208500 |     41336 |
|               CAS_COUNT_WR               | MBOX2C1 |    302358 |     61284 |
|               CAS_COUNT_RD               | MBOX3C0 |    191272 |     33146 |
|               CAS_COUNT_WR               | MBOX3C1 |    298624 |     58036 |
|               CAS_COUNT_RD               | MBOX4C0 |    189161 |     33480 |
|               CAS_COUNT_WR               | MBOX4C1 |    296786 |     57849 |
|               CAS_COUNT_RD               | MBOX5C0 |    185719 |     33001 |
|               CAS_COUNT_WR               | MBOX5C1 |    291154 |     57770 |
+------------------------------------------+---------+-----------+-----------+
+-----------------------------------------------+---------+-----------+-----------+-----------+-------------+
|                     Event                     | Counter |    Sum    |    Min    |    Max    |     Avg     |
+-----------------------------------------------+---------+-----------+-----------+-----------+-------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 593806600 | 268330000 | 325476600 |   296903300 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 487572900 | 221000700 | 266572200 |   243786450 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 667520100 | 302459900 | 365060200 |   333760050 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |   11.8359 |    5.5280 |    6.3079 |      5.9179 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |    1.6400 |    0.7761 |    0.8639 |      0.8200 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |      4800 |         0 |      4800 |        2400 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  | 136642210 |  67152790 |  69489420 |    68321105 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |      2800 |         0 |      2800 |        1400 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |     11152 |         0 |     11152 |        5576 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |    282728 |     41040 |    241688 |      141364 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |    404031 |     60902 |    343129 | 202015.5000 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |    250208 |     40871 |    209337 |      125104 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |    365406 |     61167 |    304239 |      182703 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |    249836 |     41336 |    208500 |      124918 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |    363642 |     61284 |    302358 |      181821 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |    224418 |     33146 |    191272 |      112209 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |    356660 |     58036 |    298624 |      178330 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |    222641 |     33480 |    189161 | 111320.5000 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |    354635 |     57849 |    296786 | 177317.5000 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |    218720 |     33001 |    185719 |      109360 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |    348924 |     57770 |    291154 |      174462 |
+-----------------------------------------------+---------+-----------+-----------+-----------+-------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |    0.1956 |    0.1955 |
|        Runtime unhalted [s]       |    0.0852 |    0.1028 |
|            Clock [MHz]            | 1895.2104 | 1894.0087 |
|                CPI                |    0.6790 |    0.9934 |
|             Energy [J]            |    6.3079 |    5.5280 |
|             Power [W]             |   32.2532 |   28.2731 |
|          Energy DRAM [J]          |    0.8639 |    0.7761 |
|           Power DRAM [W]          |    4.4173 |    3.9692 |
|            DP [MFLOP/s]           |  355.8712 |  343.4572 |
|          AVX DP [MFLOP/s]         |    0.5134 |         0 |
|          Packed [MUOPS/s]         |    0.0959 |         0 |
|          Scalar [MUOPS/s]         |  355.3087 |  343.4572 |
|  Memory read bandwidth [MBytes/s] |  401.0912 |   72.9538 |
|  Memory read data volume [GBytes] |    0.0784 |    0.0143 |
| Memory write bandwidth [MBytes/s] |  600.9085 |  116.8602 |
| Memory write data volume [GBytes] |    0.1175 |    0.0228 |
|    Memory bandwidth [MBytes/s]    | 1001.9997 |  189.8140 |
|    Memory data volume [GBytes]    |    0.1960 |    0.0371 |
|       Operational intensity       |    0.3552 |    1.8094 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |    0.3911 |    0.1955 |    0.1956 |    0.1956 |
|        Runtime unhalted [s] STAT       |    0.1880 |    0.0852 |    0.1028 |    0.0940 |
|            Clock [MHz] STAT            | 3789.2191 | 1894.0087 | 1895.2104 | 1894.6096 |
|                CPI STAT                |    1.6724 |    0.6790 |    0.9934 |    0.8362 |
|             Energy [J] STAT            |   11.8359 |    5.5280 |    6.3079 |    5.9179 |
|             Power [W] STAT             |   60.5263 |   28.2731 |   32.2532 |   30.2631 |
|          Energy DRAM [J] STAT          |    1.6400 |    0.7761 |    0.8639 |    0.8200 |
|           Power DRAM [W] STAT          |    8.3865 |    3.9692 |    4.4173 |    4.1932 |
|            DP [MFLOP/s] STAT           |  699.3284 |  343.4572 |  355.8712 |  349.6642 |
|          AVX DP [MFLOP/s] STAT         |    0.5134 |         0 |    0.5134 |    0.2567 |
|          Packed [MUOPS/s] STAT         |    0.0959 |         0 |    0.0959 |    0.0479 |
|          Scalar [MUOPS/s] STAT         |  698.7659 |  343.4572 |  355.3087 |  349.3829 |
|  Memory read bandwidth [MBytes/s] STAT |  474.0450 |   72.9538 |  401.0912 |  237.0225 |
|  Memory read data volume [GBytes] STAT |    0.0927 |    0.0143 |    0.0784 |    0.0464 |
| Memory write bandwidth [MBytes/s] STAT |  717.7687 |  116.8602 |  600.9085 |  358.8843 |
| Memory write data volume [GBytes] STAT |    0.1403 |    0.0228 |    0.1175 |    0.0701 |
|    Memory bandwidth [MBytes/s] STAT    | 1191.8137 |  189.8140 | 1001.9997 |  595.9068 |
|    Memory data volume [GBytes] STAT    |    0.2331 |    0.0371 |    0.1960 |    0.1166 |
|       Operational intensity STAT       |    2.1646 |    0.3552 |    1.8094 |    1.0823 |
+----------------------------------------+-----------+-----------+-----------+-----------+
--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Solve, Group 1: MEM_DP
+-------------------+----------+----------+
|    Region Info    |  Core 0  |  Core 14 |
+-------------------+----------+----------+
| RDTSC Runtime [s] | 7.583367 | 7.583348 |
|     call count    |        1 |        1 |
+-------------------+----------+----------+
+------------------------------------------+---------+-------------+-------------+
|                   Event                  | Counter |    Core 0   |   Core 14   |
+------------------------------------------+---------+-------------+-------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 32300260000 | 30919110000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 12633160000 | 14174210000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 17287860000 | 19396270000 |
|              PWR_PKG_ENERGY              |   PWR0  |    242.0807 |    214.2365 |
|              PWR_DRAM_ENERGY             |   PWR3  |     30.1427 |     31.4624 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |       22816 |           0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  8348417000 |  8344647000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       26128 |           0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       58742 |           0 |
|               CAS_COUNT_RD               | MBOX0C0 |     8660702 |    10467910 |
|               CAS_COUNT_WR               | MBOX0C1 |     2048390 |     1997874 |
|               CAS_COUNT_RD               | MBOX1C0 |     8654639 |    10492580 |
|               CAS_COUNT_WR               | MBOX1C1 |     1987660 |     1975748 |
|               CAS_COUNT_RD               | MBOX2C0 |     8642327 |    10518030 |
|               CAS_COUNT_WR               | MBOX2C1 |     1979452 |     2004093 |
|               CAS_COUNT_RD               | MBOX3C0 |     8627399 |    10437760 |
|               CAS_COUNT_WR               | MBOX3C1 |     2002679 |     2008597 |
|               CAS_COUNT_RD               | MBOX4C0 |     8629937 |    10375880 |
|               CAS_COUNT_WR               | MBOX4C1 |     2005012 |     1989597 |
|               CAS_COUNT_RD               | MBOX5C0 |     8632536 |    10438040 |
|               CAS_COUNT_WR               | MBOX5C1 |     1976647 |     1996232 |
+------------------------------------------+---------+-------------+-------------+
+-----------------------------------------------+---------+-------------+-------------+-------------+--------------+
|                     Event                     | Counter |     Sum     |     Min     |     Max     |      Avg     |
+-----------------------------------------------+---------+-------------+-------------+-------------+--------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 63219370000 | 30919110000 | 32300260000 |  31609685000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 26807370000 | 12633160000 | 14174210000 |  13403685000 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 36684130000 | 17287860000 | 19396270000 |  18342065000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    456.3172 |    214.2365 |    242.0807 |     228.1586 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     61.6051 |     30.1427 |     31.4624 |      30.8026 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |       22816 |           0 |       22816 |        11408 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  | 16693064000 |  8344647000 |  8348417000 |   8346532000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |       26128 |           0 |       26128 |        13064 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |       58742 |           0 |       58742 |        29371 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |    19128612 |     8660702 |    10467910 |      9564306 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |     4046264 |     1997874 |     2048390 |      2023132 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |    19147219 |     8654639 |    10492580 | 9.573610e+06 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |     3963408 |     1975748 |     1987660 |      1981704 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |    19160357 |     8642327 |    10518030 | 9.580178e+06 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |     3983545 |     1979452 |     2004093 | 1.991772e+06 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |    19065159 |     8627399 |    10437760 | 9.532580e+06 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |     4011276 |     2002679 |     2008597 |      2005638 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |    19005817 |     8629937 |    10375880 | 9.502908e+06 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |     3994609 |     1989597 |     2005012 | 1.997304e+06 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |    19070576 |     8632536 |    10438040 |      9535288 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |     3972879 |     1976647 |     1996232 | 1.986440e+06 |
+-----------------------------------------------+---------+-------------+-------------+-------------+--------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |    7.5834 |    7.5833 |
|        Runtime unhalted [s]       |    4.8706 |    5.4647 |
|            Clock [MHz]            | 1895.4065 | 1895.4496 |
|                CPI                |    0.3911 |    0.4584 |
|             Energy [J]            |  242.0807 |  214.2365 |
|             Power [W]             |   31.9226 |   28.2509 |
|          Energy DRAM [J]          |   30.1427 |   31.4624 |
|           Power DRAM [W]          |    3.9748 |    4.1489 |
|            DP [MFLOP/s]           | 1100.9670 | 1100.3909 |
|          AVX DP [MFLOP/s]         |    0.0758 |         0 |
|          Packed [MUOPS/s]         |    0.0142 |         0 |
|          Scalar [MUOPS/s]         | 1100.8853 | 1100.3909 |
|  Memory read bandwidth [MBytes/s] |  437.5685 |  529.4143 |
|  Memory read data volume [GBytes] |    3.3182 |    4.0147 |
| Memory write bandwidth [MBytes/s] |  101.2729 |  101.0394 |
| Memory write data volume [GBytes] |    0.7680 |    0.7662 |
|    Memory bandwidth [MBytes/s]    |  538.8414 |  630.4537 |
|    Memory data volume [GBytes]    |    4.0862 |    4.7809 |
|       Operational intensity       |    2.0432 |    1.7454 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |   15.1667 |    7.5833 |    7.5834 |    7.5834 |
|        Runtime unhalted [s] STAT       |   10.3353 |    4.8706 |    5.4647 |    5.1677 |
|            Clock [MHz] STAT            | 3790.8561 | 1895.4065 | 1895.4496 | 1895.4280 |
|                CPI STAT                |    0.8495 |    0.3911 |    0.4584 |    0.4247 |
|             Energy [J] STAT            |  456.3172 |  214.2365 |  242.0807 |  228.1586 |
|             Power [W] STAT             |   60.1735 |   28.2509 |   31.9226 |   30.0868 |
|          Energy DRAM [J] STAT          |   61.6051 |   30.1427 |   31.4624 |   30.8026 |
|           Power DRAM [W] STAT          |    8.1237 |    3.9748 |    4.1489 |    4.0618 |
|            DP [MFLOP/s] STAT           | 2201.3579 | 1100.3909 | 1100.9670 | 1100.6789 |
|          AVX DP [MFLOP/s] STAT         |    0.0758 |         0 |    0.0758 |    0.0379 |
|          Packed [MUOPS/s] STAT         |    0.0142 |         0 |    0.0142 |    0.0071 |
|          Scalar [MUOPS/s] STAT         | 2201.2762 | 1100.3909 | 1100.8853 | 1100.6381 |
|  Memory read bandwidth [MBytes/s] STAT |  966.9828 |  437.5685 |  529.4143 |  483.4914 |
|  Memory read data volume [GBytes] STAT |    7.3329 |    3.3182 |    4.0147 |    3.6665 |
| Memory write bandwidth [MBytes/s] STAT |  202.3123 |  101.0394 |  101.2729 |  101.1561 |
| Memory write data volume [GBytes] STAT |    1.5342 |    0.7662 |    0.7680 |    0.7671 |
|    Memory bandwidth [MBytes/s] STAT    | 1169.2951 |  538.8414 |  630.4537 |  584.6476 |
|    Memory data volume [GBytes] STAT    |    8.8671 |    4.0862 |    4.7809 |    4.4336 |
|       Operational intensity STAT       |    3.7886 |    1.7454 |    2.0432 |    1.8943 |
+----------------------------------------+-----------+-----------+-----------+-----------+
