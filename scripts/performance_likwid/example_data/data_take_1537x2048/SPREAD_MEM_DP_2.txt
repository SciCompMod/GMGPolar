--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+----------+
|    Region Info    |  Core 0  |  Core 14 |
+-------------------+----------+----------+
| RDTSC Runtime [s] | 0.653177 | 0.653154 |
|     call count    |        2 |        2 |
+-------------------+----------+----------+
+------------------------------------------+---------+------------+------------+
|                   Event                  | Counter |   Core 0   |   Core 14  |
+------------------------------------------+---------+------------+------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 1250704000 | 1055541000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  674966800 |  714833900 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  |  923795300 |  979209400 |
|              PWR_PKG_ENERGY              |   PWR0  |    20.8014 |    17.9576 |
|              PWR_DRAM_ENERGY             |   PWR3  |     3.1231 |     2.4900 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      23000 |          0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  277630400 |  268540500 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      15090 |          0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     141851 |          0 |
|               CAS_COUNT_RD               | MBOX0C0 |    1149202 |      79525 |
|               CAS_COUNT_WR               | MBOX0C1 |    1502954 |     115773 |
|               CAS_COUNT_RD               | MBOX1C0 |    1105570 |      79888 |
|               CAS_COUNT_WR               | MBOX1C1 |    1452919 |     116709 |
|               CAS_COUNT_RD               | MBOX2C0 |    1097734 |      80511 |
|               CAS_COUNT_WR               | MBOX2C1 |    1438430 |     116427 |
|               CAS_COUNT_RD               | MBOX3C0 |    1007230 |      67542 |
|               CAS_COUNT_WR               | MBOX3C1 |    1411028 |     109108 |
|               CAS_COUNT_RD               | MBOX4C0 |    1093999 |      67906 |
|               CAS_COUNT_WR               | MBOX4C1 |    1531899 |     109509 |
|               CAS_COUNT_RD               | MBOX5C0 |    1014289 |      67060 |
|               CAS_COUNT_WR               | MBOX5C1 |    1417799 |     109227 |
+------------------------------------------+---------+------------+------------+
+-----------------------------------------------+---------+------------+------------+------------+-------------+
|                     Event                     | Counter |     Sum    |     Min    |     Max    |     Avg     |
+-----------------------------------------------+---------+------------+------------+------------+-------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 2306245000 | 1055541000 | 1250704000 |  1153122500 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 1389800700 |  674966800 |  714833900 |   694900350 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 1903004700 |  923795300 |  979209400 |   951502350 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    38.7590 |    17.9576 |    20.8014 |     19.3795 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     5.6131 |     2.4900 |     3.1231 |      2.8066 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |      23000 |          0 |      23000 |       11500 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  546170900 |  268540500 |  277630400 |   273085450 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |      15090 |          0 |      15090 |        7545 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |     141851 |          0 |     141851 |  70925.5000 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |    1228727 |      79525 |    1149202 | 614363.5000 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |    1618727 |     115773 |    1502954 | 809363.5000 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |    1185458 |      79888 |    1105570 |      592729 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |    1569628 |     116709 |    1452919 |      784814 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |    1178245 |      80511 |    1097734 | 589122.5000 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |    1554857 |     116427 |    1438430 | 777428.5000 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |    1074772 |      67542 |    1007230 |      537386 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |    1520136 |     109108 |    1411028 |      760068 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |    1161905 |      67906 |    1093999 | 580952.5000 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |    1641408 |     109509 |    1531899 |      820704 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |    1081349 |      67060 |    1014289 | 540674.5000 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |    1527026 |     109227 |    1417799 |      763513 |
+-----------------------------------------------+---------+------------+------------+------------+-------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |    0.6532 |    0.6532 |
|        Runtime unhalted [s]       |    0.2602 |    0.2756 |
|            Clock [MHz]            | 1895.0787 | 1893.4339 |
|                CPI                |    0.5397 |    0.6772 |
|             Energy [J]            |   20.8014 |   17.9576 |
|             Power [W]             |   31.8465 |   27.4937 |
|          Energy DRAM [J]          |    3.1231 |    2.4900 |
|           Power DRAM [W]          |    4.7814 |    3.8122 |
|            DP [MFLOP/s]           |  426.9463 |  411.1445 |
|          AVX DP [MFLOP/s]         |    1.8298 |         0 |
|          Packed [MUOPS/s]         |    0.2755 |         0 |
|          Scalar [MUOPS/s]         |  425.0461 |  411.1445 |
|  Memory read bandwidth [MBytes/s] |  633.7539 |   43.3522 |
|  Memory read data volume [GBytes] |    0.4140 |    0.0283 |
| Memory write bandwidth [MBytes/s] |  857.8406 |   66.3124 |
| Memory write data volume [GBytes] |    0.5603 |    0.0433 |
|    Memory bandwidth [MBytes/s]    | 1491.5945 |  109.6646 |
|    Memory data volume [GBytes]    |    0.9743 |    0.0716 |
|       Operational intensity       |    0.2862 |    3.7491 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |    1.3064 |    0.6532 |    0.6532 |    0.6532 |
|        Runtime unhalted [s] STAT       |    0.5358 |    0.2602 |    0.2756 |    0.2679 |
|            Clock [MHz] STAT            | 3788.5126 | 1893.4339 | 1895.0787 | 1894.2563 |
|                CPI STAT                |    1.2169 |    0.5397 |    0.6772 |    0.6084 |
|             Energy [J] STAT            |   38.7590 |   17.9576 |   20.8014 |   19.3795 |
|             Power [W] STAT             |   59.3402 |   27.4937 |   31.8465 |   29.6701 |
|          Energy DRAM [J] STAT          |    5.6131 |    2.4900 |    3.1231 |    2.8066 |
|           Power DRAM [W] STAT          |    8.5936 |    3.8122 |    4.7814 |    4.2968 |
|            DP [MFLOP/s] STAT           |  838.0908 |  411.1445 |  426.9463 |  419.0454 |
|          AVX DP [MFLOP/s] STAT         |    1.8298 |         0 |    1.8298 |    0.9149 |
|          Packed [MUOPS/s] STAT         |    0.2755 |         0 |    0.2755 |    0.1378 |
|          Scalar [MUOPS/s] STAT         |  836.1906 |  411.1445 |  425.0461 |  418.0953 |
|  Memory read bandwidth [MBytes/s] STAT |  677.1061 |   43.3522 |  633.7539 |  338.5531 |
|  Memory read data volume [GBytes] STAT |    0.4423 |    0.0283 |    0.4140 |    0.2211 |
| Memory write bandwidth [MBytes/s] STAT |  924.1530 |   66.3124 |  857.8406 |  462.0765 |
| Memory write data volume [GBytes] STAT |    0.6036 |    0.0433 |    0.5603 |    0.3018 |
|    Memory bandwidth [MBytes/s] STAT    | 1601.2591 |  109.6646 | 1491.5945 |  800.6295 |
|    Memory data volume [GBytes] STAT    |    1.0459 |    0.0716 |    0.9743 |    0.5230 |
|       Operational intensity STAT       |    4.0353 |    0.2862 |    3.7491 |    2.0176 |
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
| RDTSC Runtime [s] | 26.930680 | 26.930650 |
|     call count    |         1 |         1 |
+-------------------+-----------+-----------+
+------------------------------------------+---------+--------------+--------------+
|                   Event                  | Counter |    Core 0    |    Core 14   |
+------------------------------------------+---------+--------------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 120497100000 | 118104200000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  47260570000 |  50012340000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  |  64673460000 |  68437860000 |
|              PWR_PKG_ENERGY              |   PWR0  |     865.4214 |     770.2034 |
|              PWR_DRAM_ENERGY             |   PWR3  |     125.0579 |     125.1614 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        65560 |            0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  31951170000 |  31965840000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |        83556 |            0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       417868 |            0 |
|               CAS_COUNT_RD               | MBOX0C0 |     73785000 |     70283970 |
|               CAS_COUNT_WR               | MBOX0C1 |     19513020 |     16725890 |
|               CAS_COUNT_RD               | MBOX1C0 |     73648560 |     70310430 |
|               CAS_COUNT_WR               | MBOX1C1 |     19388110 |     16800500 |
|               CAS_COUNT_RD               | MBOX2C0 |     73607850 |     70349530 |
|               CAS_COUNT_WR               | MBOX2C1 |     19377680 |     16795860 |
|               CAS_COUNT_RD               | MBOX3C0 |     73202810 |     69708880 |
|               CAS_COUNT_WR               | MBOX3C1 |     19128300 |     16625820 |
|               CAS_COUNT_RD               | MBOX4C0 |     73647060 |     69623440 |
|               CAS_COUNT_WR               | MBOX4C1 |     19683070 |     16572670 |
|               CAS_COUNT_RD               | MBOX5C0 |     73314150 |     69685420 |
|               CAS_COUNT_WR               | MBOX5C1 |     19250470 |     16666220 |
+------------------------------------------+---------+--------------+--------------+
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
|                     Event                     | Counter |      Sum     |      Min     |      Max     |      Avg     |
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 238601300000 | 118104200000 | 120497100000 | 119300650000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  |  97272910000 |  47260570000 |  50012340000 |  48636455000 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 133111320000 |  64673460000 |  68437860000 |  66555660000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    1635.6248 |     770.2034 |     865.4214 |     817.8124 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     250.2193 |     125.0579 |     125.1614 |     125.1097 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |        65560 |            0 |        65560 |        32780 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  63917010000 |  31951170000 |  31965840000 |  31958505000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |        83556 |            0 |        83556 |        41778 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |       417868 |            0 |       417868 |       208934 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |    144068970 |     70283970 |     73785000 |     72034485 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |     36238910 |     16725890 |     19513020 |     18119455 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |    143958990 |     70310430 |     73648560 |     71979495 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |     36188610 |     16800500 |     19388110 |     18094305 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |    143957380 |     70349530 |     73607850 |     71978690 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |     36173540 |     16795860 |     19377680 |     18086770 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |    142911690 |     69708880 |     73202810 |     71455845 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |     35754120 |     16625820 |     19128300 |     17877060 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |    143270500 |     69623440 |     73647060 |     71635250 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |     36255740 |     16572670 |     19683070 |     18127870 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |    142999570 |     69685420 |     73314150 |     71499785 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |     35916690 |     16666220 |     19250470 |     17958345 |
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |   26.9307 |   26.9306 |
|        Runtime unhalted [s]       |   18.2213 |   19.2822 |
|            Clock [MHz]            | 1895.3675 | 1895.4020 |
|                CPI                |    0.3922 |    0.4235 |
|             Energy [J]            |  865.4214 |  770.2034 |
|             Power [W]             |   32.1351 |   28.5995 |
|          Energy DRAM [J]          |  125.0579 |  125.1614 |
|           Power DRAM [W]          |    4.6437 |    4.6475 |
|            DP [MFLOP/s]           | 1186.5641 | 1186.9688 |
|          AVX DP [MFLOP/s]         |    0.1365 |         0 |
|          Packed [MUOPS/s]         |    0.0211 |         0 |
|          Scalar [MUOPS/s]         | 1186.4227 | 1186.9688 |
|  Memory read bandwidth [MBytes/s] | 1048.5122 |  998.0282 |
|  Memory read data volume [GBytes] |   28.2371 |   26.8775 |
| Memory write bandwidth [MBytes/s] |  276.4803 |  238.0917 |
| Memory write data volume [GBytes] |    7.4458 |    6.4120 |
|    Memory bandwidth [MBytes/s]    | 1324.9925 | 1236.1199 |
|    Memory data volume [GBytes]    |   35.6829 |   33.2895 |
|       Operational intensity       |    0.8955 |    0.9602 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |   53.8613 |   26.9306 |   26.9307 |   26.9306 |
|        Runtime unhalted [s] STAT       |   37.5035 |   18.2213 |   19.2822 |   18.7518 |
|            Clock [MHz] STAT            | 3790.7695 | 1895.3675 | 1895.4020 | 1895.3848 |
|                CPI STAT                |    0.8157 |    0.3922 |    0.4235 |    0.4078 |
|             Energy [J] STAT            | 1635.6248 |  770.2034 |  865.4214 |  817.8124 |
|             Power [W] STAT             |   60.7346 |   28.5995 |   32.1351 |   30.3673 |
|          Energy DRAM [J] STAT          |  250.2193 |  125.0579 |  125.1614 |  125.1097 |
|           Power DRAM [W] STAT          |    9.2912 |    4.6437 |    4.6475 |    4.6456 |
|            DP [MFLOP/s] STAT           | 2373.5329 | 1186.5641 | 1186.9688 | 1186.7665 |
|          AVX DP [MFLOP/s] STAT         |    0.1365 |         0 |    0.1365 |    0.0683 |
|          Packed [MUOPS/s] STAT         |    0.0211 |         0 |    0.0211 |    0.0106 |
|          Scalar [MUOPS/s] STAT         | 2373.3915 | 1186.4227 | 1186.9688 | 1186.6958 |
|  Memory read bandwidth [MBytes/s] STAT | 2046.5404 |  998.0282 | 1048.5122 | 1023.2702 |
|  Memory read data volume [GBytes] STAT |   55.1146 |   26.8775 |   28.2371 |   27.5573 |
| Memory write bandwidth [MBytes/s] STAT |  514.5720 |  238.0917 |  276.4803 |  257.2860 |
| Memory write data volume [GBytes] STAT |   13.8578 |    6.4120 |    7.4458 |    6.9289 |
|    Memory bandwidth [MBytes/s] STAT    | 2561.1124 | 1236.1199 | 1324.9925 | 1280.5562 |
|    Memory data volume [GBytes] STAT    |   68.9724 |   33.2895 |   35.6829 |   34.4862 |
|       Operational intensity STAT       |    1.8557 |    0.8955 |    0.9602 |    0.9279 |
+----------------------------------------+-----------+-----------+-----------+-----------+
