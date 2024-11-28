--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+----------+
|    Region Info    |  Core 0  |  Core 14 |
+-------------------+----------+----------+
| RDTSC Runtime [s] | 0.522996 | 0.523040 |
|     call count    |        2 |        2 |
+-------------------+----------+----------+
+------------------------------------------+---------+------------+------------+
|                   Event                  | Counter |   Core 0   |   Core 14  |
+------------------------------------------+---------+------------+------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 1342922000 | 1207888000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  694306800 |  684418700 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  |  950160200 |  936556500 |
|              PWR_PKG_ENERGY              |   PWR0  |    16.5695 |    14.5920 |
|              PWR_DRAM_ENERGY             |   PWR3  |     2.1042 |     2.0066 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      23000 |          0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  361496800 |  351336700 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      15090 |          0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     141851 |          0 |
|               CAS_COUNT_RD               | MBOX0C0 |     215287 |     122205 |
|               CAS_COUNT_WR               | MBOX0C1 |     469622 |     155621 |
|               CAS_COUNT_RD               | MBOX1C0 |     202685 |     121740 |
|               CAS_COUNT_WR               | MBOX1C1 |     451139 |     154942 |
|               CAS_COUNT_RD               | MBOX2C0 |     200606 |     121776 |
|               CAS_COUNT_WR               | MBOX2C1 |     447984 |     155008 |
|               CAS_COUNT_RD               | MBOX3C0 |     268770 |     105458 |
|               CAS_COUNT_WR               | MBOX3C1 |     575209 |     149804 |
|               CAS_COUNT_RD               | MBOX4C0 |     168285 |     106580 |
|               CAS_COUNT_WR               | MBOX4C1 |     451733 |     150302 |
|               CAS_COUNT_RD               | MBOX5C0 |     160929 |     108560 |
|               CAS_COUNT_WR               | MBOX5C1 |     439444 |     153340 |
+------------------------------------------+---------+------------+------------+
+-----------------------------------------------+---------+------------+------------+------------+-------------+
|                     Event                     | Counter |     Sum    |     Min    |     Max    |     Avg     |
+-----------------------------------------------+---------+------------+------------+------------+-------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 2550810000 | 1207888000 | 1342922000 |  1275405000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 1378725500 |  684418700 |  694306800 |   689362750 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 1886716700 |  936556500 |  950160200 |   943358350 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    31.1615 |    14.5920 |    16.5695 |     15.5808 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     4.1108 |     2.0066 |     2.1042 |      2.0554 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |      23000 |          0 |      23000 |       11500 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  712833500 |  351336700 |  361496800 |   356416750 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |      15090 |          0 |      15090 |        7545 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |     141851 |          0 |     141851 |  70925.5000 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |     337492 |     122205 |     215287 |      168746 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |     625243 |     155621 |     469622 | 312621.5000 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |     324425 |     121740 |     202685 | 162212.5000 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |     606081 |     154942 |     451139 | 303040.5000 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |     322382 |     121776 |     200606 |      161191 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |     602992 |     155008 |     447984 |      301496 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |     374228 |     105458 |     268770 |      187114 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |     725013 |     149804 |     575209 | 362506.5000 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |     274865 |     106580 |     168285 | 137432.5000 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |     602035 |     150302 |     451733 | 301017.5000 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |     269489 |     108560 |     160929 | 134744.5000 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |     592784 |     153340 |     439444 |      296392 |
+-----------------------------------------------+---------+------------+------------+------------+-------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |    0.5230 |    0.5230 |
|        Runtime unhalted [s]       |    0.2677 |    0.2639 |
|            Clock [MHz]            | 1895.2254 | 1895.3707 |
|                CPI                |    0.5170 |    0.5666 |
|             Energy [J]            |   16.5695 |   14.5920 |
|             Power [W]             |   31.6818 |   27.8984 |
|          Energy DRAM [J]          |    2.1042 |    2.0066 |
|           Power DRAM [W]          |    4.0233 |    3.8365 |
|            DP [MFLOP/s]           |  693.5764 |  671.7206 |
|          AVX DP [MFLOP/s]         |    2.2852 |         0 |
|          Packed [MUOPS/s]         |    0.3441 |         0 |
|          Scalar [MUOPS/s]         |  691.2032 |  671.7206 |
|  Memory read bandwidth [MBytes/s] |  148.8729 |   83.9791 |
|  Memory read data volume [GBytes] |    0.0779 |    0.0439 |
| Memory write bandwidth [MBytes/s] |  346.9400 |  112.4524 |
| Memory write data volume [GBytes] |    0.1814 |    0.0588 |
|    Memory bandwidth [MBytes/s]    |  495.8129 |  196.4315 |
|    Memory data volume [GBytes]    |    0.2593 |    0.1027 |
|       Operational intensity       |    1.3989 |    3.4196 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |    1.0460 |    0.5230 |    0.5230 |    0.5230 |
|        Runtime unhalted [s] STAT       |    0.5316 |    0.2639 |    0.2677 |    0.2658 |
|            Clock [MHz] STAT            | 3790.5961 | 1895.2254 | 1895.3707 | 1895.2980 |
|                CPI STAT                |    1.0836 |    0.5170 |    0.5666 |    0.5418 |
|             Energy [J] STAT            |   31.1615 |   14.5920 |   16.5695 |   15.5808 |
|             Power [W] STAT             |   59.5802 |   27.8984 |   31.6818 |   29.7901 |
|          Energy DRAM [J] STAT          |    4.1108 |    2.0066 |    2.1042 |    2.0554 |
|           Power DRAM [W] STAT          |    7.8598 |    3.8365 |    4.0233 |    3.9299 |
|            DP [MFLOP/s] STAT           | 1365.2970 |  671.7206 |  693.5764 |  682.6485 |
|          AVX DP [MFLOP/s] STAT         |    2.2852 |         0 |    2.2852 |    1.1426 |
|          Packed [MUOPS/s] STAT         |    0.3441 |         0 |    0.3441 |    0.1721 |
|          Scalar [MUOPS/s] STAT         | 1362.9238 |  671.7206 |  691.2032 |  681.4619 |
|  Memory read bandwidth [MBytes/s] STAT |  232.8520 |   83.9791 |  148.8729 |  116.4260 |
|  Memory read data volume [GBytes] STAT |    0.1218 |    0.0439 |    0.0779 |    0.0609 |
| Memory write bandwidth [MBytes/s] STAT |  459.3924 |  112.4524 |  346.9400 |  229.6962 |
| Memory write data volume [GBytes] STAT |    0.2402 |    0.0588 |    0.1814 |    0.1201 |
|    Memory bandwidth [MBytes/s] STAT    |  692.2444 |  196.4315 |  495.8129 |  346.1222 |
|    Memory data volume [GBytes] STAT    |    0.3620 |    0.1027 |    0.2593 |    0.1810 |
|       Operational intensity STAT       |    4.8185 |    1.3989 |    3.4196 |    2.4093 |
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
| RDTSC Runtime [s] | 69.955300 | 69.955280 |
|     call count    |         1 |         1 |
+-------------------+-----------+-----------+
+------------------------------------------+---------+--------------+--------------+
|                   Event                  | Counter |    Core 0    |    Core 14   |
+------------------------------------------+---------+--------------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 284180100000 | 285447300000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 131491500000 | 128457200000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 179936100000 | 175781700000 |
|              PWR_PKG_ENERGY              |   PWR0  |    2219.5780 |    1990.8430 |
|              PWR_DRAM_ENERGY             |   PWR3  |     275.1535 |     266.6727 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        65560 |            0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  84628280000 |  84346820000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |        83556 |            0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       417868 |            0 |
|               CAS_COUNT_RD               | MBOX0C0 |     52655420 |     43788870 |
|               CAS_COUNT_WR               | MBOX0C1 |     29821090 |     21495040 |
|               CAS_COUNT_RD               | MBOX1C0 |     52680000 |     43920420 |
|               CAS_COUNT_WR               | MBOX1C1 |     29836690 |     21532860 |
|               CAS_COUNT_RD               | MBOX2C0 |     52684950 |     43982330 |
|               CAS_COUNT_WR               | MBOX2C1 |     29809200 |     21620910 |
|               CAS_COUNT_RD               | MBOX3C0 |     53112340 |     43465630 |
|               CAS_COUNT_WR               | MBOX3C1 |     30698570 |     21303350 |
|               CAS_COUNT_RD               | MBOX4C0 |     52503160 |     43259910 |
|               CAS_COUNT_WR               | MBOX4C1 |     29991500 |     21161440 |
|               CAS_COUNT_RD               | MBOX5C0 |     52333280 |     43403130 |
|               CAS_COUNT_WR               | MBOX5C1 |     29705620 |     21198410 |
+------------------------------------------+---------+--------------+--------------+
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
|                     Event                     | Counter |      Sum     |      Min     |      Max     |      Avg     |
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 569627400000 | 284180100000 | 285447300000 | 284813700000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 259948700000 | 128457200000 | 131491500000 | 129974350000 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 355717800000 | 175781700000 | 179936100000 | 177858900000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    4210.4210 |    1990.8430 |    2219.5780 |    2105.2105 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     541.8262 |     266.6727 |     275.1535 |     270.9131 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |        65560 |            0 |        65560 |        32780 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  | 168975100000 |  84346820000 |  84628280000 |  84487550000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |        83556 |            0 |        83556 |        41778 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |       417868 |            0 |       417868 |       208934 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |     96444290 |     43788870 |     52655420 |     48222145 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |     51316130 |     21495040 |     29821090 |     25658065 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |     96600420 |     43920420 |     52680000 |     48300210 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |     51369550 |     21532860 |     29836690 |     25684775 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |     96667280 |     43982330 |     52684950 |     48333640 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |     51430110 |     21620910 |     29809200 |     25715055 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |     96577970 |     43465630 |     53112340 |     48288985 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |     52001920 |     21303350 |     30698570 |     26000960 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |     95763070 |     43259910 |     52503160 |     47881535 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |     51152940 |     21161440 |     29991500 |     25576470 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |     95736410 |     43403130 |     52333280 |     47868205 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |     50904030 |     21198410 |     29705620 |     25452015 |
+-----------------------------------------------+---------+--------------+--------------+--------------+--------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |   69.9553 |   69.9553 |
|        Runtime unhalted [s]       |   50.6981 |   49.5282 |
|            Clock [MHz]            | 1895.3336 | 1895.3572 |
|                CPI                |    0.4627 |    0.4500 |
|             Energy [J]            | 2219.5780 | 1990.8430 |
|             Power [W]             |   31.7285 |   28.4588 |
|          Energy DRAM [J]          |  275.1535 |  266.6727 |
|           Power DRAM [W]          |    3.9333 |    3.8120 |
|            DP [MFLOP/s]           | 1209.8024 | 1205.7249 |
|          AVX DP [MFLOP/s]         |    0.0526 |         0 |
|          Packed [MUOPS/s]         |    0.0081 |         0 |
|          Scalar [MUOPS/s]         | 1209.7479 | 1205.7249 |
|  Memory read bandwidth [MBytes/s] |  289.0707 |  239.5316 |
|  Memory read data volume [GBytes] |   20.2220 |   16.7565 |
| Memory write bandwidth [MBytes/s] |  164.5509 |  117.3888 |
| Memory write data volume [GBytes] |   11.5112 |    8.2120 |
|    Memory bandwidth [MBytes/s]    |  453.6216 |  356.9204 |
|    Memory data volume [GBytes]    |   31.7332 |   24.9685 |
|       Operational intensity       |    2.6670 |    3.3781 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |  139.9106 |   69.9553 |   69.9553 |   69.9553 |
|        Runtime unhalted [s] STAT       |  100.2263 |   49.5282 |   50.6981 |   50.1131 |
|            Clock [MHz] STAT            | 3790.6908 | 1895.3336 | 1895.3572 | 1895.3454 |
|                CPI STAT                |    0.9127 |    0.4500 |    0.4627 |    0.4564 |
|             Energy [J] STAT            | 4210.4210 | 1990.8430 | 2219.5780 | 2105.2105 |
|             Power [W] STAT             |   60.1873 |   28.4588 |   31.7285 |   30.0937 |
|          Energy DRAM [J] STAT          |  541.8262 |  266.6727 |  275.1535 |  270.9131 |
|           Power DRAM [W] STAT          |    7.7453 |    3.8120 |    3.9333 |    3.8727 |
|            DP [MFLOP/s] STAT           | 2415.5273 | 1205.7249 | 1209.8024 | 1207.7636 |
|          AVX DP [MFLOP/s] STAT         |    0.0526 |         0 |    0.0526 |    0.0263 |
|          Packed [MUOPS/s] STAT         |    0.0081 |         0 |    0.0081 |    0.0040 |
|          Scalar [MUOPS/s] STAT         | 2415.4728 | 1205.7249 | 1209.7479 | 1207.7364 |
|  Memory read bandwidth [MBytes/s] STAT |  528.6023 |  239.5316 |  289.0707 |  264.3012 |
|  Memory read data volume [GBytes] STAT |   36.9785 |   16.7565 |   20.2220 |   18.4892 |
| Memory write bandwidth [MBytes/s] STAT |  281.9397 |  117.3888 |  164.5509 |  140.9699 |
| Memory write data volume [GBytes] STAT |   19.7232 |    8.2120 |   11.5112 |    9.8616 |
|    Memory bandwidth [MBytes/s] STAT    |  810.5420 |  356.9204 |  453.6216 |  405.2710 |
|    Memory data volume [GBytes] STAT    |   56.7017 |   24.9685 |   31.7332 |   28.3509 |
|       Operational intensity STAT       |    6.0451 |    2.6670 |    3.3781 |    3.0225 |
+----------------------------------------+-----------+-----------+-----------+-----------+
