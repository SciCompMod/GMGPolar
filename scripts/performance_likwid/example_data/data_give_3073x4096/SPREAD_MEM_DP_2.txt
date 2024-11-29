--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+----------+
|    Region Info    |  Core 0  |  Core 14 |
+-------------------+----------+----------+
| RDTSC Runtime [s] | 2.057629 | 2.057758 |
|     call count    |        2 |        2 |
+-------------------+----------+----------+
+------------------------------------------+---------+------------+------------+
|                   Event                  | Counter |   Core 0   |   Core 14  |
+------------------------------------------+---------+------------+------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 5305266000 | 4819998000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 2799556000 | 2624223000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 3830557000 | 3591142000 |
|              PWR_PKG_ENERGY              |   PWR0  |    63.7252 |    63.6863 |
|              PWR_DRAM_ENERGY             |   PWR3  |     8.5299 |     7.7218 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |     122920 |          0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 1444643000 | 1406082000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      81488 |          0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |    1080242 |          0 |
|               CAS_COUNT_RD               | MBOX0C0 |    1095178 |     720381 |
|               CAS_COUNT_WR               | MBOX0C1 |    2268937 |    1013794 |
|               CAS_COUNT_RD               | MBOX1C0 |    1050557 |     722675 |
|               CAS_COUNT_WR               | MBOX1C1 |    2221399 |    1015872 |
|               CAS_COUNT_RD               | MBOX2C0 |    1250490 |     720659 |
|               CAS_COUNT_WR               | MBOX2C1 |    2469179 |    1012443 |
|               CAS_COUNT_RD               | MBOX3C0 |    1016022 |     678423 |
|               CAS_COUNT_WR               | MBOX3C1 |    2304798 |    1001079 |
|               CAS_COUNT_RD               | MBOX4C0 |    1094123 |     673282 |
|               CAS_COUNT_WR               | MBOX4C1 |    2409395 |     997774 |
|               CAS_COUNT_RD               | MBOX5C0 |     912249 |     672190 |
|               CAS_COUNT_WR               | MBOX5C1 |    2197546 |     997662 |
+------------------------------------------+---------+------------+------------+
+-----------------------------------------------+---------+-------------+------------+------------+--------------+
|                     Event                     | Counter |     Sum     |     Min    |     Max    |      Avg     |
+-----------------------------------------------+---------+-------------+------------+------------+--------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 10125264000 | 4819998000 | 5305266000 |   5062632000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  |  5423779000 | 2624223000 | 2799556000 |   2711889500 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  |  7421699000 | 3591142000 | 3830557000 |   3710849500 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    127.4115 |    63.6863 |    63.7252 |      63.7058 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     16.2517 |     7.7218 |     8.5299 |       8.1258 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |      122920 |          0 |     122920 |        61460 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  2850725000 | 1406082000 | 1444643000 |   1425362500 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |       81488 |          0 |      81488 |        40744 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |     1080242 |          0 |    1080242 |       540121 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |     1815559 |     720381 |    1095178 |  907779.5000 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |     3282731 |    1013794 |    2268937 | 1.641366e+06 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |     1773232 |     722675 |    1050557 |       886616 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |     3237271 |    1015872 |    2221399 | 1.618636e+06 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |     1971149 |     720659 |    1250490 |  985574.5000 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |     3481622 |    1012443 |    2469179 |      1740811 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |     1694445 |     678423 |    1016022 |  847222.5000 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |     3305877 |    1001079 |    2304798 | 1.652938e+06 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |     1767405 |     673282 |    1094123 |  883702.5000 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |     3407169 |     997774 |    2409395 | 1.703584e+06 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |     1584439 |     672190 |     912249 |  792219.5000 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |     3195208 |     997662 |    2197546 |      1597604 |
+-----------------------------------------------+---------+-------------+------------+------------+--------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |    2.0576 |    2.0578 |
|        Runtime unhalted [s]       |    1.0793 |    1.0117 |
|            Clock [MHz]            | 1895.6496 | 1895.3918 |
|                CPI                |    0.5277 |    0.5444 |
|             Energy [J]            |   63.7252 |   63.6863 |
|             Power [W]             |   30.9702 |   30.9494 |
|          Energy DRAM [J]          |    8.5299 |    7.7218 |
|           Power DRAM [W]          |    4.1455 |    3.7525 |
|            DP [MFLOP/s]           |  706.5689 |  683.3078 |
|          AVX DP [MFLOP/s]         |    4.3584 |         0 |
|          Packed [MUOPS/s]         |    0.6243 |         0 |
|          Scalar [MUOPS/s]         |  702.0911 |  683.3078 |
|  Memory read bandwidth [MBytes/s] |  199.6432 |  130.2423 |
|  Memory read data volume [GBytes] |    0.4108 |    0.2680 |
| Memory write bandwidth [MBytes/s] |  431.4482 |  187.8121 |
| Memory write data volume [GBytes] |    0.8878 |    0.3865 |
|    Memory bandwidth [MBytes/s]    |  631.0914 |  318.0544 |
|    Memory data volume [GBytes]    |    1.2986 |    0.6545 |
|       Operational intensity       |    1.1196 |    2.1484 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |    4.1154 |    2.0576 |    2.0578 |    2.0577 |
|        Runtime unhalted [s] STAT       |    2.0910 |    1.0117 |    1.0793 |    1.0455 |
|            Clock [MHz] STAT            | 3791.0414 | 1895.3918 | 1895.6496 | 1895.5207 |
|                CPI STAT                |    1.0721 |    0.5277 |    0.5444 |    0.5360 |
|             Energy [J] STAT            |  127.4115 |   63.6863 |   63.7252 |   63.7058 |
|             Power [W] STAT             |   61.9196 |   30.9494 |   30.9702 |   30.9598 |
|          Energy DRAM [J] STAT          |   16.2517 |    7.7218 |    8.5299 |    8.1258 |
|           Power DRAM [W] STAT          |    7.8980 |    3.7525 |    4.1455 |    3.9490 |
|            DP [MFLOP/s] STAT           | 1389.8767 |  683.3078 |  706.5689 |  694.9384 |
|          AVX DP [MFLOP/s] STAT         |    4.3584 |         0 |    4.3584 |    2.1792 |
|          Packed [MUOPS/s] STAT         |    0.6243 |         0 |    0.6243 |    0.3121 |
|          Scalar [MUOPS/s] STAT         | 1385.3989 |  683.3078 |  702.0911 |  692.6995 |
|  Memory read bandwidth [MBytes/s] STAT |  329.8855 |  130.2423 |  199.6432 |  164.9427 |
|  Memory read data volume [GBytes] STAT |    0.6788 |    0.2680 |    0.4108 |    0.3394 |
| Memory write bandwidth [MBytes/s] STAT |  619.2603 |  187.8121 |  431.4482 |  309.6301 |
| Memory write data volume [GBytes] STAT |    1.2743 |    0.3865 |    0.8878 |    0.6371 |
|    Memory bandwidth [MBytes/s] STAT    |  949.1458 |  318.0544 |  631.0914 |  474.5729 |
|    Memory data volume [GBytes] STAT    |    1.9531 |    0.6545 |    1.2986 |    0.9766 |
|       Operational intensity STAT       |    3.2680 |    1.1196 |    2.1484 |    1.6340 |
+----------------------------------------+-----------+-----------+-----------+-----------+
--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Solve, Group 1: MEM_DP
+-------------------+------------+------------+
|    Region Info    |   Core 0   |   Core 14  |
+-------------------+------------+------------+
| RDTSC Runtime [s] | 265.830100 | 265.830100 |
|     call count    |          1 |          1 |
+-------------------+------------+------------+
+------------------------------------------+---------+---------------+---------------+
|                   Event                  | Counter |     Core 0    |    Core 14    |
+------------------------------------------+---------+---------------+---------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 1083403000000 | 1085537000000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  |  499668100000 |  494584100000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  |  683775700000 |  676801200000 |
|              PWR_PKG_ENERGY              |   PWR0  |     8317.0450 |     8520.9640 |
|              PWR_DRAM_ENERGY             |   PWR3  |     1075.6880 |     1055.7080 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        218232 |             0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  322798100000 |  322311700000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |        302358 |             0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       2178750 |             0 |
|               CAS_COUNT_RD               | MBOX0C0 |     322248100 |     301901400 |
|               CAS_COUNT_WR               | MBOX0C1 |     189503800 |     160348700 |
|               CAS_COUNT_RD               | MBOX1C0 |     322181300 |     302204700 |
|               CAS_COUNT_WR               | MBOX1C1 |     189845500 |     160523200 |
|               CAS_COUNT_RD               | MBOX2C0 |     322886500 |     301970900 |
|               CAS_COUNT_WR               | MBOX2C1 |     190571200 |     160389600 |
|               CAS_COUNT_RD               | MBOX3C0 |     321366100 |     299963800 |
|               CAS_COUNT_WR               | MBOX3C1 |     190074300 |     159234000 |
|               CAS_COUNT_RD               | MBOX4C0 |     321472700 |     300102100 |
|               CAS_COUNT_WR               | MBOX4C1 |     190083100 |     159350200 |
|               CAS_COUNT_RD               | MBOX5C0 |     320998100 |     300426100 |
|               CAS_COUNT_WR               | MBOX5C1 |     189633600 |     159492500 |
+------------------------------------------+---------+---------------+---------------+
+-----------------------------------------------+---------+---------------+---------------+---------------+---------------+
|                     Event                     | Counter |      Sum      |      Min      |      Max      |      Avg      |
+-----------------------------------------------+---------+---------------+---------------+---------------+---------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 2168940000000 | 1083403000000 | 1085537000000 | 1084470000000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  |  994252200000 |  494584100000 |  499668100000 |  497126100000 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 1360576900000 |  676801200000 |  683775700000 |  680288450000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    16838.0090 |     8317.0450 |     8520.9640 |     8419.0045 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     2131.3960 |     1055.7080 |     1075.6880 |     1065.6980 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |        218232 |             0 |        218232 |        109116 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  645109800000 |  322311700000 |  322798100000 |  322554900000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |        302358 |             0 |        302358 |        151179 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |       2178750 |             0 |       2178750 |       1089375 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |     624149500 |     301901400 |     322248100 |     312074750 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |     349852500 |     160348700 |     189503800 |     174926250 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |     624386000 |     302204700 |     322181300 |     312193000 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |     350368700 |     160523200 |     189845500 |     175184350 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |     624857400 |     301970900 |     322886500 |     312428700 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |     350960800 |     160389600 |     190571200 |     175480400 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |     621329900 |     299963800 |     321366100 |     310664950 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |     349308300 |     159234000 |     190074300 |     174654150 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |     621574800 |     300102100 |     321472700 |     310787400 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |     349433300 |     159350200 |     190083100 |     174716650 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |     621424200 |     300426100 |     320998100 |     310712100 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |     349126100 |     159492500 |     189633600 |     174563050 |
+-----------------------------------------------+---------+---------------+---------------+---------------+---------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |  265.8301 |  265.8301 |
|        Runtime unhalted [s]       |  192.6419 |  190.6818 |
|            Clock [MHz]            | 1895.3909 | 1895.4392 |
|                CPI                |    0.4612 |    0.4556 |
|             Energy [J]            | 8317.0450 | 8520.9640 |
|             Power [W]             |   31.2871 |   32.0542 |
|          Energy DRAM [J]          | 1075.6880 | 1055.7080 |
|           Power DRAM [W]          |    4.0465 |    3.9714 |
|            DP [MFLOP/s]           | 1214.3741 | 1212.4726 |
|          AVX DP [MFLOP/s]         |    0.0701 |         0 |
|          Packed [MUOPS/s]         |    0.0102 |         0 |
|          Scalar [MUOPS/s]         | 1214.3023 | 1212.4726 |
|  Memory read bandwidth [MBytes/s] |  464.9352 |  434.9410 |
|  Memory read data volume [GBytes] |  123.5938 |  115.6204 |
| Memory write bandwidth [MBytes/s] |  274.3916 |  230.9657 |
| Memory write data volume [GBytes] |   72.9415 |   61.3976 |
|    Memory bandwidth [MBytes/s]    |  739.3268 |  665.9068 |
|    Memory data volume [GBytes]    |  196.5353 |  177.0181 |
|       Operational intensity       |    1.6425 |    1.8208 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+------------+-----------+-----------+-----------+
|                 Metric                 |     Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+------------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |   531.6602 |  265.8301 |  265.8301 |  265.8301 |
|        Runtime unhalted [s] STAT       |   383.3237 |  190.6818 |  192.6419 |  191.6619 |
|            Clock [MHz] STAT            |  3790.8301 | 1895.3909 | 1895.4392 | 1895.4151 |
|                CPI STAT                |     0.9168 |    0.4556 |    0.4612 |    0.4584 |
|             Energy [J] STAT            | 16838.0090 | 8317.0450 | 8520.9640 | 8419.0045 |
|             Power [W] STAT             |    63.3413 |   31.2871 |   32.0542 |   31.6707 |
|          Energy DRAM [J] STAT          |  2131.3960 | 1055.7080 | 1075.6880 | 1065.6980 |
|           Power DRAM [W] STAT          |     8.0179 |    3.9714 |    4.0465 |    4.0090 |
|            DP [MFLOP/s] STAT           |  2426.8467 | 1212.4726 | 1214.3741 | 1213.4234 |
|          AVX DP [MFLOP/s] STAT         |     0.0701 |         0 |    0.0701 |    0.0350 |
|          Packed [MUOPS/s] STAT         |     0.0102 |         0 |    0.0102 |    0.0051 |
|          Scalar [MUOPS/s] STAT         |  2426.7749 | 1212.4726 | 1214.3023 | 1213.3875 |
|  Memory read bandwidth [MBytes/s] STAT |   899.8762 |  434.9410 |  464.9352 |  449.9381 |
|  Memory read data volume [GBytes] STAT |   239.2142 |  115.6204 |  123.5938 |  119.6071 |
| Memory write bandwidth [MBytes/s] STAT |   505.3573 |  230.9657 |  274.3916 |  252.6787 |
| Memory write data volume [GBytes] STAT |   134.3391 |   61.3976 |   72.9415 |   67.1696 |
|    Memory bandwidth [MBytes/s] STAT    |  1405.2336 |  665.9068 |  739.3268 |  702.6168 |
|    Memory data volume [GBytes] STAT    |   373.5534 |  177.0181 |  196.5353 |  186.7767 |
|       Operational intensity STAT       |     3.4633 |    1.6425 |    1.8208 |    1.7317 |
+----------------------------------------+------------+-----------+-----------+-----------+
