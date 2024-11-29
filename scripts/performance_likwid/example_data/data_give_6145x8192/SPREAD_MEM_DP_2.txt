--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+----------+
|    Region Info    |  Core 0  |  Core 14 |
+-------------------+----------+----------+
| RDTSC Runtime [s] | 7.848168 | 7.848194 |
|     call count    |        2 |        2 |
+-------------------+----------+----------+
+------------------------------------------+---------+-------------+-------------+
|                   Event                  | Counter |    Core 0   |   Core 14   |
+------------------------------------------+---------+-------------+-------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 21093680000 | 19260230000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 10438780000 |  9620295000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 14285300000 | 13164700000 |
|              PWR_PKG_ENERGY              |   PWR0  |    246.4139 |    216.7717 |
|              PWR_DRAM_ENERGY             |   PWR3  |     33.1498 |     31.4436 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      761098 |           0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  5774999000 |  5626934000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      402100 |           0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     8148557 |           0 |
|               CAS_COUNT_RD               | MBOX0C0 |     4828935 |     3516342 |
|               CAS_COUNT_WR               | MBOX0C1 |     9720440 |     4682268 |
|               CAS_COUNT_RD               | MBOX1C0 |     5575262 |     3516446 |
|               CAS_COUNT_WR               | MBOX1C1 |    10575930 |     4682242 |
|               CAS_COUNT_RD               | MBOX2C0 |     4820838 |     3519059 |
|               CAS_COUNT_WR               | MBOX2C1 |     9720271 |     4684474 |
|               CAS_COUNT_RD               | MBOX3C0 |     4339612 |     3355839 |
|               CAS_COUNT_WR               | MBOX3C1 |     9782633 |     4648751 |
|               CAS_COUNT_RD               | MBOX4C0 |     4287617 |     3358448 |
|               CAS_COUNT_WR               | MBOX4C1 |     9685919 |     4652427 |
|               CAS_COUNT_RD               | MBOX5C0 |     4256381 |     3353432 |
|               CAS_COUNT_WR               | MBOX5C1 |     9620484 |     4646876 |
+------------------------------------------+---------+-------------+-------------+
+-----------------------------------------------+---------+-------------+-------------+-------------+--------------+
|                     Event                     | Counter |     Sum     |     Min     |     Max     |      Avg     |
+-----------------------------------------------+---------+-------------+-------------+-------------+--------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 40353910000 | 19260230000 | 21093680000 |  20176955000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 20059075000 |  9620295000 | 10438780000 |  10029537500 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 27450000000 | 13164700000 | 14285300000 |  13725000000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    463.1856 |    216.7717 |    246.4139 |     231.5928 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     64.5934 |     31.4436 |     33.1498 |      32.2967 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |      761098 |           0 |      761098 |       380549 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  | 11401933000 |  5626934000 |  5774999000 |   5700966500 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |      402100 |           0 |      402100 |       201050 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |     8148557 |           0 |     8148557 | 4.074278e+06 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |     8345277 |     3516342 |     4828935 | 4.172638e+06 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |    14402708 |     4682268 |     9720440 |      7201354 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |     9091708 |     3516446 |     5575262 |      4545854 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |    15258172 |     4682242 |    10575930 |      7629086 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |     8339897 |     3519059 |     4820838 | 4.169948e+06 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |    14404745 |     4684474 |     9720271 | 7.202372e+06 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |     7695451 |     3355839 |     4339612 | 3.847726e+06 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |    14431384 |     4648751 |     9782633 |      7215692 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |     7646065 |     3358448 |     4287617 | 3.823032e+06 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |    14338346 |     4652427 |     9685919 |      7169173 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |     7609813 |     3353432 |     4256381 | 3.804906e+06 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |    14267360 |     4646876 |     9620484 |      7133680 |
+-----------------------------------------------+---------+-------------+-------------+-------------+--------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |    7.8482 |    7.8482 |
|        Runtime unhalted [s]       |    4.0246 |    3.7090 |
|            Clock [MHz]            | 1895.3632 | 1895.4376 |
|                CPI                |    0.4949 |    0.4995 |
|             Energy [J]            |  246.4139 |  216.7717 |
|             Power [W]             |   31.3976 |   27.6206 |
|          Energy DRAM [J]          |   33.1498 |   31.4436 |
|           Power DRAM [W]          |    4.2239 |    4.0065 |
|            DP [MFLOP/s]           |  744.5455 |  716.9718 |
|          AVX DP [MFLOP/s]         |    8.5111 |         0 |
|          Packed [MUOPS/s]         |    1.1865 |         0 |
|          Scalar [MUOPS/s]         |  735.8404 |  716.9718 |
|  Memory read bandwidth [MBytes/s] |  229.2195 |  168.1472 |
|  Memory read data volume [GBytes] |    1.7990 |    1.3197 |
| Memory write bandwidth [MBytes/s] |  481.9932 |  228.3086 |
| Memory write data volume [GBytes] |    3.7828 |    1.7918 |
|    Memory bandwidth [MBytes/s]    |  711.2127 |  396.4559 |
|    Memory data volume [GBytes]    |    5.5817 |    3.1115 |
|       Operational intensity       |    1.0469 |    1.8085 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |   15.6964 |    7.8482 |    7.8482 |    7.8482 |
|        Runtime unhalted [s] STAT       |    7.7336 |    3.7090 |    4.0246 |    3.8668 |
|            Clock [MHz] STAT            | 3790.8008 | 1895.3632 | 1895.4376 | 1895.4004 |
|                CPI STAT                |    0.9944 |    0.4949 |    0.4995 |    0.4972 |
|             Energy [J] STAT            |  463.1856 |  216.7717 |  246.4139 |  231.5928 |
|             Power [W] STAT             |   59.0182 |   27.6206 |   31.3976 |   29.5091 |
|          Energy DRAM [J] STAT          |   64.5934 |   31.4436 |   33.1498 |   32.2967 |
|           Power DRAM [W] STAT          |    8.2304 |    4.0065 |    4.2239 |    4.1152 |
|            DP [MFLOP/s] STAT           | 1461.5173 |  716.9718 |  744.5455 |  730.7586 |
|          AVX DP [MFLOP/s] STAT         |    8.5111 |         0 |    8.5111 |    4.2556 |
|          Packed [MUOPS/s] STAT         |    1.1865 |         0 |    1.1865 |    0.5933 |
|          Scalar [MUOPS/s] STAT         | 1452.8122 |  716.9718 |  735.8404 |  726.4061 |
|  Memory read bandwidth [MBytes/s] STAT |  397.3667 |  168.1472 |  229.2195 |  198.6834 |
|  Memory read data volume [GBytes] STAT |    3.1187 |    1.3197 |    1.7990 |    1.5594 |
| Memory write bandwidth [MBytes/s] STAT |  710.3018 |  228.3086 |  481.9932 |  355.1509 |
| Memory write data volume [GBytes] STAT |    5.5746 |    1.7918 |    3.7828 |    2.7873 |
|    Memory bandwidth [MBytes/s] STAT    | 1107.6686 |  396.4559 |  711.2127 |  553.8343 |
|    Memory data volume [GBytes] STAT    |    8.6932 |    3.1115 |    5.5817 |    4.3466 |
|       Operational intensity STAT       |    2.8554 |    1.0469 |    1.8085 |    1.4277 |
+----------------------------------------+-----------+-----------+-----------+-----------+
--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Solve, Group 1: MEM_DP
+-------------------+-------------+-------------+
|    Region Info    |    Core 0   |   Core 14   |
+-------------------+-------------+-------------+
| RDTSC Runtime [s] | 1041.367000 | 1041.367000 |
|     call count    |           1 |           1 |
+-------------------+-------------+-------------+
+------------------------------------------+---------+---------------+---------------+
|                   Event                  | Counter |     Core 0    |    Core 14    |
+------------------------------------------+---------+---------------+---------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 4226209000000 | 4231077000000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 1932840000000 | 1929193000000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 2644996000000 | 2639955000000 |
|              PWR_PKG_ENERGY              |   PWR0  |    33249.6200 |    29873.7100 |
|              PWR_DRAM_ENERGY             |   PWR3  |     4396.3710 |     4371.0610 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |        714138 |             0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 1259761000000 | 1259140000000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       1138324 |             0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |      10452010 |             0 |
|               CAS_COUNT_RD               | MBOX0C0 |    1411853000 |    1358505000 |
|               CAS_COUNT_WR               | MBOX0C1 |     827646500 |     761853200 |
|               CAS_COUNT_RD               | MBOX1C0 |    1421717000 |    1358239000 |
|               CAS_COUNT_WR               | MBOX1C1 |     839293000 |     761538500 |
|               CAS_COUNT_RD               | MBOX2C0 |    1413194000 |    1358629000 |
|               CAS_COUNT_WR               | MBOX2C1 |     829360300 |     761979100 |
|               CAS_COUNT_RD               | MBOX3C0 |    1409251000 |    1353055000 |
|               CAS_COUNT_WR               | MBOX3C1 |     828888500 |     758839800 |
|               CAS_COUNT_RD               | MBOX4C0 |    1410589000 |    1353165000 |
|               CAS_COUNT_WR               | MBOX4C1 |     830228000 |     759001900 |
|               CAS_COUNT_RD               | MBOX5C0 |    1408875000 |    1352859000 |
|               CAS_COUNT_WR               | MBOX5C1 |     827822200 |     758605600 |
+------------------------------------------+---------+---------------+---------------+
+-----------------------------------------------+---------+---------------+---------------+---------------+---------------+
|                     Event                     | Counter |      Sum      |      Min      |      Max      |      Avg      |
+-----------------------------------------------+---------+---------------+---------------+---------------+---------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 8457286000000 | 4226209000000 | 4231077000000 | 4228643000000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 3862033000000 | 1929193000000 | 1932840000000 | 1931016500000 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 5284951000000 | 2639955000000 | 2644996000000 | 2642475500000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    63123.3300 |    29873.7100 |    33249.6200 |    31561.6650 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     8767.4320 |     4371.0610 |     4396.3710 |     4383.7160 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |        714138 |             0 |        714138 |        357069 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  | 2518901000000 | 1259140000000 | 1259761000000 | 1259450500000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |       1138324 |             0 |       1138324 |        569162 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |      10452010 |             0 |      10452010 |       5226005 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |    2770358000 |    1358505000 |    1411853000 |    1385179000 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |    1589499700 |     761853200 |     827646500 |     794749850 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |    2779956000 |    1358239000 |    1421717000 |    1389978000 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |    1600831500 |     761538500 |     839293000 |     800415750 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |    2771823000 |    1358629000 |    1413194000 |    1385911500 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |    1591339400 |     761979100 |     829360300 |     795669700 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |    2762306000 |    1353055000 |    1409251000 |    1381153000 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |    1587728300 |     758839800 |     828888500 |     793864150 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |    2763754000 |    1353165000 |    1410589000 |    1381877000 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |    1589229900 |     759001900 |     830228000 |     794614950 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |    2761734000 |    1352859000 |    1408875000 |    1380867000 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |    1586427800 |     758605600 |     827822200 |     793213900 |
+-----------------------------------------------+---------+---------------+---------------+---------------+---------------+
+-----------------------------------+------------+------------+
|               Metric              |   Core 0   |   Core 14  |
+-----------------------------------+------------+------------+
|        Runtime (RDTSC) [s]        |  1041.3670 |  1041.3670 |
|        Runtime unhalted [s]       |   745.1845 |   743.7785 |
|            Clock [MHz]            |  1895.4091 |  1895.4452 |
|                CPI                |     0.4573 |     0.4560 |
|             Energy [J]            | 33249.6200 | 29873.7100 |
|             Power [W]             |    31.9288 |    28.6870 |
|          Energy DRAM [J]          |  4396.3710 |  4371.0610 |
|           Power DRAM [W]          |     4.2217 |     4.1974 |
|            DP [MFLOP/s]           |  1209.8046 |  1209.1222 |
|          AVX DP [MFLOP/s]         |     0.0847 |          0 |
|          Packed [MUOPS/s]         |     0.0118 |          0 |
|          Scalar [MUOPS/s]         |  1209.7186 |  1209.1222 |
|  Memory read bandwidth [MBytes/s] |   520.8833 |   499.9245 |
|  Memory read data volume [GBytes] |   542.4307 |   520.6049 |
| Memory write bandwidth [MBytes/s] |   306.2583 |   280.3588 |
| Memory write data volume [GBytes] |   318.9273 |   291.9564 |
|    Memory bandwidth [MBytes/s]    |   827.1416 |   780.2833 |
|    Memory data volume [GBytes]    |   861.3579 |   812.5613 |
|       Operational intensity       |     1.4626 |     1.5496 |
+-----------------------------------+------------+------------+
+----------------------------------------+------------+------------+------------+------------+
|                 Metric                 |     Sum    |     Min    |     Max    |     Avg    |
+----------------------------------------+------------+------------+------------+------------+
|        Runtime (RDTSC) [s] STAT        |  2082.7340 |  1041.3670 |  1041.3670 |  1041.3670 |
|        Runtime unhalted [s] STAT       |  1488.9630 |   743.7785 |   745.1845 |   744.4815 |
|            Clock [MHz] STAT            |  3790.8543 |  1895.4091 |  1895.4452 |  1895.4271 |
|                CPI STAT                |     0.9133 |     0.4560 |     0.4573 |     0.4567 |
|             Energy [J] STAT            | 63123.3300 | 29873.7100 | 33249.6200 | 31561.6650 |
|             Power [W] STAT             |    60.6158 |    28.6870 |    31.9288 |    30.3079 |
|          Energy DRAM [J] STAT          |  8767.4320 |  4371.0610 |  4396.3710 |  4383.7160 |
|           Power DRAM [W] STAT          |     8.4191 |     4.1974 |     4.2217 |     4.2096 |
|            DP [MFLOP/s] STAT           |  2418.9268 |  1209.1222 |  1209.8046 |  1209.4634 |
|          AVX DP [MFLOP/s] STAT         |     0.0847 |          0 |     0.0847 |     0.0423 |
|          Packed [MUOPS/s] STAT         |     0.0118 |          0 |     0.0118 |     0.0059 |
|          Scalar [MUOPS/s] STAT         |  2418.8408 |  1209.1222 |  1209.7186 |  1209.4204 |
|  Memory read bandwidth [MBytes/s] STAT |  1020.8078 |   499.9245 |   520.8833 |   510.4039 |
|  Memory read data volume [GBytes] STAT |  1063.0356 |   520.6049 |   542.4307 |   531.5178 |
| Memory write bandwidth [MBytes/s] STAT |   586.6171 |   280.3588 |   306.2583 |   293.3085 |
| Memory write data volume [GBytes] STAT |   610.8837 |   291.9564 |   318.9273 |   305.4418 |
|    Memory bandwidth [MBytes/s] STAT    |  1607.4249 |   780.2833 |   827.1416 |   803.7124 |
|    Memory data volume [GBytes] STAT    |  1673.9192 |   812.5613 |   861.3579 |   836.9596 |
|       Operational intensity STAT       |     3.0122 |     1.4626 |     1.5496 |     1.5061 |
+----------------------------------------+------------+------------+------------+------------+
