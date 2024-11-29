--------------------------------------------------------------------------------
CPU name:	Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
CPU type:	Intel Skylake SP processor
CPU clock:	2.59 GHz
--------------------------------------------------------------------------------
Region Setup, Group 1: MEM_DP
+-------------------+----------+----------+
|    Region Info    |  Core 0  |  Core 14 |
+-------------------+----------+----------+
| RDTSC Runtime [s] | 0.136643 | 0.136752 |
|     call count    |        2 |        2 |
+-------------------+----------+----------+
+------------------------------------------+---------+-----------+-----------+
|                   Event                  | Counter |   Core 0  |  Core 14  |
+------------------------------------------+---------+-----------+-----------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 347081100 | 304000500 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 182132400 | 215580500 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 249252200 | 295027600 |
|              PWR_PKG_ENERGY              |   PWR0  |    4.4886 |    3.9463 |
|              PWR_DRAM_ENERGY             |   PWR3  |    0.5446 |    0.5259 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |      4800 |         0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  90640360 |  87670520 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |      2800 |         0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |     11152 |         0 |
|               CAS_COUNT_RD               | MBOX0C0 |     47143 |     19718 |
|               CAS_COUNT_WR               | MBOX0C1 |     88070 |     20687 |
|               CAS_COUNT_RD               | MBOX1C0 |     40951 |     19609 |
|               CAS_COUNT_WR               | MBOX1C1 |     80423 |     20619 |
|               CAS_COUNT_RD               | MBOX2C0 |     54214 |     19949 |
|               CAS_COUNT_WR               | MBOX2C1 |     95590 |     20977 |
|               CAS_COUNT_RD               | MBOX3C0 |     33937 |     15473 |
|               CAS_COUNT_WR               | MBOX3C1 |     80450 |     19417 |
|               CAS_COUNT_RD               | MBOX4C0 |     32677 |     15297 |
|               CAS_COUNT_WR               | MBOX4C1 |     79206 |     19494 |
|               CAS_COUNT_RD               | MBOX5C0 |     31935 |     15236 |
|               CAS_COUNT_WR               | MBOX5C1 |     77514 |     19212 |
+------------------------------------------+---------+-----------+-----------+
+-----------------------------------------------+---------+-----------+-----------+-----------+------------+
|                     Event                     | Counter |    Sum    |    Min    |    Max    |     Avg    |
+-----------------------------------------------+---------+-----------+-----------+-----------+------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 651081600 | 304000500 | 347081100 |  325540800 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  | 397712900 | 182132400 | 215580500 |  198856450 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  | 544279800 | 249252200 | 295027600 |  272139900 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    8.4349 |    3.9463 |    4.4886 |     4.2174 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |    1.0705 |    0.5259 |    0.5446 |     0.5353 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |      4800 |         0 |      4800 |       2400 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  | 178310880 |  87670520 |  90640360 |   89155440 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |      2800 |         0 |      2800 |       1400 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |     11152 |         0 |     11152 |       5576 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |     66861 |     19718 |     47143 | 33430.5000 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |    108757 |     20687 |     88070 | 54378.5000 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |     60560 |     19609 |     40951 |      30280 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |    101042 |     20619 |     80423 |      50521 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |     74163 |     19949 |     54214 | 37081.5000 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |    116567 |     20977 |     95590 | 58283.5000 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |     49410 |     15473 |     33937 |      24705 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |     99867 |     19417 |     80450 | 49933.5000 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |     47974 |     15297 |     32677 |      23987 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |     98700 |     19494 |     79206 |      49350 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |     47171 |     15236 |     31935 | 23585.5000 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |     96726 |     19212 |     77514 |      48363 |
+-----------------------------------------------+---------+-----------+-----------+-----------+------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |    0.1366 |    0.1368 |
|        Runtime unhalted [s]       |    0.0702 |    0.0831 |
|            Clock [MHz]            | 1895.3079 | 1895.3018 |
|                CPI                |    0.5248 |    0.7091 |
|             Energy [J]            |    4.4886 |    3.9463 |
|             Power [W]             |   32.8490 |   28.8577 |
|          Energy DRAM [J]          |    0.5446 |    0.5259 |
|           Power DRAM [W]          |    3.9852 |    3.8457 |
|            DP [MFLOP/s]           |  664.1422 |  641.0918 |
|          AVX DP [MFLOP/s]         |    0.7349 |         0 |
|          Packed [MUOPS/s]         |    0.1372 |         0 |
|          Scalar [MUOPS/s]         |  663.3370 |  641.0918 |
|  Memory read bandwidth [MBytes/s] |  112.8111 |   49.2721 |
|  Memory read data volume [GBytes] |    0.0154 |    0.0067 |
| Memory write bandwidth [MBytes/s] |  234.7738 |   56.3501 |
| Memory write data volume [GBytes] |    0.0321 |    0.0077 |
|    Memory bandwidth [MBytes/s]    |  347.5849 |  105.6222 |
|    Memory data volume [GBytes]    |    0.0475 |    0.0144 |
|       Operational intensity       |    1.9107 |    6.0697 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |    0.2734 |    0.1366 |    0.1368 |    0.1367 |
|        Runtime unhalted [s] STAT       |    0.1533 |    0.0702 |    0.0831 |    0.0766 |
|            Clock [MHz] STAT            | 3790.6097 | 1895.3018 | 1895.3079 | 1895.3048 |
|                CPI STAT                |    1.2339 |    0.5248 |    0.7091 |    0.6169 |
|             Energy [J] STAT            |    8.4349 |    3.9463 |    4.4886 |    4.2174 |
|             Power [W] STAT             |   61.7067 |   28.8577 |   32.8490 |   30.8533 |
|          Energy DRAM [J] STAT          |    1.0705 |    0.5259 |    0.5446 |    0.5353 |
|           Power DRAM [W] STAT          |    7.8309 |    3.8457 |    3.9852 |    3.9154 |
|            DP [MFLOP/s] STAT           | 1305.2340 |  641.0918 |  664.1422 |  652.6170 |
|          AVX DP [MFLOP/s] STAT         |    0.7349 |         0 |    0.7349 |    0.3674 |
|          Packed [MUOPS/s] STAT         |    0.1372 |         0 |    0.1372 |    0.0686 |
|          Scalar [MUOPS/s] STAT         | 1304.4288 |  641.0918 |  663.3370 |  652.2144 |
|  Memory read bandwidth [MBytes/s] STAT |  162.0832 |   49.2721 |  112.8111 |   81.0416 |
|  Memory read data volume [GBytes] STAT |    0.0221 |    0.0067 |    0.0154 |    0.0111 |
| Memory write bandwidth [MBytes/s] STAT |  291.1239 |   56.3501 |  234.7738 |  145.5619 |
| Memory write data volume [GBytes] STAT |    0.0398 |    0.0077 |    0.0321 |    0.0199 |
|    Memory bandwidth [MBytes/s] STAT    |  453.2071 |  105.6222 |  347.5849 |  226.6036 |
|    Memory data volume [GBytes] STAT    |    0.0619 |    0.0144 |    0.0475 |    0.0309 |
|       Operational intensity STAT       |    7.9804 |    1.9107 |    6.0697 |    3.9902 |
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
| RDTSC Runtime [s] | 19.672270 | 19.672240 |
|     call count    |         1 |         1 |
+-------------------+-----------+-----------+
+------------------------------------------+---------+-------------+-------------+
|                   Event                  | Counter |    Core 0   |   Core 14   |
+------------------------------------------+---------+-------------+-------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 74654430000 | 75272750000 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 37052420000 | 35478990000 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 50703370000 | 48549560000 |
|              PWR_PKG_ENERGY              |   PWR0  |    619.9255 |    554.9335 |
|              PWR_DRAM_ENERGY             |   PWR3  |     72.2812 |     72.2073 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE |   PMC0  |       22816 |           0 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  | 22177770000 | 22013120000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE |   PMC2  |       26128 |           0 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE |   PMC3  |       58742 |           0 |
|               CAS_COUNT_RD               | MBOX0C0 |     3589846 |     2980585 |
|               CAS_COUNT_WR               | MBOX0C1 |     2474371 |     1941146 |
|               CAS_COUNT_RD               | MBOX1C0 |     3555644 |     2987155 |
|               CAS_COUNT_WR               | MBOX1C1 |     2415566 |     1943785 |
|               CAS_COUNT_RD               | MBOX2C0 |     3617958 |     2999196 |
|               CAS_COUNT_WR               | MBOX2C1 |     2512148 |     1977059 |
|               CAS_COUNT_RD               | MBOX3C0 |     3555438 |     2920823 |
|               CAS_COUNT_WR               | MBOX3C1 |     2415628 |     1890843 |
|               CAS_COUNT_RD               | MBOX4C0 |     3675262 |     2926941 |
|               CAS_COUNT_WR               | MBOX4C1 |     2577971 |     1907267 |
|               CAS_COUNT_RD               | MBOX5C0 |     3582328 |     2937690 |
|               CAS_COUNT_WR               | MBOX5C1 |     2434235 |     1909535 |
+------------------------------------------+---------+-------------+-------------+
+-----------------------------------------------+---------+--------------+-------------+-------------+--------------+
|                     Event                     | Counter |      Sum     |     Min     |     Max     |      Avg     |
+-----------------------------------------------+---------+--------------+-------------+-------------+--------------+
|             INSTR_RETIRED_ANY STAT            |  FIXC0  | 149927180000 | 74654430000 | 75272750000 |  74963590000 |
|           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  |  72531410000 | 35478990000 | 37052420000 |  36265705000 |
|           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  |  99252930000 | 48549560000 | 50703370000 |  49626465000 |
|              PWR_PKG_ENERGY STAT              |   PWR0  |    1174.8590 |    554.9335 |    619.9255 |     587.4295 |
|              PWR_DRAM_ENERGY STAT             |   PWR3  |     144.4885 |     72.2073 |     72.2812 |      72.2442 |
| FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |        22816 |           0 |       22816 |        11408 |
|    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  44190890000 | 22013120000 | 22177770000 |  22095445000 |
| FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |        26128 |           0 |       26128 |        13064 |
| FP_ARITH_INST_RETIRED_512B_PACKED_DOUBLE STAT |   PMC3  |        58742 |           0 |       58742 |        29371 |
|               CAS_COUNT_RD STAT               | MBOX0C0 |      6570431 |     2980585 |     3589846 | 3.285216e+06 |
|               CAS_COUNT_WR STAT               | MBOX0C1 |      4415517 |     1941146 |     2474371 | 2.207758e+06 |
|               CAS_COUNT_RD STAT               | MBOX1C0 |      6542799 |     2987155 |     3555644 | 3.271400e+06 |
|               CAS_COUNT_WR STAT               | MBOX1C1 |      4359351 |     1943785 |     2415566 | 2.179676e+06 |
|               CAS_COUNT_RD STAT               | MBOX2C0 |      6617154 |     2999196 |     3617958 |      3308577 |
|               CAS_COUNT_WR STAT               | MBOX2C1 |      4489207 |     1977059 |     2512148 | 2.244604e+06 |
|               CAS_COUNT_RD STAT               | MBOX3C0 |      6476261 |     2920823 |     3555438 | 3.238130e+06 |
|               CAS_COUNT_WR STAT               | MBOX3C1 |      4306471 |     1890843 |     2415628 | 2.153236e+06 |
|               CAS_COUNT_RD STAT               | MBOX4C0 |      6602203 |     2926941 |     3675262 | 3.301102e+06 |
|               CAS_COUNT_WR STAT               | MBOX4C1 |      4485238 |     1907267 |     2577971 |      2242619 |
|               CAS_COUNT_RD STAT               | MBOX5C0 |      6520018 |     2937690 |     3582328 |      3260009 |
|               CAS_COUNT_WR STAT               | MBOX5C1 |      4343770 |     1909535 |     2434235 |      2171885 |
+-----------------------------------------------+---------+--------------+-------------+-------------+--------------+
+-----------------------------------+-----------+-----------+
|               Metric              |   Core 0  |  Core 14  |
+-----------------------------------+-----------+-----------+
|        Runtime (RDTSC) [s]        |   19.6723 |   19.6722 |
|        Runtime unhalted [s]       |   14.2852 |   13.6785 |
|            Clock [MHz]            | 1895.4455 | 1895.4726 |
|                CPI                |    0.4963 |    0.4713 |
|             Energy [J]            |  619.9255 |  554.9335 |
|             Power [W]             |   31.5127 |   28.2090 |
|          Energy DRAM [J]          |   72.2812 |   72.2073 |
|           Power DRAM [W]          |    3.6743 |    3.6705 |
|            DP [MFLOP/s]           | 1127.3935 | 1118.9941 |
|          AVX DP [MFLOP/s]         |    0.0292 |         0 |
|          Packed [MUOPS/s]         |    0.0055 |         0 |
|          Scalar [MUOPS/s]         | 1127.3620 | 1118.9941 |
|  Memory read bandwidth [MBytes/s] |   70.1950 |   57.7541 |
|  Memory read data volume [GBytes] |    1.3809 |    1.1362 |
| Memory write bandwidth [MBytes/s] |   48.2463 |   37.6397 |
| Memory write data volume [GBytes] |    0.9491 |    0.7405 |
|    Memory bandwidth [MBytes/s]    |  118.4413 |   95.3938 |
|    Memory data volume [GBytes]    |    2.3300 |    1.8766 |
|       Operational intensity       |    9.5186 |   11.7303 |
+-----------------------------------+-----------+-----------+
+----------------------------------------+-----------+-----------+-----------+-----------+
|                 Metric                 |    Sum    |    Min    |    Max    |    Avg    |
+----------------------------------------+-----------+-----------+-----------+-----------+
|        Runtime (RDTSC) [s] STAT        |   39.3445 |   19.6722 |   19.6723 |   19.6722 |
|        Runtime unhalted [s] STAT       |   27.9637 |   13.6785 |   14.2852 |   13.9818 |
|            Clock [MHz] STAT            | 3790.9181 | 1895.4455 | 1895.4726 | 1895.4590 |
|                CPI STAT                |    0.9676 |    0.4713 |    0.4963 |    0.4838 |
|             Energy [J] STAT            | 1174.8590 |  554.9335 |  619.9255 |  587.4295 |
|             Power [W] STAT             |   59.7217 |   28.2090 |   31.5127 |   29.8608 |
|          Energy DRAM [J] STAT          |  144.4885 |   72.2073 |   72.2812 |   72.2442 |
|           Power DRAM [W] STAT          |    7.3448 |    3.6705 |    3.6743 |    3.6724 |
|            DP [MFLOP/s] STAT           | 2246.3876 | 1118.9941 | 1127.3935 | 1123.1938 |
|          AVX DP [MFLOP/s] STAT         |    0.0292 |         0 |    0.0292 |    0.0146 |
|          Packed [MUOPS/s] STAT         |    0.0055 |         0 |    0.0055 |    0.0027 |
|          Scalar [MUOPS/s] STAT         | 2246.3561 | 1118.9941 | 1127.3620 | 1123.1780 |
|  Memory read bandwidth [MBytes/s] STAT |  127.9491 |   57.7541 |   70.1950 |   63.9745 |
|  Memory read data volume [GBytes] STAT |    2.5171 |    1.1362 |    1.3809 |    1.2586 |
| Memory write bandwidth [MBytes/s] STAT |   85.8860 |   37.6397 |   48.2463 |   42.9430 |
| Memory write data volume [GBytes] STAT |    1.6896 |    0.7405 |    0.9491 |    0.8448 |
|    Memory bandwidth [MBytes/s] STAT    |  213.8351 |   95.3938 |  118.4413 |  106.9176 |
|    Memory data volume [GBytes] STAT    |    4.2066 |    1.8766 |    2.3300 |    2.1033 |
|       Operational intensity STAT       |   21.2489 |    9.5186 |   11.7303 |   10.6244 |
+----------------------------------------+-----------+-----------+-----------+-----------+
