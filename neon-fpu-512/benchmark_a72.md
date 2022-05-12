# NEON A72

| FFT | Foward FFT | Inverse FFT
|:-------------|----------:|-----------:|
| FFT 0 |       21 |        9
| FFT 1 |       21 |        9
| FFT 2 |       25 |       20
| FFT 3 |       54 |       49
| FFT 4 |       92 |       88
| FFT 5 |      221 |      217
| FFT 6 |      540 |      543
| FFT 7 |     1155 |     1216
| FFT 8 |     2770 |     2913
| FFT 9 |     5951 |     6135
| FFT 10 |    14060 |    14705

| NTT | Foward NTT | Inverse NTT
|:-------------|----------:|-----------:|
| NTT 9 |     3561 |     3776
| NTT 10 |     7688 |     8247

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_add |        0 |       19
| poly_add |        1 |        9
| poly_add |        2 |        9
| poly_add |        3 |       19
| poly_add |        4 |       36
| poly_add |        5 |       68
| poly_add |        6 |      132
| poly_add |        7 |      260
| poly_add |        8 |      516
| poly_add |        9 |     1050
| poly_add |       10 |     2074

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_sub |        0 |       19
| poly_sub |        1 |        9
| poly_sub |        2 |        9
| poly_sub |        3 |       19
| poly_sub |        4 |       36
| poly_sub |        5 |       68
| poly_sub |        6 |      132
| poly_sub |        7 |      260
| poly_sub |        8 |      516
| poly_sub |        9 |     1050
| poly_sub |       10 |     2074

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_neg |        0 |       15
| poly_neg |        1 |        9
| poly_neg |        2 |        9
| poly_neg |        3 |       15
| poly_neg |        4 |       27
| poly_neg |        5 |       51
| poly_neg |        6 |       99
| poly_neg |        7 |      195
| poly_neg |        8 |      387
| poly_neg |        9 |      793
| poly_neg |       10 |     1561

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_adj_fft |        0 |       17
| poly_adj_fft |        1 |        9
| poly_adj_fft |        2 |        9
| poly_adj_fft |        3 |        9
| poly_adj_fft |        4 |       17
| poly_adj_fft |        5 |       29
| poly_adj_fft |        6 |       53
| poly_adj_fft |        7 |      101
| poly_adj_fft |        8 |      197
| poly_adj_fft |        9 |      407
| poly_adj_fft |       10 |      794

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_mul_fft |        0 |       10
| poly_mul_fft |        1 |       12
| poly_mul_fft |        2 |        9
| poly_mul_fft |        3 |       17
| poly_mul_fft |        4 |       48
| poly_mul_fft |        5 |       85
| poly_mul_fft |        6 |      161
| poly_mul_fft |        7 |      313
| poly_mul_fft |        8 |      617
| poly_mul_fft |        9 |     1228
| poly_mul_fft |       10 |     2459

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_invnorm2_fft |        0 |       10
| poly_invnorm2_fft |        1 |        9
| poly_invnorm2_fft |        2 |       16
| poly_invnorm2_fft |        3 |       29
| poly_invnorm2_fft |        4 |      143
| poly_invnorm2_fft |        5 |      290
| poly_invnorm2_fft |        6 |      584
| poly_invnorm2_fft |        7 |     1172
| poly_invnorm2_fft |        8 |     2348
| poly_invnorm2_fft |        9 |     4700
| poly_invnorm2_fft |       10 |     9404

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_mul_autoadj_fft |        0 |       10
| poly_mul_autoadj_fft |        1 |       10
| poly_mul_autoadj_fft |        2 |        9
| poly_mul_autoadj_fft |        3 |       13
| poly_mul_autoadj_fft |        4 |       32
| poly_mul_autoadj_fft |        5 |       59
| poly_mul_autoadj_fft |        6 |      113
| poly_mul_autoadj_fft |        7 |      221
| poly_mul_autoadj_fft |        8 |      437
| poly_mul_autoadj_fft |        9 |      891
| poly_mul_autoadj_fft |       10 |     1756

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_LDL_fft |        0 |       11
| poly_LDL_fft |        1 |       29
| poly_LDL_fft |        2 |       42
| poly_LDL_fft |        3 |       91
| poly_LDL_fft |        4 |      225
| poly_LDL_fft |        5 |      451
| poly_LDL_fft |        6 |      903
| poly_LDL_fft |        7 |     1807
| poly_LDL_fft |        8 |     3615
| poly_LDL_fft |        9 |     7260
| poly_LDL_fft |       10 |    14807

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_LDLmv_fft |        0 |       12
| poly_LDLmv_fft |        1 |       36
| poly_LDLmv_fft |        2 |       43
| poly_LDLmv_fft |        3 |       92
| poly_LDLmv_fft |        4 |      234
| poly_LDLmv_fft |        5 |      466
| poly_LDLmv_fft |        6 |      930
| poly_LDLmv_fft |        7 |     1858
| poly_LDLmv_fft |        8 |     3714
| poly_LDLmv_fft |        9 |     7436
| poly_LDLmv_fft |       10 |    15172

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_split_fft |        0 |       10
| poly_split_fft |        1 |        9
| poly_split_fft |        2 |       12
| poly_split_fft |        3 |       20
| poly_split_fft |        4 |       50
| poly_split_fft |        5 |       86
| poly_split_fft |        6 |      158
| poly_split_fft |        7 |      302
| poly_split_fft |        8 |      590
| poly_split_fft |        9 |     1185
| poly_split_fft |       10 |     2355

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_merge_fft |        0 |       10
| poly_merge_fft |        1 |       10
| poly_merge_fft |        2 |       14
| poly_merge_fft |        3 |       23
| poly_merge_fft |        4 |       40
| poly_merge_fft |        5 |       77
| poly_merge_fft |        6 |      138
| poly_merge_fft |        7 |      260
| poly_merge_fft |        8 |      504
| poly_merge_fft |        9 |      992
| poly_merge_fft |       10 |     2038


|degree|  kg(us)|  ek(us)|  sd(us)| sdc(us)|  st(us)| stc(us)|  vv(us)| vvc(us)|
| ---- |  ----- | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 512: | 18381.76 |   114.92 |   577.59 |   634.32 |   415.52 |   472.67 |    33.91 |    91.89 |
|1024: | 51452.50 |   247.69 |  1174.00 |  1289.24 |   825.50 |   939.59 |    69.65 |   183.83 |

|degree|  kg(kc)|   ek(kc)|  sd(kc)| sdc(kc)|  st(kc)| stc(kc)|  vv(kc)| vvc(kc)|
| ---- | ------ |  ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 512: | 33086.90 |   206.68 |  1039.53 |  1141.62 |   747.79 |   850.68 |    60.93 |   165.31 |
|1024: | 92605.82 |   445.66 |  2113.03 |  2320.46 |  1485.75 |  1691.10 |   125.25 |   330.79 |
