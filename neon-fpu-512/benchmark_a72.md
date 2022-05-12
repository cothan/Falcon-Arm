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
| poly_add |        0 |       18
| poly_add |        1 |        9
| poly_add |        2 |        9
| poly_add |        3 |       18
| poly_add |        4 |       35
| poly_add |        5 |       65
| poly_add |        6 |      125
| poly_add |        7 |      245
| poly_add |        8 |      485
| poly_add |        9 |      985
| poly_add |       10 |     1945

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_sub |        0 |       18
| poly_sub |        1 |        9
| poly_sub |        2 |        9
| poly_sub |        3 |       18
| poly_sub |        4 |       35
| poly_sub |        5 |       65
| poly_sub |        6 |      125
| poly_sub |        7 |      245
| poly_sub |        8 |      485
| poly_sub |        9 |      985
| poly_sub |       10 |     1945

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
| poly_adj_fft |        9 |      389
| poly_adj_fft |       10 |      794

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_mul_fft |        0 |       10
| poly_mul_fft |        1 |       11
| poly_mul_fft |        2 |        9
| poly_mul_fft |        3 |       17
| poly_mul_fft |        4 |       47
| poly_mul_fft |        5 |       83
| poly_mul_fft |        6 |      155
| poly_mul_fft |        7 |      299
| poly_mul_fft |        8 |      587
| poly_mul_fft |        9 |     1163
| poly_mul_fft |       10 |     2332

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_invnorm2_fft |        0 |       10
| poly_invnorm2_fft |        1 |        9
| poly_invnorm2_fft |        2 |       16
| poly_invnorm2_fft |        3 |       29
| poly_invnorm2_fft |        4 |      135
| poly_invnorm2_fft |        5 |      264
| poly_invnorm2_fft |        6 |      522
| poly_invnorm2_fft |        7 |     1038
| poly_invnorm2_fft |        8 |     2070
| poly_invnorm2_fft |        9 |     4134
| poly_invnorm2_fft |       10 |     8262

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_mul_autoadj_fft |        0 |       10
| poly_mul_autoadj_fft |        1 |       10
| poly_mul_autoadj_fft |        2 |        9
| poly_mul_autoadj_fft |        3 |       13
| poly_mul_autoadj_fft |        4 |       31
| poly_mul_autoadj_fft |        5 |       57
| poly_mul_autoadj_fft |        6 |      109
| poly_mul_autoadj_fft |        7 |      213
| poly_mul_autoadj_fft |        8 |      421
| poly_mul_autoadj_fft |        9 |      837
| poly_mul_autoadj_fft |       10 |     1690

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_LDL_fft |        0 |       11
| poly_LDL_fft |        1 |       30
| poly_LDL_fft |        2 |       39
| poly_LDL_fft |        3 |       78
| poly_LDL_fft |        4 |      207
| poly_LDL_fft |        5 |      415
| poly_LDL_fft |        6 |      831
| poly_LDL_fft |        7 |     1663
| poly_LDL_fft |        8 |     3327
| poly_LDL_fft |        9 |     6655
| poly_LDL_fft |       10 |    13339

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_LDLmv_fft |        0 |       12
| poly_LDLmv_fft |        1 |       31
| poly_LDLmv_fft |        2 |       42
| poly_LDLmv_fft |        3 |       80
| poly_LDLmv_fft |        4 |      208
| poly_LDLmv_fft |        5 |      416
| poly_LDLmv_fft |        6 |      832
| poly_LDLmv_fft |        7 |     1664
| poly_LDLmv_fft |        8 |     3329
| poly_LDLmv_fft |        9 |     6660
| poly_LDLmv_fft |       10 |    13334

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_split_fft |        0 |       10
| poly_split_fft |        1 |        9
| poly_split_fft |        2 |       12
| poly_split_fft |        3 |       21
| poly_split_fft |        4 |       51
| poly_split_fft |        5 |       89
| poly_split_fft |        6 |      165
| poly_split_fft |        7 |      317
| poly_split_fft |        8 |      621
| poly_split_fft |        9 |     1229
| poly_split_fft |       10 |     2641

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_merge_fft |        0 |       10
| poly_merge_fft |        1 |        9
| poly_merge_fft |        2 |       19
| poly_merge_fft |        3 |       23
| poly_merge_fft |        4 |       45
| poly_merge_fft |        5 |       78
| poly_merge_fft |        6 |      141
| poly_merge_fft |        7 |      267
| poly_merge_fft |        8 |      519
| poly_merge_fft |        9 |     1042
| poly_merge_fft |       10 |     2077


|degree|  kg(kc)|   ek(kc)|  sd(kc)| sdc(kc)|  st(kc)| stc(kc)|  vv(kc)| vvc(kc)|
| ---- |  ----- | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 512: | 18437.26 |   110.87 |   556.04 |   601.61 |   400.72 |   445.48 |    34.69 |    77.26 |
|1024: | 51469.18 |   239.48 |  1148.44 |  1236.19 |   804.57 |   892.07 |    73.17 |   151.70 |


|degree|  kg(us)|  ek(us)|  sd(us)| sdc(us)|  st(us)| stc(us)|  vv(us)| vvc(us)|
| ---- | ------ |  ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 512: | 33186.82 |   199.43 |  1000.69 |  1082.70 |   721.17 |   801.74 |    62.35 |   138.96 |
|1024: | 92644.28 |   430.95 |  2067.03 |  2224.99 |  1448.08 |  1605.58 |   131.59 |   272.92 |
