# NEON M1 

| FFT | Foward FFT | Inverse FFT
|:-------------|----------:|-----------:|
| FFT 0 |       76 |       77
| FFT 1 |       76 |       77
| FFT 2 |       79 |       80
| FFT 3 |       91 |       86
| FFT 4 |      107 |      105
| FFT 5 |      145 |      144
| FFT 6 |      234 |      220
| FFT 7 |      402 |      381
| FFT 8 |      765 |      757
| FFT 9 |     1500 |     1530
| FFT 10 |     3386 |     3412

| NTT | Foward NTT | Inverse NTT
|:-------------|----------:|-----------:|
| NTT 9 |      840 |      848
| NTT 10 |     1693 |     1756

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_add |        0 |       79
| poly_add |        1 |       76
| poly_add |        2 |       76
| poly_add |        3 |       79
| poly_add |        4 |       82
| poly_add |        5 |       87
| poly_add |        6 |      103
| poly_add |        7 |      134
| poly_add |        8 |      256
| poly_add |        9 |      369
| poly_add |       10 |      612

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_sub |        0 |       79
| poly_sub |        1 |       76
| poly_sub |        2 |       76
| poly_sub |        3 |       79
| poly_sub |        4 |       82
| poly_sub |        5 |       87
| poly_sub |        6 |      103
| poly_sub |        7 |      134
| poly_sub |        8 |      256
| poly_sub |        9 |      370
| poly_sub |       10 |      636

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_neg |        0 |       76
| poly_neg |        1 |       76
| poly_neg |        2 |       75
| poly_neg |        3 |       76
| poly_neg |        4 |       82
| poly_neg |        5 |       83
| poly_neg |        6 |       90
| poly_neg |        7 |      111
| poly_neg |        8 |      241
| poly_neg |        9 |      361
| poly_neg |       10 |      633

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_adj_fft |        0 |       77
| poly_adj_fft |        1 |       75
| poly_adj_fft |        2 |       76
| poly_adj_fft |        3 |       75
| poly_adj_fft |        4 |       77
| poly_adj_fft |        5 |       83
| poly_adj_fft |        6 |       86
| poly_adj_fft |        7 |       91
| poly_adj_fft |        8 |      181
| poly_adj_fft |        9 |      164
| poly_adj_fft |       10 |      315

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_mul_fft |        0 |       77
| poly_mul_fft |        1 |       76
| poly_mul_fft |        2 |       76
| poly_mul_fft |        3 |       77
| poly_mul_fft |        4 |       83
| poly_mul_fft |        5 |       91
| poly_mul_fft |        6 |      110
| poly_mul_fft |        7 |      153
| poly_mul_fft |        8 |      283
| poly_mul_fft |        9 |      411
| poly_mul_fft |       10 |      659

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_invnorm2_fft |        0 |       77
| poly_invnorm2_fft |        1 |       78
| poly_invnorm2_fft |        2 |       83
| poly_invnorm2_fft |        3 |       82
| poly_invnorm2_fft |        4 |       87
| poly_invnorm2_fft |        5 |       95
| poly_invnorm2_fft |        6 |      114
| poly_invnorm2_fft |        7 |      157
| poly_invnorm2_fft |        8 |      260
| poly_invnorm2_fft |        9 |      378
| poly_invnorm2_fft |       10 |      587

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_mul_autoadj_fft |        0 |       77
| poly_mul_autoadj_fft |        1 |       77
| poly_mul_autoadj_fft |        2 |       75
| poly_mul_autoadj_fft |        3 |       79
| poly_mul_autoadj_fft |        4 |       82
| poly_mul_autoadj_fft |        5 |       87
| poly_mul_autoadj_fft |        6 |       99
| poly_mul_autoadj_fft |        7 |      134
| poly_mul_autoadj_fft |        8 |      261
| poly_mul_autoadj_fft |        9 |      400
| poly_mul_autoadj_fft |       10 |      722

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_LDL_fft |        0 |       76
| poly_LDL_fft |        1 |       93
| poly_LDL_fft |        2 |       87
| poly_LDL_fft |        3 |       92
| poly_LDL_fft |        4 |       96
| poly_LDL_fft |        5 |      119
| poly_LDL_fft |        6 |      167
| poly_LDL_fft |        7 |      240
| poly_LDL_fft |        8 |      475
| poly_LDL_fft |        9 |      739
| poly_LDL_fft |       10 |     1305

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_LDLmv_fft |        0 |       77
| poly_LDLmv_fft |        1 |       92
| poly_LDLmv_fft |        2 |       88
| poly_LDLmv_fft |        3 |       93
| poly_LDLmv_fft |        4 |       95
| poly_LDLmv_fft |        5 |      120
| poly_LDLmv_fft |        6 |      174
| poly_LDLmv_fft |        7 |      244
| poly_LDLmv_fft |        8 |      449
| poly_LDLmv_fft |        9 |      741
| poly_LDLmv_fft |       10 |     1311

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_split_fft |        0 |       77
| poly_split_fft |        1 |       74
| poly_split_fft |        2 |       77
| poly_split_fft |        3 |       79
| poly_split_fft |        4 |       89
| poly_split_fft |        5 |      106
| poly_split_fft |        6 |      142
| poly_split_fft |        7 |      199
| poly_split_fft |        8 |      303
| poly_split_fft |        9 |      511
| poly_split_fft |       10 |      927

| Function | logn | cycles |
|:-------------|----------:|-----------:|
| poly_merge_fft |        0 |       77
| poly_merge_fft |        1 |       74
| poly_merge_fft |        2 |       77
| poly_merge_fft |        3 |       79
| poly_merge_fft |        4 |       87
| poly_merge_fft |        5 |      102
| poly_merge_fft |        6 |      131
| poly_merge_fft |        7 |      183
| poly_merge_fft |        8 |      274
| poly_merge_fft |        9 |      471
| poly_merge_fft |       10 |      810


|degree|  kg(kc)|   ek(kc)|  sd(kc)| sdc(kc)|  st(kc)| stc(kc)|  vv(kc)| vvc(kc)|
| ---- | ------ |  ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 512: | 15750.38 |    44.36 |   437.71 |   480.56 |   357.99 |   402.36 |    22.64 |    68.19 |
|1024: | 48666.79 |    98.22 |   882.08 |   969.50 |   710.68 |   797.54 |    42.88 |   135.46 |

|degree|  kg(us)|  ek(us)|  sd(us)| sdc(us)|  st(us)| stc(us)|  vv(us)| vvc(us)|
| ---- |  ----- | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 512: |  4923.00 |    15.17 |   138.75 |   151.79 |   113.12 |   128.04 |     8.21 |    22.33 |
|1024: | 15266.08 |    32.50 |   279.33 |   306.58 |   224.50 |   253.21 |    14.54 |    43.38 |

|degree|  kg(Ghz)|   ek(Ghz)|   sd(Ghz)|  sdc(Ghz)|   st(Ghz)|  stc(Ghz)|   vv(Ghz)|  vvc(Ghz)|
| ---- |  ----- | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 512: |     3.14 |      3.15 |     3.14  |    3.13 |     3.12 |    3.14 |     3.17  |    3.14 |
|1024: |     3.16 |      3.14 |     3.12  |    3.13 |     3.13 |    3.13 |     3.08  |    3.14 |


## Benchmark 59 bytes message

### N = 512

mlen, smlen = 59, 715

| Median   | sign | verify |
|  --  | --:  | ---: |
| kc |   477.67 |    21.68 |
| us |   151.79 |     6.79 |
| Ghz|     3.15 |     3.19 |

| Average | sign | verify |
|  --  | --:  | ---: |
| kc |   478.29 |    21.67 |
| us |   152.99 |     6.81 |
| Ghz|     3.13 |     3.18 |



### N = 1024

mlen, smlen = 59, 1332

| Median   | sign | verify |
| --   | --:  | ---: |
| kc |   892.07 |    43.33 |
| us |   280.83 |    13.75 |
| Ghz|     3.18 |     3.15 |

| Average | sign | verify |
| --   | --:  | ---: |
| kc |   892.17 |    43.33 |
| us |   282.07 |    13.78 |
| Ghz|     3.16 |     3.14 |