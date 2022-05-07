/*
 * Config
 *
 * =============================================================================
 * Copyright (c) 2021 by Cryptographic Engineering Research Group (CERG)
 * ECE Department, George Mason University
 * Fairfax, VA, U.S.A.
 * Author: Duc Tri Nguyen
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * =============================================================================
 * @author   Duc Tri Nguyen <dnguye69@gmu.edu>
 */

#ifndef CONFIG_H
#define CONFIG_H

/*
 * Enable FMA/FMS instruction, this affect accuracy but better speed
 */
#define FMA 1

/*
 * Define Table for Forward and Inverse NTT, for better caching purpose
 */
#ifndef FALCON_LOGN
#define FALCON_LOGN 9
#endif

#define FALCON_N (1 << FALCON_LOGN)
#define FALCON_Q 12289
#define FALCON_QINV (-12287)
#define FALCON_V 5461
#define FALCON_MONT 4091   // pow(2, 16, 12289)
#define FALCON_MONT2 10952 // pow(4, 16, 12289)
#define FALCON_MONT_QINV 10908

#define _APPLE_M1_ 0

#endif
