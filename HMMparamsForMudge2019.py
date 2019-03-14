#!/usr/bin/env python
# Copyright 2019 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
SampleHMMparams.py
This file has the HMM parameters used to generate the PhyloCSF Regions described in
Mudge-2019.
"""
HMMparams = {# These numbers result from EstimateHMMparams.estimate_hmm_params_for_genome
             # with 4-state hmm
             'Human' : # GENCODE v16
                (0.0018577729491349902, 56.83017481932195,
                 [0.6417673350475721, 0.2150725662530122, 0.1431600986994177],
                 [3376.1207312280044, 122153.58921094926, 328.20142626163494]),
             'Mouse' : #  GENCODE vM5
                (0.0021896341760155695, 59.329475343093236,
                 [0.643717544086626, 0.23112850841659197, 0.12515394749677808],
                 [2684.764830024852, 103708.05331414141, 291.85282506359187]),
             'Fly' : # Flybase dmel6.11
                (0.026522813874694712, 132.46982781088735,
                 [0.4881037880339513, 0.2959725369315241, 0.21592367503453075],
                 [1204.1041757499499, 13023.816348780832, 58.481616736901515]),
             'Worm' : # Wormbase cele259
                (0.04189592132902114, 67.15347014666943,
                 [0.5235551125702645, 0.3258454578240301, 0.15059942960569747],
                 [428.96501601251316, 3997.1457299207464, 37.266297382251388]),
             'Mosquito' : # Vectorbase agam3
                (0.01266864049727561, 131.33128349181644,
                 [0.4557991824202843, 0.30302892345945, 0.24117189412026765],
                 [1985.8627005310557, 31058.442254764159, 72.295158470383157])
            }
