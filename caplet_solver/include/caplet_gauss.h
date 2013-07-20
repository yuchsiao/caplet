/*
Created: Aug 12, 2010
Author : Yu-Chung Hsiao
Email  : project.caplet@gmail.com
*/

/*
This file is part of CAPLET.

CAPLET is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAPLET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef CAPLET_GAUSS_H_
#define CAPLET_GAUSS_H_

namespace gauss{

const int NN = 12; // max degrees
const int N  = (NN+1)/2;

const float p_dummy[N] = {0};
const float w_dummy[N] = {0};

const float p1[N] = {0};
const float w1[N] = {2};

const float p2[N] = {0.5773502691896257};
const float w2[N] = {1};

const float p3[N] = {0, 0.7745966692414834};
const float w3[N] = {0.888888888888889, 0.5555555555555553};

const float p4[N] = {0.3399810435848562, 0.861136311594053};
const float w4[N] = {0.6521451548625462, 0.3478548451374539};

const float p5[N] = {0, 0.5384693101056829, 0.906179845938664};
const float w5[N] = {0.5688888888888889, 0.4786286704993665, 0.2369268850561887};

const float p6[N] = {0.2386191860831968, 0.6612093864662646, 0.932469514203152};
const float w6[N] = {0.4679139345726913, 0.3607615730481379, 0.1713244923791709};

const float p7[N] = {0, 0.4058451513773971, 0.7415311855993945, 0.949107912342759};
const float w7[N] = {0.4179591836734694, 0.3818300505051188, 0.279705391489276, 0.1294849661688697};

const float p8[N] = {0.1834346424956498, 0.5255324099163289, 0.7966664774136262, 0.960289856497537};
const float w8[N] = {0.3626837833783619, 0.3137066458778874, 0.2223810344533786, 0.1012285362903738};

const float p9[N] = {0, 0.3242534234038088, 0.6133714327005908, 0.836031107326635, 0.968160239507627};
const float w9[N] = {0.3302393550012597, 0.3123470770400025, 0.2606106964029353, 0.1806481606948577, 0.0812743883615721};

const float p10[N] = {0.1488743389816312, 0.433395394129247, 0.6794095682990246, 0.865063366688984, 0.973906528517172};
const float w10[N] = {0.2955242247147529, 0.2692667193099968, 0.2190863625159832, 0.149451349150573, 0.06667134430868681};

const float p11[N] = {0, 0.2695431559523449, 0.5190961292068117, 0.73015200557405, 0.887062599768093, 0.978228658146058};
const float w11[N] = {0.2729250867779006, 0.2628045445102466, 0.2331937645919933, 0.1862902109277339, 0.1255803694649132, 0.05566856711616958};

const float p12[N] = {0.1252334085114688, 0.3678314989981804, 0.5873179542866143, 0.7699026741943177, 0.904117256370452, 0.981560634246732};
const float w12[N] = {0.2491470458134027, 0.2334925365383534, 0.2031674267230672, 0.1600783285433586, 0.1069393259953637, 0.04717533638647547};

const float* const p[][N] = {
		{p_dummy}, {p1}, {p2}, {p3}, {p4}, {p5}, {p6}, {p7}, {p8}, {p9}, {p10}, {p11}, {p12}
};
const float* const w[][N] = {
		{p_dummy}, {w1}, {w2}, {w3}, {w4}, {w5}, {w6}, {w7}, {w8}, {w9}, {w10}, {w11}, {w12}
};


}

#endif /* CAPLET_GAUSS_H_ */
