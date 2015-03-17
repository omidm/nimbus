#!/usr/bin/env python

import sys
import argparse
import os
import re
import numpy
import decimal
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


no_lb_time = [decimal.Decimal('57.917835950'), decimal.Decimal('91.446925878'),
              decimal.Decimal('121.726398468'), decimal.Decimal('153.140338182'),
              decimal.Decimal('183.183759927'), decimal.Decimal('213.767978668'),
              decimal.Decimal('244.179278373'), decimal.Decimal('274.742840290'),
              decimal.Decimal('305.279839277'), decimal.Decimal('335.560076713'),
              decimal.Decimal('366.390034914'), decimal.Decimal('397.089244842'),
              decimal.Decimal('427.908239841'), decimal.Decimal('457.833162069'),
              decimal.Decimal('488.411057472'), decimal.Decimal('527.793320179'),
              decimal.Decimal('580.275755167'), decimal.Decimal('633.495539665'),
              decimal.Decimal('687.026797533'), decimal.Decimal('740.016180753'),
              decimal.Decimal('792.802663087'), decimal.Decimal('846.483469724'),
              decimal.Decimal('899.664909362'), decimal.Decimal('953.109672069'),
              decimal.Decimal('1006.067237615'), decimal.Decimal('1060.013171911'),
              decimal.Decimal('1113.499240636'), decimal.Decimal('1166.699974060'),
              decimal.Decimal('1257.412032842'), decimal.Decimal('1309.123740434'),
              decimal.Decimal('1363.282684564'), decimal.Decimal('1417.588613987'),
              decimal.Decimal('1472.151320457'), decimal.Decimal('1526.087681770'),
              decimal.Decimal('1580.636390209'), decimal.Decimal('1635.215444564'),
              decimal.Decimal('1689.953913211'), decimal.Decimal('1744.417998790'),
              decimal.Decimal('1799.956670761'), decimal.Decimal('1854.724833965'),
              decimal.Decimal('1909.837639331'), decimal.Decimal('1964.713756561'),
              decimal.Decimal('2020.122782468'), decimal.Decimal('2074.246917009'),
              decimal.Decimal('2128.961943864'), decimal.Decimal('2182.754194736'),
              decimal.Decimal('2237.743175745'), decimal.Decimal('2294.518662452'),
              decimal.Decimal('2350.588138341'), decimal.Decimal('2405.688771724'),
              decimal.Decimal('2460.711556911'), decimal.Decimal('2515.261056184'),
              decimal.Decimal('2570.360728740'), decimal.Decimal('2627.494531154'),
              decimal.Decimal('2682.632342100'), decimal.Decimal('2737.553836822'),
              decimal.Decimal('2793.238141059'), decimal.Decimal('2847.781440973'),
              decimal.Decimal('2902.667300939'), decimal.Decimal('2957.889224052'),
              decimal.Decimal('3011.747851133'), decimal.Decimal('3067.610474824'),
              decimal.Decimal('3122.461780786'), decimal.Decimal('3177.087759733'),
              decimal.Decimal('3230.808935880'), decimal.Decimal('3324.780348300'),
              decimal.Decimal('3377.474926948'), decimal.Decimal('3433.088095188'),
              decimal.Decimal('3487.901290893'), decimal.Decimal('3542.191685915'),
              decimal.Decimal('3597.129776001'), decimal.Decimal('3651.672860145'),
              decimal.Decimal('3706.363636732'), decimal.Decimal('3761.411725759'),
              decimal.Decimal('3817.293852806'), decimal.Decimal('3872.329704999'),
              decimal.Decimal('3927.897531747'), decimal.Decimal('3983.667123794'),
              decimal.Decimal('4038.578639030'), decimal.Decimal('4093.051434755'),
              decimal.Decimal('4148.077949047'), decimal.Decimal('4201.949028015'),
              decimal.Decimal('4256.492571115'), decimal.Decimal('4310.345893144'),
              decimal.Decimal('4365.054386377'), decimal.Decimal('4419.460530281'),
              decimal.Decimal('4473.207107543'), decimal.Decimal('4527.032786130'),
              decimal.Decimal('4581.528235673'), decimal.Decimal('4637.197320938'),
              decimal.Decimal('4691.968028783'), decimal.Decimal('4745.923574924'),
              decimal.Decimal('4800.742099285'), decimal.Decimal('4856.121926069'),
              decimal.Decimal('4912.625047445'), decimal.Decimal('4967.473665952'),
              decimal.Decimal('5022.286271810'), decimal.Decimal('5076.425526618'),
              decimal.Decimal('5131.542181015'), decimal.Decimal('5186.030579328'),
              decimal.Decimal('5240.378130197'), decimal.Decimal('5294.618986129'),
              decimal.Decimal('5349.215530395'), decimal.Decimal('5403.709317207'),
              decimal.Decimal('5457.529114246'), decimal.Decimal('5510.965689420'),
              decimal.Decimal('5565.652396202'), decimal.Decimal('5620.326745986'),
              decimal.Decimal('5674.726508379'), decimal.Decimal('5729.685341596'),
              decimal.Decimal('5785.089265346')]


lb_time = [decimal.Decimal('58.369339466'), decimal.Decimal('92.347867489'),
           decimal.Decimal('122.579607964'), decimal.Decimal('154.140292883'),
           decimal.Decimal('183.845167399'), decimal.Decimal('213.608921289'),
           decimal.Decimal('244.151135683'), decimal.Decimal('274.179576397'),
           decimal.Decimal('304.982108116'), decimal.Decimal('335.830384493'),
           decimal.Decimal('366.228529930'), decimal.Decimal('396.695441484'),
           decimal.Decimal('427.238857269'), decimal.Decimal('458.005823135'),
           decimal.Decimal('489.033401251'), decimal.Decimal('541.904859304'),
           decimal.Decimal('594.783228874'), decimal.Decimal('644.649731875'),
           decimal.Decimal('691.875217438'), decimal.Decimal('741.522109509'),
           decimal.Decimal('784.138109446'), decimal.Decimal('826.678420305'),
           decimal.Decimal('870.670739651'), decimal.Decimal('913.058146238'),
           decimal.Decimal('955.096115351'), decimal.Decimal('997.727478266'),
           decimal.Decimal('1039.536827326'), decimal.Decimal('1082.406356573'),
           decimal.Decimal('1166.610794544'), decimal.Decimal('1209.104840994'),
           decimal.Decimal('1252.494080305'), decimal.Decimal('1295.055364132'),
           decimal.Decimal('1338.038009167'), decimal.Decimal('1379.949409485'),
           decimal.Decimal('1421.716799021'), decimal.Decimal('1464.554733753'),
           decimal.Decimal('1506.775873423'), decimal.Decimal('1549.008729935'),
           decimal.Decimal('1591.226543188'), decimal.Decimal('1633.135487557'),
           decimal.Decimal('1675.476188898'), decimal.Decimal('1717.775759697'),
           decimal.Decimal('1759.501713991'), decimal.Decimal('1801.998857498'),
           decimal.Decimal('1844.298150778'), decimal.Decimal('1886.262128830'),
           decimal.Decimal('1928.772400141'), decimal.Decimal('1971.532629013'),
           decimal.Decimal('2014.230008125'), decimal.Decimal('2056.730457544'),
           decimal.Decimal('2099.274816275'), decimal.Decimal('2142.480417252'),
           decimal.Decimal('2186.530419827'), decimal.Decimal('2228.944834709'),
           decimal.Decimal('2271.038614750'), decimal.Decimal('2314.822394609'),
           decimal.Decimal('2358.437840462'), decimal.Decimal('2400.924898624'),
           decimal.Decimal('2443.470822811'), decimal.Decimal('2486.063863993'),
           decimal.Decimal('2528.290863514'), decimal.Decimal('2571.237352133'),
           decimal.Decimal('2614.138123512'), decimal.Decimal('2656.451458454'),
           decimal.Decimal('2698.208684206'), decimal.Decimal('2775.904321194'),
           decimal.Decimal('2818.498358011'), decimal.Decimal('2860.811707258'),
           decimal.Decimal('2905.052045584'), decimal.Decimal('2949.332514048'),
           decimal.Decimal('2992.327894926'), decimal.Decimal('3035.878123999'),
           decimal.Decimal('3078.859702110'), decimal.Decimal('3121.628813505'),
           decimal.Decimal('3164.163958073'), decimal.Decimal('3206.525780678'),
           decimal.Decimal('3250.366357327'), decimal.Decimal('3292.067573786'),
           decimal.Decimal('3334.433430910'), decimal.Decimal('3376.738303423'),
           decimal.Decimal('3418.597016573'), decimal.Decimal('3460.502395630'),
           decimal.Decimal('3502.865504742'), decimal.Decimal('3544.784678459'),
           decimal.Decimal('3587.300551415'), decimal.Decimal('3628.555335045'),
           decimal.Decimal('3670.916536093'), decimal.Decimal('3712.800442219'),
           decimal.Decimal('3754.896687746'), decimal.Decimal('3797.741059065'),
           decimal.Decimal('3839.596092701'), decimal.Decimal('3881.809545279'),
           decimal.Decimal('3924.001489401'), decimal.Decimal('3966.298411369'),
           decimal.Decimal('4008.181195736'), decimal.Decimal('4049.769968033'),
           decimal.Decimal('4091.763573170'), decimal.Decimal('4133.467551947'),
           decimal.Decimal('4175.128409386'), decimal.Decimal('4216.798864842'),
           decimal.Decimal('4259.504122019'), decimal.Decimal('4302.296735287'),
           decimal.Decimal('4344.065224171'), decimal.Decimal('4386.549253464'),
           decimal.Decimal('4428.235793352'), decimal.Decimal('4470.769942284'),
           decimal.Decimal('4513.425975323'), decimal.Decimal('4555.752370834'),
           decimal.Decimal('4598.048165321'), decimal.Decimal('4640.171785593'),
           decimal.Decimal('4682.637504101')]



           
lb_time = [t / 60 for t in lb_time]
no_lb_time = [t / 60 for t in no_lb_time]

length = min(len(lb_time), len(no_lb_time))

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

iter_num = range(0, length)

Legends = []
Legends.append('With Load Balancing')
Legends.append('Without Load Balancing')

Parts = []

line = plt.plot(lb_time[0:length], iter_num, color='g', linewidth=3)
Parts.append(line[0])
line = plt.plot(no_lb_time[0:length], iter_num, color='r', linewidth=3)
Parts.append(line[0])

plt.xlabel('Time (minute)')
plt.ylabel('Iteration Number')
# plt.axis([40, 160, 0, 0.03])
plt.grid(True)
plt.legend(Parts, Legends, loc='upper left')

plt.show()


