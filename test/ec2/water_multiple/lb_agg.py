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

physbam_time = [decimal.Decimal('0'), decimal.Decimal('24'), decimal.Decimal('43'), decimal.Decimal('63'),
                decimal.Decimal('82'), decimal.Decimal('102'), decimal.Decimal('122'), decimal.Decimal('142'),
                decimal.Decimal('162'), decimal.Decimal('182'), decimal.Decimal('202'), decimal.Decimal('222'),
                decimal.Decimal('242'), decimal.Decimal('262'), decimal.Decimal('282'), decimal.Decimal('302'),
                decimal.Decimal('322'), decimal.Decimal('342'), decimal.Decimal('362'), decimal.Decimal('382'),
                decimal.Decimal('402'), decimal.Decimal('423'), decimal.Decimal('443'), decimal.Decimal('463'),
                decimal.Decimal('565'), decimal.Decimal('752'), decimal.Decimal('941'), decimal.Decimal('1129'),
                decimal.Decimal('1317'), decimal.Decimal('1508'), decimal.Decimal('1767'),
                decimal.Decimal('1963'), decimal.Decimal('2154'), decimal.Decimal('2326'),
                decimal.Decimal('2516'), decimal.Decimal('2707'), decimal.Decimal('2900'),
                decimal.Decimal('3091'), decimal.Decimal('3282'), decimal.Decimal('3473'),
                decimal.Decimal('3664'), decimal.Decimal('3857'), decimal.Decimal('4048'),
                decimal.Decimal('4236'), decimal.Decimal('4429'), decimal.Decimal('4618'),
                decimal.Decimal('4809'), decimal.Decimal('5000'), decimal.Decimal('5190'),
                decimal.Decimal('5379'), decimal.Decimal('5572'), decimal.Decimal('5763'),
                decimal.Decimal('5951'), decimal.Decimal('6141'), decimal.Decimal('6334'),
                decimal.Decimal('6526'), decimal.Decimal('6718'), decimal.Decimal('6909'),
                decimal.Decimal('7101'), decimal.Decimal('7292'), decimal.Decimal('7482'),
                decimal.Decimal('7673'), decimal.Decimal('7863'), decimal.Decimal('8052'),
                decimal.Decimal('8245'), decimal.Decimal('8436'), decimal.Decimal('8625'),
                decimal.Decimal('8882'), decimal.Decimal('9079'), decimal.Decimal('9267'),
                decimal.Decimal('9459'), decimal.Decimal('9650'), decimal.Decimal('9844'),
                decimal.Decimal('10034'), decimal.Decimal('10224'), decimal.Decimal('10417'),
                decimal.Decimal('10607'), decimal.Decimal('10799'), decimal.Decimal('10990'),
                decimal.Decimal('11180'), decimal.Decimal('11370'), decimal.Decimal('11560'),
                decimal.Decimal('11749'), decimal.Decimal('11936'), decimal.Decimal('12131'),
                decimal.Decimal('12318'), decimal.Decimal('12513'), decimal.Decimal('12704'),
                decimal.Decimal('12895'), decimal.Decimal('13085'), decimal.Decimal('13274'),
                decimal.Decimal('13469'), decimal.Decimal('13662'), decimal.Decimal('13857'),
                decimal.Decimal('14050'), decimal.Decimal('14242'), decimal.Decimal('14432'),
                decimal.Decimal('14625'), decimal.Decimal('14818'), decimal.Decimal('15008'),
                decimal.Decimal('15199'), decimal.Decimal('15392'), decimal.Decimal('15584'),
                decimal.Decimal('15773'), decimal.Decimal('15967'), decimal.Decimal('16158'),
                decimal.Decimal('16351'), decimal.Decimal('16541'), decimal.Decimal('16734'),
                decimal.Decimal('16923'), decimal.Decimal('17114'), decimal.Decimal('17304'),
                decimal.Decimal('17497')]


nimbus_time = [decimal.Decimal('50.618775845'), decimal.Decimal('79.250428916'),
               decimal.Decimal('105.085516215'), decimal.Decimal('132.119888545'),
               decimal.Decimal('157.538612605'), decimal.Decimal('183.246523858'),
               decimal.Decimal('208.665585518'), decimal.Decimal('233.963278533'),
               decimal.Decimal('259.615947247'), decimal.Decimal('285.539758683'),
               decimal.Decimal('310.893000603'), decimal.Decimal('336.769743920'),
               decimal.Decimal('362.568735838'), decimal.Decimal('388.047821045'),
               decimal.Decimal('413.337495089'), decimal.Decimal('438.909369708'),
               decimal.Decimal('464.364185095'), decimal.Decimal('505.840686083'),
               decimal.Decimal('699.793387414'), decimal.Decimal('861.270599366'),
               decimal.Decimal('1004.267471314'), decimal.Decimal('1140.576018572'),
               decimal.Decimal('1261.492698670'), decimal.Decimal('1368.155442477'),
               decimal.Decimal('1468.186063052'), decimal.Decimal('1575.955931426'),
               decimal.Decimal('1670.974357367'), decimal.Decimal('1765.077263594'),
               decimal.Decimal('1927.124108792'), decimal.Decimal('1967.215937138'),
               decimal.Decimal('2005.547369004'), decimal.Decimal('2044.631623507'),
               decimal.Decimal('2086.216301680'), decimal.Decimal('2127.070912123'),
               decimal.Decimal('2171.488625050'), decimal.Decimal('2211.573472500'),
               decimal.Decimal('2253.399616004'), decimal.Decimal('2298.655289650'),
               decimal.Decimal('2345.497198582'), decimal.Decimal('2389.459612847'),
               decimal.Decimal('2432.871437073'), decimal.Decimal('2480.860258580'),
               decimal.Decimal('2530.780672789'), decimal.Decimal('2578.828957797'),
               decimal.Decimal('2625.646914006'), decimal.Decimal('2674.393803120'),
               decimal.Decimal('2726.307642222'), decimal.Decimal('2777.545148612'),
               decimal.Decimal('2826.687381506'), decimal.Decimal('2876.586103440'),
               decimal.Decimal('2923.956985474'), decimal.Decimal('2977.479158640'),
               decimal.Decimal('3029.281610966'), decimal.Decimal('3070.971880913'),
               decimal.Decimal('3112.283965111'), decimal.Decimal('3152.906273365'),
               decimal.Decimal('3193.365506888'), decimal.Decimal('3233.388521910'),
               decimal.Decimal('3275.027091265'), decimal.Decimal('3315.997920275'),
               decimal.Decimal('3355.516775847'), decimal.Decimal('3393.765501500'),
               decimal.Decimal('3433.720615387'), decimal.Decimal('3472.773760796'),
               decimal.Decimal('3512.572277785'), decimal.Decimal('3614.194466830'),
               decimal.Decimal('3655.754209280'), decimal.Decimal('3698.621353388'),
               decimal.Decimal('3740.161655665'), decimal.Decimal('3781.867182255'),
               decimal.Decimal('3825.403713227'), decimal.Decimal('3866.218448163'),
               decimal.Decimal('3906.771632433'), decimal.Decimal('3946.560483695'),
               decimal.Decimal('3984.599286318'), decimal.Decimal('4024.534130335'),
               decimal.Decimal('4068.228868485'), decimal.Decimal('4110.975626708'),
               decimal.Decimal('4152.994430304'), decimal.Decimal('4195.869360924'),
               decimal.Decimal('4236.568170310'), decimal.Decimal('4277.415185214'),
               decimal.Decimal('4320.381929875'), decimal.Decimal('4359.984858990'),
               decimal.Decimal('4399.751604319'), decimal.Decimal('4435.784503937'),
               decimal.Decimal('4468.297230721'), decimal.Decimal('4504.662152529'),
               decimal.Decimal('4541.211733103'), decimal.Decimal('4579.170803786'),
               decimal.Decimal('4618.069453717'), decimal.Decimal('4657.076435566'),
               decimal.Decimal('4698.736326695'), decimal.Decimal('4739.461994410'),
               decimal.Decimal('4778.201165438'), decimal.Decimal('4817.431811572'),
               decimal.Decimal('4856.024112463'), decimal.Decimal('4894.217893839'),
               decimal.Decimal('4932.192692996'), decimal.Decimal('4971.043957949'),
               decimal.Decimal('5008.334688902'), decimal.Decimal('5046.405546666'),
               decimal.Decimal('5085.756400347'), decimal.Decimal('5126.863268137'),
               decimal.Decimal('5167.299972058'), decimal.Decimal('5205.185098649'),
               decimal.Decimal('5244.992057086'), decimal.Decimal('5284.866942645'),
               decimal.Decimal('5326.913691521'), decimal.Decimal('5369.488580466'),
               decimal.Decimal('5409.463281155')]


# no_lb_time = [decimal.Decimal('57.917835950'), decimal.Decimal('91.446925878'),
#               decimal.Decimal('121.726398468'), decimal.Decimal('153.140338182'),
#               decimal.Decimal('183.183759927'), decimal.Decimal('213.767978668'),
#               decimal.Decimal('244.179278373'), decimal.Decimal('274.742840290'),
#               decimal.Decimal('305.279839277'), decimal.Decimal('335.560076713'),
#               decimal.Decimal('366.390034914'), decimal.Decimal('397.089244842'),
#               decimal.Decimal('427.908239841'), decimal.Decimal('457.833162069'),
#               decimal.Decimal('488.411057472'), decimal.Decimal('527.793320179'),
#               decimal.Decimal('580.275755167'), decimal.Decimal('633.495539665'),
#               decimal.Decimal('687.026797533'), decimal.Decimal('740.016180753'),
#               decimal.Decimal('792.802663087'), decimal.Decimal('846.483469724'),
#               decimal.Decimal('899.664909362'), decimal.Decimal('953.109672069'),
#               decimal.Decimal('1006.067237615'), decimal.Decimal('1060.013171911'),
#               decimal.Decimal('1113.499240636'), decimal.Decimal('1166.699974060'),
#               decimal.Decimal('1257.412032842'), decimal.Decimal('1309.123740434'),
#               decimal.Decimal('1363.282684564'), decimal.Decimal('1417.588613987'),
#               decimal.Decimal('1472.151320457'), decimal.Decimal('1526.087681770'),
#               decimal.Decimal('1580.636390209'), decimal.Decimal('1635.215444564'),
#               decimal.Decimal('1689.953913211'), decimal.Decimal('1744.417998790'),
#               decimal.Decimal('1799.956670761'), decimal.Decimal('1854.724833965'),
#               decimal.Decimal('1909.837639331'), decimal.Decimal('1964.713756561'),
#               decimal.Decimal('2020.122782468'), decimal.Decimal('2074.246917009'),
#               decimal.Decimal('2128.961943864'), decimal.Decimal('2182.754194736'),
#               decimal.Decimal('2237.743175745'), decimal.Decimal('2294.518662452'),
#               decimal.Decimal('2350.588138341'), decimal.Decimal('2405.688771724'),
#               decimal.Decimal('2460.711556911'), decimal.Decimal('2515.261056184'),
#               decimal.Decimal('2570.360728740'), decimal.Decimal('2627.494531154'),
#               decimal.Decimal('2682.632342100'), decimal.Decimal('2737.553836822'),
#               decimal.Decimal('2793.238141059'), decimal.Decimal('2847.781440973'),
#               decimal.Decimal('2902.667300939'), decimal.Decimal('2957.889224052'),
#               decimal.Decimal('3011.747851133'), decimal.Decimal('3067.610474824'),
#               decimal.Decimal('3122.461780786'), decimal.Decimal('3177.087759733'),
#               decimal.Decimal('3230.808935880'), decimal.Decimal('3324.780348300'),
#               decimal.Decimal('3377.474926948'), decimal.Decimal('3433.088095188'),
#               decimal.Decimal('3487.901290893'), decimal.Decimal('3542.191685915'),
#               decimal.Decimal('3597.129776001'), decimal.Decimal('3651.672860145'),
#               decimal.Decimal('3706.363636732'), decimal.Decimal('3761.411725759'),
#               decimal.Decimal('3817.293852806'), decimal.Decimal('3872.329704999'),
#               decimal.Decimal('3927.897531747'), decimal.Decimal('3983.667123794'),
#               decimal.Decimal('4038.578639030'), decimal.Decimal('4093.051434755'),
#               decimal.Decimal('4148.077949047'), decimal.Decimal('4201.949028015'),
#               decimal.Decimal('4256.492571115'), decimal.Decimal('4310.345893144'),
#               decimal.Decimal('4365.054386377'), decimal.Decimal('4419.460530281'),
#               decimal.Decimal('4473.207107543'), decimal.Decimal('4527.032786130'),
#               decimal.Decimal('4581.528235673'), decimal.Decimal('4637.197320938'),
#               decimal.Decimal('4691.968028783'), decimal.Decimal('4745.923574924'),
#               decimal.Decimal('4800.742099285'), decimal.Decimal('4856.121926069'),
#               decimal.Decimal('4912.625047445'), decimal.Decimal('4967.473665952'),
#               decimal.Decimal('5022.286271810'), decimal.Decimal('5076.425526618'),
#               decimal.Decimal('5131.542181015'), decimal.Decimal('5186.030579328'),
#               decimal.Decimal('5240.378130197'), decimal.Decimal('5294.618986129'),
#               decimal.Decimal('5349.215530395'), decimal.Decimal('5403.709317207'),
#               decimal.Decimal('5457.529114246'), decimal.Decimal('5510.965689420'),
#               decimal.Decimal('5565.652396202'), decimal.Decimal('5620.326745986'),
#               decimal.Decimal('5674.726508379'), decimal.Decimal('5729.685341596'),
#               decimal.Decimal('5785.089265346')]
# 
# 
# lb_time = [decimal.Decimal('58.369339466'), decimal.Decimal('92.347867489'),
#            decimal.Decimal('122.579607964'), decimal.Decimal('154.140292883'),
#            decimal.Decimal('183.845167399'), decimal.Decimal('213.608921289'),
#            decimal.Decimal('244.151135683'), decimal.Decimal('274.179576397'),
#            decimal.Decimal('304.982108116'), decimal.Decimal('335.830384493'),
#            decimal.Decimal('366.228529930'), decimal.Decimal('396.695441484'),
#            decimal.Decimal('427.238857269'), decimal.Decimal('458.005823135'),
#            decimal.Decimal('489.033401251'), decimal.Decimal('541.904859304'),
#            decimal.Decimal('594.783228874'), decimal.Decimal('644.649731875'),
#            decimal.Decimal('691.875217438'), decimal.Decimal('741.522109509'),
#            decimal.Decimal('784.138109446'), decimal.Decimal('826.678420305'),
#            decimal.Decimal('870.670739651'), decimal.Decimal('913.058146238'),
#            decimal.Decimal('955.096115351'), decimal.Decimal('997.727478266'),
#            decimal.Decimal('1039.536827326'), decimal.Decimal('1082.406356573'),
#            decimal.Decimal('1166.610794544'), decimal.Decimal('1209.104840994'),
#            decimal.Decimal('1252.494080305'), decimal.Decimal('1295.055364132'),
#            decimal.Decimal('1338.038009167'), decimal.Decimal('1379.949409485'),
#            decimal.Decimal('1421.716799021'), decimal.Decimal('1464.554733753'),
#            decimal.Decimal('1506.775873423'), decimal.Decimal('1549.008729935'),
#            decimal.Decimal('1591.226543188'), decimal.Decimal('1633.135487557'),
#            decimal.Decimal('1675.476188898'), decimal.Decimal('1717.775759697'),
#            decimal.Decimal('1759.501713991'), decimal.Decimal('1801.998857498'),
#            decimal.Decimal('1844.298150778'), decimal.Decimal('1886.262128830'),
#            decimal.Decimal('1928.772400141'), decimal.Decimal('1971.532629013'),
#            decimal.Decimal('2014.230008125'), decimal.Decimal('2056.730457544'),
#            decimal.Decimal('2099.274816275'), decimal.Decimal('2142.480417252'),
#            decimal.Decimal('2186.530419827'), decimal.Decimal('2228.944834709'),
#            decimal.Decimal('2271.038614750'), decimal.Decimal('2314.822394609'),
#            decimal.Decimal('2358.437840462'), decimal.Decimal('2400.924898624'),
#            decimal.Decimal('2443.470822811'), decimal.Decimal('2486.063863993'),
#            decimal.Decimal('2528.290863514'), decimal.Decimal('2571.237352133'),
#            decimal.Decimal('2614.138123512'), decimal.Decimal('2656.451458454'),
#            decimal.Decimal('2698.208684206'), decimal.Decimal('2775.904321194'),
#            decimal.Decimal('2818.498358011'), decimal.Decimal('2860.811707258'),
#            decimal.Decimal('2905.052045584'), decimal.Decimal('2949.332514048'),
#            decimal.Decimal('2992.327894926'), decimal.Decimal('3035.878123999'),
#            decimal.Decimal('3078.859702110'), decimal.Decimal('3121.628813505'),
#            decimal.Decimal('3164.163958073'), decimal.Decimal('3206.525780678'),
#            decimal.Decimal('3250.366357327'), decimal.Decimal('3292.067573786'),
#            decimal.Decimal('3334.433430910'), decimal.Decimal('3376.738303423'),
#            decimal.Decimal('3418.597016573'), decimal.Decimal('3460.502395630'),
#            decimal.Decimal('3502.865504742'), decimal.Decimal('3544.784678459'),
#            decimal.Decimal('3587.300551415'), decimal.Decimal('3628.555335045'),
#            decimal.Decimal('3670.916536093'), decimal.Decimal('3712.800442219'),
#            decimal.Decimal('3754.896687746'), decimal.Decimal('3797.741059065'),
#            decimal.Decimal('3839.596092701'), decimal.Decimal('3881.809545279'),
#            decimal.Decimal('3924.001489401'), decimal.Decimal('3966.298411369'),
#            decimal.Decimal('4008.181195736'), decimal.Decimal('4049.769968033'),
#            decimal.Decimal('4091.763573170'), decimal.Decimal('4133.467551947'),
#            decimal.Decimal('4175.128409386'), decimal.Decimal('4216.798864842'),
#            decimal.Decimal('4259.504122019'), decimal.Decimal('4302.296735287'),
#            decimal.Decimal('4344.065224171'), decimal.Decimal('4386.549253464'),
#            decimal.Decimal('4428.235793352'), decimal.Decimal('4470.769942284'),
#            decimal.Decimal('4513.425975323'), decimal.Decimal('4555.752370834'),
#            decimal.Decimal('4598.048165321'), decimal.Decimal('4640.171785593'),
#            decimal.Decimal('4682.637504101')]



           
nimbus_time = [t / 60 for t in nimbus_time]
physbam_time = [t / 60 for t in physbam_time]

length = min(len(nimbus_time), len(physbam_time))

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

iter_num = range(0, length)

Legends = []
Legends.append('Nimbus')
Legends.append('PhysBAM')

Parts = []

line = plt.plot(nimbus_time[0:length], iter_num, color='g', linewidth=6)
Parts.append(line[0])
line = plt.plot(physbam_time[0:length], iter_num, color='r', linewidth=6)
Parts.append(line[0])

plt.xlabel('Time (minute)', size=40)
plt.ylabel('Iteration Number', size=40)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
# plt.axis([40, 160, 0, 0.03])
plt.grid(True)
plt.legend(Parts, Legends, loc='lower right', prop={'size':40})

plt.show()


