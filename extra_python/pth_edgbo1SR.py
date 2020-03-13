import numpy as np

def pth_sol(r, m, a, zeta, Ee, L):
    pth = np.sqrt(pow(r,2) + (80*pow(a,2)*zeta)/pow(r,8) - 
     (2536*pow(a,2)*zeta)/(45.*pow(r,7)) - 
     (3254*pow(a,2)*zeta)/(105.*pow(r,6)) - 
     (12371*pow(a,2)*zeta)/(735.*pow(r,5)) + 
     (12673*pow(a,2)*zeta)/(1575.*pow(r,4)) + 
     (266911*pow(a,2)*zeta)/(36750.*pow(r,3)) + 
     (2074*pow(a,2)*zeta)/(525.*pow(r,2)) + (4463*pow(a,2)*zeta)/(2625.*r))*np.sqrt(-1 + (pow(a,2)*pow(Ee,2))/
      (pow(a,2) - 2*r + pow(r,2) + (160*pow(a,4)*zeta)/pow(r,12) - 
        (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
        (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - (80*pow(a,2)*zeta)/pow(r,9) + 
        (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
        (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
        (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
        (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
        (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
        (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
        (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
        (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
        (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
        (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
        (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
        (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
        (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
        (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
        (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
        (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
        (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
        (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
        (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
        (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
        (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
        (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
        (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
        (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
        (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
        (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
        (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
        (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
        (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
        (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
        (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
        (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
        (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
        (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
        (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
        (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
        (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
        (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
        (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
        (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
        (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
        (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4))) - 
     pow(L,2)/(pow(a,2) - 2*r + pow(r,2) + (160*pow(a,4)*zeta)/pow(r,12) - 
        (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
        (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - (80*pow(a,2)*zeta)/pow(r,9) + 
        (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
        (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
        (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
        (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
        (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
        (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
        (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
        (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
        (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
        (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
        (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
        (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
        (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
        (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
        (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
        (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
        (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
        (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
        (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
        (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
        (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
        (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
        (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
        (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
        (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
        (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
        (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
        (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
        (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
        (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
        (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
        (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
        (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
        (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
        (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
        (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
        (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
        (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
        (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
        (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
        (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
        (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4))) + 
     (2*pow(a,2)*pow(Ee,2))/
      (r*(pow(a,2) - 2*r + pow(r,2) + (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (4*a*Ee*L)/(r*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (2*pow(L,2))/
      (r*(pow(a,2) - 2*r + pow(r,2) + (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (pow(Ee,2)*pow(r,2))/
      (pow(a,2) - 2*r + pow(r,2) + (160*pow(a,4)*zeta)/pow(r,12) - 
        (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
        (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - (80*pow(a,2)*zeta)/pow(r,9) + 
        (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
        (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
        (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
        (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
        (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
        (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
        (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
        (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
        (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
        (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
        (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
        (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
        (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
        (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
        (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
        (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
        (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
        (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
        (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
        (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
        (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
        (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
        (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
        (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
        (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
        (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
        (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
        (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
        (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
        (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
        (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
        (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
        (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
        (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
        (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
        (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
        (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
        (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
        (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
        (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
        (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
        (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4))) - 
     (80*pow(a,2)*pow(L,2)*zeta)/
      (pow(r,11)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (3736*pow(a,2)*pow(L,2)*zeta)/
      (45.*pow(r,10)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (382*pow(a,2)*pow(L,2)*zeta)/
      (21.*pow(r,9)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (80*pow(a,2)*pow(Ee,2)*zeta)/
      (pow(r,8)*(pow(a,2) - 2*r + pow(r,2) + (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (293*pow(a,2)*pow(L,2)*zeta)/
      (21.*pow(r,8)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (2536*pow(a,2)*pow(Ee,2)*zeta)/
      (45.*pow(r,7)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (160*a*Ee*L*zeta)/
      (3.*pow(r,7)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (80*pow(L,2)*zeta)/
      (3.*pow(r,7)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (281221*pow(a,2)*pow(L,2)*zeta)/
      (11025.*pow(r,7)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (3254*pow(a,2)*pow(Ee,2)*zeta)/
      (105.*pow(r,6)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (96*a*Ee*L*zeta)/
      (5.*pow(r,6)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (32*pow(L,2)*zeta)/
      (5.*pow(r,6)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (17511*pow(a,2)*pow(L,2)*zeta)/
      (2450.*pow(r,6)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (12371*pow(a,2)*pow(Ee,2)*zeta)/
      (735.*pow(r,5)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (12*a*Ee*L*zeta)/
      (pow(r,5)*(pow(a,2) - 2*r + pow(r,2) + (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (22*pow(L,2)*zeta)/
      (5.*pow(r,5)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (59329*pow(a,2)*pow(L,2)*zeta)/
      (18375.*pow(r,5)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (12673*pow(a,2)*pow(Ee,2)*zeta)/
      (1575.*pow(r,4)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (56*a*Ee*L*zeta)/
      (3.*pow(r,4)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (26*pow(L,2)*zeta)/
      (3.*pow(r,4)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (10588*pow(a,2)*pow(L,2)*zeta)/
      (2625.*pow(r,4)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (266911*pow(a,2)*pow(Ee,2)*zeta)/
      (36750.*pow(r,3)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (6*a*Ee*L*zeta)/
      (5.*pow(r,3)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) - 
     (pow(L,2)*zeta)/
      (3.*pow(r,3)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (3267*pow(a,2)*pow(L,2)*zeta)/
      (1750.*pow(r,3)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (2074*pow(a,2)*pow(Ee,2)*zeta)/
      (525.*pow(r,2)*(pow(a,2) - 2*r + pow(r,2) + 
          (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))) + 
     (4463*pow(a,2)*pow(Ee,2)*zeta)/
      (2625.*r*(pow(a,2) - 2*r + pow(r,2) + (160*pow(a,4)*zeta)/pow(r,12) - 
          (3872*pow(a,4)*zeta)/(45.*pow(r,11)) - 
          (37612*pow(a,4)*zeta)/(315.*pow(r,10)) - 
          (80*pow(a,2)*zeta)/pow(r,9) + (68*pow(a,4)*zeta)/(7.*pow(r,9)) + 
          (7336*pow(a,2)*zeta)/(45.*pow(r,8)) + 
          (716267*pow(a,4)*zeta)/(11025.*pow(r,8)) - 
          (20422*pow(a,2)*zeta)/(315.*pow(r,7)) + 
          (87764*pow(a,4)*zeta)/(2205.*pow(r,7)) + 
          (1917*pow(a,2)*zeta)/(245.*pow(r,6)) + 
          (25349*pow(a,4)*zeta)/(36750.*pow(r,6)) - (80*zeta)/(3.*pow(r,5)) - 
          (253756*pow(a,2)*zeta)/(11025.*pow(r,5)) - 
          (69187*pow(a,4)*zeta)/(6125.*pow(r,5)) + (32*zeta)/(5.*pow(r,4)) + 
          (838039*pow(a,2)*zeta)/(110250.*pow(r,4)) - 
          (20389*pow(a,4)*zeta)/(2625.*pow(r,4)) + (22*zeta)/(5.*pow(r,3)) - 
          (18551*pow(a,2)*zeta)/(5250.*pow(r,3)) - 
          (3267*pow(a,4)*zeta)/(1750.*pow(r,3)) + (26*zeta)/(3.*pow(r,2)) - 
          (3048*pow(a,2)*zeta)/(875.*pow(r,2)) + zeta/(3.*r) - 
          (pow(a,2)*zeta)/(6.*r) + (6400*pow(a,4)*pow(zeta,2))/pow(r,19) - 
          (100352*pow(a,4)*pow(zeta,2))/(9.*pow(r,18)) + 
          (10550272*pow(a,4)*pow(zeta,2))/(14175.*pow(r,17)) + 
          (111387328*pow(a,4)*pow(zeta,2))/(33075.*pow(r,16)) - 
          (6400*pow(a,2)*pow(zeta,2))/(3.*pow(r,15)) + 
          (127640476*pow(a,4)*pow(zeta,2))/(33075.*pow(r,15)) + 
          (73600*pow(a,2)*pow(zeta,2))/(27.*pow(r,14)) - 
          (3746962904*pow(a,4)*pow(zeta,2))/(3.472875e6*pow(r,14)) + 
          (53504*pow(a,2)*pow(zeta,2))/(175.*pow(r,13)) - 
          (12269284921*pow(a,4)*pow(zeta,2))/(5.788125e6*pow(r,13)) + 
          (5159968*pow(a,2)*pow(zeta,2))/(11025.*pow(r,12)) - 
          (13550983498*pow(a,4)*pow(zeta,2))/(1.3505625e7*pow(r,12)) - 
          (4788428*pow(a,2)*pow(zeta,2))/(3675.*pow(r,11)) + 
          (18389183644*pow(a,4)*pow(zeta,2))/(1.21550625e8*pow(r,11)) - 
          (7569622*pow(a,2)*pow(zeta,2))/(23625.*pow(r,10)) + 
          (111727950998*pow(a,4)*pow(zeta,2))/(2.02584375e8*pow(r,10)) - 
          (1719318*pow(a,2)*pow(zeta,2))/(30625.*pow(r,9)) + 
          (223785533413*pow(a,4)*pow(zeta,2))/(8.103375e8*pow(r,9)) + 
          (20121062*pow(a,2)*pow(zeta,2))/(118125.*pow(r,8)) + 
          (47728183259*pow(a,4)*pow(zeta,2))/(1.012921875e9*pow(r,8)) + 
          (86896648*pow(a,2)*pow(zeta,2))/(826875.*pow(r,7)) - 
          (2166678709*pow(a,4)*pow(zeta,2))/(4.8234375e7*pow(r,7)) + 
          (8176739*pow(a,2)*pow(zeta,2))/(183750.*pow(r,6)) - 
          (6749467699*pow(a,4)*pow(zeta,2))/(1.929375e8*pow(r,6)) + 
          (42136*pow(a,2)*pow(zeta,2))/(2625.*pow(r,5)) - 
          (14010347*pow(a,4)*pow(zeta,2))/(984375.*pow(r,5)) + 
          (4463*pow(a,2)*pow(zeta,2))/(7875.*pow(r,4)) - 
          (4860207*pow(a,4)*pow(zeta,2))/(1.53125e6*pow(r,4)))))
    return(pth)
