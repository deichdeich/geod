int edgbo1SR(double t, gsl_vector * in_state, gsl_vector * out_state){

double r, P_r, Th, P_Th, Ph, P_Ph;

r = gsl_vector_get(in_state, 0);
P_r = gsl_vector_get(in_state, 1);
Th = gsl_vector_get(in_state, 2);
P_Th = gsl_vector_get(in_state, 3);
Ph = gsl_vector_get(in_state, 4);
P_Ph = gsl_vector_get(in_state, 5);

r_d = (220500*pow(MASS,4)*Pr*pow(r,9)*pow(-2*MASS + r,3))/
   (-7350*pow(MASS,2)*(2*MASS - r)*pow(r,4)*
      (1840*pow(MASS,5)*ZETA - 48*pow(MASS,4)*r*ZETA - 15*MASS*pow(r,4)*ZETA - 
        15*pow(r,5)*ZETA - 30*pow(MASS,3)*pow(r,2)*(pow(r,4) + ZETA) + 
        5*pow(MASS,2)*(3*pow(r,7) - 52*pow(r,3)*ZETA)) + 
     pow(SPIN,2)*(-299880000*pow(MASS,10)*ZETA + 583296000*pow(MASS,9)*r*ZETA - 
        310986200*pow(MASS,8)*pow(r,2)*ZETA + 55258000*pow(MASS,7)*pow(r,3)*ZETA - 
        45095130*pow(MASS,6)*pow(r,4)*ZETA - 1018518*pow(MASS,3)*pow(r,7)*ZETA + 
        327054*pow(MASS,2)*pow(r,8)*ZETA + 132321*MASS*pow(r,9)*ZETA + 
        55125*pow(r,10)*ZETA - pow(MASS,4)*pow(r,6)*
         (110250*pow(r,4) + 3428093*ZETA) + 
        20*pow(MASS,5)*(11025*pow(r,9) + 1463804*pow(r,5)*ZETA)) + 
     3*pow(SPIN,2)*MASS*(2*MASS - r)*(149940000*pow(MASS,8)*ZETA - 
        201978000*pow(MASS,7)*r*ZETA + 101014900*pow(MASS,6)*pow(r,2)*ZETA - 
        18766650*pow(MASS,5)*pow(r,3)*ZETA - 55626*pow(MASS,2)*pow(r,6)*ZETA + 
        150696*MASS*pow(r,7)*ZETA + 187446*pow(r,8)*ZETA + 
        30*pow(MASS,4)*(2450*pow(r,8) + 394463*pow(r,4)*ZETA) - 
        5*pow(MASS,3)*(7350*pow(r,9) + 1509019*pow(r,5)*ZETA))*pow(cos(Th),2));
Th_d = (2*PTh)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2) - 
     (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
          3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
          887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
          435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*(1 + 3*cos(2*Th)))/
      (220500.*pow(MASS,3)*pow(r,8)));
Ph_d = (2*(Jz*(-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
          (ZETA*(7350*pow(MASS,2)*pow(r,4)*
                (400*pow(MASS,4) - 96*pow(MASS,3)*r - 66*pow(MASS,2)*pow(r,2) - 
                  130*MASS*pow(r,3) - 5*pow(r,4)) + 
               pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                  2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                  2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                  355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                  205821*pow(r,8)) + 
               3*pow(SPIN,2)*(8820000*pow(MASS,8) - 8173200*pow(MASS,7)*r - 
                  15803900*pow(MASS,6)*pow(r,2) + 4198950*pow(MASS,5)*pow(r,3) + 
                  4061710*pow(MASS,4)*pow(r,4) + 2275145*pow(MASS,3)*pow(r,5) - 
                  164874*pow(MASS,2)*pow(r,6) - 187446*MASS*pow(r,7) - 
                  187446*pow(r,8))*pow(cos(Th),2)))/
           (110250.*pow(MASS,3)*pow(r,11))) - 
       (SPIN*Ee*(400*pow(MASS,4)*ZETA - 144*pow(MASS,3)*r*ZETA - 140*MASS*pow(r,3)*ZETA - 
            9*pow(r,4)*ZETA + 30*pow(MASS,2)*(pow(r,6) - 3*pow(r,2)*ZETA))*
          pow(sin(Th),2))/(15.*MASS*pow(r,7))))/
   (-(pow(SPIN,2)*pow(-400*pow(MASS,4)*ZETA + 144*pow(MASS,3)*r*ZETA + 
           140*MASS*pow(r,3)*ZETA + 9*pow(r,4)*ZETA - 
           30*pow(MASS,2)*(pow(r,6) - 3*pow(r,2)*ZETA),2)*pow(sin(Th),4))/
      (225.*pow(MASS,2)*pow(r,14)) + 
     (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
        (ZETA*(7350*pow(MASS,2)*pow(r,4)*
              (400*pow(MASS,4) - 96*pow(MASS,3)*r - 66*pow(MASS,2)*pow(r,2) - 
                130*MASS*pow(r,3) - 5*pow(r,4)) + 
             pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 205821*pow(r,8))
               + 3*pow(SPIN,2)*(8820000*pow(MASS,8) - 8173200*pow(MASS,7)*r - 
                15803900*pow(MASS,6)*pow(r,2) + 4198950*pow(MASS,5)*pow(r,3) + 
                4061710*pow(MASS,4)*pow(r,4) + 2275145*pow(MASS,3)*pow(r,5) - 
                164874*pow(MASS,2)*pow(r,6) - 187446*MASS*pow(r,7) - 187446*pow(r,8))
               *pow(cos(Th),2)))/(110250.*pow(MASS,3)*pow(r,11)))*
      (pow(r,2)*pow(sin(Th),2) - 
        (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
             3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
             887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
             435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*(1 + 3*cos(2*Th))*
           pow(sin(Th),2))/(220500.*pow(MASS,3)*pow(r,8)) + 
        pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)));
P_r_d = (pow(PTh,2)*(2*r - (pow(SPIN,2)*(-6213200*pow(MASS,6) - 6833400*pow(MASS,5)*r - 
             5566950*pow(MASS,4)*pow(r,2) + 3548440*pow(MASS,3)*pow(r,3) + 
             4003665*pow(MASS,2)*pow(r,4) + 2613240*MASS*pow(r,5) + 1312122*pow(r,6))
            *ZETA*(-1 + 3*pow(cos(Th),2)))/(110250.*pow(MASS,3)*pow(r,8)) + 
        (4*pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
             3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
             887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
             435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*(-1 + 3*pow(cos(Th),2))
           )/(55125.*pow(MASS,3)*pow(r,9))))/
    pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2) - 
      (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
           3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
           887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
           435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*(-1 + 3*pow(cos(Th),2)))/
       (110250.*pow(MASS,3)*pow(r,8)),2) + 
   (pow(Pr,2)*(-(r/pow(-2*MASS + r,2)) + 1/(-2*MASS + r) + 
        (pow(SPIN,2)*(-1 + pow(cos(Th),2)))/(r*pow(-2*MASS + r,2)) - 
        (2*pow(SPIN,2)*(-r - 2*MASS*pow(cos(Th),2) + r*pow(cos(Th),2)))/
         (r*pow(-2*MASS + r,3)) - (pow(SPIN,2)*
           (-r - 2*MASS*pow(cos(Th),2) + r*pow(cos(Th),2)))/
         (pow(r,2)*pow(-2*MASS + r,2)) + 
        ZETA*((-48*pow(MASS,4) - 60*pow(MASS,3)*r - 780*pow(MASS,2)*pow(r,2) - 
              60*MASS*pow(r,3) - 75*pow(r,4))/
            (15.*pow(MASS,2)*pow(2*MASS - r,2)*pow(r,5)) - 
           (1840*pow(MASS,5) - 48*pow(MASS,4)*r - 30*pow(MASS,3)*pow(r,2) - 
              260*pow(MASS,2)*pow(r,3) - 15*MASS*pow(r,4) - 15*pow(r,5))/
            (3.*pow(MASS,2)*pow(2*MASS - r,2)*pow(r,6)) + 
           (2*(1840*pow(MASS,5) - 48*pow(MASS,4)*r - 30*pow(MASS,3)*pow(r,2) - 
                260*pow(MASS,2)*pow(r,3) - 15*MASS*pow(r,4) - 15*pow(r,5)))/
            (15.*pow(MASS,2)*pow(2*MASS - r,3)*pow(r,5)) - 
           (pow(SPIN,2)*(583296000*pow(MASS,9) - 621972400*pow(MASS,8)*r + 
                165774000*pow(MASS,7)*pow(r,2) - 180380520*pow(MASS,6)*pow(r,3) + 
                146380400*pow(MASS,5)*pow(r,4) - 20568558*pow(MASS,4)*pow(r,5) - 
                7129626*pow(MASS,3)*pow(r,6) + 2616432*pow(MASS,2)*pow(r,7) + 
                1190889*MASS*pow(r,8) + 551250*pow(r,9) - 
                1661688000*pow(MASS,9)*pow(cos(Th),2) + 
                2424046800*pow(MASS,8)*r*pow(cos(Th),2) - 
                1246933800*pow(MASS,7)*pow(r,2)*pow(cos(Th),2) + 
                509213160*pow(MASS,6)*pow(r,3)*pow(cos(Th),2) - 
                403861200*pow(MASS,5)*pow(r,4)*pow(cos(Th),2) + 
                133809174*pow(MASS,4)*pow(r,5)*pow(cos(Th),2) + 
                7497378*pow(MASS,3)*pow(r,6)*pow(cos(Th),2) + 
                5380704*pow(MASS,2)*pow(r,7)*pow(cos(Th),2) - 
                5061042*MASS*pow(r,8)*pow(cos(Th),2)))/
            (110250.*pow(MASS,4)*pow(2*MASS - r,3)*pow(r,9)) + 
           (pow(SPIN,2)*(-299880000*pow(MASS,10) + 583296000*pow(MASS,9)*r - 
                310986200*pow(MASS,8)*pow(r,2) + 55258000*pow(MASS,7)*pow(r,3) - 
                45095130*pow(MASS,6)*pow(r,4) + 29276080*pow(MASS,5)*pow(r,5) - 
                3428093*pow(MASS,4)*pow(r,6) - 1018518*pow(MASS,3)*pow(r,7) + 
                327054*pow(MASS,2)*pow(r,8) + 132321*MASS*pow(r,9) + 
                55125*pow(r,10) + 899640000*pow(MASS,10)*pow(cos(Th),2) - 
                1661688000*pow(MASS,9)*r*pow(cos(Th),2) + 
                1212023400*pow(MASS,8)*pow(r,2)*pow(cos(Th),2) - 
                415644600*pow(MASS,7)*pow(r,3)*pow(cos(Th),2) + 
                127303290*pow(MASS,6)*pow(r,4)*pow(cos(Th),2) - 
                80772240*pow(MASS,5)*pow(r,5)*pow(cos(Th),2) + 
                22301529*pow(MASS,4)*pow(r,6)*pow(cos(Th),2) + 
                1071054*pow(MASS,3)*pow(r,7)*pow(cos(Th),2) + 
                672588*pow(MASS,2)*pow(r,8)*pow(cos(Th),2) - 
                562338*MASS*pow(r,9)*pow(cos(Th),2)))/
            (12250.*pow(MASS,4)*pow(2*MASS - r,3)*pow(r,10)) - 
           (pow(SPIN,2)*(-299880000*pow(MASS,10) + 583296000*pow(MASS,9)*r - 
                310986200*pow(MASS,8)*pow(r,2) + 55258000*pow(MASS,7)*pow(r,3) - 
                45095130*pow(MASS,6)*pow(r,4) + 29276080*pow(MASS,5)*pow(r,5) - 
                3428093*pow(MASS,4)*pow(r,6) - 1018518*pow(MASS,3)*pow(r,7) + 
                327054*pow(MASS,2)*pow(r,8) + 132321*MASS*pow(r,9) + 
                55125*pow(r,10) + 899640000*pow(MASS,10)*pow(cos(Th),2) - 
                1661688000*pow(MASS,9)*r*pow(cos(Th),2) + 
                1212023400*pow(MASS,8)*pow(r,2)*pow(cos(Th),2) - 
                415644600*pow(MASS,7)*pow(r,3)*pow(cos(Th),2) + 
                127303290*pow(MASS,6)*pow(r,4)*pow(cos(Th),2) - 
                80772240*pow(MASS,5)*pow(r,5)*pow(cos(Th),2) + 
                22301529*pow(MASS,4)*pow(r,6)*pow(cos(Th),2) + 
                1071054*pow(MASS,3)*pow(r,7)*pow(cos(Th),2) + 
                672588*pow(MASS,2)*pow(r,8)*pow(cos(Th),2) - 
                562338*MASS*pow(r,9)*pow(cos(Th),2)))/
            (36750.*pow(MASS,4)*pow(2*MASS - r,4)*pow(r,9)))))/
    pow(r/(-2*MASS + r) + (pow(SPIN,2)*
         (-r - 2*MASS*pow(cos(Th),2) + r*pow(cos(Th),2)))/(r*pow(-2*MASS + r,2))\
       + ZETA*((1840*pow(MASS,5) - 48*pow(MASS,4)*r - 30*pow(MASS,3)*pow(r,2) - 
            260*pow(MASS,2)*pow(r,3) - 15*MASS*pow(r,4) - 15*pow(r,5))/
          (15.*pow(MASS,2)*pow(2*MASS - r,2)*pow(r,5)) - 
         (pow(SPIN,2)*(-299880000*pow(MASS,10) + 583296000*pow(MASS,9)*r - 
              310986200*pow(MASS,8)*pow(r,2) + 55258000*pow(MASS,7)*pow(r,3) - 
              45095130*pow(MASS,6)*pow(r,4) + 29276080*pow(MASS,5)*pow(r,5) - 
              3428093*pow(MASS,4)*pow(r,6) - 1018518*pow(MASS,3)*pow(r,7) + 
              327054*pow(MASS,2)*pow(r,8) + 132321*MASS*pow(r,9) + 55125*pow(r,10) + 
              899640000*pow(MASS,10)*pow(cos(Th),2) - 
              1661688000*pow(MASS,9)*r*pow(cos(Th),2) + 
              1212023400*pow(MASS,8)*pow(r,2)*pow(cos(Th),2) - 
              415644600*pow(MASS,7)*pow(r,3)*pow(cos(Th),2) + 
              127303290*pow(MASS,6)*pow(r,4)*pow(cos(Th),2) - 
              80772240*pow(MASS,5)*pow(r,5)*pow(cos(Th),2) + 
              22301529*pow(MASS,4)*pow(r,6)*pow(cos(Th),2) + 
              1071054*pow(MASS,3)*pow(r,7)*pow(cos(Th),2) + 
              672588*pow(MASS,2)*pow(r,8)*pow(cos(Th),2) - 
              562338*MASS*pow(r,9)*pow(cos(Th),2)))/
          (110250.*pow(MASS,4)*pow(2*MASS - r,3)*pow(r,9))),2) - 
   Jz*(-((Jz*(-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
             ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
                (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                     2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                     2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                     355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                     205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                     24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                     47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                     12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                     12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                     6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                     494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                     562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                     562338*pow(r,8)*pow(cos(Th),2)))/
                 (110250.*pow(MASS,3)*pow(r,11))))*
           (-2*((2*a*MASS*pow(sin(Th),2))/pow(r,2) + 
                (SPIN*(144*pow(MASS,3) + 180*pow(MASS,2)*r + 420*MASS*pow(r,2) + 
                     36*pow(r,3))*ZETA*pow(sin(Th),2))/(15.*MASS*pow(r,7)) - 
                (7*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                     140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
                 (15.*MASS*pow(r,8)))*
              ((-2*a*MASS*pow(sin(Th),2))/r + 
                (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                     140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
                 (15.*MASS*pow(r,7))) + 
             (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
                ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                   (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                   26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                   (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                        2005500*pow(MASS,6)*pow(r,2) - 
                        1538250*pow(MASS,5)*pow(r,3) - 
                        2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                        355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                        205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                        24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                        47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                        12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                        12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                        6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                        494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                        562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                        562338*pow(r,8)*pow(cos(Th),2)))/
                    (110250.*pow(MASS,3)*pow(r,11))))*
              (2*r*pow(sin(Th),2) - 
                (pow(SPIN,2)*(-6213200*pow(MASS,6) - 6833400*pow(MASS,5)*r - 
                     5566950*pow(MASS,4)*pow(r,2) + 3548440*pow(MASS,3)*pow(r,3) + 
                     4003665*pow(MASS,2)*pow(r,4) + 2613240*MASS*pow(r,5) + 
                     1312122*pow(r,6))*ZETA*(-1 + 3*pow(cos(Th),2))*
                   pow(sin(Th),2))/(110250.*pow(MASS,3)*pow(r,8)) + 
                (4*pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                     3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                     887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                     435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                   (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
                 (55125.*pow(MASS,3)*pow(r,9)) - 
                (2*pow(SPIN,2)*MASS*pow(sin(Th),4))/pow(r,2)) + 
             ((-2*MASS)/pow(r,2) + (6*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,4) + 
                ZETA*((-560*pow(MASS,3))/(3.*pow(r,8)) + 
                   (192*pow(MASS,2))/(5.*pow(r,7)) + (22*MASS)/pow(r,6) + 
                   104/(3.*pow(r,5)) + 1/(MASS*pow(r,4)) + 
                   (pow(SPIN,2)*(9153200*pow(MASS,7) + 4011000*pow(MASS,6)*r - 
                        4614750*pow(MASS,5)*pow(r,2) - 
                        11248840*pow(MASS,4)*pow(r,3) - 
                        3939975*pow(MASS,3)*pow(r,4) + 
                        2135844*pow(MASS,2)*pow(r,5) + 3112872*MASS*pow(r,6) + 
                        1646568*pow(r,7) - 24519600*pow(MASS,7)*pow(cos(Th),2) - 
                        94823400*pow(MASS,6)*r*pow(cos(Th),2) + 
                        37790550*pow(MASS,5)*pow(r,2)*pow(cos(Th),2) + 
                        48740520*pow(MASS,4)*pow(r,3)*pow(cos(Th),2) + 
                        34127175*pow(MASS,3)*pow(r,4)*pow(cos(Th),2) - 
                        2967732*pow(MASS,2)*pow(r,5)*pow(cos(Th),2) - 
                        3936366*MASS*pow(r,6)*pow(cos(Th),2) - 
                        4498704*pow(r,7)*pow(cos(Th),2)))/
                    (110250.*pow(MASS,3)*pow(r,11)) - 
                   (11*pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                        2005500*pow(MASS,6)*pow(r,2) - 
                        1538250*pow(MASS,5)*pow(r,3) - 
                        2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                        355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                        205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                        24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                        47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                        12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                        12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                        6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                        494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                        562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                        562338*pow(r,8)*pow(cos(Th),2)))/
                    (110250.*pow(MASS,3)*pow(r,12))))*
              (pow(r,2)*pow(sin(Th),2) - 
                (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                     3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                     887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                     435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                   (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
                 (110250.*pow(MASS,3)*pow(r,8)) + 
                pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))))/
         pow(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7)),2) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)),2)) - 
      (Ee*((-2*a*MASS*pow(sin(Th),2))/r + 
           (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
            (15.*MASS*pow(r,7)))*(-2*
            ((2*a*MASS*pow(sin(Th),2))/pow(r,2) + 
              (SPIN*(144*pow(MASS,3) + 180*pow(MASS,2)*r + 420*MASS*pow(r,2) + 
                   36*pow(r,3))*ZETA*pow(sin(Th),2))/(15.*MASS*pow(r,7)) - 
              (7*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,8)))*
            ((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7))) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (2*r*pow(sin(Th),2) - 
              (pow(SPIN,2)*(-6213200*pow(MASS,6) - 6833400*pow(MASS,5)*r - 
                   5566950*pow(MASS,4)*pow(r,2) + 3548440*pow(MASS,3)*pow(r,3) + 
                   4003665*pow(MASS,2)*pow(r,4) + 2613240*MASS*pow(r,5) + 
                   1312122*pow(r,6))*ZETA*(-1 + 3*pow(cos(Th),2))*
                 pow(sin(Th),2))/(110250.*pow(MASS,3)*pow(r,8)) + 
              (4*pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (55125.*pow(MASS,3)*pow(r,9)) - 
              (2*pow(SPIN,2)*MASS*pow(sin(Th),4))/pow(r,2)) + 
           ((-2*MASS)/pow(r,2) + (6*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,4) + 
              ZETA*((-560*pow(MASS,3))/(3.*pow(r,8)) + 
                 (192*pow(MASS,2))/(5.*pow(r,7)) + (22*MASS)/pow(r,6) + 
                 104/(3.*pow(r,5)) + 1/(MASS*pow(r,4)) + 
                 (pow(SPIN,2)*(9153200*pow(MASS,7) + 4011000*pow(MASS,6)*r - 
                      4614750*pow(MASS,5)*pow(r,2) - 11248840*pow(MASS,4)*pow(r,3) - 
                      3939975*pow(MASS,3)*pow(r,4) + 2135844*pow(MASS,2)*pow(r,5) + 
                      3112872*MASS*pow(r,6) + 1646568*pow(r,7) - 
                      24519600*pow(MASS,7)*pow(cos(Th),2) - 
                      94823400*pow(MASS,6)*r*pow(cos(Th),2) + 
                      37790550*pow(MASS,5)*pow(r,2)*pow(cos(Th),2) + 
                      48740520*pow(MASS,4)*pow(r,3)*pow(cos(Th),2) + 
                      34127175*pow(MASS,3)*pow(r,4)*pow(cos(Th),2) - 
                      2967732*pow(MASS,2)*pow(r,5)*pow(cos(Th),2) - 
                      3936366*MASS*pow(r,6)*pow(cos(Th),2) - 
                      4498704*pow(r,7)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11)) - 
                 (11*pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,12))))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))))/
       pow(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)),2) + 
      (Jz*((-2*MASS)/pow(r,2) + (6*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,4) + 
           ZETA*((-560*pow(MASS,3))/(3.*pow(r,8)) + 
              (192*pow(MASS,2))/(5.*pow(r,7)) + (22*MASS)/pow(r,6) + 
              104/(3.*pow(r,5)) + 1/(MASS*pow(r,4)) + 
              (pow(SPIN,2)*(9153200*pow(MASS,7) + 4011000*pow(MASS,6)*r - 
                   4614750*pow(MASS,5)*pow(r,2) - 11248840*pow(MASS,4)*pow(r,3) - 
                   3939975*pow(MASS,3)*pow(r,4) + 2135844*pow(MASS,2)*pow(r,5) + 
                   3112872*MASS*pow(r,6) + 1646568*pow(r,7) - 
                   24519600*pow(MASS,7)*pow(cos(Th),2) - 
                   94823400*pow(MASS,6)*r*pow(cos(Th),2) + 
                   37790550*pow(MASS,5)*pow(r,2)*pow(cos(Th),2) + 
                   48740520*pow(MASS,4)*pow(r,3)*pow(cos(Th),2) + 
                   34127175*pow(MASS,3)*pow(r,4)*pow(cos(Th),2) - 
                   2967732*pow(MASS,2)*pow(r,5)*pow(cos(Th),2) - 
                   3936366*MASS*pow(r,6)*pow(cos(Th),2) - 
                   4498704*pow(r,7)*pow(cos(Th),2)))/
               (110250.*pow(MASS,3)*pow(r,11)) - 
              (11*pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                   2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                   2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                   355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                   205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                   24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                   47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                   12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                   12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                   6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                   494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                   562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                   562338*pow(r,8)*pow(cos(Th),2)))/
               (110250.*pow(MASS,3)*pow(r,12)))))/
       (-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))) + 
      (Ee*((2*a*MASS*pow(sin(Th),2))/pow(r,2) + 
           (SPIN*(144*pow(MASS,3) + 180*pow(MASS,2)*r + 420*MASS*pow(r,2) + 36*pow(r,3))*
              ZETA*pow(sin(Th),2))/(15.*MASS*pow(r,7)) - 
           (7*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
            (15.*MASS*pow(r,8))))/
       (-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)))) + 
   Ee*((Jz*((-2*a*MASS*pow(sin(Th),2))/r + 
           (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
            (15.*MASS*pow(r,7)))*(-2*
            ((2*a*MASS*pow(sin(Th),2))/pow(r,2) + 
              (SPIN*(144*pow(MASS,3) + 180*pow(MASS,2)*r + 420*MASS*pow(r,2) + 
                   36*pow(r,3))*ZETA*pow(sin(Th),2))/(15.*MASS*pow(r,7)) - 
              (7*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,8)))*
            ((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7))) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (2*r*pow(sin(Th),2) - 
              (pow(SPIN,2)*(-6213200*pow(MASS,6) - 6833400*pow(MASS,5)*r - 
                   5566950*pow(MASS,4)*pow(r,2) + 3548440*pow(MASS,3)*pow(r,3) + 
                   4003665*pow(MASS,2)*pow(r,4) + 2613240*MASS*pow(r,5) + 
                   1312122*pow(r,6))*ZETA*(-1 + 3*pow(cos(Th),2))*
                 pow(sin(Th),2))/(110250.*pow(MASS,3)*pow(r,8)) + 
              (4*pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (55125.*pow(MASS,3)*pow(r,9)) - 
              (2*pow(SPIN,2)*MASS*pow(sin(Th),4))/pow(r,2)) + 
           ((-2*MASS)/pow(r,2) + (6*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,4) + 
              ZETA*((-560*pow(MASS,3))/(3.*pow(r,8)) + 
                 (192*pow(MASS,2))/(5.*pow(r,7)) + (22*MASS)/pow(r,6) + 
                 104/(3.*pow(r,5)) + 1/(MASS*pow(r,4)) + 
                 (pow(SPIN,2)*(9153200*pow(MASS,7) + 4011000*pow(MASS,6)*r - 
                      4614750*pow(MASS,5)*pow(r,2) - 11248840*pow(MASS,4)*pow(r,3) - 
                      3939975*pow(MASS,3)*pow(r,4) + 2135844*pow(MASS,2)*pow(r,5) + 
                      3112872*MASS*pow(r,6) + 1646568*pow(r,7) - 
                      24519600*pow(MASS,7)*pow(cos(Th),2) - 
                      94823400*pow(MASS,6)*r*pow(cos(Th),2) + 
                      37790550*pow(MASS,5)*pow(r,2)*pow(cos(Th),2) + 
                      48740520*pow(MASS,4)*pow(r,3)*pow(cos(Th),2) + 
                      34127175*pow(MASS,3)*pow(r,4)*pow(cos(Th),2) - 
                      2967732*pow(MASS,2)*pow(r,5)*pow(cos(Th),2) - 
                      3936366*MASS*pow(r,6)*pow(cos(Th),2) - 
                      4498704*pow(r,7)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11)) - 
                 (11*pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,12))))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))))/
       pow(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)),2) + 
      (Ee*(pow(r,2)*pow(sin(Th),2) - 
           (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
              (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
            (110250.*pow(MASS,3)*pow(r,8)) + 
           pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))*
         (-2*((2*a*MASS*pow(sin(Th),2))/pow(r,2) + 
              (SPIN*(144*pow(MASS,3) + 180*pow(MASS,2)*r + 420*MASS*pow(r,2) + 
                   36*pow(r,3))*ZETA*pow(sin(Th),2))/(15.*MASS*pow(r,7)) - 
              (7*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,8)))*
            ((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7))) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (2*r*pow(sin(Th),2) - 
              (pow(SPIN,2)*(-6213200*pow(MASS,6) - 6833400*pow(MASS,5)*r - 
                   5566950*pow(MASS,4)*pow(r,2) + 3548440*pow(MASS,3)*pow(r,3) + 
                   4003665*pow(MASS,2)*pow(r,4) + 2613240*MASS*pow(r,5) + 
                   1312122*pow(r,6))*ZETA*(-1 + 3*pow(cos(Th),2))*
                 pow(sin(Th),2))/(110250.*pow(MASS,3)*pow(r,8)) + 
              (4*pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (55125.*pow(MASS,3)*pow(r,9)) - 
              (2*pow(SPIN,2)*MASS*pow(sin(Th),4))/pow(r,2)) + 
           ((-2*MASS)/pow(r,2) + (6*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,4) + 
              ZETA*((-560*pow(MASS,3))/(3.*pow(r,8)) + 
                 (192*pow(MASS,2))/(5.*pow(r,7)) + (22*MASS)/pow(r,6) + 
                 104/(3.*pow(r,5)) + 1/(MASS*pow(r,4)) + 
                 (pow(SPIN,2)*(9153200*pow(MASS,7) + 4011000*pow(MASS,6)*r - 
                      4614750*pow(MASS,5)*pow(r,2) - 11248840*pow(MASS,4)*pow(r,3) - 
                      3939975*pow(MASS,3)*pow(r,4) + 2135844*pow(MASS,2)*pow(r,5) + 
                      3112872*MASS*pow(r,6) + 1646568*pow(r,7) - 
                      24519600*pow(MASS,7)*pow(cos(Th),2) - 
                      94823400*pow(MASS,6)*r*pow(cos(Th),2) + 
                      37790550*pow(MASS,5)*pow(r,2)*pow(cos(Th),2) + 
                      48740520*pow(MASS,4)*pow(r,3)*pow(cos(Th),2) + 
                      34127175*pow(MASS,3)*pow(r,4)*pow(cos(Th),2) - 
                      2967732*pow(MASS,2)*pow(r,5)*pow(cos(Th),2) - 
                      3936366*MASS*pow(r,6)*pow(cos(Th),2) - 
                      4498704*pow(r,7)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11)) - 
                 (11*pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,12))))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))))/
       pow(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)),2) - 
      (Jz*((2*a*MASS*pow(sin(Th),2))/pow(r,2) + 
           (SPIN*(144*pow(MASS,3) + 180*pow(MASS,2)*r + 420*MASS*pow(r,2) + 36*pow(r,3))*
              ZETA*pow(sin(Th),2))/(15.*MASS*pow(r,7)) - 
           (7*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
            (15.*MASS*pow(r,8))))/
       (-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))) - 
      (Ee*(2*r*pow(sin(Th),2) - 
           (pow(SPIN,2)*(-6213200*pow(MASS,6) - 6833400*pow(MASS,5)*r - 
                5566950*pow(MASS,4)*pow(r,2) + 3548440*pow(MASS,3)*pow(r,3) + 
                4003665*pow(MASS,2)*pow(r,4) + 2613240*MASS*pow(r,5) + 
                1312122*pow(r,6))*ZETA*(-1 + 3*pow(cos(Th),2))*
              pow(sin(Th),2))/(110250.*pow(MASS,3)*pow(r,8)) + 
           (4*pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
              (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
            (55125.*pow(MASS,3)*pow(r,9)) - 
           (2*pow(SPIN,2)*MASS*pow(sin(Th),4))/pow(r,2)))/
       (-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))));
P_Th_d = (pow(PTh,2)*(-2*pow(SPIN,2)*cos(Th)*sin(Th) + 
        (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
             3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
             887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
             435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*sin(Th))/
         (18375.*pow(MASS,3)*pow(r,8))))/
    pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2) - 
      (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
           3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
           887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
           435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*(-1 + 3*pow(cos(Th),2)))/
       (110250.*pow(MASS,3)*pow(r,8)),2) + 
   (pow(Pr,2)*((pow(SPIN,2)*(4*MASS*cos(Th)*sin(Th) - 2*r*cos(Th)*sin(Th)))/
         (r*pow(-2*MASS + r,2)) - (pow(SPIN,2)*ZETA*
           (-1799280000*pow(MASS,10)*cos(Th)*sin(Th) + 
             3323376000*pow(MASS,9)*r*cos(Th)*sin(Th) - 
             2424046800*pow(MASS,8)*pow(r,2)*cos(Th)*sin(Th) + 
             831289200*pow(MASS,7)*pow(r,3)*cos(Th)*sin(Th) - 
             254606580*pow(MASS,6)*pow(r,4)*cos(Th)*sin(Th) + 
             161544480*pow(MASS,5)*pow(r,5)*cos(Th)*sin(Th) - 
             44603058*pow(MASS,4)*pow(r,6)*cos(Th)*sin(Th) - 
             2142108*pow(MASS,3)*pow(r,7)*cos(Th)*sin(Th) - 
             1345176*pow(MASS,2)*pow(r,8)*cos(Th)*sin(Th) + 
             1124676*MASS*pow(r,9)*cos(Th)*sin(Th)))/
         (110250.*pow(MASS,4)*pow(2*MASS - r,3)*pow(r,9))))/
    pow(r/(-2*MASS + r) + (pow(SPIN,2)*(-r - 2*MASS*pow(cos(Th),2) + r*pow(cos(Th),2)))/
       (r*pow(-2*MASS + r,2)) + ZETA*
       ((1840*pow(MASS,5) - 48*pow(MASS,4)*r - 30*pow(MASS,3)*pow(r,2) - 
            260*pow(MASS,2)*pow(r,3) - 15*MASS*pow(r,4) - 15*pow(r,5))/
          (15.*pow(MASS,2)*pow(2*MASS - r,2)*pow(r,5)) - 
         (pow(SPIN,2)*(-299880000*pow(MASS,10) + 583296000*pow(MASS,9)*r - 
              310986200*pow(MASS,8)*pow(r,2) + 55258000*pow(MASS,7)*pow(r,3) - 
              45095130*pow(MASS,6)*pow(r,4) + 29276080*pow(MASS,5)*pow(r,5) - 
              3428093*pow(MASS,4)*pow(r,6) - 1018518*pow(MASS,3)*pow(r,7) + 
              327054*pow(MASS,2)*pow(r,8) + 132321*MASS*pow(r,9) + 55125*pow(r,10) + 
              899640000*pow(MASS,10)*pow(cos(Th),2) - 
              1661688000*pow(MASS,9)*r*pow(cos(Th),2) + 
              1212023400*pow(MASS,8)*pow(r,2)*pow(cos(Th),2) - 
              415644600*pow(MASS,7)*pow(r,3)*pow(cos(Th),2) + 
              127303290*pow(MASS,6)*pow(r,4)*pow(cos(Th),2) - 
              80772240*pow(MASS,5)*pow(r,5)*pow(cos(Th),2) + 
              22301529*pow(MASS,4)*pow(r,6)*pow(cos(Th),2) + 
              1071054*pow(MASS,3)*pow(r,7)*pow(cos(Th),2) + 
              672588*pow(MASS,2)*pow(r,8)*pow(cos(Th),2) - 
              562338*MASS*pow(r,9)*pow(cos(Th),2)))/
          (110250.*pow(MASS,4)*pow(2*MASS - r,3)*pow(r,9))),2) - 
   Jz*((Ee*((-4*a*MASS*cos(Th)*sin(Th))/r + 
           (2*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*cos(Th)*sin(Th))/
            (15.*MASS*pow(r,7))))/
       (-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))) + 
      (Jz*((4*pow(SPIN,2)*MASS*cos(Th)*sin(Th))/pow(r,3) + 
           (pow(SPIN,2)*ZETA*(-52920000*pow(MASS,8)*cos(Th)*sin(Th) + 
                49039200*pow(MASS,7)*r*cos(Th)*sin(Th) + 
                94823400*pow(MASS,6)*pow(r,2)*cos(Th)*sin(Th) - 
                25193700*pow(MASS,5)*pow(r,3)*cos(Th)*sin(Th) - 
                24370260*pow(MASS,4)*pow(r,4)*cos(Th)*sin(Th) - 
                13650870*pow(MASS,3)*pow(r,5)*cos(Th)*sin(Th) + 
                989244*pow(MASS,2)*pow(r,6)*cos(Th)*sin(Th) + 
                1124676*MASS*pow(r,7)*cos(Th)*sin(Th) + 
                1124676*pow(r,8)*cos(Th)*sin(Th)))/(110250.*pow(MASS,3)*pow(r,11))))
        /(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))) - 
      (Jz*(-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
           ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
              (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
              (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                   2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                   2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                   355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                   205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                   24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                   47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                   12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                   12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                   6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                   494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                   562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                   562338*pow(r,8)*pow(cos(Th),2)))/
               (110250.*pow(MASS,3)*pow(r,11))))*
         (-2*((-4*a*MASS*cos(Th)*sin(Th))/r + 
              (2*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*cos(Th)*sin(Th))/
               (15.*MASS*pow(r,7)))*
            ((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7))) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (2*pow(r,2)*cos(Th)*sin(Th) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
                 (-1 + 3*pow(cos(Th),2))*sin(Th))/(55125.*pow(MASS,3)*pow(r,8)) + 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
                 pow(sin(Th),3))/(18375.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(2*cos(Th)*sin(Th) + (8*MASS*cos(Th)*pow(sin(Th),3))/r)) + 
           ((4*pow(SPIN,2)*MASS*cos(Th)*sin(Th))/pow(r,3) + 
              (pow(SPIN,2)*ZETA*(-52920000*pow(MASS,8)*cos(Th)*sin(Th) + 
                   49039200*pow(MASS,7)*r*cos(Th)*sin(Th) + 
                   94823400*pow(MASS,6)*pow(r,2)*cos(Th)*sin(Th) - 
                   25193700*pow(MASS,5)*pow(r,3)*cos(Th)*sin(Th) - 
                   24370260*pow(MASS,4)*pow(r,4)*cos(Th)*sin(Th) - 
                   13650870*pow(MASS,3)*pow(r,5)*cos(Th)*sin(Th) + 
                   989244*pow(MASS,2)*pow(r,6)*cos(Th)*sin(Th) + 
                   1124676*MASS*pow(r,7)*cos(Th)*sin(Th) + 
                   1124676*pow(r,8)*cos(Th)*sin(Th)))/
               (110250.*pow(MASS,3)*pow(r,11)))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))))/
       pow(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)),2) - 
      (Ee*((-2*a*MASS*pow(sin(Th),2))/r + 
           (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
            (15.*MASS*pow(r,7)))*(-2*
            ((-4*a*MASS*cos(Th)*sin(Th))/r + 
              (2*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*cos(Th)*sin(Th))/
               (15.*MASS*pow(r,7)))*
            ((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7))) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (2*pow(r,2)*cos(Th)*sin(Th) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
                 (-1 + 3*pow(cos(Th),2))*sin(Th))/(55125.*pow(MASS,3)*pow(r,8)) + 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
                 pow(sin(Th),3))/(18375.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(2*cos(Th)*sin(Th) + (8*MASS*cos(Th)*pow(sin(Th),3))/r)) + 
           ((4*pow(SPIN,2)*MASS*cos(Th)*sin(Th))/pow(r,3) + 
              (pow(SPIN,2)*ZETA*(-52920000*pow(MASS,8)*cos(Th)*sin(Th) + 
                   49039200*pow(MASS,7)*r*cos(Th)*sin(Th) + 
                   94823400*pow(MASS,6)*pow(r,2)*cos(Th)*sin(Th) - 
                   25193700*pow(MASS,5)*pow(r,3)*cos(Th)*sin(Th) - 
                   24370260*pow(MASS,4)*pow(r,4)*cos(Th)*sin(Th) - 
                   13650870*pow(MASS,3)*pow(r,5)*cos(Th)*sin(Th) + 
                   989244*pow(MASS,2)*pow(r,6)*cos(Th)*sin(Th) + 
                   1124676*MASS*pow(r,7)*cos(Th)*sin(Th) + 
                   1124676*pow(r,8)*cos(Th)*sin(Th)))/
               (110250.*pow(MASS,3)*pow(r,11)))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))))/
       pow(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)),2)) + 
   Ee*(-((Jz*((-4*a*MASS*cos(Th)*sin(Th))/r + 
             (2*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                  140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*cos(Th)*sin(Th))/
              (15.*MASS*pow(r,7))))/
         (-pow((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7)),2) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)))) - 
      (Ee*(2*pow(r,2)*cos(Th)*sin(Th) - 
           (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
              (-1 + 3*pow(cos(Th),2))*sin(Th))/(55125.*pow(MASS,3)*pow(r,8)) + 
           (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*pow(sin(Th),3))
             /(18375.*pow(MASS,3)*pow(r,8)) + 
           pow(SPIN,2)*(2*cos(Th)*sin(Th) + (8*MASS*cos(Th)*pow(sin(Th),3))/r)))/
       (-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))) + 
      (Jz*((-2*a*MASS*pow(sin(Th),2))/r + 
           (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
            (15.*MASS*pow(r,7)))*(-2*
            ((-4*a*MASS*cos(Th)*sin(Th))/r + 
              (2*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*cos(Th)*sin(Th))/
               (15.*MASS*pow(r,7)))*
            ((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7))) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (2*pow(r,2)*cos(Th)*sin(Th) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
                 (-1 + 3*pow(cos(Th),2))*sin(Th))/(55125.*pow(MASS,3)*pow(r,8)) + 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
                 pow(sin(Th),3))/(18375.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(2*cos(Th)*sin(Th) + (8*MASS*cos(Th)*pow(sin(Th),3))/r)) + 
           ((4*pow(SPIN,2)*MASS*cos(Th)*sin(Th))/pow(r,3) + 
              (pow(SPIN,2)*ZETA*(-52920000*pow(MASS,8)*cos(Th)*sin(Th) + 
                   49039200*pow(MASS,7)*r*cos(Th)*sin(Th) + 
                   94823400*pow(MASS,6)*pow(r,2)*cos(Th)*sin(Th) - 
                   25193700*pow(MASS,5)*pow(r,3)*cos(Th)*sin(Th) - 
                   24370260*pow(MASS,4)*pow(r,4)*cos(Th)*sin(Th) - 
                   13650870*pow(MASS,3)*pow(r,5)*cos(Th)*sin(Th) + 
                   989244*pow(MASS,2)*pow(r,6)*cos(Th)*sin(Th) + 
                   1124676*MASS*pow(r,7)*cos(Th)*sin(Th) + 
                   1124676*pow(r,8)*cos(Th)*sin(Th)))/
               (110250.*pow(MASS,3)*pow(r,11)))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))))/
       pow(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)),2) + 
      (Ee*(pow(r,2)*pow(sin(Th),2) - 
           (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
              (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
            (110250.*pow(MASS,3)*pow(r,8)) + 
           pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))*
         (-2*((-4*a*MASS*cos(Th)*sin(Th))/r + 
              (2*a*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*cos(Th)*sin(Th))/
               (15.*MASS*pow(r,7)))*
            ((-2*a*MASS*pow(sin(Th),2))/r + 
              (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                   140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
               (15.*MASS*pow(r,7))) + 
           (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
              ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - 
                 (32*pow(MASS,2))/(5.*pow(r,6)) - (22*MASS)/(5.*pow(r,5)) - 
                 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
                 (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                      2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                      2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                      355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                      205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                      24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                      47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                      12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                      12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                      6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                      494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                      562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                      562338*pow(r,8)*pow(cos(Th),2)))/
                  (110250.*pow(MASS,3)*pow(r,11))))*
            (2*pow(r,2)*cos(Th)*sin(Th) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
                 (-1 + 3*pow(cos(Th),2))*sin(Th))/(55125.*pow(MASS,3)*pow(r,8)) + 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*cos(Th)*
                 pow(sin(Th),3))/(18375.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(2*cos(Th)*sin(Th) + (8*MASS*cos(Th)*pow(sin(Th),3))/r)) + 
           ((4*pow(SPIN,2)*MASS*cos(Th)*sin(Th))/pow(r,3) + 
              (pow(SPIN,2)*ZETA*(-52920000*pow(MASS,8)*cos(Th)*sin(Th) + 
                   49039200*pow(MASS,7)*r*cos(Th)*sin(Th) + 
                   94823400*pow(MASS,6)*pow(r,2)*cos(Th)*sin(Th) - 
                   25193700*pow(MASS,5)*pow(r,3)*cos(Th)*sin(Th) - 
                   24370260*pow(MASS,4)*pow(r,4)*cos(Th)*sin(Th) - 
                   13650870*pow(MASS,3)*pow(r,5)*cos(Th)*sin(Th) + 
                   989244*pow(MASS,2)*pow(r,6)*cos(Th)*sin(Th) + 
                   1124676*MASS*pow(r,7)*cos(Th)*sin(Th) + 
                   1124676*pow(r,8)*cos(Th)*sin(Th)))/
               (110250.*pow(MASS,3)*pow(r,11)))*
            (pow(r,2)*pow(sin(Th),2) - 
              (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                   3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                   887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                   435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
                 (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (110250.*pow(MASS,3)*pow(r,8)) + 
              pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r))))/
       pow(-pow((-2*a*MASS*pow(sin(Th),2))/r + 
            (SPIN*(-400*pow(MASS,4) + 144*pow(MASS,3)*r + 90*pow(MASS,2)*pow(r,2) + 
                 140*MASS*pow(r,3) + 9*pow(r,4))*ZETA*pow(sin(Th),2))/
             (15.*MASS*pow(r,7)),2) + 
         (-1 + (2*MASS)/r - (2*pow(SPIN,2)*MASS*pow(cos(Th),2))/pow(r,3) + 
            ZETA*((80*pow(MASS,3))/(3.*pow(r,7)) - (32*pow(MASS,2))/(5.*pow(r,6)) - 
               (22*MASS)/(5.*pow(r,5)) - 26/(3.*pow(r,4)) - 1/(3.*MASS*pow(r,3)) + 
               (pow(SPIN,2)*(-8820000*pow(MASS,8) + 9153200*pow(MASS,7)*r + 
                    2005500*pow(MASS,6)*pow(r,2) - 1538250*pow(MASS,5)*pow(r,3) - 
                    2812210*pow(MASS,4)*pow(r,4) - 787995*pow(MASS,3)*pow(r,5) + 
                    355974*pow(MASS,2)*pow(r,6) + 444696*MASS*pow(r,7) + 
                    205821*pow(r,8) + 26460000*pow(MASS,8)*pow(cos(Th),2) - 
                    24519600*pow(MASS,7)*r*pow(cos(Th),2) - 
                    47411700*pow(MASS,6)*pow(r,2)*pow(cos(Th),2) + 
                    12596850*pow(MASS,5)*pow(r,3)*pow(cos(Th),2) + 
                    12185130*pow(MASS,4)*pow(r,4)*pow(cos(Th),2) + 
                    6825435*pow(MASS,3)*pow(r,5)*pow(cos(Th),2) - 
                    494622*pow(MASS,2)*pow(r,6)*pow(cos(Th),2) - 
                    562338*MASS*pow(r,7)*pow(cos(Th),2) - 
                    562338*pow(r,8)*pow(cos(Th),2)))/
                (110250.*pow(MASS,3)*pow(r,11))))*
          (pow(r,2)*pow(sin(Th),2) - 
            (pow(SPIN,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                 887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                 435540*MASS*pow(r,6) + 187446*pow(r,7))*ZETA*
               (-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (110250.*pow(MASS,3)*pow(r,8)) + 
            pow(SPIN,2)*(pow(sin(Th),2) + (2*MASS*pow(sin(Th),4))/r)),2));
P_Ph_d = 0;

gsl_vector_set(out_state, 0, r_d/2);
gsl_vector_set(out_state, 1, P_r_d/2);
gsl_vector_set(out_state, 2, Th_d/2);
gsl_vector_set(out_state, 3, P_Th_d/2);
gsl_vector_set(out_state, 4, Ph_d/2);
gsl_vector_set(out_state, 5, P_Ph_d/2);
return 0;
}
