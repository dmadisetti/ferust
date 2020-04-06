function [P] = local_projection_4(X, Y, nodes, u, l, local_nodes)
  %local_stiffness Calculates the local projection vector for a given element
  % and its type, using an appropriate integration method.
  % NOTE: THIS IS GENERATED CODE. REFER TO local_projection.m.tmpl and generate.py
  if local_nodes == 4
    [P] = crunch_quad(X(1), X(2), X(3), X(4), ...
                      Y(1), Y(2), Y(3), Y(4), ...
                      nodes(1, 1), nodes(2, 1), nodes(3, 1), nodes(4, 1), ...
                      nodes(1, 2), nodes(2, 2), nodes(3, 2), nodes(4, 2), ...
                      u, l);
  elseif local_nodes == 8
    [P] = crunch_serindipity(X(1), X(2), X(3), X(4), ...
                           X(5), X(6), X(7), X(8), ...
                           Y(1), Y(2), Y(3), Y(4), ...
                           Y(5), Y(6), Y(7), Y(8), ...
                           nodes(1, 1), nodes(2, 1), nodes(3, 1), nodes(4, 1), ...
                           nodes(5, 1), nodes(6, 1), nodes(7, 1), nodes(8, 1), ...
                           nodes(1, 2), nodes(2, 2), nodes(3, 2), nodes(4, 2), ...
                           nodes(5, 2), nodes(6, 2), nodes(7, 2), nodes(8, 2), ...
                           u, l);
  else
    [P] = crunch_bubble(X(1), X(2), X(3), X(4), ...
                           X(5), X(6), X(7), X(8), X(9), ...
                           Y(1), Y(2), Y(3), Y(4), X(5), ...
                           Y(6), Y(7), Y(8), Y(9), ...
                           nodes(1, 1), nodes(2, 1), nodes(3, 1), nodes(4, 1), ...
                           nodes(5, 1), nodes(6, 1), nodes(7, 1), nodes(8, 1), nodes(9, 1), ...
                           nodes(1, 2), nodes(2, 2), nodes(3, 2), nodes(4, 2), ...
                           nodes(5, 2), nodes(6, 2), nodes(7, 2), nodes(8, 2), nodes(9, 2), ...
                           u, l);
  end
end

function [P] = crunch_quad(A, B, C, D, ...
                           R, S, T, U, ...
                           d0X, d1X, d2X, d3X, ...
                           d0Y, d1Y, d2Y, d3Y, ...
                           u, l)
  
    x0=sqrt(3);x1=x0 + 3;x2=D*x1;x3=3 - x0;x4=B*x3 ...
    ;x5=-A*x1 + C*x3;x6=x2 - x4 + x5;x7=x3*x6;x8=B*x1 ...
    ;x9=D*x3;x10=x5 + x8 - x9;x11=x1*x10;x12=x11 + x7;x13=d3Y* ...
    l;x14=l + 2*u;x15=U*x1;x16=S*x3;x17=-R*x1 + T*x3; ...
    x18=x15 - x16 + x17;x19=x1*x18;x20=-x19;x21=U*x3;x22=S*x1;x23= ...
    x17 - x21 + x22;x24=x1*x23;x25=x20 + x24;x26=d0X*x25;x27=x23*x3;x28 ...
    =x19 + x27;x29=d1X*x28;x30=x1*x6;x31=-x11;x32=x30 + x31;x33=d0Y*x32 ...
    ;x34=-x7;x35=x10*x3;x36=x34 + x35;x37=d2Y*l;x38=-x30;x39=-x35;x40 ...
    =x38 + x39;x41=d1Y*l;x42=x18*x3;x43=-x27;x44=x42 + x43;x45=d2X*x14 ...
    ;x46=-x42;x47=-x24;x48=x46 + x47;x49=d3X*x48;x50=l*x33 + x12*x13  ...
    + x14*x26 + x14*x29 + x14*x49 + x36*x37 + x40*x41 + x44*x45;x51 ...
    =x1^2;x52=-A*x3 + C*x1;x53=x52 - x8 + x9;x54=x3*x53;x55 ...
    =-x2 + x4 + x52;x56=x3*x55;x57=-x56;x58=x54 + x57;x59=d0Y*l; ...
    x60=-R*x3 + T*x1;x61=x21 - x22 + x60;x62=x3*x61;x63=-x15 + ...
     x16 + x60;x64=x1*x63;x65=x62 + x64;x66=d1X*x65;x67=-x54;x68=x1*x55 ...
    ;x69=-x68;x70=x67 + x69;x71=x1*x53;x72=x56 + x71;x73=-x71;x74=x68  ...
    + x73;x75=x1*x61;x76=-x64;x77=x75 + x76;x78=-x62;x79=x3*x63;x80=x78 ...
     + x79;x81=d0X*x80;x82=-x75;x83=-x79;x84=x82 + x83;x85=d3X*x84;x86= ...
    x13*x72 + x14*x66 + x14*x81 + x14*x85 + x37*x74 + x41*x70 + x45 ...
    *x77 + x58*x59;x87=x3^2;x88=x46 + x64;x89=d0X*x88;x90=x69 + x7; ...
    x91=x34 + x57;x92=x42 + x79;x93=d1X*x92;x94=x30 + x68;x95=x38 + x56 ...
    ;x96=x20 + x76;x97=d3X*x96;x98=x19 + x83;x99=x13*x94 + x14*x89 +  ...
    x14*x93 + x14*x97 + x37*x95 + x41*x91 + x45*x98 + x59*x90;x100=x1 ...
    *x3;x101=x27 + x82;x102=d0X*x101;x103=x39 + x71;x104=x31 + x73;x105=x24  ...
    + x75;x106=d1X*x105;x107=x35 + x54;x108=x11 + x67;x109=x47 + x62;x110=x43 ...
     + x78;x111=d3X*x110;x112=x102*x14 + x103*x59 + x104*x41 + x106*x14 +  ...
    x107*x13 + x108*x37 + x109*x45 + x111*x14;x113=x100*x112 + x100*x99;x114=d1Y ...
    *x14;x115=d2X*l;x116=d3Y*x14;x117=d2Y*x14;x118=l*x26 + l*x29 + l* ...
    x49 + x114*x40 + x115*x44 + x116*x12 + x117*x36 + x14*x33;x119=d0Y*x14 ...
    ;x120=l*x66 + l*x81 + l*x85 + x114*x70 + x115*x77 + x116*x72  ...
    + x117*x74 + x119*x58;x121=l*x102 + l*x106 + l*x111 + x103*x119 + ...
     x104*x114 + x107*x116 + x108*x117 + x109*x115;x122=l*x89 + l*x93 +  ...
    l*x97 + x114*x91 + x115*x98 + x116*x94 + x117*x95 + x119*x90;x123=x100 ...
    *x121 + x100*x122;x124=d0Y*u;x125=d0X*u;x126=d2Y*u;x127=d1Y*u;x128=d3X* ...
    u;x129=d2X*u;x130=d1X*u;x131=d3Y*u;x132=x12*x128 + x124*x25 + x125*x32 ...
     + x126*x44 + x127*x28 + x129*x36 + x130*x40 + x131*x48;x133=x124*x80  ...
    + x125*x58 + x126*x77 + x127*x65 + x128*x72 + x129*x74 + x130*x70 + ...
     x131*x84;x134=x124*x88 + x125*x90 + x126*x98 + x127*x92 + x128*x94 +  ...
    x129*x95 + x130*x91 + x131*x96;x135=x101*x124 + x103*x125 + x104*x130 + x105 ...
    *x127 + x107*x128 + x108*x129 + x109*x126 + x110*x131;x136=x100*x134 + x100* ...
    x135;x137=x100*x50 + x100*x86;x138=x100*x118 + x100*x120;x139=x100*x132 + x100*x133;
    P=[x113/5184 + x50*x51/5184 + x86*x87/5184, x118*x51/5184 + x120*x87/ ...
    5184 + x123/5184, x132*x51/5184 + x133*x87/5184 + x136/5184, x112*x51/5184 + ...
     x137/5184 + x87*x99/5184, x121*x51/5184 + x122*x87/5184 + x138/5184, x134* ...
    x87/5184 + x135*x51/5184 + x139/5184, x113/5184 + x50*x87/5184 + x51*x86/ ...
    5184, x118*x87/5184 + x120*x51/5184 + x123/5184, x132*x87/5184 + x133*x51/5184 ...
     + x136/5184, x112*x87/5184 + x137/5184 + x51*x99/5184, x121*x87/5184 +  ...
    x122*x51/5184 + x138/5184, x134*x51/5184 + x135*x87/5184 + x139/5184];
end

function [P] = crunch_serindipity(A, B, C, D, E, F, G, H, ...
                                  R, S, T, U, V, W, X, Y, ...
                                  d0X, d1X, d2X, d3X, d4X, d5X, d6X, d7X, ...
                                  d0Y, d1Y, d2Y, d3Y, d4Y, d5Y, d6Y, d7Y, ...
                                  u, l)
  
    x0=sqrt(3);x1=x0 + 3;x2=x1^2 - 4*x1;x3=3 - x0 ...
    ;x4=x0*x3;x5=F*x4;x6=x0*x1;x7=H*x6;x8=-3*x3;x9=-2*x4  ...
    + x8 + 6;x10=2*x6;x11=3*x0;x12=x11 + 3;x13=-x10 + x12; ...
    x14=-3*x1;x15=x14 - 2*x6 + 6;x16=2*x4;x17=3 - x11;x18=- ...
    x16 + x17;x19=A*x15 + C*x18;x20=12*G;x21=-12*E + x20;x22=B ...
    *x9 + D*x13 + x19 + x21 + 4*x5 + 4*x7;x23=E*x6;x24= ...
    3*F - 3*H;x25=4*G;x26=x25*x4;x27=B*x13 + D*x9 + x19 ...
     + 4*x23 + 4*x24 + x26;x28=-x15*x22 + x15*x27;x29=d0Y*x28;x30 ...
    =3*x27;x31=-x22*x4 + x30;x32=4*l;x33=d6Y*x32;x34=-x13*x22 +  ...
    x27*x9;x35=d1Y*x34;x36=Y*x6;x37=W*x4;x38=3*X;x39=-3*V + x38 ...
    ;x40=R*x15 + T*x18;x41=S*x9 + U*x13 + 4*x36 + 4*x37  ...
    + 4*x39 + x40;x42=V*x6;x43=X*x4;x44=-3*Y;x45=3*W;x46=x44 ...
     + x45;x47=S*x13 + U*x9 + x40 + 4*x42 + 4*x43 + 4* ...
    x46;x48=3*x47;x49=x4*x41 - x48;x50=l + 2*u;x51=4*x50;x52=d6X ...
    *x51;x53=x13*x41 - x47*x9;x54=d1X*x53;x55=x15*x41 - x15*x47;x56=d0X* ...
    x50;x57=3*x22;x58=x27*x4 - x57;x59=d5Y*x32;x60=x18*x41 - x18*x47;x61 ...
    =d2X*x60;x62=3*x41;x63=-x4*x47 + x62;x64=d5X*x63;x65=-x22*x6 -  ...
    x30;x66=d4Y*x65;x67=x13*x27 - x22*x9;x68=d3Y*x67;x69=x41*x6 + x48;x70 ...
    =d4X*x51;x71=-x18*x22 + x18*x27;x72=d2Y*x71;x73=-x13*x47 + x41*x9; ...
    x74=d3X*x73;x75=-x47*x6 - x62;x76=d7X*x75;x77=x27*x6 + x57;x78=d7Y* ...
    x32;x79=l*x29 + l*x35 + l*x68 + l*x72 + x31*x33 + x32*x66 ...
     + x49*x52 + x50*x54 + x50*x61 + x50*x74 + x51*x64 + x51*x76  ...
    + x55*x56 + x58*x59 + x69*x70 + x77*x78;x80=x3^2 - 4*x3;x81 ...
    =F*x6;x82=H*x4;x83=x16 + x17;x84=x14 + 2*x6 + 6;x85=x10  ...
    + x12;x86=2*x4 + x8 + 6;x87=A*x86 + C*x85;x88=B*x84 + ...
     D*x83 - 12*E + x20 - 4*x81 - 4*x82 + x87;x89=-12*H; ...
    x90=x25*x6;x91=E*x4;x92=B*x83 + D*x84 + 12*F + x87 + x89 ...
     - x90 - 4*x91;x93=3*x92;x94=x6*x88 + x93;x95=x4*x88 - x93; ...
    x96=d4Y*x32;x97=x83*x92 - x84*x88;x98=d3Y*x97;x99=W*x6;x100=Y*x4;x101 ...
    =-12*V;x102=R*x86 + T*x85;x103=S*x84 + U*x83 - 4*x100 +  ...
    x101 + x102 + 4*x38 - 4*x99;x104=3*x103;x105=V*x4;x106=X*x6;x107 ...
    =S*x83 + U*x84 + x102 - 4*x105 - 4*x106 + 4*x46;x108=x104  ...
    + x107*x6;x109=d5X*x108;x110=3*x88;x111=x110 - x4*x92;x112=-x110 - x6* ...
    x92;x113=x103*x85 - x107*x85;x114=d2X*x113;x115=-x104 + x107*x4;x116=d7X*x115; ...
    x117=-x85*x88 + x85*x92;x118=d2Y*x117;x119=x103*x84 - x107*x83;x120=d3X*x119; ...
    x121=-x86*x88 + x86*x92;x122=d0Y*x121;x123=3*x107;x124=-x103*x6 - x123;x125 ...
    =x103*x83 - x107*x84;x126=d1X*x125;x127=-x83*x88 + x84*x92;x128=d1Y*x127;x129 ...
    =-x103*x4 + x123;x130=x103*x86 - x107*x86;x131=d0X*x130;x132=l*x118 + l ...
    *x122 + l*x128 + l*x98 + x109*x51 + x111*x78 + x112*x59 + x114* ...
    x50 + x116*x51 + x120*x50 + x124*x52 + x126*x50 + x129*x70 + x131*x50 ...
     + x33*x94 + x95*x96;x133=A*x84 + C*x83;x134=B*x86 + D*x85  ...
    + x133 + x21 - 4*x5 - 4*x7;x135=A*x9 + C*x13;x136=B*x18 ...
     + D*x15 + 12*F + x135 + x89 + x90 + 4*x91;x137=3*x136; ...
    x138=-x134*x6 + x137;x139=R*x84 + T*x83;x140=S*x86 + U*x85 +  ...
    x139 - 4*x36 - 4*x37 + 4*x39;x141=R*x9 + T*x13;x142=S*x18 ...
     + U*x15 + 12*W + 4*x105 + 4*x106 + x141 + 4*x44;x143= ...
    x13*x140 - x142*x83;x144=d2X*x143;x145=3*x134;x146=-x136*x4 - x145;x147=-x136 ...
    *x6 + x145;x148=3*x142;x149=x140*x6 - x148;x150=d6X*x149;x151=x140*x18 -  ...
    x142*x86;x152=d1X*x151;x153=x140*x4 + x148;x154=d4X*x153;x155=3*x140;x156=x142*x6 ...
     - x155;x157=d7X*x156;x158=-x134*x18 + x136*x86;x159=d1Y*l;x160=x140*x9 - ...
     x142*x84;x161=-x134*x15 + x136*x85;x162=d3Y*x161;x163=-x13*x134 + x136*x83; ...
    x164=d2Y*x163;x165=-x134*x9 + x136*x84;x166=d0Y*x165;x167=x142*x4 + x155;x168= ...
    d5X*x167;x169=x140*x15 - x142*x85;x170=d3X*x50;x171=-x134*x4 - x137;x172=d4Y* ...
    x171;x173=l*x162 + l*x164 + l*x166 + x138*x33 + x144*x50 + x146*x59 ...
     + x147*x78 + x150*x51 + x152*x50 + x154*x51 + x157*x51 + x158*x159  ...
    + x160*x56 + x168*x51 + x169*x170 + x172*x32;x174=x1*x3 - 12;x175=S ...
    *x15 + U*x18 + 4*x100 + x101 + x141 + 4*x38 + 4*x99;x176= ...
    S*x85 + U*x86 - 12*Y + x139 - 4*x42 - 4*x43 + 4*x45 ...
    ;x177=3*x176;x178=-x175*x6 + x177;x179=-x15*x176 + x175*x85;x180=d1X*x179; ...
    x181=-x13*x176 + x175*x83;x182=d2X*x181;x183=B*x15 + D*x18 + x135 +  ...
    x21 + 4*x81 + 4*x82;x184=B*x85 + D*x86 + x133 - 4*x23 + ...
     4*x24 - x26;x185=3*x184;x186=x183*x4 + x185;x187=-x183*x84 + x184*x9 ...
    ;x188=d0Y*x187;x189=3*x175;x190=-x176*x4 - x189;x191=d7X*x190;x192=-x176*x6  ...
    + x189;x193=d5X*x192;x194=x175*x84 - x176*x9;x195=x175*x86 - x176*x18;x196=x15 ...
    *x184 - x183*x85;x197=x183*x6 - x185;x198=x13*x184 - x183*x83;x199=d2Y*x198; ...
    x200=3*x183;x201=x184*x4 + x200;x202=x184*x6 - x200;x203=x18*x184 - x183*x86 ...
    ;x204=d3Y*x203;x205=-x175*x4 - x177;x206=l*x188 + l*x199 + l*x204 + ...
     x159*x196 + x170*x195 + x178*x70 + x180*x50 + x182*x50 + x186*x33 +  ...
    x191*x51 + x193*x51 + x194*x56 + x197*x96 + x201*x78 + x202*x59 + x205 ...
    *x52;x207=x173*x174 + x174*x206;x208=d6Y*x51;x209=d0X*l;x210=d5Y*x51;x211=d4X* ...
    x32;x212=d7Y*x51;x213=d6X*x32;x214=l*x54 + l*x61 + l*x74 + x208*x31 ...
     + x209*x55 + x210*x58 + x211*x69 + x212*x77 + x213*x49 + x29*x50  ...
    + x32*x64 + x32*x76 + x35*x50 + x50*x68 + x50*x72 + x51*x66;x215 ...
    =d4Y*x51;x216=l*x114 + l*x120 + l*x126 + l*x131 + x109*x32 +  ...
    x111*x212 + x112*x210 + x116*x32 + x118*x50 + x122*x50 + x124*x213 + x128 ...
    *x50 + x129*x211 + x208*x94 + x215*x95 + x50*x98;x217=d1Y*x50;x218=d3X* ...
    l;x219=l*x180 + l*x182 + x178*x211 + x186*x208 + x188*x50 + x191*x32 ...
     + x193*x32 + x194*x209 + x195*x218 + x196*x217 + x197*x215 + x199*x50  ...
    + x201*x212 + x202*x210 + x204*x50 + x205*x213;x220=l*x144 + l*x152 + ...
     4*l*x154 + x138*x208 + x146*x210 + x147*x212 + x150*x32 + x157*x32  ...
    + x158*x217 + x160*x209 + x162*x50 + x164*x50 + x166*x50 + x168*x32 + ...
     x169*x218 + 4*x172*x50;x221=x174*x219 + x174*x220;x222=4*u;x223=d4Y*x222; ...
    x224=d5Y*x222;x225=d1Y*u;x226=d2X*u;x227=d7Y*x222;x228=d7X*x222;x229=d3X*u;x230 ...
    =d0X*u;x231=d3Y*u;x232=d6X*x222;x233=d2Y*u;x234=d5X*x222;x235=d1X*u;x236= ...
    d0Y*u;x237=d6Y*x222;x238=d4X*x222;x239=x223*x69 + x224*x63 + x225*x53 + x226 ...
    *x71 + x227*x75 + x228*x77 + x229*x67 + x230*x28 + x231*x73 + x232* ...
    x31 + x233*x60 + x234*x58 + x235*x34 + x236*x55 + x237*x49 + x238*x65 ...
    ;x240=x108*x224 + x111*x228 + x112*x234 + x113*x233 + x115*x227 + x117*x226  ...
    + x119*x231 + x121*x230 + x124*x237 + x125*x225 + x127*x235 + x129*x223 + ...
     x130*x236 + x229*x97 + x232*x94 + x238*x95;x241=x138*x232 + x143*x233 +  ...
    x146*x234 + x147*x228 + x149*x237 + x151*x225 + x153*x223 + x156*x227 + x158 ...
    *x235 + x160*x236 + x161*x229 + x163*x226 + x165*x230 + x167*x224 + x169* ...
    x231 + x171*x238;x242=x178*x223 + x179*x225 + x181*x233 + x186*x232 + x187*x230 ...
     + x190*x227 + x192*x224 + x194*x236 + x195*x231 + x196*x235 + x197*x238  ...
    + x198*x226 + x201*x228 + x202*x234 + x203*x229 + x205*x237;x243=x174*x241 + ...
     x174*x242;x244=x132*x174 + x174*x79;x245=x174*x214 + x174*x216;x246=x174*x239 +  ...
    x174*x240;x247=x1*x206 + x173*x3;x248=x1*x79 + x132*x3;x249=x1*x219 + x220 ...
    *x3;x250=x1*x214 + x216*x3;x251=x1*x242 + x241*x3;x252=x1*x239 + x240* ...
    x3;x253=x1*x132 + x3*x79;x254=x1*x216 + x214*x3;x255=x1*x240 + x239*x3 ...
    ;x256=x1*x173 + x206*x3;x257=x1*x220 + x219*x3;x258=x1*x241 + x242*x3;
    P=[x132*x80/46656 + x2*x79/46656 + x207/46656, x2*x214/46656 + x216*x80/ ...
    46656 + x221/46656, x2*x239/46656 + x240*x80/46656 + x243/46656, x173*x80/46656 + ...
     x2*x206/46656 + x244/46656, x2*x219/46656 + x220*x80/46656 + x245/46656, x2* ...
    x242/46656 + x241*x80/46656 + x246/46656, x132*x2/46656 + x207/46656 + x79*x80/ ...
    46656, x2*x216/46656 + x214*x80/46656 + x221/46656, x2*x240/46656 + x239*x80/46656 ...
     + x243/46656, x173*x2/46656 + x206*x80/46656 + x244/46656, x2*x220/46656 +  ...
    x219*x80/46656 + x245/46656, x2*x241/46656 + x242*x80/46656 + x246/46656, x247/11664 ...
     + x248/11664, x249/11664 + x250/11664, x251/11664 + x252/11664, x247/11664 + x253 ...
    /11664, x249/11664 + x254/11664, x251/11664 + x255/11664, x253/11664 + x256/11664,  ...
    x254/11664 + x257/11664, x255/11664 + x258/11664, x248/11664 + x256/11664, x250/11664  ...
    + x257/11664, x252/11664 + x258/11664];
end

function [P] = crunch_bubble(A, B, C, D, E, F, G, H, I, ...
                             R, S, T, U, V, W, X, Y, Z, ...
                             d0X, d1X, d2X, d3X, d4X, d5X, d6X, d7X, d8X, ...
                             d0Y, d1Y, d2Y, d3Y, d4Y, d5Y, d6Y, d7Y, d8Y, ...
                             u, l)
  
    x0=sqrt(3);x1=x0 + 3;x2=x1^2 - 4*x1 + 4;x3=16 ...
    *I*x0;x4=2*x0;x5=-x4;x6=-x0;x7=x6 + 3;x8=x0*x7;x9= ...
    x5 + x8;x10=4*x9;x11=x0*x1;x12=2*x11;x13=7*x0;x14=-x12 +  ...
    x13 + 3;x15=x4 - x8;x16=-3*x7;x17=2*x15 + x16 + 6;x18=x11 ...
     + x5;x19=4*H;x20=x5 + 3;x21=G*x20;x22=x5 - 3;x23=E* ...
    x22;x24=2*x8;x25=x1 - x24;x26=-x11 + x4;x27=-3*x1;x28=2*x26  ...
    + x27 + 6;x29=A*x28 + C*x25;x30=B*x17 + D*x14 + F*x10 ...
     + x18*x19 + 4*x21 + 4*x23 + x29 + x3;x31=H*x22;x32=4* ...
    E;x33=4*F;x34=x20*x33 + x3;x35=B*x14 + D*x17 + G*x10 + ...
     x18*x32 + x29 + 4*x31 + x34;x36=-x0*x30 + x0*x35;x37=16*d8Y ...
    ;x38=x36*x37;x39=-x28*x30 + x28*x35;x40=d0Y*x39;x41=Z*x0;x42=4*x41 ...
    ;x43=4*Y;x44=R*x28 + T*x25;x45=4*X;x46=4*V;x47=x20*x45  ...
    + x22*x46;x48=S*x17 + U*x14 + 4*W*x9 + x18*x43 + 4*x42 ...
     + x44 + x47;x49=16*x41;x50=4*W;x51=x20*x50 + x22*x43 + x49; ...
    x52=S*x14 + U*x17 + X*x10 + x18*x46 + x44 + x51;x53=x0*x48 ...
     - x0*x52;x54=2*u;x55=l + x54;x56=16*d8X;x57=x55*x56;x58=x20* ...
    x48 - x52*x9;x59=4*x55;x60=d5X*x59;x61=4*l;x62=x18*x35 - x22*x30 ...
    ;x63=d7Y*x62;x64=-x14*x52 + x17*x48;x65=d3X*x55;x66=-x18*x52 + x22* ...
    x48;x67=d7X*x66;x68=x14*x48 - x17*x52;x69=d1X*x68;x70=-x25*x30 + x25* ...
    x35;x71=d2Y*l;x72=-x18*x30 + x22*x35;x73=d4Y*x72;x74=x28*x48 - x28* ...
    x52;x75=d0X*x55;x76=x18*x48 - x22*x52;x77=d4X*x59;x78=-x20*x30 + x35* ...
    x9;x79=d5Y*x78;x80=x25*x48 - x25*x52;x81=d2X*x80;x82=-x14*x30 + x17* ...
    x35;x83=d1Y*x82;x84=-x20*x52 + x48*x9;x85=d6X*x84;x86=x20*x35 - x30* ...
    x9;x87=d6Y*x61;x88=x14*x35 - x17*x30;x89=d3Y*l;x90=l*x38 + l*x40 ...
     + l*x83 + x53*x57 + x55*x69 + x55*x81 + x58*x60 + x59*x67  ...
    + x59*x85 + x61*x63 + x61*x73 + x61*x79 + x64*x65 + x70*x71 + ...
     x74*x75 + x76*x77 + x86*x87 + x88*x89;x91=x7^2 - 4*x7 +  ...
    4;x92=x4 + 3;x93=G*x92;x94=x4 - 3;x95=E*x94;x96=2*x18 + ...
     x27 + 6;x97=-x13 + x24 + 3;x98=x12 + x7;x99=-x3;x100=x16  ...
    + 2*x9 + 6;x101=A*x100 + C*x98 + x99;x102=B*x96 + D*x97 ...
     + 4*F*x26 + 4*H*x15 + x101 + 4*x93 + 4*x95;x103=H* ...
    x94;x104=x33*x92;x105=B*x97 + D*x96 + 4*E*x15 + 4*G*x26 + ...
     x101 + 4*x103 + x104;x106=x0*x102 - x0*x105;x107=l*x37;x108=-x102*x26 ...
     + x105*x92;x109=-x49;x110=x45*x92 + x46*x94;x111=R*x100 + T*x98;x112 ...
    =S*x96 + U*x97 + x109 + x110 + x111 + x15*x43 + x26*x50;x113= ...
    x109 + x43*x94 + x50*x92;x114=S*x97 + U*x96 + x111 + x113 + x15 ...
    *x46 + x26*x45;x115=-x0*x112 + x0*x114;x116=x112*x15 - x114*x94;x117=d4X ...
    *x116;x118=x100*x112 - x100*x114;x119=d0X*x118;x120=x112*x26 - x114*x92;x121=d6X* ...
    x120;x122=-x102*x98 + x105*x98;x123=x112*x97 - x114*x96;x124=d1X*x123;x125=-x100 ...
    *x102 + x100*x105;x126=d0Y*x125;x127=x112*x96 - x114*x97;x128=-x102*x97 + x105 ...
    *x96;x129=d1Y*l;x130=x112*x92 - x114*x26;x131=-x102*x15 + x105*x94;x132=d4Y ...
    *x131;x133=-x102*x94 + x105*x15;x134=d7Y*x133;x135=-x102*x92 + x105*x26;x136= ...
    d5Y*x135;x137=-x102*x96 + x105*x97;x138=x112*x98 - x114*x98;x139=d2X*x138;x140= ...
    x112*x94 - x114*x15;x141=d7X*x140;x142=l*x126 + x106*x107 + x108*x87 + x115 ...
    *x57 + x117*x59 + x119*x55 + x121*x59 + x122*x71 + x124*x55 + x127* ...
    x65 + x128*x129 + x130*x60 + x132*x61 + x134*x61 + x136*x61 + x137*x89 ...
     + x139*x55 + x141*x59;x143=x1*x7 - 8;x144=4*x15;x145=A*x96 +  ...
    C*x97 + x99;x146=B*x100 + D*x98 + F*x144 + x145 + x19*x26 + ...
     4*x93 + 4*x95;x147=A*x17 + C*x14;x148=B*x25 + D*x28 +  ...
    4*E*x9 + 4*G*x18 + x147 + 4*x31 + x34;x149=-x146*x22 +  ...
    x148*x26;x150=d7Y*x149;x151=R*x96 + T*x97;x152=S*x100 + U*x98 + 4 ...
    *W*x15 + x110 + x151 + x26*x43 - 4*x42;x153=R*x17 + T*x14; ...
    x154=S*x25 + U*x28 + V*x10 + x153 + x18*x45 + x51;x155=x14*x152 ...
     - x154*x97;x156=d2X*x155;x157=x152*x9 - x154*x94;x158=x0*x152 + x0*x154; ...
    x159=-x146*x17 + x148*x96;x160=d0Y*x159;x161=-x100*x154 + x152*x25;x162=d1X*x161 ...
    ;x163=x152*x22 - x154*x26;x164=d7X*x163;x165=x152*x18 - x154*x92;x166=d6X*x165; ...
    x167=-x14*x146 + x148*x97;x168=x100*x148 - x146*x25;x169=-x146*x9 + x148*x94 ...
    ;x170=d4Y*x169;x171=x152*x17 - x154*x96;x172=-x0*x146 - x0*x148;x173=-x146* ...
    x20 + x148*x15;x174=d5Y*x173;x175=-x15*x154 + x152*x20;x176=x152*x28 - x154* ...
    x98;x177=-x146*x18 + x148*x92;x178=-x146*x28 + x148*x98;x179=l*x160 + x107 ...
    *x172 + x129*x168 + x150*x61 + x156*x55 + x157*x77 + x158*x57 + x162* ...
    x55 + x164*x59 + x166*x59 + x167*x71 + x170*x61 + x171*x75 + x174*x61 ...
     + x175*x60 + x176*x65 + x177*x87 + x178*x89;x180=B*x28 + D*x25  ...
    + 4*F*x18 + H*x10 + x147 + 4*x21 + 4*x23 + x3;x181=B ...
    *x98 + D*x100 + G*x144 + 4*x103 + x104 + x145 + x26*x32;x182= ...
    -x180*x98 + x181*x28;x183=d1Y*x182;x184=x0*x180 + x0*x181;x185=x184*x37;x186= ...
    S*x28 + U*x25 + Y*x10 + x153 + x18*x50 + x47 + x49;x187=S ...
    *x98 + U*x100 + x113 + x15*x45 + x151 + x26*x46;x188=-x17*x187 + ...
     x186*x96;x189=d0X*x188;x190=-x0*x186 - x0*x187;x191=x190*x56;x192=x186*x98 - ...
     x187*x28;x193=d1X*x192;x194=x186*x26 - x187*x22;x195=x100*x186 - x187*x25;x196= ...
    -x14*x187 + x186*x97;x197=d2X*x196;x198=x14*x181 - x180*x97;x199=x18*x181 -  ...
    x180*x92;x200=d5Y*x199;x201=x17*x181 - x180*x96;x202=d0Y*x201;x203=-x180*x94 +  ...
    x181*x9;x204=d7Y*x203;x205=-x18*x187 + x186*x92;x206=x186*x94 - x187*x9;x207= ...
    d7X*x206;x208=-x100*x180 + x181*x25;x209=d3Y*x208;x210=-x15*x180 + x181*x20;x211 ...
    =-x180*x26 + x181*x22;x212=d4Y*x211;x213=x15*x186 - x187*x20;x214=d6X*x213;x215 ...
    =l*x183 + l*x185 + l*x202 + l*x209 + x189*x55 + x191*x55 +  ...
    x193*x55 + x194*x77 + x195*x65 + x197*x55 + x198*x71 + x200*x61 + x204 ...
    *x61 + x205*x60 + x207*x59 + x210*x87 + x212*x61 + x214*x59;x216=x143* ...
    x179 + x143*x215;x217=d3X*l;x218=d5X*x61;x219=l*x56;x220=d2Y*x55;x221=d6Y*x59 ...
    ;x222=2*l;x223=d0X*l;x224=d3Y*x55;x225=2*d4X*x222*x76 + l*x69 +  ...
    l*x81 + x217*x64 + x218*x58 + x219*x53 + x220*x70 + x221*x86 + 2 ...
    *x222*x67 + x223*x74 + x224*x88 + x38*x55 + x40*x55 + 4*x55*x73  ...
    + x55*x83 + x59*x63 + x59*x79 + x61*x85;x226=x37*x55;x227=d1Y*x55;x228 ...
    =l*x119 + l*x124 + l*x139 + x106*x226 + x108*x221 + x115*x219 +  ...
    2*x117*x222 + x121*x61 + x122*x220 + x126*x55 + x127*x217 + x128*x227 + ...
     x130*x218 + x132*x59 + x134*x59 + x136*x59 + x137*x224 + 2*x141*x222; ...
    x229=d4X*x61;x230=l*x156 + l*x162 + x150*x59 + x157*x229 + x158*x219 + ...
     x160*x55 + x164*x61 + x166*x61 + x167*x220 + x168*x227 + x170*x59 +  ...
    x171*x223 + x172*x226 + x174*x59 + x175*x218 + x176*x217 + x177*x221 + x178 ...
    *x224;x231=l*x189 + l*x191 + l*x193 + l*x197 + x183*x55 + x185* ...
    x55 + x194*x229 + x195*x217 + x198*x220 + x200*x59 + x202*x55 + x204*x59 ...
     + x205*x218 + x207*x61 + x209*x55 + x210*x221 + x212*x59 + x214*x61; ...
    x232=x143*x230 + x143*x231;x233=4*u;x234=d4X*x233;x235=d3Y*u;x236=d1Y*u;x237 ...
    =d7Y*x233;x238=d6X*x233;x239=u*x56;x240=d5X*x233;x241=d0Y*u;x242=u*x37;x243= ...
    d2X*u;x244=d0X*u;x245=d3X*u;x246=d1X*u;x247=d6Y*x233;x248=d4Y*x54;x249=d5Y ...
    *x54;x250=d2Y*u;x251=d7X*x233;x252=x234*x72 + x235*x64 + x236*x68 + x237* ...
    x66 + x238*x86 + x239*x36 + x240*x78 + x241*x74 + x242*x53 + x243*x70 ...
     + x244*x39 + x245*x88 + x246*x82 + x247*x84 + 2*x248*x76 + 2* ...
    x249*x58 + x250*x80 + x251*x62;x253=d4X*x54;x254=x106*x239 + x108*x238 + x115 ...
    *x242 + 2*x116*x248 + x118*x241 + x120*x247 + x122*x243 + x123*x236 +  ...
    x125*x244 + x127*x235 + x128*x246 + 2*x130*x249 + 2*x131*x253 + x133*x251 ...
     + x135*x240 + x137*x245 + x138*x250 + x140*x237;x255=d4Y*x157*x233 + d5Y* ...
    x175*x233 + x149*x251 + x155*x250 + x158*x242 + x159*x244 + x161*x236 + x163 ...
    *x237 + x165*x247 + x167*x243 + x168*x246 + x169*x234 + x171*x241 + x172* ...
    x239 + x173*x240 + x176*x235 + x177*x238 + x178*x245;x256=x182*x246 + x184*x239 ...
     + x188*x241 + x190*x242 + x192*x236 + 2*x194*x248 + x195*x235 + x196* ...
    x250 + x198*x243 + x199*x240 + x201*x244 + x203*x251 + 2*x205*x249 + x206 ...
    *x237 + x208*x245 + x210*x238 + 2*x211*x253 + x213*x247;x257=x143*x255 +  ...
    x143*x256;x258=x142*x143 + x143*x90;x259=x143*x225 + x143*x228;x260=x143*x252 + x143 ...
    *x254;x261=x0 + 1;x262=x6 + 1;x263=x142*x262 + x261*x90;x264=x179*x262  ...
    + x215*x261;x265=x230*x262 + x231*x261;x266=x225*x261 + x228*x262;x267=x255*x262 + ...
     x256*x261;x268=x252*x261 + x254*x262;x269=x142*x261 + x262*x90;x270=x225*x262 +  ...
    x228*x261;x271=x252*x262 + x254*x261;x272=x179*x261 + x215*x262;x273=x230*x261 + x231 ...
    *x262;x274=x255*x261 + x256*x262;
    P=[x142*x91/46656 + x2*x90/46656 + x216/46656, x2*x225/46656 + x228*x91/ ...
    46656 + x232/46656, x2*x252/46656 + x254*x91/46656 + x257/46656, x179*x91/46656 + ...
     x2*x215/46656 + x258/46656, x2*x231/46656 + x230*x91/46656 + x259/46656, x2* ...
    x256/46656 + x255*x91/46656 + x260/46656, x142*x2/46656 + x216/46656 + x90*x91/ ...
    46656, x2*x228/46656 + x225*x91/46656 + x232/46656, x2*x254/46656 + x252*x91/46656 ...
     + x257/46656, x179*x2/46656 + x215*x91/46656 + x258/46656, x2*x230/46656 +  ...
    x231*x91/46656 + x259/46656, x2*x255/46656 + x256*x91/46656 + x260/46656, x263/11664 ...
     + x264/11664, x265/11664 + x266/11664, x267/11664 + x268/11664, x264/11664 + x269 ...
    /11664, x265/11664 + x270/11664, x267/11664 + x271/11664, x269/11664 + x272/11664,  ...
    x270/11664 + x273/11664, x271/11664 + x274/11664, x263/11664 + x272/11664, x266/11664  ...
    + x273/11664, x268/11664 + x274/11664, x142/2916 + x179/2916 + x215/2916 + x90 ...
    /2916, x225/2916 + x228/2916 + x230/2916 + x231/2916, x252/2916 + x254/2916  ...
    + x255/2916 + x256/2916];
end
