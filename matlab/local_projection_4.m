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
                           Y(1), Y(2), Y(3), Y(4), Y(5), ...
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
  
    x0=l + 2*u;x1=sqrt(3)/3;x2=x1 + 1;x3=x2/4;x4 ...
    =D*x3;x5=1 - x1;x6=x5/4;x7=B*x6;x8=-A*x6 + C*x3 ...
    ;x9=-x4 + x7 + x8;x10=U*x6;x11=S*x3;x12=-R*x6 + T* ...
    x3;x13=x10 - x11 + x12;x14=D*x6;x15=B*x3;x16=x14 - x15 + x8 ...
    ;x17=U*x3;x18=S*x6;x19=x12 - x17 + x18;x20=x13*x9 - x16*x19; ...
    x21=1./x20;x22=x0*x21;x23=x13*x6;x24=-x23;x25=x19*x6;x26=x24 + x25; ...
    x27=d0X*x26;x28=l*x21;x29=x16*x6;x30=x6*x9;x31=-x30;x32=x29 + x31; ...
    x33=d0Y*x32;x34=x19*x3;x35=x23 + x34;x36=d1X*x35;x37=x16*x3;x38=x30 + ...
     x37;x39=d3Y*x38;x40=-x29;x41=x3*x9;x42=-x41;x43=x40 + x42;x44=d1Y* ...
    x43;x45=-x37;x46=x41 + x45;x47=d2Y*x46;x48=x13*x3;x49=-x34;x50=x48 + ...
     x49;x51=d2X*x50;x52=-x48;x53=-x25;x54=x52 + x53;x55=d3X*x54;x56=x22* ...
    x27 + x22*x36 + x22*x51 + x22*x55 + x28*x33 + x28*x39 + x28*x44 ...
     + x28*x47;x57=x5^2/4;x58=x20*x57;x59=x2^2/4;x60=-A*x3 + ...
     C*x6;x61=-x14 + x15 + x60;x62=-R*x3 + T*x6;x63=x17 -  ...
    x18 + x62;x64=x4 + x60 - x7;x65=-x10 + x11 + x62;x66=x61*x63  ...
    - x64*x65;x67=x3*x63;x68=-x67;x69=x3*x65;x70=x68 + x69;x71=d0X*x70; ...
    x72=1./x66;x73=x0*x72;x74=l*x72;x75=x3*x64;x76=x3*x61;x77=-x76;x78= ...
    x75 + x77;x79=d0Y*x78;x80=-x75;x81=x6*x61;x82=-x81;x83=x80 + x82;x84 ...
    =d1Y*x83;x85=x6*x65;x86=x67 + x85;x87=d1X*x86;x88=x6*x64;x89=-x88;x90 ...
    =x81 + x89;x91=d2Y*x90;x92=x6*x63;x93=-x85;x94=x92 + x93;x95=d2X*x94 ...
    ;x96=x76 + x88;x97=d3Y*x96;x98=-x92;x99=-x69;x100=x72*(x98 + x99); ...
    x101=d3X*x100;x102=x66*(x0*x101 + x71*x73 + x73*x87 + x73*x95 + x74* ...
    x79 + x74*x84 + x74*x91 + x74*x97);x103=x13*x61 - x16*x65;x104=1./ ...
    x103;x105=x0*x104;x106=x52 + x85;x107=d0X*x106;x108=x48 + x69;x109=d1X*x108;x110 ...
    =l*x104;x111=x37 + x82;x112=d0Y*x111;x113=x29 + x81;x114=d3Y*x113;x115=x45  ...
    + x77;x116=d1Y*x115;x117=x40 + x76;x118=d2Y*x117;x119=x23 + x99;x120=d2X*x119 ...
    ;x121=x24 + x93;x122=d3X*x121;x123=x103*(x105*x107 + x105*x109 + x105*x120 + ...
     x105*x122 + x110*x112 + x110*x114 + x110*x116 + x110*x118);x124=x3*x5;x125 ...
    =-x19*x64 + x63*x9;x126=1./x125;x127=x0*x126;x128=x34 + x98;x129=d0X*x128 ...
    ;x130=x31 + x89;x131=d1Y*x130;x132=l*x126;x133=x30 + x80;x134=d2Y*x133;x135= ...
    x25 + x92;x136=d1X*x135;x137=x42 + x88;x138=d0Y*x137;x139=x53 + x67;x140=d2X ...
    *x139;x141=x41 + x75;x142=d3Y*x141;x143=x49 + x68;x144=d3X*x143;x145=x127*x129  ...
    + x127*x136 + x127*x140 + x127*x144 + x131*x132 + x132*x134 + x132*x138 + ...
     x132*x142;x146=x124*x125;x147=x123*x124 + x145*x146;x148=x22*x33 + x22*x39 +  ...
    x22*x44 + x22*x47 + x27*x28 + x28*x36 + x28*x51 + x28*x55;x149=x66 ...
    *(l*x101 + x71*x74 + x73*x79 + x73*x84 + x73*x91 + x73*x97 + ...
     x74*x87 + x74*x95);x150=x127*x131 + x127*x134 + x127*x138 + x127*x142 + ...
     x129*x132 + x132*x136 + x132*x140 + x132*x144;x151=x103*(x105*x112 + x105*x114 ...
     + x105*x116 + x105*x118 + x107*x110 + x109*x110 + x110*x120 + x110*x122) ...
    ;x152=x124*x151 + x146*x150;x153=u*x21;x154=d0X*x153*x32 + d0Y*x153*x26 +  ...
    d1X*x153*x43 + d1Y*x153*x35 + d2X*x153*x46 + d2Y*x153*x50 + d3X*x153*x38 ...
     + d3Y*x153*x54;x155=u*x72;x156=d3Y*u;x157=x66*(d0X*x155*x78 + d0Y*x155 ...
    *x70 + d1X*x155*x83 + d1Y*x155*x86 + d2X*x155*x90 + d2Y*x155*x94 +  ...
    d3X*x155*x96 + x100*x156);x158=u*x104;x159=x103*(d0X*x111*x158 + d0Y*x106*x158 ...
     + d1X*x115*x158 + d1Y*x108*x158 + d2X*x117*x158 + d2Y*x119*x158 + d3X* ...
    x113*x158 + x104*x121*x156);x160=u*x126;x161=d0X*x137*x160 + d0Y*x128*x160 +  ...
    d1X*x130*x160 + d1Y*x135*x160 + d2X*x133*x160 + d2Y*x139*x160 + d3X*x141*x160 ...
     + x126*x143*x156;x162=x124*x159 + x146*x161;x163=x125*x57;x164=x124*x20;x165=x102* ...
    x124 + x164*x56;x166=x124*x149 + x148*x164;x167=x124*x157 + x154*x164;x168=x20*x59 ...
    ;x169=x125*x59;
    P=[x102*x59 + x147 + x56*x58, x148*x58 + x149*x59 + x152, x154*x58 ...
     + x157*x59 + x162, x123*x59 + x145*x163 + x165, x150*x163 + x151*x59  ...
    + x166, x159*x59 + x161*x163 + x167, x102*x57 + x147 + x168*x56, x148* ...
    x168 + x149*x57 + x152, x154*x168 + x157*x57 + x162, x123*x57 + x145*x169 ...
     + x165, x150*x169 + x151*x57 + x166, x159*x57 + x161*x169 + x167];
end

function [P] = crunch_serindipity(A, B, C, D, E, F, G, H, ...
                                  R, S, T, U, V, W, X, Y, ...
                                  d0X, d1X, d2X, d3X, d4X, d5X, d6X, d7X, ...
                                  d0Y, d1Y, d2Y, d3Y, d4Y, d5Y, d6Y, d7Y, ...
                                  u, l)
  
    x0=sqrt(3);x1=x0/3;x2=1 - x1;x3=x2/3;x4=x0*x3;x5 ...
    =Y*x4;x6=x1 + 1;x7=x6/3;x8=x0*x7;x9=W*x8;x10=x0/6; ...
    x11=x10*x6;x12=x6/4;x13=1/6 - x12;x14=x11 + x13;x15=x10*x2;x16 ...
    =x0/12;x17=1/12 - x16;x18=x15 + x17;x19=x16 + 1/12;x20=x11  ...
    + x19;x21=1/6 - x2/4;x22=x15 + x21;x23=R*x22 + T*x20;x24 ...
    =-V/3 + X/3;x25=S*x14 + U*x18 + x23 + x24 - x5 - ...
     x9;x26=G*x8;x27=E*x4;x28=F/3 - H/3;x29=A*x22 + C* ...
    x20;x30=B*x18 + D*x14 - x26 - x27 + x28 + x29;x31=X*x8;x32 ...
    =V*x4;x33=W/3 - Y/3;x34=S*x18 + U*x14 + x23 - x31  ...
    - x32 + x33;x35=H*x4;x36=F*x8;x37=-E/3 + G/3;x38=B* ...
    x14 + D*x18 + x29 - x35 - x36 + x37;x39=x25*x30 - x34*x38;x40 ...
    =1./x39;x41=l*x40;x42=-x14*x38 + x18*x30;x43=d3Y*x42;x44=x20*x25 - ...
     x20*x34;x45=l + 2*u;x46=x40*x45;x47=x14*x25 - x18*x34;x48=x38/ ...
    3;x49=x0*x30;x50=-x3*x49 + x48;x51=d7Y*x50;x52=-x48 - x49*x7;x53 ...
    =d5Y*x52;x54=x22*x25 - x22*x34;x55=d0X*x54;x56=x0*x48;x57=x30/3;x58= ...
    x56*x6 + x57;x59=d6Y*x58;x60=x25/3;x61=x34*x4 - x60;x62=d7X*x61;x63 ...
    =x2*x56 - x57;x64=d4Y*x63;x65=-x14*x34 + x18*x25;x66=d1X*x65;x67=x34 ...
    *x8 + x60;x68=d5X*x67;x69=x0*x25;x70=x34/3;x71=-x69*x7 - x70;x72 ...
    =d6X*x71;x73=x20*x30 - x20*x38;x74=d2Y*x73;x75=x14*x30 - x18*x38;x76= ...
    d1Y*x75;x77=x22*x30 - x22*x38;x78=d0Y*x77;x79=-x3*x69 + x70;x80=d4X* ...
    x79;x81=d2X*x44*x46 + d3X*x46*x47 + x41*x43 + x41*x51 + x41*x53 + ...
     x41*x59 + x41*x64 + x41*x74 + x41*x76 + x41*x78 + x46*x55 +  ...
    x46*x62 + x46*x66 + x46*x68 + x46*x72 + x46*x80;x82=x2^2/4 - ...
     x3;x83=x39*x82;x84=x6^2/4 - x7;x85=W*x4;x86=Y*x8;x87=-x15 ...
    ;x88=x21 + x87;x89=-x11;x90=x19 + x89;x91=x17 + x87;x92=x13 + x89 ...
    ;x93=R*x92 + T*x91;x94=S*x88 + U*x90 + x24 + x85 + x86  ...
    + x93;x95=E*x8;x96=G*x4;x97=A*x92 + C*x91;x98=B*x90 + D ...
    *x88 + x28 + x95 + x96 + x97;x99=X*x4;x100=V*x8;x101=S*x90  ...
    + U*x88 + x100 + x33 + x93 + x99;x102=H*x8;x103=F*x4;x104=B ...
    *x88 + D*x90 + x102 + x103 + x37 + x97;x105=-x101*x104 + x94*x98 ...
    ;x106=1./x105;x107=l*x106;x108=x104/3;x109=x0*x6;x110=x98/3;x111=-x108*x109 ...
     - x110;x112=d4Y*x111;x113=x106*x45;x114=x101/3;x115=-x114 + x4*x94;x116=d6X ...
    *x115;x117=-x104*x88 + x90*x98;x118=d3Y*x117;x119=-x101*x91 + x91*x94;x120= ...
    d2X*x119;x121=-x104*x91 + x91*x98;x122=d2Y*x121;x123=d3X*x45;x124=-x101*x90 + ...
     x88*x94;x125=x106*x124;x126=x0*x2;x127=-x108*x126 + x110;x128=d6Y*x127;x129=- ...
    x101*x92 + x92*x94;x130=d0X*x129;x131=-x108 + x110*x126;x132=d5Y*x131;x133=x94/ ...
    3;x134=-x101*x8 - x133;x135=d7X*x134;x136=-x104*x90 + x88*x98;x137=d1Y*x136 ...
    ;x138=-x104*x92 + x92*x98;x139=d0Y*x138;x140=x108 + x109*x110;x141=d7Y*x140;x142 ...
    =x114 + x8*x94;x143=x106*x142;x144=-x101*x88 + x90*x94;x145=d1X*x144;x146=- ...
    x101*x4 + x133;x147=d5X*x146;x148=x105*(d4X*x143*x45 + x107*x112 + x107*x118  ...
    + x107*x122 + x107*x128 + x107*x132 + x107*x137 + x107*x139 + x107*x141 + ...
     x113*x116 + x113*x120 + x113*x130 + x113*x135 + x113*x145 + x113*x147 +  ...
    x123*x125);x149=R*x88 + T*x90;x150=S*x92 + U*x91 + x149 + x24  ...
    + x5 + x9;x151=A*x14 + C*x18;x152=B*x20 + D*x22 + x151 + ...
     x28 - x95 - x96;x153=R*x14 + T*x18;x154=S*x20 + U*x22 -  ...
    x100 + x153 + x33 - x99;x155=A*x88 + C*x90;x156=B*x92 + D*x91 ...
     + x155 + x35 + x36 + x37;x157=x150*x152 - x154*x156;x158=1./x157;x159= ...
    l*x158;x160=x0*x156;x161=x152/3;x162=x160*x7 - x161;x163=d4Y*x162;x164=x158*x45 ...
    ;x165=x150/3;x166=x154/3;x167=-x126*x165 - x166;x168=d6X*x167;x169=x158*(x150* ...
    x22 - x154*x91);x170=x152*x90 - x156*x18;x171=d2Y*x170;x172=x152*x91 - x156* ...
    x22;x173=d3Y*x172;x174=x156/3;x175=x126*x161 + x174;x176=d7Y*x175;x177=-x154*x4  ...
    - x165;x178=d7X*x177;x179=-x14*x156 + x152*x88;x180=d0Y*x179;x181=x160*x3 +  ...
    x161;x182=d6Y*x181;x183=x109*x161 - x174;x184=d5Y*x183;x185=x150*x18 - x154*x90;x186 ...
    =x14*x150 - x154*x88;x187=d0X*x186;x188=-x154*x8 + x165;x189=d5X*x188;x190=x152 ...
    *x92 - x156*x20;x191=d1Y*x190;x192=-x109*x165 + x166;x193=d4X*x192;x194=x150*x20 ...
     - x154*x92;x195=d1X*x194;x196=d2X*x164*x185 + x123*x169 + x159*x163 + x159* ...
    x171 + x159*x173 + x159*x176 + x159*x180 + x159*x182 + x159*x184 + x159*x191 ...
     + x164*x168 + x164*x178 + x164*x187 + x164*x189 + x164*x193 + x164*x195; ...
    x197=x12*x2 - x3/2 - x7/2;x198=x157*x197;x199=S*x22 + U*x20 + ...
     x153 + x24 - x85 - x86;x200=B*x91 + D*x92 + x155 + x26 +  ...
    x27 + x28;x201=S*x91 + U*x92 + x149 + x31 + x32 + x33;x202=B ...
    *x22 + D*x20 - x102 - x103 + x151 + x37;x203=x199*x200 - x201*x202; ...
    x204=1./x203;x205=x204*x45;x206=x201/3;x207=x199*x8 - x206;x208=d6X*x207;x209=l ...
    *x204;x210=x20*x200 - x202*x92;x211=d3Y*x210;x212=x199*x92 - x20*x201;x213=x18* ...
    x200 - x202*x90;x214=d2Y*x213;x215=x202/3;x216=x0*x215;x217=x200/3;x218=-x216* ...
    x6 + x217;x219=d6Y*x218;x220=x0*x200;x221=-x215 - x220*x3;x222=d5Y*x221;x223= ...
    -x18*x201 + x199*x90;x224=x204*x223;x225=x215 - x220*x7;x226=d7Y*x225;x227=x14* ...
    x200 - x202*x88;x228=d0Y*x227;x229=x199/3;x230=x201*x8 - x229;x231=d7X*x230;x232 ...
    =-x2*x216 - x217;x233=d4Y*x232;x234=x200*x22 - x202*x91;x235=d1Y*x234;x236=x204 ...
    *(x199*x4 + x206);x237=d4X*x236;x238=x199*x91 - x201*x22;x239=d1X*x238;x240= ...
    -x14*x201 + x199*x88;x241=d0X*x240;x242=x201*x4 + x229;x243=d5X*x242;x244=x203* ...
    (d2X*x224*x45 + d3X*x205*x212 + x205*x208 + x205*x231 + x205*x239 + x205* ...
    x241 + x205*x243 + x209*x211 + x209*x214 + x209*x219 + x209*x222 + x209*x226 ...
     + x209*x228 + x209*x233 + x209*x235 + x237*x45);x245=x196*x198 + x197*x244 ...
    ;x246=d3X*l;x247=d2X*l;x248=x246*x40*x47 + x247*x40*x44 + x41*x55 +  ...
    x41*x62 + x41*x66 + x41*x68 + x41*x72 + x41*x80 + x43*x46 + x46 ...
    *x51 + x46*x53 + x46*x59 + x46*x64 + x46*x74 + x46*x76 + x46* ...
    x78;x249=x105*(d3X*x107*x124 + d4X*x107*x142 + x107*x116 + x107*x120 + x107* ...
    x130 + x107*x135 + x107*x145 + x107*x147 + x112*x113 + x113*x118 + x113*x122 ...
     + x113*x128 + x113*x132 + x113*x137 + x113*x139 + x113*x141);x250=x158*x185 ...
    *x247 + x159*x168 + x159*x178 + x159*x187 + x159*x189 + x159*x193 + x159* ...
    x195 + x163*x164 + x164*x171 + x164*x173 + x164*x176 + x164*x180 + x164*x182 ...
     + x164*x184 + x164*x191 + x169*x246;x251=x203*(l*x237 + x204*x212*x246 + ...
     x205*x211 + x205*x214 + x205*x219 + x205*x222 + x205*x226 + x205*x228 +  ...
    x205*x233 + x205*x235 + x208*x209 + x209*x231 + x209*x239 + x209*x241 + x209 ...
    *x243 + x224*x247);x252=x197*x251 + x198*x250;x253=u*x40;x254=d0X*x253*x77 + ...
     d0Y*x253*x54 + d1X*x253*x75 + d1Y*x253*x65 + d2X*x253*x73 + d2Y*x253* ...
    x44 + d3X*x253*x42 + d3Y*x253*x47 + d4X*x253*x63 + d4Y*x253*x79 + d5X ...
    *x253*x52 + d5Y*x253*x67 + d6X*x253*x58 + d6Y*x253*x71 + d7X*x253*x50  ...
    + d7Y*x253*x61;x255=u*x106;x256=d3Y*u;x257=d4Y*u;x258=x105*(d0X*x138*x255  ...
    + d0Y*x129*x255 + d1X*x136*x255 + d1Y*x144*x255 + d2X*x121*x255 + d2Y*x119 ...
    *x255 + d3X*x117*x255 + d4X*x111*x255 + d5X*x131*x255 + d5Y*x146*x255 +  ...
    d6X*x127*x255 + d6Y*x115*x255 + d7X*x140*x255 + d7Y*x134*x255 + x125*x256 + ...
     x143*x257);x259=u*x158;x260=d0X*x179*x259 + d0Y*x186*x259 + d1X*x190*x259 + ...
     d1Y*x194*x259 + d2X*x170*x259 + d2Y*x185*x259 + d3X*x172*x259 + d4X*x162* ...
    x259 + d5X*x183*x259 + d5Y*x188*x259 + d6X*x181*x259 + d6Y*x167*x259 + d7X ...
    *x175*x259 + d7Y*x177*x259 + x158*x192*x257 + x169*x256;x261=u*x204;x262=x203* ...
    (d0X*x227*x261 + d0Y*x240*x261 + d1X*x234*x261 + d1Y*x238*x261 + d2X*x213* ...
    x261 + d2Y*x223*x261 + d3X*x210*x261 + d3Y*x212*x261 + d4X*x232*x261 + d5X ...
    *x221*x261 + d5Y*x242*x261 + d6X*x218*x261 + d6Y*x207*x261 + d7X*x225*x261  ...
    + d7Y*x230*x261 + x236*x257);x263=x197*x262 + x198*x260;x264=x157*x84;x265=x197* ...
    x39;x266=x148*x197 + x265*x81;x267=x197*x249 + x248*x265;x268=x197*x258 + x254*x265 ...
    ;x269=x39*x84;x270=x157*x82;x271=x39*x81;x272=x148*x7 + x271*x3;x273=x157*x196; ...
    x274=x244*x3 + x273*x7;x275=x248*x39;x276=x249*x7 + x275*x3;x277=x157*x250;x278 ...
    =x251*x3 + x277*x7;x279=x254*x39;x280=x258*x7 + x279*x3;x281=x157*x260;x282= ...
    x262*x3 + x281*x7;x283=x148*x3 + x271*x7;x284=x249*x3 + x275*x7;x285=x258 ...
    *x3 + x279*x7;x286=x244*x7 + x273*x3;x287=x251*x7 + x277*x3;x288=x262* ...
    x7 + x281*x3;
    P=[x148*x84 + x245 + x81*x83, x248*x83 + x249*x84 + x252, x254*x83 ...
     + x258*x84 + x263, x196*x264 + x244*x82 + x266, x250*x264 + x251*x82  ...
    + x267, x260*x264 + x262*x82 + x268, x148*x82 + x245 + x269*x81, x248* ...
    x269 + x249*x82 + x252, x254*x269 + x258*x82 + x263, x196*x270 + x244*x84 ...
     + x266, x250*x270 + x251*x84 + x267, x260*x270 + x262*x84 + x268, x272 ...
     + x274, x276 + x278, x280 + x282, x274 + x283, x278 + x284, x282 + ...
     x285, x283 + x286, x284 + x287, x285 + x288, x272 + x286, x276 + x287 ...
    , x280 + x288];
end

function [P] = crunch_bubble(A, B, C, D, E, F, G, H, I, ...
                             R, S, T, U, V, W, X, Y, Z, ...
                             d0X, d1X, d2X, d3X, d4X, d5X, d6X, d7X, d8X, ...
                             d0Y, d1Y, d2Y, d3Y, d4Y, d5Y, d6Y, d7Y, d8Y, ...
                             u, l)
  
    x0=sqrt(3);x1=x0/3;x2=1 - x1;x3=x2^2/4 - x2/3 ...
     + 1/9;x4=2*x0/9;x5=x1 + 1;x6=x1*x5;x7=x4 - x6; ...
    x8=x1*x2;x9=x4 - x8;x10=x0/6;x11=x10*x5;x12=x0/9;x13=-x12; ...
    x14=1/6 - x5/4;x15=x11 + x13 + x14;x16=x10*x2;x17=7*x0/36 ...
    ;x18=x16 - x17 + 1/12;x19=x0/36;x20=x11 - x19 + 1/12;x21= ...
    4*x0/9;x22=I*x21;x23=-x22;x24=x2/4;x25=1/6 - x24;x26=x13  ...
    + x16 + x25;x27=A*x26 + C*x20 + x23;x28=x4 - 1/3;x29=x4 ...
     + 1/3;x30=F*x29 + H*x28;x31=B*x18 + D*x15 + E*x9  ...
    + G*x7 + x27 + x30;x32=V*x28 + X*x29;x33=Z*x21;x34=-x33; ...
    x35=R*x26 + T*x20 + x34;x36=S*x15 + U*x18 + W*x7 + Y ...
    *x9 + x32 + x35;x37=E*x28 + G*x29;x38=B*x15 + D*x18 +  ...
    F*x7 + H*x9 + x27 + x37;x39=W*x29 + Y*x28;x40=S*x18 + ...
     U*x15 + V*x9 + X*x7 + x35 + x39;x41=x31*x36 - x38*x40; ...
    x42=1./x41;x43=l*x42;x44=-x15*x38 + x18*x31;x45=d3Y*x44;x46=l +  ...
    2*u;x47=x42*x46;x48=x15*x36 - x18*x40;x49=d3X*x48;x50=-x29*x38 +  ...
    x31*x7;x51=d5Y*x50;x52=-x15*x40 + x18*x36;x53=d1X*x52;x54=-x28*x40 + ...
     x36*x9;x55=d4X*x54;x56=x15*x31 - x18*x38;x57=d1Y*x56;x58=-x21*x36 + ...
     x21*x40;x59=d8X*x58;x60=-x28*x38 + x31*x9;x61=d7Y*x60;x62=x29*x31 - ...
     x38*x7;x63=d6Y*x62;x64=x20*x31 - x20*x38;x65=d2Y*x64;x66=x29*x36 -  ...
    x40*x7;x67=d5X*x66;x68=x26*x36 - x26*x40;x69=d0X*x68;x70=-x21*x31 +  ...
    x21*x38;x71=d8Y*x70;x72=-x29*x40 + x36*x7;x73=d6X*x72;x74=x26*x31 -  ...
    x26*x38;x75=d0Y*x74;x76=x28*x36 - x40*x9;x77=d7X*x76;x78=x20*x36 - x20 ...
    *x40;x79=d2X*x78;x80=x28*x31 - x38*x9;x81=d4Y*x80;x82=x41*(x43*x45 + ...
     x43*x51 + x43*x57 + x43*x61 + x43*x63 + x43*x65 + x43*x71 +  ...
    x43*x75 + x43*x81 + x47*x49 + x47*x53 + x47*x55 + x47*x59 + x47 ...
    *x67 + x47*x69 + x47*x73 + x47*x77 + x47*x79);x83=x5^2/4 - ...
     x5/3 + 1/9;x84=-x16;x85=x12 + x25 + x84;x86=-x4;x87=x6  ...
    + x86;x88=x8 + x86;x89=-x11;x90=x17 + x89 + 1/12;x91=x19 +  ...
    x84 + 1/12;x92=x12 + x14 + x89;x93=A*x92 + C*x91;x94=x86 - ...
     1/3;x95=x86 + 1/3;x96=F*x95 + H*x94 + x22;x97=B*x90  ...
    + D*x85 + E*x87 + G*x88 + x93 + x96;x98=V*x94 + X*x95 ...
     + x33;x99=R*x92 + T*x91;x100=S*x85 + U*x90 + W*x88 +  ...
    Y*x87 + x98 + x99;x101=E*x94 + G*x95 + x22;x102=B*x85 + D ...
    *x90 + F*x88 + H*x87 + x101 + x93;x103=W*x95 + Y*x94 +  ...
    x33;x104=S*x90 + U*x85 + V*x87 + X*x88 + x103 + x99;x105=x100 ...
    *x97 - x102*x104;x106=1./x105;x107=x106*x46;x108=x100*x85 - x104*x90;x109=d3X* ...
    x108;x110=l*x106;x111=-x102*x91 + x91*x97;x112=d2Y*x111;x113=x100*x90 - x104* ...
    x85;x114=d1X*x113;x115=x100*x87 - x104*x94;x116=d4X*x115;x117=-x102*x90 + x85* ...
    x97;x118=d1Y*x117;x119=-x102*x94 + x87*x97;x120=d7Y*x119;x121=-x102*x85 + x90 ...
    *x97;x122=d3Y*x121;x123=-x102*x88 + x95*x97;x124=d6Y*x123;x125=x100*x94 - x104 ...
    *x87;x126=d7X*x125;x127=x100*x88 - x104*x95;x128=-x102*x95 + x88*x97;x129=d5Y ...
    *x128;x130=-x102*x21 + x21*x97;x131=d8Y*x130;x132=x100*x91 - x104*x91;x133=d2X ...
    *x132;x134=x100*x21 - x104*x21;x135=d8X*x134;x136=x100*x95 - x104*x88;x137=d5X* ...
    x136;x138=-x102*x92 + x92*x97;x139=d0Y*x138;x140=x100*x92 - x104*x92;x141=d0X* ...
    x140;x142=-x102*x87 + x94*x97;x143=d4Y*x142;x144=x105*(d6X*x107*x127 + x107*x109 ...
     + x107*x114 + x107*x116 + x107*x126 + x107*x133 + x107*x135 + x107*x137  ...
    + x107*x141 + x110*x112 + x110*x118 + x110*x120 + x110*x122 + x110*x124 + ...
     x110*x129 + x110*x131 + x110*x139 + x110*x143);x145=x5/3;x146=x2/3;x147 ...
    =-x145/2 - x146/2 + x24*x5 + 1/9;x148=A*x15 + C*x18 + ...
     x23;x149=B*x20 + D*x26 + E*x7 + G*x9 + x148 + x30;x150= ...
    R*x85 + T*x90;x151=S*x92 + U*x91 + W*x87 + Y*x88 + x150 ...
     + x98;x152=A*x85 + C*x90;x153=B*x92 + D*x91 + F*x87 +  ...
    H*x88 + x101 + x152;x154=R*x15 + T*x18 + x34;x155=S*x20 + U ...
    *x26 + V*x7 + X*x9 + x154 + x39;x156=x149*x151 - x153*x155;x157= ...
    x149*x91 - x153*x26;x158=d3Y*x157;x159=1./x156;x160=l*x159;x161=x159*x46;x162=x151 ...
    *x26 - x155*x91;x163=d3X*x162;x164=x149*x90 - x153*x18;x165=d2Y*x164;x166=x151* ...
    x28 - x155*x88;x167=d7X*x166;x168=x149*x92 - x153*x20;x169=d1Y*x168;x170=x149*x88 ...
     - x153*x28;x171=d7Y*x170;x172=-x151*x21 - x155*x21;x173=d8X*x172;x174=x149*x94 ...
     - x153*x7;x175=d4Y*x174;x176=x151*x20 - x155*x92;x177=d1X*x176;x178=x149*x21  ...
    + x153*x21;x179=d8Y*x178;x180=x159*(x151*x9 - x155*x95);x181=x149*x95 - x153 ...
    *x9;x182=d6Y*x181;x183=x151*x29 - x155*x87;x184=x149*x85 - x15*x153;x185=d0Y* ...
    x184;x186=x149*x87 - x153*x29;x187=d5Y*x186;x188=x151*x7 - x155*x94;x189=d4X*x188 ...
    ;x190=x15*x151 - x155*x85;x191=d0X*x190;x192=x151*x18 - x155*x90;x193=d2X*x192; ...
    x194=x156*(d5X*x161*x183 + d6X*x180*x46 + x158*x160 + x160*x165 + x160*x169  ...
    + x160*x171 + x160*x175 + x160*x179 + x160*x182 + x160*x185 + x160*x187 + ...
     x161*x163 + x161*x167 + x161*x173 + x161*x177 + x161*x189 + x161*x191 +  ...
    x161*x193);x195=B*x91 + D*x92 + E*x88 + G*x87 + x152 + x96; ...
    x196=S*x26 + U*x20 + W*x9 + Y*x7 + x154 + x32;x197=B*x26 ...
     + D*x20 + F*x9 + H*x7 + x148 + x37;x198=S*x91 + U* ...
    x92 + V*x88 + X*x87 + x103 + x150;x199=x195*x196 - x197*x198;x200=1. ...
    /x199;x201=l*x200;x202=x195*x20 - x197*x92;x203=d3Y*x202;x204=x200*x46;x205=x196* ...
    x92 - x198*x20;x206=d3X*x205;x207=x195*x9 - x197*x95;x208=d5Y*x207;x209=x196*x91 ...
     - x198*x26;x210=d1X*x209;x211=x195*x26 - x197*x91;x212=d1Y*x211;x213=x196*x94  ...
    - x198*x7;x214=d7X*x213;x215=x195*x7 - x197*x94;x216=d7Y*x215;x217=x18*x195 - ...
     x197*x90;x218=d2Y*x217;x219=-x195*x21 - x197*x21;x220=d8Y*x219;x221=x196*x87 - ...
     x198*x29;x222=x196*x88 - x198*x28;x223=d4X*x222;x224=x196*x21 + x198*x21;x225= ...
    d8X*x224;x226=x196*x95 - x198*x9;x227=d5X*x226;x228=x195*x29 - x197*x87;x229=d6Y ...
    *x228;x230=x15*x195 - x197*x85;x231=d0Y*x230;x232=-x15*x198 + x196*x85;x233=d0X ...
    *x232;x234=-x18*x198 + x196*x90;x235=d2X*x234;x236=x195*x28 - x197*x88;x237=d4Y ...
    *x236;x238=x199*(d6X*x204*x221 + x201*x203 + x201*x208 + x201*x212 + x201*x216 ...
     + x201*x218 + x201*x220 + x201*x229 + x201*x231 + x201*x237 + x204*x206  ...
    + x204*x210 + x204*x214 + x204*x223 + x204*x225 + x204*x227 + x204*x233 + ...
     x204*x235);x239=x147*x194 + x147*x238;x240=x41*(x43*x49 + x43*x53 + x43* ...
    x55 + x43*x59 + x43*x67 + x43*x69 + x43*x73 + x43*x77 + x43*x79 ...
     + x45*x47 + x47*x51 + x47*x57 + x47*x61 + x47*x63 + x47*x65  ...
    + x47*x71 + x47*x75 + x47*x81);x241=d6X*l;x242=x106*x127;x243=x105*(x107 ...
    *x112 + x107*x118 + x107*x120 + x107*x122 + x107*x124 + x107*x129 + x107* ...
    x131 + x107*x139 + x107*x143 + x109*x110 + x110*x114 + x110*x116 + x110*x126 ...
     + x110*x133 + x110*x135 + x110*x137 + x110*x141 + x241*x242);x244=x159*x183 ...
    ;x245=x156*(d5X*l*x244 + x158*x161 + x160*x163 + x160*x167 + x160*x173 + ...
     x160*x177 + x160*x189 + x160*x191 + x160*x193 + x161*x165 + x161*x169 +  ...
    x161*x171 + x161*x175 + x161*x179 + x161*x182 + x161*x185 + x161*x187 + x180 ...
    *x241);x246=x200*x221;x247=x199*(x201*x206 + x201*x210 + x201*x214 + x201*x223  ...
    + x201*x225 + x201*x227 + x201*x233 + x201*x235 + x203*x204 + x204*x208 + ...
     x204*x212 + x204*x216 + x204*x218 + x204*x220 + x204*x229 + x204*x231 +  ...
    x204*x237 + x241*x246);x248=x147*x245 + x147*x247;x249=u*x42;x250=d6Y*u;x251= ...
    d5Y*u;x252=x41*(d0X*x249*x74 + d0Y*x249*x68 + d1X*x249*x56 + d1Y*x249* ...
    x52 + d2X*x249*x64 + d2Y*x249*x78 + d3X*x249*x44 + d3Y*x249*x48 + d4X ...
    *x249*x80 + d4Y*x249*x54 + d5X*x249*x50 + d6X*x249*x62 + d7X*x249*x60  ...
    + d7Y*x249*x76 + d8X*x249*x70 + d8Y*x249*x58 + x250*x42*x72 + x251*x42 ...
    *x66);x253=u*x106;x254=x105*(d0X*x138*x253 + d0Y*x140*x253 + d1X*x117*x253  ...
    + d1Y*x113*x253 + d2X*x111*x253 + d2Y*x132*x253 + d3X*x121*x253 + d3Y*x108 ...
    *x253 + d4X*x142*x253 + d4Y*x115*x253 + d5X*x128*x253 + d6X*x123*x253 +  ...
    d7X*x119*x253 + d7Y*x125*x253 + d8X*x130*x253 + d8Y*x134*x253 + x106*x136*x251 ...
     + x242*x250);x255=u*x159;x256=x156*(d0X*x184*x255 + d0Y*x190*x255 + d1X* ...
    x168*x255 + d1Y*x176*x255 + d2X*x164*x255 + d2Y*x192*x255 + d3X*x157*x255 + ...
     d3Y*x162*x255 + d4X*x174*x255 + d4Y*x188*x255 + d5X*x186*x255 + d6X*x181* ...
    x255 + d7X*x170*x255 + d7Y*x166*x255 + d8X*x178*x255 + d8Y*x172*x255 + x180 ...
    *x250 + x244*x251);x257=u*x200;x258=x199*(d0X*x230*x257 + d0Y*x232*x257 +  ...
    d1X*x211*x257 + d1Y*x209*x257 + d2X*x217*x257 + d2Y*x234*x257 + d3X*x202*x257 ...
     + d3Y*x205*x257 + d4X*x236*x257 + d4Y*x222*x257 + d5X*x207*x257 + d6X* ...
    x228*x257 + d7X*x215*x257 + d7Y*x213*x257 + d8X*x219*x257 + d8Y*x224*x257 + ...
     x200*x226*x251 + x246*x250);x259=x147*x256 + x147*x258;x260=x144*x147 + x147*x82 ...
    ;x261=x147*x240 + x147*x243;x262=x147*x252 + x147*x254;x263=x146 - 2/9;x264= ...
    x145 - 2/9;x265=x144*x264 + x263*x82;x266=x194*x264 + x238*x263;x267=x245*x264 ...
     + x247*x263;x268=x240*x263 + x243*x264;x269=x252*x263 + x254*x264;x270=x256*x264  ...
    + x258*x263;x271=x144*x263 + x264*x82;x272=x240*x264 + x243*x263;x273=x252*x264 + ...
     x254*x263;x274=x194*x263 + x238*x264;x275=x245*x263 + x247*x264;x276=x256*x263 +  ...
    x258*x264;
    P=[x144*x83 + x239 + x3*x82, x240*x3 + x243*x83 + x248, x252*x3 ...
     + x254*x83 + x259, x194*x83 + x238*x3 + x260, x245*x83 + x247*x3  ...
    + x261, x256*x83 + x258*x3 + x262, x144*x3 + x239 + x82*x83, x240* ...
    x83 + x243*x3 + x248, x252*x83 + x254*x3 + x259, x194*x3 + x238*x83 ...
     + x260, x245*x3 + x247*x83 + x261, x256*x3 + x258*x83 + x262, x265 ...
     + x266, x267 + x268, x269 + x270, x266 + x271, x267 + x272, x270 + ...
     x273, x271 + x274, x272 + x275, x273 + x276, x265 + x274, x268 + x275 ...
    , x269 + x276, 4*x144/9 + 4*x194/9 + 4*x238/9 + 4*x82/ ...
    9, 4*x240/9 + 4*x243/9 + 4*x245/9 + 4*x247/9, 4*x252 ...
    /9 + 4*x254/9 + 4*x256/9 + 4*x258/9];
end

