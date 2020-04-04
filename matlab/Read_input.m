
function [node, element, elemType, nel, nen, nIntPts,nnd,ps, nu, E, Force_Node, bforce, disp_BC]=Read_input(filename)
fid=fopen(filename, 'r');

textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
E=C{1,1};

textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nu=C{1,1};

textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
ps=C{1,1};

textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nnd=C{1,1};
C=textscan(fid,'%n %n','MultipleDelimsAsOne', 1);
X=C{1};
Y=C{2};
node=[X Y];

textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nel=C{1,1};
C=textscan(fid,'%n',1);
nen=C{1,1};
if (nen==4)
    elemType='Q4';
    C=textscan(fid,'%u %u %u %u','MultipleDelimsAsOne', 1);
    element=cell2mat(C(1:4));
elseif (nen==9)
    elemType='Q9';
    C=textscan(fid,'%u %u %u %u %u %u %u %u %u','MultipleDelimsAsOne', 1);
    element=cell2mat(C(1:9));
elseif (nen==8)
    elemType='Q8';
    C=textscan(fid,'%u %u %u %u %u %u %u %u','MultipleDelimsAsOne', 1);
    element=cell2mat(C(1:8));
end

textscan(fid,'%s',1);
C=textscan(fid,'%n',1);
nIntPts=C{1,1};

textscan(fid,'%s',7);
C=textscan(fid,'%n %n %n','MultipleDelimsAsOne', 1);
clear disp_BC
disp_BC=[C{1} C{2} C{3}];

textscan(fid,'%s',4);
C=textscan(fid,'%u','MultipleDelimsAsOne', 1);
Force_Node=C{1};
textscan(fid,'%s',2);
C=textscan(fid,'%n',1);
bforce=C{1,1};
fclose(fid);