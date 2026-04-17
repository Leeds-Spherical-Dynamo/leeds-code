% Diffusivity model for Boussinesq case.
% Outer core has nn mesh points inner core nic mesh points.
% nic=1 is no inner core
nn=64;
nic=16;
%
% Print out the required quatities at the meshpoints
% Radius ratio is rr, usually rr=0.35 for Earth 
rr=0.35;
ri =rr/(1-rr); 
r=zeros(1,nn);
for i=1:nn
    cch=0.5*(1.0+cos(pi*(nn-i)/(nn-1.0)));
    r(i)=ri+cch;
end
ric=zeros(1,nic);
n2ic=2*nic;
for i = nic+1: n2ic
    ric(i-nic) =  ri * cos(pi*(n2ic - i)/(n2ic-1));
end
r
ric
fid = fopen('diffusivity.txt','w');
% Enter desired diffusivity model here
%
diff=zeros(1,nn);
for i=1:nn
    diff(i) = 1.0 + 0.1*(r(nn)-r(i));
end
diffspl=spline(r,diff);
diffspld=mmppder(diffspl);
for i=1:nn 
    ddiffdr=ppval(diffspld,r(i));
    fprintf(fid,'%12.8e %12.8e \n',diff(i),ddiffdr);
end
diffic=zeros(1,nic);       
for i=1:nic
    diffic(i)=1.1+0.1*(ri-ric(i));
end
splic=spline(ric,diffic);
difficspld=mmppder(splic);
for i=1:nic
    difficd=ppval(difficspld,ric(i));
    fprintf(fid,'%12.8e %12.8e \n',diffic(i),difficd);
end
status=fclose(fid);


