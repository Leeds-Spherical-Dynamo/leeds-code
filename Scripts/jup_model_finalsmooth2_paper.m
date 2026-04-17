% Smoothed Jupiter model, based on French et al 2012.
% If rtil is r/rjup, form of xi taken as 
% xi = 25.2rtil^3 - 31.5rtil^2 + 11.4rtil - 0.68 - 1.14/(0.9964-rtil)
% then by analyitic integration
% log rho =6.3rtil^4-10.5rtil^3+5.7rtil^2-0.68rtil+1.14log(0.9964-rtil)+8.53
% where the constant of integration is taken as 8.53 giving 
% rho=4.42e03 kg/m^3 at rtil=0.0923, the inner core boundary. 
% These formulae valid for the French et al range 0.0923 < rtil < 0.994.
% Temperature gradient now a cubic polynomial plus another term
% of the form quadratic/(1.09-rtil)*((0.9319-rtil)^2+0.0025)
% parameters chosen so that the model is close to
% the Jupiter model of
% French et al Ap.J Supp series 202:5  2012 model J11-8a
%global rinp rhoinp
%
%  Data from French et al Table 1
%
%  Note the rho data at 0.629 is the midpoint of the two values
%  given by French et al. This point not used in te best fit procedure
%
%global rinp Tinp np rtil1
format long
Tinp=[1000 2000 2500 3000 3500 4000 4400 4500 4600 4700 4800 ...
    5000 6000 8000 10000 11000 12000 14000 16000 18000 19500];
rinp=[0.994 0.986 0.980 0.972 0.964 0.952 0.940 0.930 0.918 0.908 0.900...
    0.890 0.852 0.770 0.680 0.629 0.584 0.478 0.350 0.196 0.0923];
rhoinp=[0.0113 0.0497 0.0848 0.1320 0.1770 0.2400 0.3120 0.3800 0.4640 0.5310 0.5810...
    0.638 0.824 1.2000 1.6000 1.9200 2.2600 2.8000 3.3900 3.9900 4.4200];
pi=4.0*atan(1.0);
rjup=6.9894e07;
ri=0.0923*rjup
rr=rjup*rinp;
temp=spline(rr,Tinp);
tempd=mmppder(temp);
rhos=spline(rr,rhoinp*1e03);
rhosd=mmppder(rhos);
rhosdd=mmppder(rhosd);
np=200;
for i=1:np
    rtil2(i)=0.0923+(i-1)*(1-0.0923)/(np-1);
    rho2(i)=ppval(rhos,rtil2(i)*rjup);
    xi2(i)=rjup*ppval(rhosd,rtil2(i)*rjup)/rho2(i);
    dxidr2(i)=rjup*rjup*ppval(rhosdd,rtil2(i)*rjup)/rho2(i)-xi2(i)^2;
    temp2(i)=ppval(temp,rtil2(i)*rjup);
    dtempdr2(i)=rjup*ppval(tempd,rtil2(i)*rjup);
%    poltrop(i)=4482.0*sin(pi*rtil2(i))/(pi*rtil2(i));
end
for i=1:21
    dtemdr(i)=rjup*ppval(tempd,rinp(i)*rjup);
    xidat(i)=rjup*ppval(rhosd,rinp(i)*rjup)/ppval(rhos,rinp(i)*rjup);
end
%
% Form of dT/dr = te1*rtil^3 + te2*rtil^2 + te3*rtil +te4
%  + tac*((taa-rtil)^2+tab)/((tb1-rrtil)*((tb2-rtil)^2+tb3))
% Constants te1 etc found by least squares fit to the Tinp data
% from the French model 
t1=2.4198e05
t2=-2.9725e05
t3=1.1094e05
t4=-0.1786e05
%te4=-0.0520
% rational polynomial part is
% +tac*((taa-rtil^2+tab)/((tb1-rtil)((tb2-rtil)^2+tb3)
a1=0.9136;
a2=sqrt(0.0002517);
c1=-0.5987e04;
%
% Rational part can be written 
% tac*(rtil^2+ta1*rtil+ta2)/((tb1-rtil)((tb2-rtil)^2+tb3)
%ta1=-2*taa;
%ta2=taa*taa+tab;
b1=1.0248
b2=0.9097;
b3=sqrt(0.0005081);
for i=1:np
    ctemp=c1*((rtil2(i)-a1)^2+a2^2);
    q3(i)=t1*rtil2(i)^3+t2*rtil2(i)^2+t3*rtil2(i)+t4+ ...
        ctemp/((b1-rtil2(i))*(((b2-rtil2(i))^2)+b3^2));
end
%
% Split rational part into partial fractions
% tac*(rtil^2+ta1*rtil+ta2)/((tb1-rtil)((tb2-rtil)^2+tb3) =
% tac*( tc1/(tb2-rtil) + (tc2*rtil+tc3)/((tb2-rtil)^2 + tb3) )
%
c2=c1*((a1-b1)^2+a2^2)/(b3^2+(b2-b1)^2)
c3=(2*c1*a1+b1*(c2-c1)-b2*(c1+c2))/b3
q4=[t1 t2 t3 t4]
% Uncomment this section and use with tempmin11 to do the least squares fit
%
%for jj=1:10
%v=[te1/1e05 te2/1e05 te3/1e05 te4/1e05 taa tab tac/1e04 tb1 tb2 tb3]
%options=optimset('MaxFunEvals',30000,'MaxIter',20000,'TolX',1e-9,'TolFun',1e-9,'Display','final');
%[tea,fval,exitflag,output] = fminsearch(@tempmin11,v,options)
%fval
%taas=tea(5)*100
%tabs=tea(6)*1000
%tb1s=tea(8)*100
%tb2s=tea(9)*100
%tb3s=tea(10)*1000
%te1=tea(1)*1e05;
%te2=tea(2)*1e05;
%te3=tea(3)*1e05;
%te4=tea(4)*1e05;
%q4=[tea(1)*1e05 tea(2)*1e05 tea(3)*1e05 tea(4)*1e05]
%taa=tea(5);
%tab=tea(6);
%tac=tea(7)*1e04;
%tb1=tea(8);
%tb2=tea(9);
%tb3=tea(10);
%ta1=-2*taa;
%ta2=taa*taa+tab;
%tc1=(tb1*tb1+ta1*tb1+ta2)/(tb3+(tb2-tb1)^2)
%tc2=tc1-1
%tc3=-ta1-2*tb2*tc1+tb1*tc2
%end
roots(q4)
for i=1:np
    ctemp=c1*((rtil2(i)-a1)^2+a2^2);
    q3(i)=polyval(q4,rtil2(i)) + ...
     ctemp/((b1-rtil2(i))*(((b2-rtil2(i))^2)+b3^2));
end
% Uncomment two lines below to see Temperature Gradient plot
plot(rtil2,dtempdr2,'r',rtil2,q3,'b',rinp,dtemdr,'xk')
axis([0.0923 1 -60000 0]);
% 
%  Integrate dT/dr
%
q5=polyint(q4)
for i=1:np
    r=rtil2(i);
    tempii2(i)=polyval(q5,r)-c2*log(b1-r)+0.5*(c2-c1)*log((r-b2)^2+b3^2) ...
        +c3*atan((r-b2)/b3);
end
Tcint=-tempii2(1)+19500.0
for i=1:np
    tempii(i)=tempii2(i)+Tcint;
end
%
% Evaluate least squares fit value: with these parameters in the formula
% its 1.63868e+04
%
sumsq=0;
for i=1:21
     r=rinp(i);
    tempii3(i)=polyval(q5,r)-c2*log(b1-r)+0.5*(c2-c1)*log((r-b2)^2+b3^2) ...
        +c3*atan((r-b2)/b3)+Tcint;
end
for i=1:21
    sumsq=sumsq+(Tinp(i)-tempii3(i))^2;
end
sumsq
% Uncomment here to get temperature plot, with or without the spline fit
%plot(rtil2,temp2,'r',rtil2,tempii,'b',rinp,Tinp,'xk')
%plot(rtil2,tempii,'b',rinp,Tinp,'xk')
%axis([0.0923 1 0 20000]);
%
% Do density structure
% Approximation type for xi is the same as the approximation type for
% dT/dr that is cubic + quadratic/cubic
%
r1=29.6363
r2=-29.0005
r3=10.6323
r4=0.239974
rhoq4=[r1 r2 r3 r4]
d1=0.646554 
d2=sqrt(0.0010387)
f1=-1.9152
%a1=-2*aa;
%a2=aa*aa+ab;
e1=1.005638
e2=0.6462
e3=sqrt(0.0007259);
f2=f1*((d1-e1)^2+d2^2)/(e3^2+(e1-e2)^2)
f3=(2*f1*d1+e1*(f2-f1)-e2*(f2+f1))/e3
% Uncomment to do minimisation over rho using rhomin11
% Note that the globals need to be set for either rho or t minimisation
%for jj=1:1
%v=[e1 e2 e3 e4 aa ab ac b1 b2 b3]
%%v=[te1 te2 te3]
%options=optimset('MaxFunEvals',30000,'MaxIter',20000,'TolX',1e-9,'TolFun',1e-9,'Display','final');
%[tea,fval,exitflag,output] = fminsearch(@rhomin11,v,options)
%fval
%e1=tea(1)
%e2=tea(2)
%e3=tea(3)
%e4=tea(4)
%rhoq4=[tea(1) tea(2) tea(3) tea(4)]
%aa=tea(5);
%ab=tea(6);
%ac=tea(7);
%b1=tea(8);
%b2=tea(9);
%b3=tea(10);
%end
%a1=-2*aa;
%a2=aa*aa+ab;
%c1=(b1*b1+a1*b1+a2)/(b3+(b2-b1)^2)
%c2=c1-1
%c3=-a1-2*b2*c1+b1*c2
for i=1:np
    ctemp=f1*((rtil2(i)-d1)^2+d2^2);
    xii(i)=polyval(rhoq4,rtil2(i)) + ...
     ctemp/((e1-rtil2(i))*(((e2-rtil2(i))^2)+e3^2));
end
% Uncomment to plot xi
%plot(rtil2,xi2,'r',rtil2,xii,'b',rinp,xidat,'xk')
%axis([0.0923 0.994 -50 0]);
%
% Integrate to get log(rho)
%
rhoq5=polyint(rhoq4)
for i=1:np
    r=rtil2(i);
    lrhoii2(i)=polyval(rhoq5,r)-f2*log(e1-r)+0.5*(f2-f1)*log((r-e2)^2+e3^2) ...
        +f3*atan((r-e2)/e3);
end
rhocint=-lrhoii2(1)+log(4420)
for i=1:np
    lrhoii(i)=lrhoii2(i)+rhocint;
end
sumsq=0;
for i=1:21
     r=rinp(i);
    lrhoii3(i)=polyval(rhoq5,r)-f2*log(e1-r)+0.5*(f2-f1)*log((r-e2)^2+e3^2) ...
        +f3*atan((r-e2)/e3)+rhocint;
end
for i=1:15
    sumsq=sumsq+(1000*rhoinp(i)-exp(lrhoii3(i)))^2;
end
for i=17:21
    sumsq=sumsq+(1000*rhoinp(i)-exp(lrhoii3(i)))^2;
end
% least squares parameter currently 72.119
sumsq
%Uncomment to get rho with or without spline fit
%plot(rtil2,rho2,'r',rtil2,exp(lrhoii),'b',rinp,1000*rhoinp,'xk')
for i=1:15
    rhopic1(i)=rhoinp(i)*1000;
    rpic1(i)=rinp(i);
end
for i=17:21
    rhopic2(i-16)=rhoinp(i)*1000;
    rpic2(i-16)=rinp(i);
end
rhopic3(1)=1000*1.82;
rhopic4(1)=1000*2.02;
rpic3(1)=rinp(16);
rpic4(1)=rinp(16);
%plot(rtil2,exp(lrhoii),'b',rpic1,rhopic1,'xk',rpic2,rhopic2,'xk', ...
%    rpic3,rhopic3,'xk',rpic4,rhopic4,'xk')
%axis([0.0923 1.0 0 5000]);
%
% Differentiate to get dxi/dr
%
rhoq4d=polyder(rhoq4);
for i=1:np
    r=rtil2(i);
    dxiidr(i)=polyval(rhoq4d,r)+f2/((e1-r)^2) ...
        +(f2-f1)/(((r-e2)^2+e3^2)) ...
     +2*(e2-r)*((f2-f1)*r+f1*(2*d1-e1)+f2*(e1-2*e2))/(((r-e2)^2+e3^2)^2);
end
% Uncomment to plot dxi/dr
%plot(rtil2,dxidr2,'r',rtil2,dxiidr,'b')
%axis([0.0923 0.994 -700 60]);
%
% Introduce the model cut-off radius and evaluate the midpoint values 
% eta is the radius ratio used in the program
%
rcut_off=6.7e07
eta=ri/rcut_off
rho_cut_off=ppval(rhos,rcut_off)
r=rcut_off/rjup;
rho_cut_off=exp(polyval(rhoq5,r)-f2*log(e1-r)+0.5*(f2-f1)*log((r-e2)^2+e3^2) ...
        +f3*atan((r-e2)/e3)+rhocint)
rmid=0.5*(rcut_off+ri)
d=rcut_off-ri
rmidtil=rmid/rjup
lrhomid=polyval(rhoq5,rmidtil)-f2*log(e1-rmidtil)+ ...
  0.5*(f2-f1)*log((rmidtil-e2)^2+e3^2) ...
    +f3*atan((rmidtil-e2)/e3)+rhocint;
rhomid=exp(lrhomid)
Tmid=polyval(q5,rmidtil)-c2*log(b1-rmidtil)+ ...
  0.5*(c2-c1)*log((rmidtil-b2)^2+b3^2) ...
  +c3*atan((rmidtil-b2)/b3)+Tcint
%
% Put in the electrical conductivity here from table 2 of French et al
%
%rinp2=[0.0923 0.196 0.350 0.478 0.584 0.680 0.770 0.852 0.890 0.930 0.952 0.972 0.980 0.986];
%sigmainp=[3.39e06 3.05e06 2.48e06 1.97e06 1.48e06 1.09e06 7.20e05 3.60e05 ...
%    1.46e05 1.5e03 6.60 3.5e-02 3.5e-04 1e-07];
%logsig=log(sigmainp);
%cond=spline(rinp2*rjup,logsig);
%conddiff=mmppder(cond);
%sigmid=exp(ppval(cond,rmid))
%diffusivitymid=1/(4*pi*1e-07*sigmid)
%
% Hyperbolic diffusivity model
%
daa=-4.2791e-06;
dbb=274.9;
%  dbb = 180.0 givesv a Saturn-like diffusivity
% Also, the diffusivity gets very big near surface, so if it was greater than 10^6
% I set it to 10^6, and derivative I cut off at 10^8.
dcc=-2.544e-08;
ddd=1.801;
dee=20.28;
bbb=(daa+dcc)*rmid+dbb+ddd;
ccc=(daa*rmid+dbb)*(dcc*rmid+ddd)-dee;
ldiffmid=(-bbb+sqrt(bbb*bbb-4*ccc))/2.0;
diffusivitymid=exp(ldiffmid)
%
%  Set the number of mesh points nn used for the run here
%
nn=160;
%
% Print out the required quatities at the meshpoints
%
for i=1:nn
    cch=0.5*(1.0+cos(pi*(nn-i)/(nn-1.0)));
    chebmesh(i)=ri+d*cch;
end
fid = fopen('planet_model.txt.in','w');
for i=1:nn
    rr=chebmesh(i);
    bbb=(daa+dcc)*rr+dbb+ddd;
    ccc=(daa*rr+dbb)*(dcc*rr+ddd)-dee;
    ldiff=(-bbb+sqrt(bbb*bbb-4*ccc))/2.0;
    ldiffdr=-0.5*(daa+dcc)+(0.5*bbb*(daa+dcc)-2*daa*dcc*rr-dbb*dcc-daa*ddd)/sqrt(bbb*bbb-4*ccc);
    rtil=chebmesh(i)/rjup;
    T=polyval(q5,rtil)-c2*log(b1-rtil)+ ...
      0.5*(c2-c1)*log((rtil-b2)^2+b3^2) ...
      +c3*atan((rtil-b2)/b3)+Tcint;
    ctemp=c1*((rtil-a1)^2+a2^2);
    dtdrtil=    polyval(q4,rtil) + ...
           ctemp/((b1-rtil)*(((b2-rtil)^2)+b3^2));
    dtdr=dtdrtil/rjup;
    diffusivity=exp(ldiff);
    ddiffdr=ldiffdr*diffusivity;
    lrho1=polyval(rhoq5,rtil)-f2*log(e1-rtil)+ ...
        0.5*(f2-f1)*log((rtil-e2)^2+e3^2) ...
        +f3*atan((rtil-e2)/e3)+rhocint;
    rho1=exp(lrho1);
    rhodim=rho1/rhomid;
    ctemp=f2*((rtil-d1)^2+d2^2);
    xidim=(polyval(rhoq4,rtil) + ...
        ctemp/((e1-rtil)*(((e2-rtil)^2)+e3^2)))*(d/rjup);
    dxidrdim= d*d*(polyval(rhoq4d,rtil)+f2/((e1-rtil)^2) ...
        +(f2-f1)/(((rtil-e2)^2+e3^2)) ...
      +2*(e2-rtil)*((f2-f1)*rtil+f1*(2*d1-e1) ...
      +f2*(e1-2*e2))/(((rtil-e2)^2+e3^2)^2))/(rjup^2);
    Tdim=T/Tmid;
    dtdrdim=d*dtdr/Tmid;
    diffusivitydim=diffusivity/diffusivitymid;
    ddiffdrdim=d*ddiffdr/diffusivitymid;
    fprintf(fid,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8e %12.8e \n',rhodim ...
        ,xidim,dxidrdim,Tdim,dtdrdim,diffusivitydim,ddiffdrdim);
end
%aa=eta/(1.0-eta);
%bb=1.0/(1.0-eta);
%tol=1d-06;
%flux=quad(@fun_int,aa,bb,tol)
status=fclose(fid);


