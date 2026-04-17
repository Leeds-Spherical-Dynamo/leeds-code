% Anelasticbenchmark model

%
%  Set the number of mesh points nn used for the run here
%
nn=128;
%
% Print out the required quatities at the meshpoints
%
beta=0.35
ri=beta/(1-beta)
for i=1:nn
    cch=0.5*(1.0+cos(pi*(nn-i)/(nn-1.0)));
    chebmesh(i)=ri+cch;
end
fid = fopen('planet_model.txt.in','w');
% The formulae below used to generate the reference state are detailed in
% Jones et al. 2011, Icarus vol 216, pp 120-135. 
for i=1:nn
    rr=chebmesh(i);
    r(i)=rr; 
    c0=-0.459754;
    c1=1.515898;
    n=2;
    zeta(i)=c0+c1/rr;
    rho(i)=zeta(i)^n;
    xi(i)=-n*c1/((rr^2)*zeta(i));
    dxidr(i)=2*n*c1/((rr^3)*zeta(i))-n*c1*c1/((rr^4)*zeta(i)*zeta(i));
    T(i)=zeta(i)/c1;
    dTdr(i)=-1/(rr*rr);
    eta(i)=1;
    detadr(i)=0;
    fprintf(fid,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8e %12.8e \n',rho(i) ...
        ,xi(i),dxidr(i),T(i),dTdr(i),eta(i),detadr(i));
end
r
rho
xi
dxidr
T
%aa=eta/(1.0-eta);
%bb=1.0/(1.0-eta);
%tol=1d-06;
%flux=quad(@fun_int,aa,bb,tol)
status=fclose(fid);


