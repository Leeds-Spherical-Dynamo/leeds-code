% This script outputs poloidal field coefficients into the file
%  mag_PolBC.txt.in which can be read into the unified dynamo model
% in routine io_load_magPolBC in io.F90
% The format of mag_PolBC.txt.in is that each row mustave 3 integers
% and  two reals. The first integer is 1 for inner boundary or 2 for
% the outer boundary. The next two integers are the degree l and order m
% of the spherical harmonic coefficient, starting at l=1. The two
% real numbers are the real and imaginary parts f the poloidal harmonic
% coefficients. If there is no l,m entry, the  harmonic coefficients default 
%  zero.   
%
% Currently set up to generate the poloidal magnetic field coefficients
% required to give the Earth's field (year 2000) at the CMB.
% earthgauss_2000.txt contains the Gauss coefficients in nanotesla up to
% degree 12.
%
% Read in the Gauss coefficients
%
out = [];
fid = fopen('earthgauss_2000.txt','rt');
%assert(fid>=3,msg)
while ~feof(fid)
    str = fgetl(fid);
    vec = sscanf(str,'%f',[1,Inf]);
    num = numel(vec);
    if num
        out(end+1,1:num) = vec;
    end
end  
fclose(fid);
% Convert the Gauss coefficients to poloidal harmonic coefficients as
% defined in the unified dynamo code. The coefficients are scaled by an
% arbitrary dnorm factor, which can be varied according to how strng an
% imposed field is desired.
%
fid2 = fopen('mag_PolBC.txt.in','w')
rCMB=3486;
rS=6371;
fac=rCMB/rS;
eta=0.35;
r=0;
dnorm=2e05;
for n=1:13
    np1=n+1;
    for mp1=1:np1
        % n,m degree and order
        m=mp1-1;
        r=r+1;
        %  bdy=2 for coefficients suplliedat outer boundary only
        %  Add in some bdy=1 entries if an imposed inner boundary field
        % is desired. 
        bdy=2;
%  a cos m phi + b sin m phi corresponds to c exp i m phi with
%  Re{c} = a/2, Im(c)=-b/2 if m not equal to 0
%  Re{c} = a, Im{c}=0 if m=0. (Ashley Willis notes page 10)
        if m==0
            fac2=1.0
        else
            fac2=0.5
        end
% Re{B_r}_nm = q_nm = (n+1)g_nm / fac^(n+2)
% q_nm = n(n+1) P_nm / r where P_nm is poloidal 
% Note - sign for imaginary part because Im(c) = -b/2. r=1/(1-eta). 
%  ((-1)^m) factor because geomagnetic field is usually plotted with
% longitude 0 at the centre, whereas hammer plots usually have phi=pi
% at the centre. ((-1)^m) factor rotates phi through angle pi. 
        repart=((-1)^m)*fac2*(1-eta)*out(r,3)/((n*dnorm)*(fac^(n+2)));
        impart=-((-1)^m)*(1-eta)*fac2*out(r,4)/((n*dnorm)*(fac^(n+2)));
        fprintf(fid2,'%4d %4d %4d %16.8e %16.8e \n',bdy,n,m,repart,impart);
    end
end
fclose(fid2)


% Enter desired diffusivity model here
%
%for i=1:nn 
%    ddiffdr=ppval(diffspld,r(i));
%    fprintf(fid,'%12.8e %12.8e \n',diff(i),ddiffdr);
%end


