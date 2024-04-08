
function [output]=FieldB_plasma(r_evaluate,Params)
%load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
%load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))
rho=r_evaluate(1); z=r_evaluate(3);
R=Params.R; Rp=Params.Rp; nodes=Params.MRWn; weights=Params.MRWw;



if Params.test
    subplot(1,2,1)
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
hold on; axis equal; plot(rho,z,'o')
%pause()
end

if ((rho-R)^2+z^2)<=Rp^2
    %Psimin=0; Psimax=2*pi;
    Intx=zeros(size(nodes)); funx=Intx; Intsumx=0;
    Intz=zeros(size(nodes)); funz=Intz; Intsumz=0;
    for i=1:length(nodes)
         [~,Dmax]=LimsDInt(rho,z,nodes(i)*2*pi,Params);

        for j=1:length(nodes)
            [funx(j),funz(j)]=IntegrandB(rho,z,nodes(j)*Dmax,nodes(i)*2*pi,Params);
            Intx(j)=2*pi*Dmax*funx(j)*weights(i)*weights(j);
            Intz(j)=2*pi*Dmax*funz(j)*weights(i)*weights(j);
        end  
        Intx(isnan(Intx)|isinf(Intx))=0;
        Intsumx=Intsumx+sum(Intx);
        Intz(isnan(Intz)|isinf(Intz))=0;
        Intsumz=Intsumz+sum(Intz);
        if Params.test
            subplot(1,2,2)
            plot(nodes,funx)
            shg
            pause();
        end
    end

else
 [Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params);
    Intx=zeros(size(nodes)); Intsumx=0; funx=Intx;
    Intz=zeros(size(nodes)); Intsumz=0; funz=Intz;
    for i=1:length(nodes)
        Psi=Psimin+nodes(i)*(Psimax-Psimin);
        [Dmin,Dmax]=LimsDExt(Psi,Dc,Psic,Params);
        for j=1:length(nodes)
            D=Dmin+nodes(j)*(Dmax-Dmin);
            [funx(j),funz(j)]=IntegrandB(rho,z,D,Psi,Params);
            Intx(j)=(Psimax-Psimin)*(Dmax-Dmin)*weights(i)*weights(j)*funx(j);
            Intz(j)=(Psimax-Psimin)*(Dmax-Dmin)*weights(i)*weights(j)*funz(j);

        end       
        Intx(isnan(Intx)|isinf(Intx))=0;
        Intsumx=Intsumx+sum(Intx);
        Intz(isnan(Intz)|isinf(Intz))=0;
        Intsumz=Intsumz+sum(Intz);
        if Params.test
            subplot(1,2,2)
            plot(nodes,funx)
            pause();
        end
    end


end
Bx=Intsumx; Bz=Intsumz;

%B=struct('Bx',Bx,'Bz',Bz);
output=struct();
output.B=[Bx,0,Bz]';
output.mod_B=norm(output.B);
end

function [Dmin,Dmax]=LimsDExt(Psi,Dc,Psic,Params)
Rp=Params.Rp;
Dmax=Dc*abs(cos(Psi-Psic))+sqrt(Rp^2-Dc^2*sin(Psi-Psic)^2);
Dmin=Dc*abs(cos(Psi-Psic))-sqrt(Rp^2-Dc^2*sin(Psi-Psic)^2);

end

function [Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params)
R=Params.R; Rp=Params.Rp;
Psic=atan2(-z,(R-rho));
Dc=sqrt((R-rho)^2+z^2);

Psimax=Psic+asin(Rp/Dc);
Psimin=Psic-asin(Rp/Dc);

end


function [Dmin,Dmax]=LimsDInt(rho,z,Psi,Params)
R=Params.R; Rp=Params.Rp;

Dmin=0;
Dmax=((R-rho)*cos(Psi)-z*sin(Psi))+sqrt(((R-rho)*cos(Psi)-z*sin(Psi))^2+(Rp^2-z^2-(R-rho)^2));

end


function [H1,H2]=IntegrandB(rho,z,D,Psi,Params)
% por hacer

rhop=rho+D*cos(Psi); zp=z+D*sin(Psi);
if Params.test
    subplot(1,2,1)
plot(rhop,zp,'x'); shg
%pause();
end

xi=(rho^2+rhop^2+(z-zp)^2)/(2*rho*rhop);

DerQhalf=DerQ12(xi);
Qhalf=Q12(xi);

Jphi=Params.J(rhop,zp);

H1=-Params.mu/(2*pi)*(z-zp)/(sqrt(rhop)*rho^(3/2))*DerQhalf*Jphi*D;

H2=Params.mu/(2*pi)*(sqrt(rhop)/(2*rho^(3/2))*Qhalf+...
    ((rho^2-rhop^2-(z-zp)^2)/(2*rho^(5/2)*sqrt(rhop)))*DerQhalf)*Jphi*D;

end


function DerQhalf= DerQ12(x)
% if imag(x)~=0
%     a=1;
% end

%[K,E]=ellipke(sqrt(2./(x+1)));
try
[K,E]=ellipke(2/(x+1));
DerQhalf=0.5*sqrt(2/(x+1))*K-x*(2*(x+1))^0.5/(2*(x^2-1))*E;
catch
DerQhalf=0; errorc=1;
end
end
function [Qhalf,errorc]= Q12(x)

% if imag(x)~=0
%     a=1;
% end
M=(2/(x+1));



try
%[K,E]=ellipke(sqrt(2./(x+1)));
[K,E]=ellipke(M);
Qhalf=x*sqrt(2/(x+1))*K-sqrt(2*(x+1))*E;
errorc=0;
catch
    %warning('error en q1/2')
    errorc=1; Qhalf=0;
end

end
