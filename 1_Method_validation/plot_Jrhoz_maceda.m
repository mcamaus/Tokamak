Jrhoz=@(rho,z) 1*(((rho-R).^2+(z).^2)<=Htheta(atan2(z,rho-R))^2);

rhovec=linspace(R-2*Rp,R+2*Rp,100);
zvec=linspace(-2*Rp,2*Rp,100);
J=zeros(length(rhovec),length(zvec));

for i=1:length(rhovec)
    for j=1:length(zvec)
        angu=atan2(zvec(j),rhovec(i));
        H=Htheta(angu);
        J(i,j)=Jrhoz(rhovec(i),zvec(j));
        hold on
        %scatter(R+H*cos(angu),H*sin(angu),15)
    end
end

[rhomat,zmat]=ndgrid(rhovec,zvec);

surf(rhomat,zmat,J)

shg
view(2)
axis equal
shading interp
