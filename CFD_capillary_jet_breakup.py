%%Numerical Simulation of Capillary Jet Breakup
Re =200 ;
k = 0.45;
epsilon = 0.05 ;
rho = 1 ;
sigma = 1 ;
r = 0.5 ;
CFL=0.001;
tc= sqrt((rho*r^3)/sigma);
lamda = (2*pi*r)/k ; %% the Z domain is defined by r and k
mu = sqrt ( rho * sigma * r)/ Re;

%%creating staggered mesh
Lx=lamda; %influenced by r and k
zmin=1;
zmax=60;
z(zmin:zmax)= linspace (0,Lx,zmax);
dz_z= z(zmin+1)-z(zmin);
zm(zmin:zmax-1)=0.5*(z(zmin:zmax-1)+z(zmin+1:zmax));
dz_zm= zm(zmin+1)-zm(zmin);

dt_viscous = (dz_zm^2 * rho)/(mu * 2) ; % influenced by Re, r , dz_zm
dt_tension = sqrt (( rho* (dz_zm^3) )/ sigma ) ;

% initialization of f value across z
f_z(zmin:zmax)= (r*(1+epsilon*cos((2*pi*z(zmin:zmax))/lamda))).^2; % f=h^2
f_zm(zmin:zmax-1)= (r*(1+epsilon*cos((2*pi*zm(zmin:zmax-1))/lamda))).^2; % f=h^2

%initialization of velocity
u = zeros(zmin,zmax-1);

dt=0.01;
dt_courant=0.06;
fmin=r^2;
T=0;
%THE LOOP
while fmin > ((r^2)*1e-6)
    A= [ dt_viscous ; dt_tension ; dt_courant ];
    if dt <= min(A)
        dt=min(A)*CFL;
        T=T+dt/tc
        %U_avg (we count minus or not at the ends)??
        u_avg = zeros(zmin,zmax);
        for n=zmin:zmax
            if n==zmin
                u_avg(n)=0.5*(u(n)-u(n)); %u(n-1)=-u(n)
            elseif n==zmax
                u_avg(n)=0.5*(u(n-1)-u(n-1)); %u(n)=-u(n-1)
            else
                u_avg(n)=0.5*(u(n)+u(n-1)); % at f points not u
            end
        end     
        %find out new f_new based on new u_new new f_z
        f_new = zeros(zmin,zmax);
       for w=zmin:zmax
            if w==zmin
                f_new(w)=f_z(w)-f_z(w)*dt*((u(w)+u(w))/dz_z)-(u_avg(w)*dt*(f_z(w+1)-f_z(w+1))/(2*dz_z));
            elseif w==zmax
                f_new(w)=f_z(w)-f_z(w)*dt*((-u(w-1)-u(w-1))/dz_z)-(u_avg(w)*dt*(f_z(w-1)-f_z(w-1))/(2*dz_z));
            else
                f_new(w)=f_z(w)-f_z(w)*dt*((u(w)-u(w-1))/dz_z)-(u_avg(w)*dt*(f_z(w+1)-f_z(w-1))/(2*dz_z));
            end
       end
        fmin=min(f_new);
        %str = sprintf('Drop Evolution Time=%f (Re=%03.1f,k=%03.2f)',(T/0.353553391),Re,k);
        fz = zeros(zmin,zmax);
        for i=zmin:zmax
            if i==zmin
                fz(i)= 0 ;
            elseif i==zmax
                fz(i)= 0 ;
            else
                fz(i)= (f_new(i+1)-f_new(i-1))/(2*dz_z) ;
            end
        end
        
        %fzz
        fzz = zeros(zmin,zmax);
        for j=zmin:zmax
            if j==zmin
                fzz(j)= (f_new(j+1)-2*f_new(j)+f_new(j+1))/(dz_z^2) ;
            elseif j==zmax
                fzz(j)= (f_new(j-1)-2*f_new(j)+f_new(j-1))/(dz_z^2) ;
            else
                fzz(j)= (f_new(j+1)-2*f_new(j)+f_new(j-1))/(dz_z^2) ;
            end
        end
        
        %H
        H = zeros(zmin,zmax);
        for q=zmin:zmax
            H(q)=0.25*((f_z(q)*(2-fzz(q))+(fz(q))^2)/((f_z(q)+0.25*(fz(q))^2)^(3/2)));
        end
        
        %Si(x points) there are zmax points of H but we r taking zmax-1
        Si= zeros(zmin,zmax-1);
        for m=zmin:zmax-1
            Si(m)=(((-2)*sigma)/rho)*((H(m+1)-H(m))/dz_z);
        end
        
        %f_avg (x points)(perfecto)
        f_avg = zeros(zmin,zmax-1);
        for n=zmin:zmax-1
            f_avg(n)=0.5*(f_new(n)+f_new(n+1));
        end
          
        %Bi (x points)(perfecto)
        Bi = zeros(zmin,zmax-1);
        for p=zmin:zmax-1
            if p==zmin
                Bi(p)= ((3*mu)/(rho*f_avg(p)))*((f_new(p+1)-f_new(p))/dz_zm)*((u(p+1)+u(p))/(2*dz_zm))+ (((3*mu)/rho)*((u(p+1)-2*u(p)-u(p))/(dz_zm^2)));
            elseif p==zmax-1
                Bi(p)= ((3*mu)/(rho*f_avg(p)))*((f_new(p+1)-f_new(p))/dz_zm)*((-u(p)-u(p-1))/(2*dz_zm))+ (((3*mu)/rho)*((-u(p)-2*u(p)+u(p-1))/(dz_zm^2)));
            else
                Bi(p)= ((3*mu)/(rho*f_avg(p)))*((f_new(p+1)-f_new(p))/dz_zm)*((u(p+1)-u(p-1))/(2*dz_zm))+ (((3*mu)/rho)*((u(p+1)-2*u(p)+u(p-1))/(dz_zm^2)));
            end
        end
        u_new = zeros(zmin,zmax-1);
        %THE LOOP OF U
        for s=zmin:zmax-1
            if s==zmin
                u_new(s) = u(s)-(u(s)*(u(s+1)+u(s))*dt)/(2*dz_zm)+Bi(s)*dt+Si(s)*dt; % -U1=U1
            elseif s==zmax-1
                u_new(s) = u(s)-(u(s)*(-u(s)-u(s-1))*dt)/(2*dz_zm)+Bi(s)*dt+Si(s)*dt; % Un-1=-Un+1
            else
                u_new(s) = u(s)-(u(s)*(u(s+1)-u(s-1))*dt)/(2*dz_zm)+Bi(s)*dt+Si(s)*dt; 
            end
        end
        umax=max(u_new);
        dt_courant = dz_zm/umax;
        u=u_new;
        f_z=f_new;
    else
        dt=dt*0.5;
    end
end
        figure(2)
        plot(z,f_new);
        str = sprintf('Drop Evolution Time=%f (Re=%03.1f,k=%03.2f)',(T),Re,k);
        title(str);
        hold on
        plot(z,-f_new);
