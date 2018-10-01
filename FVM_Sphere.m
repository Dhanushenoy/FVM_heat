%FVM method for a hemi-sphere
clear all;
clc;

%Parameters----------------------
a=1;            %radius
n=100;          %angular discretization
m=100;          %radial discretization
alpha=0.1;      %diffusivity (k/rho*Cp)*Vc
Tinf=0;         %ambient Temperature
dt=0.001;       %Time Step C<0.5
tmax=50.0;      %max Time
k=20;           %conductivity 

%discretization------------------
dr=a/m;     
dtheta=pi/n;

%initialize----------------------
T=zeros(n+1,m+1);       %Temperature array time=m
Tnew=zeros(n+1,m+1);    %Temperature at time=m+1
T(:,:)=100.0;           %initial value

time=0.0;               %beginning of time march

%--------------------------------
%       Calculation
%--------------------------------
while time<=tmax
for i=1:n+1         %angular
    for j=1:m+1     %radial
     
%------------Control volume faces------------------------
        N(i,j)=2*pi*sin(dtheta*(i-1))*(dr*(j-1)+dr/2);
        S(i,j)=2*pi*sin(dtheta*(i-1))*(dr*(j-1)-dr/2);
        E(i,j)=2*pi*(dr*(j-1))*sin(dtheta*(i-1)+dtheta/2);
        W(i,j)=2*pi*(dr*(j-1))*sin(dtheta*(i-1)-dtheta/2);
%--------------------------------------------------------

%------------Boundary Conditions-------------------------
        %region at the center
        if (j==1) && (i==1)
                
        Tnew(i,j)=T(i,j)+alpha*dt*((((T(i,j+1)-T(i,j))/dr)*N(i,j))...
                    +((1/(dr*j)*(T(i+1,j)-T(i,j))/dtheta)*E(i,j)));

        elseif (j==1) && (i==n+1)
        Tnew(i,j)=T(i,j)+alpha*dt*(((T(i,j+1)-T(i,j))/dr)*N(i,j)...
                    -((1/(dr*j)*(T(i,j)-T(i-1,j))/dtheta)*W(i,j)));

        elseif(j==1) && (i<n+1) && (i>1)
        Tnew(i,j)=T(i,j)+alpha*dt*(((T(i,j+1)-T(i,j))/dr)*N(i,j)...
                    +((1/(dr*j)*(T(i+1,j)-T(i,j))/dtheta)*E(i,j))...
                    -((1/(dr*j)*(T(i,j)-T(i-1,j))/dtheta)*W(i,j)));
           
        
        %at angle 0 and pi
        elseif(i==1) && j>1 && j<m+1
        Tnew(i,j)=T(i,j)+alpha*dt*(((T(i,j+1)-T(i,j))/dr)*N(i,j)...
                    -((T(i,j)-T(i,j-1))/dr*S(i,j))...
                    +((1/(dr*j)*(T(i+1,j)-T(i,j))/dtheta)*E(i,j)));

        elseif(i==n+1) && j>1 && j<m+1
        Tnew(i,j)=T(i,j)+alpha*dt*(((T(i,j+1)-T(i,j))/dr)*N(i,j)...
                    -((T(i,j)-T(i,j-1))/dr*S(i,j))...
                    -((1/(dr*j)*(T(i,j)-T(i-1,j))/dtheta)*W(i,j)));
     
        
        elseif(j==m+1) && i>1 && i<n+1
        Tnew(i,j)=T(i,j)+alpha*dt*((k*sin(dtheta*j)*(T(i,j)-Tinf)*N(i,j))...
                    -((T(i,j)-T(i,j-1))/dr*S(i,j))...
                    +((1/(dr*j)*(T(i+1,j)-T(i,j))/dtheta)*E(i,j))...
                    -((1/(dr*j)*(T(i,j)-T(i-1,j))/dtheta)*W(i,j)));

        elseif(j==m+1) && i==1
        Tnew(i,j)=T(i,j)+alpha*dt*((k*sin(dtheta*j)*(T(i,j)-Tinf)*N(i,j))...
                    -((T(i,j)-T(i,j-1))/dr*S(i,j))...
                    +((1/(dr*j)*(T(i+1,j)-T(i,j))/dtheta)*E(i,j)));

        elseif(j==m+1) && i==n+1
        Tnew(i,j)=T(i,j)+alpha*dt*((k*sin(dtheta*j)*(T(i,j)-Tinf)*N(i,j))...
                    -((T(i,j)-T(i,j-1))/dr*S(i,j))...
                    -((1/(dr*j)*(T(i,j)-T(i-1,j))/dtheta)*W(i,j)));
%---------------------------------------------------------------
        else
            
        Tnew(i,j)=T(i,j)+alpha*dt*(((T(i,j+1)-T(i,j))/dr)*N(i,j)...
                    -((T(i,j)-T(i,j-1))/dr*S(i,j))...
                    +((1/(dr*j)*(T(i+1,j)-T(i,j))/dtheta)*E(i,j))...
                    -((1/(dr*j)*(T(i,j)-T(i-1,j))/dtheta)*W(i,j)));
        end
%---------------------------------------------------------------
    end
end
T(:,:)=Tnew(:,:);   %Replacing with new values
time=time+dt;       %Time stepping
end

%Corrdinate transform
for i=1:n+1
for j=1:m+1
    X(i,j)=(dr*(j-1))*cos(dtheta*(i-1));
    Y(i,j)=(dr*(j-1))*sin(dtheta*(i-1));     
end
end

%plot
set(gcf, 'Position', [100, 100, 1000, 500])
colormap(jet)
pcolor(X,Y,T)