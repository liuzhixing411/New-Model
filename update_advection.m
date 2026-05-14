function alpha_new=update_advection(Fluid,grid)
Ny=grid.Ny;
Nx=grid.Nx;
start=grid.start;
endy=grid.endy;
endx=grid.endx;
alpha=Fluid.alpha;
alpha_star=zeros(Ny,Nx);
U=Fluid.u;
V=Fluid.v;
h=grid.h;
dt=grid.dt;
ghostnum=grid.ghostnum;
%% 第一步：x方向扫描，得到alpha_star
for j=start:endy
    for i=start:endx
        eps=1e-6;
        d0=1/10;d1=3/5;d2=3/10;
        
        % 在i+1/2面上的alpha值 (右面)
        if(U(j,i+1)>=0)
            beta0=13/12*(alpha(j,i-2)-2*alpha(j,i-1)+alpha(j,i))^2+1/4*(alpha(j,i-2)-4*alpha(j,i-1)+3*alpha(j,i))^2;
            beta1=13/12*(alpha(j,i-1)-2*alpha(j,i)+alpha(j,i+1))^2+1/4*(alpha(j,i-1)-alpha(j,i+1))^2;
            beta2=13/12*(alpha(j,i)-2*alpha(j,i+1)+alpha(j,i+2))^2+1/4*(3*alpha(j,i)-4*alpha(j,i+1)+alpha(j,i+2))^2;
            f0=1/3*alpha(j,i-2)-7/6*alpha(j,i-1)+11/6*alpha(j,i);
            f1=-1/6*alpha(j,i-1)+5/6*alpha(j,i)+1/3*alpha(j,i+1);
            f2=1/3*alpha(j,i)+5/6*alpha(j,i+1)-1/6*alpha(j,i+2);
        else
            beta0=13/12*(alpha(j,i+3)-2*alpha(j,i+2)+alpha(j,i+1))^2+1/4*(alpha(j,i+3)-4*alpha(j,i+2)+3*alpha(j,i+1))^2;
            beta1=13/12*(alpha(j,i+2)-2*alpha(j,i+1)+alpha(j,i))^2+1/4*(alpha(j,i+2)-alpha(j,i))^2;
            beta2=13/12*(alpha(j,i+1)-2*alpha(j,i)+alpha(j,i-1))^2+1/4*(3*alpha(j,i+1)-4*alpha(j,i)+alpha(j,i-1))^2;
            f0=1/3*alpha(j,i+3)-7/6*alpha(j,i+2)+11/6*alpha(j,i+1);
            f1=-1/6*alpha(j,i+2)+5/6*alpha(j,i+1)+1/3*alpha(j,i);
            f2=1/3*alpha(j,i+1)+5/6*alpha(j,i)-1/6*alpha(j,i-1);
        end
        alpha0=d0/(beta0+eps)^2;alpha1=d1/(beta1+eps)^2;alpha2=d2/(beta2+eps)^2;
        sum_alpha=alpha0+alpha1+alpha2;
        omega0=alpha0/sum_alpha;omega1=alpha1/sum_alpha;omega2=alpha2/sum_alpha;
        alphastar_phalf=omega0*f0+omega1*f1+omega2*f2;
        % 在i-1/2面上的alpha值 (左面)
        if(U(j,i)>=0)
            beta0=13/12*(alpha(j,i-3)-2*alpha(j,i-2)+alpha(j,i-1))^2+1/4*(alpha(j,i-3)-4*alpha(j,i-2)+3*alpha(j,i-1))^2;
            beta1=13/12*(alpha(j,i-2)-2*alpha(j,i-1)+alpha(j,i))^2+1/4*(alpha(j,i-2)-alpha(j,i))^2;
            beta2=13/12*(alpha(j,i-1)-2*alpha(j,i)+alpha(j,i+1))^2+1/4*(3*alpha(j,i-1)-4*alpha(j,i)+alpha(j,i+1))^2;
            f0=1/3*alpha(j,i-3)-7/6*alpha(j,i-2)+11/6*alpha(j,i-1);
            f1=-1/6*alpha(j,i-2)+5/6*alpha(j,i-1)+1/3*alpha(j,i);
            f2=1/3*alpha(j,i-1)+5/6*alpha(j,i)-1/6*alpha(j,i+1);
        else
            beta0=13/12*(alpha(j,i+2)-2*alpha(j,i+1)+alpha(j,i))^2+1/4*(alpha(j,i+2)-4*alpha(j,i+1)+3*alpha(j,i))^2;
            beta1=13/12*(alpha(j,i+1)-2*alpha(j,i)+alpha(j,i-1))^2+1/4*(alpha(j,i+1)-alpha(j,i-1))^2;
            beta2=13/12*(alpha(j,i)-2*alpha(j,i-1)+alpha(j,i-2))^2+1/4*(3*alpha(j,i)-4*alpha(j,i-1)+alpha(j,i-2))^2;
            f0=1/3*alpha(j,i+2)-7/6*alpha(j,i+1)+11/6*alpha(j,i);
            f1=-1/6*alpha(j,i+1)+5/6*alpha(j,i)+1/3*alpha(j,i-1);
            f2=1/3*alpha(j,i)+5/6*alpha(j,i-1)-1/6*alpha(j,i-2);
        end
        
        alpha0=d0/(beta0+eps)^2;alpha1=d1/(beta1+eps)^2;alpha2=d2/(beta2+eps)^2;
        sum_alpha=alpha0+alpha1+alpha2;
        omega0=alpha0/sum_alpha;omega1=alpha1/sum_alpha;omega2=alpha2/sum_alpha;
        alpha_mhalf=omega0*f0+omega1*f1+omega2*f2;
        % x方向更新：dα/dt + d(αu)/dx = 0
        term1=(alphastar_phalf*U(j,i+1)-alpha_mhalf*U(j,i))/h;
        term2=(U(j,i+1)-U(j,i))/h*alpha(j,i);
        alpha_star(j,i)=alpha(j,i)+dt*(-term1+term2);
        alpha_star(j,i) = min(max(alpha_star(j,i),0),1);
    end
end
% 对alpha_star应用边界条件
alpha_star = apply_boundary_conditions('vof', grid, alpha_star);
%% 第二步：y方向扫描，得到alpha_new（使用alpha_star作为输入）
alpha_new=zeros(Ny,Nx);
for j=start:endy
    for i=start:endx
        eps=1e-6;
        d0=1/10;d1=3/5;d2=3/10;
        
        % 在j+1/2面上的alpha值 (上面) - **修正：使用alpha_star**
        if(V(j+1,i)>=0)
            beta0=13/12*(alpha_star(j-2,i)-2*alpha_star(j-1,i)+alpha_star(j,i))^2+1/4*(alpha_star(j-2,i)-4*alpha_star(j-1,i)+3*alpha_star(j,i))^2;
            beta1=13/12*(alpha_star(j-1,i)-2*alpha_star(j,i)+alpha_star(j+1,i))^2+1/4*(alpha_star(j-1,i)-alpha_star(j+1,i))^2;
            beta2=13/12*(alpha_star(j,i)-2*alpha_star(j+1,i)+alpha_star(j+2,i))^2+1/4*(3*alpha_star(j,i)-4*alpha_star(j+1,i)+alpha_star(j+2,i))^2;
            f0=1/3*alpha_star(j-2,i)-7/6*alpha_star(j-1,i)+11/6*alpha_star(j,i);
            f1=-1/6*alpha_star(j-1,i)+5/6*alpha_star(j,i)+1/3*alpha_star(j+1,i);
            f2=1/3*alpha_star(j,i)+5/6*alpha_star(j+1,i)-1/6*alpha_star(j+2,i);
        else
            beta0=13/12*(alpha_star(j+3,i)-2*alpha_star(j+2,i)+alpha_star(j+1,i))^2+1/4*(alpha_star(j+3,i)-4*alpha_star(j+2,i)+3*alpha_star(j+1,i))^2;
            beta1=13/12*(alpha_star(j+2,i)-2*alpha_star(j+1,i)+alpha_star(j,i))^2+1/4*(alpha_star(j+2,i)-alpha_star(j,i))^2;
            beta2=13/12*(alpha_star(j+1,i)-2*alpha_star(j,i)+alpha_star(j-1,i))^2+1/4*(3*alpha_star(j+1,i)-4*alpha_star(j,i)+alpha_star(j-1,i))^2;
            f0=1/3*alpha_star(j+3,i)-7/6*alpha_star(j+2,i)+11/6*alpha_star(j+1,i);
            f1=-1/6*alpha_star(j+2,i)+5/6*alpha_star(j+1,i)+1/3*alpha_star(j,i);
            f2=1/3*alpha_star(j+1,i)+5/6*alpha_star(j,i)-1/6*alpha_star(j-1,i);
        end
        alpha0=d0/(beta0+eps)^2;alpha1=d1/(beta1+eps)^2;alpha2=d2/(beta2+eps)^2;
        sum_alpha=alpha0+alpha1+alpha2;
        omega0=alpha0/sum_alpha;omega1=alpha1/sum_alpha;omega2=alpha2/sum_alpha;
        alphastar_phalf=omega0*f0+omega1*f1+omega2*f2;
        
        if(V(j,i)>=0)
            beta0=13/12*(alpha_star(j-3,i)-2*alpha_star(j-2,i)+alpha_star(j-1,i))^2+1/4*(alpha_star(j-3,i)-4*alpha_star(j-2,i)+3*alpha_star(j-1,i))^2;
            beta1=13/12*(alpha_star(j-2,i)-2*alpha_star(j-1,i)+alpha_star(j,i))^2+1/4*(alpha_star(j-2,i)-alpha_star(j,i))^2;
            beta2=13/12*(alpha_star(j-1,i)-2*alpha_star(j,i)+alpha_star(j+1,i))^2+1/4*(3*alpha_star(j-1,i)-4*alpha_star(j,i)+alpha_star(j+1,i))^2;
            f0=1/3*alpha_star(j-3,i)-7/6*alpha_star(j-2,i)+11/6*alpha_star(j-1,i);
            f1=-1/6*alpha_star(j-2,i)+5/6*alpha_star(j-1,i)+1/3*alpha_star(j,i);
            f2=1/3*alpha_star(j-1,i)+5/6*alpha_star(j,i)-1/6*alpha_star(j+1,i);
        else
            beta0=13/12*(alpha_star(j+2,i)-2*alpha_star(j+1,i)+alpha_star(j,i))^2+1/4*(alpha_star(j+2,i)-4*alpha_star(j+1,i)+3*alpha_star(j,i))^2;
            beta1=13/12*(alpha_star(j+1,i)-2*alpha_star(j,i)+alpha_star(j-1,i))^2+1/4*(alpha_star(j+1,i)-alpha_star(j-1,i))^2;
            beta2=13/12*(alpha_star(j,i)-2*alpha_star(j-1,i)+alpha_star(j-2,i))^2+1/4*(3*alpha_star(j,i)-4*alpha_star(j-1,i)+alpha_star(j-2,i))^2;
            f0=1/3*alpha_star(j+2,i)-7/6*alpha_star(j+1,i)+11/6*alpha_star(j,i);
            f1=-1/6*alpha_star(j+1,i)+5/6*alpha_star(j,i)+1/3*alpha_star(j-1,i);
            f2=1/3*alpha_star(j,i)+5/6*alpha_star(j-1,i)-1/6*alpha_star(j-2,i);
        end
        
        alpha0=d0/(beta0+eps)^2;alpha1=d1/(beta1+eps)^2;alpha2=d2/(beta2+eps)^2;
        sum_alpha=alpha0+alpha1+alpha2;
        omega0=alpha0/sum_alpha;omega1=alpha1/sum_alpha;omega2=alpha2/sum_alpha;
        alphastar_mhalf=omega0*f0+omega1*f1+omega2*f2;
        % y方向更新
        termstar1=(alphastar_phalf*V(j+1,i)-alphastar_mhalf*V(j,i))/h;
        termstar2=(V(j+1,i)-V(j,i))/h*alpha_star(j,i);
        alpha_new(j,i)=alpha_star(j,i)+dt*(-termstar1+termstar2);
        alpha_new(j,i) = min(max(alpha_new(j,i),0),1);
    end
end

% 边界条件
alpha_new  = apply_boundary_conditions('vof', grid, alpha_new);
end

