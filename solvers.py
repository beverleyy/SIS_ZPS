global np
import numpy as np

import warnings
warnings.filterwarnings("error")

def build_velocity(b,rho,dt,u,v,dx,dy):
    b[1:-1, 1:-1] = (rho*(1/dt*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+(v[2:,1:-1]-v[0:-2,1:-1])/(2*dy)) 
                          -((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))**2 
                          -2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy)*(v[1:-1,2:]-v[1:-1,0:-2])/(2*dx))
                          -((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))**2))
    return b

def build_pressure(p,dx,dy,b,nit):
    pn = np.empty_like(p)
    pn = p.copy()
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = (((pn[1:-1,2:]+pn[1:-1,0:-2])*dy**2+(pn[2:,1:-1]+pn[0:-2,1:-1])*dx**2)/
                       (2*(dx**2+dy**2)) - dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1])
        
        #boundary conditions
        p[-1,:] = 0       #top edge, p=0
        p[0,:] = p[1,:]   #bottom edge, dp/dy=0
        p[:,0] = p[:,1]   #left edge, dp/dx=0
        p[:,-1] = p[:,-2] #right edge, dp/dx=0
        
    return p

def calculate_cfl(u,v,dx,dy,dt):
    return np.max(np.sqrt(np.square(u*dt/dx/dx)+np.square(v*dt/dy/dy)))

def solve_cavity(u_lid,u,v,dt,dx,dy,p,rho,nu,nx,ny,nit):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.zeros((ny,nx))
    
    vel = np.empty_like(u)
    cfl = calculate_cfl(u,v,dx,dy,dt)
    
    try:
        un = u.copy()
        vn = v.copy()

        b = build_velocity(b,rho,dt,u,v,dx,dy)
        p = build_pressure(p,dx,dy,b,nit)

        #u-momentum
        u[1:-1,1:-1] = (un[1:-1,1:-1]-un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[1:-1,0:-2])
                       -vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[0:-2,1:-1])
                       -dt/(2*rho*dx)*(p[1:-1,2:]-p[1:-1,0:-2])
                        +nu*(dt/(dx**2)*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])
                        +dt/(dy**2)*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])))
        #v-momentum
        v[1:-1,1:-1] = (vn[1:-1,1:-1]-un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[1:-1,0:-2])
                       -vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[0:-2,1:-1])
                       -dt/(2*rho*dy)*(p[2:,1:-1]-p[0:-2,1:-1])
                        +nu*(dt/(dx**2)*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])
                        +dt/(dy**2)*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])))

        #no-slip
        u[0,:] = 0
        u[:,0] = 0
        u[:,-1] = 0
        u[-1,:] = u_lid
        v[0,:] = 0
        v[-1,:] = 0
        v[:,0] = 0
        v[:,-1] = 0

        vel = np.sqrt(np.square(u)+np.square(v))
        cfl = calculate_cfl(u,v,dx,dy,dt)
        run=True
            
    except RuntimeWarning:
        print(f"Solution diverged. CFL_max={cfl:.3}...")
        run=False
        
    return u,v,vel,p,run
