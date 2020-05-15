## Calculating Potentials

## This function Generates a Projected Potential using potentials from Kirkland suitable for TEM,  i.e intensity is proportional to Z^1/3
## Usage: ``V = stempot(xmax,ymax,nx,ny,potfile)``
##Inputs: xmax, ymax are the size of the slice in angstroms
##        nx,ny are the number of pixels in the x and y directions

def stempot(xmax,ymax,nx,ny,potfile):
    
    zed =2 #0.67;  # zed=2 for rutherford scattering of the nucleus, less for screening
    ix = np.arange(1,nx+1)
    iy = np.arange(1,ny+1)
    
    dx = xmax/(nx-1)
    dy = ymax/(ny-1)
    rx = np.arange(0,(xmax-dx),dx)
    ry = np.arange(0,(ymax-dy),dy)
    
    # Loading atomic positions from potfile
    M = np.loadtxt(potfile, delimiter=',')
    Zatom = M[:,0]
    ax= M[:,1]
    ay=M[:,2]
    az =M[:,3]
    wt=M[:,4]
    tds=0
    amax = len(Zatom)
    
    #find boundaries of the axes
    axmin = ax.min(0)
    axmax = ax.max(0)
    aymin = ay.min(0)
    aymax = ay.max(0)
    
    #shift coords to fit in box
    ax = ax - axmin
    ay = ay - aymin
    
    V= np.zeros((nx,ny))
    # map x and y coords of the atoms to the nearest grid points
    # A fraction of the atom must be assigned to the closest gridpoints
    # to avoid sum and difference frequencies appearing in the image
    
    #grid point to the left of the atom
    iax = np.maximum(1,np.minimum(np.floor(ax/dx)+1,nx)).astype(int)
    iax = iax -1
    
    #create periodic boundary conditions
    ibx = np.mod(iax,nx) #+ 1
    
    #fraction of atom at iax 
    fax = 1-np.mod((ax/dx),1)
    
    #grid point above the atom
    iay = np.maximum(1,np.minimum(np.floor(ay/dy)+1,ny)).astype(int) 
    iay = iay-1
    
    #create periodic boundary conditions
    iby = np.mod(iay,ny) #+1
    
    #fraction of atom at iay
    fay = 1-np.mod((ay/dy),1)
    
    #Add each atom to the potential grid j is too large to makegrid(iax,iay) which would allow us to vectorize V
    V1 = fax*fay*(Zatom**zed)
    V2 = (1-fax)*fay*(Zatom**zed)
    V3 = fax*(1-fay)*(Zatom**zed)
    V4 = (1-fax)*(1-fay)*(Zatom**zed);
    
    for a0 in range(0,amax):
        V[iax[a0],iay[a0]] = V[iax[a0],iay[a0]] + V1[a0]
        V[ibx[a0],iay[a0]] = V[ibx[a0],iay[a0]] + V2[a0]
        V[iax[a0],iby[a0]] = V[iax[a0],iby[a0]] + V3[a0]
        V[ibx[a0],iby[a0]] = V[ibx[a0],iby[a0]] + V4[a0]
        
    
    return V
    