# this function is to plot truncation for specific airfoil with indicated 
# indices, mid points(control points), unit normal and tangential vector 
def truncation_plotter(af):

    indices = af.indices_generate()

    n=len(indices)
    indices = np.array([[indices[i,0],indices[i,1],0] for i in range(0,n)])

   # define normal and tangential directions for each panel
    indices_2 = np.concatenate((indices,np.array([indices[0]])),axis=0)

    length_panel = np.linalg.norm(np.diff(indices_2,axis=0),axis=1)
    t_hat = np.dot(np.diag(1/length_panel),np.diff(indices_2,axis=0))

    z_hat = np.array([0,0,1])
    z_hat = np.tile(z_hat,(n,1))

    n_hat = np.cross(z_hat,t_hat)
   # define midpoints for each panel
    midpt = af.midpts(indices)
    return indices,midpt,n_hat,t_hat,length_panel

