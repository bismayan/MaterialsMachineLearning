from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt

def get_sim(dt_frame,n_rows=2000,plt_flag=False,sort_flag=True,out_file="sim.png",plot_every=1):
    if sort_flag:
        dist=euclidean_distances(dt_frame.values,dt_frame.values[0])
        dt_temp=dt_frame.copy()
        dt_temp["dist"]=dist
        dt_sort=dt_temp.sort("dist").drop("dist",axis=1)
    else:
        dt_sort=dt_frame.copy()
    dist_full=euclidean_distances(dt_sort[0:n_rows].values)
    plt.figure()
    plt.imshow(dist_full[::plot_every,::plot_every],extent=(0,n_rows,n_rows,0))
    plt.colorbar()
    plt.savefig(out_file)
    if plt_flag:
        plt.show()
