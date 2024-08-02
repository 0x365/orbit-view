import matplotlib.pyplot as plt
import numpy as np

from propagator import prop

plot_1_orbit = False

num_loops = 20
per_loop = 1000



def get_orbit(v_1,v_2, loops):
    v_1 = float(v_1)
    v_2 = float(v_2)
    m1,m2,m3 = 1,1,1
    r1x, r1y, r2x, r2y, r3x, r3y = -1, 0, 1, 0, 0, 0
    v1x, v1y, v2x, v2y, v3x, v3y = v_1, v_2, v_1, v_2, -2*v_1/m3, -2*v_2/m3
    initial = [m1,m2,m3,r1x,r1y,r2x,r2y,r3x,r3y,v1x,v1y,v2x,v2y,v3x,v3y]
    time = loops
    timestep = loops*per_loop
    return prop(initial, time=time, timestep=timestep)


def onclick(event):

    if event.xdata != None and event.ydata != None:
        counters = np.linspace(-1,1, np.shape(grid)[0])
        ii = np.argsort(np.square(counters - event.ydata))[0]
        jj = np.argsort(np.square(counters - event.xdata))[0]

        print("Plot", ii, jj)        

        fig2 = plt.figure(2, figsize=(8,5))
        fig2.clf()

        if np.isnan(period[ii,jj]) or np.isinf(period[ii,jj]) or not plot_1_orbit:
            ender = -1
            orbit_data = get_orbit(counters[ii],counters[jj], num_loops)
        else:
            ender = int(period[ii,jj]*per_loop)
            orbit_data = get_orbit(counters[ii],counters[jj], period[ii,jj])

        t = orbit_data["t"]
        r1x, r1y, r2x, r2y, r3x, r3y, _, _, _, _, _, _ = orbit_data["y"]

        plt.plot(r1x[:ender],r1y[:ender])
        plt.plot(r2x[:ender],r2y[:ender])
        plt.plot(r3x[:ender],r3y[:ender])

        plt.xlabel("X")
        plt.ylabel("Y")

        if len(t) != total_steps:
            collisional = "True"
            per = np.nan
            std = np.nan
        else:
            collisional = "False"
            if np.isnan(period[ii,jj]):
                per = np.inf
                std = np.inf
            else:
                per = period[ii,jj]
                std = 1/((1-grid[ii,jj])**2)
        
        textstr = ""
        textstr += "V1=%.9f\n"%(counters[ii])
        textstr += "V2=%.9f\n"%(counters[jj])
        textstr += "\n"
        textstr += "Predicted Period=%.3f\n"%(per)
        textstr += "STD in Euclidean Distance=%.0f\n"%(std)
        textstr += "\n"
        textstr += "Collisional Orbit="+collisional+"\n"
        textstr += "Total ODE Steps=%.0f\n"%(len(t))

        plt.text(-0.5, 0.5, textstr, va='center', ha="center", fontsize=10, transform=plt.gca().transAxes)
        plt.subplots_adjust(left=0.44)
        plt.title("Propagated Orbit for 1 of the predicted period")

        plt.axis('equal')
        fig2.canvas.draw()
        fig2.canvas.flush_events()
        plt.show()


def main_map():

    fig, axs = plt.subplots(1,1, figsize=(12,10))
    colors = plt.imshow(grid, extent=[-1,1,1,-1], cmap="BuPu_r")
    axs.invert_yaxis()

    ticker = np.arange(0, np.nanmax(grid), 0.01)
    cbar = plt.colorbar(colors, ticks=[*ticker, np.nanmax(grid)])
    cbar.ax.set_yticklabels([*np.array(np.around((1/(np.square(1-ticker))),decimals=-2), dtype=str), "inf"])

    plt.axis('equal')

    plt.xlabel("V2")
    plt.ylabel("V1")

    # plt.savefig("stability.png",dpi=1000)

    fig.canvas.mpl_connect('button_press_event',lambda event: onclick(event))

    plt.show()
    plt.clf()
    plt.close()



if __name__ == '__main__':
    total_steps = num_loops*per_loop

    grid = np.load("data/orbit_deviation_std.npy")
    grid[grid == np.inf] = np.nan
    grid = 1-np.sqrt(1/grid)

    period = np.load("data/periods.npy")

    main_map()
