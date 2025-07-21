import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox
from scipy.ndimage import gaussian_filter
from intensities2detint import intensities2detint

cif_path = '/epics/crystals/SiO2/EntryWithCollCode176.cif'
hkl_path = '/epics/crystals/SiO2/EntryWithCollCode176.hkl'

gauss_sig = 5
gamma_axis = [0,0,1]
delta_axis = [0,-1,0]

R = 70
wavelength = 1.486
cyl_center = (0,0)
ray_origin = np.array([0.0, 0.0, 0.0])
min_intensity = 1
window_width = 120

angle_deg = 12
zmax = R * np.tan(np.deg2rad(angle_deg))
zmin = -zmax
y_range = 2*zmax

det_angle_deg = 7.5
det_zmax = R * np.tan(np.deg2rad(det_angle_deg))
det_zmin = -det_zmax
det_height = det_zmax-det_zmin

lst = intensities2detint(cif_path, hkl_path, wavelength, min_intensity, R,
    cyl_center, ray_origin, zmin, zmax, gamma_axis, delta_axis)

if lst != []:
    data = np.array(lst)
    theta, z, intensity, h, k, l, mu, omega, chi, phi, gamma, delta = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8], data[:,9], data[:,10], data[:,11]
    data_dict =  {'theta': theta, 'z': z, 'intensity':intensity, 'h': h, 'k':k, 'l':l, 'mu':mu, 'omega':omega, 'chi':chi, 'phi':phi, 'gamma':gamma, 'delta':delta}
    #df = pd.DataFrame(data_dict)
    #df.to_csv('refls.csv')
else:
    print("NO DATA")

nx, ny = 500, 200
theta_grid = np.linspace(0, 360, nx, endpoint=False)
z_grid = np.linspace(zmin, zmax, ny)

initial_y_offset = 0
initial_mu = -180
initial_mu_center = initial_mu + 60

def compute_heatmap(threshold, cur_mu): #TODO pass in theta, mu
    heatmap = np.zeros((ny, nx))
    for t, zz, inten, m in zip(theta, z, intensity, mu):
        if (inten > threshold) and ((cur_mu-2)<=m<=(cur_mu+2)):
            i = int(nx*t/360) % nx
            j = np.searchsorted(z_grid, zz)
            if 0 <= j < ny:
                heatmap[j, i] += inten
    blurred = gaussian_filter(heatmap, sigma=gauss_sig)
    if blurred.max() != 0:
        blurred /= blurred.max()
    return heatmap, blurred

heatmap, heatmap_blurred = compute_heatmap(min_intensity, initial_mu_center)

fig, ax = plt.subplots(figsize=(10, 5))
plt.subplots_adjust(bottom=0.4)


def slice_window(hmap, mu_start, y_offset, det_zmin, det_zmax,
                 nx, ny, y_range, window_width, det_height):
    z_start = det_zmin + y_offset
    z_end = z_start + (det_zmax - det_zmin)
    #i_start = int(nx * omega_start / 360)
    #i_end   = (i_start + int(nx * window_width / 360)) % nx
    j_start = int(ny*z_start/y_range) + int(ny/2)
    j_end   = j_start + int(ny * det_height / y_range)
    hmap_window = hmap[j_start:j_end, :]
    theta_extent = [0, 120] # detector doesn't rotate
    z_extent = [z_start, z_end]
    return hmap_window, theta_extent, z_extent

heatmap_window, theta_extent, z_extent = slice_window(heatmap_blurred,
    initial_mu, initial_y_offset, det_zmin, det_zmax, nx, ny, y_range,
    window_width, det_height)

heatmap_im = ax.imshow(heatmap_window,
            extent=[theta_extent[0], theta_extent[1], z_extent[0], z_extent[1]],
            origin='lower', aspect='auto', cmap='viridis', vmin=0, vmax=1)

plt.colorbar(heatmap_im, ax=ax, label="Normalized Intensity")
ax.set_xlabel("theta (°)")
ax.set_ylabel("z")

axcolor = 'lightgoldenrodyellow'
ax_slider_mu = plt.axes([0.2, 0.25, 0.65, 0.03], facecolor=axcolor)
slider_mu = Slider(ax_slider_mu, 'theta start (°)', -180, 60, valinit=initial_mu)

ax_slider_y = plt.axes([0.2, 0.2, 0.65, 0.03], facecolor=axcolor)
slider_y = Slider(ax_slider_y, 'z offset', -5, 5, valinit=initial_y_offset)

ax_text_minint = plt.axes([0.2, 0.12, 0.1, 0.03], facecolor=axcolor)
text_minint = TextBox(ax_text_minint, 'min Intensity', initial=str(min_intensity))

scatter = None  # placeholder for scatter plot

annot = ax.annotate("", xy=(0,0), xytext=(20,20),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)

def update(val):
    global heatmap_blurred, scatter
    mu_start = slider_mu.val
    cur_mu_center = mu_start + 60
    y_offset = slider_y.val
    try:
        current_min_intensity = float(text_minint.text)
    except ValueError:
        current_min_intensity = min_intensity

    heatmap, heatmap_blurred = compute_heatmap(current_min_intensity, cur_mu_center)

    hmap_win, theta_extent, z_extent = slice_window(heatmap_blurred,
        mu_start, y_offset, det_zmin, det_zmax, nx, ny, y_range,
        window_width, det_height)
    heatmap_im.set_data(hmap_win)
    heatmap_im.set_extent([theta_extent[0], theta_extent[1], z_extent[0], z_extent[1]])
    ax.set_xlim(*theta_extent)
    ax.set_ylim(*z_extent)

    # update scatter points for hover-over info
    #if scatter is not None:
    #    scatter.remove()
    #mask = ((z >= z_extent[0]) & (z <= z_extent[1]) & \
    #    (theta >= theta_extent[0]) & (theta <= theta_extent[1]))
    #    #(intensity >= current_min_intensity) &
    #    #(o >= cur_omega_center-2) & (o <= cur_omega_center+2) &
    #print(heatmap)
    #scatter_data = heatmap#np.column_stack([theta[mask], z[mask]])
    #scatter.set_offsets(scatter_data) if scatter else None
    #scatter = ax.scatter(scatter_data[:,0], scatter_data[:,1], s=20, facecolors='none', edgecolors='red', picker=True)

    fig.canvas.draw_idle()

slider_mu.on_changed(update)
slider_y.on_changed(update)
text_minint.on_submit(lambda text: update(None))
update(0)

#def hover(event):
#    if event.inaxes != ax or scatter is None:
#        annot.set_visible(False)
#        fig.canvas.draw_idle()
#        return
#    cont, ind = scatter.contains(event)
#    if cont:
#        index = ind["ind"][0]
#        x, y = scatter.get_offsets()[index]
#        annot.xy = (x, y)
#        # find the hkl info for this point
#        idx = np.where((theta==x) & (z==y))[0]
#        if len(idx):
#            i = idx[0]
#            annot.set_text(
#                f"hkl=({int(h[i])},{int(k[i])},{int(l[i])})\nω={theta[i]:.1f}° z={z[i]:.2f}"
#            )
#            annot.set_visible(True)
#    else:
#        annot.set_visible(False)
#    fig.canvas.draw_idle()

#fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()

