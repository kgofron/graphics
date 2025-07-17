import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox
from scipy.ndimage import gaussian_filter
from intensities2detint import intensities2detint
import sys

#TODO filter shown reflections by what is actually visible based on omega position 
# will need to treat omega and tth values separately, assuming omega is sample
# and tth is the point detector

cif_path = '/epics/crystals/SiO2/EntryWithCollCode176.cif'
hkl_path = '/epics/crystals/SiO2/EntryWithCollCode176.hkl'

R = 70 #cm
wavelength = 1.486 # Angs
cyl_center = (0,0)
ray_origin = np.array([0.0, 0.0, 0.0])
min_intensity = 100
window_width = 120  # degrees (horizontally from sample)
sigma_blur = 2

angle_deg = 12 # degrees (vertically from sample)
zmax = R * np.tan(np.deg2rad(angle_deg))
zmin = -zmax
y_range = 2*zmax

det_angle_deg = 7.5
det_zmax = R * np.tan(np.deg2rad(det_angle_deg))
det_zmin = -det_zmax
det_height = det_zmax-det_zmin

lst = intensities2detint(
    cif_path, hkl_path, wavelength, min_intensity,
    R, cyl_center, ray_origin, zmin, zmax
)

data = np.array(lst)
theta, z, intensity, h, k, l = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

nx, ny = 500, 200
theta_grid = np.linspace(0, 360, nx, endpoint=False)
z_grid = np.linspace(zmin, zmax, ny)

def compute_heatmap(threshold):
    heatmap = np.zeros((ny, nx))
    for t, zz, inten in zip(theta, z, intensity):
        if inten < threshold:
            continue
        i = int(nx*t/360) % nx
        j = np.searchsorted(z_grid, zz)
        if 0 <= j < ny:
            heatmap[j, i] += inten
    blurred = gaussian_filter(heatmap, sigma=sigma_blur)
    blurred /= blurred.max()
    return blurred

heatmap_blurred = compute_heatmap(min_intensity)

#plt.imshow(heatmap_blurred, cmap='hot')
#plt.savefig('test.png')
#sys.exit()

# interactive plot
fig, ax = plt.subplots(figsize=(10, 5))
plt.subplots_adjust(bottom=0.4)

initial_omega = 0
initial_y_offset = 0

def slice_window(hmap, omega_start, y_offset, det_zmin, det_zmax,
                 nx, ny, y_range, window_width, det_height):
    z_start = det_zmin + y_offset
    z_end = z_start + (det_zmax - det_zmin)
    i_start = int(nx * omega_start / 360)
    i_end   = (i_start + int(nx * window_width / 360)) % nx
    j_start = int(ny*z_start/y_range) + int(ny/2)
    j_end   = j_start + int(ny * det_height / y_range)
    if i_start < i_end:
        hmap_theta = hmap[:, i_start:i_end]
    else:
        hmap_theta = np.hstack((hmap[:, i_start:], hmap[:, :i_end]))
    hmap_window = hmap_theta[j_start:j_end, :]
    theta_extent = [omega_start, omega_start + window_width]
    z_extent = [z_start, z_end]
    return hmap_window, theta_extent, z_extent

heatmap_window, theta_extent, z_extent = slice_window(heatmap_blurred, \
    initial_omega, initial_y_offset, det_zmin, det_zmax, nx, ny, y_range, \
    window_width, det_height)

heatmap_im = ax.imshow(heatmap_window, \
            extent=[theta_extent[0], theta_extent[1], z_extent[0], z_extent[1]], \
            origin='lower', aspect='auto', cmap='viridis', vmin=0, vmax=1)

plt.colorbar(heatmap_im, ax=ax, label="Normalized Intensity")
ax.set_xlabel("ω (degrees)")
ax.set_ylabel("z")
#title = ax.set_title("Ω window: 0°–120°")

# sliders
axcolor = 'lightgoldenrodyellow'
ax_slider_omega = plt.axes([0.2, 0.25, 0.65, 0.03], facecolor=axcolor)
slider_omega = Slider(ax_slider_omega, 'omega start (°)', 0, 240, valinit=initial_omega)

ax_slider_y = plt.axes([0.2, 0.2, 0.65, 0.03], facecolor=axcolor)
slider_y = Slider(ax_slider_y, 'z offset', -5, 5, valinit=initial_y_offset)

ax_text_minint = plt.axes([0.2, 0.12, 0.1, 0.03], facecolor=axcolor)
text_minint = TextBox(ax_text_minint, 'min Intensity', initial=str(min_intensity))

def update(val):
    omega_start = slider_omega.val
    y_offset = slider_y.val
    try:
        current_min_intensity = float(text_minint.text)
    except ValueError:
        current_min_intensity = min_intensity

    global heatmap_blurred
    heatmap_blurred = compute_heatmap(current_min_intensity)

    #hmap_win, theta_extent, z_extent = slice_window(heatmap_blurred, omega_start, y_offset, det_zmin, det_zmax)

    hmap_win, theta_extent, z_extent = slice_window(heatmap_blurred, \
        omega_start, y_offset, det_zmin, det_zmax, nx, ny, y_range, \
        window_width, det_height)
    heatmap_im.set_data(hmap_win)
    heatmap_im.set_extent([theta_extent[0], theta_extent[1], z_extent[0], z_extent[1]])
    ax.set_xlim(*theta_extent)
    ax.set_ylim(*z_extent)
    fig.canvas.draw_idle()

slider_omega.on_changed(update)
slider_y.on_changed(update)
text_minint.on_submit(lambda text: update(None))
update(0)
# hover tooltip, info on all nearby peaks
annot = ax.annotate("", xy=(0,0), xytext=(20,20),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)

def on_hover(event):
    if event.inaxes != ax:
        annot.set_visible(False)
        fig.canvas.draw_idle()
        return

    x, y = event.xdata, event.ydata
    if x is None or y is None:
        annot.set_visible(False)
        fig.canvas.draw_idle()
        return

    omega_start = slider_omega.val
    y_offset = slider_y.val
    try:
        current_min_intensity = float(text_minint.text)
    except ValueError:
        current_min_intensity = min_intensity

    omega_end = (omega_start + window_width) % 360
    y_min = y_offset + det_zmin
    y_max = y_offset + det_zmax

    # mask by current window & min_intensity
    mask = (intensity >= current_min_intensity) & (z >= y_min) & (z <= y_max)

    # handle omega wrap-around
    if omega_start < omega_end:
        mask &= (theta >= omega_start) & (theta <= omega_end)
    else:
        mask &= (theta >= omega_start) | (theta <= omega_end)

    if not np.any(mask):
        annot.set_visible(False)
        fig.canvas.draw_idle()
        return

    theta_visible = theta[mask]
    z_visible = z[mask]
    h_visible = h[mask]
    k_visible = k[mask]
    l_visible = l[mask]

    distances = np.sqrt(((theta_visible - (x%360))/window_width)**2 + ((z_visible - y)/y_range)**2)
    close_idx = np.where(distances < 0.05)[0]  # hovering threshold

    if len(close_idx) > 0:
        annot.xy = (theta_visible[close_idx[0]], z_visible[close_idx[0]])
        lines = []
        for i in close_idx:
            lines.append(
                f"hkl=({int(h_visible[i])},{int(k_visible[i])},{int(l_visible[i])})"
                f" ω={theta_visible[i]:.1f}° z={z_visible[i]:.2f}"
            )
        annot.set_text('\n'.join(lines))
        annot.set_visible(True)
    else:
        annot.set_visible(False)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", on_hover)
plt.show()

