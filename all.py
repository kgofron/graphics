# in shell: source /epics/iocs/ioc-hkl/iochkl/bin/activate
# export GI_TYPELIB_PATH=/usr/local/lib/girepository-1.0

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
from matplotlib import cm
from intensities2detint import intensities2detint

def rotation_matrix(axis, angle_rad):
    axis = np.array(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    C = 1 - c
    R = np.array([
        [c + x*x*C,     x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, c + y*y*C,     y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, c + z*z*C]
    ])
    return R

def partial_cylinder_intersection(ray_origin, ray_dir, R, theta_min, theta_max, height, cyl_center=(0,0)):
    ox, oy, oz = ray_origin
    dx, dy, dz = ray_dir
    x_c, y_c = cyl_center

    ox_shift = ox - x_c
    oy_shift = oy - y_c

    A = dx**2 + dy**2
    B = 2 * (ox_shift * dx + oy_shift * dy)
    C = ox_shift**2 + oy_shift**2 - R**2

    discriminant = B**2 - 4*A*C
    if discriminant < 0:
        return None

    t1 = (-B - np.sqrt(discriminant)) / (2*A)
    t2 = (-B + np.sqrt(discriminant)) / (2*A)

    t_candidates = [t for t in [t1, t2] if t > 0]
    if not t_candidates:
        return None

    t = max(t_candidates)
    hit_point = ray_origin + t * ray_dir

    z = hit_point[2]
    theta = np.arctan2(hit_point[1] - y_c, hit_point[0] - x_c)
    if theta < 0:
        theta += 2*np.pi

    theta_min = np.deg2rad(theta_min) % (2*np.pi)
    theta_max = np.deg2rad(theta_max) % (2*np.pi)
    if theta_min < theta_max:
        if (theta_min <= theta <= theta_max) and (-height <= z <= height):
            return hit_point
    else:
        if (theta >= theta_min or theta <= theta_max) and (-height <= z <= height):
            return hit_point
    return None

def main():
    #TODO feed values in from EPICS
    scatter_holder = {'scatter': None}
    hkl_path = '/epics/crystals/SiO2/EntryWithCollCode176.hkl'
    min_intensity = 1
    R = 10.0
    cyl_center = (4, 0)
    ray_origin = np.array([0.0, 0.0, 0.0])
    detect_height = 5
    detect_width = 10

    omega_axis = (0, 0, 1)
    chi_axis   = (1, 0, 0)
    phi_axis   = (0, 1, 0)
    theta_axis = (0, 0, 1)

    ray_dir_0 = np.array([1.0, 0.0, 0.0])

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(bottom=0.25)

    #TODO would need to change ax to ax3d, ax2d for flattened plot shown in window
    #fig = plt.figure(figsize=(10, 10))
    #ax3d = fig.add_subplot(211, projection='3d')  # top: 3D plot
    #ax2d = fig.add_subplot(212)                  # bottom: 2D intensity
    #plt.subplots_adjust(bottom=0.25, hspace=0.4)

    x_c, y_c = cyl_center
    height = detect_height
    box_size = 10

    theta_min_deg = -45
    theta_max_deg = 45

    theta_min = np.deg2rad(theta_min_deg)
    theta_max = np.deg2rad(theta_max_deg)

    # TODO calculate theta values from width value and radius
    xmin = 5
    ymin = y_c-detect_width/2 #TODO replace with theta
    ymax = y_c+detect_width/2 #TODO replace with theta
    zmin = -height
    zmax = height

    z = np.linspace(-height, height, 100)
    theta_grid = np.linspace(theta_min, theta_max, 50)
    theta_grid, z_grid = np.meshgrid(theta_grid, z)
    x_cyl = x_c + R * np.cos(theta_grid)
    y_cyl = y_c + R * np.sin(theta_grid)

    text1 = ax.text2D(-0.5, 1, f"beam direction: [{ray_dir_0[0]:.2f}, {ray_dir_0[1]:.2f}, {ray_dir_0[2]:.2f}]", transform=ax.transAxes, fontsize=10, color='blue')
    text2 = ax.text2D(-0.5, 0.95, f"sample position: {ray_origin}", transform=ax.transAxes, fontsize=10, color='blue')
    text3 = ax.text2D(-0.5, 0.85, f"detector radius: {R} \ndetector cylinder center: {cyl_center}", transform=ax.transAxes, fontsize=10, color='blue')
    text4 = ax.text2D(-0.5, 0.65, f"Rotation axes: \n omega {omega_axis} \n chi       {chi_axis} \n phi       {phi_axis} \n tth       {theta_axis}", transform=ax.transAxes, fontsize=10, color='blue')
    intersection_text = ax.text2D(-0.5, 0.5, "detector intersection: \nNone", transform=ax.transAxes, fontsize=10, color='black')

    ax.plot_surface(x_cyl, y_cyl, z_grid, alpha=0.3, color='cyan', edgecolor='none')

    ray_line, = ax.plot([], [], [], color='red', linewidth=2, label='Ray')
    intersection_dot = ax.scatter([], [], [], color='black', s=50, label='Intersection')

    ax.set_xlim(-box_size, box_size)
    ax.set_ylim(-box_size, box_size)
    ax.set_zlim(-box_size, box_size)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.title("reflection-detector intersection")

    orig_elev, orig_azim = ax.elev, ax.azim

    axcolor = 'lightgoldenrodyellow'
    ax_omega = plt.axes([0.2, 0.14, 0.65, 0.03], facecolor=axcolor)
    ax_chi   = plt.axes([0.2, 0.11, 0.65, 0.03], facecolor=axcolor)
    ax_phi   = plt.axes([0.2, 0.08, 0.65, 0.03], facecolor=axcolor)
    ax_theta = plt.axes([0.2, 0.05, 0.65, 0.03], facecolor=axcolor)

    s_omega = Slider(ax_omega, 'Omega', -180, 180, valinit=0)
    s_chi   = Slider(ax_chi,   'Chi',   -180, 180, valinit=45)
    s_phi   = Slider(ax_phi,   'Phi',   -180, 180, valinit=45)
    s_theta = Slider(ax_theta, 'Theta', -180, 180, valinit=0)

    def update(val):
        omega = np.deg2rad(s_omega.val)
        chi   = np.deg2rad(s_chi.val)
        phi   = np.deg2rad(s_phi.val)
        theta = np.deg2rad(s_theta.val)

        R_omega = rotation_matrix(omega_axis, omega)
        R_chi   = rotation_matrix(chi_axis, chi)
        R_phi   = rotation_matrix(phi_axis, phi)
        R_theta = rotation_matrix(theta_axis, theta)

        total_R = R_omega @ R_chi @ R_phi @ R_theta
        ray_dir_rotated = total_R @ ray_dir_0
        ray_dir_rotated /= np.linalg.norm(ray_dir_rotated)

        hit_point = partial_cylinder_intersection(ray_origin, ray_dir_rotated, R, cyl_center=cyl_center, theta_min=theta_min_deg, theta_max=theta_max_deg, height=height)

        t_vals = np.linspace(0, 10, 100)
        ray_points = ray_origin[:, None] + ray_dir_rotated[:, None] * t_vals
        ray_line.set_data(ray_points[0], ray_points[1])
        ray_line.set_3d_properties(ray_points[2])

        if hit_point is not None:
            intersection_dot._offsets3d = ([hit_point[0]], [hit_point[1]], [hit_point[2]])
            intersection_text.set_text(f"detector intersection: \n[{hit_point[0]:.2f}, {hit_point[1]:.2f}, {hit_point[2]:.2f}]")
        else:
            intersection_dot._offsets3d = ([], [], [])
            intersection_text.set_text("detector intersection: \nNone")
        fig.canvas.draw_idle()

    s_omega.on_changed(update)
    s_chi.on_changed(update)
    s_phi.on_changed(update)
    s_theta.on_changed(update)

    button_width = 0.12
    button_height = 0.05
    button_left = 0.85
    button_bottom_start = 0.7
    button_spacing = 0.07

    ax_orig = plt.axes([button_left, button_bottom_start + 3*button_spacing, button_width, button_height])
    ax_z    = plt.axes([button_left, button_bottom_start + 2*button_spacing, button_width, button_height])
    ax_x    = plt.axes([button_left, button_bottom_start + 1*button_spacing, button_width, button_height])
    ax_y    = plt.axes([button_left, button_bottom_start + 0*button_spacing, button_width, button_height])
    ax_run_intensities = plt.axes([button_left, button_bottom_start - 1*button_spacing, button_width, button_height])
    ax_save_png = plt.axes([button_left, button_bottom_start - 2*button_spacing, button_width, button_height])
    #ax_show_2d = plt.axes([button_left, button_bottom_start - 3*button_spacing, button_width, button_height])


    btn_orig = Button(ax_orig, 'Original View')
    btn_z = Button(ax_z, 'View Z')
    btn_x = Button(ax_x, 'View X')
    btn_y = Button(ax_y, 'View Y')
    btn_run_intensities = Button(ax_run_intensities, 'Run Intensities')
    btn_save_png = Button(ax_save_png, 'Save PNG')
    #btn_show_2d = Button(ax_show_2d, 'Show 2D Map')

    def view_z(event): ax.view_init(elev=90, azim=-90); plt.draw()
    def view_x(event): ax.view_init(elev=0, azim=0); plt.draw()
    def view_y(event): ax.view_init(elev=0, azim=90); plt.draw()
    def view_original(event): ax.view_init(elev=orig_elev, azim=orig_azim); plt.draw()

    btn_orig.on_clicked(view_original)
    btn_z.on_clicked(view_z)
    btn_x.on_clicked(view_x)
    btn_y.on_clicked(view_y)
    intensity_scatter = None

    def run_intensities(event):
        nonlocal intensity_scatter
        intensity_data = intensities2detint(hkl_path, min_intensity, R, cyl_center, ray_origin,xmin,ymin,ymax,zmin,zmax)
        xs, ys, zs, intensities = zip(*intensity_data)
        xs, ys, zs, intensities = map(np.array, (xs, ys, zs, intensities))
        if intensity_scatter is not None:
            intensity_scatter.remove()
        #sizes = 50 + 200 * (intensities / intensities.max())
        intensity_scatter = ax.scatter(
            xs, ys, zs,
            c=intensities,
            cmap='viridis',
            marker='o',
        #    s=50,
            alpha=0.8
        )
        if scatter_holder['scatter'] is not None:
            scatter_holder['scatter'].remove()
        scatter_holder['scatter'] = intensity_scatter
        fig.colorbar(intensity_scatter, ax=ax, shrink=0.5, aspect=10, label='Intensity')
        plt.draw()
    btn_run_intensities.on_clicked(run_intensities)
    
    def save_png(event):
        if scatter_holder['scatter'] is None:
            print("no scatter to save")
            return
        intensity_scatter = scatter_holder['scatter']
        intensities = intensity_scatter.get_array()
        xs, ys, zs = intensity_scatter._offsets3d
        thetas = np.arctan2(ys - y_c, xs - x_c)
        thetas = np.where(thetas < theta_min, thetas + 2*np.pi, thetas)
        n_theta = 360
        n_z = 200
        theta_grid = np.linspace(theta_min, theta_max, n_theta)
        z_grid = np.linspace(zmin, zmax, n_z)
        theta_idx = np.searchsorted(theta_grid, thetas)
        z_idx = np.searchsorted(z_grid, zs)
        img = np.zeros((n_z, n_theta))
        counts = np.zeros((n_z, n_theta))
        for i in range(len(intensities)):
            t_idx = np.clip(theta_idx[i], 0, n_theta-1)
            z_i = np.clip(z_idx[i], 0, n_z-1)
            img[z_i, t_idx] += intensities[i]
            counts[z_i, t_idx] += 1
        counts[counts==0] = 1
        img /= counts
        plt.figure(figsize=(8,4))
        plt.imshow(img, extent=[np.rad2deg(theta_min), np.rad2deg(theta_max), zmin, zmax], aspect='auto', origin='lower', cmap='viridis')
        plt.xlabel(r'$\theta$ (degrees)')
        plt.ylabel('z')
        plt.title('Detector intensities unwrapped')
        plt.colorbar(label='Intensity')
        plt.tight_layout()
        plt.savefig('detector_intensities.png', dpi=300)
        plt.close()
        print("saved: detector_intensities.png")
    btn_save_png.on_clicked(save_png)

    #def show_2d_map(event):
    #    nonlocal intensity_scatter
    #    if intensity_scatter is None:
    #        print("No intensities to display.")
    #        return
    #    xs, ys, zs = intensity_scatter._offsets3d
    #    xs = np.array(xs)
    #    ys = np.array(ys)
    #    zs = np.array(zs)
    #    thetas = np.arctan2(ys - y_c, xs - x_c)
    #    thetas = (thetas + 2*np.pi) % (2*np.pi)
    #    theta_grid = np.linspace(theta_min, theta_max, n_theta)
    #    z_grid = np.linspace(zmin, zmax, n_z)
    #    theta_idx = np.searchsorted(theta_grid, thetas)
    #    z_idx = np.searchsorted(z_grid, zs)
    #    img = np.zeros((n_z, n_theta))
    #    counts = np.zeros((n_z, n_theta))
    #    for i in range(len(thetas)):
    #        t_idx = np.clip(theta_idx[i], 0, n_theta-1)
    #        z_i = np.clip(z_idx[i], 0, n_z-1)
    #        img[z_i, t_idx] += intensities[i]
    #        counts[z_i, t_idx] += 1
    #    counts[counts==0] = 1
    #    img /= counts
    #    ax2d.clear()
    #    cax = ax2d.imshow(
    #        img,
    #        extent=[np.rad2deg(theta_min), np.rad2deg(theta_max), zmin, zmax],
    #        aspect='auto',
    #        origin='lower',
    #        cmap='viridis'
    #    )
    #    ax2d.set_xlabel(r'$\theta$ (degrees)')
    #    ax2d.set_ylabel('z')
    #    ax2d.set_title('Detector intensities unwrapped')
    #    fig.colorbar(cax, ax=ax2d, orientation='vertical', label='Intensity')
    #    plt.draw()
    #btn_show_2d.on_clicked(show_2d_map)
    update(None)
    plt.show()


if __name__ == "__main__":
    main()
