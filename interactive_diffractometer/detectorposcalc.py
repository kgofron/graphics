import numpy as np

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

def cylinder_intersection(ray_origin, ray_dir, R, cyl_center=(0,0)):
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
    return ray_origin + t * ray_dir

def real2det(s_omega, s_chi, s_phi, s_theta, R, cyl_center, ray_origin):
    # Parameters
    #R = 3.0
    #cyl_center = (4, 0)
    #ray_origin = np.array([0.0, 0.0, 0.0])

    # Rotation axes
    omega_axis = (0, 0, 1)
    chi_axis   = (1, 0, 0)
    phi_axis   = (0, 1, 0)
    theta_axis = (0, 0, 1)

    ray_dir_0 = np.array([1.0, 0.0, 0.0])

    x_c, y_c = cyl_center
    height = 5
    box_size = 5

    # Cylinder surface
    z = np.linspace(-height, height, 100)
    theta_grid = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid = np.meshgrid(theta_grid, z)
    x_cyl = x_c + R * np.cos(theta_grid)
    y_cyl = y_c + R * np.sin(theta_grid)

    omega = np.deg2rad(s_omega)
    chi   = np.deg2rad(s_chi)
    phi   = np.deg2rad(s_phi)
    theta = np.deg2rad(s_theta)

    R_omega = rotation_matrix(omega_axis, omega)
    R_chi   = rotation_matrix(chi_axis, chi)
    R_phi   = rotation_matrix(phi_axis, phi)
    R_theta = rotation_matrix(theta_axis, theta)

    #total_R = R_theta @ R_phi @ R_chi @ R_omega
    total_R = R_omega @ R_chi @ R_phi @ R_theta #TODO order affects result, 
    # a theta rotation is only correct if not applied first (last in order)
    # equivalent to total_R = R_omega @ (R_chi @ (R_phi @ (R_theta(beam))))
    # Now if omega has a value, chi, phi, theta are 0, rotating chi does nothing
    # does chi rotation need phi value?
    # Can I just multiple this to UB matrix?
    # Note, this is all for a made up orientation of omega, chi, phi, tth
    ray_dir_rotated = total_R @ ray_dir_0
    ray_dir_rotated /= np.linalg.norm(ray_dir_rotated)
    hit_point = cylinder_intersection(ray_origin, ray_dir_rotated, R, cyl_center=cyl_center)
    if hit_point is not None:
            #return([hit_point[0]], [hit_point[1]], [hit_point[2]])
            return([hit_point[0], hit_point[1], hit_point[2]])
    else:
            return None


