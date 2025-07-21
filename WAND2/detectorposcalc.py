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

def real2det(gamma_axis, delta_axis, s_gamma, s_delta, R, cyl_center, ray_origin):
    # Rotation axes
    #gamma_axis   = (0, 0, 1)
    #delta_axis = (0, -1, 0)

    ray_dir_0 = np.array([1.0, 0.0, 0.0])
    x_c, y_c = cyl_center

    gamma = np.deg2rad(s_gamma)
    delta = np.deg2rad(s_delta)

    R_gamma = rotation_matrix(gamma_axis, gamma)
    R_delta = rotation_matrix(delta_axis, delta)

    total_R = R_delta @ R_gamma
    ray_dir_rotated = total_R @ ray_dir_0
    ray_dir_rotated /= np.linalg.norm(ray_dir_rotated)

    hit_point = cylinder_intersection(ray_origin, ray_dir_rotated, R, cyl_center=cyl_center)

    if hit_point is not None:
        x_hit, y_hit, z_hit = hit_point

        # compute unfolded coordinates:
        theta_rad = np.arctan2(y_hit - y_c, x_hit - x_c)
        if theta_rad < 0:
            theta_rad += 2*np.pi

        theta_deg = np.degrees(theta_rad)

        # return unfolded (theta_deg, z)
        return [theta_deg, z_hit]
    else:
        return None

