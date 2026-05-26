import numpy as np
import matplotlib.pyplot as plt



# Needed functions:
print(np.pi)
# T_from_DH_row

def T_from_DH_row(DH_row):
    # This function takes in a DH row (array with 4 numerical values) and returns the corresponding 4x4 transformation matrix.
    d, theta, a, alpha = DH_row
    ct, st = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([
        [ct, -st*ca,  st*sa, a*ct],
        [st,  ct*ca, -ct*sa, a*st],
        [0,   sa,     ca,    d   ],
        [0,   0,      0,     1   ],
    ])

def n_T_n_from_DH_table(DH_table, n_start, n_end):
    # This function takes in a DH table (2D array where each row is a DH row) and returns the overall transformation matrix from the base to the end-effector.
    T = np.eye(4)
    for i in range(n_start, n_end):
        T = T @ T_from_DH_row(DH_table[i])
    return T


def geometric_jacobian_column_from_DH_table(DH_table, n):
    # This function takes in a DH table and returns the geometric Jacobian column corresponding to the joint n.
    T_prev = n_T_n_from_DH_table(DH_table, 0, n)
    T_end  = n_T_n_from_DH_table(DH_table, 0, len(DH_table))

    z_prev = T_prev[:3, 2]
    p_prev = T_prev[:3, 3]
    p_end  = T_end[:3, 3]

    J_v = np.cross(z_prev, p_end - p_prev)
    J_w = z_prev
    return np.concatenate([J_v, J_w])

def geometric_jacobian_from_DH_table(DH_table):
    # This function takes in a DH table and returns the full geometric Jacobian matrix for the end-effector. It uses the geometric_jacobian_column_from_DH_table function to compute each column of the Jacobian.
    J = np.zeros((6, len(DH_table)))
    for n in range(len(DH_table)):
        J[:, n] = geometric_jacobian_column_from_DH_table(DH_table, n)
    return J

def plot_robot_from_DH_table(DH_table):
    # Frame origins: 0 = base, N = end-effector. Intermediate ones are the joint positions.
    origins = np.array([n_T_n_from_DH_table(DH_table, 0, i)[:3, 3]
                        for i in range(len(DH_table) + 1)])

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Robot arm: lines between consecutive frame origins.
    ax.plot(origins[:, 0], origins[:, 1], origins[:, 2],
            '-', color='black', linewidth=2, label='Links')

    # Joints (everything between base and end-effector) with a number label.
    for i in range(1, len(origins) - 1):
        ax.scatter(*origins[i], color='steelblue', s=60)
        ax.text(*origins[i], f' {i}', fontsize=11, color='steelblue')

    # Base and end-effector markers.
    ax.scatter(*origins[0],  color='black',   marker='s', s=120, label='Base')
    ax.scatter(*origins[-1], color='magenta', marker='*', s=200, label='End-effector')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Robot arm from DH table')
    ax.legend(loc='upper left')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()

def CLIK(p_d, R_d, omega_d, p_d_dot, q, q_dot, J, DH_func,
         K_P=10.0, K_O=10.0, K_N=0.5, dt=0.01, tol=1e-3, max_iter=5000):
    # Closed-Loop Inverse Kinematics, following the schematic block-for-block:
    # the desired feedforward velocities (p_d_dot, omega_d) are summed with
    # proportional corrections (K_P * e_P, K_O * e_O) on the position and
    # orientation errors, then mapped to joint velocities by the Jacobian
    # pseudo-inverse and integrated to update q.
    #
    # J        : callable q -> 6xN Jacobian matrix
    # DH_func  : callable q -> DH table (needed for p(q), R(q))
    q = np.asarray(q, dtype=float).copy()
    q_history = [q.copy()]

    for _ in range(max_iter):
        # Forward kinematics: end-effector position p(q) and orientation R(q).
        T_e = n_T_n_from_DH_table(DH_func(q), 0, len(q))
        p_e = T_e[:3, 3]
        R_e = T_e[:3, :3]

        # Position error e_P = p_d - p(q).
        e_P = p_d - p_e

        # Orientation error from R(theta, r) = R_d @ R_e^T, then e_O = r * theta.
        R_err = R_d @ R_e.T
        angle = np.arccos(np.clip((np.trace(R_err) - 1) / 2, -1, 1))
        if abs(angle) < 1e-9:
            e_O = np.zeros(3)
        else:
            axis = np.array([R_err[2, 1] - R_err[1, 2],
                             R_err[0, 2] - R_err[2, 0],
                             R_err[1, 0] - R_err[0, 1]]) / (2 * np.sin(angle))
            e_O = axis * angle

        # Reference velocity = feedforward + proportional correction.
        v_ref = np.concatenate([p_d_dot + K_P * e_P,
                                omega_d + K_O * e_O])

        # Joint velocity via Jacobian pseudo-inverse + null-space term that
        # pulls q toward zero in directions that don't affect the end-effector.
        J_q    = J(q)
        J_pinv = np.linalg.pinv(J_q)
        N      = np.eye(len(q)) - J_pinv @ J_q   # null-space projector
        q_dot  = J_pinv @ v_ref + N @ (-K_N * q)
        q = q + q_dot * dt
        q_history.append(q.copy())

        if np.linalg.norm(np.concatenate([e_P, e_O])) < tol:
            break

    return q, q_dot, np.array(q_history)


def plot_CLIK_convergence(q_history, DH_func):
    # Plots the robot arm at every iteration of CLIK: start in blue, final in
    # black with the magenta end-effector marker, all intermediates in gray.
    def origins(q):
        DH_table = DH_func(q)
        return np.array([n_T_n_from_DH_table(DH_table, 0, i)[:3, 3]
                         for i in range(len(DH_table) + 1)])

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Intermediate iterations in gray.
    for q in q_history[1:-1]:
        o = origins(q)
        ax.plot(o[:, 0], o[:, 1], o[:, 2], '-', color='lightgray',
                linewidth=1, alpha=0.4)

    # Start state.
    o_start = origins(q_history[0])
    ax.plot(o_start[:, 0], o_start[:, 1], o_start[:, 2],
            '-', color='steelblue', linewidth=2, label='Start')

    # Final state.
    o_end = origins(q_history[-1])
    ax.plot(o_end[:, 0], o_end[:, 1], o_end[:, 2],
            '-', color='black', linewidth=2, label='Final')

    # Base and end-effector markers (use final config for the base; it doesn't move).
    ax.scatter(*o_end[0],  color='black',   marker='s', s=120, label='Base')
    ax.scatter(*o_end[-1], color='magenta', marker='*', s=200, label='End-effector')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('CLIK convergence')
    ax.legend(loc='upper left')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()
