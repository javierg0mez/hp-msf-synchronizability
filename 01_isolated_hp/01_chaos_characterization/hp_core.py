import numpy as np
from numba import njit
from scipy.stats import t as student_t


# Numerical tools for the isolated Hastings-Powell system
# parameters = [a1, a2, b1, b2, d1, d2]


# HP equations

@njit
def hp_rhs(state, parameters):

    x, y, z = state
    a1, a2, b1, b2, d1, d2 = parameters

    dx = x * (
        1.0
        - x
        - (a1 * y) / (1.0 + b1 * x)
    )

    dy = y * (
        (a1 * x) / (1.0 + b1 * x)
        - (a2 * z) / (1.0 + b2 * y)
        - d1
    )

    dz = z * (
        (a2 * y) / (1.0 + b2 * y)
        - d2
    )

    return np.array([dx, dy, dz])


# Jacobian of the isolated system

@njit
def hp_jacobian(state, parameters):

    x, y, z = state
    a1, a2, b1, b2, d1, d2 = parameters

    A = 1.0 + b1 * x
    B = 1.0 + b2 * y

    J11 = 1.0 - 2.0*x - (a1*y) / A**2
    J12 = -(a1*x) / A
    J13 = 0.0

    J21 = (a1*y) / A**2
    J22 = (a1*x) / A - d1 - (a2*z) / B**2
    J23 = -(a2*y) / B

    J31 = 0.0
    J32 = (a2*z) / B**2
    J33 = (a2*y) / B - d2

    return np.array([
        [J11, J12, J13],
        [J21, J22, J23],
        [J31, J32, J33]
    ])


# One RK4 step for the HP system

@njit
def rk4_step(state, parameters, dt):

    k1 = hp_rhs(state, parameters)

    k2 = hp_rhs(
        state + 0.5 * dt * k1,
        parameters
    )

    k3 = hp_rhs(
        state + 0.5 * dt * k2,
        parameters
    )

    k4 = hp_rhs(
        state + dt * k3,
        parameters
    )

    return state + (dt / 6.0) * (
        k1 + 2.0*k2 + 2.0*k3 + k4
    )


# Tangent matrix dynamics

@njit
def tangent_matrix_rhs(state, tangent_matrix, parameters):

    J = hp_jacobian(state, parameters)

    return J @ tangent_matrix


@njit
def rk4_tangent_matrix_step(
    state,
    tangent_matrix,
    parameters,
    dt
):

    k1_state = hp_rhs(state, parameters)

    k1_tangent = tangent_matrix_rhs(
        state,
        tangent_matrix,
        parameters
    )


    state_k2 = state + 0.5 * dt * k1_state
    tangent_k2 = tangent_matrix + 0.5 * dt * k1_tangent

    k2_state = hp_rhs(
        state_k2,
        parameters
    )

    k2_tangent = tangent_matrix_rhs(
        state_k2,
        tangent_k2,
        parameters
    )


    state_k3 = state + 0.5 * dt * k2_state
    tangent_k3 = tangent_matrix + 0.5 * dt * k2_tangent

    k3_state = hp_rhs(
        state_k3,
        parameters
    )

    k3_tangent = tangent_matrix_rhs(
        state_k3,
        tangent_k3,
        parameters
    )


    state_k4 = state + dt * k3_state
    tangent_k4 = tangent_matrix + dt * k3_tangent

    k4_state = hp_rhs(
        state_k4,
        parameters
    )

    k4_tangent = tangent_matrix_rhs(
        state_k4,
        tangent_k4,
        parameters
    )


    new_state = state + (dt / 6.0) * (
        k1_state
        + 2.0*k2_state
        + 2.0*k3_state
        + k4_state
    )

    new_tangent_matrix = tangent_matrix + (dt / 6.0) * (
        k1_tangent
        + 2.0*k2_tangent
        + 2.0*k3_tangent
        + k4_tangent
    )

    return new_state, new_tangent_matrix


# Gram-Schmidt for the three tangent directions

@njit
def gram_schmidt(tangent_matrix):

    q1 = tangent_matrix[:, 0].copy()

    norm1 = np.linalg.norm(q1)
    q1 = q1 / norm1


    q2 = tangent_matrix[:, 1].copy()

    q2 = q2 - np.dot(q1, q2) * q1

    norm2 = np.linalg.norm(q2)
    q2 = q2 / norm2


    q3 = tangent_matrix[:, 2].copy()

    q3 = q3 - np.dot(q1, q3) * q1
    q3 = q3 - np.dot(q2, q3) * q2

    norm3 = np.linalg.norm(q3)
    q3 = q3 / norm3


    Q = np.column_stack((
        q1,
        q2,
        q3
    ))

    growth = np.array([
        norm1,
        norm2,
        norm3
    ])

    return Q, growth


# Lyapunov spectrum computed in blocks

@njit
def lyapunov_spectrum_blocks(
    parameters,
    initial_state,
    dt,
    transient_time,
    alignment_time,
    average_time,
    renormalization_time,
    block_time
):

    state = initial_state.copy()


    # Remove the transient

    n_transient = int(
        transient_time / dt
    )

    for _ in range(n_transient):

        state = rk4_step(
            state,
            parameters,
            dt
        )


    # Start with three independent tangent directions

    tangent_matrix = np.eye(3)

    steps_per_renormalization = int(
        renormalization_time / dt
    )


    # Let the tangent directions align

    n_alignment = int(
        alignment_time / renormalization_time
    )

    for _ in range(n_alignment):

        for _ in range(
            steps_per_renormalization
        ):

            state, tangent_matrix = rk4_tangent_matrix_step(
                state,
                tangent_matrix,
                parameters,
                dt
            )

        tangent_matrix, growth = gram_schmidt(
            tangent_matrix
        )


    # Compute the spectrum in blocks

    n_blocks = int(
        average_time / block_time
    )

    renormalizations_per_block = int(
        block_time / renormalization_time
    )

    block_exponents = np.zeros(
        (n_blocks, 3)
    )

    trace_sum = 0.0
    trace_count = 0


    for block in range(n_blocks):

        log_growth = np.zeros(3)

        for _ in range(
            renormalizations_per_block
        ):

            for _ in range(
                steps_per_renormalization
            ):

                J = hp_jacobian(
                    state,
                    parameters
                )

                trace_sum += (
                    J[0, 0]
                    + J[1, 1]
                    + J[2, 2]
                )

                trace_count += 1


                state, tangent_matrix = rk4_tangent_matrix_step(
                    state,
                    tangent_matrix,
                    parameters,
                    dt
                )


            tangent_matrix, growth = gram_schmidt(
                tangent_matrix
            )

            log_growth += np.log(
                growth
            )


        block_exponents[block] = (
            log_growth / block_time
        )


    mean_trace = (
        trace_sum / trace_count
    )

    return (
        block_exponents,
        state,
        mean_trace
    )


# Keep the last successive maxima of z

@njit
def z_maxima(
    parameters,
    initial_state,
    dt,
    transient_time,
    record_time,
    max_peaks
):

    state = initial_state.copy()


    # Remove the transient

    n_transient = int(
        transient_time / dt
    )

    for _ in range(n_transient):

        state = rk4_step(
            state,
            parameters,
            dt
        )


    # Start reading z

    z_old = state[2]

    state = rk4_step(
        state,
        parameters,
        dt
    )

    z_current = state[2]


    peaks = np.empty(max_peaks)
    peak_count = 0

    n_record = int(
        record_time / dt
    )


    for _ in range(n_record):

        state = rk4_step(
            state,
            parameters,
            dt
        )

        z_new = state[2]


        if (
            z_old < z_current
            and z_current >= z_new
        ):

            peaks[
                peak_count % max_peaks
            ] = z_current

            peak_count += 1


        z_old = z_current
        z_current = z_new


    # Return the last maxima in chronological order

    if peak_count <= max_peaks:

        return (
            peaks[:peak_count].copy(),
            state
        )


    first = peak_count % max_peaks

    ordered_peaks = np.concatenate((
        peaks[first:],
        peaks[:first]
    ))

    return (
        ordered_peaks,
        state
    )


# Statistics for the Lyapunov spectrum

def spectrum_statistics(block_exponents):

    spectrum = np.mean(block_exponents, axis=0)
    error_estimates = []

    # Use blocks grouped by 1, 2 and 4
    for scale in (1, 2, 4):

        n_groups = len(block_exponents) // scale

        if n_groups < 8:
            continue

        grouped = block_exponents[:n_groups * scale].reshape(
            n_groups, scale, 3
        ).mean(axis=1)

        standard_error = np.std(
            grouped,
            axis=0,
            ddof=1
        ) / np.sqrt(n_groups)

        error = student_t.ppf(
            0.975,
            n_groups - 1
        ) * standard_error

        error_estimates.append(error)

    errors = np.max(
        np.array(error_estimates),
        axis=0
    )

    # Simple stationarity check
    half = len(block_exponents) // 2

    first_half = np.mean(
        block_exponents[:half],
        axis=0
    )

    second_half = np.mean(
        block_exponents[half:],
        axis=0
    )

    drift = np.abs(
        first_half - second_half
    )

    stationary = drift <= errors

    return (
        spectrum,
        errors,
        first_half,
        second_half,
        drift,
        stationary
    )