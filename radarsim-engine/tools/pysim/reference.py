"""Reference point-target baseband simulator (pure Python/NumPy).

This is the executable specification for radarsim-engine's PointSimulator.
The model was recovered from:
  - the original pure-Python simulator in radarsimpy history (src/simpy.py
    @4fe7b70): amplitude chain and waveform-phase beat structure
  - the golden references in tests/test_module_sim_radar_ideal.py (current
    C++ engine behavior: linear extrapolation of the waveform beyond the
    sampled grid, per-leg free-space spreading, center-frequency wavelength)

Signal model per (frame, tx m, rx n, pulse p, sample s):
  T     = timestamp[ch, p, s]                 (absolute, s)
  pos   = target position at T                (static + speed*T or arrays)
  R_tx  = |pos - tx_pos_world(T)|,  R_rx = |pos - rx_pos_world(T)|
  tau   = (R_tx + R_rx) / c
  u     = pulse-local waveform time = T - frame_start - pulse_start[p] - delay[m]
  beat  = phi_p(u) - phi_p(u - tau)           (cycles; phi integrates f_offset[p] + f(u))
  amp   = radar-equation voltage (see _amplitude), patterns in dB
  bb   += amp * mod(u) * pulse_mod[p] * exp(j*(2*pi*beat + phase_t))
"""

import numpy as np

C = 299792458.0


def _waveform_phase_fn(f, t):
    """Cumulative phase (cycles) of the sampled waveform f(t).

    Linear interpolation of f between samples -> quadratic phase.
    Beyond the grid, f is extended linearly with the end-segment slope
    (required to match the engine's golden references; the historical
    Python implementation clamped instead).
    """
    f = np.asarray(f, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    k = np.diff(f) / np.diff(t)  # per-segment frequency slope (Hz/s)
    dt = np.diff(t)
    # cumulative cycles at grid points; phi(t[0]) = 0
    phi_grid = np.concatenate(([0.0], np.cumsum(f[:-1] * dt + 0.5 * k * dt**2)))

    def phi(u):
        u = np.asarray(u, dtype=np.float64)
        out = np.empty_like(u)
        flat_u = u.ravel()
        flat_out = out.ravel()
        # segment index: last segment for u beyond the end, first for u < t[0]
        idx = np.searchsorted(t, flat_u, side="right") - 1
        idx = np.clip(idx, 0, len(t) - 2)
        du = flat_u - t[idx]
        flat_out[:] = (
            phi_grid[idx] + f[idx] * du + 0.5 * k[idx] * du**2
        )
        return out

    return phi


def _rot_matrix(yaw_pitch_roll_deg):
    """R = Rz(yaw) Ry(-pitch) Rx(roll), angles in degrees (radar convention)."""
    yaw, pitch, roll = np.radians(np.asarray(yaw_pitch_roll_deg, dtype=float))
    cy, sy = np.cos(yaw), np.sin(yaw)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cr, sr = np.cos(roll), np.sin(roll)
    return np.array(
        [
            [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
            [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
            [-sp, cp * sr, cp * cr],
        ]
    )


def _pattern_gain_db(az_deg, el_deg, ch):
    """Antenna gain (dB): NEAREST table entry on both axes + gain.

    Golden references (test_simc_tx_az_pattern / test_simc_tx_el_pattern)
    rule out linear interpolation: at 45 deg with a [-46, 0, 46] table the
    engine applies the -46 entry, i.e. nearest-neighbor lookup.
    """
    az_patterns = np.asarray(ch["az_patterns"])
    el_patterns = np.asarray(ch["el_patterns"])
    ia = np.argmin(np.abs(np.asarray(ch["az_angles"]) - az_deg))
    ie = np.argmin(np.abs(np.asarray(ch["el_angles"]) - el_deg))
    return az_patterns[ia] + el_patterns[ie] + ch["antenna_gain"]


def _channel_view_angles(direction_body):
    """Azimuth/elevation (deg) of a unit direction in the channel frame.

    az: 0 at +x, +90 at +y. el: above the x-y plane (el = 90 - theta).
    """
    az = np.degrees(np.arctan2(direction_body[..., 1], direction_body[..., 0]))
    el = np.degrees(np.arcsin(np.clip(direction_body[..., 2], -1.0, 1.0)))
    return az, el


def _platform_pose(radar, T):
    """Platform location and rotation matrix at absolute time T."""
    rp = radar.radar_prop
    loc = np.asarray(rp["location"], dtype=float)
    rot = np.asarray(rp["rotation"], dtype=float)
    speed = np.asarray(rp["speed"], dtype=float)
    rate = np.asarray(rp["rotation_rate"], dtype=float)
    if loc.ndim == 1:
        return loc + speed * T, _rot_matrix(rot + rate * T)
    raise NotImplementedError("time-varying platform motion not yet supported")


def _target_state(target, T, ch_idx):
    """Target position (m), rcs (dBsm), phase (deg) at time T.

    Time-varying arrays follow the timestamp shape [channels, pulses, samples].
    """
    loc = np.asarray(target["location"], dtype=float)
    speed = np.asarray(target.get("speed", (0, 0, 0)), dtype=float)
    if loc.ndim == 1:
        pos = loc + speed * T
    else:
        pos = loc[ch_idx]  # caller indexes [ch, p, s] externally
    rcs = np.asarray(target.get("rcs", 0.0), dtype=float)
    phase = np.asarray(target.get("phase", 0.0), dtype=float)
    return pos, rcs, phase


def sim_interference_reference(radar, interf):
    """Radar-to-radar interference (InterferenceSimulator reference).

    Recovered from test_simc_interference and gen_docs/user_guide/interference.rst:
      - one-way Friis path with an extra 1/(4*pi): P_r = P_t G lambda^2 / ((4 pi)^3 R^2)
      - the victim demodulates with its own LO: phase = phi_v(u_v) - phi_i(u_i)
        (same "local minus emission" structure as the target beat)
      - hard passband gate: |f_v(u_v) - f_i(u_i)| > noise_bandwidth -> no contribution
      - the interferer's waveform/pulse modulation applies; the victim's does not gate
    """
    tx_v = radar.radar_prop["transmitter"]
    rx_v = radar.radar_prop["receiver"]
    tx_i = interf.radar_prop["transmitter"]
    bb = rx_v.bb_prop

    f_v = np.asarray(tx_v.waveform_prop["f"], dtype=float)
    t_v = np.asarray(tx_v.waveform_prop["t"], dtype=float)
    foff_v = np.asarray(tx_v.waveform_prop["f_offset"], dtype=float)
    pstart_v = np.asarray(tx_v.waveform_prop["pulse_start_time"], dtype=float)

    f_i = np.asarray(tx_i.waveform_prop["f"], dtype=float)
    t_i = np.asarray(tx_i.waveform_prop["t"], dtype=float)
    foff_i = np.asarray(tx_i.waveform_prop["f_offset"], dtype=float)
    pstart_i = np.asarray(tx_i.waveform_prop["pulse_start_time"], dtype=float)
    pulses_i = int(tx_i.waveform_prop["pulses"])
    t_i_end = t_i[-1]

    phi_v_base = _waveform_phase_fn(f_v, t_v)
    phi_i_base = _waveform_phase_fn(f_i, t_i)
    fi_of_u = _freq_fn(f_i, t_i)
    fv_of_u = _freq_fn(f_v, t_v)

    fc_i = 0.5 * (f_i.min() + f_i.max())

    timestamp = radar.time_prop["timestamp"]
    n_ch_total, pulses, samples = timestamp.shape
    n_frames = np.atleast_1d(radar.time_prop["frame_start_time"]).size
    n_ch = n_ch_total // n_frames
    n_tx_v = int(tx_v.txchannel_prop["size"])
    n_rx_v = int(rx_v.rxchannel_prop["size"])
    n_tx_i = int(tx_i.txchannel_prop["size"])
    frame_start = np.atleast_1d(
        np.asarray(radar.time_prop["frame_start_time"], dtype=float)
    )
    frame_start_i = np.atleast_1d(
        np.asarray(interf.time_prop["frame_start_time"], dtype=float)
    )

    noise_bw = bb["noise_bandwidth"]
    out = np.zeros_like(timestamp, dtype=np.complex128)

    for ch in range(n_ch_total):
        f_idx = ch // n_ch
        ch_in_frame = ch % n_ch
        m = ch_in_frame // n_rx_v
        n = ch_in_frame % n_rx_v
        for mi in range(n_tx_i):
            mod_i = tx_i.txchannel_prop["waveform_mod"][mi]
            for p in range(pulses):
                for s in range(samples):
                    T = timestamp[ch, p, s]
                    u_v = (
                        T
                        - frame_start[f_idx]
                        - pstart_v[p]
                        - tx_v.txchannel_prop["delay"][m]
                    )
                    plat_v, rot_v = _platform_pose(radar, T)
                    rx_pos = plat_v + rot_v @ np.asarray(
                        rx_v.rxchannel_prop["locations"][n], dtype=float
                    )
                    for q in range(pulses_i):
                        plat_i, rot_i = _platform_pose(interf, T)
                        tx_pos_i = plat_i + rot_i @ np.asarray(
                            tx_i.txchannel_prop["locations"][mi], dtype=float
                        )
                        R = np.linalg.norm(tx_pos_i - rx_pos)
                        if R == 0.0:
                            continue
                        tau = R / C
                        u_i = (
                            T
                            - tau
                            - frame_start_i[min(f_idx, frame_start_i.size - 1)]
                            - pstart_i[q]
                            - tx_i.txchannel_prop["delay"][mi]
                        )
                        if u_i < t_i[0] or u_i > t_i_end:
                            continue
                        beat_f = (fv_of_u(u_v) + foff_v[p]) - (
                            fi_of_u(u_i) + foff_i[q]
                        )
                        if abs(beat_f) > noise_bw:
                            continue

                        # antenna gains: interferer tx toward victim, victim rx
                        # toward interferer, each in its own body frame
                        d_i = rx_pos - tx_pos_i
                        az_i, el_i = _channel_view_angles(rot_i.T @ (d_i / R))
                        d_v = tx_pos_i - rx_pos
                        az_r, el_r = _channel_view_angles(rot_v.T @ (d_v / R))
                        g_db = _pattern_gain_db(
                            az_i, el_i, _ch(tx_i.txchannel_prop, mi)
                        ) + _pattern_gain_db(
                            az_r, el_r, _ch(rx_v.rxchannel_prop, n)
                        )

                        p_w = (
                            1e-3
                            * 10 ** (tx_i.rf_prop["tx_power"] / 10)
                            * 10 ** (g_db / 10)
                            * (C / (fc_i + foff_i[q])) ** 2
                            / ((4 * np.pi) ** 3 * R**2)
                        )
                        amp = (
                            np.sqrt(2 * p_w * bb["load_resistor"])
                            * 10 ** (rx_v.rf_prop["rf_gain"] / 20)
                            * 10 ** (bb["baseband_gain"] / 20)
                        )
                        ph = 2 * np.pi * (
                            (foff_v[p] * u_v + phi_v_base(u_v))
                            - (foff_i[q] * u_i + phi_i_base(u_i))
                        )
                        pol = np.abs(
                            np.vdot(
                                np.asarray(rx_v.rxchannel_prop["polarization"][n]),
                                np.asarray(
                                    tx_i.txchannel_prop["polarization"][mi]
                                ),
                            )
                        )
                        val = amp * pol * np.exp(1j * ph)
                        if mod_i["enabled"]:
                            val *= _mod_zoh(mod_i, u_i)
                        val *= tx_i.txchannel_prop["pulse_mod"][mi][q]
                        out[ch, p, s] += val
    return out


def sim_radar_reference(radar, targets, interf=None, **_kwargs):
    """Reference implementation of sim_radar (point targets only)."""
    tx = radar.radar_prop["transmitter"]
    rx = radar.radar_prop["receiver"]
    wf = tx.waveform_prop
    bb = rx.bb_prop

    f = np.asarray(wf["f"], dtype=float)
    t = np.asarray(wf["t"], dtype=float)
    f_offset = np.asarray(wf["f_offset"], dtype=float)
    pulse_start = np.asarray(wf["pulse_start_time"], dtype=float)
    pulses = int(wf["pulses"])

    fs = float(bb["fs"])
    frame_start = np.atleast_1d(
        np.asarray(radar.time_prop["frame_start_time"], dtype=float)
    )
    n_frames = frame_start.size

    n_tx = int(tx.txchannel_prop["size"])
    n_rx = int(rx.rxchannel_prop["size"])
    n_ch = n_tx * n_rx
    samples = int(radar.sample_prop["samples_per_pulse"])

    timestamp = radar.time_prop["timestamp"]  # [F*M*N, pulses, samples]
    baseband = np.zeros((n_frames * n_ch, pulses, samples), dtype=np.complex128)

    fc_base = 0.5 * (f.min() + f.max())

    phi_base = _waveform_phase_fn(f, t)

    txch = tx.txchannel_prop
    rxch = rx.rxchannel_prop

    for f_idx in range(n_frames):
        for m in range(n_tx):
            tx_mod = txch["waveform_mod"][m]
            tx_pol = np.asarray(txch["polarization"][m])
            for n in range(n_rx):
                ch = f_idx * n_ch + m * n_rx + n
                rx_pol = np.asarray(rxch["polarization"][n])
                pol_factor = np.abs(np.vdot(rx_pol, tx_pol))

                for p in range(pulses):
                    f_off = f_offset[p]

                    def phi(u, _f_off=f_off):
                        return _f_off * np.asarray(u, dtype=float) + phi_base(u)

                    for s in range(samples):
                        T = timestamp[ch, p, s]
                        u = (
                            T
                            - frame_start[f_idx]
                            - pulse_start[p]
                            - txch["delay"][m]
                        )
                        plat_loc, plat_rot = _platform_pose(radar, T)
                        tx_pos = plat_loc + plat_rot @ np.asarray(
                            txch["locations"][m], dtype=float
                        )
                        rx_pos = plat_loc + plat_rot @ np.asarray(
                            rxch["locations"][n], dtype=float
                        )

                        for tgt in targets:
                            pos, rcs_db, phs_deg = _target_state_at(
                                tgt, T, ch, p, s
                            )
                            d_tx = pos - tx_pos
                            d_rx = pos - rx_pos
                            R_tx = np.linalg.norm(d_tx)
                            R_rx = np.linalg.norm(d_rx)
                            if R_tx == 0.0 or R_rx == 0.0:
                                continue
                            tau = (R_tx + R_rx) / C

                            # antenna gains (dB) in each channel's body frame
                            rot_inv = plat_rot.T
                            az_tx, el_tx = _channel_view_angles(
                                rot_inv @ (d_tx / R_tx)
                            )
                            az_rx, el_rx = _channel_view_angles(
                                rot_inv @ (d_rx / R_rx)
                            )
                            g_db = (
                                _pattern_gain_db(az_tx, el_tx, _ch(txch, m))
                                + _pattern_gain_db(az_rx, el_rx, _ch(rxch, n))
                            )

                            # wavelength rides the per-pulse frequency offset
                            fc = fc_base + f_off
                            pr_dbm = (
                                tx.rf_prop["tx_power"]
                                + g_db
                                - 10 * np.log10(4 * np.pi * R_tx**2)
                                + rcs_db
                                - 10 * np.log10(4 * np.pi * R_rx**2)
                                + 10 * np.log10((C / fc) ** 2 / (4 * np.pi))
                                + rx.rf_prop["rf_gain"]
                            )
                            amp = (
                                np.sqrt(
                                    1e-3 * 10 ** (pr_dbm / 10)
                                    * bb["load_resistor"]
                                )
                                * 10 ** (bb["baseband_gain"] / 20)
                                * np.sqrt(2)
                            )

                            beat = phi(u) - phi(u - tau)
                            val = (
                                amp
                                * pol_factor
                                * np.exp(
                                    1j * (2 * np.pi * beat + np.radians(phs_deg))
                                )
                            )
                            if tx_mod["enabled"]:
                                val *= _mod_zoh(tx_mod, u - tau)
                            val *= txch["pulse_mod"][m][p]
                            baseband[ch, p, s] += val

    interference = None
    if interf is not None:
        interference = sim_interference_reference(radar, interf)

    return {
        "baseband": baseband,
        "noise": None,
        "timestamp": timestamp,
        "interference": interference,
    }


def _freq_fn(f, t):
    """Linear-interpolated instantaneous frequency, clamped to endpoints."""
    kseg = np.diff(f) / np.diff(t)

    def f_of_u(u):
        u = float(u)
        idx = int(np.clip(np.searchsorted(t, u, side="right") - 1, 0, len(t) - 2))
        return f[idx] + kseg[idx] * (u - t[idx])

    return f_of_u


def _mod_zoh(mod, u_echo):
    """Waveform modulation: zero-order hold over the mod_t table, evaluated at
    the echo's transmit time.

    Recovered from the pulsed-radar and waveform-modulation goldens plus the
    upstream "modulation index fix" commit: the table index is
    ``floor((u_echo - mod_t[0]) / step) + 1`` with ``step = mod_t[1] -
    mod_t[0]``, and out-of-range indices contribute zero.
    """
    mt = np.asarray(mod["t"], dtype=float)
    var = np.asarray(mod["var"])
    step = mt[1] - mt[0]
    idx = int(np.floor((u_echo - mt[0]) / step)) + 1
    if 0 <= idx < len(var):
        return var[idx]
    return 0.0 + 0.0j


def _ch(channel_prop, idx):
    return {
        "az_angles": channel_prop["az_angles"][idx],
        "az_patterns": channel_prop["az_patterns"][idx],
        "el_angles": channel_prop["el_angles"][idx],
        "el_patterns": channel_prop["el_patterns"][idx],
        "antenna_gain": channel_prop["antenna_gains"][idx],
    }


def _target_state_at(target, T, ch, p, s):
    loc = np.asarray(target["location"], dtype=float)
    speed = np.asarray(target.get("speed", (0, 0, 0)), dtype=float)
    rcs = np.asarray(target.get("rcs", 0.0), dtype=float)
    phase = np.asarray(target.get("phase", 0.0), dtype=float)

    if loc.ndim == 1:
        pos = loc + speed * T
    else:
        pos = loc[ch, p, s]
    if rcs.ndim > 0:
        rcs = rcs[ch, p, s] if rcs.ndim == 3 else float(rcs)
    if phase.ndim > 0:
        phase = phase[ch, p, s] if phase.ndim == 3 else float(phase)
    return pos, float(rcs), float(phase)
