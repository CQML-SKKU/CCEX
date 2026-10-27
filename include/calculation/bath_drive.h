#ifndef __CCEX_BATH_DRIVE_H_
#define __CCEX_BATH_DRIVE_H_

/**
 * @file bath_drive.h
 * @brief Continuous driving of bath spins in the cCCE (pulse.bath_drive): rotating frame, rotating-wave approximation.
 *
 * Ported from pyCCEX ccex/simulator/drive.py (feature/spin-squeezing b3ce43e), the setting of Garcia Olivan et al.,
 * arXiv 2512.06948: the qubit gets instantaneous pi pulses (Ramsey, Hahn, CPMG), electron-like spin-1/2 bath spins on
 * the lines of a Defect table are driven by one or more tones. docs/methods/bath-drive.md has the model.
 */

#include <string>
#include <vector>

/** One drive tone: a field of Rabi frequency rabi on the listed lines, at line + offset, with phase. */
struct BathDriveTone {
    std::vector<double> lines;   ///< lines it targets [rad*kHz] (input MHz)
    double rabi = 0.0;           ///< Rabi frequency W [rad*kHz] (input MHz); H = W Jx gives Rabi oscillation at W
    double offset = 0.0;         ///< carrier offset from the line [rad*kHz] (input MHz)
    double phase = 0.0;          ///< phase of the field in the frame [rad] (input degrees)
};

/** pulse.bath_drive (NULL in Pulse: no drive). */
struct BathDrive {
    std::vector<BathDriveTone> tones;
    bool shifts = false;              ///< Bloch-Siegert, AC Stark, Van Vleck and drive-assisted cross-line flip-flop
    bool qubit_transverse = true;     ///< the qubit coupling A_zx, A_zy kept exactly (tilted branch frames)
    double line_tolerance = 0.0;      ///< a spin is on a line when closer than this [rad*kHz] (input MHz, default 1e-3)
};

#endif
