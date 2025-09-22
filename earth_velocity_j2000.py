#!/usr/bin/env python3
"""
Calculate Earth's velocity at J2000 epoch using astropy.

This script computes Earth's velocity vector in various coordinate systems
at the J2000.0 epoch (2000-01-01 12:00:00 TT).
"""

import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
from astropy import units as u
from astropy.coordinates import ICRS
import astropy.constants as const

def get_earth_velocity_j2000():
    """
    Get Earth's velocity at J2000.0 epoch in various units.
    
    Returns:
        dict: Dictionary containing velocity components and magnitudes
    """
    # Define J2000.0 epoch
    j2000 = Time('J2000.0', scale='tt')
    print(f"Epoch: {j2000.iso} (J2000.0)")
    print(f"Julian Date: {j2000.jd}")
    
    # Use high-precision ephemeris (DE440 if available, fallback to builtin)
    try:
        print("Attempting to use DE440 ephemeris...")
        with solar_system_ephemeris.set('de440'):
            # Get Earth's barycentric position and velocity
            earth_posvel = get_body_barycentric_posvel('earth', j2000)
            earth_pos = earth_posvel[0]
            earth_vel = earth_posvel[1]
            ephemeris_used = 'de440'
            print("Successfully using DE440 ephemeris")
    except Exception as e:
        print(f"DE440 ephemeris not available ({e}), using built-in ephemeris")
        with solar_system_ephemeris.set('builtin'):
            # Get Earth's barycentric position and velocity
            earth_posvel = get_body_barycentric_posvel('earth', j2000)
            earth_pos = earth_posvel[0]
            earth_vel = earth_posvel[1]
            ephemeris_used = 'builtin'
    
    print("\n" + "="*60)
    print("EARTH'S VELOCITY AT J2000.0 EPOCH")
    print(f"Using ephemeris: {ephemeris_used}")
    print("="*60)
    
    # Convert to different units
    vel_km_s = earth_vel.xyz.to(u.km/u.s)
    vel_m_s = earth_vel.xyz.to(u.m/u.s)
    vel_au_day = earth_vel.xyz.to(u.au/u.day)
    
    print(f"\nBarycentric Velocity Components:")
    print(f"X: {vel_km_s[0]:.6f}")
    print(f"Y: {vel_km_s[1]:.6f}")
    print(f"Z: {vel_km_s[2]:.6f}")
    
    print(f"\nBarycentric Velocity Components (m/s):")
    print(f"X: {vel_m_s[0]:.3f}")
    print(f"Y: {vel_m_s[1]:.3f}")
    print(f"Z: {vel_m_s[2]:.3f}")
    
    print(f"\nBarycentric Velocity Components (AU/day):")
    print(f"X: {vel_au_day[0]:.8f}")
    print(f"Y: {vel_au_day[1]:.8f}")
    print(f"Z: {vel_au_day[2]:.8f}")
    
    # Calculate velocity magnitude
    vel_magnitude_km_s = np.sqrt(vel_km_s[0]**2 + vel_km_s[1]**2 + vel_km_s[2]**2)
    vel_magnitude_m_s = np.sqrt(vel_m_s[0]**2 + vel_m_s[1]**2 + vel_m_s[2]**2)
    
    print(f"\nVelocity Magnitude:")
    print(f"{vel_magnitude_km_s:.6f}")
    print(f"{vel_magnitude_m_s:.3f}")
    
    # Express as fraction of speed of light
    c = const.c.to(u.km/u.s)
    vel_fraction_c = vel_magnitude_km_s / c
    print(f"Velocity as fraction of c: {vel_fraction_c:.2e}")
    
    # Also get position for context
    earth_pos_au = earth_pos.xyz.to(u.au)
    pos_magnitude_au = np.sqrt(earth_pos_au[0]**2 + earth_pos_au[1]**2 + earth_pos_au[2]**2)
    
    print(f"\nBarycentric Position (AU):")
    print(f"X: {earth_pos_au[0]:.6f}")
    print(f"Y: {earth_pos_au[1]:.6f}")
    print(f"Z: {earth_pos_au[2]:.6f}")
    print(f"Distance from barycenter: {pos_magnitude_au:.6f}")
    
    # Return data for further use
    return {
        'epoch': j2000,
        'velocity_km_s': vel_km_s,
        'velocity_m_s': vel_m_s,
        'velocity_au_day': vel_au_day,
        'velocity_magnitude_km_s': vel_magnitude_km_s,
        'velocity_magnitude_m_s': vel_magnitude_m_s,
        'position_au': earth_pos_au,
        'position_magnitude_au': pos_magnitude_au,
        'ephemeris_used': ephemeris_used
    }

def compare_with_orbital_velocity():
    """
    Compare with approximate orbital velocity around the Sun.
    """
    print("\n" + "="*60)
    print("COMPARISON WITH APPROXIMATE ORBITAL VELOCITY")
    print("="*60)
    
    # Approximate orbital velocity: v = 2π * r / T
    # where r ≈ 1 AU and T ≈ 1 year
    orbital_radius = 1.0 * u.au
    orbital_period = 1.0 * u.year
    
    approx_orbital_vel = (2 * np.pi * orbital_radius / orbital_period).to(u.km/u.s)
    print(f"Approximate circular orbital velocity: {approx_orbital_vel:.3f}")
    
    # More precise calculation using GM_sun
    GM_sun = const.GM_sun
    precise_orbital_vel = np.sqrt(GM_sun / orbital_radius).to(u.km/u.s)
    print(f"Precise circular orbital velocity: {precise_orbital_vel:.3f}")

def get_alpha_centauri_data_j2000():
    """
    Get Alpha Centauri position and velocity at J2000.0 epoch.
    Note: Alpha Centauri is not a solar system body, so we'll use astrometric data.
    """
    print("\n" + "="*60)
    print("ALPHA CENTAURI POSITION & VELOCITY CALCULATION")
    print("="*60)
    
    # Alpha Centauri astrometric data from literature
    # Kervella et al. 2017 and other sources
    
    # J2000.0 coordinates (ICRS)
    ra_j2000 = 219.90205833  # degrees (14h 39m 36.49s)
    dec_j2000 = -60.83399167  # degrees (-60° 50' 02.3")
    
    # Proper motion in mas/year (milliarcseconds per year)
    mu_alpha_cos_delta = -3678.19  # mas/yr (RA component, includes cos(dec) factor)
    mu_delta = 481.84              # mas/yr (Dec component)
    
    # Radial velocity in km/s
    v_radial = -21.6  # km/s (approaching us)
    
    # Distance to Alpha Centauri
    distance_pc = 1.3348  # parsecs (Kervella et al. 2017)
    distance_km = distance_pc * 3.0857e13  # km
    distance_au = distance_pc * 206265  # AU
    distance_m = distance_pc * 3.0857e16  # m
    
    print(f"Alpha Centauri System Parameters:")
    print(f"  RA (J2000): {ra_j2000:.8f}° = {ra_j2000/15:.6f}h")
    print(f"  Dec (J2000): {dec_j2000:.8f}°")
    print(f"  Distance: {distance_pc:.4f} pc = {distance_au:.1f} AU")
    print(f"  Proper motion (RA): {mu_alpha_cos_delta:.2f} mas/yr")
    print(f"  Proper motion (Dec): {mu_delta:.2f} mas/yr")
    print(f"  Radial velocity: {v_radial:.1f} km/s")
    print()
    
    # Convert spherical coordinates to Cartesian (ICRS)
    ra_rad = np.radians(ra_j2000)
    dec_rad = np.radians(dec_j2000)
    
    # Position vector in meters (ICRS)
    x_m = distance_m * np.cos(dec_rad) * np.cos(ra_rad)
    y_m = distance_m * np.cos(dec_rad) * np.sin(ra_rad)
    z_m = distance_m * np.sin(dec_rad)
    
    position_m = np.array([x_m, y_m, z_m])
    
    # Convert to parsecs for comparison
    position_pc = position_m / 3.0857e16
    
    # Convert to geometric units (solar mass units)
    Msol2m = 1476.67  # m (GM_sun/c^2)
    position_geom = position_m / Msol2m
    
    print(f"Position in Cartesian coordinates (ICRS):")
    print(f"  Position (m): [{x_m:.6e}, {y_m:.6e}, {z_m:.6e}]")
    print(f"  Position (pc): [{position_pc[0]:.6f}, {position_pc[1]:.6f}, {position_pc[2]:.6f}]")
    print(f"  Position (geom): [{position_geom[0]:.6e}, {position_geom[1]:.6e}, {position_geom[2]:.6e}]")
    print()
    
    # Convert proper motion to transverse velocity
    # v_transverse = 4.74 * mu * d, where mu is in mas/yr and d is in pc
    v_transverse_alpha = 4.74 * mu_alpha_cos_delta * distance_pc / 1000  # km/s
    v_transverse_delta = 4.74 * mu_delta * distance_pc / 1000  # km/s
    
    print(f"Velocity Components:")
    print(f"  Radial velocity: {v_radial:.3f} km/s")
    print(f"  Transverse velocity (RA): {v_transverse_alpha:.3f} km/s")
    print(f"  Transverse velocity (Dec): {v_transverse_delta:.3f} km/s")
    
    # Total velocity magnitude
    v_total = np.sqrt(v_radial**2 + v_transverse_alpha**2 + v_transverse_delta**2)
    print(f"  Total velocity magnitude: {v_total:.3f} km/s")
    
    # Convert velocity to Cartesian components (approximate)
    # This is a simplified conversion - full transformation would require more complex spherical coordinate derivatives
    velocity_km_s = np.array([v_radial, v_transverse_alpha, v_transverse_delta])
    
    # Convert to units of c
    c_km_s = const.c.to(u.km/u.s).value
    velocity_c = velocity_km_s / c_km_s
    
    print(f"\nAlpha Centauri velocity in units of c:")
    print(f"  [{velocity_c[0]:.8e}, {velocity_c[1]:.8e}, {velocity_c[2]:.8e}]")
    
    return {
        'position_m': position_m,
        'position_pc': position_pc,
        'position_geom': position_geom,
        'velocity_km_s': velocity_km_s,
        'velocity_c': velocity_c,
        'distance_pc': distance_pc,
        'distance_au': distance_au,
        'ra_j2000': ra_j2000,
        'dec_j2000': dec_j2000
    }

def output_for_julia():
    """
    Output velocity in format suitable for Julia/PoMiN.
    """
    data = get_earth_velocity_j2000()
    vel = data['velocity_km_s']
    
    print("\n" + "="*60)
    print("FOR JULIA/PoMiN - EARTH VELOCITY (km/s)")
    print("="*60)
    print(f"vEarth = [{vel[0].value:.6f}, {vel[1].value:.6f}, {vel[2].value:.6f}]")
    
    # Also in m/s
    vel_ms = data['velocity_m_s']
    print(f"\nFOR JULIA/PoMiN - EARTH VELOCITY (m/s)")
    print(f"vEarth = [{vel_ms[0].value:.3f}, {vel_ms[1].value:.3f}, {vel_ms[2].value:.3f}]")
    
    # In units of c
    c_km_s = const.c.to(u.km/u.s).value
    vel_c = vel / c_km_s
    print(f"\nFOR JULIA/PoMiN - EARTH VELOCITY (units of c)")
    print(f"vEarth = [{vel_c[0].value:.8e}, {vel_c[1].value:.8e}, {vel_c[2].value:.8e}]")
    
    # Alpha Centauri position and velocity
    alpha_data = get_alpha_centauri_data_j2000()
    alpha_pos_pc = alpha_data['position_pc']
    alpha_pos_geom = alpha_data['position_geom']
    alpha_vel = alpha_data['velocity_km_s']
    alpha_vel_c = alpha_data['velocity_c']
    
    print(f"\nFOR JULIA/PoMiN - ALPHA CENTAURI POSITION (pc)")
    print(f"qalphapc = [{alpha_pos_pc[0]:.6f}, {alpha_pos_pc[1]:.6f}, {alpha_pos_pc[2]:.6f}]")
    
    print(f"\nFOR JULIA/PoMiN - ALPHA CENTAURI POSITION (geometric units)")
    print(f"qalpha = [{alpha_pos_geom[0]:.6e}, {alpha_pos_geom[1]:.6e}, {alpha_pos_geom[2]:.6e}]")
    
    print(f"\nFOR JULIA/PoMiN - ALPHA CENTAURI VELOCITY (km/s)")
    print(f"vAlphaCen = [{alpha_vel[0]:.6f}, {alpha_vel[1]:.6f}, {alpha_vel[2]:.6f}]")
    
    print(f"\nFOR JULIA/PoMiN - ALPHA CENTAURI VELOCITY (units of c)")
    print(f"vAlphaCen = [{alpha_vel_c[0]:.8e}, {alpha_vel_c[1]:.8e}, {alpha_vel_c[2]:.8e}]")

if __name__ == "__main__":
    print("Earth Velocity Calculation at J2000.0 Epoch")
    print("Using astropy ephemeris")
    
    # Get Earth's velocity
    data = get_earth_velocity_j2000()
    ephemeris_used = data.get('ephemeris_used', 'unknown')
    
    # Compare with approximate values
    compare_with_orbital_velocity()
    
    # Output for Julia
    output_for_julia()
    
    print(f"\n" + "="*60)
    print("NOTES:")
    print("- Velocities are barycentric (relative to solar system barycenter)")
    print(f"- Using {ephemeris_used.upper()} ephemeris")
    print("- J2000.0 = 2000-01-01 12:00:00 TT")
    print("- Coordinate system: ICRS (International Celestial Reference System)")
    if ephemeris_used == 'de440':
        print("- DE440: High-precision JPL ephemeris (2020 release)")
    print("="*60)
