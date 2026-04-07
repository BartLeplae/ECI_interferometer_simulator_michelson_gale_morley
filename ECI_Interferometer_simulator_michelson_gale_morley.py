"""
3D ECI Interferometer Kinematics Simulator
==========================================

This script provides a high-precision, 3D Earth-Centered Inertial (ECI) kinematic 
simulation of classical interferometer experiments, specifically the Michelson-Morley 
and Michelson-Gale experiments. 

Concepts:
- Light travels at the constant speed 'c' relative to the Earth-Centered Inertial 
  (ECI) reference frame. 
- The simulator allows you to add a hypothetical drift vector to test the 
  Michelson-Morley hypothesis.
- The ECI frame is centered on the Earth, aligned with its axis of rotation, and 
  maintains a fixed orientation relative to the distant stars.
- The mirrors of the interferometer travel dynamically with the surface of the 
  rotating Earth while the light is in transit from one mirror to the next.

It models the exact time-of-flight of light traveling between mirrors on the surface 
of a rotating Earth, accounting for:
- The Sagnac effect (Earth's rotation).
- Arbitrary geographical latitudes.
- Custom arm lengths and intersection angles.
- Spherical Earth curvature (using exact 3D spatial chord adjustments).
- Hypothetical "Drift" velocity vectors (e.g., ether wind).

The simulator reads experiment configurations from a `scenarios.yaml` file, 
processes the relativistic/kinematic time differences down to 50 decimal places 
to capture sub-nanometer fringe shifts, and outputs a compiled PDF report 
containing the numerical results, setup schematics, and 360-degree rotation graphs.

Dependencies:
    - pyyaml
    - matplotlib
    - fpdf

Usage:
    Ensure `scenarios.yaml` is in the same directory, then run:
    $ python ECI_Interferometer_simulator_michelson_gale_morley.py
"""

import math
import yaml
import os
import matplotlib.pyplot as plt
from decimal import Decimal, getcontext
from fpdf import FPDF

# ---------------------------------------------------------
# 1. SETUP & HIGH PRECISION CONFIGURATION
# ---------------------------------------------------------
# Interference phase differences require extreme precision to avoid floating-point 
# truncation. 50 decimal places are used to ensure nanosecond-scale accuracy.
getcontext().prec = 50

def dec_sin(x):
    """
    Computes the sine of x using a Taylor series expansion to maintain 
    Decimal precision, bypassing the standard math.sin() float conversion.
    """
    getcontext().prec += 5
    x = Decimal(x)
    t = x
    sum_sin = x
    n = 1
    while True:
        t = t * (-x**2) / Decimal((2*n) * (2*n + 1))
        if abs(t) < Decimal('1e-55'):
            break
        sum_sin += t
        n += 1
    getcontext().prec -= 5
    return sum_sin + Decimal('0')

def dec_cos(x):
    """
    Computes the cosine of x using a Taylor series expansion to maintain 
    Decimal precision, bypassing the standard math.cos() float conversion.
    """
    getcontext().prec += 5
    x = Decimal(x)
    t = Decimal('1')
    sum_cos = t
    n = 1
    while True:
        t = t * (-x**2) / Decimal((2*n - 1) * (2*n))
        if abs(t) < Decimal('1e-55'):
            break
        sum_cos += t
        n += 1
    getcontext().prec -= 5
    return sum_cos + Decimal('0')

class Vec3:
    """Simple 3D Vector class supporting Decimal math operations."""
    def __init__(self, x, y, z):
        self.x = Decimal(x)
        self.y = Decimal(y)
        self.z = Decimal(z)
        
    def __sub__(self, other):
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)
        
    def __add__(self, other):
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)
        
    def __mul__(self, scalar):
        s = Decimal(scalar)
        return Vec3(self.x * s, self.y * s, self.z * s)
        
    def magnitude(self):
        return (self.x**2 + self.y**2 + self.z**2).sqrt()

# ---------------------------------------------------------
# 2. CORE PHYSICS ENGINE
# ---------------------------------------------------------
def run_interferometer_physics(lat_deg, d_AD, d_AF, angle_AD_AF, orientation_deg, wl, earth_rotation, drift_v, drift_angle_deg):
    """
    Executes the 3D kinematic time-of-flight simulation for a given setup.

    Args:
        lat_deg (float): Latitude of the experiment.
        d_AD (Decimal): Length of the primary arm in meters.
        d_AF (Decimal): Length of the secondary arm in meters.
        angle_AD_AF (float): Internal angle between the two arms.
        orientation_deg (float): Azimuthal orientation of the apparatus.
        wl (Decimal): Wavelength of the light source in meters.
        earth_rotation (bool): Enables/disables Earth angular velocity.
        drift_v (float): Magnitude of the external drift vector in m/s.
        drift_angle_deg (float): Direction of the external drift vector.

    Returns:
        dict: Time differences (in seconds) for Michelson-Gale and Michelson-Morley paths.
    """
    c = Decimal('299792458')
    R = Decimal('6371000') # Earth radius in meters
    omega = Decimal('7.292115e-5') if earth_rotation else Decimal('0')
    pi = Decimal(math.pi)

    lat_center_rad = Decimal(lat_deg) * pi / Decimal('180')
    
    orient_rad = Decimal(orientation_deg) * pi / Decimal('180')
    ang_rad = Decimal(angle_AD_AF) * pi / Decimal('180')

    # Construct initial geometry in a local flat 2D tangent plane
    Ax_0, Ay_0 = Decimal('0'), Decimal('0')
    
    Dx_0 = d_AD * dec_sin(orient_rad)
    Dy_0 = d_AD * dec_cos(orient_rad)
    
    Fx_0 = d_AF * dec_sin(orient_rad + ang_rad)
    Fy_0 = d_AF * dec_cos(orient_rad + ang_rad)
    
    Ex_0 = Dx_0 + Fx_0
    Ey_0 = Dy_0 + Fy_0

    # Center the coordinates around the centroid to minimize spherical distortion
    cx = (Ax_0 + Dx_0 + Ex_0 + Fx_0) / Decimal('4')
    cy = (Ay_0 + Dy_0 + Ey_0 + Fy_0) / Decimal('4')

    Ax, Ay = Ax_0 - cx, Ay_0 - cy
    Dx_init, Dy_init = Dx_0 - cx, Dy_0 - cy
    Fx_init, Fy_init = Fx_0 - cx, Fy_0 - cy
    Ex_init, Ey_init = Ex_0 - cx, Ey_0 - cy

    def local_to_ecef(x, y):
        """Projects local 2D tangent plane coordinates onto the 3D spherical Earth (ECEF)."""
        lat_point = lat_center_rad + (y / R)
        lon_point = x / (R * dec_cos(lat_point))
        X = R * dec_cos(lat_point) * dec_cos(lon_point)
        Y = R * dec_cos(lat_point) * dec_sin(lon_point)
        Z = R * dec_sin(lat_point)
        return Vec3(X, Y, Z)

    def micro_adjust(start_local, end_local, start_ecef, target_dist):
        """
        Iterative Newton-Raphson approximation to ensure the 3D ECEF straight-line
        chord distance exactly matches the defined macroscopic parameter distance.
        """
        sx, sy = start_local
        ex, ey = end_local
        dx, dy = ex - sx, ey - sy
        s = Decimal('1')
        for _ in range(30):
            test_x = sx + dx * s
            test_y = sy + dy * s
            test_ecef = local_to_ecef(test_x, test_y)
            current_dist = (test_ecef - start_ecef).magnitude()
            if abs(target_dist - current_dist) < Decimal('1e-45'):
                break
            s = s * (target_dist / current_dist)
        return local_to_ecef(sx + dx * s, sy + dy * s)

    # Lock node coordinates based on exact 3D distances
    node_A_ecef = local_to_ecef(Ax, Ay)
    node_D_ecef = micro_adjust((Ax, Ay), (Dx_init, Dy_init), node_A_ecef, d_AD)
    node_F_ecef = micro_adjust((Ax, Ay), (Fx_init, Fy_init), node_A_ecef, d_AF)
    node_E_ecef = micro_adjust((Dx_init, Dy_init), (Ex_init, Ey_init), node_D_ecef, d_AF)

    nodes_ecef = {'A': node_A_ecef, 'D': node_D_ecef, 'E': node_E_ecef, 'F': node_F_ecef}

    # Map the external drift vector to the Earth-Centered Inertial (ECI) frame
    drift_ang_rad = Decimal(drift_angle_deg) * pi / Decimal('180')
    V_mag = Decimal(drift_v)
    
    V_x_local = V_mag * dec_cos(drift_ang_rad) 
    V_y_local = V_mag * dec_sin(drift_ang_rad) 

    sin_lat_c = dec_sin(lat_center_rad)
    cos_lat_c = dec_cos(lat_center_rad)
    V_drift_eci = Vec3(-V_y_local * sin_lat_c, V_x_local, V_y_local * cos_lat_c)

    def get_eci_pos(node_id, t):
        """Calculates dynamic 3D ECI position of a node at time t."""
        ecef = nodes_ecef[node_id]
        theta = omega * t
        cos_t = dec_cos(theta)
        sin_t = dec_sin(theta)
        x_eci = ecef.x * cos_t - ecef.y * sin_t
        y_eci = ecef.x * sin_t + ecef.y * cos_t
        z_eci = ecef.z
        return Vec3(x_eci, y_eci, z_eci) + (V_drift_eci * t)

    def time_of_flight(start_node, end_node, t_start):
        """Iteratively solves for the exact time a photon intersects the target node."""
        p_start = get_eci_pos(start_node, t_start)
        dt = Decimal('0')
        for _ in range(25):
            t_end = t_start + dt
            p_end = get_eci_pos(end_node, t_end)
            dt_new = (p_end - p_start).magnitude() / c
            if abs(dt_new - dt) < Decimal('1e-45'):
                break
            dt = dt_new
        return dt

    # --- Execute Michelson-Gale Loops ---
    t = Decimal('0')
    t += time_of_flight('A', 'D', t)
    t += time_of_flight('D', 'E', t)
    t += time_of_flight('E', 'F', t)
    t += time_of_flight('F', 'A', t)
    t_CW = t

    t = Decimal('0')
    t += time_of_flight('A', 'F', t)
    t += time_of_flight('F', 'E', t)
    t += time_of_flight('E', 'D', t)
    t += time_of_flight('D', 'A', t)
    t_CCW = t

    # --- Execute Michelson-Morley Arms ---
    t_A_D = time_of_flight('A', 'D', Decimal('0'))
    t_ADA = t_A_D + time_of_flight('D', 'A', t_A_D)

    t_A_F = time_of_flight('A', 'F', Decimal('0'))
    t_AFA = t_A_F + time_of_flight('F', 'A', t_A_F)

    return {
        'MG_time': t_CCW - t_CW,
        'MM_time': t_AFA - t_ADA
    }

# ---------------------------------------------------------
# 3. SCENARIO ORCHESTRATOR & PLOTTING
# ---------------------------------------------------------
def process_scenario(config, scenario_index):
    """
    Parses configuration, invokes the physics engine, formats results, 
    and generates matplotlib plots if requested in the YAML configuration.
    """
    lat_deg = config.get("lat_deg", 41.76666666666667)
    dist_AD = config.get("dist_AD", 339.2424)
    dist_AF = config.get("dist_AF", 612.648)
    angle_AD_AF = config.get("angle_AD_AF", 90.0)
    wavelength = config.get("wavelength", 570e-9)
    earth_rotation = config.get("earth_rotation", True)
    drift_vector = config.get("drift_vector", [0.0, 0.0])
    print_mode = config.get("print_mode", "both").lower()
    
    generate_graph = config.get("generate_graph", False)
    generate_setup_plot = config.get("generate_setup_plot", False)
    setup_plot_type = config.get("setup_plot_type", "mg").lower()
    
    drift_v, drift_angle_deg = drift_vector
    wl = Decimal(wavelength)
    c = Decimal('299792458')

    # Baseline calculation at 0-degree orientation
    res = run_interferometer_physics(
        lat_deg, Decimal(dist_AD), Decimal(dist_AF), angle_AD_AF, 0.0, 
        wl, earth_rotation, drift_v, drift_angle_deg
    )

    status_rot = "ON" if earth_rotation else "OFF"
    output_text = [
        f"Location:       {lat_deg} deg N",
        f"Loop Shape:     A-D: {dist_AD}m, A-F: {dist_AF}m @ Angle: {angle_AD_AF} deg",
        f"Wavelength:     {wavelength} m",
        f"Kinematics:     Earth Rot: {status_rot} | Drift Vector: {drift_v} m/s @ {drift_angle_deg} deg"
    ]

    if print_mode in ["both", "mg"]:
        delta_d_MG = res['MG_time'] * c
        output_text.extend([
            "\n--- Michelson-Gale (0 deg Orientation) ---",
            f"Time Diff:      {res['MG_time'] * Decimal('1e9'):.10f} ns",
            f"Distance Diff:  {delta_d_MG * Decimal('1e9'):.10f} nm",
            f"Fringe Shift:   {delta_d_MG / wl:.10f} fringes"
        ])

    if print_mode in ["both", "mm"]:
        delta_d_MM = res['MM_time'] * c
        output_text.extend([
            "\n--- Michelson-Morley (0 deg Orientation) ---",
            f"Time Diff:      {res['MM_time'] * Decimal('1e9'):.10f} ns",
            f"Distance Diff:  {delta_d_MM * Decimal('1e9'):.10f} nm",
            f"Fringe Shift:   {delta_d_MM / wl:.10f} fringes"
        ])

    setup_filename = None
    graph_filename = None

    # Process Geometric Apparatus Setup Plot
    if generate_setup_plot:
        print(f"Generating geometric setup plot for Scenario {scenario_index}...")
        plt.figure(figsize=(5, 5))
        
        ax, ay = 0.0, 0.0
        dx, dy = 0.0, float(dist_AD)
        ang_rad = math.radians(float(angle_AD_AF))
        fx, fy = float(dist_AF) * math.sin(ang_rad), float(dist_AF) * math.cos(ang_rad)
        ex, ey = dx + fx, dy + fy
        
        node_style = dict(boxstyle="circle,pad=0.3", facecolor="white", edgecolor="black", linewidth=1.5)

        if setup_plot_type == "mm":
            plt.plot([ax, dx], [ay, dy], linestyle='-', color='#d62728', linewidth=2.5, alpha=0.8, zorder=1)
            plt.plot([ax, fx], [ay, fy], linestyle='-', color='#1f77b4', linewidth=2.5, alpha=0.8, zorder=1)
            
            plt.text(ax, ay, 'A', fontsize=12, fontweight='bold', ha='center', va='center', bbox=node_style, zorder=5)
            plt.text(dx, dy, 'D', fontsize=12, fontweight='bold', ha='center', va='center', bbox=node_style, zorder=5)
            plt.text(fx, fy, 'F', fontsize=12, fontweight='bold', ha='center', va='center', bbox=node_style, zorder=5)
            
            plt.title(f'Michelson-Morley Setup (Angle: {angle_AD_AF}°)')

        else:
            xs = [ax, dx, ex, fx, ax]
            ys = [ay, dy, ey, fy, ay]
            plt.plot(xs, ys, linestyle='-', color='#2ca02c', linewidth=2.5, alpha=0.8, zorder=1)
            
            plt.text(ax, ay, 'A', fontsize=12, fontweight='bold', ha='center', va='center', bbox=node_style, zorder=5)
            plt.text(dx, dy, 'D', fontsize=12, fontweight='bold', ha='center', va='center', bbox=node_style, zorder=5)
            plt.text(ex, ey, 'E', fontsize=12, fontweight='bold', ha='center', va='center', bbox=node_style, zorder=5)
            plt.text(fx, fy, 'F', fontsize=12, fontweight='bold', ha='center', va='center', bbox=node_style, zorder=5)
            
            plt.title(f'Michelson-Gale Setup (Angle: {angle_AD_AF}°)')

        plt.xlabel('East (meters)')
        plt.ylabel('North (meters)')
        
        plt.margins(0.15)
        plt.axis('equal') 
        plt.grid(True, linestyle=':', alpha=0.7, zorder=0)

        setup_filename = f"setup_scenario_{scenario_index}.png"
        plt.savefig(setup_filename, bbox_inches='tight', dpi=150)
        plt.close()

    # Process 360-Degree Rotation Kinematic Plot
    if generate_graph:
        print(f"Generating 360-degree rotation plot for Scenario {scenario_index}...")
        angles = range(0, 365, 5)
        mm_fringe_shifts = []

        for angle in angles:
            rot_res = run_interferometer_physics(
                lat_deg, Decimal(dist_AD), Decimal(dist_AF), angle_AD_AF, angle, 
                wl, earth_rotation, drift_v, drift_angle_deg
            )
            fringe_shift = float((rot_res['MM_time'] * c) / wl)
            mm_fringe_shifts.append(fringe_shift)

        max_shift = max(mm_fringe_shifts)
        min_shift = min(mm_fringe_shifts)
        total_shift_diff = max_shift - min_shift

        plt.figure(figsize=(8, 4))
        plt.plot(angles, mm_fringe_shifts, marker='o', linestyle='-', color='#1f77b4', markersize=3)
        plt.title('Michelson-Morley Fringe Shift vs Apparatus Rotation')
        plt.xlabel('Apparatus Orientation (Degrees)')
        plt.ylabel('Fringe Shift')
        
        textstr = f'Max: {max_shift:.4f} | Min: {min_shift:.4f}\nTotal Difference (\u0394N): {total_shift_diff:.4f} fringes'
        props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
        plt.text(180, 0.0, textstr, fontsize=10, verticalalignment='center', horizontalalignment='center', bbox=props, zorder=5)
        plt.axhline(0, color='black', linewidth=0.8, alpha=0.3, zorder=1)

        plt.grid(True, linestyle='--', alpha=0.7, zorder=0)
        plt.xlim(0, 360)
        plt.xticks(range(0, 361, 45))
        
        graph_filename = f"plot_scenario_{scenario_index}.png"
        plt.savefig(graph_filename, bbox_inches='tight', dpi=150)
        plt.close()

    return "\n".join(output_text), graph_filename, setup_filename

# ---------------------------------------------------------
# 4. PDF GENERATOR CLASS
# ---------------------------------------------------------
class PDFReport(FPDF):
    """Custom FPDF class to compile text output and generated plots into a report."""
    def header(self):
        self.set_font('Courier', 'B', 14)
        self.cell(0, 10, '3D ECI Interferometer Simulation Report', 0, 1, 'C')
        self.ln(5)

    def footer(self):
        self.set_y(-15)
        self.set_font('Courier', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

    def add_scenario(self, title, results_text, image_path, setup_path):
        self.set_font('Courier', 'B', 12)
        self.set_fill_color(220, 220, 220)
        self.cell(0, 8, title, 0, 1, 'L', 1)
        self.ln(2)
        
        if setup_path and os.path.exists(setup_path):
            self.image(setup_path, w=110)
            self.ln(5)

        self.set_font('Courier', '', 10)
        self.multi_cell(0, 5, results_text)
        self.ln(5)

        if image_path and os.path.exists(image_path):
            self.image(image_path, w=150)
            self.ln(10)


# =========================================================
# MAIN EXECUTION BLOCK
# =========================================================
if __name__ == "__main__":
    
    # Load configuration from YAML
    try:
        with open("scenarios.yaml", "r") as file:
            data = yaml.safe_load(file)
            scenarios = data.get("scenarios", [])
    except FileNotFoundError:
        print("Error: 'scenarios.yaml' file not found.")
        print("Please create it in the same directory as this script.")
        exit()

    if not scenarios:
        print("No scenarios found in the YAML file.")
        exit()

    # Initialize FPDF report
    pdf = PDFReport()
    pdf.add_page()
    generated_images = []

    # Process dynamically loaded scenarios
    for idx, config in enumerate(scenarios, start=1):
        title = config.get("name", f"SCENARIO {idx}")
        
        print("\n" + "="*50)
        print(f"Processing: {title}")
        print("="*50)
        
        results_text, graph_path, setup_path = process_scenario(config, idx)
        
        print(results_text)
        pdf.add_scenario(title, results_text, graph_path, setup_path)

        # Track temporal image artifacts for cleanup
        if graph_path: generated_images.append(graph_path)
        if setup_path: generated_images.append(setup_path)

    # Save artifact
    pdf_filename = "Interferometer_Report.pdf"
    pdf.output(pdf_filename)
    print(f"\nSuccess! Open '{pdf_filename}' to view your formatted report.")

    # Cleanup temp images
    for img in generated_images:
        try:
            os.remove(img)
        except OSError:
            pass