#!/usr/bin/env python3
"""
Generate a Vesta-shaped mesh, run the C ray tracer, and create a
rotating solar flux GIF for one full revolution.

Vesta parameters (from Dawn mission):
  Semi-axes: 286.3 x 278.6 x 223.2 km
  Obliquity: ~27.5 deg
  Rotation period: 5.342 hours
"""
import os
import sys
import subprocess
import struct
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.animation as animation

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "build")

# ── Vesta parameters ──
A_KM, B_KM, C_KM = 286.3, 278.6, 223.2   # triaxial semi-axes [km]
ZA_DEG = 27.5       # obliquity
TT = 48             # time steps per rotation (48 -> 7.5 deg each)
NMONTE = 80         # MC samples per face

# ── Mesh generation: subdivided icosphere → triaxial ellipsoid ──

def icosphere(subdivisions=3):
    """Create an icosphere by subdividing an icosahedron."""
    # Golden ratio
    t = (1.0 + np.sqrt(5.0)) / 2.0

    # 12 vertices of icosahedron
    verts = np.array([
        [-1,  t,  0], [ 1,  t,  0], [-1, -t,  0], [ 1, -t,  0],
        [ 0, -1,  t], [ 0,  1,  t], [ 0, -1, -t], [ 0,  1, -t],
        [ t,  0, -1], [ t,  0,  1], [-t,  0, -1], [-t,  0,  1],
    ], dtype=float)

    # Normalize to unit sphere
    verts /= np.linalg.norm(verts, axis=1, keepdims=True)

    # 20 triangular faces
    faces = [
        [0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],
        [1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],
        [3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],
        [4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1],
    ]
    faces = np.array(faces, dtype=int)

    # Subdivide
    for _ in range(subdivisions):
        edge_midpoints = {}
        new_faces = []

        def get_midpoint(i0, i1):
            key = (min(i0,i1), max(i0,i1))
            if key in edge_midpoints:
                return edge_midpoints[key]
            mid = (verts[i0] + verts[i1]) / 2.0
            mid /= np.linalg.norm(mid)
            idx = len(verts)
            # Can't append to numpy, collect and extend later
            edge_midpoints[key] = (idx, mid)
            return (idx, mid)

        new_verts = list(verts)
        edge_midpoints_resolved = {}

        def get_mid_idx(i0, i1):
            key = (min(i0,i1), max(i0,i1))
            if key in edge_midpoints_resolved:
                return edge_midpoints_resolved[key]
            mid = (verts[i0] + verts[i1]) / 2.0
            mid /= np.linalg.norm(mid)
            idx = len(new_verts)
            new_verts.append(mid)
            edge_midpoints_resolved[key] = idx
            return idx

        for f in faces:
            a, b, c = f
            ab = get_mid_idx(a, b)
            bc = get_mid_idx(b, c)
            ca = get_mid_idx(c, a)
            new_faces.append([a, ab, ca])
            new_faces.append([b, bc, ab])
            new_faces.append([c, ca, bc])
            new_faces.append([ab, bc, ca])

        verts = np.array(new_verts)
        faces = np.array(new_faces, dtype=int)

    return verts, faces


def make_vesta_mesh(subdivisions=3):
    """Create a Vesta-shaped mesh (triaxial ellipsoid + surface roughness)."""
    verts, faces = icosphere(subdivisions)
    npoint = len(verts)
    nface = len(faces)

    # Scale to Vesta triaxial ellipsoid (normalized to unit scale)
    # We normalize so max axis = 1.0
    scale = max(A_KM, B_KM, C_KM)
    verts[:, 0] *= A_KM / scale
    verts[:, 1] *= B_KM / scale
    verts[:, 2] *= C_KM / scale

    # Add some surface roughness (gentle perturbation with spherical harmonics-like noise)
    np.random.seed(42)
    for i in range(npoint):
        r = np.linalg.norm(verts[i])
        # Low-frequency bumps (Vesta has a big south pole basin ~Rheasilvia)
        lat = np.arcsin(verts[i, 2] / r)
        lon = np.arctan2(verts[i, 1], verts[i, 0])
        bump = 0.02 * np.cos(2*lat) * np.sin(lon)  # degree-2 perturbation
        bump += 0.015 * np.sin(3*lat + 1.0)         # asymmetry
        # South pole depression (Rheasilvia)
        if lat < -1.0:
            bump -= 0.04 * np.cos(lat + 1.0)
        verts[i] *= (1.0 + bump)

    print(f"Vesta mesh: {npoint} vertices, {nface} faces")
    return verts, faces, npoint, nface


def save_mesh(verts, faces, vert_file, face_file):
    """Save mesh in the format expected by the C ray tracer."""
    with open(vert_file, "w") as f:
        for v in verts:
            f.write(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n")

    with open(face_file, "w") as f:
        for fc in faces:
            # C code converts from 1-based to 0-based
            f.write(f"{fc[0]+1} {fc[1]+1} {fc[2]+1}\n")


def read_flux(filename):
    """Read flux output: header nface, tt, then nface*tt values row-major."""
    with open(filename, "r") as f:
        nface = int(f.readline())
        tt = int(f.readline())
        data = np.array([float(f.readline()) for _ in range(nface * tt)])
    return data.reshape(nface, tt), nface, tt


def create_rotation_gif(verts, faces, flux, tt, output_gif):
    """Create a GIF of the asteroid rotating with solar flux coloring."""
    nface = len(faces)

    # Precompute face centroids for sorting
    centroids = np.mean(verts[faces], axis=1)  # (nface, 3)

    # Color setup
    vmin, vmax = 0.0, np.max(flux)
    cmap = plt.cm.hot

    frames = []
    fig = plt.figure(figsize=(7, 7), facecolor='black')
    ax = fig.add_subplot(111, projection='3d', facecolor='black')

    def make_frame(t):
        ax.clear()
        ax.set_facecolor('black')
        fig.patch.set_facecolor('black')

        # Rotation angle for this frame
        angle = t / tt * 360.0

        # Get flux for this time step
        face_flux = flux[:, t]

        # Rotate vertices for viewing
        theta = np.radians(angle)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        R = np.array([[cos_t, -sin_t, 0],
                       [sin_t,  cos_t, 0],
                       [0,      0,     1]])
        rot_verts = verts @ R.T

        # Build triangles
        tri_verts = rot_verts[faces]  # (nface, 3, 3)

        # Depth sort (painter's algorithm)
        rot_centroids = centroids @ R.T
        order = np.argsort(rot_centroids[:, 1])  # sort by depth (y)

        # Build polygon collection
        sorted_tris = tri_verts[order]
        sorted_flux = face_flux[order]

        # Normalize flux to [0, 1] for colormap
        if vmax > 0:
            normed = sorted_flux / vmax
        else:
            normed = np.zeros(len(sorted_flux))
        colors = cmap(normed)
        # Darken unlit faces
        colors[sorted_flux <= 0] = [0.05, 0.05, 0.08, 1.0]

        poly3d = Poly3DCollection(sorted_tris, linewidths=0.1,
                                   edgecolors=(0.15, 0.15, 0.15, 0.3))
        poly3d.set_facecolor(colors)
        ax.add_collection3d(poly3d)

        lim = 1.2
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_axis_off()

        # Fixed camera angle (slightly tilted to show obliquity)
        ax.view_init(elev=25, azim=0)

        # Title
        local_time = t / tt * 24.0
        ax.set_title(f"Vesta Solar Flux  |  Local Time {local_time:.1f}h  |  "
                     f"Obliquity {ZA_DEG}°",
                     color='white', fontsize=11, pad=10)

        # Sun indicator
        fi = (t + 0.5) / tt * 2 * np.pi + np.pi
        za_rad = np.radians(ZA_DEG)
        sun_dir = np.array([np.cos(-za_rad)*np.sin(fi),
                            np.cos(-za_rad)*np.cos(fi),
                            np.sin(-za_rad)])
        sun_pos = sun_dir * 1.5
        sun_rot = sun_pos @ R.T
        ax.scatter([sun_rot[0]], [sun_rot[1]], [sun_rot[2]],
                   color='yellow', s=80, marker='*', zorder=10)

        return []

    print(f"Rendering {tt} frames...")
    ani = animation.FuncAnimation(fig, make_frame, frames=tt,
                                   interval=120, blit=False)
    ani.save(output_gif, writer='pillow', fps=8, dpi=100,
             savefig_kwargs={'facecolor': 'black'})
    plt.close(fig)
    print(f"GIF saved: {output_gif}")


def main():
    # 1. Generate Vesta mesh
    print("=" * 50)
    print("Vesta Ray Tracing + Rotation GIF")
    print("=" * 50)

    verts, faces, npoint, nface = make_vesta_mesh(subdivisions=3)

    vert_file = os.path.join(BUILD_DIR, "vesta_vertices.dat")
    face_file = os.path.join(BUILD_DIR, "vesta_faces.dat")
    flux_file = os.path.join(BUILD_DIR, "vesta_flux.dat")
    output_gif = os.path.join(BUILD_DIR, "vesta_rotation.gif")

    save_mesh(verts, faces, vert_file, face_file)
    print(f"Mesh saved: {vert_file}, {face_file}")

    # 2. Run C ray tracer
    dsmc_exe = os.path.join(BUILD_DIR, "dsmc_suite")
    if not os.path.exists(dsmc_exe):
        print(f"ERROR: {dsmc_exe} not found. Build first.")
        sys.exit(1)

    cmd = [
        dsmc_exe, "raytrace",
        "-m", vert_file,
        "-f", face_file,
        "--npoint", str(npoint),
        "--nface", str(nface),
        "--za", str(ZA_DEG),
        "--tt", str(TT),
        "--nmonte", str(NMONTE),
        "-o", os.path.join(BUILD_DIR, "vesta"),
        "-t", "8",
    ]
    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("STDERR:", result.stderr)
        sys.exit(1)

    # 3. Read flux results
    flux_path = os.path.join(BUILD_DIR, "vesta_flux.dat")
    flux, nf_read, tt_read = read_flux(flux_path)
    print(f"Flux loaded: {nf_read} faces x {tt_read} time steps")
    print(f"  Max flux: {flux.max():.4f}")
    print(f"  Mean illuminated flux: {flux[flux>0].mean():.4f}")

    # 4. Create rotation GIF
    create_rotation_gif(verts, faces, flux, TT, output_gif)

    print("\nDone!")


if __name__ == "__main__":
    main()
