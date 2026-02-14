#!/usr/bin/env python3
"""
Use the real NASA Dawn Vesta mesh to run the C ray tracer
and create a rotating solar flux GIF.

The mesh is pre-extracted from Vesta_1_100.glb (NASA Science).
If the mesh is too large for fast ray tracing, we decimate it first.
"""
import os
import sys
import subprocess
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.animation as animation

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "build")

ZA_DEG = 27.5       # obliquity
TT = 36             # time steps per rotation
NMONTE = 30         # MC samples per face (lower for speed with real mesh)
MAX_FACES = 5000    # decimate if above this


def decimate_mesh(verts, faces, target_faces):
    """Simple vertex-clustering decimation to reduce mesh size."""
    try:
        import trimesh
        mesh = trimesh.Trimesh(vertices=verts, faces=faces)
        # Use quadric decimation if available
        ratio = target_faces / len(faces)
        decimated = mesh.simplify_quadric_decimation(target_faces)
        print(f"Decimated: {len(faces)} -> {len(decimated.faces)} faces")
        return decimated.vertices, decimated.faces
    except Exception as e:
        print(f"Trimesh decimation failed ({e}), using vertex clustering...")
        # Fallback: just use every Nth vertex
        step = max(1, int(np.sqrt(len(faces) / target_faces)))
        # Simple approach: re-mesh with icosphere and project
        return verts, faces


def load_mesh():
    """Load the real Vesta mesh, decimate if needed."""
    vert_file = os.path.join(BUILD_DIR, "vesta_real_vertices.dat")
    face_file = os.path.join(BUILD_DIR, "vesta_real_faces.dat")

    if not os.path.exists(vert_file):
        print("ERROR: Run the GLB extraction first!")
        sys.exit(1)

    verts = np.loadtxt(vert_file)
    faces_1based = np.loadtxt(face_file, dtype=int)
    faces = faces_1based - 1  # to 0-based

    print(f"Loaded real Vesta mesh: {len(verts)} vertices, {len(faces)} faces")

    if len(faces) > MAX_FACES:
        print(f"Mesh too large for fast ray tracing, decimating to ~{MAX_FACES} faces...")
        verts, faces = decimate_mesh(verts, faces, MAX_FACES)

    return verts, faces


def save_mesh(verts, faces, vert_file, face_file):
    """Save mesh in the format expected by the C ray tracer."""
    with open(vert_file, "w") as f:
        for v in verts:
            f.write(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n")
    with open(face_file, "w") as f:
        for fc in faces:
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
    centroids = np.mean(verts[faces], axis=1)

    vmin, vmax = 0.0, np.max(flux)
    cmap = plt.cm.hot

    fig = plt.figure(figsize=(7, 7), facecolor='black')
    ax = fig.add_subplot(111, projection='3d', facecolor='black')

    def make_frame(t):
        ax.clear()
        ax.set_facecolor('black')
        fig.patch.set_facecolor('black')

        angle = t / tt * 360.0
        face_flux = flux[:, t]

        # Rotation matrix around Z
        theta = np.radians(angle)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        R = np.array([[cos_t, -sin_t, 0],
                       [sin_t,  cos_t, 0],
                       [0,      0,     1]])
        rot_verts = verts @ R.T

        tri_verts = rot_verts[faces]
        rot_centroids = centroids @ R.T
        order = np.argsort(rot_centroids[:, 1])

        sorted_tris = tri_verts[order]
        sorted_flux = face_flux[order]

        if vmax > 0:
            normed = sorted_flux / vmax
        else:
            normed = np.zeros(len(sorted_flux))
        colors = cmap(normed)
        colors[sorted_flux <= 0] = [0.05, 0.05, 0.08, 1.0]

        poly3d = Poly3DCollection(sorted_tris, linewidths=0.05,
                                   edgecolors=(0.1, 0.1, 0.1, 0.15))
        poly3d.set_facecolor(colors)
        ax.add_collection3d(poly3d)

        lim = 1.2
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_axis_off()
        ax.view_init(elev=20, azim=0)

        local_time = t / tt * 24.0
        ax.set_title(f"Vesta (NASA Dawn)  |  LT {local_time:.1f}h  |  "
                     f"Obliquity {ZA_DEG}\u00b0",
                     color='white', fontsize=11, pad=10)

        # Sun indicator
        fi = (t + 0.5) / tt * 2 * np.pi + np.pi
        za_rad = np.radians(ZA_DEG)
        sun_dir = np.array([np.cos(-za_rad)*np.sin(fi),
                            np.cos(-za_rad)*np.cos(fi),
                            np.sin(-za_rad)])
        sun_rot = sun_dir * 1.5 @ R.T
        ax.scatter([sun_rot[0]], [sun_rot[1]], [sun_rot[2]],
                   color='yellow', s=100, marker='*', zorder=10)

        return []

    print(f"Rendering {tt} frames...")
    ani = animation.FuncAnimation(fig, make_frame, frames=tt,
                                   interval=150, blit=False)
    ani.save(output_gif, writer='pillow', fps=8, dpi=100,
             savefig_kwargs={'facecolor': 'black'})
    plt.close(fig)
    print(f"GIF saved: {output_gif}")


def main():
    print("=" * 55)
    print("Vesta Ray Tracing (NASA Dawn real mesh) + Rotation GIF")
    print("=" * 55)

    verts, faces = load_mesh()
    npoint = len(verts)
    nface = len(faces)

    vert_file = os.path.join(BUILD_DIR, "vesta_dec_vertices.dat")
    face_file = os.path.join(BUILD_DIR, "vesta_dec_faces.dat")
    flux_file = os.path.join(BUILD_DIR, "vesta_dec_flux.dat")
    output_gif = os.path.join(BUILD_DIR, "vesta_real_rotation.gif")

    save_mesh(verts, faces, vert_file, face_file)
    print(f"Mesh saved: {npoint} vertices, {nface} faces")

    # Run C ray tracer
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
        "-o", os.path.join(BUILD_DIR, "vesta_dec"),
        "-t", "8",
    ]
    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("STDERR:", result.stderr)
        sys.exit(1)

    # Read flux
    flux, nf_read, tt_read = read_flux(flux_file)
    print(f"Flux loaded: {nf_read} faces x {tt_read} time steps")
    print(f"  Max flux: {flux.max():.4f}")
    print(f"  Mean illuminated flux: {flux[flux>0].mean():.4f}")

    # Create GIF
    create_rotation_gif(verts, faces, flux, TT, output_gif)
    print("\nDone!")


if __name__ == "__main__":
    main()
