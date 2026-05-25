# Approximate Link Inertias — KUKA LBR Med 7 R800

## Method

Each of the 7 links is approximated as a **single solid cylinder** sized to the
arm's outer envelope. Lengths come from the published iiwa kinematic offsets
(0.34 / 0.40 / 0.40 / 0.126 m for the limb segments; the shoulder, elbow and
wrist joint-pairs are co-located, so their housings are modelled as short
cylinders). Radii are read off the photograph, tapering from ~65 mm at the base
to ~35 mm at the flange.

A real cobot link is a hollow aluminium shell containing a motor, gearbox and
torque sensor. Instead of modelling those internals, the solid cylinder is given
a single **effective density** ρ = 1565 kg/m³ (~58% of solid aluminium). That
value is calibrated so the seven cylinder masses sum to the catalogue robot mass
of ~24 kg. Using the larger cylinders for the proximal links automatically puts
more mass near the base, as in the real arm.

Link frame *i* sits at joint *i* with its **z-axis along the joint axis**. Since
each cylinder is axisymmetric about that axis, every COM lies on z and each
inertia tensor is **diagonal** in the link frame (products of inertia ≈ 0).

## Formulas

For a cylinder of mass m = ρ·πr²L, length L along its axis, radius r:

- axial inertia: **Izz = ½ m r²**
- transverse inertia: **Ixx = Iyy = (1/12) m (3r² + L²)**

both taken about the COM, located at the cylinder's mid-length.

## Per-link results

| Link | Calculation (L, r → m; Izz; Ixx=Iyy) |
|---|---|
| 1 base/shoulder | L=0.34, r=0.065 → m=7.06 kg; Izz=0.0149; Ixx=Iyy=0.0755; com_z=0.150 |
| 2 shoulder pitch | L=0.18, r=0.060 → m=3.19 kg; Izz=0.00573; Ixx=Iyy=0.0115; com_z=0 |
| 3 upper arm | L=0.40, r=0.055 → m=5.95 kg; Izz=0.00900; Ixx=Iyy=0.0838; com_z=0.180 |
| 4 elbow pitch | L=0.16, r=0.050 → m=1.97 kg; Izz=0.00246; Ixx=Iyy=0.00542; com_z=0 |
| 5 forearm | L=0.40, r=0.045 → m=3.98 kg; Izz=0.00403; Ixx=Iyy=0.0551; com_z=0.170 |
| 6 wrist pitch | L=0.14, r=0.040 → m=1.10 kg; Izz=0.000881; Ixx=Iyy=0.00224; com_z=0 |
| 7 flange/tool | L=0.126, r=0.035 → m=0.76 kg; Izz=0.000465; Ixx=Iyy=0.00124; com_z=0.050 |

Units: kg, m, kg·m². Total mass = 24.0 kg.