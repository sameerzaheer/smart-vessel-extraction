# SmartVesselExtraction

A C++ research tool, written in 2012–2013 as part of a Master's thesis, that
extracts **blood-vessel centrelines from 3D medical images**.

Starting from an image volume and a set of initial vessel lines/seeds, it
tracks each vessel through the volume:

1. **Hessian analysis** at each point gives the local vessel direction
   (eigenvectors of the Hessian of the image).
2. A **cross-section** perpendicular to that direction is sampled, and the
   point is re-centred on the vessel's intensity maximum.
3. **Rays are cast** across the cross-section (8 directions) to estimate the
   vessel **radius**.
4. Tracking steps forward, smoothing the path, until termination checks
   (direction reversal, intensity/eigenvalue tests) stop it.

The output is a VTK polyline file of the extracted centrelines (plus
`radiusData.txt` with per-point radii).

## Usage

```
SmartVesselExtraction <image.mhd> <lines.vtk>
```

`Readme.txt` holds the original working notes with the commands used in each
experiment.

## Code layout

| File | What it does |
|---|---|
| `itkVesselExtraction.h/.hxx` | Main extraction class: `Extraction_MAIN`, tracking loop |
| `itkVesselExtractionTools.h/.hxx` | Hessian, cross-section, recentring, radius and I/O helpers |
| `itkSandBox.h/.hxx` | Experimental copy of the extraction class (2013) |
| `VesselInfoNew.h` | Per-vessel data (seed, direction, length, …) |
| `MainSmartVesselExtraction.cxx`, `TestVesselExtraction.cxx` | Program entry point / driver |
| `CmdVesselExtraction.cxx` | Alternative command-line driver (not in the CMake target) |
| `CodeRefs/` | Earlier prototypes (`VesselTree9`, `VesselInfo`) kept for reference |

## Dependencies (as originally built)

| Dependency | Version |
|---|---|
| [ITK](https://itk.org) | 3.20.1 |
| [VTK](https://vtk.org) | 5.8.0 (needs `VTK_USE_RENDERING`) |
| CMake | 2.8.7 |
| Compiler / IDE | Visual Studio 2008 (2012 snapshots), Visual Studio 2010 (2013 snapshot) |
| OS | Windows (32-bit) |

`CMakeLists.txt` hard-codes `ITK_DIR` / `VTK_DIR` to the original Windows
paths (`C:/VITK_10/...`) and links VTK 5.x `.lib` files by name, including
`MapReduceMPI.lib` and `mpistubs.lib`. Expect to adjust it before building
with a modern toolchain; the code has not been built since 2013.

## Data

No image data is included. The original experiments used clinical and
animal (pig lung) image volumes that cannot be shared. File names such as
`EHPS3R3Crop.mhd` in `Readme.txt` refer to those datasets.

## History

History reconstructed in 2026 from dated snapshot archives (2012–2013). Each
backdated commit corresponds to one archive and is dated by the
modification times of the files inside it. The project had no version
control at the time, so changes between snapshots appear as single commits.
Build output, IDE files and data were excluded.

Some source files carry copyright headers from the lab where the work was
done (CIGITI, The Hospital for Sick Children) or ITK template headers. They
are left as they were.
