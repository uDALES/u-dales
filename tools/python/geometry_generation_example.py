"""
Geometry Generation Examples for uDALES

This script demonstrates how to use the geometry generation module to create
synthetic geometries for urban canopy simulations.

Copyright (C) 2024 the uDALES Team.
"""

import numpy as np
import sys
from pathlib import Path

# Add the tools/python directory to the path
sys.path.insert(0, str(Path(__file__).parent))

try:
    from geometry_generation import (
        create_flat_surface,
        create_cubes,
        create_canyons,
        create_realistic
    )
    print("✓ Geometry generation module imported successfully\n")
except ImportError as e:
    print(f"✗ Error importing geometry generation module: {e}")
    sys.exit(1)


def example_1_flat_surface():
    """Example 1: Create a simple flat surface"""
    print("=" * 60)
    print("Example 1: Creating Flat Surface")
    print("=" * 60)
    
    # Create 100m x 100m flat surface with 10m edge length
    geom = create_flat_surface(xsize=100, ysize=100, edgelength=10)
    
    print(f"Created flat surface:")
    print(f"  - Number of faces: {geom.n_faces}")
    print(f"  - Number of vertices: {geom.n_vertices}")
    print(f"  - Total area: {geom.total_area:.2f} m²")
    print(f"  - Bounds: {geom.bounds}")
    
    # Save to file
    output_path = Path(__file__).parent / "flat_surface.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}\n")
    
    return geom


def example_2_single_cube():
    """Example 2: Create a single cube"""
    print("=" * 60)
    print("Example 2: Creating Single Cube")
    print("=" * 60)
    
    # Create 10x10x15m cube centered in 100x100m domain
    geom = create_cubes(
        xsize=100, ysize=100,
        Hx=10, Hy=10, Hz=15,
        geom_option='S',
        add_ground=True
    )
    
    print(f"Created single cube with ground:")
    print(f"  - Number of faces: {geom.n_faces}")
    print(f"  - Cube dimensions: 10m x 10m x 15m")
    print(f"  - Domain: 100m x 100m")
    print(f"  - Is watertight: {geom.is_watertight}")
    
    # Save to file
    output_path = Path(__file__).parent / "single_cube.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}\n")
    
    return geom


def example_3_aligned_cube_array():
    """Example 3: Create aligned cube array"""
    print("=" * 60)
    print("Example 3: Creating Aligned Cube Array")
    print("=" * 60)
    
    # Create array of 10x10x15m cubes with 10m spacing
    geom = create_cubes(
        xsize=100, ysize=100,
        Hx=10, Hy=10, Hz=15,
        Cx=10, Cy=10,
        geom_option='AC',
        add_ground=True
    )
    
    print(f"Created aligned cube array:")
    print(f"  - Cube dimensions: 10m x 10m x 15m")
    print(f"  - Spacing: 10m x 10m")
    print(f"  - Number of faces: {geom.n_faces}")
    print(f"  - Total surface area: {geom.total_area:.2f} m²")
    
    # Calculate volume
    print(f"  - Total volume: {geom.volume:.2f} m³")
    
    # Save to file
    output_path = Path(__file__).parent / "aligned_cubes.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}\n")
    
    return geom


def example_4_staggered_cube_array():
    """Example 4: Create staggered cube array"""
    print("=" * 60)
    print("Example 4: Creating Staggered Cube Array")
    print("=" * 60)
    
    # Create staggered array of cubes
    geom = create_cubes(
        xsize=100, ysize=100,
        Hx=10, Hy=10, Hz=15,
        Cx=10, Cy=10,
        geom_option='SC',
        add_ground=True
    )
    
    print(f"Created staggered cube array:")
    print(f"  - Cube dimensions: 10m x 10m x 15m")
    print(f"  - Spacing: 10m x 10m")
    print(f"  - Pattern: Staggered (offset every other row)")
    print(f"  - Number of faces: {geom.n_faces}")
    
    # Save to file
    output_path = Path(__file__).parent / "staggered_cubes.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}\n")
    
    return geom


def example_5_street_canyon_x():
    """Example 5: Create street canyons along x-axis"""
    print("=" * 60)
    print("Example 5: Creating Street Canyons (x-orientation)")
    print("=" * 60)
    
    # Create canyons aligned along x-axis
    geom = create_canyons(
        xsize=200, ysize=100,
        Hx=200, Hy=10, Hz=20,
        Cx=0, Cy=10,
        orientation='x',
        add_ground=True
    )
    
    print(f"Created street canyons along x-axis:")
    print(f"  - Building dimensions: 200m x 10m x 20m")
    print(f"  - Street width: 10m")
    print(f"  - Domain: 200m x 100m")
    print(f"  - Number of faces: {geom.n_faces}")
    
    # Save to file
    output_path = Path(__file__).parent / "canyon_x.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}\n")
    
    return geom


def example_6_street_canyon_y():
    """Example 6: Create street canyons along y-axis"""
    print("=" * 60)
    print("Example 6: Creating Street Canyons (y-orientation)")
    print("=" * 60)
    
    # Create canyons aligned along y-axis
    geom = create_canyons(
        xsize=100, ysize=200,
        Hx=10, Hy=200, Hz=20,
        Cx=10, Cy=0,
        orientation='y',
        add_ground=True
    )
    
    print(f"Created street canyons along y-axis:")
    print(f"  - Building dimensions: 10m x 200m x 20m")
    print(f"  - Street width: 10m")
    print(f"  - Domain: 100m x 200m")
    print(f"  - Number of faces: {geom.n_faces}")
    
    # Save to file
    output_path = Path(__file__).parent / "canyon_y.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}\n")
    
    return geom


def example_7_realistic_simple():
    """Example 7: Create simple realistic building layout"""
    print("=" * 60)
    print("Example 7: Creating Simple Realistic Layout")
    print("=" * 60)
    
    # Define buildings with different sizes
    buildings = [
        {'position': (25, 25), 'size': (20, 20, 30)},
        {'position': (75, 25), 'size': (15, 25, 25)},
        {'position': (25, 75), 'size': (18, 18, 35)},
        {'position': (75, 75), 'size': (22, 16, 28)},
    ]
    
    geom = create_realistic(
        xsize=100, ysize=100,
        building_configs=buildings,
        add_ground=True
    )
    
    print(f"Created realistic building layout:")
    print(f"  - Number of buildings: {len(buildings)}")
    print(f"  - Domain: 100m x 100m")
    print(f"  - Number of faces: {geom.n_faces}")
    print(f"  - Total volume: {geom.volume:.2f} m³")
    
    # Print individual building info
    print("\nBuilding details:")
    for i, b in enumerate(buildings, 1):
        pos = b['position']
        size = b['size']
        print(f"  Building {i}: {size[0]}x{size[1]}x{size[2]}m at ({pos[0]}, {pos[1]})")
    
    # Save to file
    output_path = Path(__file__).parent / "realistic_simple.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}\n")
    
    return geom


def example_8_realistic_rotated():
    """Example 8: Create realistic layout with rotated buildings"""
    print("=" * 60)
    print("Example 8: Creating Realistic Layout with Rotations")
    print("=" * 60)
    
    # Define buildings with rotations
    buildings = [
        {'position': (30, 30), 'size': (20, 15, 25), 'rotation': 0},
        {'position': (70, 30), 'size': (18, 18, 30), 'rotation': 15},
        {'position': (30, 70), 'size': (22, 12, 28), 'rotation': -20},
        {'position': (70, 70), 'size': (16, 20, 32), 'rotation': 30},
    ]
    
    geom = create_realistic(
        xsize=100, ysize=100,
        building_configs=buildings,
        add_ground=True
    )
    
    print(f"Created realistic layout with rotations:")
    print(f"  - Number of buildings: {len(buildings)}")
    print(f"  - Domain: 100m x 100m")
    print(f"  - Number of faces: {geom.n_faces}")
    
    # Print individual building info
    print("\nBuilding details:")
    for i, b in enumerate(buildings, 1):
        pos = b['position']
        size = b['size']
        rot = b.get('rotation', 0)
        print(f"  Building {i}: {size[0]}x{size[1]}x{size[2]}m at ({pos[0]}, {pos[1]}), "
              f"rotation={rot}°")
    
    # Save to file
    output_path = Path(__file__).parent / "realistic_rotated.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}\n")
    
    return geom


def example_9_comparison():
    """Example 9: Compare different geometry types"""
    print("=" * 60)
    print("Example 9: Comparing Geometry Properties")
    print("=" * 60)
    
    # Create different geometries
    print("Creating geometries...")
    
    flat = create_flat_surface(100, 100, edgelength=10)
    single = create_cubes(100, 100, 10, 10, 15, geom_option='S', add_ground=False)
    aligned = create_cubes(100, 100, 10, 10, 15, Cx=10, Cy=10, 
                          geom_option='AC', add_ground=False)
    canyon = create_canyons(100, 100, 10, 100, 20, Cx=10, Cy=0,
                           orientation='y', add_ground=False)
    
    # Print comparison table
    print("\nGeometry Comparison:")
    print("-" * 80)
    print(f"{'Type':<20} {'Faces':>10} {'Vertices':>10} {'Volume (m³)':>15} {'Area (m²)':>15}")
    print("-" * 80)
    
    geometries = [
        ('Flat Surface', flat),
        ('Single Cube', single),
        ('Aligned Array', aligned),
        ('Street Canyon', canyon)
    ]
    
    for name, geom in geometries:
        print(f"{name:<20} {geom.n_faces:>10} {geom.n_vertices:>10} "
              f"{geom.volume:>15.2f} {geom.total_area:>15.2f}")
    
    print("-" * 80)
    print()


def example_10_visualization():
    """Example 10: Visualize generated geometries"""
    print("=" * 60)
    print("Example 10: Visualizing Generated Geometries")
    print("=" * 60)
    
    print("This example demonstrates visualization of generated geometries.")
    print("Uncomment the show() calls to display 3D visualizations.\n")
    
    # Create a simple geometry
    geom = create_cubes(
        xsize=50, ysize=50,
        Hx=10, Hy=10, Hz=15,
        Cx=5, Cy=5,
        geom_option='AC',
        add_ground=True
    )
    
    print(f"Created aligned cube array:")
    print(f"  - Domain: 50m x 50m")
    print(f"  - Cube size: 10m x 10m x 15m")
    print(f"  - Spacing: 5m")
    print(f"  - Total faces: {geom.n_faces}")
    
    # Save
    output_path = Path(__file__).parent / "visualization_demo.stl"
    geom.save(str(output_path))
    print(f"✓ Saved to {output_path}")
    
    # Visualization (uncomment to display)
    # print("\nDisplaying geometry...")
    # geom.show()
    
    print("\nTo visualize: ")
    print("  geom.show()  # Opens interactive 3D viewer")
    print()


def main():
    """Run all examples"""
    print("\n" + "=" * 60)
    print("UDALES Geometry Generation Examples")
    print("=" * 60 + "\n")
    
    examples = [
        example_1_flat_surface,
        example_2_single_cube,
        example_3_aligned_cube_array,
        example_4_staggered_cube_array,
        example_5_street_canyon_x,
        example_6_street_canyon_y,
        example_7_realistic_simple,
        example_8_realistic_rotated,
        example_9_comparison,
        example_10_visualization,
    ]
    
    for i, example in enumerate(examples, 1):
        try:
            example()
        except Exception as e:
            print(f"✗ Example {i} failed: {e}\n")
            import traceback
            traceback.print_exc()
            continue
    
    print("=" * 60)
    print("All examples completed!")
    print("=" * 60)
    print("\nGenerated STL files can be found in:")
    print(f"  {Path(__file__).parent}")
    print("\nThese geometries can be used as input for uDALES simulations.")


if __name__ == "__main__":
    main()
