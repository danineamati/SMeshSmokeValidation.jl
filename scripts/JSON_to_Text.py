#!/usr/bin/env python3

import json
import glob
import os
import sys
from shapely.geometry import shape, Point
import pyproj

def determine_utm_crs(longitude, latitude):
    """
    Given a longitude/latitude in WGS84, returns a pyproj CRS for the
    corresponding UTM zone.
    """
    zone_number = int((longitude + 180) // 6) + 1
    epsg_code = f"326{zone_number:02d}" if latitude >= 0 else f"327{zone_number:02d}"
    return pyproj.CRS.from_epsg(int(epsg_code))

def process_geojson_files(folder_path, output_text_path):
    """
    Processes multiple GeoJSON files from the given folder and writes all
    transformed UTM coordinates into a single text file.
    """
    input_geojson_paths = glob.glob(os.path.join(folder_path, "*.json")) + glob.glob(os.path.join(folder_path, "*.geojson"))

    if not input_geojson_paths:
        print(f"No GeoJSON files found in {folder_path}")
        sys.exit(1)

    with open(output_text_path, 'w', encoding='utf-8') as txtfile:
        # Write header
        txtfile.write("file_name\tfeature_index\tgeometry_type\tpart_index\tcoordinates\n")

        for input_geojson_path in input_geojson_paths:
            with open(input_geojson_path, 'r', encoding='utf-8') as f:
                geojson_data = json.load(f)

            # Ensure it's a FeatureCollection
            if geojson_data.get('type') != 'FeatureCollection':
                print(f"Skipping {input_geojson_path}: Not a FeatureCollection.")
                continue

            features = geojson_data.get('features', [])
            if not features:
                print(f"Skipping {input_geojson_path}: No features found.")
                continue

            # Determine UTM zone from centroid of first feature
            first_geom_data = features[0]['geometry']
            coords = first_geom_data['coordinates']

            # Ensure first feature has valid coordinates
            if isinstance(coords, list) and len(coords) > 2:
                coords = coords[:2]  # Keep only (lon, lat)

            try:
                first_geom = shape({'type': first_geom_data['type'], 'coordinates': coords})
            except ValueError as e:
                print(f"Skipping {input_geojson_path}: Invalid first geometry - {e}")
                continue

            centroid = first_geom.centroid
            wgs84_lon, wgs84_lat = centroid.x, centroid.y
            utm_crs = determine_utm_crs(wgs84_lon, wgs84_lat)

            # Prepare transformer from WGS84 -> UTM
            wgs84 = pyproj.CRS("EPSG:4326")
            transformer = pyproj.Transformer.from_crs(wgs84, utm_crs, always_xy=True)

            for feature_index, feature in enumerate(features):
                geom_data = feature['geometry']
                coords = geom_data['coordinates']

                # Ensure valid Point coordinates
                if geom_data['type'] == "Point" and isinstance(coords, list) and len(coords) > 2:
                    coords = coords[:2]  # Only take (lon, lat)

                try:
                    geom = shape({'type': geom_data['type'], 'coordinates': coords})
                except ValueError as e:
                    print(f"Skipping feature {feature_index} in {input_geojson_path}: {e}")
                    continue

                geom_type = geom.geom_type

                def transform_coords(coords):
                    return [f"[{transformer.transform(x, y)[0]}, {transformer.transform(x, y)[1]}]" for x, y in coords]

                if geom_type == 'Point':
                    coords_str = f"[{transformer.transform(geom.x, geom.y)[0]}, {transformer.transform(geom.x, geom.y)[1]}]"
                    txtfile.write(f"{os.path.basename(input_geojson_path)}\t{feature_index}\t{geom_type}\t0\t{coords_str}\n")

                elif geom_type in ("LineString", "LinearRing"):
                    coords_str = " ".join(transform_coords(geom.coords))
                    txtfile.write(f"{os.path.basename(input_geojson_path)}\t{feature_index}\t{geom_type}\t0\t{coords_str}\n")

                elif geom_type == "Polygon":
                    coords_str = " ".join(transform_coords(geom.exterior.coords))
                    txtfile.write(f"{os.path.basename(input_geojson_path)}\t{feature_index}\t{geom_type}_ring0\t0\t{coords_str}\n")

                elif geom_type.startswith("Multi"):
                    for part_index, part in enumerate(geom.geoms):
                        coords_str = " ".join(transform_coords(part.coords))
                        txtfile.write(f"{os.path.basename(input_geojson_path)}\t{feature_index}\t{geom_type}\t{part_index}\t{coords_str}\n")

                else:
                    print(f"Skipping unsupported geometry type in {input_geojson_path}: {geom_type}")

    print(f"Transformed coordinates from {len(input_geojson_paths)} files saved to: {output_text_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <folder_path> <output_text_file>")
        sys.exit(1)

    folder_path = sys.argv[1]
    output_text_path = sys.argv[2]

    # Ensure the output directory exists
    output_dir = os.path.dirname(output_text_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    process_geojson_files(folder_path, output_text_path)
