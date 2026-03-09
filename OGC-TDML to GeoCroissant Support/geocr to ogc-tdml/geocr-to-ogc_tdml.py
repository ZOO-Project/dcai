import argparse
import json
import os
import re
from pathlib import Path
import pytdml
import pytdml.io
from pytdml.type.extended_types import EOTrainingDataset, AI_EOTrainingData, AI_PixelLabel, MD_Band, AI_EOTask
from pytdml.type.basic_types import NamedValue, MD_Identifier

def convert_geocroissant_to_tdml_objectmodel(geocroissant_path, tdml_output_path):
    # Load GeoCroissant JSON directly
    with open(geocroissant_path) as f:
        croissant = json.load(f)

    # Extract standard schema.org properties defined in croissant.ttl
    identifier = croissant.get("identifier", "")
    name = croissant.get("name", "")
    description = croissant.get("description", "No description provided.")
    license_ = croissant.get("license", "")
    if isinstance(license_, list):
        license_ = license_[0] if license_ else ""
    
    # Handle creator/providers (schema:creator)
    creator = croissant.get("creator", {})
    if isinstance(creator, dict):
        providers = [creator.get("name", "Unknown")]
    elif isinstance(creator, list):
        providers = [c.get("name", "Unknown") if isinstance(c, dict) else str(c) for c in creator]
    else:
        providers = []
    
    # Handle dates (schema:dateCreated, schema:dateModified, schema:datePublished)
    created_time = croissant.get("dateCreated", "") or croissant.get("datePublished", "2025-07-17")
    updated_time = croissant.get("dateModified", "") or created_time
    
    # Handle version (schema:version)
    version = croissant.get("version", "")
    
    # Handle spatialCoverage (schema:spatialCoverage with schema:GeoShape)
    spatial_coverage = croissant.get("spatialCoverage")
    extent = None
    if spatial_coverage and isinstance(spatial_coverage, dict):
        geo = spatial_coverage.get("geo")
        if geo and isinstance(geo, dict):
            # Handle box format: "minLat minLon maxLat maxLon"
            box = geo.get("box")
            if box:
                coords = [float(x) for x in box.split()]
                if len(coords) == 4:
                    # box format is: minLat minLon maxLat maxLon
                    # extent format is: [minLon, minLat, maxLon, maxLat]
                    extent = [coords[1], coords[0], coords[3], coords[2]]
            # Handle coordinates format (polygon)
            elif "coordinates" in geo:
                coords = geo.get("coordinates", [[]])[0]
                if coords and len(coords) > 0:
                    min_x = min(c[0] for c in coords)
                    min_y = min(c[1] for c in coords)
                    max_x = max(c[0] for c in coords)
                    max_y = max(c[1] for c in coords)
                    extent = [min_x, min_y, max_x, max_y]
    
    # Extract bands from geocr:spectralBandMetadata
    bands = []
    spectral_bands = croissant.get("geocr:spectralBandMetadata", [])
    for band in spectral_bands:
        if not isinstance(band, dict):
            continue
        band_name = band.get("name", "")
        # Extract center wavelength
        center_wl = band.get("geocr:centerWavelength", {})
        if isinstance(center_wl, dict):
            wl_value = center_wl.get("value", "")
            wl_unit = center_wl.get("unitText", "")
            unit_str = f"{wl_value}{wl_unit}" if wl_value and wl_unit else ""
            if unit_str:
                bands.append(MD_Band(
                    name=[MD_Identifier(code=band_name)],
                    description=band_name,
                    units=unit_str
                ))
            else:
                bands.append(MD_Band(
                    name=[MD_Identifier(code=band_name)],
                    description=band_name
                ))
        else:
            bands.append(MD_Band(
                name=[MD_Identifier(code=band_name)],
                description=band_name
            ))
    
    # Extract classes from keywords or create defaults for segmentation
    classes = []
    keywords = croissant.get("keywords", [])
    
    # Look for class information in the dataset name or keywords
    dataset_name_lower = name.lower()
    if "burn" in dataset_name_lower or any("burn" in str(k).lower() for k in keywords):
        # Burn scar segmentation classes
        classes = [
            NamedValue(key="background", name="background", value=0),
            NamedValue(key="burn_scar", name="burn_scar", value=1)
        ]
    else:
        # Generic segmentation classes
        classes = [
            NamedValue(key="background", name="background", value=0),
            NamedValue(key="foreground", name="foreground", value=1)
        ]
    
    # Parse recordSet to find image/mask field pairs
    record_sets = croissant.get("recordSet", [])
    image_field = None
    mask_field = None
    
    for rs in record_sets:
        if not isinstance(rs, dict):
            continue
        fields = rs.get("field", [])
        for field in fields:
            if not isinstance(field, dict):
                continue
            field_name = field.get("name", "").lower()
            field_id = field.get("@id", "").lower()
            
            # Identify image and mask fields
            if "image" in field_name or "image" in field_id:
                image_field = field
            elif "mask" in field_name or "mask" in field_id or "label" in field_name:
                mask_field = field
    
    # Build data entries from distribution
    # For proper GeoCroissant, distribution contains FileObject/FileSet definitions
    # The actual data references are in the recordSet fields via source
    distribution = croissant.get("distribution", [])
    
    # Parse FileSet to get file patterns (if needed for TDML)
    # For now, we'll create placeholder data entries since the actual file enumeration
    # would require accessing the file system
    data = []
    
    # Check if there are actual file URLs in distribution (legacy format)
    has_legacy_distribution = False
    for dist in distribution:
        if isinstance(dist, dict) and dist.get("@type") not in ["cr:FileObject", "cr:FileSet"]:
            has_legacy_distribution = True
            break
    
    if has_legacy_distribution:
        # Legacy format: distribution contains direct file URLs
        for i in range(0, len(distribution), 2):
            img_entry = distribution[i]
            mask_entry = distribution[i+1] if i+1 < len(distribution) else None
            
            if not isinstance(img_entry, dict):
                continue
                
            img_url = img_entry.get("contentUrl", "")
            data_id = f"data_{i//2}"
            labels = []
            
            if mask_entry and isinstance(mask_entry, dict):
                mask_url = mask_entry.get("contentUrl", "")
                mask_format = mask_entry.get("encodingFormat", "image/tiff")
                labels = [AI_PixelLabel(
                    type="AI_PixelLabel",
                    image_url=[mask_url],
                    image_format=[mask_format],
                    class_=""
                )]
            
            data.append(AI_EOTrainingData(
                type="AI_EOTrainingData",
                id=data_id,
                data_url=[img_url],
                labels=labels
            ))
    else:
        # Proper Croissant format: Files are referenced via FileSet
        # Enumerate files from the file system based on patterns
        print("Note: GeoCroissant uses FileSet structure. Enumerating files from filesystem...")
        
        # Get base path from FileObject distribution
        base_path = ""
        for dist in distribution:
            if isinstance(dist, dict) and dist.get("@type") == "cr:FileObject":
                base_path = dist.get("contentUrl", "")
                break
        
        # Extract regex patterns from image and mask fields
        image_pattern = None
        mask_pattern = None
        
        if image_field:
            source = image_field.get("source", {})
            transform = source.get("transform", {})
            image_pattern = transform.get("regex", ".*_merged\\.tif$")
        
        if mask_field:
            source = mask_field.get("source", {})
            transform = source.get("transform", {})
            mask_pattern = transform.get("regex", ".*\\.mask\\.tif$")
        
        # Build absolute base path
        if base_path.startswith("./") or base_path.startswith("../"):
            # Relative path - resolve relative to the croissant file location
            croissant_dir = os.path.dirname(os.path.abspath(geocroissant_path))
            base_path = os.path.normpath(os.path.join(croissant_dir, base_path))
        
        # Enumerate files matching the patterns
        image_files = []
        mask_files = []
        
        if os.path.exists(base_path):
            # Walk through the directory
            for root, dirs, files in os.walk(base_path):
                for filename in files:
                    file_path = os.path.join(root, filename)
                    rel_path = os.path.relpath(file_path, base_path)
                    
                    # Check if file matches image pattern
                    if image_pattern and re.search(image_pattern, filename):
                        image_files.append(rel_path)
                    
                    # Check if file matches mask pattern
                    if mask_pattern and re.search(mask_pattern, filename):
                        mask_files.append(rel_path)
            
            # Sort to ensure consistent ordering
            image_files.sort()
            mask_files.sort()
            
            print(f"Found {len(image_files)} image files and {len(mask_files)} mask files")
            
            # Create data entries by pairing images with masks
            # Assuming paired naming: image ends with _merged.tif, mask ends with .mask.tif
            # Extract base names and match them
            image_base_map = {}
            for img in image_files:
                # Extract base name (remove _merged.tif suffix)
                base_name = re.sub(r'_merged\.tif$', '', os.path.basename(img))
                image_base_map[base_name] = img
            
            mask_base_map = {}
            for msk in mask_files:
                # Extract base name (remove .mask.tif suffix)
                base_name = re.sub(r'\.mask\.tif$', '', os.path.basename(msk))
                mask_base_map[base_name] = msk
            
            # Match images with masks
            for idx, (base_name, img_rel_path) in enumerate(image_base_map.items()):
                mask_rel_path = mask_base_map.get(base_name)
                
                if mask_rel_path:
                    data_id = f"data_{idx}"
                    img_full_path = os.path.join(base_path, img_rel_path)
                    mask_full_path = os.path.join(base_path, mask_rel_path)
                    
                    labels = [AI_PixelLabel(
                        type="AI_PixelLabel",
                        image_url=[mask_full_path],
                        image_format=["image/tiff"],
                        class_=""
                    )]
                    
                    data.append(AI_EOTrainingData(
                        type="AI_EOTrainingData",
                        id=data_id,
                        data_url=[img_full_path],
                        labels=labels
                    ))
        else:
            print(f"Warning: Base path '{base_path}' does not exist. Creating placeholder entry.")
            # Create minimal data entry as placeholder
            data.append(AI_EOTrainingData(
                type="AI_EOTrainingData",
                id="data_0",
                data_url=[base_path or "./data"],
                labels=[AI_PixelLabel(
                    type="AI_PixelLabel",
                    image_url=["./labels"],
                    image_format=["image/tiff"],
                    class_=""
                )]
            ))

    # Ensure data is not empty
    if not data:
        raise ValueError("No data entries found. The distribution field may be empty or not parsed correctly.")

    # Build tasks
    tasks = [AI_EOTask(
        type="AI_EOTask",
        id="task_0",
        name="Segmentation Task",
        description="Semantic segmentation task for geospatial imagery.",
        input_type="image",
        output_type="mask",
        taskType="segmentation"
    )]

    # Create EOTrainingDataset object
    tdml_obj = EOTrainingDataset(
        type="AI_EOTrainingDataset",
        id=identifier,
        name=name,
        description=description,
        license=license_,
        providers=providers,
        created_time=created_time,
        updated_time=updated_time,
        version=version,
        tasks=tasks,
        classes=classes,
        bands=bands,
        data=data,
        extent=extent,
        amount_of_training_data=len(data),
        number_of_classes=len(classes)
    )

    # Write TDML JSON using pytdml
    pytdml.io.write_to_json(tdml_obj, tdml_output_path)
    print(f"TDML file written to {tdml_output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GeoCroissant JSON to TDML JSON using pytdml object model.")
    parser.add_argument("geocroissant_path", help="Path to input GeoCroissant JSON")
    parser.add_argument("tdml_output_path", help="Path to output TDML JSON")
    args = parser.parse_args()
    convert_geocroissant_to_tdml_objectmodel(args.geocroissant_path, args.tdml_output_path) 
