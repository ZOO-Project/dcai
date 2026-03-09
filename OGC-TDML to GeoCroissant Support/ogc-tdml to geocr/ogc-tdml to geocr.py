import json
import re
import os
from datetime import datetime
import pytdml.io

def tdml_to_geocroissant(tdml_path, output_path, dataset_url=None):
    """
    Convert OGC-TDML JSON to GeoCroissant JSON-LD format.
    Fully compliant with Croissant 1.1 and GeoCroissant 1.0 specifications.
    
    Args:
        tdml_path: Path to input TDML JSON file
        output_path: Path to output GeoCroissant JSON file
        dataset_url: Optional URL for the dataset landing page
    """
    tdml = pytdml.io.read_from_json(tdml_path)
    
    # Sanitize the name for forbidden characters
    sanitized_name = re.sub(r'[^A-Za-z0-9_-]', '_', tdml.name)
    
    # Build proper @context with all required namespaces
    context = {
        "@language": "en",
        "@vocab": "https://schema.org/",
        "citeAs": "cr:citeAs",
        "column": "cr:column",
        "conformsTo": "dct:conformsTo",
        "cr": "http://mlcommons.org/croissant/",
        "geocr": "http://mlcommons.org/croissant/geo/",
        "rai": "http://mlcommons.org/croissant/RAI/",
        "dct": "http://purl.org/dc/terms/",
        "sc": "https://schema.org/",
        "data": {"@id": "cr:data", "@type": "@json"},
        "examples": {"@id": "cr:examples", "@type": "@json"},
        "dataBiases": "cr:dataBiases",
        "dataCollection": "cr:dataCollection",
        "dataType": {"@id": "cr:dataType", "@type": "@vocab"},
        "extract": "cr:extract",
        "field": "cr:field",
        "fileProperty": "cr:fileProperty",
        "fileObject": "cr:fileObject",
        "fileSet": "cr:fileSet",
        "format": "cr:format",
        "includes": "cr:includes",
        "isLiveDataset": "cr:isLiveDataset",
        "jsonPath": "cr:jsonPath",
        "key": "cr:key",
        "md5": "cr:md5",
        "parentField": "cr:parentField",
        "path": "cr:path",
        "personalSensitiveInformation": "cr:personalSensitiveInformation",
        "recordSet": "cr:recordSet",
        "references": "cr:references",
        "regex": "cr:regex",
        "repeated": "cr:repeated",
        "replace": "cr:replace",
        "samplingRate": "cr:samplingRate",
        "separator": "cr:separator",
        "source": "cr:source",
        "subField": "cr:subField",
        "transform": "cr:transform"
    }
    
    # Extract keywords from task types and classes
    keywords = [tdml.name]
    if hasattr(tdml, "tasks") and tdml.tasks:
        for task in tdml.tasks:
            if hasattr(task, "taskType"):
                keywords.append(task.taskType)
    if hasattr(tdml, "classes") and tdml.classes:
        keywords.extend([c.key for c in tdml.classes[:5]])  # Add first 5 class names
    
    # Build spectralBandMetadata from TDML bands
    spectral_bands = []
    if hasattr(tdml, "bands") and tdml.bands:
        for band in tdml.bands:
            band_entry = {
                "@type": "geocr:SpectralBand",
                "name": band.description if hasattr(band, "description") else band.name[0].code if band.name else "Unknown"
            }
            if hasattr(band, "units") and band.units:
                # Parse units like "865nm" -> value: 865, unitText: "nm"
                units_str = str(band.units)
                match = re.match(r'(\d+\.?\d*)([a-zA-Z]+)', units_str)
                if match:
                    band_entry["geocr:centerWavelength"] = {
                        "@type": "QuantitativeValue",
                        "value": float(match.group(1)),
                        "unitText": match.group(2)
                    }
            spectral_bands.append(band_entry)
    
    # Build spatialCoverage from extent
    spatial_coverage = None
    if hasattr(tdml, "extent") and tdml.extent:
        extent = tdml.extent
        spatial_coverage = {
            "@type": "Place",
            "geo": {
                "@type": "GeoShape",
                "box": f"{extent[1]} {extent[0]} {extent[3]} {extent[2]}"  # minLat minLon maxLat maxLon
            }
        }
    
    # Determine if we have local files or URLs
    has_local_files = False
    if hasattr(tdml, "data") and len(tdml.data) > 0:
        first_url = tdml.data[0].data_url[0] if tdml.data[0].data_url else ""
        has_local_files = not first_url.startswith(('http://', 'https://'))
    
    # Build distribution with proper FileObject/FileSet structure
    distribution = []
    
    if has_local_files:
        # For local files, use FileObject + FileSet structure
        # Extract base directory from first file path
        if tdml.data and tdml.data[0].data_url:
            first_path = tdml.data[0].data_url[0]
            # Find common parent directory
            base_dir = os.path.dirname(os.path.dirname(first_path))  # Go up to parent of training/validation
            
            distribution.append({
                "@type": "cr:FileObject",
                "@id": "data_repo",
                "name": "data_repo",
                "description": "Directory containing the dataset files",
                "contentUrl": base_dir,
                "encodingFormat": "local_directory",
                "sha256": "placeholder_hash_for_directory"
            })
            
            distribution.append({
                "@type": "cr:FileSet",
                "@id": "tiff-files",
                "name": "tiff-files",
                "description": "All TIFF files (images and masks).",
                "containedIn": {"@id": "data_repo"},
                "encodingFormat": "image/tiff",
                "includes": "**/*.tif"
            })
    else:
        # For remote URLs, list individual FileObjects
        for idx, d in enumerate(tdml.data[:100]):  # Limit to first 100 to avoid huge files
            if d.data_url:
                distribution.append({
                    "@type": "cr:FileObject",
                    "@id": f"image_{idx}",
                    "name": os.path.basename(d.data_url[0]),
                    "contentUrl": d.data_url[0],
                    "encodingFormat": "image/tiff"
                })
            if d.labels and hasattr(d.labels[0], "image_url"):
                mask_url = d.labels[0].image_url[0]
                distribution.append({
                    "@type": "cr:FileObject",
                    "@id": f"mask_{idx}",
                    "name": os.path.basename(mask_url),
                    "contentUrl": mask_url,
                    "encodingFormat": d.labels[0].image_format[0] if hasattr(d.labels[0], "image_format") else "image/tiff"
                })
    
    # Build recordSet with field definitions
    record_set = []
    
    # Determine band names
    band_names = []
    if spectral_bands:
        band_names = [band.get("name", f"Band_{i+1}") for i, band in enumerate(spectral_bands)]
    
    fields = []
    
    # Image field
    image_field = {
        "@type": "cr:Field",
        "@id": f"{sanitized_name}/image",
        "name": f"{sanitized_name}/image",
        "description": "Satellite imagery with multiple spectral bands",
        "dataType": "sc:Text"
    }
    
    if has_local_files:
        image_field["source"] = {
            "fileSet": {"@id": "tiff-files"},
            "extract": {"fileProperty": "fullpath"},
            "transform": {"regex": ".*_merged\\.tif$"}
        }
    
    if band_names:
        image_field["geocr:bandConfiguration"] = {
            "@type": "geocr:BandConfiguration",
            "geocr:totalBands": len(band_names),
            "geocr:bandNameList": band_names
        }
    
    fields.append(image_field)
    
    # Mask/label field
    mask_field = {
        "@type": "cr:Field",
        "@id": f"{sanitized_name}/mask",
        "name": f"{sanitized_name}/mask",
        "description": "Mask annotations with values representing different classes",
        "dataType": "sc:Text"
    }
    
    if has_local_files:
        mask_field["source"] = {
            "fileSet": {"@id": "tiff-files"},
            "extract": {"fileProperty": "fullpath"},
            "transform": {"regex": ".*\\.mask\\.tif$"}
        }
    
    mask_field["geocr:bandConfiguration"] = {
        "@type": "geocr:BandConfiguration",
        "geocr:totalBands": 1,
        "geocr:bandNameList": ["mask"]
    }
    
    fields.append(mask_field)
    
    record_set.append({
        "@type": "cr:RecordSet",
        "@id": sanitized_name,
        "name": sanitized_name,
        "description": f"{tdml.description}",
        "field": fields
    })
    
    # Build creator structure
    creator = {
        "@type": "Organization",
        "name": tdml.providers[0] if tdml.providers else "Unknown"
    }
    if dataset_url:
        creator["url"] = dataset_url
    
    # Build citeAs (BibTeX format)
    cite_as = f"""@dataset{{{sanitized_name},
    title={{{tdml.name}}},
    author={{{tdml.providers[0] if tdml.providers else "Unknown"}}},
    year={{{tdml.created_time[:4] if tdml.created_time else datetime.now().year}}},
    publisher={{TDML Dataset}}"""
    
    if tdml.id:
        cite_as += f",\n    doi={{{tdml.id}}}"
    cite_as += "\n}"
    
    # Construct the GeoCroissant document
    geocroissant = {
        "@context": context,
        "@type": "sc:Dataset",
        "name": tdml.name,
        "description": tdml.description,
        "conformsTo": [
            "http://mlcommons.org/croissant/1.1",
            "http://mlcommons.org/croissant/geo/1.0"
        ],
        "identifier": tdml.id,
        "license": tdml.license or "Unknown",
        "url": dataset_url or f"file://./{sanitized_name}",
        "creator": creator,
        "datePublished": tdml.created_time or datetime.now().isoformat()[:10],
        "version": getattr(tdml, "version", "1.0"),
        "keywords": keywords,
        "citeAs": cite_as,
        "distribution": distribution,
        "recordSet": record_set
    }
    
    # Add optional GeoCroissant properties
    if spectral_bands:
        geocroissant["geocr:spectralBandMetadata"] = spectral_bands
    
    if spatial_coverage:
        geocroissant["spatialCoverage"] = spatial_coverage
    
    if hasattr(tdml, "created_time"):
        geocroissant["dateCreated"] = tdml.created_time
    
    if hasattr(tdml, "updated_time"):
        geocroissant["dateModified"] = tdml.updated_time
    
    # Write output
    with open(output_path, "w") as f:
        json.dump(geocroissant, f, indent=2, ensure_ascii=False)

# Example usage:
tdml_to_geocroissant("hls_burn_scars_tdml.json", "hls_burn_scars_geo-croissant.json")