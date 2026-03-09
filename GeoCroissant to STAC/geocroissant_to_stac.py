import json
from datetime import datetime
from pystac import Item, Asset, MediaType
from pystac.extensions.projection import ProjectionExtension

def geocroissant_to_stac(geocroissant_data):
    """Convert GeoCroissant metadata to STAC Item."""
    
    # Extract basic metadata
    item_id = geocroissant_data.get("name", "unknown").replace(" ", "_")
    title = geocroissant_data.get("name", "")
    description = geocroissant_data.get("description", "")
    license_info = geocroissant_data.get("license", "proprietary")
    keywords = geocroissant_data.get("keywords", [])
    
    # Parse spatial coverage (GeoCroissant format: "south west north east")
    spatial_coverage = geocroissant_data.get("spatialCoverage", {})
    geo_info = spatial_coverage.get("geo", {}) if isinstance(spatial_coverage, dict) else {}
    bbox_string = geo_info.get("box", "") if isinstance(geo_info, dict) else ""
    
    if bbox_string:
        coords = [float(x) for x in bbox_string.split()]
        south, west, north, east = coords
        bbox = [west, south, east, north]  # STAC format: [west, south, east, north]
        geometry = {
            "type": "Polygon",
            "coordinates": [[
                [west, south], [west, north], [east, north], [east, south], [west, south]
            ]]
        }
    else:
        # Default to global extent if no spatial coverage provided
        bbox = [-180, -90, 180, 90]
        geometry = {
            "type": "Polygon",
            "coordinates": [[
                [-180, -90], [-180, 90], [180, 90], [180, -90], [-180, -90]
            ]]
        }
    
    # Parse temporal coverage (ISO 8601: "start/end")
    temporal_coverage = geocroissant_data.get("temporalCoverage", "")
    if temporal_coverage and "/" in temporal_coverage:
        start_str, end_str = temporal_coverage.split("/")
        start_dt = datetime.fromisoformat(start_str)
        end_dt = datetime.fromisoformat(end_str)
        midpoint_dt = start_dt + (end_dt - start_dt) / 2
    else:
        start_dt = end_dt = None
        midpoint_dt = datetime.now()
    
    # Build STAC properties
    properties = {
        "title": title,
        "description": description,
        "license": license_info,
        "keywords": keywords,
    }
    
    if start_dt and end_dt:
        properties["start_datetime"] = start_dt.isoformat() + "Z"
        properties["end_datetime"] = end_dt.isoformat() + "Z"
    
    # Add GeoCroissant metadata
    crs = geocroissant_data.get("geocr:coordinateReferenceSystem")
    
    spatial_res = geocroissant_data.get("geocr:spatialResolution", {})
    if isinstance(spatial_res, dict) and spatial_res.get("value"):
        properties["gsd"] = float(spatial_res["value"])
    
    temporal_res = geocroissant_data.get("geocr:temporalResolution", {})
    if isinstance(temporal_res, dict) and temporal_res.get("value"):
        properties["geocr:temporalResolution"] = f"{temporal_res['value']} {temporal_res.get('unitText', '')}"
    
    sampling_strategy = geocroissant_data.get("geocr:samplingStrategy")
    if sampling_strategy:
        properties["geocr:samplingStrategy"] = sampling_strategy
    
    conforms_to = geocroissant_data.get("conformsTo", [])
    if conforms_to:
        properties["conformsTo"] = conforms_to
    
    # Create STAC Item
    item = Item(
        id=item_id,
        geometry=geometry,
        bbox=bbox,
        datetime=midpoint_dt,
        properties=properties
    )
    
    # Add projection extension if CRS present
    if crs and "EPSG:" in crs:
        proj_ext = ProjectionExtension.ext(item, add_if_missing=True)
        proj_ext.epsg = int(crs.replace("EPSG:", ""))
    
    # Process distribution to add assets
    distribution = geocroissant_data.get("distribution", [])
    for dist_item in distribution:
        item_type = dist_item.get("@type", "")
        content_url = dist_item.get("contentUrl", "")
        
        # Skip directory entries and file:// URLs
        if not content_url or content_url.startswith("file://"):
            continue
        if "directory" in dist_item.get("encodingFormat", "").lower():
            continue
        
        asset_id = dist_item.get("@id", dist_item.get("name", "asset")).replace(" ", "_").lower()
        encoding_format = dist_item.get("encodingFormat", "")
        
        # Determine media type
        if "tiff" in encoding_format.lower() or "tif" in encoding_format.lower():
            media_type = MediaType.GEOTIFF
        elif "json" in encoding_format.lower():
            media_type = MediaType.JSON
        elif "parquet" in encoding_format.lower():
            media_type = MediaType.PARQUET
        else:
            media_type = encoding_format
        
        # Determine roles
        roles = ["data"]
        if "FileSet" in item_type:
            roles.append("collection")
        
        asset = Asset(
            href=content_url,
            media_type=media_type,
            title=dist_item.get("description", dist_item.get("name", "")),
            roles=roles
        )
        
        # Add file pattern for FileSets
        includes = dist_item.get("includes")
        if includes:
            asset.extra_fields["file_pattern"] = includes
        
        item.add_asset(asset_id, asset)
    
    # Add spectral band metadata to GEOTIFF assets if present
    spectral_bands = geocroissant_data.get("geocr:spectralBandMetadata", [])
    if spectral_bands:
        raster_bands = []
        for band_info in spectral_bands:
            raster_band = {"name": band_info.get("name", "")}
            
            center_wl = band_info.get("geocr:centerWavelength", {})
            if isinstance(center_wl, dict) and center_wl.get("value"):
                raster_band["center_wavelength"] = float(center_wl["value"])
            
            bandwidth = band_info.get("geocr:bandwidth", {})
            if isinstance(bandwidth, dict) and bandwidth.get("value"):
                raster_band["bandwidth"] = float(bandwidth["value"])
            
            raster_bands.append(raster_band)
        
        # Apply to GEOTIFF assets
        for asset_key, asset in item.assets.items():
            if asset.media_type in [MediaType.GEOTIFF, MediaType.COG]:
                asset.extra_fields["raster:bands"] = raster_bands
    
    return item.to_dict()
