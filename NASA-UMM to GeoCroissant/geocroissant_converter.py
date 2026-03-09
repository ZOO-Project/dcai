#!/usr/bin/env python3
"""
NASA UMM-G to GeoCroissant Converter

This script converts NASA UMM-G JSON to GeoCroissant format
"""

import json
import re
from typing import Dict, List, Any, Optional
from datetime import datetime

class CompleteNASAUMMGToGeoCroissantConverter:
    """Converter that uses only properties defined in the official TTL schemas."""
    
    def __init__(self):
        self.setup_context()
    
    def setup_context(self):
        """Setup the JSON-LD context for GeoCroissant."""
        self.context = {
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
            "data": {
                "@id": "cr:data",
                "@type": "@json"
            },
            "examples": {
                "@id": "cr:examples",
                "@type": "@json"
            },
            "dataBiases": "cr:dataBiases",
            "dataCollection": "cr:dataCollection",
            "dataType": {
                "@id": "cr:dataType",
                "@type": "@vocab"
            },
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
    
    def create_dataset_structure(self, meta: Dict[str, Any], umm: Dict[str, Any]) -> Dict[str, Any]:
        """Create the main Dataset structure using only TTL-defined properties."""
        concept_id = meta.get('concept-id')
        granule_ur = umm.get('GranuleUR', 'HLS_Sentinel2_Dataset')
        
        dataset = {
            "@context": self.context,
            "@type": "sc:Dataset",
            "name": granule_ur,
            "description": umm.get('CollectionReference', {}).get('EntryTitle', 'HLS Sentinel-2 satellite imagery dataset'),
            "url": f"https://cmr.earthdata.nasa.gov/search/concepts/{concept_id}.html",
            "datePublished": meta.get('revision-date'),
            "version": str(meta.get('revision-id', '1.0')),
            "license": "https://creativecommons.org/publicdomain/zero/1.0/",
            "citeAs": f"{granule_ur}. NASA EOSDIS Land Processes Distributed Active Archive Center. https://cmr.earthdata.nasa.gov/search/concepts/{concept_id}.html",
            "conformsTo": [
                "http://mlcommons.org/croissant/1.1",
                "http://mlcommons.org/croissant/geo/1.0"
            ]
        }
        
        # Add GeoCroissant properties at dataset level
        self.add_geospatial_properties(dataset, umm)
        self.add_temporal_properties(dataset, umm)
        self.add_band_properties(dataset, umm)
        self.add_instrument_properties(dataset, umm)
        self.add_sampling_properties(dataset, umm)
        
        # Add distribution
        self.add_distribution(dataset, umm)
        
        # Add record set
        dataset["recordSet"] = [self.create_record(meta, umm)]
        
        return dataset
    
    def create_record(self, meta: Dict[str, Any], umm: Dict[str, Any]) -> Dict[str, Any]:
        """Create a single record within the RecordSet."""
        record = {
            "@type": "cr:RecordSet",
            "@id": meta.get('concept-id'),
            "name": umm.get('GranuleUR'),
            "description": umm.get('CollectionReference', {}).get('EntryTitle')
        }
        
        return record
    
    def add_geospatial_properties(self, dataset: Dict[str, Any], umm: Dict[str, Any]):
        """Add geospatial properties using only TTL-defined properties."""
        additional_attrs = umm.get('AdditionalAttributes', [])
        
        # coordinateReferenceSystem (valid property)
        crs_code = self.find_additional_attribute(additional_attrs, 'HORIZONTAL_CS_CODE')
        if crs_code:
            dataset["geocr:coordinateReferenceSystem"] = crs_code
        
        # spatialResolution (valid property - as QuantitativeValue)
        spatial_resolution = self.find_additional_attribute(additional_attrs, 'SPATIAL_RESOLUTION')
        if spatial_resolution:
            dataset["geocr:spatialResolution"] = {
                "@type": "sc:QuantitativeValue",
                "value": float(spatial_resolution),
                "unitText": "m"
            }
        
        # spatialCoverage (schema.org property)
        spatial_extent = umm.get('SpatialExtent', {})
        if spatial_extent:
            horizontal_domain = spatial_extent.get('HorizontalSpatialDomain', {})
            geometry = horizontal_domain.get('Geometry', {})
            polygons = geometry.get('GPolygons', [])
            
            if polygons:
                points = polygons[0].get('Boundary', {}).get('Points', [])
                if points:
                    bbox = self.calculate_bounding_box(points)
                    if bbox:
                        dataset["spatialCoverage"] = {
                            "@type": "sc:Place",
                            "geo": {
                                "@type": "sc:GeoShape",
                                "box": f"{bbox['south']} {bbox['west']} {bbox['north']} {bbox['east']}"
                            }
                        }
    
    def add_temporal_properties(self, dataset: Dict[str, Any], umm: Dict[str, Any]):
        """Add temporal properties using only TTL-defined properties."""
        temporal_extent = umm.get('TemporalExtent', {})
        if temporal_extent:
            range_datetime = temporal_extent.get('RangeDateTime', {})
            if range_datetime:
                start = range_datetime.get('BeginningDateTime')
                end = range_datetime.get('EndingDateTime')
                if start and end:
                    # temporalCoverage (schema.org property)
                    dataset["temporalCoverage"] = f"{start}/{end}"
    
    def add_band_properties(self, dataset: Dict[str, Any], umm: Dict[str, Any]):
        """Add band configuration and spectral metadata using TTL-defined properties."""
        # bandConfiguration (valid property)
        dataset["geocr:bandConfiguration"] = {
            "@type": "geocr:BandConfiguration",
            "geocr:totalBands": 13,
            "geocr:bandNamesList": ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B10", "B11", "B12"]
        }
        
        # spectralBandMetadata (valid property)
        band_info = {
            "B01": {"wavelength": 443, "bandwidth": 65, "description": "Coastal aerosol"},
            "B02": {"wavelength": 490, "bandwidth": 65, "description": "Blue"},
            "B03": {"wavelength": 560, "bandwidth": 60, "description": "Green"},
            "B04": {"wavelength": 665, "bandwidth": 30, "description": "Red"},
            "B05": {"wavelength": 705, "bandwidth": 15, "description": "Red edge 1"},
            "B06": {"wavelength": 740, "bandwidth": 15, "description": "Red edge 2"},
            "B07": {"wavelength": 783, "bandwidth": 20, "description": "Red edge 3"},
            "B08": {"wavelength": 842, "bandwidth": 115, "description": "NIR"},
            "B8A": {"wavelength": 865, "bandwidth": 20, "description": "NIR narrow"},
            "B09": {"wavelength": 945, "bandwidth": 20, "description": "Water vapour"},
            "B10": {"wavelength": 1375, "bandwidth": 30, "description": "SWIR cirrus"},
            "B11": {"wavelength": 1610, "bandwidth": 90, "description": "SWIR 1"},
            "B12": {"wavelength": 2190, "bandwidth": 180, "description": "SWIR 2"}
        }
        
        spectral_bands = []
        for band_name, info in band_info.items():
            spectral_bands.append({
                "@type": "geocr:SpectralBand",
                "name": band_name,
                "description": info["description"],
                "geocr:centerWavelength": {
                    "@type": "sc:QuantitativeValue",
                    "value": info["wavelength"],
                    "unitText": "nm"
                },
                "geocr:bandwidth": {
                    "@type": "sc:QuantitativeValue",
                    "value": info["bandwidth"],
                    "unitText": "nm"
                }
            })
        
        dataset["geocr:spectralBandMetadata"] = spectral_bands
    
    def add_instrument_properties(self, dataset: Dict[str, Any], umm: Dict[str, Any]):
        """Add instrument and observatory properties using TTL-defined properties."""
        platforms = umm.get('Platforms', [])
        if platforms:
            platform = platforms[0]
            instruments = platform.get('Instruments', [])
            
            if instruments:
                # solarInstrumentCharacteristics (valid for space weather, but can be used here)
                dataset["geocr:solarInstrumentCharacteristics"] = {
                    "@type": "geocr:SolarInstrumentCharacteristics",
                    "geocr:observatory": platform.get('ShortName'),
                    "geocr:instrument": instruments[0].get('ShortName')
                }
    
    def add_sampling_properties(self, dataset: Dict[str, Any], umm: Dict[str, Any]):
        """Add sampling strategy using TTL-defined properties."""
        additional_attrs = umm.get('AdditionalAttributes', [])
        
        # samplingStrategy (valid property)
        spatial_coverage = self.find_additional_attribute(additional_attrs, 'SPATIAL_COVERAGE')
        if spatial_coverage:
            dataset["geocr:samplingStrategy"] = f"Spatial coverage: {spatial_coverage}%"
    
    def add_distribution(self, dataset: Dict[str, Any], umm: Dict[str, Any]):
        """Add distribution information using schema.org properties."""
        distributions = []
        
        # Get all related URLs
        related_urls = umm.get('RelatedUrls', [])
        
        for url_info in related_urls:
            url = url_info.get('URL', '')
            url_type = url_info.get('Type', '')
            subtype = url_info.get('Subtype', '')
            description = url_info.get('Description', '')
            
            # Determine encoding format based on URL or type
            encoding_format = self.determine_encoding_format(url, url_type, subtype)
            
            distribution = {
                "@type": "cr:FileObject",
                "name": url.split('/')[-1] or "data_file",
                "contentUrl": url,
                "encodingFormat": encoding_format,
                "sha256": "https://github.com/mlcommons/croissant/issues/80",
                "description": description or f"Download {url.split('/')[-1]}"
            }
            
            distributions.append(distribution)
        
        if distributions:
            dataset["distribution"] = distributions
    
    def convert_polygon_to_wkt(self, points: List[Dict[str, float]]) -> str:
        """Convert polygon points to WKT format."""
        if not points:
            return ""
        
        coords = []
        for point in points:
            lon = point.get('Longitude', 0)
            lat = point.get('Latitude', 0)
            coords.append(f"{lon} {lat}")
        
        if coords and coords[0] != coords[-1]:
            coords.append(coords[0])
        
        return f"POLYGON(({', '.join(coords)}))"
    
    def calculate_bounding_box(self, points: List[Dict[str, float]]) -> Dict[str, float]:
        """Calculate bounding box from polygon points."""
        if not points:
            return {}
        
        lons = [p.get('Longitude', 0) for p in points]
        lats = [p.get('Latitude', 0) for p in points]
        
        return {
            "west": min(lons),
            "south": min(lats),
            "east": max(lons),
            "north": max(lats)
        }
    
    def find_additional_attribute(self, attributes: List[Dict], name: str) -> Optional[str]:
        """Find value of additional attribute by name."""
        for attr in attributes:
            if attr.get('Name') == name:
                values = attr.get('Values', [])
                return values[0] if values else None
        return None
    
    def find_additional_attribute_values(self, attributes: List[Dict], name: str) -> List[str]:
        """Find all values of additional attribute by name."""
        for attr in attributes:
            if attr.get('Name') == name:
                return attr.get('Values', [])
        return []
    
    def determine_encoding_format(self, url: str, url_type: str, subtype: str) -> str:
        """Determine the encoding format based on URL and type."""
        if url.endswith('.tif') or url.endswith('.tiff'):
            return "image/tiff"
        elif url.endswith('.jpg') or url.endswith('.jpeg'):
            return "image/jpeg"
        elif url.endswith('.json'):
            return "application/json"
        elif url.endswith('.xml'):
            return "application/xml"
        elif url.endswith('.hdf') or url.endswith('.h5'):
            return "application/x-hdf"
        elif url.endswith('.nc'):
            return "application/x-netcdf"
        elif 's3credentials' in url:
            return "application/octet-stream"
        else:
            return "application/octet-stream"
    
    def convert_to_complete_geocroissant(self, ummg_data: Dict[str, Any]) -> Dict[str, Any]:
        """Main conversion method - using only TTL-defined properties."""
        # Extract main sections
        meta = ummg_data.get('meta', {})
        umm = ummg_data.get('umm', {})
        
        # Create the complete GeoCroissant structure
        return self.create_dataset_structure(meta, umm)

def main():
    """Main function to demonstrate conversion using only TTL-defined properties."""
    
    # Load the NASA UMM-G JSON
    with open('nasa_ummg_h.json', 'r') as f:
        ummg_data = json.load(f)
    
    # Convert to GeoCroissant using only TTL properties
    converter = CompleteNASAUMMGToGeoCroissantConverter()
    geocroissant_data = converter.convert_to_complete_geocroissant(ummg_data)
    
    # Save the converted data
    with open('geocroissant_output.json', 'w') as f:
        json.dump(geocroissant_data, f, indent=2)
    
    print("Conversion completed using only TTL-defined properties!")
    print(f"Input: nasa_ummg_h.json")
    print(f"Output: geocroissant_output.json")
    
    # Print statistics
    print("\nGeoCroissant Properties Used (from TTL):")
    ttl_properties = [
        "geocr:coordinateReferenceSystem",
        "geocr:spatialResolution",
        "geocr:temporalResolution",
        "geocr:bandConfiguration",
        "geocr:spectralBandMetadata",
        "geocr:samplingStrategy",
        "geocr:solarInstrumentCharacteristics"
    ]
    
    for prop in ttl_properties:
        if prop in str(geocroissant_data):
            print(f"  ✓ {prop}")
    
    print("\nSchema.org Properties Used:")
    schema_properties = [
        "spatialCoverage",
        "temporalCoverage",
        "distribution"
    ]
    
    for prop in schema_properties:
        if prop in str(geocroissant_data):
            print(f"  ✓ {prop}")
    

if __name__ == "__main__":
    main()
