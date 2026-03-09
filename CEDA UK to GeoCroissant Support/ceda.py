import json

def cloud_product_to_geocroissant(products, stac_item):
    """Convert CEDA cloud products to valid GeoCroissant format"""
    # Get properties from STAC item
    properties = stac_item.stac_attributes.get('properties', {})
    bbox = stac_item.bbox
    geometry = stac_item.stac_attributes.get('geometry', {})
    item_id = stac_item.id
    
    # Extract CMIP6 metadata
    variable_name = properties.get('cmip6:variable_long_name', 'Unknown')
    variable_id = properties.get('cmip6:variable_id', 'tas')
    variable_units = properties.get('cmip6:variable_units', 'K')

    # Format bounding box for spatialCoverage
    if bbox and len(bbox) >= 4:
        bbox_str = f"{bbox[1]} {bbox[0]} {bbox[3]} {bbox[2]}"  # south west north east
    else:
        bbox_str = "-90 -180 90 180"

    # Build TTL-compliant GeoCroissant metadata
    croissant_metadata = {
        "@context": {
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
        },
        "@type": "sc:Dataset",
        "name": properties.get('title', item_id),
        "description": f"CMIP6 dataset for {variable_name}",
        "version": "1.0.0",
        "license": "CC-BY-4.0",
        "conformsTo": [
            "http://mlcommons.org/croissant/1.1",
            "http://mlcommons.org/croissant/geo/1.0"
        ],
        "citeAs": f"@dataset{{ceda_cmip6_{variable_id}, title={{CEDA CMIP6 {variable_name}}}, year={{2024}}, url={{https://catalogue.ceda.ac.uk/}}}}",
        "datePublished": "2024-01-01",
        
        # Standard spatial coverage using schema.org
        "spatialCoverage": {
            "@type": "Place",
            "geo": {
                "@type": "GeoShape",
                "box": bbox_str
            }
        },
        
        # GeoCroissant properties
        "geocr:coordinateReferenceSystem": "EPSG:4326",
        "geocr:spatialResolution": {
            "@type": "QuantitativeValue",
            "value": 1.0,
            "unitText": "degrees"
        },
        "geocr:temporalResolution": {
            "@type": "QuantitativeValue", 
            "value": 1,
            "unitText": "month"
        },
        
        # Temporal coverage
        "temporalCoverage": f"{properties.get('start_datetime', '2015-01-01')}/{properties.get('end_datetime', '2100-12-31')}",
        
        "keywords": [variable_id, "cmip6", "climate", "temperature", "ceda"],
        
        "distribution": []
    }
    
    # Add products as distribution items
    for product in products:
        if not hasattr(product, 'href'):
            continue
            
        # Extract asset name from product ID
        asset_name = product.id.split('-')[-1] if '-' in product.id else product.id
        
        # Determine encoding from URL
        product_url = product.href
        if product_url.endswith('.json'):
            encoding = 'application/json'
        elif product_url.endswith(('.nc', '.netcdf')):
            encoding = 'application/netcdf'
        elif product_url.endswith('.zarr'):
            encoding = 'application/zarr'
        else:
            encoding = "application/octet-stream"
        
        file_obj = {
            "@type": "cr:FileObject",
            "@id": asset_name,
            "name": asset_name,
            "contentUrl": product_url,
            "encodingFormat": encoding,
            "sha256": "placeholder"
        }
        croissant_metadata["distribution"].append(file_obj)
    
    # Add recordSet
    if croissant_metadata["distribution"]:
        croissant_metadata["recordSet"] = [
                {
                    "@type": "cr:RecordSet",
                    "name": "climate_data",
                    "field": [
                        {
                            "@type": "cr:Field",
                            "name": variable_id,
                            "description": variable_name,
                            "dataType": "sc:Float",
                            "source": {
                                "fileObject": {
                                    "@id": croissant_metadata["distribution"][0]["@id"]
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "latitude",
                            "description": "Latitude coordinate",
                            "dataType": "sc:Float", 
                            "source": {
                                "fileObject": {
                                    "@id": croissant_metadata["distribution"][0]["@id"]
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "longitude", 
                            "description": "Longitude coordinate",
                            "dataType": "sc:Float",
                            "source": {
                                "fileObject": {
                                    "@id": croissant_metadata["distribution"][0]["@id"]
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "time",
                            "description": "Time coordinate",
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": croissant_metadata["distribution"][0]["@id"]
                                }
                            }
                        }
                    ]
                }
            ]
    
    return croissant_metadata