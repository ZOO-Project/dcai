import json
from datetime import datetime

# Create a proper GeoCroissant JSON-LD document according to the schema
geocroissant_json = {
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
    },
    "@type": "sc:Dataset",
    "name": "Landslide4sense",
    "description": "The Landslide4Sense dataset contains satellite imagery for landslide detection, consisting of 3799 training, 245 validation, and 800 test image patches. Each patch is a 128x128 pixel composite of 14 bands including Sentinel-2 multispectral data (B1-B12), slope data, and digital elevation model (DEM) from ALOS PALSAR. All bands are resampled to ~10m resolution and labeled pixel-wise for landslide presence.",
    "url": "https://huggingface.co/datasets/ibm-nasa-geospatial/Landslide4sense",
    "citeAs": "@dataset{Landslide4sense, title={Landslide4Sense: Reference Benchmark Data and Deep Learning Models for Landslide Detection}, author={Ghorbanzadeh, Omid and Xu, Yongjun and Zhao, Haokun and Wang, Junxi and Zhong, Yanfei and Zhao, Dong and Zang, Qiushi and Wang, Shimin and Zhang, Fang and Shi, Yiliang and others}, year={2022}, url={https://github.com/iarai/Landslide4Sense-2022}}",
    "datePublished": "2022-01-01",
    "version": "1.0",
    "license": "Unknown",
    "conformsTo": [
        "http://mlcommons.org/croissant/1.1",
        "http://mlcommons.org/croissant/geo/1.0"
    ],
    "identifier": "ibm-nasa-geospatial/Landslide4sense",
    "alternateName": ["Landslide4Sense-2022", "L4S"],
    "creator": {
        "@type": "Organization",
        "name": "IARAI - Institute for Advanced Research in Artificial Intelligence",
        "url": "https://github.com/iarai/Landslide4Sense-2022"
    },
    "keywords": [
        "Landslide4sense",
        "landslide detection",
        "natural disasters",
        "remote sensing",
        "satellite imagery",
        "Sentinel-2",
        "ALOS PALSAR",
        "DEM",
        "slope",
        "geospatial",
        "semantic segmentation",
        "English",
        "1K - 10K",
        "Image",
        "HDF5",
        "Datasets",
        "Croissant"
    ],
    "temporalCoverage": "2015-01-01/2022-01-01",
    "geocr:coordinateReferenceSystem": "EPSG:4326",
    "spatialCoverage": {
        "@type": "Place",
        "geo": {
            "@type": "GeoShape",
            "description": "Global coverage with focus on landslide-prone regions"
        }
    },
    "geocr:spatialResolution": {
        "@type": "QuantitativeValue",
        "value": 10.0,
        "unitText": "m"
    },
    "geocr:samplingStrategy": "Subsetted to 128x128 pixel windows covering landslide and non-landslide areas from global landslide-prone regions",
    "geocr:spectralBandMetadata": [
        {
            "@type": "geocr:SpectralBand",
            "name": "B1_Coastal_Aerosol",
            "description": "Sentinel-2 Band 1 - Coastal aerosol",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 443,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 20,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B2_Blue",
            "description": "Sentinel-2 Band 2 - Blue",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 490,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 65,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B3_Green",
            "description": "Sentinel-2 Band 3 - Green",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 560,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 35,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B4_Red",
            "description": "Sentinel-2 Band 4 - Red",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 665,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 30,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B5_Red_Edge_1",
            "description": "Sentinel-2 Band 5 - Vegetation Red Edge",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 705,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 15,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B6_Red_Edge_2",
            "description": "Sentinel-2 Band 6 - Vegetation Red Edge",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 740,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 15,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B7_Red_Edge_3",
            "description": "Sentinel-2 Band 7 - Vegetation Red Edge",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 783,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 20,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B8_NIR",
            "description": "Sentinel-2 Band 8 - Near Infrared",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 842,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 115,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B9_Water_Vapour",
            "description": "Sentinel-2 Band 9 - Water vapour",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 945,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 20,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B10_SWIR_Cirrus",
            "description": "Sentinel-2 Band 10 - SWIR - Cirrus",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 1375,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 30,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B11_SWIR_1",
            "description": "Sentinel-2 Band 11 - Short Wave Infrared 1",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 1610,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 90,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B12_SWIR_2",
            "description": "Sentinel-2 Band 12 - Short Wave Infrared 2",
            "geocr:centerWavelength": {
                "@type": "QuantitativeValue",
                "value": 2190,
                "unitText": "nm"
            },
            "geocr:bandwidth": {
                "@type": "QuantitativeValue",
                "value": 180,
                "unitText": "nm"
            }
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B13_Slope",
            "description": "ALOS PALSAR - Slope data derived from radar",
            "dataType": "Topographic"
        },
        {
            "@type": "geocr:SpectralBand",
            "name": "B14_DEM",
            "description": "ALOS PALSAR - Digital Elevation Model",
            "dataType": "Elevation"
        }
    ],
    "distribution": [
        {
            "@type": "cr:FileObject",
            "@id": "data_repo",
            "name": "data_repo",
            "description": "Landslide4sense dataset directory containing HDF5 files",
            "contentUrl": "./Landslide4sense",
            "encodingFormat": "local_directory",
            "sha256": "placeholder_checksum_for_directory"
        },
        {
            "@type": "cr:FileSet",
            "@id": "h5-files",
            "name": "h5-files",
            "description": "All HDF5 files containing images and masks",
            "containedIn": {
                "@id": "data_repo"
            },
            "encodingFormat": "application/x-hdf5",
            "includes": "**/*.h5"
        }
    ],
    "recordSet": [
        {
            "@type": "cr:RecordSet",
            "@id": "landslide4sense",
            "name": "landslide4sense",
            "description": "Landslide4sense dataset with 14-band satellite imagery and binary mask annotations for landslide detection.",
            "field": [
                {
                    "@type": "cr:Field",
                    "@id": "landslide4sense/image",
                    "name": "landslide4sense/image",
                    "description": "File path to HDF5 image file containing 14-band satellite imagery (Sentinel-2 B1-B12, Slope, DEM). Each image is 128x128 pixels with bands stacked as (height, width, bands). Data type: float64. Dataset name within HDF5: 'img'.",
                    "dataType": "sc:Text",
                    "source": {
                        "fileSet": {
                            "@id": "h5-files"
                        },
                        "extract": {
                            "fileProperty": "fullpath"
                        },
                        "transform": {
                            "regex": "images/.*/image_.*\\.h5$"
                        }
                    },
                    "geocr:bandConfiguration": {
                        "@type": "geocr:BandConfiguration",
                        "geocr:totalBands": 14,
                        "geocr:bandNamesList": [
                            "B1_Coastal_Aerosol", 
                            "B2_Blue", 
                            "B3_Green", 
                            "B4_Red", 
                            "B5_Red_Edge_1", 
                            "B6_Red_Edge_2", 
                            "B7_Red_Edge_3", 
                            "B8_NIR", 
                            "B9_Water_Vapour", 
                            "B10_SWIR_Cirrus", 
                            "B11_SWIR_1", 
                            "B12_SWIR_2", 
                            "B13_Slope", 
                            "B14_DEM"
                        ]
                    }
                },
                {
                    "@type": "cr:Field",
                    "@id": "landslide4sense/mask",
                    "name": "landslide4sense/mask",
                    "description": "File path to HDF5 mask file containing binary landslide annotations. Each mask is 128x128 pixels with values: 0 (non-landslide) and 1 (landslide). Data type: uint8. Dataset name within HDF5: 'mask'.",
                    "dataType": "sc:Text",
                    "source": {
                        "fileSet": {
                            "@id": "h5-files"
                        },
                        "extract": {
                            "fileProperty": "fullpath"
                        },
                        "transform": {
                            "regex": "annotations/.*/mask_.*\\.h5$"
                        }
                    },
                    "geocr:bandConfiguration": {
                        "@type": "geocr:BandConfiguration",
                        "geocr:totalBands": 1,
                        "geocr:bandNamesList": ["mask"]
                    }
                },
                {
                    "@type": "cr:Field",
                    "@id": "landslide4sense/split",
                    "name": "landslide4sense/split",
                    "description": "Dataset split indicator extracted from file path (train/validation/test)",
                    "dataType": "sc:Text",
                    "source": {
                        "fileSet": {
                            "@id": "h5-files"
                        },
                        "extract": {
                            "fileProperty": "fullpath"
                        },
                        "transform": {
                            "regex": "images/(train|validation|test)/"
                        }
                    }
                }
            ]
        }
    ]
}

# Write the GeoCroissant JSON-LD to file
with open("landslide4sense_geocroissant.json", "w") as f:
    json.dump(geocroissant_json, f, indent=2)
