#!/usr/bin/env python3
"""
Dynamic GeoCroissant metadata generator for SDO NetCDF datasets.
Automatically extracts metadata from NetCDF files and generates compliant GeoCroissant JSON-LD.

Conforms to:
- Croissant 1.1: http://mlcommons.org/croissant/1.1
- GeoCroissant 1.0: http://mlcommons.org/croissant/geo/1.0
"""

import json
import netCDF4
import glob
import os
from datetime import datetime
from pathlib import Path


class SDOGeoCroissantConverter:
    """Convert SDO NetCDF data to GeoCroissant metadata format."""
    
    # AIA wavelength metadata (from SDO documentation)
    AIA_BANDS = {
        "aia94": {"wavelength": 94, "bandwidth": 3, "ion": "Fe XVIII", "temp": "~6.3 MK", "region": "corona and flare plasma"},
        "aia131": {"wavelength": 131, "bandwidth": 7, "ion": "Fe VIII/XXI", "temp": "~0.4 MK and ~11 MK", "region": "cool and hot plasma"},
        "aia171": {"wavelength": 171, "bandwidth": 6, "ion": "Fe IX", "temp": "~0.6 MK", "region": "quiet corona"},
        "aia193": {"wavelength": 193, "bandwidth": 6, "ion": "Fe XII/XXIV", "temp": "~1.5 MK and ~20 MK", "region": "corona and hot flare plasma"},
        "aia211": {"wavelength": 211, "bandwidth": 6, "ion": "Fe XIV", "temp": "~2 MK", "region": "active regions"},
        "aia304": {"wavelength": 304, "bandwidth": 4, "ion": "He II", "temp": "~0.05 MK", "region": "chromosphere and transition region"},
        "aia335": {"wavelength": 335, "bandwidth": 4, "ion": "Fe XVI", "temp": "~2.5 MK", "region": "active region corona"},
        "aia1600": {"wavelength": 1600, "bandwidth": 55, "ion": "C IV continuum", "temp": "~0.005 MK", "region": "upper photosphere and transition region"}
    }
    
    def __init__(self, data_dir):
        """Initialize converter with data directory."""
        self.data_dir = Path(data_dir)
        self.nc_files = sorted(glob.glob(str(self.data_dir / "*.nc")))
        
        if not self.nc_files:
            raise ValueError(f"No NetCDF files found in {data_dir}")
        
        # Read metadata from first file as representative
        self.sample_nc = netCDF4.Dataset(self.nc_files[0], 'r')
        self.variables = self._extract_variables()
        
    def _extract_variables(self):
        """Extract variable information from NetCDF file."""
        variables = {}
        
        for var_name, var in self.sample_nc.variables.items():
            var_info = {
                "name": var_name,
                "shape": var.shape,
                "dtype": str(var.dtype),
                "dimensions": var.dimensions,
                "attributes": {}
            }
            
            # Extract important attributes
            for attr in ['unit', 'description', 't_obs', 'qflag']:
                if hasattr(var, attr):
                    var_info['attributes'][attr] = str(getattr(var, attr))
            
            variables[var_name] = var_info
        
        return variables
    
    def _get_temporal_coverage(self):
        """Extract temporal coverage from filenames."""
        # Extract dates from filenames like 20110120_0100.nc
        dates = []
        for f in self.nc_files:
            basename = os.path.basename(f)
            date_str = basename.split('_')[0]
            if len(date_str) == 8:  # YYYYMMDD
                dates.append(f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}")
        
        if dates:
            return f"{min(dates)}/{max(dates)}"
        return "Unknown"
    
    def _generate_spectral_band_metadata(self):
        """Generate spectral band metadata for AIA channels."""
        bands = []
        
        for var_name in self.variables:
            if var_name.startswith('aia') and var_name in self.AIA_BANDS:
                band_info = self.AIA_BANDS[var_name]
                bands.append({
                    "@type": "geocr:SpectralBand",
                    "name": f"AIA_{band_info['wavelength']}",
                    "description": f"AIA {band_info['wavelength']}Å - {band_info['ion']} emission, {band_info['temp']} {band_info['region']}",
                    "geocr:centerWavelength": {
                        "@type": "QuantitativeValue",
                        "value": band_info['wavelength'],
                        "unitText": "Angstrom"
                    },
                    "geocr:bandwidth": {
                        "@type": "QuantitativeValue",
                        "value": band_info['bandwidth'],
                        "unitText": "Angstrom"
                    }
                })
        
        return bands
    
    def _generate_field_definitions(self):
        """Generate field definitions for all variables."""
        fields = []
        
        for var_name, var_info in self.variables.items():
            # Get description from attributes or generate one
            description = var_info['attributes'].get('description', '')
            unit = var_info['attributes'].get('unit', '')
            
            # Enhanced descriptions
            if var_name.startswith('aia'):
                if var_name in self.AIA_BANDS:
                    band_info = self.AIA_BANDS[var_name]
                    description = f"AIA {band_info['wavelength']}Å channel - {band_info['ion']} emission for {band_info['region']} at {band_info['temp']}."
            elif var_name == 'hmi_m':
                description = "HMI line-of-sight magnetogram - measures the magnetic field component along the line of sight from observer to Sun."
            elif var_name == 'hmi_bx':
                description = "HMI vector magnetic field X-component (east-west direction in heliographic coordinates)."
            elif var_name == 'hmi_by':
                description = "HMI vector magnetic field Y-component (north-south direction in heliographic coordinates)."
            elif var_name == 'hmi_bz':
                description = "HMI vector magnetic field Z-component (radial direction normal to solar surface)."
            elif var_name == 'hmi_v':
                description = "HMI line-of-sight Doppler velocity - measures plasma motion via Doppler shift of spectral lines. Positive values indicate motion away from observer."
            
            # Add units and dimensions
            if unit:
                description += f" Data units: {unit}."
            if len(var_info['shape']) == 2:
                description += f" Image dimensions: {var_info['shape'][0]}x{var_info['shape'][1]} pixels."
            
            # Dynamically determine number of bands from variable shape
            # Shape: (y, x) = 1 band (2D array)
            # Shape: (y, x, bands) = multiple bands (3D array)
            # Shape: (bands, y, x) = multiple bands (3D array with bands first)
            num_bands = 1
            band_names = [var_name]
            
            if len(var_info['shape']) == 3:
                # 3D array - detect if last or first dimension is bands
                # Typically bands is the smallest dimension
                dims = var_info['shape']
                if dims[0] < dims[1] and dims[0] < dims[2]:  # bands first: (bands, y, x)
                    num_bands = dims[0]
                elif dims[2] < dims[0] and dims[2] < dims[1]:  # bands last: (y, x, bands)
                    num_bands = dims[2]
                
                # Generate band names if multiple bands
                if num_bands > 1:
                    band_names = [f"{var_name}_band{i+1}" for i in range(num_bands)]
            
            field_def = {
                "@type": "cr:Field",
                "@id": f"sdo/{var_name}",
                "name": var_name,
                "description": description.strip(),
                "dataType": "sc:Float",
                "source": {
                    "fileSet": {
                        "@id": "nc-files"
                    },
                    "extract": {
                        "column": var_name
                    }
                },
                "geocr:bandConfiguration": {
                    "@type": "geocr:BandConfiguration",
                    "geocr:totalBands": num_bands,
                    "geocr:bandNamesList": band_names
                }
            }
            
            fields.append(field_def)
        
        return fields
    
    def generate_geocroissant(self, output_file="sdo_geocroissant.json"):
        """Generate complete GeoCroissant metadata compliant with Croissant 1.1 and GeoCroissant 1.0."""
        
        # Get AIA wavelength list
        aia_wavelengths = [f"{self.AIA_BANDS[v]['wavelength']}Å" 
                          for v in self.variables if v in self.AIA_BANDS]
        
        geocroissant = {
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
            },
            "@type": "sc:Dataset",
            "name": "SDO Multi-Instrument Solar Observations",
            "description": f"Solar Dynamics Observatory (SDO) multi-instrument dataset containing synchronized observations from the Atmospheric Imaging Assembly (AIA) and Helioseismic and Magnetic Imager (HMI). Each timestep includes {len([v for v in self.variables if v.startswith('aia')])} AIA extreme ultraviolet (EUV) wavelength channels ({', '.join(aia_wavelengths)}) and {len([v for v in self.variables if v.startswith('hmi')])} HMI magnetic field measurements. All data are provided as {self.sample_nc.dimensions['x'].size}x{self.sample_nc.dimensions['y'].size} pixel full-disk images in NetCDF4 format, capturing the dynamic solar atmosphere and magnetic field evolution for space weather research and solar physics studies.",
            "url": "https://sdo.gsfc.nasa.gov/",
            "citeAs": "@article{Pesnell2012, title={The Solar Dynamics Observatory (SDO)}, author={Pesnell, W. Dean and Thompson, B. J. and Chamberlin, P. C.}, journal={Solar Physics}, volume={275}, pages={3--15}, year={2012}, doi={10.1007/s11207-011-9841-3}}",
            "datePublished": "2010-02-11",
            "version": "1.0",
            "license": "CC0-1.0",
            "conformsTo": [
                "http://mlcommons.org/croissant/1.1",
                "http://mlcommons.org/croissant/geo/1.0"
            ],
            "identifier": "nasa-sdo/aia-hmi-core",
            "alternateName": ["SDO", "Solar Dynamics Observatory"],
            "creator": {
                "@type": "Organization",
                "name": "NASA Solar Dynamics Observatory",
                "url": "https://sdo.gsfc.nasa.gov/"
            },
            "keywords": [
                "SDO", "Solar Dynamics Observatory", "AIA", "HMI",
                "space weather", "heliophysics", "solar physics",
                "magnetic field", "extreme ultraviolet", "EUV",
                "magnetogram", "Doppler velocity", "solar atmosphere",
                "NetCDF", "time series"
            ],
            "temporalCoverage": self._get_temporal_coverage(),
            "geocr:temporalResolution": {
                "@type": "QuantitativeValue",
                "value": 12,
                "unitText": "minutes"
            },
            "geocr:coordinateReferenceSystem": "Helioprojective-Cartesian (HPC)",
            "spatialCoverage": {
                "@type": "Place",
                "geo": {
                    "@type": "GeoShape",
                    "description": "Full solar disk coverage from Sun-Earth L1 Lagrange point"
                }
            },
            "geocr:spatialResolution": {
                "@type": "QuantitativeValue",
                "value": 0.6,
                "unitText": "arcsec/pixel"
            },
            "geocr:samplingStrategy": "Full-disk synchronized observations from AIA and HMI instruments, temporally aligned to 12-minute cadence snapshots",
            "geocr:solarInstrumentCharacteristics": {
                "@type": "geocr:SolarInstrumentCharacteristics",
                "geocr:observatory": "Solar Dynamics Observatory (SDO)",
                "geocr:instrument": "AIA (Atmospheric Imaging Assembly) and HMI (Helioseismic and Magnetic Imager)"
            },
            "geocr:multiWavelengthConfiguration": {
                "@type": "geocr:MultiWavelengthConfiguration",
                "geocr:channelList": aia_wavelengths,
                "description": "AIA multi-wavelength EUV channels sampling different temperature regimes of the solar atmosphere"
            },
            "geocr:spectralBandMetadata": self._generate_spectral_band_metadata(),
            "distribution": [
                {
                    "@type": "cr:FileObject",
                    "@id": "data_repo",
                    "name": "data_repo",
                    "description": "SDO dataset directory containing NetCDF4 files",
                    "contentUrl": str(self.data_dir),
                    "encodingFormat": "local_directory",
                    "sha256": "placeholder_checksum_for_directory"
                },
                {
                    "@type": "cr:FileSet",
                    "@id": "nc-files",
                    "name": "nc-files",
                    "description": f"All NetCDF4 files containing synchronized AIA and HMI observations ({len(self.nc_files)} files)",
                    "containedIn": {"@id": "data_repo"},
                    "encodingFormat": "application/x-netcdf",
                    "includes": "*.nc"
                }
            ],
            "recordSet": [
                {
                    "@type": "cr:RecordSet",
                    "@id": "sdo_observations",
                    "name": "sdo_observations",
                    "description": f"SDO full-disk multi-instrument observations with {len([v for v in self.variables if v.startswith('aia')])} AIA EUV channels and {len([v for v in self.variables if v.startswith('hmi')])} HMI magnetic field variables",
                    "field": self._generate_field_definitions()
                }
            ]
        }
        
        # Write output
        with open(output_file, 'w') as f:
            json.dump(geocroissant, f, indent=2)
        
        self.sample_nc.close()
        return geocroissant
