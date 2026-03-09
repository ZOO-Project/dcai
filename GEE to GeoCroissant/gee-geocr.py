import ee
from google.oauth2 import service_account
from google.auth.transport.requests import Request
import json
import datetime

# 1. Authenticate to Earth Engine
SERVICE_ACCOUNT_FILE = "code-earthengine.json"
ASSET_ID = "COPERNICUS/S2/20170430T190351_20170430T190351_T10SEG"

creds = service_account.Credentials.from_service_account_file(
    SERVICE_ACCOUNT_FILE, scopes=["https://www.googleapis.com/auth/earthengine"]
)
creds.refresh(Request())
ee.Initialize(credentials=creds)
token = creds.token

# 2. Fetch asset metadata
meta = ee.data.getAsset(ASSET_ID)
props = meta["properties"]

# 3. Compute bounding box
coords = meta["geometry"]["coordinates"][0]
lons, lats = zip(*coords)
bbox = f"{min(lats)} {min(lons)} {max(lats)} {max(lons)}"
wkt = "POLYGON((" + ", ".join(f"{x} {y}" for x, y in coords) + "))"

# 4. Build per-band assets
assets = {}
for band in meta["bands"]:
    band_id = band["id"]
    res = band["grid"]["affineTransform"]["scaleX"]
    assets[band_id] = {
        "href": f"ee://{ASSET_ID}/{band_id}",
        "type": "image/tiff",
        "roles": ["data"],
        "band_name": band_id,
        "data_type": band["dataType"]["precision"].lower(),
        "spatial_resolution": res,
        "description": f"Sentinel-2 band {band_id} of image {ASSET_ID}",
    }

# 5. Convert bbox to array format
bbox_coords = [min(lons), min(lats), max(lons), max(lats)]

# 6. Build fileSet for the bands
fileSet_id = f"sentinel2-bands-{ASSET_ID.replace('/', '-')}"
band_files = []
for band in meta["bands"]:
    band_files.append(
        {
            "name": f"{band['id']}.tif",
            "path": f"ee://{ASSET_ID}/{band['id']}",
            "contentSize": (
                band.get("dimensions", {}).get("width", 0)
                * band.get("dimensions", {}).get("height", 0)
                * 2
            ),  # Approximate size
        }
    )

# 7. Assemble Geo-Croissant JSON-LD (using correct prefixes & geocr IRIs)
geo_croissant = {
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
        "transform": "cr:transform",
    },
    "@type": "sc:Dataset",
    "name": ASSET_ID.replace("/", "_"),
    "alternateName": [
        ASSET_ID.replace("/", "-"),
        f"Sentinel-2-{props.get('MGRS_TILE', '')}",
    ],
    "description": (
        f"Sentinel-2 Level-1C image over MGRS tile {props.get('MGRS_TILE','')} acquired"
        f" on {meta['startTime'][:10]}. This dataset contains"
        f" {len(meta['bands'])} spectral bands with spatial resolutions ranging from"
        " 10m to 60m."
    ),
    "conformsTo": [
        "http://mlcommons.org/croissant/1.0",
        "http://mlcommons.org/croissant/geo/1.0"
    ],
    "version": "1.0.0",
    "creator": {
        "@type": "Organization",
        "name": "European Space Agency (ESA)",
        "url": "https://www.esa.int/",
    },
    "url": f"https://earthengine.googleapis.com/v1alpha/projects/earthengine-public/assets/{ASSET_ID}",
    "keywords": [
        "Sentinel-2",
        "satellite imagery",
        "remote sensing",
        "multispectral",
        "Earth observation",
        f"MGRS-{props.get('MGRS_TILE', '')}",
        "Level-1C",
        "ESA",
        "Copernicus",
    ],
    "citeAs": f"https://earthengine.googleapis.com/v1alpha/projects/earthengine-public/assets/{ASSET_ID}",
    "datePublished": meta["startTime"][:10],
    "license": "https://creativecommons.org/licenses/by/4.0/",
    "spatialCoverage": {
        "@type": "Place",
        "geo": {
            "@type": "GeoShape",
            "box": f"{min(lats)} {min(lons)} {max(lats)} {max(lons)}"
        }
    },
    "temporalCoverage": f"{meta['startTime']}/{meta['endTime']}",
    "geocr:spatialResolution": "10-60m",
    "geocr:coordinateReferenceSystem": "EPSG:4326",
    "variableMeasured": [
        {
            "@type": "sc:PropertyValue",
            "sc:name": "Cloudy pixel percentage",
            "sc:value": props.get("CLOUDY_PIXEL_PERCENTAGE", 0),
        },
        {
            "@type": "sc:PropertyValue",
            "sc:name": "Cloud coverage assessment",
            "sc:value": props.get("CLOUD_COVERAGE_ASSESSMENT", 0),
        },
    ],
    "recordSet": [
        {
            "@type": "cr:RecordSet",
            "@id": f"sentinel2_bands_{ASSET_ID.replace('/', '_')}",
            "name": f"sentinel2_bands_{ASSET_ID.replace('/', '_')}",
            "description": f"Spectral bands for Sentinel-2 image {ASSET_ID}",
            "field": [
                {
                    "@type": "cr:Field",
                    "@id": f"{ASSET_ID.replace('/', '_')}/asset_id",
                    "name": f"{ASSET_ID.replace('/', '_')}/asset_id",
                    "description": "Asset identifier",
                    "dataType": "sc:Text",
                    "source": {
                        "fileSet": {"@id": fileSet_id},
                        "extract": {"fileProperty": "fullpath"},
                    },
                },
                {
                    "@type": "cr:Field",
                    "@id": f"{ASSET_ID.replace('/', '_')}/image_data",
                    "name": f"{ASSET_ID.replace('/', '_')}/image_data",
                    "description": "Satellite imagery data",
                    "dataType": "sc:ImageObject",
                    "source": {
                        "fileSet": {"@id": fileSet_id},
                        "extract": {"fileProperty": "fullpath"},
                    },
                    "geocr:bandConfiguration": {
                        "@type": "geocr:BandConfiguration",
                        "geocr:totalBands": len(meta["bands"]),
                        "geocr:bandNameList": [band["id"] for band in meta["bands"]],
                    },
                },
            ],
        }
    ],
    "fileSet": [
        {
            "@type": "cr:FileSet",
            "@id": fileSet_id,
            "name": fileSet_id,
            "description": f"Sentinel-2 spectral bands for {ASSET_ID}",
            "includes": "*.tif",
            "encodingFormat": "image/tiff",
            "fileObject": [
                {
                    "@type": "cr:FileObject",
                    "name": f"{band['id']}.tif",
                    "contentUrl": f"https://earthengine.googleapis.com/v1alpha/projects/earthengine-public/assets/{ASSET_ID}/{band['id']}",
                    "encodingFormat": "image/tiff",
                    "contentSize": (
                        band.get("dimensions", {}).get("width", 0)
                        * band.get("dimensions", {}).get("height", 0)
                        * 2
                    ),
                }
                for band in meta["bands"]
            ],
        }
    ],
    "distribution": [
        {
            "@type": "cr:FileObject",
            "@id": f"sentinel2-bands-{ASSET_ID.replace('/', '-')}",
            "name": f"Sentinel-2 Bands for {ASSET_ID}",
            "description": f"Downloadable Sentinel-2 spectral bands for {ASSET_ID}",
            "contentUrl": f"https://earthengine.googleapis.com/v1alpha/projects/earthengine-public/assets/{ASSET_ID}",
            "encodingFormat": "application/json",
        }
    ],
}

# 8. Write to file
with open("gee.json", "w") as f:
    json.dump(geo_croissant, f, indent=2)

print("Geo-Croissant JSON-LD saved to gee.json")