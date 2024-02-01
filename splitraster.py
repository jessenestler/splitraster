import argparse
import math
from pathlib import Path

import geopandas as gpd
from osgeo import gdal  # Import before rasterio to avoid DLL errors
import rasterio as rio
from rasterio.windows import Window
from shapely.geometry import box


class TileGenerator:
    def __init__(self, in_raster: Path, out_folder: Path,
                 tile_area_sqm: int, nodata: int, file_name: str):
        self.in_raster = in_raster
        self.out_folder = out_folder
        self.tile_area_sqm = tile_area_sqm
        self.nodata = nodata
        self.file_name = file_name
        self._raster = self.raster
        self._metadata = self.metadata
        self._res_x = self.res_x
        self._res_y = self.res_y
        self._tile_dims = self.tile_dims
        self._colormap = self.colormap

    @property
    def raster(self) -> rio.DatasetReader:
        return rio.open(self.in_raster, 'r')

    @property
    def metadata(self) -> dict:
        """Returns the metadata of the CHM, with nodata values filled in.

        Returns
        -------
        dict
            Metadata understandable by rasterio functions and methods
        """
        # Create a dict with each rasterio data type's nodata value
        null_values = {
            'uint8': 255,
            'uint16': 65535,
            'uint32': 4294967295,
            'int8': -128,
            'int16': -32768,
            'int32': -2147483648,
            'float32': -9999.0,
            'float64': -9999.0
        }

        # Save the metadata
        meta = self._raster.profile

        # Alter nodata if it's null
        if self.nodata:
            meta['nodata'] = self.nodata
        else:
            try:
                meta['nodata'] = null_values[meta['dtype']]
            except KeyError:
                pass

        return meta

    @property
    def res_x(self) -> int:
        """The CHM raster resolution in the x direction"""
        return self._raster.transform[0]

    @property
    def res_y(self) -> int:
        """The CHM raster resolution in the y direction"""
        return -self._raster.transform[4]

    @property
    def tile_dims(self) -> tuple:
        """The tile dimensions in number of pixels in the (x, y) directions."""
        # Obtain number of linear units per m
        factor = self._raster.crs.linear_units_factor[1]

        # Determine side length of a square with spec'd tile area
        side_length_m = math.sqrt(self.tile_area_sqm)

        # Perform unit conversion to get number of pixels
        tile_width = math.ceil(side_length_m / factor / self._res_x)
        tile_height = math.ceil(side_length_m / factor / self._res_y)

        return (tile_width, tile_height)

    @property
    def colormap(self) -> dict:
        """The color map for the raster."""
        try:
            color = self._raster.colormap(1)
        except ValueError:
            color = None
        return color

    def generate_tiles(self):
        """
        Generate tiles from the input raster based on the specified tile area.
        """
        # Define widths and heights
        tile_width, tile_height = self._tile_dims
        total_width = self._raster.width
        total_height = self._raster.height

        # Iterate through the raster using these widths and heights
        for r_num, row in enumerate(range(0, total_height, tile_height)):
            for c_num, col in enumerate(range(0, total_width, tile_width)):
                # Get data within the window of interest
                window = Window(col, row, tile_width, tile_height)
                tile_data = self._raster.read(window=window)

                # Skip creating the tile if it contains only NoData values
                meta = self._metadata
                if not tile_data[tile_data != meta['nodata']].any():
                    continue

                # Update metadata
                meta['width'], meta['height'] = tile_width, tile_height
                meta['transform'] = rio.windows.transform(
                    window, self._raster.transform)

                tile_name = f'R{r_num}C{c_num}{self.file_name}.tif'
                out_tile = self.out_folder / tile_name
                with rio.open(out_tile, 'w', **meta) as dst:
                    dst.write(tile_data)
                    if self._colormap:
                        dst.write_colormap(1, self._colormap)

    def create_tile_index(self):
        """
        Create a shapefile index for the generated tiles and export it.
        """
        # Initialize empty list of dataframes
        gdfs = []

        # Iterate through each tile
        for tile_file in self.out_folder.glob('*.tif'):
            # Grab the tile name
            name = tile_file.stem.split('_')[0]

            # Grab tile boundaries and create gdf
            with rio.open(tile_file) as src:
                xmin, ymin, xmax, ymax = src.bounds
                tile_geometry = box(xmin, ymin, xmax, ymax)
                gdf = gpd.GeoDataFrame(geometry=[tile_geometry],
                                       crs=src.crs)
                gdf['tile'] = name
                gdfs.append(gdf)

        # Concat all gdfs into one
        tiles_gdf = gpd.pd.concat(gdfs, ignore_index=True)

        # Export
        out_index = self.out_folder / f'index{self.file_name}.shp'
        tiles_gdf.to_file(out_index)


# Create parser with a description and arguments
description = (
    "Splits a large raster into a series of tiles, and generates a tile index."
)
parser = argparse.ArgumentParser(
    description=description
)

# Add arguments
parser.add_argument(
    "raster",
    type=Path,
    help="The raster that needs to be split into tiles."
)

parser.add_argument(
    "folder",
    type=Path,
    help="The output location for the tiles."
)

parser.add_argument(
    "-a",
    "--area",
    default=1000000,
    type=int,
    help=(
        "The area of the tiles in square meters. Defaults to 1,000,000 "
        "square meters (e.g. 1 square kilometer)."
    )
)

parser.add_argument(
    "-n",
    "--nodata",
    default=None,
    type=int,
    help=(
        "The desired no data value, if it differs from the raster's "
        "standard one."
    )
)

parser.add_argument(
    "-f",
    "--file-name",
    default='_lulc',
    type=str,
    help=(
        "The string to append to the end of the tile file names. The tool "
        "creates file names according to which row and column the tile "
        "belongs to, and this will append a string to the end of every file "
        "(e.g. R0C4_lulc.tif or R0C4_ppa.tif, etc). Defaults to '_lulc'."
    )
)

if __name__ == '__main__':
    # Parse arguments
    args = parser.parse_args()

    # Initialize a TileGenerator object
    tile_gen = TileGenerator(args.raster, args.folder,
                             args.area, args.nodata, args.file_name)

    # Generate tiles
    tile_gen.generate_tiles()

    # Create tile index
    tile_gen.create_tile_index()
